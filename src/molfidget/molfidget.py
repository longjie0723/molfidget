import argparse
import dataclasses
import os
import trimesh
import yaml
import numpy as np
from collections import OrderedDict as Orderdict
from dataclasses import dataclass

atom_radius_table = {
    "C": 1.7,
    "O": 1.52,
    "H": 1.2,
    "S": 1.8,
    "N": 1.55,
}

atom_color_table = {
    "C": [100, 100, 100, 255],
    "O": [255, 0, 0, 255],
    "H": [255, 255, 255, 255],
    "S": [255, 255, 0, 255],
    "N": [0, 0, 255, 255],
}

# Bond distance table in Angstroms for single, double, and triple bonds
bond_distance_table = {
    frozenset(["C", "C"]): (1.54, 1.33, 1.20),
    frozenset(["C", "O"]): (1.43, 1.23, 1.17),
    frozenset(["C", "H"]): (1.09, 0.00, 0.00),
    frozenset(["O", "O"]): (1.48, 1.2075, 1.28),
    frozenset(["H", "O"]): (0.96, 0.00, 0.00),
    frozenset(["C", "S"]): (1.81, 1.50, 1.40),
    frozenset(["O", "S"]): (1.70, 1.45, 1.35),
    frozenset(["H", "S"]): (1.34, 0.00, 0.00),
    frozenset(["C", "N"]): (1.47, 1.34, 1.20),
    frozenset(["N", "O"]): (1.40, 1.20, 1.15),
    frozenset(["H", "N"]): (1.01, 0.00, 0.00),
    frozenset(["N", "N"]): (1.45, 1.25, 1.10),
}

@dataclass
class ShapeConfig:
    scale: float = 10.0 # Scale factor for the whole model
    vdw_scale: float = 0.8 # Scale factor for van der Waals radius
    shaft_radius: float = 0.3 # Radius of the shaft [Angstrom]
    shaft_length: float = 0.3 # Length of the shaft [Angstrom]
    stopper_radius: float = 0.4 # Radius of the stopper [Angstrom]
    stopper_length: float = 0.2 # Length of the stopper [Angstrom]
    hole_radius: float = 0.3 # Radius of the hole [Angstrom]
    hole_length: float = 0.3 # Length of the hole [Angstrom]
    chamfer_length: float = 0.1 # Length of the chamfer [Angstrom]
    wall_thickness: float = 0.1 # Thickness of the wall [Angstrom]
    shaft_gap: float = 0.03 # Gap between the shaft and the hole [Angstrom]
    bond_gap: float = 0.0 # Gap between the bond plane [Angstrom]

    slice_distance: float = 0.0 # Distance from the atom to the bond plane, calculated in Bond class

class Atom:
    def __init__(self, id: int, name: str, x: float, y: float, z: float):
        self.id = id
        self.name = name[0]
        if self.name not in atom_radius_table:
            raise ValueError(f"Unknown atom {id}: name: {self.name}")
        self.radius = atom_radius_table[self.name]
        self.x = x
        self.y = y
        self.z = z
        self.pairs = {}

    def __repr__(self):
        return f"{self.name}_{self.id}({self.x}, {self.y}, {self.z})"
    
    def create_trimesh_model(self, config: ShapeConfig):
        # Shpere mesh for the atom
        mesh = trimesh.primitives.Sphere(radius=config.vdw_scale*self.radius, center=[0, 0, 0])
        # Scrupt the sphere to represent bonds
        for pair in self.pairs.values():
            pair.update_slice_distance(config)
            mesh = pair.sculpt_trimesh_model(mesh, config)
        mesh.apply_translation([self.x, self.y, self.z])
        mesh.visual.vertex_colors = atom_color_table[self.name]
        return mesh

class Bond:
    def __init__(self, atom1:Atom, atom2:Atom, type:str, shaft: bool):
        '''
        結合を表すクラス
        atom1(Atom): 原子1(分割したとき軸側になる)
        atom2(Atom): 原子2(分割したとき穴側になる)
        type(str): 結合の種類(single, double, triple)
        shaft(bool): True: 軸側 / False 穴側
        '''
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type
        self.shaft = shaft
        self.atom_distance = np.linalg.norm(np.array([atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z]))
        self.vector = np.array([atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z])
        self.vector = self.vector / np.linalg.norm(self.vector)
    
    def __repr__(self):
        return f"Bond({self.atom1.id}, {self.atom2}, type={self.type})"
    
    def update_slice_distance(self, config: ShapeConfig):
        # Update the slice distance based on the configuration
        r1 = config.vdw_scale * self.atom1.radius
        r2 = config.vdw_scale * self.atom2.radius
        self.slice_distance = (r1**2 - r2**2 + self.atom_distance**2)/(2*self.atom_distance)
    
    def sculpt_trimesh_model(self, mesh: trimesh.Trimesh, config: ShapeConfig):
        mesh = self.slice_by_bond_plane(mesh, config)

        if self.type == "none":
            return mesh
        
        if self.shaft:
            if self.type == "single":
                # Create the cavity
                cavity = self.create_cavity_shape(config)
                cavity.apply_translation([0, 0, self.slice_distance])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                cavity.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.difference([mesh, cavity], check_volume=False)
                # Create the shaft
                shaft = self.create_rotate_shaft(config)
                shaft.apply_translation([0, 0, self.slice_distance])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                shaft.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.union([mesh, shaft], check_volume=False)
            else:
                # Create the fixed shaft
                shaft = self.create_fixed_shaft_shape(config)
                shaft.apply_translation([0, 0, self.slice_distance - config.bond_gap])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                shaft.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.union([mesh, shaft], check_volume=False)
        else:
            hole = self.create_hole_shape(config)
            hole.apply_translation([0, 0, self.slice_distance])
            rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
            hole.apply_transform(rotation_matrix)
            mesh = trimesh.boolean.difference([mesh, hole], check_volume=False)
        return mesh

    def slice_by_bond_plane(self, mesh: trimesh.Trimesh, config: ShapeConfig):
        box = trimesh.primitives.Box(extents=[self.atom1.radius*2, self.atom1.radius*2, self.atom1.radius*2])
        box.apply_translation([0, 0, self.slice_distance-self.atom1.radius-config.bond_gap])
        z_axis = np.array([0, 0, 1])
        rotation_matrix = trimesh.geometry.align_vectors(z_axis, self.vector)
        box.apply_transform(rotation_matrix)
        mesh = trimesh.boolean.intersection([mesh, box], check_volume=False)
        #mesh = trimesh.boolean.union([mesh, box], check_volume=False)
        return mesh
    
    def create_rotate_shaft(self, config):
        # Create a shaft
        # d1: Shaft length including the wall thickness and gap without chamfer
        d1 = config.shaft_length + config.wall_thickness + config.shaft_gap
        cylinder1 = trimesh.creation.cylinder(radius=config.shaft_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1/2])
        # Create the chamfer on the shaft
        cylinder3 = trimesh.creation.cylinder(radius=config.shaft_radius, height=config.chamfer_length)
        cylinder3.apply_translation([0, 0, config.chamfer_length/2])
        cone1 = trimesh.creation.cone(radius=config.shaft_radius, height=2*config.shaft_radius, sections=32)
        cone1 = trimesh.boolean.intersection([cone1, cylinder3], check_volume=False)
        cylinder1 = trimesh.boolean.union([cylinder1, cone1], check_volume=False)
        cylinder1.apply_translation([0, 0, config.shaft_length - config.chamfer_length])
        # Create the stopper
        cylinder2 = trimesh.creation.cylinder(radius=config.stopper_radius, height=config.stopper_length)
        cylinder2.apply_translation([0, 0, - config.stopper_length/2 - config.wall_thickness - config.shaft_gap])
        mesh =  trimesh.boolean.union([cylinder1, cylinder2], check_volume=False)
        return mesh
    
    def create_cavity_shape(self, config):
        eps = 0.01  # Small epsilon to avoid numerical issues
        # Create the cavity shape for the shaft
        d1 = config.wall_thickness + eps
        cylinder1 = trimesh.creation.cylinder(radius=config.shaft_radius+config.shaft_gap, height=d1)
        cylinder1.apply_translation([0, 0, -d1/2 + eps])
        # Create the cavity for the stopper
        d2 = config.stopper_length + 2*config.shaft_gap
        cylinder2 = trimesh.creation.cylinder(radius=config.stopper_radius+config.shaft_gap, height=d2)
        cylinder2.apply_translation([0, 0, -d2/2 - d1 + eps])
        mesh =  trimesh.boolean.union([cylinder1, cylinder2], check_volume=False)
        return mesh
    
    def create_fixed_shaft_shape(self, config):
        eps = 0.01  # Small epsilon to avoid numerical issues
        # Create a fixed shaft shape
        d1 = config.shaft_length + config.bond_gap - config.chamfer_length
        cylinder1 = trimesh.creation.cylinder(radius=config.shaft_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1/2])
        # Create the chamfer on the shaft
        cylinder2 = trimesh.creation.cylinder(radius=config.shaft_radius, height=config.chamfer_length)
        cylinder2.apply_translation([0, 0, config.chamfer_length/2])
        cone1 = trimesh.creation.cone(radius=config.shaft_radius, height=2*config.shaft_radius, sections=32)
        cone1 = trimesh.boolean.intersection([cone1, cylinder2], check_volume=False)
        cylinder1 = trimesh.boolean.union([cylinder1, cone1], check_volume=False)
        cylinder1.apply_translation([0, 0, config.shaft_length + config.bond_gap - config.chamfer_length - eps])
        return cylinder1
    
    def create_hole_shape(self, config):
        # Create a hole shape for the shaft
        eps = 0.01  # Small epsilon to avoid numerical issues
        d1 = config.hole_length + eps
        cylinder1 = trimesh.creation.cylinder(radius=config.hole_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1/2 + eps])
        return cylinder1

def atom_distance(atom1: Atom, atom2: Atom):
    # Calculate the distance between two atoms
    return float(np.linalg.norm(np.array([atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z])))

class Molecule:
    def __init__(self):
        # Dictionary to hold atoms by their ids
        self.atoms = Orderdict()
        # Group of atoms that can be merged
        self.atom_groups = Orderdict()

    def load_mol_file(self, file_name):
        # Load a MOL file and populate the molecule with atoms and bonds
        with open(file_name, 'r') as file:
            lines = file.readlines()
        # Fist line contains the name of the molecule
        self.name = lines[0].strip()
        counts = lines[3].split()
        atom_count = int(counts[0])
        bond_count = int(counts[1])
        # Parse the atom lines
        for i in range(atom_count):
            data = lines[4 + i].strip().split()
            x = float(data[0])
            y = float(data[1])
            z = float(data[2])
            name = data[3]
            self.atoms[i+1] = Atom(i+1, name, x, y, z)
            print(f"Loaded atom: {self.atoms[i+1]}")
        # Create the pairs of atoms based on the distances
        self.create_pairs()
        # Parse the bond lines
        for i in range(bond_count):
            data = lines[4 + atom_count + i].strip().split()
            id1 = int(data[0])
            id2 = int(data[1])
            type = int(data[2])
            if type == 1:
                bond_type = "single"
            elif type == 2:
                bond_type = "double"
            elif type == 3:
                bond_type = "triple"
            elif type == 4:
                bond_type = "1.5"
            else:
                raise ValueError(f"Unknown bond type: {type}")
            self.atoms[id1].pairs[id1, id2] = Bond(self.atoms[id1], self.atoms[id2], type=bond_type, shaft=id1 < id2)
            self.atoms[id2].pairs[id1, id2] = Bond(self.atoms[id2], self.atoms[id1], type=bond_type, shaft=id2 < id1)
            print(f"Loaded bond: {self.atoms[id1].name}_{id1} - {self.atoms[id2].name}_{id2}, type={bond_type}")

    def load_pdb_file(self, file_name):
        # Load a PDB file and populate the molecule with atoms and bonds
        with open(file_name, 'r') as file:
            for line in file:
                # Do we have the name of the molecule somewhere?
                if line.startswith("COMPND"):
                    self.name = line[10:].strip()
                # Parse ATOM and HETATM lines to extract atom information
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    id = int(line[6:11].strip())
                    name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    self.atoms[id] = Atom(id, name, x, y, z)
        # Create the atom pairs based on the atom distances
        self.create_pairs()
        # Find the bonds based on the atom pairs
        self.find_bonds()

    def create_pairs(self):
        # Update the pairs of atoms based on the distance between them
        for id1, atom1 in self.atoms.items():
            for id2, atom2 in self.atoms.items():
                if id1 == id2:
                    continue
                if atom_distance(atom1, atom2) < atom1.radius + atom2.radius:
                    atom1.pairs[id1, id2] = Bond(atom1, atom2, type="none", shaft=id1 < id2)

    def find_bonds(self):
        # Update the bonds based on the distance between atoms
        for id1, atom1 in self.atoms.items():
            for id2, atom2 in self.atoms.items():
                if id1 == id2:
                    continue
                if frozenset([atom1.name, atom2.name]) not in bond_distance_table:
                    continue
                # Check if the atoms are within triple bond distance
                if atom_distance(atom1, atom2) < 1.05*bond_distance_table[frozenset([atom1.name, atom2.name])][2]:
                    atom1.pairs[id1, id2] = Bond(atom1, atom2, type="triple", shaft=id1 < id2)
                # Check if the atoms are within double bond distance
                elif atom_distance(atom1, atom2) < 1.05*bond_distance_table[frozenset([atom1.name, atom2.name])][1]:
                    atom1.pairs[id1, id2] = Bond(atom1, atom2, type="double", shaft=id1 < id2)
                # Check if the atoms are within single bond distance
                elif atom_distance(atom1, atom2) < 1.05*bond_distance_table[frozenset([atom1.name, atom2.name])][0]:
                    atom1.pairs[id1, id2] = Bond(atom1, atom2, type="single", shaft=id1 < id2)

    def __repr__(self):
        return f"Molecule({self.name}, {len(self.atoms)} atoms)"
    
    # acces the center of the molecule
    @property
    def center(self):
        # Get the center of the molecule
        x = np.mean([atom.x for atom in self.atoms.values()])
        y = np.mean([atom.y for atom in self.atoms.values()])
        z = np.mean([atom.z for atom in self.atoms.values()])
        return np.array([x, y, z])
    
    def create_trimesh_scene(self, config: ShapeConfig):
        # Create a trimesh model for the molecule
        scene = trimesh.Scene()
        for atom in self.atoms.values():
            mesh = atom.create_trimesh_model(config)
            scene.add_geometry(mesh)
        # center the scene
        scene.apply_translation(-self.center)
        scene.apply_scale(config.scale)
        return scene

    def save_stl_files(self, config: ShapeConfig, output_dir: str ='output'):
        os.makedirs(output_dir, exist_ok=True)
        for atom in self.atoms.values():
            mesh = atom.create_trimesh_model(config)
            mesh.apply_scale(config.scale)
            mesh.export(os.path.join(output_dir, f"{atom.name}_{atom.id}.stl"))

    def merge_atoms(self):
        counter = 0
        for atom in self.atoms.values():
            for pair in atom.pairs.values():
                if pair.type == "none":
                    continue
                if pair.atom1.name != pair.atom2.name:
                    continue
                # Search group containing atom1 or atom2
                group = next((g for g in self.atom_groups.values() if pair.atom1.id in g or pair.atom2.id in g), None)
                if group is None:
                    # Create a new group if not found
                    self.atom_groups[f"group_{counter}"] = set()
                    self.atom_groups[f"group_{counter}"].add(pair.atom1.id)
                    self.atom_groups[f"group_{counter}"].add(pair.atom2.id)
                    counter += 1
                else:
                    # Add the atoms to the existing group
                    group.add(pair.atom1.id)
                    group.add(pair.atom2.id)

        print(f"Merged atoms into {len(self.atom_groups)} groups")
        print("Groups:", self.atom_groups)

    def save_group_stl_files(self, config: ShapeConfig, output_dir: str ='output'):
        os.makedirs(output_dir, exist_ok=True)
        for group_name, group in self.atom_groups.items():
            # Merge the atoms and save as a single file
            meshes = [self.atoms[id].create_trimesh_model(config) for id in group]
            merged_mesh = trimesh.util.concatenate(meshes)
            merged_mesh.apply_scale(config.scale)
            merged_mesh.export(os.path.join(output_dir, f"{group_name}.stl"))

def main():
    config = ShapeConfig()
    # command line interface
    parser = argparse.ArgumentParser(description="Molecule visualization and manipulation")
    parser.add_argument("file_name", type=str, help="PDB or MOL file to load")
    parser.add_argument('--config-file', type=str, help="Configuration yaml file")
    parser.add_argument("--scale", type=float, default=config.scale, help="Scale of the molecule (default: %(default)s)")
    parser.add_argument('--vdw-radius-scale', type=float, default=config.vdw_scale, help="Scale factor for van der Waals radius (default: %(default)s)")
    parser.add_argument('--shaft-radius', type=float, default=config.shaft_radius, help="Radius of the shaft [Angstrom] (default: %(default)s)")
    parser.add_argument('--shaft-length', type=float, default=config.shaft_length, help="Length of the shaft [Angstrom] (default: %(default)s)")
    parser.add_argument('--stopper-radius', type=float, default=config.stopper_radius, help="Radius of the stopper [Angstrom] (default: %(default)s)")
    parser.add_argument('--stopper-length', type=float, default=config.stopper_length, help="Length of the stopper [Angstrom] (default: %(default)s)")
    parser.add_argument('--hole-radius', type=float, default=config.hole_radius, help="Radius of the hole [Angstrom] (default: %(default)s)")
    parser.add_argument('--hole-length', type=float, default=config.hole_length, help="Length of the hole [Angstrom] (default: %(default)s)")
    parser.add_argument('--chamfer-length', type=float, default=config.chamfer_length, help="Length of the chamfer [Angstrom] (default: %(default)s)")
    parser.add_argument('--wall-thickness', type=float, default=config.wall_thickness, help="Thickness of the wall [Angstrom] (default: %(default)s)")
    parser.add_argument('--shaft-gap', type=float, default=0.35, help="Gap between the shaft and the hole [mm] (default: %(default)s)")
    parser.add_argument('--bond-gap', type=float, default=0.0, help="Gap between the bond plane [mm] (default: %(default)s)") 
    parser.add_argument('--output-dir', type=str, default='output', help="Output directory for STL files (default: %(default)s)")

    args = parser.parse_args()

    # Load configuration from a YAML file if provided
    if args.config_file:
        with open(args.config_file, 'r') as file:
            config_data = yaml.safe_load(file)
            config = ShapeConfig(**config_data)

    config.scale = args.scale
    config.vdw_scale = args.vdw_radius_scale
    config.shaft_radius = args.shaft_radius
    config.shaft_length = args.shaft_length
    config.stopper_radius = args.stopper_radius
    config.stopper_length = args.stopper_length
    config.hole_radius = args.hole_radius
    config.hole_length = args.hole_length
    config.chamfer_length = args.chamfer_length
    config.wall_thickness = args.wall_thickness
    config.shaft_gap = min(0.05, args.shaft_gap / args.scale)
    config.bond_gap = min(0.05, args.bond_gap / args.scale)

    molecule = Molecule()

    if args.file_name.endswith(".pdb"):
        molecule.load_pdb_file(args.file_name)
    elif args.file_name.endswith(".mol"):
        molecule.load_mol_file(args.file_name)
    else:
        exit(f"Unsupported file format: {args.file_name}. Please provide a .pdb or .mol file.")

    scene = molecule.create_trimesh_scene(config)
    scene.show()
    # Save the molecule as STL files
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    scene.export(os.path.join(args.output_dir, "molecule.stl"))
    scene.export(os.path.join(args.output_dir, "molecule.ply"))
    
    molecule.save_stl_files(config, output_dir=args.output_dir)
    print(f"Loaded {len(molecule.atoms)} atoms from {args.file_name}")

    molecule.merge_atoms()
    molecule.save_group_stl_files(config, output_dir=args.output_dir)

    # Save the configuration to a YAML file
    config_data = dataclasses.asdict(config)
    config_data["file_name"] = args.file_name
    with open(os.path.join(args.output_dir, "config.yaml"), 'w') as file:
        yaml.dump(config_data, file, default_flow_style=False)

if __name__ == "__main__":
    main()