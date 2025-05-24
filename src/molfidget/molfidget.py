import scipy
import trimesh
import numpy as np
from collections import OrderedDict as Orderdict

g_check_volume = False
g_vdw_scale = 0.8  # Scale factor for van der Waals radius

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
    ("C", "C"): (1.54, 1.33, 1.20),
    ("C", "O"): (1.43, 1.23, 1.17),
    ("C", "H"): (1.09, 0.00, 0.00),
    ("O", "O"): (1.48, 1.2075, 1.28),
    ("H", "O"): (0.96, 0.00, 0.00),
    ("S", "C"): (1.81, 1.50, 1.40),
    ("S", "O"): (1.70, 1.45, 1.35),
    ("S", "H"): (1.34, 0.00, 0.00),
    ("O", "H"): (0.96, 0.00, 0.00),
    ("N", "C"): (1.47, 1.34, 1.20),
    ("N", "O"): (1.40, 1.20, 1.15),
    ("N", "H"): (1.01, 0.00, 0.00),
    ("N", "N"): (1.45, 1.25, 1.10),
}

class Atom:
    def __init__(self, id, name, x, y, z):
        self.id = id
        self.name = name[0]
        if self.name not in atom_radius_table:
            raise ValueError(f"Unknown atom {id}: name: {self.name}")

        # ファンデルワールス半径のままだと干渉するのでちょっと小さくする
        self.radius = g_vdw_scale * atom_radius_table[self.name]
        self.x = x
        self.y = y
        self.z = z
        self.bonds = []

    def __repr__(self):
        return f"{self.name}_{self.id}({self.x}, {self.y}, {self.z})"
    
    def create_trimesh_model(self):
        # Trimeshモデルを生成するためのコードをここに記述
        mesh = trimesh.primitives.Sphere(radius=self.radius, center=[0, 0, 0])
        mesh.visual.vertex_colors = atom_color_table[self.name]
        # Scrupt the sphere to represent bonds
        for bond in self.bonds:
            mesh = bond.sculpt_trimesh_model(mesh)
        mesh.apply_translation([self.x, self.y, self.z])
        return mesh
    
class Bond:
    # Size of the shaft and the hole
    wall_thickness = 0.2
    shaft_r1 = 0.3
    shaft_r2 = 0.4
    shaft_r3 = 0.01
    shaft_d1 = 0.6
    shaft_d2 = 0.3
    shaft_c1 = 0.1 # Edge chanfer length [Angstrom]

    def __init__(self, atom1:Atom, atom2:Atom, type:str, shaft: bool, shaft_gap: float = 0.0, bond_gap: float = 0.0):
        '''
        結合を表すクラス
        atom1(Atom): 原子1(分割したとき軸側になる)
        atom2(Atom): 原子2(分割したとき穴側になる)
        type(str): 結合の種類(single, double, triple)
        shaft(bool): 軸側を可動にするかどうか
        shaft_gap(float): 軸と穴の隙間[Angstrom]
        '''
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type
        self.shaft = shaft
        self.shaft_gap = shaft_gap
        self.bond_gap = bond_gap
        self.vector = np.array([atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z])
        self.atom_distance = np.linalg.norm(np.array([atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z]))
        self.vector = self.vector / np.linalg.norm(self.vector)
        self.slice_distance = (atom1.radius**2 - atom2.radius**2 + self.atom_distance**2)/(2*self.atom_distance)
    
    def __repr__(self):
        return f"Bond({self.atom1.id}, {self.atom2}, type={self.type})"
    
    def sculpt_trimesh_model(self, mesh):
        mesh = self.slice_by_bond_plane(mesh)
        if self.shaft:
            if self.type and self.type == "single":
                # Create the cavity
                cavity = self.create_cavity_shape(self.shaft_r1+2*self.shaft_gap, self.shaft_d1, self.shaft_r2+2*self.shaft_gap, self.shaft_d2+2*self.shaft_gap)
                cavity.apply_translation([0, 0, self.slice_distance - self.wall_thickness])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                cavity.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.difference([mesh, cavity], check_volume=g_check_volume)
                # Create the shaft
                shaft = self.create_shaft_shape(self.shaft_r1, self.shaft_d1, self.shaft_r2, self.shaft_d2, self.shaft_r3, 2*self.shaft_gap, self.shaft_c1)
                shaft.apply_translation([0, 0, self.slice_distance - self.wall_thickness - self.shaft_gap])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                shaft.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.union([mesh, shaft], check_volume=g_check_volume)
            else:
                tool = trimesh.primitives.Cylinder(radius=0.3, height=0.5)
                tool.apply_translation([0, 0, self.slice_distance+0.2])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                tool.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.union([mesh, tool], check_volume=g_check_volume)
        else:
            tool = trimesh.primitives.Cylinder(radius=0.3, height=0.5)
            tool.apply_translation([0, 0, self.slice_distance-0.24])
            rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
            tool.apply_transform(rotation_matrix)
            mesh = trimesh.boolean.difference([mesh, tool], check_volume=g_check_volume)
        mesh.visual.vertex_colors = atom_color_table[self.atom1.name]
        return mesh

    def slice_by_bond_plane(self, mesh):
        box = trimesh.primitives.Box(extents=[self.atom1.radius*2, self.atom1.radius*2, self.atom1.radius*2])
        box.apply_translation([0, 0, self.slice_distance-self.atom1.radius-self.bond_gap])
        z_axis = np.array([0, 0, 1])
        rotation_matrix = trimesh.geometry.align_vectors(z_axis, self.vector)
        box.apply_transform(rotation_matrix)
        # create the slice
        slice = trimesh.boolean.intersection([mesh, box], check_volume=g_check_volume)
        return slice
    
    def create_shaft_shape(self, r1, d1, r2, d2, r3, d3, c1):
        # Subtract the chamfer length from the height of the first cylinder
        d1 = d1 - c1
        cylinder1 = trimesh.creation.cylinder(radius=r1, height=d1)
        # Create the chamfer on the shaft
        if c1 > 0:
            cone1 = trimesh.creation.cone(radius=r1, height=2*r1, sections=32)
            cylinder3 = trimesh.creation.cylinder(radius=r1, height=c1)
            cylinder3.apply_translation([0, 0, c1/2])
            cone1 = trimesh.boolean.intersection([cone1, cylinder3], check_volume=g_check_volume)
            cone1.apply_translation([0, 0, d1/2])
            cylinder1 = trimesh.boolean.union([cylinder1, cone1], check_volume=g_check_volume)
        cylinder1.apply_translation([0, 0, d1/2])
        cylinder2 = trimesh.creation.cylinder(radius=r2, height=d2)
        cylinder2.apply_translation([0, 0, -d2/2])
        cylinder3 = trimesh.creation.cylinder(radius=r3, height=d3)
        cylinder3.apply_translation([0, 0, -d2-d3/2])
        mesh =  trimesh.boolean.union([cylinder1, cylinder2, cylinder3], check_volume=g_check_volume)
        return mesh
    
    def create_cavity_shape(self, r1, d1, r2, d2):
        cylinder1 = trimesh.creation.cylinder(radius=r1, height=d1)
        cylinder1.apply_translation([0, 0, d1/2])
        cylinder2 = trimesh.creation.cylinder(radius=r2, height=d2)
        cylinder2.apply_translation([0, 0, -d2/2])
        mesh =  trimesh.boolean.union([cylinder1, cylinder2], check_volume=g_check_volume)
        return mesh
    
def atom_distance(atom1: Atom, atom2: Atom):
    # Calculate the distance between two atoms
    return float(np.linalg.norm(np.array([atom1.x - atom2.x, atom1.y - atom2.y, atom1.z - atom2.z])))

class Molecule:
    def __init__(self, scale=1.0, shaft_gap_mm=0.2, bond_gap_mm=0.0):
        # Dictionary to hold atoms by their ids
        self.atoms = Orderdict()
        # Dictionary to hold bonds by their ids
        self.bonds = {}
        # Scale of the molecule
        self.scale = scale
        # Gap between the shaft and the cavity [angstrom]
        self.shaft_gap = shaft_gap_mm / self.scale
        # Too much gap is not good for sculpting
        if self.shaft_gap > 0.05:
            self.shaft_gap = 0.05
        # Gap between the bond plane[angstrom]
        self.bond_gap = bond_gap_mm / self.scale
        # Too much gap is not good for sculpting
        if self.bond_gap > 0.05:
            self.bond_gap = 0.05

    def add_atom(self, atom):
        # Add an atom to the molecule
        self.atoms[atom.id] = atom

    def load_pdb_file(self, file_name):
        # Load a PDB file and populate the molecule with atoms and bonds
        with open(file_name, 'r') as file:
            for line in file:
                if line.startswith("COMPND"):
                    self.name = line[10:].strip()
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    id = int(line[6:11].strip())
                    name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    self.atoms[id] = Atom(id, name, x, y, z)
        self.update_bond()


    def update_bond(self):
        # Update the bonds based on the distance between atoms
        for id1, atom1 in self.atoms.items():
            for id2, atom2 in self.atoms.items():
                if id1 == id2:
                    continue
                if tuple(sorted((atom1.name, atom2.name))) not in bond_distance_table:
                    continue
                # Check if the atoms are within triple bond distance
                if atom_distance(atom1, atom2) < 1.05*bond_distance_table[tuple(sorted((atom1.name, atom2.name)))][2]:
                    atom1.bonds.append(Bond(atom1, atom2, type="triple", shaft=id1 < id2, shaft_gap=self.shaft_gap, bond_gap=self.bond_gap))
                # Check if the atoms are within double bond distance
                elif atom_distance(atom1, atom2) < 1.05*bond_distance_table[tuple(sorted((atom1.name, atom2.name)))][1]:
                    atom1.bonds.append(Bond(atom1, atom2, type="double", shaft=id1 < id2, shaft_gap=self.shaft_gap, bond_gap=self.bond_gap))
                # Check if the atoms are within single bond distance
                elif atom_distance(atom1, atom2) < 1.05*bond_distance_table[tuple(sorted((atom1.name, atom2.name)))][0]:
                    atom1.bonds.append(Bond(atom1, atom2, type="single", shaft=id1 < id2, shaft_gap=self.shaft_gap, bond_gap=self.bond_gap))  

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
    
    def create_trimesh_scene(self):
        # Create a trimesh model for the molecule
        scene = trimesh.Scene()
        for atom in self.atoms.values():
            mesh = atom.create_trimesh_model()
            scene.add_geometry(mesh)
        # center the scene
        scene.apply_translation(-self.center)
        scene.apply_scale(self.scale)
        return scene

    def save_stl_files(self):
        for atom in self.atoms.values():
            mesh = atom.create_trimesh_model()
            mesh.apply_scale(self.scale)
            mesh.export(f"{atom.name}_{atom.id}.stl")

def main():
    # command line interface
    import argparse
    parser = argparse.ArgumentParser(description="Molecule visualization and manipulation")
    parser.add_argument("pdb_file", type=str, help="PDB file to load")
    parser.add_argument("--scale", type=float, default=1.0, help="Scale of the molecule (default: %(default)s)")
    parser.add_argument('--rotate', nargs=2, action='append', type=int, metavar=('id1', 'id2'), help="Atom id pair to make a joint")
    parser.add_argument('--shaft-gap', type=float, default=0.2, help="Gap of the shaft and cavity [mm] (default: %(default)s)")
    parser.add_argument('--bond-gap', type=float, default=0.0, help="Gap of the bond plane [mm] (default: %(default)s)")
    parser.add_argument('--check-volume', type=bool, default=False, help="Check volume of the mesh (default: %(default)s)")
    parser.add_argument('--vdw-scale', type=float, default=0.8, help="Scale factor for van der Waals radius (default: %(default)s)")

    args = parser.parse_args()

    global g_check_volume, g_vdw_scale
    g_check_volume = args.check_volume
    g_vdw_scale = args.vdw_scale

    molecule = Molecule(scale=args.scale, shaft_gap_mm=args.shaft_gap, bond_gap_mm=args.bond_gap)
    molecule.load_pdb_file(args.pdb_file)
    scene = molecule.create_trimesh_scene()
    scene.show()
    # Save the molecule as STL files
    scene.export(f"molecule.stl")
    
    molecule.save_stl_files()
    print(f"Loaded {len(molecule.atoms)} atoms and {len(molecule.bonds)} bonds from {args.pdb_file}")

if __name__ == "__main__":
    main()