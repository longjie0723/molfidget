import os
from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
import trimesh
import yaml

from molfidget.atom import Atom
from molfidget.bond import Bond
from molfidget.config import (
    AtomConfig,
    BondConfig,
    DefaultAtomConfig,
    DefaultBondConfig,
    MoleculeConfig,
)
from molfidget.constants import atom_color_table, atom_radius_table

class Molecule:
    def __init__(self, config: MoleculeConfig):
        # Group of atoms that can be merged
        self.atom_groups = OrderedDict()

        # Load the configuration for the molecule
        print(f"Loading configuration for molecule: {config.name}")
        self.name = config.name
        self.scale = config.scale if config.scale is not None else 1.0
        # Load individual atom configurations
        self.atoms = {}
        for atom_config in config.atoms:
            name = atom_config.name
            self.atoms[name] = Atom(atom_config, config.default.atom)
        # Load individual bond configurations
        self.bonds = {}
        for bond_config in config.bonds:
            atom1_name, atom2_name = bond_config.atom_pair
            self.bonds[atom1_name, atom2_name] = Bond(bond_config, config.default.bond)
        # Update atom pairs based on bonds
        for atom in self.atoms.values():
            atom.update_bonds(self.bonds)
        for bond in self.bonds.values():
            bond.update_atoms(self.atoms)

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
        return f"Molecule: {self.name}, {len(self.atoms)} atoms"

    # Access the center of the molecule
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
            atom.create_trimesh_model()
        for bond in self.bonds.values():
            bond.sculpt_atoms2()
        for atom in self.atoms.values():
            atom.mesh.apply_translation([atom.x, atom.y, atom.z])
            atom.mesh.visual.vertex_colors = atom.color
            atom.mesh.visual.face_colors = atom.color
            scene.add_geometry(atom.mesh, geom_name=f"{atom.name}")
        return scene

    def save_stl_files(self, scale, output_dir: str = "output"):
        for atom in self.atoms.values():
            mesh = atom.mesh.copy()
            mesh.apply_scale(scale)
            mesh.export(os.path.join(output_dir, f"{atom.name}.stl"))    

    def merge_atoms(self):
        counter = 0
        for atom in self.atoms.values():
            for pair in atom.pairs.values():
                if pair.bond_type == "none":
                    continue
                if pair.atom1.elem != pair.atom2.elem:
                    continue
                # Search group containing atom1 or atom2
                group = next((g for g in self.atom_groups.values() if pair.atom1 in g or pair.atom2 in g), None)
                if group is None:
                    # Create a new group if not found
                    self.atom_groups[f"group_{counter}"] = set()
                    self.atom_groups[f"group_{counter}"].add(pair.atom1)
                    self.atom_groups[f"group_{counter}"].add(pair.atom2)
                    counter += 1
                else:
                    # Add the atoms to the existing group
                    group.add(pair.atom1)
                    group.add(pair.atom2)

        print(f"Merged atoms into {len(self.atom_groups)} groups")
        print("Groups:", self.atom_groups)

    def save_group_stl_files(self, scale, output_dir: str ='output'):
        os.makedirs(output_dir, exist_ok=True)
        for group_name, group in self.atom_groups.items():
            # Merge the atoms and save as a single file
            meshes = [atom.mesh.copy().apply_scale(scale) for atom in group]
            merged_mesh = trimesh.util.concatenate(meshes)
            atom_list = sorted(group, key=lambda a: a.id)
            file_name = f"{atom_list[0].elem}_" + "_".join([atom.id for atom in atom_list]) + ".stl"
            merged_mesh.export(os.path.join(output_dir, file_name))


def load_molfidget_file(file_path: str) -> MoleculeConfig:
    import yaml
    from yaml.loader import SafeLoader

    with open(file_path, 'r') as file:
        data = yaml.load(file, Loader=SafeLoader)

    molecule_config = MoleculeConfig(**data.get('molecule', {}))
    molecule_config.default_atom = DefaultAtomConfig(**data['molecule']['default_atom'])
    molecule_config.default_bond = DefaultBondConfig(**data['molecule']['default_bond'])
    molecule_config.atoms = [AtomConfig(**atom) for atom in data['molecule'].get('atoms', [])]
    molecule_config.bonds = [BondConfig(**bond) for bond in data['molecule'].get('bonds', [])]

    return molecule_config
