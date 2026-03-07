import os
from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
import trimesh
import yaml

from molfidget.atom import Atom
from molfidget.bond import Bond
from molfidget.config import (
    BondConfig,
    DefaultConfig,
    MoleculeConfig,
)
from molfidget.constants import atom_color_table, atom_radius_table

class Molecule:
    def __init__(self, config: MoleculeConfig, default: DefaultConfig):
        self.atom_groups = OrderedDict()

        print(f"Loading configuration for molecule: {config.name}")
        self.name = config.name
        self.scale = config.scale if config.scale is not None else 1.0
        # Load individual atom configurations
        self.atoms = {}
        for atom_config in config.atoms:
            name = atom_config.name
            self.atoms[name] = Atom(atom_config, default.atom)
        # 干渉する原子ペアにplaneボンドを自動追加
        self.bonds = {}
        atom_names = list(self.atoms.keys())
        for i in range(len(atom_names)):
            for j in range(i + 1, len(atom_names)):
                a1, a2 = self.atoms[atom_names[i]], self.atoms[atom_names[j]]
                distance = np.linalg.norm(np.array([a2.x - a1.x, a2.y - a1.y, a2.z - a1.z]))
                if distance < a1.shape_radius + a2.shape_radius:
                    bond_config = BondConfig(
                        atom_pair=[atom_names[i], atom_names[j]],
                        bond_type="plane",
                    )
                    bond = Bond(bond_config, default.bond, default.shape, self.scale)
                    bond.index = 0
                    self.bonds[atom_names[i], atom_names[j]] = bond
                    print(f"Auto plane bond: {atom_names[i]} - {atom_names[j]} (distance={distance:.3f}, r1+r2={a1.shape_radius + a2.shape_radius:.3f})")
        # Load individual bond configurations (自動生成のplaneボンドを上書き)
        for idx, bond_config in enumerate(config.bonds, start=1):
            atom1_name, atom2_name = bond_config.atom_pair
            if (atom1_name, atom2_name) in self.bonds:
                del self.bonds[atom1_name, atom2_name]
            if (atom2_name, atom1_name) in self.bonds:
                del self.bonds[atom2_name, atom1_name]
            bond = Bond(bond_config, default.bond, default.shape, self.scale)
            bond.index = idx
            self.bonds[atom1_name, atom2_name] = bond
        # Update atom pairs based on bonds
        for atom in self.atoms.values():
            atom.update_bonds(self.bonds)
        for bond in self.bonds.values():
            bond.update_atoms(self.atoms)

    def __repr__(self):
        return f"Molecule: {self.name}, {len(self.atoms)} atoms"

    def create_trimesh_scene(self):
        # Create a trimesh model for the molecule
        scene = trimesh.Scene()
        for atom in self.atoms.values():
            atom.create_trimesh_model()
        for bond in self.bonds.values():
            bond.sculpt_atoms()
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
        self.atom_groups.clear()
        counter = 0
        for atom in self.atoms.values():
            for pair in atom.pairs.values():
                if pair.bond_type in ("none", "plane"): # 近接原子に自動で割り当てられるplaneも除外する
                    continue
                if pair.atom1.elem != pair.atom2.elem:
                    continue
                # Find groups containing atom1 / atom2 independently.
                group1_name = next((name for name, g in self.atom_groups.items() if pair.atom1 in g), None)
                group2_name = next((name for name, g in self.atom_groups.items() if pair.atom2 in g), None)

                if group1_name is None and group2_name is None:
                    # Create a new group if not found
                    self.atom_groups[f"group_{counter}"] = set()
                    self.atom_groups[f"group_{counter}"].add(pair.atom1)
                    self.atom_groups[f"group_{counter}"].add(pair.atom2)
                    counter += 1
                    continue

                if group1_name is not None and group2_name is not None:
                    if group1_name != group2_name:
                        # Merge two existing groups when the pair bridges them.
                        self.atom_groups[group1_name].update(self.atom_groups[group2_name])
                        del self.atom_groups[group2_name]
                    self.atom_groups[group1_name].add(pair.atom1)
                    self.atom_groups[group1_name].add(pair.atom2)
                    continue

                target_name = group1_name if group1_name is not None else group2_name
                self.atom_groups[target_name].add(pair.atom1)
                self.atom_groups[target_name].add(pair.atom2)

        print(f"Merged atoms into {len(self.atom_groups)} groups")
        print("Groups:", self.atom_groups)

    def save_group_stl_files(self, scale, output_dir: str ='output'):
        os.makedirs(output_dir, exist_ok=True)
        for group_name, group in self.atom_groups.items():
            # Merge the atoms and save as a single file
            meshes = [atom.mesh.copy().apply_scale(scale) for atom in group]
            merged_mesh = trimesh.util.concatenate(meshes)
            atom_list = sorted(group, key=lambda a: a.id)
            base_name = f"{atom_list[0].elem}_" + "_".join([atom.id for atom in atom_list])
            # 32文字制限
            if len(base_name) > 32:
                base_name = base_name[:32]
            file_name = base_name + ".stl"
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
