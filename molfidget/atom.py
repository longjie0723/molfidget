from dataclasses import dataclass

import numpy as np
import trimesh

from molfidget.config import AtomConfig, DefaultAtomConfig
from molfidget.constants import atom_radius_table, atom_color_table

class Atom:
    def __init__(self, config: AtomConfig, default: DefaultAtomConfig):
        self.name = config.name
        self.elem, self.id = self.name.split('_')
        if self.elem not in atom_radius_table:
            raise ValueError(f"Unknown atom {self.elem}: name: {self.name}")
        self.radius = atom_radius_table[self.elem]
        self.scale = config.scale if config.scale is not None else default.scale
        self.shape_radius = self.scale * self.radius
        self.x, self.y, self.z = config.position
        self.position = config.position
        self.color = config.color if config.color is not None else atom_color_table[self.elem]
        print(f"self.color:", self.color)
        self.pairs = {}

    def update_bonds(self, bonds: dict):
        self.pairs = {}
        for bond in bonds.values():
            if bond.atom1_name == self.name or bond.atom2_name == self.name:
                self.pairs[(bond.atom1_name, bond.atom2_name)] = bond

    def __repr__(self):
        return f"{self.name}: ({self.x}, {self.y}, {self.z})"

    def create_trimesh_model(self):
        # Sphere mesh for the atom
        self.mesh = trimesh.primitives.Sphere(
            radius = self.shape_radius, center=[0, 0, 0]
        )
