import numpy as np
import trimesh

from molfidget.atom import Atom
from molfidget.shape import Shape
from molfidget.config import BondConfig, DefaultBondConfig, ShapeConfig


class Bond:
    def __init__(self, config: BondConfig, default: DefaultBondConfig):
        self.atom1_name = config.atom_pair[0]
        self.atom2_name = config.atom_pair[1]
        self.atom_name = config.atom_pair
        
        self.shape_type = config.shape_type
        self.bond_type = config.bond_type if config.bond_type else default.bond_type
        self.bond_gap_mm = config.bond_gap_mm if config.bond_gap_mm is not None else default.bond_gap_mm
        self.bond_gap = self.bond_gap_mm / 10.0  # Convert mm to angstrom
        self.chamfer_length = config.chamfer_length if config.chamfer_length is not None else default.chamfer_length
        self.hole_length = config.hole_length if config.hole_length is not None else default.hole_length
        self.hole_radius = config.hole_radius if config.hole_radius is not None else default.hole_radius
        self.shaft_radius = config.shaft_radius if config.shaft_radius is not None else default.shaft_radius
        self.shaft_gap = config.shaft_gap if config.shaft_gap is not None else default.shaft_gap
        self.stopper_radius = config.stopper_radius if config.stopper_radius is not None else default.stopper_radius
        self.stopper_length = config.stopper_length if config.stopper_length is not None else default.stopper_length
        self.shaft_length = config.shaft_length if config.shaft_length is not None else default.shaft_length
        self.wall_thickness = config.wall_thickness if config.wall_thickness is not None else default.wall_thickness

        if config.bond_type == "single":
            config.shape_pair[0].shape_type = "spin"
            config.shape_pair[1].shape_type = "hole"
        elif config.bond_type == "double":
            config.shape_pair[0].shape_type = "fixed"
            config.shape_pair[1].shape_type = "hole"
        elif config.bond_type == "triple":
            config.shape_pair[0].shape_type = "fixed"
            config.shape_pair[1].shape_type = "hole"

        self.shape_pair = [Shape(self.atom1_name, config.shape_pair[0], default), Shape(self.atom2_name, config.shape_pair[1], default)]

    def update_atoms(self, atoms: dict):
        self.atom1 = atoms[self.atom1_name]
        self.atom2 = atoms[self.atom2_name]
        self.bond_distance = np.linalg.norm(np.array([self.atom2.x - self.atom1.x, self.atom2.y - self.atom1.y, self.atom2.z - self.atom1.z]))
        self.vector = np.array([self.atom2.x - self.atom1.x, self.atom2.y - self.atom1.y, self.atom2.z - self.atom1.z])
        self.atom_distance = np.linalg.norm(self.vector)
        self.vector /= np.linalg.norm(self.vector)

        self.shape_pair[0].update_atom(self.atom1, self.atom2)
        self.shape_pair[1].update_atom(self.atom2, self.atom1)

    def __repr__(self):
        return f"Bond({self.atom1.name}, {self.atom2.name})"

    def update_slice_distance(self):
        # Update the slice distance based on the configuration
        r1 = self.atom1.scale * self.atom1.radius
        r2 = self.atom2.scale * self.atom2.radius
        self.slice_distance1 = (r1**2 - r2**2 + self.atom_distance**2) / (2 * self.atom_distance)
        self.slice_distance2 = (r2**2 - r1**2 + self.atom_distance**2) / (2 * self.atom_distance)

    def sculpt_atoms2(self):
        print(f"Sculpting bond between {self.atom1.name} and {self.atom2.name}")
        self.slice_atoms_by_bond_plane()
        for shape in self.shape_pair:
            if shape.taper_radius_scale is not None and shape.taper_angle_deg is not None:
                shape.sculpt_trimesh_by_taper()
            if shape.shape_type == "spin":
                shape.sculpt_trimesh_by_spin()
            elif shape.shape_type == "fixed":
                shape.sculpt_trimesh_by_fixed()
            elif shape.shape_type == "hole":
                shape.sculpt_trimesh_by_hole()

    def sculpt_trimesh_model(self, mesh: trimesh.Trimesh):
        # mesh = self.slice_by_bond_plane(mesh)

        #if self.shaft_types[0] == "taper" or self.shaft_types[1] == "taper":
        #    return mesh
        #mesh = self.slice_by_taper(mesh)
        if True:
            if self.type == "single":
                # Create the cavity
                cavity = self.create_cavity_shape()
                cavity.apply_translation([0, 0, self.slice_distance])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                cavity.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.difference([mesh, cavity], check_volume=False)
                # Create the shaft
                shaft = self.create_rotate_shaft()
                shaft.apply_translation([0, 0, self.slice_distance])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                shaft.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.union([mesh, shaft], check_volume=False)
            else:
                # Create the fixed shaft
                shaft = self.create_fixed_shaft_shape()
                shaft.apply_translation([0, 0, self.slice_distance - self.bond_gap])
                rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
                shaft.apply_transform(rotation_matrix)
                mesh = trimesh.boolean.union([mesh, shaft], check_volume=False)
        else:
            hole = self.create_hole_shape()
            hole.apply_translation([0, 0, self.slice_distance])
            rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
            hole.apply_transform(rotation_matrix)
            mesh = trimesh.boolean.difference([mesh, hole], check_volume=False)
        return mesh
    
    def sculpt_trimesh_by_spin(self, mesh: trimesh.Trimesh):
        # Create the cavity
        cavity = self.create_cavity_shape()
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

    def slice_by_taper(self, mesh: trimesh.Trimesh):
        import math
        taper_distance = 0.3
        if taper_distance > self.slice_distance:
            taper_distance = self.slice_distance
        atom_radius = self.atom1.radius * self.atom1.scale
        # r2 は上円錐の半径
        r2 = math.sqrt(atom_radius**2 - self.slice_distance**2)
        # r1 は下円錐の半径
        r1 = math.sqrt(atom_radius**2 - (self.slice_distance -taper_distance)**2)
        # h1 下円錐の高さ
        h1 = r1 * taper_distance / (r1 - r2)
        # h 円錐の高さ
        h = h1 + (self.slice_distance - taper_distance) + atom_radius
        # r 大円錐の半径
        r = h * r1 / h1
        print(f"Creating taper cone with r1={r1}, r2={r2}, h1={h1}, h={h}, r={r}")
        cone = trimesh.creation.cone(radius=r, height=h)
        cone.apply_translation([0, 0, - atom_radius])
        z_axis = np.array([0, 0, 1])
        rotation_matrix = trimesh.geometry.align_vectors(z_axis, self.vector)
        cone.apply_transform(rotation_matrix)
        mesh = trimesh.boolean.intersection([mesh, cone], check_volume=False)
        # mesh = trimesh.boolean.union([mesh, cone], check_volume=False)

        return mesh

    def slice_atoms_by_bond_plane(self):
        for atom1, atom2 in [(self.atom1, self.atom2), (self.atom2, self.atom1)]:
            vector = np.array([atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z])
            distance = np.linalg.norm(vector)
            vector /= distance
            r1 = atom1.scale * atom1.radius
            r2 = atom2.scale * atom2.radius
            slice_distance = (r1**2 - r2**2 + distance**2) / (2 * distance)
            # atom1をbond平面でスライスする
            box = trimesh.primitives.Box(
                extents=[atom1.radius * 2, atom1.radius * 2, atom1.radius * 2]
            )
            box.apply_translation([0, 0, slice_distance - atom1.radius - self.bond_gap / 2])
            z_axis = np.array([0, 0, 1])
            rotation_matrix = trimesh.geometry.align_vectors(z_axis, vector)
            box.apply_transform(rotation_matrix)
            atom1.mesh = trimesh.boolean.intersection([atom1.mesh, box], check_volume=False)

    def create_rotate_shaft(self):
        # Create a shaft
        # d1: Shaft length including the wall thickness and gap without chamfer
        d1 = self.shaft_length + self.wall_thickness + self.shaft_gap
        cylinder1 = trimesh.creation.cylinder(radius=self.shaft_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1 / 2])
        # Create the chamfer on the shaft
        cylinder3 = trimesh.creation.cylinder(
            radius=self.shaft_radius, height=self.chamfer_length
        )
        cylinder3.apply_translation([0, 0, self.chamfer_length / 2])
        cone1 = trimesh.creation.cone(
            radius=self.shaft_radius, height=2 * self.shaft_radius, sections=32
        )
        cone1 = trimesh.boolean.intersection([cone1, cylinder3], check_volume=False)
        cylinder1 = trimesh.boolean.union([cylinder1, cone1], check_volume=False)
        cylinder1.apply_translation([0, 0, self.shaft_length - self.chamfer_length])
        # Create the stopper
        cylinder2 = trimesh.creation.cylinder(
            radius=self.stopper_radius, height=self.stopper_length
        )
        cylinder2.apply_translation(
            [0, 0, -self.stopper_length / 2 - self.wall_thickness - self.shaft_gap]
        )
        mesh = trimesh.boolean.union([cylinder1, cylinder2], check_volume=False)
        return mesh

    def create_cavity_shape(self):
        eps = 0.01  # Small epsilon to avoid numerical issues
        # Create the cavity shape for the shaft
        d1 = self.wall_thickness + eps
        cylinder1 = trimesh.creation.cylinder(
            radius=self.shaft_radius + self.shaft_gap, height=d1
        )
        cylinder1.apply_translation([0, 0, -d1 / 2 + eps])
        # Create the cavity for the stopper
        d2 = self.stopper_length + 2 * self.shaft_gap
        cylinder2 = trimesh.creation.cylinder(
            radius=self.stopper_radius + self.shaft_gap, height=d2
        )
        cylinder2.apply_translation([0, 0, -d2 / 2 - d1 + eps])
        mesh = trimesh.boolean.union([cylinder1, cylinder2], check_volume=False)
        return mesh

    def create_fixed_shaft_shape(self, config):
        eps = 0.01  # Small epsilon to avoid numerical issues
        # Create a fixed shaft shape
        d1 = self.shaft_length + self.bond_gap - self.chamfer_length
        cylinder1 = trimesh.creation.cylinder(radius=self.shaft_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1 / 2])
        # Create the chamfer on the shaft
        cylinder2 = trimesh.creation.cylinder(
            radius=self.shaft_radius, height=self.chamfer_length
        )
        cylinder2.apply_translation([0, 0, self.chamfer_length / 2])
        cone1 = trimesh.creation.cone(
            radius=self.shaft_radius, height=2 * self.shaft_radius, sections=32
        )
        cone1 = trimesh.boolean.intersection([cone1, cylinder2], check_volume=False)
        cylinder1 = trimesh.boolean.union([cylinder1, cone1], check_volume=False)
        cylinder1.apply_translation(
            [0, 0, self.shaft_length + self.bond_gap - self.chamfer_length - eps]
        )
        return cylinder1

    def create_hole_shape(self):
        # Create a hole shape for the shaft
        eps = 0.01  # Small epsilon to avoid numerical issues
        d1 = self.hole_length + eps
        cylinder1 = trimesh.creation.cylinder(radius=self.hole_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1 / 2 + eps])
        return cylinder1