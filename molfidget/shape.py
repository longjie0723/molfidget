import trimesh
from molfidget.atom import Atom
from molfidget.config import ShapeConfig, MolfidgetConfig
import numpy as np


class Shape:
    def __init__(self, atom_name: str, shape_config: ShapeConfig, molfidget_config: MolfidgetConfig, scale: float):
        default = molfidget_config.default

        self.atom_name = atom_name
        self.atom = None
        self.slice_distance = None
        self.slice_radius = None
        self.vector = None

        self.shape_type = shape_config.shape_type if shape_config.shape_type is not None else None
        # bond_gap_mmはShapeConfigにはあるがDefaultShapeConfigにはないので、デフォルト値を直接指定
        self.bond_gap_mm = shape_config.bond_gap_mm if shape_config.bond_gap_mm is not None else 0.0
        self.bond_gap = self.bond_gap_mm / scale  # Convert mm to angstrom
        self.chamfer_length = shape_config.chamfer_length if shape_config.chamfer_length is not None else default.shape.chamfer_length

        # shaft_length: mm単位が指定されていればそちらを優先
        if shape_config.shaft_length_mm is not None:
            self.shaft_length = shape_config.shaft_length_mm / scale  # mm to Angstrom
        elif default.shape.shaft_length_mm is not None:
            self.shaft_length = default.shape.shaft_length_mm / scale
        else:
            self.shaft_length = shape_config.shaft_length if shape_config.shaft_length is not None else default.shape.shaft_length

        # shaft_radius: mm単位が指定されていればそちらを優先
        if shape_config.shaft_radius_mm is not None:
            self.shaft_radius = shape_config.shaft_radius_mm / scale  # mm to Angstrom
        elif default.shape.shaft_radius_mm is not None:
            self.shaft_radius = default.shape.shaft_radius_mm / scale
        else:
            self.shaft_radius = shape_config.shaft_radius if shape_config.shaft_radius is not None else default.shape.shaft_radius

        # hole_length: mm単位が指定されていればそちらを優先
        if shape_config.hole_length_mm is not None:
            self.hole_length = shape_config.hole_length_mm / scale  # mm to Angstrom
        elif default.shape.hole_length_mm is not None:
            self.hole_length = default.shape.hole_length_mm / scale
        else:
            self.hole_length = shape_config.hole_length if shape_config.hole_length is not None else default.shape.hole_length

        # hole_radius: mm単位が指定されていればそちらを優先
        if shape_config.hole_radius_mm is not None:
            self.hole_radius = shape_config.hole_radius_mm / scale  # mm to Angstrom
        elif default.shape.hole_radius_mm is not None:
            self.hole_radius = default.shape.hole_radius_mm / scale
        else:
            self.hole_radius = shape_config.hole_radius if shape_config.hole_radius is not None else default.shape.hole_radius

        self.shaft_gap = shape_config.shaft_gap if shape_config.shaft_gap is not None else default.shape.shaft_gap
        self.stopper_radius = shape_config.stopper_radius if shape_config.stopper_radius is not None else default.shape.stopper_radius
        self.stopper_length = shape_config.stopper_length if shape_config.stopper_length is not None else default.shape.stopper_length
        self.wall_thickness = shape_config.wall_thickness if shape_config.wall_thickness is not None else default.shape.wall_thickness
        self.taper_radius_scale = shape_config.taper_radius_scale if shape_config.taper_radius_scale is not None else default.shape.taper_radius_scale
        self.taper_angle_deg = shape_config.taper_angle_deg if shape_config.taper_angle_deg is not None else default.shape.taper_angle_deg

        print(f"shape_type: {self.shape_type}")
        print(f"hole_length: {self.hole_length}")
        print(f"shaft_radius: {self.shaft_radius}")

    def __str__(self):
        return f"Shape(atom: {self.atom_name}, type: {self.shape_type})"

    def update_atom(self, atom1: Atom, atom2: Atom):
        self.atom = atom1
        self.bond_distance = np.linalg.norm(np.array([atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z]))
        self.vector = np.array([atom2.x - atom1.x, atom2.y - atom1.y, atom2.z - atom1.z])
        self.atom_distance = np.linalg.norm(self.vector)
        self.vector /= np.linalg.norm(self.vector)

        # Update the slice distance based on the configuration
        r1 = atom1.vdw_scale * atom1.radius
        r2 = atom2.vdw_scale * atom2.radius
        self.slice_distance = (r1**2 - r2**2 + self.atom_distance**2) / (2 * self.atom_distance)
        self.slice_radius = np.sqrt(r1**2 - self.slice_distance**2)
        #self.slice_distance = (r2**2 - r1**2 + self.atom_distance**2) / (2 * self.atom_distance)


    def sculpt_trimesh_by_spin(self):
        # Create the cavity
        cavity = self.create_cavity_shape()
        cavity.apply_translation([0, 0, self.slice_distance])
        rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
        cavity.apply_transform(rotation_matrix)
        self.atom.mesh = trimesh.boolean.difference([self.atom.mesh, cavity], check_volume=False)
        # Create the shaft
        shaft = self.create_rotate_shaft()
        shaft.apply_translation([0, 0, self.slice_distance])
        rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
        shaft.apply_transform(rotation_matrix)
        self.atom.mesh = trimesh.boolean.union([self.atom.mesh, shaft], check_volume=False)

    def sculpt_trimesh_by_shaft(self):
        """固定軸（Dカット無し）を作成"""
        shaft = self.create_shaft_shape()
        shaft.apply_translation([0, 0, self.slice_distance - self.bond_gap])
        rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
        shaft.apply_transform(rotation_matrix)
        self.atom.mesh = trimesh.boolean.union([self.atom.mesh, shaft], check_volume=False)

    def sculpt_trimesh_by_shaft_dcut(self):
        """固定軸（Dカット有り）を作成"""
        shaft = self.create_shaft_dcut_shape()
        shaft.apply_translation([0, 0, self.slice_distance - self.bond_gap])
        rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
        shaft.apply_transform(rotation_matrix)
        self.atom.mesh = trimesh.boolean.union([self.atom.mesh, shaft], check_volume=False)

    def sculpt_trimesh_by_hole(self):
        """丸穴（Dカット無し）を作成"""
        hole = self.create_hole_shape()
        hole.apply_translation([0, 0, self.slice_distance])
        rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
        hole.apply_transform(rotation_matrix)
        self.atom.mesh = trimesh.boolean.difference([self.atom.mesh, hole], check_volume=False)

    def sculpt_trimesh_by_hole_dcut(self):
        """Dカット穴を作成"""
        hole = self.create_hole_dcut_shape()
        hole.apply_translation([0, 0, self.slice_distance])
        rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
        hole.apply_transform(rotation_matrix)
        self.atom.mesh = trimesh.boolean.difference([self.atom.mesh, hole], check_volume=False)

    def sculpt_trimesh_by_taper(self):
        if self.taper_angle_deg <= 0.0:
            return
        taper = self.create_taper_shape()
        taper.apply_translation([0, 0, -self.atom.shape_radius])
        rotation_matrix = trimesh.geometry.align_vectors([0, 0, 1], self.vector)
        taper.apply_transform(rotation_matrix)
        self.atom.mesh = trimesh.boolean.intersection([taper, self.atom.mesh], check_volume=False)
        #self.atom.mesh = trimesh.boolean.union([self.atom.mesh, taper], check_volume=False)
    
    def create_rotate_shaft(self):
        """回転軸（Dカット無し、ストッパー付き）を作成"""
        # Create a shaft
        # d1: Shaft length including the wall thickness and gap without chamfer
        d1 = self.shaft_length + self.wall_thickness + self.shaft_gap - self.chamfer_length + self.bond_gap / 2
        cylinder1 = trimesh.creation.cylinder(radius=self.shaft_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1 / 2])
        # Create the chamfer on the shaft
        cylinder3 = trimesh.creation.cylinder(
            radius=self.shaft_radius, height=self.chamfer_length
        )
        cylinder3.apply_translation([0, 0, self.chamfer_length / 2])
        cone1 = trimesh.creation.cone(
            radius=self.shaft_radius, height=2*self.shaft_radius, sections=32
        )
        cone1 = trimesh.boolean.intersection([cone1, cylinder3], check_volume=False)
        cylinder1 = trimesh.boolean.union([cylinder1, cone1], check_volume=False)
        # 面取り後の位置に移動
        cylinder1.apply_translation(
            [0, 0, self.shaft_length - self.chamfer_length])
        # Create the stopper
        cylinder2 = trimesh.creation.cylinder(
           radius=self.stopper_radius, height=self.stopper_length
        )
        cylinder2.apply_translation(
            [0, 0, -self.stopper_length / 2 - self.wall_thickness - self.shaft_gap]
        )
        mesh = trimesh.boolean.union([cylinder1, cylinder2], check_volume=False)
        return mesh

    def create_rotate_shaft_dcut(self):
        """回転軸（Dカット有り、ストッパー付き）を作成"""
        # Create a shaft
        # d1: Shaft length including the wall thickness and gap without chamfer
        d1 = self.shaft_length + self.wall_thickness + self.shaft_gap - self.chamfer_length + self.bond_gap / 2
        cylinder1 = trimesh.creation.cylinder(radius=self.shaft_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1 / 2])
        # Create the chamfer on the shaft
        cylinder3 = trimesh.creation.cylinder(
            radius=self.shaft_radius, height=self.chamfer_length
        )
        cylinder3.apply_translation([0, 0, self.chamfer_length / 2])
        cone1 = trimesh.creation.cone(
            radius=self.shaft_radius, height=2*self.shaft_radius, sections=32
        )
        cone1 = trimesh.boolean.intersection([cone1, cylinder3], check_volume=False)
        cylinder1 = trimesh.boolean.union([cylinder1, cone1], check_volume=False)
        # D-cut the shaft
        d2 = d1 + self.chamfer_length
        box1 = trimesh.creation.box(
            extents=[2 * self.shaft_radius, 2 * self.shaft_radius, d2])
        box1.apply_translation([0, 0.3*self.shaft_radius, -d2 / 2 + self.chamfer_length])
        cylinder1 = trimesh.boolean.intersection(
            [cylinder1, box1], check_volume=False)
        cylinder1.apply_translation(
            [0, 0, self.shaft_length - self.chamfer_length])
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

    def create_shaft_shape(self):
        """固定軸（Dカット無し）を作成"""
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
        cylinder1.apply_translation([0, 0, d1])
        return cylinder1

    def create_shaft_dcut_shape(self):
        """固定軸（Dカット有り）を作成"""
        # まずDカット無し軸を作成
        cylinder1 = self.create_shaft_shape()
        # D-cut the shaft
        box1 = trimesh.creation.box(
            extents=[2 * self.shaft_radius, 2 * self.shaft_radius, self.shaft_length + self.bond_gap])
        box1.apply_translation([0, 0.3*self.shaft_radius, (self.shaft_length + self.bond_gap) / 2])
        cylinder1 = trimesh.boolean.intersection(
            [cylinder1, box1], check_volume=False)
        return cylinder1

    def create_hole_shape(self):
        """丸穴（Dカット無し）を作成"""
        d1 = self.hole_length
        cylinder1 = trimesh.creation.cylinder(radius=self.hole_radius, height=d1)
        cylinder1.apply_translation([0, 0, -d1 / 2])
        return cylinder1

    def create_hole_dcut_shape(self):
        """Dカット穴を作成"""
        d1 = self.hole_length
        cylinder1 = trimesh.creation.cylinder(radius=self.hole_radius, height=d1)
        # D-cut the hole
        box1 = trimesh.creation.box(
            extents=[2 * self.hole_radius, 2 * self.hole_radius, d1])
        box1.apply_translation([0, 0.3*self.hole_radius, 0])
        cylinder1 = trimesh.boolean.intersection(
           [cylinder1, box1], check_volume=False)
        cylinder1.apply_translation([0, 0, -d1 / 2])
        return cylinder1
    
    def create_taper_shape(self):
        r = self.slice_radius * self.taper_radius_scale
        a = self.taper_angle_deg * np.pi / 180.0
        d = r / np.cos(a)
        h = d * np.sin(a)
        H = h + self.slice_distance + self.atom.shape_radius
        R = r*H/h
        cone1 = trimesh.creation.cone(radius = R, height = H, sections = 32)
        return cone1
