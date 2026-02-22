from collections import OrderedDict
from dataclasses import asdict, dataclass, field

from dacite import from_dict
from ruamel.yaml import YAML, CommentedMap, CommentedSeq
from typing import List
import os

from Bio.PDB import PDBParser, NeighborSearch

from molfidget.constants import bond_distance_table, hydrogen_bond_params
import numpy as np


@dataclass
class DefaultAtomConfig:
    scale: float = 0.8  # Scale factor for van der Waals radius
    color: List[int] = field(default_factory=lambda: [200, 200, 200, 255])


@dataclass
class DefaultBondConfig:
    bond_gap_mm: float = 0.1  # Gap between the bond plane [mm]
    bond_gap: float = 0.0  # Gap between the bond plane [Angstrom]
    shaft_radius: float = 0.3  # Radius of the shaft [Angstrom]
    shaft_length: float = 0.3  # Length of the shaft [Angstrom]
    shaft_radius_mm: float = None  # Radius of the shaft [mm] (priority over shaft_radius)
    shaft_length_mm: float = None  # Length of the shaft [mm] (priority over shaft_length)
    stopper_radius: float = 0.4  # Radius of the stopper [Angstrom]
    stopper_length: float = 0.2  # Length of the stopper [Angstrom]
    hole_radius: float = 0.3  # Radius of the hole [Angstrom]
    hole_length: float = 0.3  # Length of the hole [Angstrom]
    hole_radius_mm: float = None  # Radius of the hole [mm] (priority over hole_radius)
    hole_length_mm: float = None  # Length of the hole [mm] (priority over hole_length)
    chamfer_length: float = 0.1  # Length of the chamfer [Angstrom]
    wall_thickness: float = 0.1  # Thickness of the wall [Angstrom]
    shaft_gap: float = 0.03  # Gap between the shaft and the cavity [Angstrom]
    taper_angle_deg: float = 0.0 # Taper angle in degrees
    taper_radius_scale: float = 1.0 # Scale factor for the taper radius
    magnetic_hole_radius_mm: float = 3.525  # Radius of the magnetic hole [mm]
    magnetic_hole_length_mm: float = 2.0  # Length of the magnetic hole [mm]
    bond_marker: str = "hetero-only-except-ch" # Option for bond number marker display (e.g., on, off, hetero-only, hetero-only-except-ch)
    bond_marker_size_mm: float = 1.5 # Size of bond number marker [mm] (inherits default when None)
    bond_marker_depth_mm: float = 0.5 # Depth of bond number marker [mm] (inherits default when None)
@dataclass
class DefaultConfig:
    atom: DefaultAtomConfig = field(default_factory=lambda: DefaultAtomConfig())
    bond: DefaultBondConfig = field(default_factory=lambda: DefaultBondConfig())


@dataclass
class AtomConfig:
    name: str = field(default=None)  # Identical name of the atom (e.g., C_1, H_2)
    position: List[float] = field(
        default=None
    )  # Position of the atom [x, y, z] in Angstrom
    scale: float = None  # Scale factor for van der Waals radius
    color: List[int] = None  # Color of the atom in RGBA format [R, G, B, A]


@dataclass
class ShapeConfig:
    shape_type: str = None  # Type of the shape (e.g., spin, fixed, hole)
    shaft_radius: float = None  # Radius of the shaft [Angstrom]
    shaft_length: float = None  # Length of the shaft [Angstrom]
    shaft_radius_mm: float = None  # Radius of the shaft [mm]
    shaft_length_mm: float = None  # Length of the shaft [mm]
    stopper_radius: float = None  # Radius of the stopper [Angstrom]
    stopper_length: float = None  # Length of the stopper [Angstrom]
    hole_radius: float = None  # Radius of the hole [Angstrom]
    hole_length: float = None  # Length of the hole [Angstrom]
    hole_radius_mm: float = None  # Radius of the hole [mm]
    hole_length_mm: float = None  # Length of the hole [mm]
    chamfer_length: float = None  # Length of the chamfer [Angstrom]
    wall_thickness: float = None  # Thickness of the wall [Angstrom]
    shaft_gap: float = None  # Gap between the shaft and the cavity [Angstrom]
    shaft_gap_mm: float = None  # Gap between the shaft and the cavity [mm]
    taper_radius_scale: float = None  # Scale factor for the taper radius
    taper_angle_deg: float = None  # Taper angle in degrees
    bond_gap: float = None  # Gap between the bond plane [Angstrom]
    bond_gap_mm: float = None  # Gap between the bond plane [mm]

    taper_distance: float = None  # Distance for tapering [Angstrom]
    taper_height: float = None  # Height for tapering [Angstrom]


@dataclass
class BondConfig:
    atom_pair: List[str] = field(
        default_factory=lambda: ["None", "None"]
    )  # Names of the two atoms forming the bond
    shape_pair: List[ShapeConfig] = field(
        default_factory=lambda: [ShapeConfig(), ShapeConfig()]
    )
    bond_type: str = "none"  # Type of the bond (e.g., single, double, triple)
    shaft_types: List[str] = None  # Types of the two shafts (e.g., spin, fixed, hole)
    shape_type: List[str] = None  # Types of the two shapes (e.g., spin, fixed, hole)
    shaft_radius: float = None  # Radius of the shaft [Angstrom]
    shaft_length: float = None  # Length of the shaft [Angstrom]
    stopper_radius: float = None  # Radius of the stopper [Angstrom]
    stopper_length: float = None  # Length of the stopper [Angstrom]
    hole_radius: float = None  # Radius of the hole [Angstrom]
    hole_length: float = None  # Length of the hole [Angstrom]
    chamfer_length: float = None  # Length of the chamfer [Angstrom]
    wall_thickness: float = None  # Thickness of the wall [Angstrom]
    shaft_gap: float = None  # Gap between the shaft and the cavity [Angstrom]
    shaft_gap_mm: float = None  # Gap between the shaft and the cavity [mm]
    bond_gap: float = None  # Gap between the bond plane [Angstrom]
    bond_gap_mm: float = None  # Gap between the bond plane [mm]
    taper_angle_deg: List[float] = None  # Taper angle at the two ends in degrees
    taper_radius_scale: List[float] = None  # Scale factor for the taper radius at the two ends
    magnetic_hole_radius_mm: float = None  # Radius of the magnetic hole [mm]
    magnetic_hole_length_mm: float = None  # Length of the magnetic hole [mm]
    bond_marker: str = None # Option for bond number marker display (inherits default when None)
    bond_marker_size_mm: float = None # Size of bond number marker [mm] (inherits default when None)
    bond_marker_depth_mm: float = None # Depth of bond number marker [mm] (inherits default when None)

@dataclass
class MoleculeConfig:
    name: str  # Name of the molecule
    scale: float  # Scale factor for the whole model
    default: DefaultConfig
    atoms: List[AtomConfig]
    bonds: List[BondConfig]


def load_molfidget_config(file_path: str) -> MoleculeConfig:
    yaml = YAML()

    with open(file_path, "r") as file:
        data = yaml.load(file)

    molecule_config = from_dict(data_class=MoleculeConfig, data=data["molecule"])

    return molecule_config


def molecule_config_representer(dumper, data):
    data_dict = data.__dict__
    return dumper.represent_mapping("tag:yaml.org,2002:map", {"molecule": data_dict})


def default_config_representer(dumper, data):
    cmap = CommentedMap()
    cmap["atom"] = data.atom
    cmap["bond"] = data.bond
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def default_atom_config_representer(dumper, data):
    data_dict = OrderedDict(asdict(data))
    print("default_atom_config:", data_dict)
    cmap = CommentedMap()
    for key, value in data_dict.items():
        print(key, value)
        if value is None:
            continue
        if key in ("color",):
            cmap[key] = CommentedSeq(value)
            cmap[key].fa.set_flow_style()
        else:
            cmap[key] = value
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def default_bond_config_representer(dumper, data):
    data_dict = OrderedDict(asdict(data))
    cmap = CommentedMap()
    for key, value in data_dict.items():
        if value is None:
            continue
        else:
            cmap[key] = value
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def atom_config_representer(dumper, data):
    """AtomConfigをOrderedDictとして表現し、フィールド順序を保持"""
    field_dict = OrderedDict(asdict(data))
    cmap = CommentedMap()

    for key, value in field_dict.items():
        if value is None:
            continue
        if key == "position":
            cmap[key] = CommentedSeq(value)
            cmap[key].fa.set_flow_style()  # Set flow style for position
        else:
            cmap[key] = value
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)

def shape_config_representer(dumper, data):
    """ShapeConfigをOrderedDictとして表現し、フィールド順序を保持"""
    print('shapeconfigrepresenter:', data)
    field_dict = OrderedDict(asdict(data))
    cmap = CommentedMap()

    for key, value in field_dict.items():
        print(key, value)
        if value is None:
            continue
        else:
            cmap[key] = value

    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def bond_config_representer(dumper, data, default_config: DefaultBondConfig = None):
    """BondConfigをOrderedDictとして表現し、フィールド順序を保持"""
    # shape_pairを除外してasdict()を実行
    field_dict = OrderedDict()
    for key in data.__dataclass_fields__.keys():
        value = getattr(data, key)
        if key == "shape_pair":
            # shape_pairはオブジェクトのまま保持（辞書に変換しない）
            field_dict[key] = value
        else:
            field_dict[key] = value

    cmap = CommentedMap()

    for key, value in field_dict.items():
        if value is None:
            continue
        if key in ("atom_pair", "shaft_types", "taper_height", "taper_distance", "taper_angle_deg", "taper_radius_scale"):
            cmap[key] = CommentedSeq(value)
            cmap[key].fa.set_flow_style()
        else:
            cmap[key] = value

    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def save_molfidget_config(config: MoleculeConfig, file_path: str):
    yaml = YAML()
    yaml.representer.add_representer(MoleculeConfig, molecule_config_representer)
    yaml.representer.add_representer(DefaultConfig, default_config_representer)
    yaml.representer.add_representer(DefaultAtomConfig, default_atom_config_representer)
    yaml.representer.add_representer(DefaultBondConfig, default_bond_config_representer)
    yaml.representer.add_representer(AtomConfig, atom_config_representer)
    yaml.representer.add_representer(BondConfig, bond_config_representer)
    yaml.representer.add_representer(ShapeConfig, shape_config_representer)

    if not file_path:
        # 標準出力へ
        import sys
        yaml.dump(config, sys.stdout)
        return
    with open(file_path, "w") as file:
        yaml.dump(config, file)


def load_mol_file(file_name: str) -> MoleculeConfig:
    # Load a MOL file and populate the molecule with atoms and bonds
    with open(file_name, "r") as file:
        lines = file.readlines()
    # First line contains the name of the molecule
    name = lines[0].strip()
    atom_count = int(lines[3][0:3])
    bond_count = int(lines[3][3:6])
    # Parse the atom lines
    atoms = []
    for i in range(atom_count):
        data = lines[4 + i].strip().split()
        x = float(data[0])
        y = float(data[1])
        z = float(data[2])
        name = data[3] + f"_{i + 1}"
        atoms.append(AtomConfig(name=name, position=[x, y, z]))
    # Parse the bond lines
    bonds = []
    for i in range(bond_count):
        line = lines[4 + atom_count + i]
        id1 = int(line[0:3])
        id2 = int(line[3:6])
        type = int(line[6:9])
        if type == 1:
            bond_type = "single"
        elif type == 2:
            bond_type = "double"
        elif type == 3:
            bond_type = "triple"
        elif type == 4:
            bond_type = "aromatic"
        elif type == 9:
            bond_type = "magnetic"
        else:
            raise ValueError(f"Unknown bond type: {type}")
        bonds.append(
            BondConfig(
                atom_pair=[atoms[id1 - 1].name, atoms[id2 - 1].name],
                bond_type=bond_type,
            )
        )
    # file_nameから拡張子を除いた名前を設定
    molecle_config = MoleculeConfig(
        name=file_name.split("/")[-1].split(".")[0], scale=10.0, default=DefaultConfig(), atoms=atoms, bonds=bonds
    )
    return molecle_config


def load_pdb_file(file_name: str) -> MoleculeConfig:
    # Load a PDB file and populate the molecule with atoms and bonds
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb", file_name)

    # Determine molecule name from header or fallback to filename
    header = getattr(structure, "header", {}) or {}
    name = header.get("name") or None
    if not name:
        compound = header.get("compound")
        if isinstance(compound, dict):
            for key in ("1", "2", "3"):
                if key in compound and isinstance(compound[key], dict):
                    name = compound[key].get("molecule") or compound[key].get("MOLECULE")
                    if name:
                        break
            if not name and isinstance(compound.get("molecule"), str):
                name = compound.get("molecule")
        elif isinstance(compound, list):
            for entry in compound:
                if isinstance(entry, dict):
                    name = entry.get("molecule") or entry.get("MOLECULE")
                    if name:
                        break
    if not name:
        name = os.path.basename(file_name).split(".")[0]

    # Build atom list and mappings
    atoms_by_serial = {}
    atom_config_by_atom = {}
    element_by_atom = {}
    elem_counts = {}

    for atom in structure.get_atoms():
        serial = atom.get_serial_number()
        if serial is None:
            continue
        residue = atom.get_parent()
        if residue is not None: # 水分子をスキップ
            resname = residue.get_resname()
            if isinstance(resname, str) and resname.strip().upper() == "HOH":
                continue

        elem = getattr(atom, "element", "") or ""
        elem = str(elem).strip()
        if elem.lower().endswith("new"): #　reduceでつけられたHの接尾辞を削除
            elem = elem[:-3].strip()

        if not elem:
            # Fallback to atom name if element is missing
            atom_name = atom.get_name().strip()
            if atom_name.lower().endswith("new"):
                atom_name = atom_name[:-3].strip()
            atom_name = atom_name.lstrip("0123456789")
            elem = atom_name[:2].strip().capitalize()
            if len(elem) == 2 and elem[1].isalpha():
                elem = elem[0].upper() + elem[1].lower()
            else:
                elem = elem[0].upper()

        # Name atoms as Element_#
        elem_counts[elem] = elem_counts.get(elem, 0) + 1
        x, y, z = atom.get_coord().tolist()
        atom_config = AtomConfig(name=f"{elem}_{elem_counts[elem]}", position=[x, y, z])

        # Map by atom object for distance-based bonds
        atom_config_by_atom[atom] = atom_config
        element_by_atom[atom] = elem

        # Map by serial if possible (for CONECT)
        if int(serial) not in atoms_by_serial:
            atoms_by_serial[int(serial)] = atom_config

    atoms = [atom_config_by_atom[a] for a in atom_config_by_atom.keys()]

    # Build bonds by distance (fallback to CONECT when present)
    bonds = []
    seen = set()
    # 1) Try distance-based bonding using Bio.PDB NeighborSearch
    atoms_list = list(atom_config_by_atom.keys())
    if atoms_list:
        ns = NeighborSearch(atoms_list)
        # Global max search radius (covers typical covalent bond lengths)
        max_radius = 2.2
        for atom in atoms_list:
            id1 = atom.get_serial_number()
            if id1 is None:
                continue
            elem1 = element_by_atom.get(atom)
            if not elem1:
                continue
            for nb in ns.search(atom.get_coord(), max_radius, level="A"):
                id2 = nb.get_serial_number()
                if id2 is None:
                    continue
                if id1 == id2:
                    continue
                pair = tuple(sorted((id(atom), id(nb))))
                if pair in seen:
                    continue
                elem2 = element_by_atom.get(nb)
                if not elem2:
                    continue
                dist = atom - nb
                # Determine cutoff and bond type
                key = frozenset([elem1, elem2])
                bond_type = "single"
                if key in bond_distance_table:
                    single_d, double_d, triple_d = bond_distance_table[key]
                    cutoff = 1.1 * single_d
                    if triple_d > 0 and dist <= 1.05 * triple_d:
                        bond_type = "triple"
                    elif double_d > 0 and dist <= 1.05 * double_d:
                        bond_type = "double"
                    else:
                        bond_type = "single"
                elif "H" in (elem1, elem2):
                    cutoff = 1.2
                else:
                    cutoff = 1.9
                if dist <= cutoff:
                    a1 = atom_config_by_atom.get(atom)
                    a2 = atom_config_by_atom.get(nb)
                    if a1 is None or a2 is None:
                        continue
                    seen.add(pair)
                    bonds.append(
                        BondConfig(
                            atom_pair=[a1.name, a2.name],
                            bond_type=bond_type,
                        )
                    )

    # 2) Parse CONECT records and add any missing bonds
    with open(file_name, "r") as file:
        for line in file:
            if not line.startswith("CONECT"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                id1 = int(parts[1])
            except ValueError:
                continue
            for id2_str in parts[2:]:
                try:
                    id2 = int(id2_str)
                except ValueError:
                    continue
                if id1 == id2:
                    continue
                pair = (min(id1, id2), max(id1, id2))
                if pair in seen:
                    continue
                a1 = atoms_by_serial.get(pair[0])
                a2 = atoms_by_serial.get(pair[1])
                if a1 is None or a2 is None:
                    continue
                seen.add(pair)
                bonds.append(
                    BondConfig(
                        atom_pair=[a1.name, a2.name],
                        bond_type="single",
                    )
                )

    # Detect hydrogen bonds and add as magnetic bonds
    name_to_config = {cfg.name: cfg for cfg in atom_config_by_atom.values()}
    elem_by_name = {cfg.name: element_by_atom[atom] for atom, cfg in atom_config_by_atom.items()}
    adjacency = {name: set() for name in name_to_config.keys()}
    for bond in bonds:
        if bond.bond_type == "magnetic":
            continue
        a1, a2 = bond.atom_pair
        if a1 in adjacency and a2 in adjacency:
            adjacency[a1].add(a2)
            adjacency[a2].add(a1)

    hydrogen_names = [n for n, e in elem_by_name.items() if e == "H"]
    acceptor_names = [n for n, e in elem_by_name.items() if e in ("O", "N", "S")]
    magnetic_seen = set()

    for h_name in hydrogen_names:
        h_cfg = name_to_config.get(h_name)
        if h_cfg is None:
            continue
        h_pos = np.array(h_cfg.position, dtype=float)
        donors = [n for n in adjacency.get(h_name, []) if elem_by_name.get(n) in ("O", "N")]
        if not donors:
            continue
        for d_name in donors:
            d_cfg = name_to_config.get(d_name)
            if d_cfg is None:
                continue
            d_pos = np.array(d_cfg.position, dtype=float)
            dh_vec = d_pos - h_pos
            dh_norm = np.linalg.norm(dh_vec)
            if dh_norm == 0.0:
                continue
            for a_name in acceptor_names:
                if a_name == h_name or a_name == d_name:
                    continue
                if a_name in adjacency.get(h_name, set()):
                    continue
                a_cfg = name_to_config.get(a_name)
                if a_cfg is None:
                    continue
                a_pos = np.array(a_cfg.position, dtype=float)
                ha_vec = a_pos - h_pos
                da_vec = a_pos - d_pos
                ha_dist = np.linalg.norm(ha_vec)
                da_dist = np.linalg.norm(da_vec)
                if ha_dist > hydrogen_bond_params["HA_max"] or da_dist > hydrogen_bond_params["DA_max"]:
                    continue
                ha_norm = ha_dist
                if ha_norm == 0.0:
                    continue
                cos_angle = float(np.dot(dh_vec, ha_vec) / (dh_norm * ha_norm))
                cos_angle = max(-1.0, min(1.0, cos_angle))
                angle_deg = np.degrees(np.arccos(cos_angle))
                if angle_deg < hydrogen_bond_params["DHA_min_deg"]:
                    continue
                pair = frozenset([h_name, a_name])
                if pair in magnetic_seen:
                    continue
                magnetic_seen.add(pair)
                bonds.append(
                    BondConfig(
                        atom_pair=[h_name, a_name],
                        bond_type="magnetic",
                    )
                )

    molecle_config = MoleculeConfig(
        name=name, scale=8.0, default=DefaultConfig(), atoms=atoms, bonds=bonds
    )
    return molecle_config
