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
    vdw_scale: float = 0.8  # Scale factor for van der Waals radius
    color: List[int] = field(default_factory=lambda: [200, 200, 200, 255])


@dataclass
class BondTypeShapeConfig:
    shaft_radius: float = None
    shaft_length: float = None
    shaft_radius_mm: float = None
    shaft_length_mm: float = None
    stopper_radius: float = None
    stopper_length: float = None
    hole_radius: float = None
    hole_length: float = None
    hole_radius_mm: float = None
    hole_length_mm: float = None
    chamfer_length: float = None
    wall_thickness: float = None
    shaft_gap_mm: float = None
    bond_gap_mm: float = None


@dataclass
class DefaultBondConfig:
    # 共通パラメータ
    bond_gap_mm: float = 0.1
    bond_marker: str = "hetero-only-except-ch"
    bond_marker_size_mm: float = 1.5
    bond_marker_depth_mm: float = 0.5
    taper_angle_deg: float = 0.0
    taper_radius_scale: float = 1.0
    # bond_typeごとのデフォルト
    spin: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        shaft_radius=0.3, shaft_length=0.3, hole_radius=0.3, hole_length=0.3,
        chamfer_length=0.1, wall_thickness=0.1, shaft_gap_mm=0.3,
        stopper_radius=0.4, stopper_length=0.2))
    normal: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        shaft_radius=0.3, shaft_length=0.3, hole_radius=0.3, hole_length=0.3,
        chamfer_length=0.1))
    fixed: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        shaft_radius=0.3, shaft_length=0.3, hole_radius=0.3, hole_length=0.3,
        chamfer_length=0.1))
    gapped: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        shaft_radius=0.8, shaft_length=0.3, hole_radius=0.8, hole_length=0.3,
        chamfer_length=0.1))
    short: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        shaft_radius=0.3, shaft_length=0.2, hole_radius=0.3, hole_length=0.3,
        chamfer_length=0.1))
    holes: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        hole_radius_mm=3.525, hole_length_mm=2.0, bond_gap_mm=0.0))
    notch_2: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        shaft_radius=0.3, shaft_length=0.3, hole_radius=0.3, hole_length=0.3,
        chamfer_length=0.1, wall_thickness=0.1, shaft_gap_mm=0.3,
        stopper_radius=0.4, stopper_length=0.2))
    notch_3: BondTypeShapeConfig = field(default_factory=lambda: BondTypeShapeConfig(
        shaft_radius=0.3, shaft_length=0.3, hole_radius=0.3, hole_length=0.3,
        chamfer_length=0.1, wall_thickness=0.1, shaft_gap_mm=0.3,
        stopper_radius=0.4, stopper_length=0.2))


@dataclass
class DefaultConfig:
    atom: DefaultAtomConfig = field(default_factory=DefaultAtomConfig)
    bond: DefaultBondConfig = field(default_factory=DefaultBondConfig)


@dataclass
class AtomConfig:
    name: str = field(default=None)  # Identical name of the atom (e.g., C_1, H_2)
    position: List[float] = field(
        default=None
    )  # Position of the atom [x, y, z] in Angstrom
    vdw_scale: float = None  # Scale factor for van der Waals radius
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
    shaft_gap_mm: float = None  # Gap between the shaft and the cavity [mm]
    taper_radius_scale: float = None  # Scale factor for the taper radius
    taper_angle_deg: float = None  # Taper angle in degrees
    bond_gap_mm: float = None  # Gap between the bond plane [mm]

    taper_distance: float = None  # Distance for tapering [Angstrom]
    taper_height: float = None  # Height for tapering [Angstrom]


@dataclass
class BondConfig:
    atom_pair: List[str] = field(
        default_factory=lambda: ["None", "None"]
    )
    shape_pair: List[ShapeConfig] = field(
        default_factory=lambda: [ShapeConfig(), ShapeConfig()]
    )
    bond_type: str = "none"
    bond_gap_mm: float = None
    bond_marker: str = None
    bond_marker_size_mm: float = None
    bond_marker_depth_mm: float = None


@dataclass
class MoleculeConfig:
    name: str
    scale: float
    atoms: List[AtomConfig]
    bonds: List[BondConfig]


@dataclass
class MolfidgetConfig:
    default: DefaultConfig
    molecule: MoleculeConfig


OLD_BOND_TYPES = {
    "single": "spin", "double": "fixed", "triple": "fixed",
    "aromatic": "gapped", "magnetic": "holes", "1.5": "gapped",
}


def load_molfidget_config(file_path: str) -> MolfidgetConfig:
    yaml = YAML()

    with open(file_path, "r") as file:
        data = yaml.load(file)

    # 旧形式の検出
    if "shape" in data.get("default", {}):
        raise ValueError(
            "旧形式のMLFファイルです。default.shape は廃止されました。"
            "default.bond.{bond_type} サブセクションに移行してください。"
        )
    for bond_data in data.get("molecule", {}).get("bonds", []):
        bt = str(bond_data.get("bond_type", ""))
        if bt in OLD_BOND_TYPES:
            raise ValueError(
                f"bond_type '{bt}' は廃止されました。"
                f"代わりに '{OLD_BOND_TYPES[bt]}' を使用してください。"
            )

    deprecated_shape_keys = {
        "shaft_gap": "shaft_gap_mm",
        "bond_gap": "bond_gap_mm",
    }

    def check_deprecated_shape_keys(node, path="root"):
        if isinstance(node, dict):
            for key, replacement in deprecated_shape_keys.items():
                if key in node:
                    raise ValueError(
                        f"{path}.{key} は廃止されました。"
                        f" 代わりに {replacement} を使用してください。"
                    )
            for child_key, child_value in node.items():
                check_deprecated_shape_keys(child_value, f"{path}.{child_key}")
        elif isinstance(node, list):
            for index, item in enumerate(node):
                check_deprecated_shape_keys(item, f"{path}[{index}]")

    check_deprecated_shape_keys(data)

    default_config = from_dict(
        data_class=DefaultConfig, data=data.get("default", {})
    )
    molecule_config = from_dict(
        data_class=MoleculeConfig, data=data["molecule"]
    )

    return MolfidgetConfig(default=default_config, molecule=molecule_config)


def molfidget_config_representer(dumper, data):
    cmap = CommentedMap()
    cmap["default"] = data.default
    cmap["molecule"] = data.molecule
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def molecule_config_representer(dumper, data):
    cmap = CommentedMap()
    cmap["name"] = data.name
    cmap["scale"] = data.scale
    cmap["atoms"] = data.atoms
    cmap["bonds"] = data.bonds
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def default_config_representer(dumper, data):
    cmap = CommentedMap()
    cmap["atom"] = data.atom
    cmap["bond"] = data.bond
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def default_atom_config_representer(dumper, data):
    data_dict = OrderedDict(asdict(data))
    cmap = CommentedMap()
    for key, value in data_dict.items():
        if value is None:
            continue
        if key in ("color",):
            cmap[key] = CommentedSeq(value)
            cmap[key].fa.set_flow_style()
        else:
            cmap[key] = value
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


BOND_TYPE_FIELDS = ("spin", "normal", "fixed", "gapped", "short", "holes", "notch_2", "notch_3")


def bond_type_shape_config_representer(dumper, data):
    data_dict = OrderedDict(asdict(data))
    cmap = CommentedMap()
    for key, value in data_dict.items():
        if value is None:
            continue
        cmap[key] = value
    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def default_bond_config_representer(dumper, data):
    data_dict = OrderedDict(asdict(data))
    cmap = CommentedMap()
    # 共通パラメータを先に出力
    for key, value in data_dict.items():
        if key in BOND_TYPE_FIELDS:
            continue
        if value is None:
            continue
        cmap[key] = value
    # bond_typeサブセクションを出力
    for bt_name in BOND_TYPE_FIELDS:
        bt_config = getattr(data, bt_name, None)
        if bt_config is not None:
            cmap[bt_name] = bt_config
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
    field_dict = OrderedDict(asdict(data))
    cmap = CommentedMap()

    for key, value in field_dict.items():
        if value is None:
            continue
        else:
            cmap[key] = value

    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def bond_config_representer(dumper, data):
    """BondConfigをOrderedDictとして表現し、フィールド順序を保持"""
    field_dict = OrderedDict()
    for key in data.__dataclass_fields__.keys():
        value = getattr(data, key)
        if key == "shape_pair":
            field_dict[key] = value
        else:
            field_dict[key] = value

    cmap = CommentedMap()
    for key, value in field_dict.items():
        if value is None:
            continue
        if key == "atom_pair":
            cmap[key] = CommentedSeq(value)
            cmap[key].fa.set_flow_style()
        else:
            cmap[key] = value

    return dumper.represent_mapping("tag:yaml.org,2002:map", cmap)


def save_molfidget_config(config: MolfidgetConfig, file_path: str):
    yaml = YAML()
    yaml.representer.add_representer(MolfidgetConfig, molfidget_config_representer)
    yaml.representer.add_representer(MoleculeConfig, molecule_config_representer)
    yaml.representer.add_representer(DefaultConfig, default_config_representer)
    yaml.representer.add_representer(DefaultAtomConfig, default_atom_config_representer)
    yaml.representer.add_representer(DefaultBondConfig, default_bond_config_representer)
    yaml.representer.add_representer(BondTypeShapeConfig, bond_type_shape_config_representer)
    yaml.representer.add_representer(AtomConfig, atom_config_representer)
    yaml.representer.add_representer(BondConfig, bond_config_representer)
    yaml.representer.add_representer(ShapeConfig, shape_config_representer)

    if not file_path:
        import sys
        yaml.dump(config, sys.stdout)
        return
    with open(file_path, "w") as file:
        yaml.dump(config, file)


def load_mol_file(file_name: str) -> MolfidgetConfig:
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
            bond_type = "spin"
        elif type == 2:
            bond_type = "fixed"
        elif type == 3:
            bond_type = "fixed"
        elif type == 4:
            bond_type = "gapped"
        elif type == 9:
            bond_type = "holes"
        else:
            raise ValueError(f"Unknown bond type: {type}")
        bonds.append(
            BondConfig(
                atom_pair=[atoms[id1 - 1].name, atoms[id2 - 1].name],
                bond_type=bond_type,
            )
        )
    molecle_config = MoleculeConfig(
        name=file_name.split("/")[-1].split(".")[0],
        scale=10.0, atoms=atoms, bonds=bonds,
    )
    return MolfidgetConfig(default=DefaultConfig(), molecule=molecle_config)


def load_pdb_file(file_name: str) -> MolfidgetConfig:
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
                bond_type = "spin"
                if key in bond_distance_table:
                    single_d, double_d, triple_d = bond_distance_table[key]
                    cutoff = 1.1 * single_d
                    if triple_d > 0 and dist <= 1.05 * triple_d:
                        bond_type = "fixed"
                    elif double_d > 0 and dist <= 1.05 * double_d:
                        bond_type = "fixed"
                    else:
                        bond_type = "spin"
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
                        bond_type="spin",
                    )
                )

    # Detect hydrogen bonds and add as magnetic bonds
    name_to_config = {cfg.name: cfg for cfg in atom_config_by_atom.values()}
    elem_by_name = {cfg.name: element_by_atom[atom] for atom, cfg in atom_config_by_atom.items()}
    adjacency = {name: set() for name in name_to_config.keys()}
    for bond in bonds:
        if bond.bond_type == "holes":
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
                        bond_type="holes",
                    )
                )

    molecle_config = MoleculeConfig(
        name=name, scale=8.0, atoms=atoms, bonds=bonds,
    )
    return MolfidgetConfig(default=DefaultConfig(), molecule=molecle_config)
