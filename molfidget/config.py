from collections import OrderedDict
from dataclasses import asdict, dataclass, field

from dacite import from_dict
from ruamel.yaml import YAML, CommentedMap, CommentedSeq
from typing import List


@dataclass
class DefaultAtomConfig:
    scale: float = 1.0  # Scale factor for van der Waals radius
    color: List[int] = field(default_factory=lambda: [200, 200, 200, 255])


@dataclass
class DefaultBondConfig:
    bond_gap_mm: float = 0.1  # Gap between the bond plane [mm]
    bond_gap: float = 0.0  # Gap between the bond plane [Angstrom]
    shaft_radius: float = 0.3  # Radius of the shaft [Angstrom]
    shaft_length: float = 0.3  # Length of the shaft [Angstrom]
    stopper_radius: float = 0.4  # Radius of the stopper [Angstrom]
    stopper_length: float = 0.2  # Length of the stopper [Angstrom]
    hole_radius: float = 0.3  # Radius of the hole [Angstrom]
    hole_length: float = 0.3  # Length of the hole [Angstrom]
    chamfer_length: float = 0.1  # Length of the chamfer [Angstrom]
    wall_thickness: float = 0.1  # Thickness of the wall [Angstrom]
    shaft_gap: float = 0.03  # Gap between the shaft and the cavity [Angstrom]
    taper_angle_deg: float = 0.0 # Taper angle in degrees
    taper_radius_scale: float = 1.0 # Scale factor for the taper radius


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
    stopper_radius: float = None  # Radius of the stopper [Angstrom]
    stopper_length: float = None  # Length of the stopper [Angstrom]
    hole_radius: float = None  # Radius of the hole [Angstrom]
    hole_length: float = None  # Length of the hole [Angstrom]
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
            bond_type = "1.5"
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
        name=file_name.split("/")[-1].split(".")[0], scale=1.0, default=DefaultConfig(), atoms=atoms, bonds=bonds
    )
    return molecle_config


def load_pdb_file(file_name: str) -> MoleculeConfig:
    # Load a PDB file and populate the molecule with atoms and bonds
    with open(file_name, "r") as file:
        lines = file.readlines()

    atoms = []
    bonds = []
    for line in lines:
        # Do we have the name of the molecule somewhere?
        if line.startswith("COMPND"):
            name = line[10:].strip()
        # Parse ATOM and HETATM lines to extract atom information
        if line.startswith("ATOM") or line.startswith("HETATM"):
            id = int(line[6:11].strip())
            name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atoms.append(AtomConfig(name=name + f"_{id}", position=[x, y, z]))
        if line.startswith("CONECT"):
            parts = line.split()
            id1 = int(parts[1])
            for id2_str in parts[2:]:
                id2 = int(id2_str)
                # Avoid duplicate bonds
                if id1 < id2:
                    bonds.append(
                        BondConfig(
                            atom_pair=[atoms[id1 - 1].name, atoms[id2 - 1].name],
                            bond_type="single",
                        )
                    )
    molecle_config = MoleculeConfig(
        name="file_name", scale=1.0, default=DefaultConfig(), atoms=atoms, bonds=bonds
    )
    return molecle_config
