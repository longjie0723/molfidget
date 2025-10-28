import argparse
import pyglet
from molfidget.config import (
    load_pdb_file,
    load_mol_file,
    save_molfidget_config,
    load_molfidget_config,
)
from molfidget.labeled_scene_viewer import LabeledSceneViewer
from molfidget.molecule import Molecule


def setup_argparse():
    parser = argparse.ArgumentParser(description="Molecule visualization and manipulation")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    convert_parser = subparsers.add_parser("convert", help="Convert molecule file formats to molfidget format")
    convert_parser.add_argument("input_file", type=str, help="Input PDB or MOL file to convert")
    convert_parser.add_argument("output_file", type=str, help="Output molfidget YAML file")

    preview_parser = subparsers.add_parser("preview", help="Preview molecule from molfidget file")
    preview_parser.add_argument("molfidget_file", type=str, help="Input molfidget YAML file to preview")

    generate_parser = subparsers.add_parser("generate", help="Generate STL files from molfidget file")
    generate_parser.add_argument("molfidget_file", type=str, help="Input molfidget YAML file to generate STL files from")

    return parser


def exec_convert(args):
    # ファイル名が .pdb か .mol かで処理を分ける
    if args.input_file.endswith('.pdb'):
        molecule_config = load_pdb_file(args.input_file)
    elif args.input_file.endswith('.mol'):
        molecule_config = load_mol_file(args.input_file)
    else:
        raise ValueError("Input file must be a PDB or MOL file")
    save_molfidget_config(molecule_config, args.output_file)


def exec_preview(args):
    molecule_config = load_molfidget_config(args.molfidget_file)
    molecule = Molecule(molecule_config)
    print(f"Molecule: {molecule.name}")

    # Create trimesh model
    scene = molecule.create_trimesh_scene()

    # Show the model
    viewer = LabeledSceneViewer(scene)
    pyglet.app.run()


def exec_generate(args):
    config = load_molfidget_config(args.molfidget_file)
    molecule = Molecule(config)
    molecule.save_stl_files()


def main():
    parser = setup_argparse()
    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        return
    if args.command == "convert":
        exec_convert(args)
    elif args.command == "preview":
        exec_preview(args)
    elif args.command == "generate":
        exec_generate(args)


if __name__ == "__main__":
    main()