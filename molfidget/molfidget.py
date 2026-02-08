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
import os

try:
    import argcomplete
    ARGCOMPLETE_AVAILABLE = True
except ImportError:
    ARGCOMPLETE_AVAILABLE = False

import numpy as np

def _normalize_rgba_u8(rgba):
    rgba = np.asarray(rgba).reshape(-1)
    if rgba.size == 3:
        rgba = np.concatenate([rgba, [255]])
    rgba = rgba[:4]
    if rgba.dtype != np.uint8:
        if np.issubdtype(rgba.dtype, np.floating) and rgba.max() <= 1.0:
            rgba = (rgba * 255.0).round()
        rgba = rgba.astype(np.uint8)
    return rgba

def export_scene_as_colored_3mf(scene, out_path, libpath=None, debug=False):
    # import
    try:
        import lib3mf as Lib3MF
        wrapper = Lib3MF.Wrapper()
    except Exception:
        import Lib3MF
        if libpath is None:
            raise RuntimeError("Lib3MF SDK版を使うなら libpath（lib3mf共有ライブラリの場所）が必要です")
        wrapper = Lib3MF.Wrapper(libpath)

    def rgba_u8_to_lib3mf_color(rgba_u8):
        c = Lib3MF.Color()
        c.Red   = int(rgba_u8[0])
        c.Green = int(rgba_u8[1])
        c.Blue  = int(rgba_u8[2])
        c.Alpha = int(rgba_u8[3])
        return c

    def make_triangle_props(resource_id, property_id):
        tp = Lib3MF.TriangleProperties()
        try:
            for k in range(3):
                tp.ResourceID[k] = int(resource_id)
                tp.PropertyID[k] = int(property_id)
        except Exception:
            tp.ResourceID1 = int(resource_id); tp.PropertyID1 = int(property_id)
            tp.ResourceID2 = int(resource_id); tp.PropertyID2 = int(property_id)
            tp.ResourceID3 = int(resource_id); tp.PropertyID3 = int(property_id)
        return tp

    model = wrapper.CreateModel()

    # mm単位
    try:
        model.SetUnit(Lib3MF.ModelUnit.Millimeter)
    except Exception:
        pass

    # ColorGroup（色定義）
    color_group = model.AddColorGroup()
    color_group_uid = color_group.GetUniqueResourceID()

    # RGBA -> ColorGroup内のPropertyID
    color2pid = {}

    if debug:
        print("ColorGroup UID:", color_group_uid)

    for geom_name, tri in scene.geometry.items():
        tri2 = tri.copy()

        # Scene graph transform を焼く（サイズ/配置対策）
        try:
            mat, _ = scene.graph.get(geom_name)
            tri2.apply_transform(mat)
        except Exception as e:
            if debug:
                print("scene.graph.get failed:", geom_name, e)

        vertices = np.asarray(tri2.vertices, dtype=float)
        faces    = np.asarray(tri2.faces, dtype=np.int64)

        mesh_obj = model.AddMeshObject()

        # vertices
        for v in vertices:
            p = Lib3MF.Position()
            p.Coordinates[0] = float(v[0])
            p.Coordinates[1] = float(v[1])
            p.Coordinates[2] = float(v[2])
            mesh_obj.AddVertex(p)

        # triangles
        for f in faces:
            t = Lib3MF.Triangle()
            t.Indices[0] = int(f[0])
            t.Indices[1] = int(f[1])
            t.Indices[2] = int(f[2])
            mesh_obj.AddTriangle(t)

        # 色抽出（単色）
        vc = getattr(tri.visual, "vertex_colors", None)
        fc = getattr(tri.visual, "face_colors", None)
        if fc is not None and len(fc) > 0:
            rgba = _normalize_rgba_u8(fc[0])
        elif vc is not None and len(vc) > 0:
            rgba = _normalize_rgba_u8(vc[0])
        else:
            rgba = np.array([200, 200, 200, 255], dtype=np.uint8)

        key = tuple(int(x) for x in rgba[:4])
        if key not in color2pid:
            pid = color_group.AddColor(rgba_u8_to_lib3mf_color(rgba))
            color2pid[key] = pid

        pid_main = color2pid[key]

        # 1) triangle properties（試す：writerが出してくれればラッキー）
        tp_main = make_triangle_props(color_group_uid, pid_main)
        try:
            mesh_obj.SetAllTriangleProperties(tp_main)
        except TypeError:
            n = mesh_obj.GetTriangleCount()
            for i in range(n):
                mesh_obj.SetTriangleProperties(i, tp_main)

        # 2) ★重要：object-level も必ず付ける（pid/pindex を確実に出す）
        # Bambuが triangle-level を読まない/出力されない場合の保険
        mesh_obj.SetObjectLevelProperty(color_group_uid, pid_main)

        model.AddBuildItem(mesh_obj, Lib3MF.Transform())

    writer = model.QueryWriter("3mf")
    writer.WriteToFile(out_path)

def setup_argparse():
    parser = argparse.ArgumentParser(description="Molecule visualization and manipulation")
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    convert_parser = subparsers.add_parser("convert", help="Convert molecule file formats to molfidget format")
    convert_parser.add_argument("input_file", type=str, help="Input PDB or MOL file to convert")
    convert_parser.add_argument("output_file", type=str, nargs="?", default=None, help="Output molfidget YAML file (optional, defaults to stdout)")

    preview_parser = subparsers.add_parser("preview", help="Preview molecule from molfidget file")
    preview_parser.add_argument("molfidget_file", type=str, help="Input molfidget YAML file to preview")

    generate_parser = subparsers.add_parser("generate", help="Generate STL files from molfidget file")
    # add parameters for generation
    # scale parameter
    generate_parser.add_argument("--scale", type=float, help="Scale factor for the output model")
    generate_parser.add_argument("--output-dir", type=str, default="output", help="Output directory for STL files")
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
    molfidget_config = load_molfidget_config(args.molfidget_file)
    # Determine scale: command line > molecule config > default config
    if args.scale is not None:
        scale = args.scale
    elif molfidget_config.molecule.scale is not None:
        scale = molfidget_config.molecule.scale
    else:
        scale = molfidget_config.default.molecule.scale
    molfidget_config.molecule.scale = scale
    print(f"Using scale factor: {scale}")

    molecule = Molecule(molfidget_config)
    print(f"Molecule: {molecule.name}")
    # Create trimesh model
    scene = molecule.create_trimesh_scene()
    # Apply scale
    scene.apply_scale(scale)

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    # export the entire molecule as 3MF
    # scene.export(os.path.join(output_dir, f"{molecule.name}.3mf"))
    export_scene_as_colored_3mf(scene, os.path.join(output_dir, f"{molecule.name}.3mf"), libpath=None, debug=False)

    # Save STL files for each component
    molecule.save_stl_files(scale=scale, output_dir=output_dir)
    # Merge atoms into groups and save group STL files
    molecule.merge_atoms()
    molecule.save_group_stl_files(scale, output_dir=output_dir)

def main():
    parser = setup_argparse()

    # Enable argcomplete if available
    if ARGCOMPLETE_AVAILABLE:
        argcomplete.autocomplete(parser)

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