# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要

**molfidget** は、分子構造データ（MOL/PDB形式）を3Dプリンタで造形可能なSTL/3MFファイルに変換するPythonツール。単結合を回転軸として可動させられる物理的な分子模型を生成する。Python 3.12以上が必要。

## 開発環境

```bash
# 必要なシステムパッケージ
sudo apt install git libglut-dev python3-poetry

# 依存関係のインストール
poetry install

# ビルド
poetry build
```

## 主要コマンド

```bash
# PDB/MOLファイルをMLF形式に変換
poetry run molfidget convert <input.pdb|input.mol> [output.mlf]

# MLFファイルを3Dプレビュー（'q'で終了、'w'でワイヤフレーム）
poetry run molfidget preview <input.mlf>

# MLFファイルからSTL/3MFファイルを生成
poetry run molfidget generate --scale 10.0 --output-dir output <input.mlf>
```

テストフレームワークは未導入。動作確認は `data/mlf/` 配下のサンプルファイルでpreview/generateを実行して行う。

## アーキテクチャ

### データフロー
```
PDB/MOL File → (config.py) → MoleculeConfig → (save) → MLF File (YAML)
MLF File → (load) → MoleculeConfig → Molecule → create_trimesh_scene() → STL/3MF
```

### 主要モジュールと役割

- **molfidget.py**: CLIエントリーポイント。3つのサブコマンド（convert/preview/generate）と3MFエクスポート機能
- **config.py**: データクラス群（`MoleculeConfig`, `AtomConfig`, `BondConfig`, `ShapeConfig`）、PDB/MOLパーサー、MLF（YAML）の読み書き。`dacite.from_dict`でYAML→dataclass変換、`ruamel.yaml`でフロースタイル保持した書き出し
- **molecule.py**: `Molecule`クラス。Atom/Bondの初期化、trimeshシーン生成、STL/3MF出力、元素グループへのマージ
- **bond.py**: `Bond`クラス。bond_typeに応じたshape_pairの自動設定、結合平面でのスライス、ノッチ機構（notch_1/2/3）、ボンドマーカー刻印
- **shape.py**: `Shape`クラス。軸（shaft/shaft_spin/shaft_dcut）と穴（hole/hole_dcut）の3Dジオメトリ生成。テーパー・面取り・Dカット処理
- **atom.py**: `Atom`クラス。ファン・デル・ワールス半径の球体メッシュ生成。bondによって形状が切り出される
- **constants.py**: 元素ごとのファン・デル・ワールス半径、CPKカラー、結合距離テーブル
- **labeled_scene_viewer.py**: pygletベースの3Dプレビューアー

### 重要な設計パターン

**パラメータ解決の優先順位**: shape_pair個別値 > BondConfig値 > DefaultBondConfig値。mm単位のパラメータ（`*_mm`）が指定されている場合、Å単位より常に優先される。

**形状生成の流れ**: `Molecule.create_trimesh_scene()` → 各Atomの球体メッシュ生成 → 各Bondの `sculpt_atoms2()` で結合平面スライス＋ジョイント形状の切り出し（trimeshブーリアン演算）→ シーンに追加

**bond_typeとデフォルト形状の対応**:
| bond_type | atom1側 | atom2側 |
|-----------|---------|---------|
| single, notch_N | shaft_spin | hole |
| double, triple | shaft_dcut | hole_dcut |
| aromatic | shaft | hole |
| magnetic | hole（磁石穴） | hole（磁石穴） |
| plane | none | none |

## 単位系

- **内部座標**: Angstrom（Å）
- **ギャップ・クリアランス系パラメータ**: mm（`bond_gap_mm`, `shaft_radius_mm`等）
- **STL出力**: scaleパラメータでÅ→mm変換（`scale=10.0`推奨）
- mm→Å変換は `value_mm / scale` で行われる（shape.py, bond.py参照）

## MLFファイル形式

YAML形式の中間表現（拡張子: `.mlf`）。詳細仕様は [docs/mlf_file.md](docs/mlf_file.md) を参照。

```yaml
molecule:
  name: 分子名
  scale: 10.0  # [mm/Å]
  default:
    atom:
      vdw_scale: 0.8
    bond:
      bond_gap_mm: 0.1
      shaft_radius: 0.3
      # ...
  atoms:
    - name: C_1  # 元素名_連番
      position: [x, y, z]  # Å単位
  bonds:
    - atom_pair: [C_1, H_1]
      bond_type: single
      shape_pair:  # オプション：個別パラメータ上書き
        - shape_type: shaft_spin
        - shape_type: hole
```

## サンプルデータ

テスト用の小分子: `data/mlf/hydrogen.mlf`, `data/mlf/h2o.mlf`, `data/mlf/ethanol.mlf`
複雑な分子: `data/mlf/cyclohexane.mlf`, `data/pdb/ibuprofen.pdb`

## 配布

PyPIパッケージ名: `molfidget`（バージョン: 1.0.0）
GitHub Actionsで自動公開: `.github/workflows/pypi-publish.yml`
