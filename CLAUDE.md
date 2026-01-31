# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要

**molfidget** は、分子構造データ（MOL形式・PDB形式）を3Dプリンタで造形可能なSTLファイルに変換するツールです。単結合を回転軸として可動させることができる物理的な分子模型を作成できます。

## 開発環境のセットアップ

### 必要なツール
```bash
sudo apt install git libglut-dev python3-poetry
```

### 依存関係のインストール
```bash
poetry install
```

## 主要コマンド

### CLIサブコマンド
```bash
# PDB/MOLファイルをMLF形式に変換
poetry run molfidget convert <input.pdb|input.mol> [output.mlf]

# MLFファイルを3Dプレビュー（'q'キーで終了）
poetry run molfidget preview <input.mlf>

# MLFファイルからSTLファイルを生成
poetry run molfidget generate --scale 10.0 --output-dir output <input.mlf>
```

### パッケージビルド
```bash
# ビルド（dist/にwheel/tar.gzが生成される）
poetry build

# バージョン更新
poetry version patch|minor|major
```

## アーキテクチャ

### データフロー
```
PDB/MOL File
  ↓ (config.py: load_pdb_file/load_mol_file)
MoleculeConfig (内部表現)
  ↓ (save_molfidget_config)
MLF File (YAML)
  ↓ (load_molfidget_config)
Molecule インスタンス
  ↓ (create_trimesh_scene)
Trimesh Scene (3Dジオメトリ)
  ↓
STL/3MF出力
```

### 主要モジュール

#### [molfidget.py](molfidget/molfidget.py)
- CLIエントリーポイント
- 3つのサブコマンド: `convert`, `preview`, `generate`
- argcompleteによるシェル補完サポート

#### [molecule.py](molfidget/molecule.py)
- `Molecule`クラス: 分子全体を管理
- trimeshシーンの生成と原子のグルーピング
- STL/3MFファイルの出力

#### [atom.py](molfidget/atom.py)
- `Atom`クラス: 個別の原子を表現
- ファン・デル・ワールス半径に基づく球体メッシュ生成
- CPKカラーリング対応

#### [bond.py](molfidget/bond.py)
- `Bond`クラス: 原子間の結合と機械的ジョイントを表現
- ジョイントタイプ:
  - `shaft_spin`: 回転軸（単結合用、Dカット無し、ストッパー付き）
  - `shaft`: 固定軸（Dカット無し）
  - `shaft_dcut`: 固定軸（Dカット有り、二重・三重結合のデフォルト）
  - `hole`: 丸穴
  - `hole_dcut`: Dカット穴（二重・三重結合のデフォルト）
  - `none`: 形状なし
- `sculpt_atoms2()`: 原子形状にジョイント形状を切り出す

#### [shape.py](molfidget/shape.py)
- `Shape`クラス: 機械的ジョイントの幾何形状を定義
- 原子間の交差平面計算
- テーパー、面取り、Dカットなどのパラメータ処理

#### [config.py](molfidget/config.py)
- データクラス定義: `MoleculeConfig`, `AtomConfig`, `BondConfig`
- PDB/MOLファイルパーサー
- MLF（YAML）形式の読み書き

#### [constants.py](molfidget/constants.py)
- 元素ごとのファン・デル・ワールス半径テーブル
- 元素ごとのCPKカラーテーブル
- 結合距離テーブル（単結合/二重結合/三重結合）

#### [labeled_scene_viewer.py](molfidget/labeled_scene_viewer.py)
- pygletベースの3Dビューアー
- 原子名ラベル表示機能

## MLF（Molfidget）ファイル形式

MLFはYAML形式の中間表現ファイル（拡張子: `.mlf`）。

### 主要構造
```yaml
name: 分子名
scale: スケール [mm/Angstrom]
default:
  atom:
    scale: 原子スケール
  bond:
    bond_gap_mm: 結合ギャップ [mm]
    shaft_radius: 軸半径 [Angstrom]
    # その他のデフォルトパラメータ
atoms:
  - name: 元素名_ID  # 例: C_1, H_2
    position: [x, y, z]  # Angstrom単位
    scale: 原子スケール（オプション）
    color: [r, g, b]  # オプション
bonds:
  - atom_pair: [atom1, atom2]
    bond_type: single|double|triple
    shape_types: [shaft_spin, hole]  # オプション
    shape_pair:  # 個別パラメータ（オプション）
      bond_gap_mm: 0.2
      shaft_radius: 0.3
```

### 形状タイプ
- `shaft_spin`: 回転軸（Dカット無し、ストッパー付き）
- `shaft`: 固定軸（Dカット無し）
- `shaft_dcut`: 固定軸（Dカット有り）
- `hole`: 丸穴
- `hole_dcut`: Dカット穴
- `none`: 形状なし

## 重要パラメータ

### スケール
- 元データはÅ（オングストローム）単位
- `--scale 1.0`だとSTLモデルが極小になる
- 3Dプリント用には`--scale 10.0`が推奨

### 軸ギャップ
- `--shaft-gap`: 可動軸の軸と穴のクリアランス（mm単位）
- デフォルト: 0.2mm
- 一般的な3Dプリンタでは0.2〜0.3mmが適切

### ジョイントパラメータ（MLFファイル内）
- `bond_gap_mm`: 結合平面の隙間 [mm]
- `shaft_radius`: 軸半径 [Angstrom]
- `shaft_length`: 軸長 [Angstrom]
- `shaft_chamfer`: 面取り長 [Angstrom]
- `shaft_spin_stopper_radius`: ストッパー半径 [Angstrom]
- `stopper_length`: ストッパー長 [Angstrom]
- `hole_radius`: 穴半径 [Angstrom]
- `hole_depth`: 穴深さ [Angstrom]

詳細は[docs/images/modeling-2.png](docs/images/modeling-2.png)を参照。

## 出力ファイル

### generate コマンドの出力
```
output/
├── <分子名>.3mf              # 全体の3MFファイル
├── <原子名>.stl              # 個別原子のSTL（例: C_1.stl）
└── <元素名>_group.stl        # 元素ごとのマージ済みSTL（例: C_group.stl）
```

## コード規約

### 単位系
- **ファイル内部の座標**: Angstrom（Å）
- **ジョイントギャップ系**: mm
- **最終STL出力**: scaleパラメータで変換

### デフォルト動作
- 単結合（single）: 最初の原子に`shaft_spin`、次の原子に`hole`
- 二重・三重結合（double/triple）: 最初の原子に`shaft_dcut`、次の原子に`hole_dcut`
- これらはMLFファイルの`shape_pair`内の`shape_type`で上書き可能

### trimeshとmanifold3d
- `trimesh`: 3Dメッシュ操作ライブラリ（v4.6.9）
- `manifold3d`: ブーリアン演算用（v3.1.0）
- ジョイント形状の切り出しに使用

## サンプルデータ

### 小分子（テスト用）
- [data/mlf/hydrogen.mlf](data/mlf/hydrogen.mlf): 水素分子（H₂）
- [data/mlf/h2o.mlf](data/mlf/h2o.mlf): 水分子
- [data/mlf/carbon_dioxide.mlf](data/mlf/carbon_dioxide.mlf): 二酸化炭素
- [data/mlf/ethanol.mlf](data/mlf/ethanol.mlf): エタノール

### やや複雑な分子
- [data/mlf/cyclohexane.mlf](data/mlf/cyclohexane.mlf): シクロヘキサン
- [data/pdb/ibuprofen.pdb](data/pdb/ibuprofen.pdb): イブプロフェン
- [data/pdb/18-Crown-6.pdb](data/pdb/18-Crown-6.pdb): クラウンエーテル

## 配布

### PyPI公開
- GitHub Actionsワークフロー: [.github/workflows/pypi-publish.yml](.github/workflows/pypi-publish.yml)
- PyPIパッケージ名: `molfidget`
- 現在のバージョン: 0.1.1（[pyproject.toml](pyproject.toml)を参照）

### インストール
```bash
pip install molfidget
```

## ドキュメント

すべてのドキュメントは日本語:
- [README.md](README.md): クイックスタート
- [docs/index.md](docs/index.md): プロジェクト概要
- [docs/mlf_file.md](docs/mlf_file.md): MLF形式仕様
- [docs/files.md](docs/files.md): ファイル構造とワークフロー
