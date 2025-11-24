# ファイル構成

## フォルダ構成


### molfidget/

pythonコード

molfidget.py コマンド本体

molecule.py 分子のクラス定義
atom.py 原子のクラス定義
bond.py 結合のクラス定義
config.py 各オブジェクト作成時の設定クラスを定義
constants.py 確定数を定義
labeled_scene_viewer.py ラベル付きのViewerクラス

### data/

サンプルデータ用フォルダ

data/mlf/ molfidgetファイル

data/pdb/ pdbファイル

data/mol/ molファイル

data/cml/ cmlファイル

### doc

各種ドキュメント


###

手順
- 各オブジェクトの生成と初期化

1 Atomのループ
 - 球の生成
 - モデルはAtomが保持

2 Bondのループ
 - 2つのAtomの　切断面の整形
 - テーパの整形

3 Bondのループ
 - atomの軸、穴形状の整形