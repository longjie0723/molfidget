# molfidget

## Quick start

### 実行環境

* Ubuntu
* もしくは Windows WSL2 

### 必要なツール

```
sudo apt install git libglut-dev python3-poetry
```

### コードを取得

```
cd ~/
git clone https://github.com/longjie0723/molfidget.git
```

### パッケージをインストール

```
cd ~/molfidget
poetry install
```

### 実行

```
poetry run molfidget --scale 10.0 --shaft-gap 0.2 pdb/ethanol.pdb
```

* プレビューは'q'キーで抜ける
* 対応するstlファイルが直下に生成される

### その他

* 元のモデルは単位がオングストロームになっているのでscale=1.0だとSTLファイルのモデルはすごく小さくなってしまう
* なのでscale=10.0とかにするとプリント可能なSTLになる
* 可動軸の軸と穴のギャップは `--shaft-gap`で指定する。単位はmmで、デフォルトは0.2になっている
* 通常の3Dプリンタだとたぶん0.2~0.3ぐらいの間で問題ない

### モデリング

軸、穴部分の形状とパラメータの名前はこのようになっている。
![モデリング図解](image/modeling-1.png)