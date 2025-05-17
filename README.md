# molfidget

## Quick start

### コードを取得

```
cd ~/
git clone https://github.com/longjie0723/molfidget.git
```

### poetryをインストール

```
sudo apt install python3-poetry
```

### パッケージをインストール

```
cd ~/molfidget
poetry install
```

### 実行

```
poetry run molfidget --scale 10.0 pdb/ethanol.pdb
```

* プレビューは'q'キーで抜ける
* 対応するstlファイルが直下に生成される