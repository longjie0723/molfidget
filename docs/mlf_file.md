# Molfidget File

- 拡張子は mlf

## molecule
- 分子全体を表すオブジェクト
- 複数のatomを持つ
- 複数のbondを持つ

name: 分子の名前
scale: 3Dモデルのスケール [mm/Angstrom]

### default
各セクションのデフォルト値

## atoms
- atomの配列を持つ

# atom
- 1つの原子を表すオブジェクト

name: 原子の名前: 元素名 + '_' + ID
IDはファイルを通して一意の自然数の番号である必要がある

scale: 原子のファン・デル・ワールス半径に対する実際の形状のスケール

position: 原子の位置: [x, y, z] 座標 [単位: Å]

color: 原子の色。デフォルト以外にしたい場合に指定

## bonds
- bondの配列を持つ

## bond
- 2つのatom同士の間の結合を表すオブジェクト

### atom_pair:

結合する原子のペアを名前で指定

- デフォルトで、前に来たatomが軸形状、後ろに来たatomが穴側
- ただし、shape_typesで変更できる
- どちらも軸、どちらも穴とかも変更可能

### bond_type
single, double, tripleなどの化学結合を表す
これらは化学結合を示すので形状を指定するものではない
デフォルトの形状が決まっている

singleの場合、最初に来たatom側にshaft_spin、後に来たatom側にholeを作成
double, tripleの場合、前に来たatom側にshaft_dcut、後に来たatom側にhole_dcutを作成

これらのデフォルトの形状はshape_typeなどで上書きされる

### shape_pair

atom_pairに個別に形状のパラメータを設定する場合のオプション。

- 形状に関係ないパラメータは無視される

おしりにmmとあるのはmm単位。何もないのはÅ単位。

bond_gap_mm: # 結合平面の隙間 [mm]

wall_thickness: 0.1 # 壁の厚さ [Angstrom] どこ？
shaft_radius: 0.3   # シャフトの半径
shaft_length: 0.3   # シャフトの長さ
shaft_chamfer: 0.1 # 軸の面取り
shaft_spin_gap_mm: 0.3   # スピン軸の隙間 [mm]
shaft_key_thickness: キー付き軸のキーの厚み

shaft_spin_stopper_radius: 0.4 # ストッパーの半径
stopper_length: 0.2 # ストッパーの長さ
hole_radius: 軸穴の半径
hole_depth: 軸穴の深さ

shaft_ball_stopper_radius: ボールの半径
shaft_ball_gap_mm: ボールと空間の隙間 [mm]

shape_d_cut: [False, False]: Dカットの抑制

### shape_types
形状のタイプ

- shaft_spin（回転軸、Dカット無し）
- shaft（固定軸、Dカット無し）
- shaft_dcut（固定軸、Dカット有り）
- hole（丸穴）
- hole_dcut（Dカット穴）
- none（形状なし）

軸の組み立て方向を指定したい場合は、shaft_dcutとhole_dcutを使用する。
