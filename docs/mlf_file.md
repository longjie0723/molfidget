# MLF (Molfidget) ファイル仕様

- 拡張子: `.mlf`
- 形式: YAML
- 分子構造データの中間表現。PDB/MOLから変換され、3Dモデル生成に使用される。

## トップレベル構造

```yaml
default:   # 各種デフォルト値
  atom: ...
  bond: ...
molecule:  # 分子データ本体
  name: ...
  scale: ...
  atoms: [...]
  bonds: [...]
```

---

## default

各セクションのデフォルト値を定義する。個別のatom/bondで指定がない場合にこの値が使われる。

### default.atom

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| vdw_scale | float | 0.8 | ファン・デル・ワールス半径に対する形状のスケール |

### default.bond

共通パラメータと、bond_typeごとのデフォルト値を持つ。

#### 共通パラメータ

すべてのbond_typeに適用される。bond_typeごとのセクション内で上書き可能。

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| bond_gap_mm | float | 0.1 | 結合平面の隙間 [mm] |
| bond_marker | str | "hetero-only-except-ch" | ボンドマーカーの刻印モード |
| bond_marker_size_mm | float | 1.5 | ボンドマーカーのサイズ [mm] |
| bond_marker_depth_mm | float | 0.5 | ボンドマーカーの刻印深さ [mm] |
| taper_angle_deg | float | 0.0 | テーパー角度 [度] |
| taper_radius_scale | float | 1.0 | テーパー半径のスケール |

**bond_marker の値:**
- `off` — マーカーなし
- `on` — すべての結合にマーカー
- `hetero-only` — 異なる元素間の結合のみ
- `hetero-only-except-ch` — 異なる元素間の結合のみ（C-H結合を除く）

#### bond_typeごとのデフォルト値

各bond_typeのセクション内に、そのタイプ固有の形状パラメータのデフォルト値を定義する。

##### spin（回転軸 + 丸穴）

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| shaft_radius | float | 0.3 | シャフトの半径 [Å] |
| shaft_length | float | 0.3 | シャフトの長さ [Å] |
| hole_radius | float | 0.3 | 穴の半径 [Å] |
| hole_length | float | 0.3 | 穴の深さ [Å] |
| chamfer_length | float | 0.1 | 面取りの長さ [Å] |
| wall_thickness | float | 0.1 | 壁の厚さ [Å] |
| shaft_gap | float | 0.03 | シャフトとキャビティの隙間 [Å] |
| stopper_radius | float | 0.4 | ストッパーの半径 [Å] |
| stopper_length | float | 0.2 | ストッパーの長さ [Å] |

##### normal（丸軸 + 丸穴）

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| shaft_radius | float | 0.3 | シャフトの半径 [Å] |
| shaft_length | float | 0.3 | シャフトの長さ [Å] |
| hole_radius | float | 0.3 | 穴の半径 [Å] |
| hole_length | float | 0.3 | 穴の深さ [Å] |
| chamfer_length | float | 0.1 | 面取りの長さ [Å] |

##### fixed（Dカット軸 + Dカット穴）

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| shaft_radius | float | 0.3 | シャフトの半径 [Å] |
| shaft_length | float | 0.3 | シャフトの長さ [Å] |
| hole_radius | float | 0.3 | 穴の半径 [Å] |
| hole_length | float | 0.3 | 穴の深さ [Å] |
| chamfer_length | float | 0.1 | 面取りの長さ [Å] |

##### gapped（丸軸 + 丸穴、一体成型用に軸を太くする）

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| shaft_radius | float | 0.8 | シャフトの半径 [Å] |
| shaft_length | float | 0.3 | シャフトの長さ [Å] |
| hole_radius | float | 0.8 | 穴の半径 [Å] |
| hole_length | float | 0.3 | 穴の深さ [Å] |
| chamfer_length | float | 0.1 | 面取りの長さ [Å] |

##### short（短い丸軸 + 丸穴、はめ込みパーツ用）

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| shaft_radius | float | 0.3 | シャフトの半径 [Å] |
| shaft_length | float | 0.2 | シャフトの長さ [Å] |
| hole_radius | float | 0.3 | 穴の半径 [Å] |
| hole_length | float | 0.3 | 穴の深さ [Å] |
| chamfer_length | float | 0.1 | 面取りの長さ [Å] |

##### holes（丸穴 + 丸穴、磁石接続用）

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| hole_radius_mm | float | 3.525 | 穴の半径 [mm] |
| hole_length_mm | float | 2.0 | 穴の深さ [mm] |
| bond_gap_mm | float | 0.0 | 結合平面の隙間 [mm] |

##### notch_2 / notch_3（回転軸 + 丸穴 + ノッチ付き）

spinと同じデフォルト値。ノッチの組数が名前で決まる（notch_2=2組、notch_3=3組）。

| パラメータ | 型 | デフォルト | 説明 |
|---|---|---|---|
| shaft_radius | float | 0.3 | シャフトの半径 [Å] |
| shaft_length | float | 0.3 | シャフトの長さ [Å] |
| hole_radius | float | 0.3 | 穴の半径 [Å] |
| hole_length | float | 0.3 | 穴の深さ [Å] |
| chamfer_length | float | 0.1 | 面取りの長さ [Å] |
| wall_thickness | float | 0.1 | 壁の厚さ [Å] |
| shaft_gap | float | 0.03 | シャフトとキャビティの隙間 [Å] |
| stopper_radius | float | 0.4 | ストッパーの半径 [Å] |
| stopper_length | float | 0.2 | ストッパーの長さ [Å] |

##### plane（平面スライスのみ、形状なし）

固有のパラメータなし。共通パラメータのbond_gap_mmのみ使用。

#### default.bond の YAML 例

```yaml
default:
  bond:
    # 共通パラメータ
    bond_gap_mm: 0.1
    bond_marker: hetero-only-except-ch
    bond_marker_size_mm: 1.5
    bond_marker_depth_mm: 0.5
    taper_angle_deg: 0.0
    taper_radius_scale: 1.0
    # bond_typeごとのデフォルト
    spin:
      shaft_radius: 0.3
      shaft_length: 0.3
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
      wall_thickness: 0.1
      shaft_gap: 0.03
      stopper_radius: 0.4
      stopper_length: 0.2
    normal:
      shaft_radius: 0.3
      shaft_length: 0.3
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
    fixed:
      shaft_radius: 0.3
      shaft_length: 0.3
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
    gapped:
      shaft_radius: 0.8
      shaft_length: 0.3
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
    short:
      shaft_radius: 0.3
      shaft_length: 0.2
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
    holes:
      hole_radius_mm: 3.525
      hole_length_mm: 2.0
      bond_gap_mm: 0.0
    notch_2:
      shaft_radius: 0.3
      shaft_length: 0.3
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
      wall_thickness: 0.1
      shaft_gap: 0.03
      stopper_radius: 0.4
      stopper_length: 0.2
    notch_3:
      shaft_radius: 0.3
      shaft_length: 0.3
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
      wall_thickness: 0.1
      shaft_gap: 0.03
      stopper_radius: 0.4
      stopper_length: 0.2
```

---

## molecule

分子全体を表すオブジェクト。

| パラメータ | 型 | 説明 |
|---|---|---|
| name | str | 分子の名前 |
| scale | float | 3Dモデルのスケール [mm/Å]（推奨: 10.0） |

---

## atoms

atomの配列。

### atom

1つの原子を表すオブジェクト。

| パラメータ | 型 | 必須 | 説明 |
|---|---|---|---|
| name | str | Yes | 原子の名前: `元素名_連番`（例: `C_1`, `H_2`）。ファイル内で一意 |
| position | [x,y,z] | Yes | 原子の位置 [Å] |
| vdw_scale | float | No | ファン・デル・ワールス半径スケール（default.atom.vdw_scaleを上書き） |
| color | [R,G,B,A] | No | 原子の色（default.atom.colorを上書き） |

---

## bonds

bondの配列。

### bond

2つのatom間の結合を表すオブジェクト。

| パラメータ | 型 | 必須 | 説明 |
|---|---|---|---|
| atom_pair | [str, str] | Yes | 結合する原子名のペア。前=軸側、後=穴側（デフォルト） |
| bond_type | str | Yes | 結合タイプ |
| bond_gap_mm | float | No | 結合平面の隙間 [mm] |
| bond_marker | str | No | ボンドマーカーモード |
| bond_marker_size_mm | float | No | ボンドマーカーサイズ [mm] |
| bond_marker_depth_mm | float | No | ボンドマーカー深さ [mm] |
| shape_pair | [ShapeConfig, ShapeConfig] | No | 個別の形状パラメータ（後述） |

### bond_type

結合のタイプを指定する。各タイプにデフォルトの形状（shape_type）が割り当てられる。

| bond_type | atom1側（前） | atom2側（後） | 説明 |
|---|---|---|---|
| spin | shaft_spin | hole | 回転軸 + 丸穴 |
| normal | shaft | hole | 丸軸 + 丸穴 |
| fixed | shaft_dcut | hole_dcut | Dカット軸 + Dカット穴 |
| gapped | shaft | hole | 丸軸 + 丸穴（一体成型用、軸太め） |
| short | shaft | hole | 短い丸軸 + 丸穴（はめ込み用） |
| holes | hole | hole | 丸穴 + 丸穴（磁石接続用） |
| notch_2 | shaft_spin | hole | 回転軸 + 丸穴 + ノッチ2組 |
| notch_3 | shaft_spin | hole | 回転軸 + 丸穴 + ノッチ3組 |
| plane | none | none | 平面スライスのみ（形状なし） |

### shape_pair

atom_pairの各原子に対して個別に形状パラメータを設定するオプション。2要素の配列で、それぞれ以下のパラメータを持つ。bond_typeのデフォルト値を上書きする。

**mm単位のパラメータが指定されている場合、Å単位より常に優先されます。**

| パラメータ | 型 | 説明 |
|---|---|---|
| shape_type | str | 形状タイプ（bond_typeのデフォルトを上書き） |
| shaft_radius | float | シャフトの半径 [Å] |
| shaft_length | float | シャフトの長さ [Å] |
| shaft_radius_mm | float | シャフトの半径 [mm] |
| shaft_length_mm | float | シャフトの長さ [mm] |
| hole_radius | float | 穴の半径 [Å] |
| hole_length | float | 穴の深さ [Å] |
| hole_radius_mm | float | 穴の半径 [mm] |
| hole_length_mm | float | 穴の深さ [mm] |
| chamfer_length | float | 面取りの長さ [Å] |
| wall_thickness | float | 壁の厚さ [Å] |
| shaft_gap | float | シャフトとキャビティの隙間 [Å] |
| shaft_gap_mm | float | シャフトとキャビティの隙間 [mm] |
| stopper_radius | float | ストッパーの半径 [Å] |
| stopper_length | float | ストッパーの長さ [Å] |
| bond_gap_mm | float | 結合平面の隙間 [mm] |
| taper_radius_scale | float | テーパー半径のスケール |
| taper_angle_deg | float | テーパー角度 [度] |
| taper_distance | float | テーパーの距離 [Å] |
| taper_height | float | テーパーの高さ [Å] |

### shape_type 一覧

| shape_type | 説明 |
|---|---|
| shaft_spin | 回転軸（Dカットなし） |
| shaft | 丸軸（Dカットなし） |
| shaft_dcut | Dカット軸 |
| hole | 丸穴 |
| hole_dcut | Dカット穴 |
| none | 形状なし |

---

## パラメータ解決の優先順位

```
shape_pair個別値 > bond個別値 > default.bond.{bond_type}値 > default.bond共通値
```

mm単位（`*_mm`）が指定されている場合、同名のÅ単位パラメータより常に優先される。

---

## ファイル例

```yaml
default:
  atom:
    vdw_scale: 0.8
  bond:
    bond_gap_mm: 0.1
    bond_marker: hetero-only-except-ch
    taper_angle_deg: 0.0
    taper_radius_scale: 1.0
    spin:
      shaft_radius: 0.3
      shaft_length: 0.3
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
      wall_thickness: 0.1
      shaft_gap: 0.03
      stopper_radius: 0.4
      stopper_length: 0.2
    normal:
      shaft_radius: 0.3
      shaft_length: 0.2
      hole_radius: 0.3
      hole_length: 0.3
      chamfer_length: 0.1
    holes:
      hole_radius_mm: 3.525
      hole_length_mm: 2.0
      bond_gap_mm: 0.0
molecule:
  name: ethanol
  scale: 8.0
  atoms:
  - name: C_1
    position: [0.445, -2.202, -5.28]
  - name: O_4
    position: [1.295, -3.437, -5.31]
  - name: H_5
    position: [1.85, -3.479, -4.479]
  bonds:
  - atom_pair: [C_1, O_4]
    bond_type: spin
  - atom_pair: [O_4, H_5]
    bond_type: spin
    bond_gap_mm: 0.3
    shape_pair:
    - taper_radius_scale: 0.3
      taper_angle_deg: 5.0
    - taper_radius_scale: 0.3
      taper_angle_deg: 5.0
```
