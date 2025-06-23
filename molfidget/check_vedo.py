from vedo import Cube, Text3D, show, Sphere

# オブジェクトを作成
#obj = Cube().pos(0, 0, 0).c("cyan")

s = Sphere().pos(0, 0.5, 0).c("cyan")

# ラベルをオブジェクトの横に配置（x方向にオフセット）
#label_pos = obj.center_of_mass() + [0, 0.5, 0]  # x方向に1.2だけずらす
#label = Text3D("Cube", pos=label_pos, s=0.2, c="black")

s_pos = s.pos() + [0, 1, 0]  # x方向に1.2だけずらす
label_s = Text3D("Sphere", pos=s_pos, s=0.2, c="black")
label_s = label_s.follow_camera(True)  # カメラの動きに合わせてラベルを回転
# 表示
show([s, label_s], axes=1)

