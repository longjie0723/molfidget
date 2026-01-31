import numpy as np
import pyglet
import trimesh
from pyglet.gl import *
from pyglet.text import Label
from trimesh.viewer import SceneViewer


class LabeledSceneViewer(SceneViewer):
    def __init__(self, scene):
        """
        geometryの名前をラベルとしてtrimesh.sceneを表示するビューアー
        """
        super().__init__(scene, start_loop=False)

        # Pygletのラベル初期化
        self.labels = {}
        for name, geom in self.scene.geometry.items():
            label = Label(
                name,
                font_name="Arial",
                font_size=16,
                color=(0, 0, 0, 255),  # 黒、不透明
                anchor_x="center",
                anchor_y="center",
            )
            # 各ジオメトリの中心位置をラベルの位置に設定
            # position = geom.bounding_box.centroid
            position = geom.centroid
            self.labels[name] = label, position

    def world_to_screen(self, point):
        """
        3Dワールド座標 point をスクリーン座標(x, y)に変換
        """
        mv = (GLdouble * 16)()
        glGetDoublev(GL_MODELVIEW_MATRIX, mv)
        pj = (GLdouble * 16)()
        glGetDoublev(GL_PROJECTION_MATRIX, pj)
        vp = (GLint * 4)()
        glGetIntegerv(GL_VIEWPORT, vp)
        win_x = GLdouble()
        win_y = GLdouble()
        win_z = GLdouble()
        # 出力変数を渡す
        gluProject(
            GLdouble(point[0]),
            GLdouble(point[1]),
            GLdouble(point[2]),
            mv,
            pj,
            vp,
            win_x,
            win_y,
            win_z,
        )
        return win_x.value, win_y.value, win_z.value

    def draw_labels(self):
        # 2D描画モードに切り替える前にラベルのスクリーン座標を計算
        label_positions = {}
        for name, (label, position) in self.labels.items():
            screen_x, screen_y, _ = self.world_to_screen(position)
            label.x = screen_x
            label.y = screen_y

        # 2D描画モードに切り替え
        glMatrixMode(GL_PROJECTION)
        glPushMatrix()
        glLoadIdentity()
        gluOrtho2D(0, self.width, 0, self.height)
        glMatrixMode(GL_MODELVIEW)
        glPushMatrix()
        glLoadIdentity()
        # ラベル描画
        for name, (label, _) in self.labels.items():
            label.draw()

        # 3D描画モードに戻す
        glMatrixMode(GL_PROJECTION)
        glPopMatrix()
        glMatrixMode(GL_MODELVIEW)
        glPopMatrix()

    def on_draw(self):
        # wireframeの線の太さを細くしたい
        glLineWidth(2.0)
        super().on_draw()
        self.draw_labels()


# 使用例
if __name__ == "__main__":
    # サンプルシーン作成
    scene = trimesh.Scene()

    # サンプルオブジェクト（立方体）
    cube = trimesh.creation.box(extents=(1, 1, 1))
    cube.apply_translation([0, 0, 0])  # 中心に配置
    scene.add_geometry(cube, geom_name="Cube")

    sphere = trimesh.creation.icosphere(subdivisions=2, radius=0.5)
    sphere.apply_translation([1.5, 0, 0])  # 立方体の右側に配置
    scene.add_geometry(sphere, geom_name="Sphere")

    cylinder = trimesh.creation.cylinder(radius=0.2, height=1.0)
    cylinder.apply_translation([-1.5, 0, 0])  # 立方体の左側に配置
    scene.add_geometry(cylinder, geom_name="Cylinder")

    # ビューアー起動
    viewer = LabeledSceneViewer(scene)
    pyglet.app.run()