from manim import *
from classes import Cell
from manim import color
import numpy as np


cellcentroids = [
    [-1, 1, 0],
    [0, 1, 0],
    [1, 1, 0],
    [-1, 0, 0],
    [0, 0, 0],
    [1, 0, 0],
    [-1, -1, 0],
    [0, -1, 0],
    [1, -1, 0],
]
internalfaces = {5, 6, 8, 9, 10, 12, 13, 15, 16, 17}
bfaces = [1, 2, 3, 4, 7, 11, 14, 18, 21, 22, 23, 24]
cellfacemapping = {
    0: [1, 4, 5, 8],
    1: [2, 5, 6, 9],
    2: [3, 6, 7, 10],
    3: [8, 11, 12, 15],
    4: [9, 12, 13, 16],
    5: [10, 13, 14, 17],
    6: [15, 18, 19, 22],
    7: [16, 19, 20, 23],
    8: [17, 20, 21, 24],
}
numintmapping = {0: 2, 1: 3, 2: 2, 3: 3, 4: 4, 5: 3, 6: 2, 7: 3, 8: 2}
meshcells = [
    Cell(
        i,
        numIntFaces=numintmapping[i],
        iFaces=cellfacemapping[i],
        centroid=np.array(cellcentroids[i]),
    )
    for i in range(9)
]

cell_faces = [
    [[-1.5, 1.5, 0.0], [-0.5, 1.5, 0.0]],
    [[-0.5, 1.5, 0.0], [0.5, 1.5, 0.0]],
    [[0.5, 1.5, 0.0], [1.5, 1.5, 0.0]],
    [[-1.5, 1.5, 0.0], [-1.5, 0.5, 0.0]],
    [[-0.5, 1.5, 0.0], [-0.5, 0.5, 0.0]],
    [[0.5, 1.5, 0.0], [0.5, 0.5, 0.0]],
    [[1.5, 1.5, 0.0], [1.5, 0.5, 0.0]],
    [[-1.5, 0.5, 0.0], [-0.5, 0.5, 0.0]],
    [[-0.5, 0.5, 0.0], [0.5, 0.5, 0.0]],
    [[0.5, 0.5, 0.0], [1.5, 0.5, 0.0]],
    [[-1.5, 0.5, 0.0], [-1.5, -0.5, 0.0]],
    [[-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0]],
    [[0.5, 0.5, 0.0], [0.5, -0.5, 0.0]],
    [[1.5, 0.5, 0.0], [1.5, -0.5, 0.0]],
    [[-1.5, -0.5, 0.0], [-0.5, -0.5, 0.0]],
    [[-0.5, -0.5, 0.0], [0.5, -0.5, 0.0]],
    [[0.5, -0.5, 0.0], [1.5, -0.5, 0.0]],
    [[-1.5, -0.5, 0.0], [-1.5, -1.5, 0.0]],
    [[-0.5, -0.5, 0.0], [-0.5, -1.5, 0.0]],
    [[0.5, -0.5, 0.0], [0.5, -1.5, 0.0]],
    [[1.5, -0.5, 0.0], [1.5, -1.5, 0.0]],
    [[-1.5, -1.5, 0.0], [-0.5, -1.5, 0.0]],
    [[-0.5, -1.5, 0.0], [0.5, -1.5, 0.0]],
    [[0.5, -1.5, 0.0], [1.5, -1.5, 0.0]],
]


class FaceBased(Scene):
    def construct(self):
        faces = VGroup()
        faceid = 0
        faceMapping: dict[int, LabeledLine] = dict()
        for i, pair in enumerate(cell_faces):
            line = LabeledLine(
                label=f"{i}",
                label_position=0.5,
                label_config={"font_size": 20},
                start=pair[0],
                end=pair[1],
            )
            faceMapping[i] = line
            faces.add(line)
            faceid += 1
        self.add(faces)
        cells = VGroup()
        for i, centroid in enumerate(cellcentroids):
            t = MathTex(rf"C_{i}").move_to(np.array(centroid))
            cells.add(t)
        self.add(cells)
        all = VGroup(faces, cells)
        all.scale(2.0)

        for iFace in range(len(cell_faces)):
            self.play(
                # faceMapping[last][0].animate.set_color(WHITE),
                # faceMapping[last][1][1].animate.set_color(WHITE),
                # faceMapping[last][1][0].animate.set_color(BLACK),
                faceMapping[iFace][0].animate.set_color(GREEN),
                faceMapping[iFace][1][1].animate.set_color(WHITE),
            )



class CellBased(Scene):
    def construct(self):
        faces = VGroup()
        faceid = 0
        faceMapping: dict[int, LabeledLine] = dict()
        for i, pair in enumerate(cell_faces):
            line = LabeledLine(
                label=f"{i}",
                label_position=0.5,
                label_config={"font_size": 20},
                start=pair[0],
                end=pair[1],
            )
            faceMapping[i] = line
            faces.add(line)
            faceid += 1
        self.add(faces)
        cells = VGroup()
        squares = VGroup()
        cellsquares: dict[int, Square] = dict()
        for i, centroid in enumerate(cellcentroids):
            s = Square(side_length=0.98)
            s.move_to(np.array(centroid))
            s.set_opacity(0.0)
            s.set_fill(ORANGE)

            s.set_z_index(-1.0)
            t = MathTex(rf"C_{i}").move_to(np.array(centroid))
            # t.set_color(BLACK)
            cells.add(t)
            squares.add(s)
            cellsquares[i] = s
        self.add(cells)
        self.add(squares)
        all = VGroup(faces, cells, squares)
        all.scale(2.0)
        for line in faceMapping.values():
            line.set_z_index(1000)
        for i, cell in enumerate(meshcells):
            numints = cell.numIntFaces
            cellsquare = cellsquares[i]
            cellsquare.set_opacity(0.4)

            for iFace in cell.iFaces:
                if iFace in bfaces:
                    continue
                self.play(
                    faceMapping[iFace-1].animate.set_color(BLUE),
                    faceMapping[iFace-1][1][1].animate.set_color(WHITE),
                )
            for iFace in cell.iFaces:
                if iFace not in bfaces:
                    continue
                self.play(
                    faceMapping[iFace-1].animate.set_color(GREEN),
                    faceMapping[iFace-1][1][1].animate.set_color(WHITE),
                )    


class BatchedFaceBased(Scene):
    def construct(self):
        faces = VGroup()
        faceid = 0
        faceMapping: dict[int, LabeledLine] = dict()
        for i, pair in enumerate(cell_faces):
            line = LabeledLine(
                label=f"{i}",
                label_position=0.5,
                label_config={"font_size": 20},
                start=pair[0],
                end=pair[1],
            )
            faceMapping[i] = line
            faces.add(line)
            faceid += 1
        self.add(faces)
        cells = VGroup()
        for i, centroid in enumerate(cellcentroids):
            t = MathTex(rf"C_{i}").move_to(np.array(centroid))
            cells.add(t)
        self.add(cells)
        all = VGroup(faces, cells)
        all.scale(2.0)

        for iFace in range(len(cell_faces)):
            if iFace +1 in bfaces:
                continue
            self.play(
                faceMapping[iFace][0].animate.set_color(GREEN),
                faceMapping[iFace][1][1].animate.set_color(WHITE),
            )

        for iFace in range(len(cell_faces)):
            if iFace +1 not in bfaces:
                continue
            self.play(
                faceMapping[iFace][0].animate.set_color(BLUE),
                faceMapping[iFace][1][1].animate.set_color(WHITE),
            )
