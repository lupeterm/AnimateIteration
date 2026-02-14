from dataclasses import dataclass, field
import numpy as np

@dataclass
class Face:
    index: int
    iNodes: list[int]
    iOwner: int
    iNeighbor: int = -1
    centroid: np.ndarray = field(default_factory=lambda: np.zeros(3))
    Sf: np.ndarray = field(default_factory=lambda: np.zeros(3))
    area: float = 0.0
    gDiff: float = 0.0
    patchIndex: int = -1


@dataclass
class Boundary:
    name: str
    type: str
    nFaces: int
    startFace: int
    index: int


@dataclass
class Cell:
    index: int = 0
    numIntFaces: int = 0
    volume: float = 0.0
    iFaces: list[int] = field(default_factory=lambda: [])
    iNeighbors: list[int]  = field(default_factory=lambda: [])
    faceSigns: list[int]  = field(default_factory=lambda: [])
    centroid: np.ndarray = field(default_factory=lambda: np.zeros(3))

    # iNodes:list[int]
    # numNeighbors:int
    # oldVolume:float


@dataclass
class Mesh:
    nodes: list[np.ndarray] = field(default_factory=lambda: [])
    faces: list[Face] = field(default_factory=lambda: [])
    boundaries: list[Boundary]= field(default_factory=lambda: [])
    numCells: int = 0
    cells: list[Cell] = field(default_factory=lambda: [])
    numInteriorFaces: int = 0
    numBoundaryCells: int = 0
    numBoundaryFaces: int = 0


@dataclass
class Field:
    values: list[list[float]]


@dataclass
class BoundaryField:
    name: str
    nFaces: int
    values: list[list[float]]
    type: str


@dataclass
class LdcMatrixAssemblyInput:
    mesh: Mesh
    nu: list[float]
    U: tuple[list[BoundaryField], Field]
