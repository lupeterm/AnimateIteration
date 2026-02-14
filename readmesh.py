import os
from classes import Mesh, Face, Boundary, Cell
import numpy as np
import assembly


def getAmountAndStart(file: str) -> tuple[int, int]:
    with open(file) as f:
        for j, line in enumerate(f.readlines()):
            try:
                amount = int(line)
                return amount, j + 2
            except ValueError:
                ...
    return 0, 0


def readPointsFile(polyMeshDir: str) -> list[np.ndarray]:
    pointsFileName = os.path.join(polyMeshDir, "points")
    if not os.path.isfile(pointsFileName):
        print(f"Case Directory '{pointsFileName}' does not exist")
        return []
    _, start = getAmountAndStart(pointsFileName)
    centroids: list[np.ndarray] = []
    i = 0
    with open(pointsFileName) as pointsFile:
        for j, line in enumerate(pointsFile.readlines()):
            if j < start:
                continue
            if line == ")\n":
                break
            s = line[1:-2].split(" ")
            centroids.append(np.array([float(ss) for ss in s]))
            i += 1
    return centroids


def readOwnersFile(polyMeshDir: str) -> list[int]:
    ownersFileName = os.path.join(polyMeshDir, "owner")
    if not os.path.isfile(ownersFileName):
        print(f"Owners file '{ownersFileName}' does not exist")
        return []
    nOwners, start = getAmountAndStart(ownersFileName)
    owners: list[int] = [0] * nOwners
    i = 0
    with open(ownersFileName) as ownersFile:
        for j, line in enumerate(ownersFile.readlines()):
            if j < start:
                continue
            if line == ")\n":
                break
            owners[i] = int(line)
            i += 1
    return owners


def readFacesFile(polyMeshDir: str, owners: list[int]) -> list[Face]:
    facesFileName = os.path.join(polyMeshDir, "faces")
    if not os.path.isfile(facesFileName):
        print(f"faces file '{facesFileName}' does not exist")
        return []
    _, start = getAmountAndStart(facesFileName)
    faces: list[Face] = []
    i = 0
    with open(facesFileName) as file:
        for j, line in enumerate(file.readlines()):
            if j < start:
                continue
            if line == ")\n":
                break
            iNodes = [int(ss) for ss in line[2:-2].split(" ")]
            faces.append(Face(index=i, iNodes=iNodes, iOwner=owners[i]))
            i += 1
    return faces


def readNeighborsFile(polyMeshDir: str, faces: list[Face]) -> tuple[int, list[Face]]:
    neighborsFileName = os.path.join(polyMeshDir, "neighbour")
    if not os.path.isfile(neighborsFileName):
        print(f"nb file '{neighborsFileName}' does not exist")
        return 0, faces
    numNeighbors, start = getAmountAndStart(neighborsFileName)
    i = 0
    with open(neighborsFileName) as file:
        for j, line in enumerate(file.readlines()):
            if j < start:
                continue
            if line == ")\n":
                break
            faces[i].iNeighbor = int(line)
            i += 1
    return numNeighbors, faces


def readBoundaryFile(polyMeshDir: str) -> list[Boundary]:
    boundariesFileName = os.path.join(polyMeshDir, "boundary")
    if not os.path.isfile(boundariesFileName):
        print(f"b file '{boundariesFileName}' does not exist")
        return []
    numBoundaries, _ = getAmountAndStart(boundariesFileName)
    i = 16
    boundaries: list[Boundary] = []
    with open(boundariesFileName) as file:
        lines = file.readlines()
        while lines[i] != "\n" or lines[i].startswith("//"):
            i += 1
        lines = lines[i + 1 : -4]
        l = "".join(lines).split("}")
        for iBound in range(numBoundaries):
            name = ""
            type = ""
            nFaces: int = 0
            startface: int = 0
            spl = l[iBound].split("{")
            name = spl[0].strip()
            arrayofpairs = spl[1].split(";")
            for pair in arrayofpairs:
                if pair.strip() == "":
                    continue
                pair = pair.lstrip().split()
                left = pair[0]
                right = pair[1]
                if left == "type":
                    type = right
                elif left == "inGroups":
                    ...
                elif left == "nFaces":
                    nFaces = int(right)
                elif left == "startFace":
                    startface = int(right)
                else:
                    print(f"unknown boundary key {left}")
            boundary = Boundary(name, type, nFaces, startface, iBound)
            boundaries.append(boundary)
    return boundaries


def constructCells(
    nodes: list[np.ndarray],
    boundaries: list[Boundary],
    faces: list[Face],
    numCells: int,
    numInteriorFaces: int,
) -> Mesh:
    cells = [Cell(index=i) for i in range(numCells + 1)]
    numInteriors = numInteriorFaces
    for InteriorFace in range(numInteriorFaces):
        iOwner = faces[InteriorFace].iOwner
        iNeighbor = faces[InteriorFace].iNeighbor
        if iNeighbor <= 0:
            # println("settings numInteriors to $numInteriors")
            numInteriors = InteriorFace - 1
            break
        cells[iOwner].iFaces.append(faces[InteriorFace].index)
        cells[iOwner].iNeighbors.append(iNeighbor)
        cells[iOwner].faceSigns.append(1)
        cells[iNeighbor].iFaces.append(faces[InteriorFace].index)
        cells[iNeighbor].iNeighbors.append(iOwner)
        cells[iNeighbor].faceSigns.append(-1)
        cells[iOwner].numIntFaces += 1
        cells[iNeighbor].numIntFaces += 1
    for boundaryFace in range(numInteriors + 1, len(faces)):
        owner = faces[boundaryFace].iOwner
        cells[owner].iFaces.append(boundaryFace)
        cells[owner].faceSigns.append(1)
    numBoundaryCells = len(faces) - numInteriors
    numBoundaryFaces = len(faces) - numInteriors
    return Mesh(
        nodes,
        faces,
        boundaries,
        numCells,
        cells,
        numInteriors,
        numBoundaryCells,
        numBoundaryFaces,
    )


def processBasicFaceGeometry(mesh: Mesh) -> Mesh:
    for face in mesh.faces:
        area = 0.0
        # special case: triangle
        if len(face.iNodes) == 3:
            ...
            # sum x,y,z and divide by 3
            # triangleNodes = [mesh.nodes[i] for i in face.iNodes]
            # face.centroid = sum(triangleNodes) / 3
            # face.Sf .= 0.5 * cross(triangleNodes[2] - triangleNodes[1], triangleNodes[3] - triangleNodes[1])
            # area = magnitude(face.Sf)
        else:  # general case, polygon is not a triangle
            nodes = np.array([mesh.nodes[i] for i in face.iNodes])
            center = sum(nodes) / len(nodes)
            # Using the center to compute the area and centroid of virtual
            # triangles based on the center and the face nodes
            triangleNode1 = center
            triangleNode3 = [0.0, 0.0, 0.0]
            # for iNodeIndex, iNode in enumerate(face.iNodes):
            for i, iNode in enumerate(face.iNodes):
                if i < len(face.iNodes)-1:
                    triangleNode3 = mesh.nodes[face.iNodes[i+1]]
                else:
                    triangleNode3 = mesh.nodes[face.iNodes[0]]
                # Calculate the centroid of a given subtriangle
                localCentroid = (triangleNode1 + mesh.nodes[iNode] + triangleNode3) / 3
                # Calculate the surface area vector of a given subtriangle by cross product
                localSf = 0.5 * np.cross(
                    mesh.nodes[iNode] - triangleNode1, triangleNode3 - triangleNode1
                )
                # Calculate the surface area of a given subtriangle
                localArea = np.sqrt(np.dot(localSf, localSf))
                face.centroid += localArea * localCentroid
                face.Sf += localSf
            area = float(np.linalg.norm(face.Sf))
            # Compute centroid of the polygon
            face.centroid /= area
        face.area = area
    return mesh


def computeElementVolumeAndCentroid(mesh: Mesh) -> Mesh:
    for iElement in range(len(mesh.cells)):
        iFaces = mesh.cells[iElement].iFaces
        elementCenter = np.zeros(3)
        for iFace in iFaces:
            elementCenter += mesh.faces[iFace].centroid
        elementCenter /= len(iFaces)
        localVolumeCentroidSum = np.zeros(3)
        localVolumeSum = 0.0
        for iFace in range(len(iFaces)):
            localFace = mesh.faces[iFaces[iFace]]
            localFaceSign = mesh.cells[iElement].faceSigns[iFace]
            Sf = localFaceSign * localFace.Sf
            d_Gf = localFace.centroid - elementCenter

            localVolume = (Sf[0] * d_Gf[0] + Sf[1] * d_Gf[1] + Sf[2] * d_Gf[2]) / 3.0
            localVolumeSum += localVolume

            localCentroid = 0.75 * localFace.centroid + 0.25 * elementCenter
            localVolumeCentroidSum += localCentroid * localVolume
        mesh.cells[iElement].centroid = (1 / localVolumeSum) * localVolumeCentroidSum
        mesh.cells[iElement].volume = localVolumeSum
        # mesh.cells[iElement].oldVolume = localVolumeSum
    return mesh


def processSecondaryFaceGeometry(mesh: Mesh) -> Mesh:
    for iFace in range(mesh.numInteriorFaces):
        theFace = mesh.faces[iFace]
        ownerElement = mesh.cells[theFace.iOwner]
        neighborElement = mesh.cells[theFace.iNeighbor]

        CN = neighborElement.centroid - ownerElement.centroid
        magCN = np.linalg.norm(CN)
        eCN = CN / magCN

        E = theFace.area * eCN
        theFace.gDiff = float(np.linalg.norm(E) / magCN)

    for iBFace in range(mesh.numInteriorFaces + 1, len(mesh.faces)):
        theBFace = mesh.faces[iBFace]
        ownerElement = mesh.cells[theBFace.iOwner]
        cn = theBFace.centroid - ownerElement.centroid
        mesh.faces[iBFace].gDiff = (
            theBFace.area * theBFace.area / np.dot(cn, theBFace.Sf)
        )
    return mesh


def labelBoundaryFaces(mesh: Mesh) -> Mesh:
    for iBoundary, boundary in enumerate(mesh.boundaries):
        startFace = boundary.startFace  # +1
        nBFaces = boundary.nFaces
        for iFace in range(startFace, startFace + nBFaces - 1):
            mesh.faces[iFace].patchIndex = iBoundary
    return mesh


def processOpenFoamMesh(mesh: Mesh) -> Mesh:
    mesh = processBasicFaceGeometry(mesh)
    mesh = computeElementVolumeAndCentroid(mesh)
    mesh = processSecondaryFaceGeometry(mesh)
    mesh = labelBoundaryFaces(mesh)
    return mesh


def readOpenFoamMesh(caseDir: str) -> Mesh:
    if not os.path.isdir(caseDir):
        print(f"Case Directory '{caseDir}' does not exist")
        return Mesh()
    polymeshDir: str = os.path.join(caseDir, "constant/polyMesh")
    if not os.path.isdir(polymeshDir):
        print(f"PolyMesh Directory '{caseDir}' does not exist")
        return Mesh()
    nodes = readPointsFile(polymeshDir)
    owner = readOwnersFile(polymeshDir)
    faces = readFacesFile(polymeshDir, owner)
    numNeighbors, faces = readNeighborsFile(polymeshDir, faces)
    numCells = max(owner)
    boundaries = readBoundaryFile(polymeshDir)
    mesh = constructCells(nodes, boundaries, faces, numCells, numNeighbors)
    mesh = processOpenFoamMesh(mesh)
    return mesh


if __name__ == "__main__":
    case = "/home/peter/Documents/uni/FVM-Prototyping/cases/square-3x3"
    mesh = readOpenFoamMesh(case)
    rows = assembly.CellBasedAssembly(mesh)
