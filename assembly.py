from classes import Mesh
import numpy as np
def CellBasedAssembly(mesh:Mesh) -> None:
    nCells = len(mesh.cells)
    idx = 0
    entriesNeeded = nCells + 2 * mesh.numInteriorFaces
    rows = np.zeros(entriesNeeded)
    cols = np.zeros(entriesNeeded)
    cellIdx = -1
    for theElement in mesh.cells:
        iElement = theElement.index
        isset = False
        numFaces = len(theElement.iFaces)
        for iFace in range(theElement.numIntFaces):
            iFaceIndex = theElement.iFaces[iFace]
            theFace = mesh.faces[iFaceIndex]
            if theElement.iNeighbors[iFace] > theElement.index and not isset:
                cellIdx = idx
                idx += 1
                isset = True
            cols[idx] = iElement
            rows[idx] = theElement.iNeighbors[iFace]
            idx += 1
        if not isset:
            cellIdx = idx
        for iFace in range(theElement.numIntFaces+1,numFaces):
            iFaceIndex = theElement.iFaces[iFace]
            # boundaryType = velocity_boundary[iBoundary].type
            # if boundaryType != "fixedValue"
            #     continue
            # end
            # theFace = mesh.faces[iFaceIndex]
            # fluxCb = nu[theFace.iOwner] * theFace.gDiff
            # relativeFaceIndex = iFaceIndex - mesh.boundaries[iBoundary].startFace
            # fluxVb::Vector{Float32} = velocity_boundary[iBoundary].values[relativeFaceIndex] .* -fluxCb
            # RHS[iElement] -= fluxVb[1]
            # RHS[iElement+nCells] -= fluxVb[2]
            # RHS[iElement+nCells+nCells] -= fluxVb[3]
            # diagx += fluxCb
            # diagy += fluxCb
            # diagz += fluxCb
        cols[cellIdx] = iElement
        rows[cellIdx] = iElement
    print(rows)
    return rows