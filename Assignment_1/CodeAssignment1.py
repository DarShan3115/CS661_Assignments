from vtk import *
import numpy as np
import vtkmodules.util.numpy_support as np_support

# Your input file path inside double quotes. Don't add more quotes in the start or end if there remove them.
filepath = r""

reader = vtkXMLImageDataReader()
reader.SetFileName(filepath)
reader.Update()

data = reader.GetOutput()

dim = data.GetDimensions()
spc = data.GetSpacing()
org = data.GetOrigin()

ptData = data.GetPointData()
scalers = ptData.GetScalars()

scalersArray = np_support.vtk_to_numpy(scalers)
scalersArray = scalersArray.reshape(250,250)
minV = scalersArray.min()
maxV = scalersArray.max()

def LERP(p1,p2,C):
    v1 = scalersArray[p1]
    v2 = scalersArray[p2]
    if v1 == v2 == C :
        return p1
    
    r = (C - v1)/(v2 - v1)

    y = p1[0] + (p2[0] - p1[0])*r
    x = p1[1] + (p2[1] - p1[1])*r

    return (x,y,0.0)

# spiral traversal and marching square
def isoLine(C):
    isoLineData = vtkPolyData()
    isoLinePoints = vtkPoints()
    isoLineLines = vtkCellArray()


    T = 0 
    B = dim[0] - 2 
    L = 0 
    R = dim[1] - 2

    def Insertion(p0,p1,p2,p3) :
        cellPoints = [] # to hold cell spiral connection points
        # every edge LERP if needed
        if (C-scalersArray[p0])*(C-scalersArray[p1]) <= 0:
            p = LERP(p0,p1,C)
            cellPoints.append(p)
        if (C-scalersArray[p1])*(C-scalersArray[p2]) <= 0:
            p = LERP(p1,p2,C)
            cellPoints.append(p)
        if (C-scalersArray[p2])*(C-scalersArray[p3]) <= 0:
            p = LERP(p2,p3,C)
            cellPoints.append(p)
        if (C-scalersArray[p3])*(C-scalersArray[p0]) <= 0:
            p = LERP(p3,p0,C)
            cellPoints.append(p)
        
        if cellPoints:
            n = len(cellPoints)
            startId = isoLinePoints.GetNumberOfPoints()

            for p in cellPoints :
                isoLinePoints.InsertNextPoint(p)

            
            polyLine = vtkPolyLine()

            polyLine.GetPointIds().SetNumberOfIds(n + 1)

            for j in range(n) :
                id = startId + j
                polyLine.GetPointIds().SetId(j , id)
            
            polyLine.GetPointIds().SetId(n , startId)

            isoLineLines.InsertNextCell(polyLine)

    while T <= B and L <= R :
        for i in range(L,R+1) :
            p0 = (T,i)
            p1 = (T,i+1)
            p2 = (T+1,i+1)
            p3 = (T+1,i)
            Insertion(p0,p1,p2,p3)
        T += 1
        for i in range(T,B+1) :
            p0 = (i,R)
            p1 = (i,R+1)
            p2 = (i+1,R+1)
            p3 = (i+1,R)
            Insertion(p0,p1,p2,p3)
        R -= 1
        if T <= B :
            for i in range(R,L-1,-1) :
                p0 = (B,i)
                p1 = (B,i+1)
                p2 = (B+1,i+1)
                p3 = (B+1,i)
                Insertion(p0,p1,p2,p3)
            B -= 1
        if L <= R :
            for i in range(B,T-1,-1) :
                p0 = (i,L)
                p1 = (i,L+1)
                p2 = (i+1,L+1)
                p3 = (i+1,L)
                Insertion(p0,p1,p2,p3)
            L += 1
    
    isoLineData.SetPoints(isoLinePoints)
    isoLineData.SetLines(isoLineLines)

    return isoLineData
    




if __name__ == "__main__":
   
    C = float(input("Enter the isovalue:\t"))
    while C < minV or C > maxV :
        print("Invalid! Out of Range!")
        C = float(input("Enter the isovalue:\t"))

    outputData = isoLine(C)

    writer = vtkXMLPolyDataWriter()
    writer.SetFileName(f"Isoline_{C}.vtp")
    writer.SetInputData(outputData)

    writer.Write()
