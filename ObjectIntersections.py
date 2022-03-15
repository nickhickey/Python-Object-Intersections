#Author: Nicholas Hickey

import maya.cmds as cmds
import math
ZV = 0.000000000000000000001

#Delete the previous interface window and replace with new one only if first is open
if('MyWin' == globals()):
    if cmds.window(MyWin, exists=True):
        cmds.deleteUI(MyWin, window=True)

#Setting height and width of the user interface
MyWin = cmds.window(title='IMD3002_Assignment01', menuBar=True,widthHeight=(500,400))

#creating a find intersection button in the user interface which calls the findIntersect() fucntion which is defined later
cmds.columnLayout( columnAttach = ('left', 5), rowSpacing=10, columnWidth=500)
cmds.button( label='Find the Intersection', command='findIntersect()')
cmds.setParent("..")
cmds.paneLayout()
cmds.textScrollList('uiPointList', numberOfRows=50, allowMultiSelection=False)
cmds.showWindow(MyWin)

#matrix multiplier to convert to a homogenous system with the w value = 1
def matrixMultiplier(Matrix, Point):
   
    PointOut = [0.0, 0.0, 0.0, 0.0]
    PointIn = [Point[0], Point[1], Point[2], 1] 
   
    PointOut[0] = (Matrix[0]*PointIn[0])+(Matrix[4]*PointIn[1])+(Matrix[8]*PointIn[2])+(Matrix[12]*PointIn[3])
    PointOut[1] = (Matrix[1]*PointIn[0])+(Matrix[5]*PointIn[1])+(Matrix[9]*PointIn[2])+(Matrix[13]*PointIn[3])
    PointOut[2] = (Matrix[2]*PointIn[0])+(Matrix[6]*PointIn[1])+(Matrix[10]*PointIn[2])+(Matrix[14]*PointIn[3])
    PointOut[3] = (Matrix[3]*PointIn[0])+(Matrix[7]*PointIn[1])+(Matrix[11]*PointIn[2])+(Matrix[15]*PointIn[3])
    
    return(PointOut)
    
#define get magnitude of point to intersection to the sphere vertex to solve for distance
def getMagnitude(sphereVtx, intVtx):
   
    magnitude = 0.0
    magnitude = (((sphereVtx[0] - intVtx[0])**2) +  ((sphereVtx[1] - intVtx[1])**2) +  ((sphereVtx[2] - intVtx[2])**2)) ** 0.5  
    return(magnitude)
    
#get the normal vector of the three cube vertex making up the facet
def getNormalVector(vtxA, vtxB, vtxC):
   
    crossProduct = [0.0, 0.0, 0.0]
    normVec = [0.0, 0.0, 0.0]
    vec1 = [0.0, 0.0, 0.0]
    vec2 = [0.0, 0.0, 0.0]
    magnitude = 0.0
   
    vec1[0] = (vtxA[0] - vtxB[0])
    vec1[1] = (vtxA[1] - vtxB[1])
    vec1[2] = (vtxA[2] - vtxB[2])
    
    vec2[0] = (vtxB[0] - vtxC[0])
    vec2[1] = (vtxB[1] - vtxC[1])
    vec2[2] = (vtxB[2] - vtxC[2])
    
    
    crossProduct[0] = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1])         
    crossProduct[1] = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2])
    crossProduct[2] = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0])
    
    magnitude = (((crossProduct[0])**2) +  ((crossProduct[1])**2) +  ((crossProduct[2])**2))**0.5
    
    normVec[0] = (crossProduct[0] / magnitude)
    normVec[1] = (crossProduct[1] / magnitude)
    normVec[2] = (crossProduct[2] / magnitude)
    
    return(normVec)

#Cross product fucntion to be used to complete cross product on two vectors
def crossProduct(vecA, vecB, vecC):
    vec1[0] = (vtxA[0] - vtxB[0])
    vec1[1] = (vtxA[1] - vtxB[1])
    vec1[2] = (vtxA[2] - vtxB[2])
    
    vec2[0] = (vtxB[0] - vtxC[0])
    vec2[1] = (vtxB[1] - vtxC[1])
    vec2[2] = (vtxB[2] - vtxC[2])
    cross = [0.0, 0.0, 0.0]
    cross[0] = (vec1[1] * vec2[2]) - (vec1[2] * vec2[1])         
    cross[1] = (vec1[2] * vec2[0]) - (vec1[0] * vec2[2])
    cross[2] = (vec1[0] * vec2[1]) - (vec1[1] * vec2[0])
    
    return(cross)
    
    
#Area of intersecting facet: magnitude of (axb) is also equal to the area of the parallelogram formed
def areaOfFacet(vecA, vecB, vecC):
    magnitude = 0.0
    area = 0.0
    crosspro = [0,0,0]

    vecA1 = [(vecB[0]-vecA[0]) , (vecB[1]-vecA[1]) , (vecB[2]-vecA[2])]
    vecA2 = [(vecC[0]-vecA[0]) , (vecC[1]-vecA[1]) , (vecC[2]-vecA[2])]
    
    crosspro = [((vecA1[1] * vecA2[2]) - (vecA1[2] * vecA2[1])),((vecA1[2] * vecA2[0]) - (vecA1[0] * vecA2[2])),((vecA1[0] * vecA2[1]) - (vecA1[1] * vecA2[0]))]
    

    area = ( ((crosspro[0])**2) + ((crosspro[1])**2) + ((crosspro[2])**2) )**0.5
    
    return(area)
    
#gets the equation of the plane using the normal   
def getNormalToPlane(normal, point):
    planeEq = [0.0, 0.0, 0.0, 0.0]
    
    d = -((normal[0] * point[0]) + (normal[1] * point[1]) +(normal[2] * point[2])) 
    
    planeEq = [normal[0], normal[1], normal[2], d]
    
    return planeEq
    
#dot product function to be used to find dot product
def dotProduct(vtx1, vtx2):
    dot = 0.0
    
    dot = ((vtx1[0] * vtx2[0]) + (vtx1[1] * vtx2[1]) + (vtx1[2] * vtx2[2]))
 
    return dot

#finding the angle to the xy grid 
def angleToGrid(PointA, PointB):
    
    angleGrid = 0.0
    dotPointABtoGrid = 0.0
    dotB = 0.0
    magPointAB = 0.0
    magPointAG = 0.0
    vecPointAB = [0,0,0]
    vecPointAG = [0,0,0]
    radToDeg = 180/3.14159265359
    
    vecPointAB = [(PointB[0]-PointA[0]) , (PointB[1]-PointA[1]) , (PointB[2]-PointA[2])]
    vecPointAG = [(PointB[0]-PointA[0]) , (PointB[1]) , (PointB[2]-PointA[2])]
     
    magPointAB = (((vecPointAB[0])**2) +  ((vecPointAB[1])**2) +  ((vecPointAB[2])**2)) ** 0.5
    magPointAG = (((vecPointAG[0])**2) +  ((vecPointAG[1])**2) +  ((vecPointAG[2])**2)) ** 0.5
  
    dotPointABtoGrid = ((vecPointAB[0] * vecPointAG[0]) + (vecPointAB[1] * vecPointAG[1]) + (vecPointAB[2] * vecPointAG[2]))
    angleGrid = math.acos(dotPointABtoGrid/(magPointAB*magPointAG))
    angleGrid = angleGrid * radToDeg
    
    #to always get the smaller angle the lines are making with the XY grid
    if (angleGrid > 90):
        angleGrid = 180 - angleGrid
 
    return angleGrid  

#check intersection fucntion to see where the point lies in space, ie to remove unwanted points and leave the point on the face remaining by comparing dot products  
def intersectionCheck(vtxA, vtxB, vtxC, vtxD, intPoint):
    
    #intiialize vectors
    vecAB = [0,0,0]
    vecAC = [0,0,0]
    vecAI = [0,0,0]
    vecDA = [0,0,0]
    vecDB = [0,0,0]
    vecDI = [0,0,0] 
    vecBA = [0,0,0]
    vecBC = [0,0,0]
    vecBI = [0,0,0]
    vecCB = [0,0,0]
    vecCA = [0,0,0]
    vecCI = [0,0,0]
    vecAD = [0,0,0]
    
    #initialize magnitudes
    magAB = 0.0
    magAC = 0.0
    magAI =  0.0
    magDA = 0.0
    magDB = 0.0
    magDI = 0.0
    magBA = 0.0
    magBC =  0.0
    magBI = 0.0
    magCB = 0.0
    magCA =  0.0
    magCI =  0.0
    
    #create the vectors 
    vecAB = [(vtxB[0]-vtxA[0]) , (vtxB[1]-vtxA[1]) , (vtxB[2]-vtxA[2])]
    vecAC = [(vtxC[0]-vtxA[0]) , (vtxC[1]-vtxA[1]) , (vtxC[2]-vtxA[2])]
    vecAI = [(intPoint[0]-vtxA[0]) , (intPoint[1]-vtxA[1]) , (intPoint[2]-vtxA[2])] 

    vecBD = [(vtxD[0]-vtxB[0]) , (vtxD[1]-vtxB[1]) , (vtxD[2]-vtxB[2])]
     
    vecDB = [(vtxB[0]-vtxD[0]) , (vtxB[1]-vtxD[1]) , (vtxB[2]-vtxD[2])]
    vecDC = [(vtxC[0]-vtxD[0]) , (vtxC[1]-vtxD[1]) , (vtxC[2]-vtxD[2])]
    
    vecCD = [(vtxD[0]-vtxC[0]) , (vtxD[1]-vtxC[1]) , (vtxD[2]-vtxC[2])]
    vecDI = [(intPoint[0]-vtxD[0]) , (intPoint[1]-vtxD[1]) , (intPoint[2]-vtxD[2])]
    
    vecBA = [(vtxD[0]-vtxB[0]) , (vtxD[1]-vtxB[1]) , (vtxD[2]-vtxB[2])]
    vecBC = [(vtxC[0]-vtxB[0]) , (vtxC[1]-vtxB[1]) , (vtxC[2]-vtxB[2])]
    vecBI = [(intPoint[0]-vtxB[0]) , (intPoint[1]-vtxB[1]) , (intPoint[2]-vtxB[2])]

    vecCB = [(vtxB[0]-vtxC[0]) , (vtxB[1]-vtxC[1]) , (vtxB[2]-vtxC[2])]
    vecCA = [(vtxA[0]-vtxC[0]) , (vtxA[1]-vtxC[1]) , (vtxA[2]-vtxC[2])]
    vecCI = [(intPoint[0]-vtxC[0]) , (intPoint[1]-vtxC[1]) , (intPoint[2]-vtxC[2])]   

   #get the respective magnitudes
    magAB = (((vecAB[0])**2) +  ((vecAB[1])**2) +  ((vecAB[2])**2))**0.5 
    magAC = (((vecAC[0])**2) +  ((vecAC[1])**2) +  ((vecAC[2])**2))**0.5 
    magAI = (((vecAI[0])**2) +  ((vecAI[1])**2) +  ((vecAI[2])**2))**0.5 
     
    magDB = (((vecDB[0])**2) +  ((vecDB[1])**2) +  ((vecDB[2])**2))**0.5 
    magDC = (((vecDC[0])**2) +  ((vecDC[1])**2) +  ((vecDC[2])**2))**0.5 
    magDI = (((vecDI[0])**2) +  ((vecDI[1])**2) +  ((vecDI[2])**2))**0.5 
    
    
    magBD = (((vecBD[0])**2) +  ((vecBD[1])**2) +  ((vecBD[2])**2))**0.5 
    magBA = (((vecBA[0])**2) +  ((vecBA[1])**2) +  ((vecBA[2])**2))**0.5 
    magBC = (((vecBC[0])**2) +  ((vecBC[1])**2) +  ((vecBC[2])**2))**0.5 
    magBI = (((vecBI[0])**2) +  ((vecBI[1])**2) +  ((vecBI[2])**2))**0.5    
    
    magCB = (((vecCB[0])**2) +  ((vecCB[1])**2) +  ((vecCB[2])**2))**0.5 
    magCA = (((vecCA[0])**2) +  ((vecCA[1])**2) +  ((vecCA[2])**2))**0.5 
    magCI = (((vecCI[0])**2) +  ((vecCI[1])**2) +  ((vecCI[2])**2))**0.5 
    magCD = (((vecCD[0])**2) +  ((vecCD[1])**2) +  ((vecCD[2])**2))**0.5 


    #intialize all the dot products
    dotA  = 0.0
    dotB  = 0.0
    dotC  = 0.0
    dotD  = 0.0
    dotD2 = 0.0
    dotAI = 0.0
    dotBI = 0.0
    dotCI = 0.0
    dotDI = 0.0
    dotCD = 0.0
 
    #get the dot product value
    dotA = ((vecAB[0] * vecAC[0]) + (vecAB[1] * vecAC[1]) + (vecAB[2] * vecAC[2]))
    dotAI = ((vecAB[0] * vecAI[0]) + (vecAB[1] * vecAI[1]) + (vecAB[2] * vecAI[2]))
    dotB = ((vecBA[0] * vecBC[0]) + (vecBA[1] * vecBC[1]) + (vecBA[2] * vecBC[2]))
    dotBI = ((vecBA[0] * vecBI[0]) + (vecBA[1] * vecBI[1]) + (vecBA[2] * vecBI[2]))
    dotC = ((vecCB[0] * vecCA[0]) + (vecCB[1] * vecCA[1]) + (vecCB[2] * vecCA[2]))  
    dotCI = ((vecCB[0] * vecCI[0]) + (vecCB[1] * vecCI[1]) + (vecCB[2] * vecCI[2])) 
    dotD = ((vecDB[0] * vecDC[0]) + (vecDB[1] * vecDC[1]) + (vecDB[2] * vecDC[2]))
    dotDI = ((vecDB[0] * vecDI[0]) + (vecDB[1] * vecDI[1]) + (vecDB[2] * vecDI[2]))
    dotBD = ((vecBC[0] * vecBD[0]) + (vecBC[1] * vecBD[1]) + (vecBC[2] * vecBD[2]))
    dotCD = ((vecCD[0] * vecCB[0]) + (vecCD[1] * vecCB[1]) + (vecCD[2] * vecCB[2]))
    
    #divide the dot products by the respective sum of the multiplying magnitudes to normalise the values
    dotA = dotA / (magAB * magAC)
    dotB = dotB / (magBA * magBC)
    dotC = dotC / (magCB * magCA)
    dotD = dotD / (magDB * magDC)
    dotCD = dotCD / (magCD * magCA)
    dotBD = dotBD / (magBC * magBD)
    dotAI = dotAI / (magAB * magAI)
    dotBI = dotBI / (magBA * magBI)
    dotCI = dotCI / (magCB * magCI)
    dotDI = dotDI / (magDB * magDI)  
     
    
    #only return true if the statement suffices all the dot quantities. ie if dotAI value is larger than dotA value essentially closer to 1.
    if (( (dotAI > dotA) and (dotBI > dotB) and (dotCI > dotC)) or ((dotDI > dotD) and (dotBI > dotBD) and (dotCI > dotCD))):
        return True

#equation used to find the t value to be used in the intersection. Returns a value from 0 to 1 
def getTValue(planeEq, PointA, PointB):
    denominator = 0.0
    numerator = 0.0
    
    denominator=(planeEq[0]*(PointA[0]-PointB[0]))+(planeEq[1]*(PointA[1]-PointB[1]))+(planeEq[2]*(PointA[2]-PointB[2]))
    
    if(abs(denominator) < ZV):
        print "the demoniator is zero"
        return False
    
    numerator = (planeEq[0] * PointA[0]) + (planeEq[1] * PointA[1]) + (planeEq[2] * PointA[2]) + planeEq[3]
    
    return(numerator/denominator)


#defining the findIntersect fucntion which is called in the interface
def findIntersect():

    selectedShapes = cmds.ls(selection=True)
    shapeCount = 0
    meshList = []
    intersectCount = 0
    distanceCount  = 0
    areaCount      = 0
    areaQuadCount  = 0 
    angleCount     = 0
    normalCount    = 0
    
    #Which shape is selected first becomes the          
    for shape in selectedShapes:
        if(cmds.objectType(shape) == 'transform'):
            childShape = cmds.listRelatives(shape, fullPath=True, shapes=True)
            shapeCount += 1     
            
        if(cmds.objectType(childShape) == 'mesh'):
            meshList.append(childShape)
             
    #if the user has selected anything under two shapes they will get a print out telling them to select two         
    if(shapeCount < 2):
        print "Please select two shapes"
        return False

    #set the variables to the respective selected shapes
    referenceObj = selectedShapes[0]
    targetObj = selectedShapes[1]
    
    #translation points of the referenece object
    translationPoint = cmds.xform(referenceObj, query=True, translation=True, worldSpace=True)
    meshXFormT = cmds.xform(targetObj, query=True, matrix=True, worldSpace=True)
    meshXFormR = cmds.xform(referenceObj, query=True, matrix=True, worldSpace=True)
   
    facetCount = cmds.polyEvaluate(targetObj, face=True)
    vertCount = cmds.polyEvaluate(referenceObj, vertex=True)
    
    #loop through all the vertex in the range of the evaluated vertex count
    for vrt in range(0, vertCount):

        vtxA1 = cmds.getAttr(referenceObj + ".vt[" + str(vrt) + "]")       
        vtxNewA1 = matrixMultiplier(meshXFormR, list(vtxA1[0]))
        
        #Pointa and Pointb values are the speheres exterior vertex and the center of the sphere
        PointA = vtxNewA1
        PointB = translationPoint
        cmds.curve(d = 1, p=[(vtxNewA1[0],vtxNewA1[1],vtxNewA1[2]),(translationPoint[0], translationPoint[1], translationPoint[2])] )        

        #loops through the faces in the facetcount range of the target object ie cube in this case
        for face in range(0, facetCount):
            
            vtxLst = cmds.polyInfo((targetObj + ".f[" + str(face) + "]"), faceToVertex=True)
            vtxIdx = str(vtxLst[0]).split()
            
            vtxA = cmds.getAttr(targetObj + ".vt[" + vtxIdx[2] + "]")
            vtxB = cmds.getAttr(targetObj + ".vt[" + vtxIdx[3] + "]")
            vtxC = cmds.getAttr(targetObj + ".vt[" + vtxIdx[4] + "]") 
            vtxD = cmds.getAttr(targetObj + ".vt[" + vtxIdx[5] + "]")                                   
                                    
            vtxNewA = matrixMultiplier(meshXFormT, list(vtxA[0]))
            vtxNewB = matrixMultiplier(meshXFormT, list(vtxB[0]))
            vtxNewC = matrixMultiplier(meshXFormT, list(vtxC[0])) 
            vtxNewD = matrixMultiplier(meshXFormT, list(vtxD[0]))  
            
            #call the normal equation using the spheres vertexs
            normal  = getNormalVector(vtxNewA, vtxNewB, vtxNewC)      
            #use the noraml in the get normal to plane equation 
            planeEq = getNormalToPlane(normal, vtxNewA)

            sValueA = (planeEq[0]*PointA[0])+(planeEq[1]*PointA[1])+(planeEq[2]*PointA[2])+planeEq[3]
            sValueB = (planeEq[0]*PointB[0])+(planeEq[1]*PointB[1])+(planeEq[2]*PointB[2])+planeEq[3]
                
            if(((sValueA>0.0) and (sValueB<0.0)) or ((sValueA<0.0) and (sValueB>0.0))):
                tValue = getTValue(planeEq, PointA, PointB)
                PointI = [0.0, 0.0, 0.0]
                PointI[0] = PointA[0] + (tValue * (PointB[0] - PointA[0]))
                PointI[1] = PointA[1] + (tValue * (PointB[1] - PointA[1]))
                PointI[2] = PointA[2] + (tValue * (PointB[2] - PointA[2]))
                
                #if the check of the intersection is true it then can draw the cubes at the intersection points and print the respective data 
                check = intersectionCheck(vtxNewA, vtxNewB, vtxNewC, vtxNewD, PointI)
                if (check == True):
                        
                    cube = cmds.polyCube(h=0.1, w = 0.1, d = 0.1)
                    cmds.move(PointI[0], PointI[1], PointI[2], cube)
                    
                    #Text outputs for the user interface window
                    sphereToInt = getMagnitude(PointA, PointI)
                    distanceText = "Distance of Intersection [" + str(distanceCount) + "]: " + str(round(sphereToInt,2))    
                    cmds.textScrollList('uiPointList', edit=True, append=[distanceText])
                    distanceCount += 1 

                    facetArea = areaOfFacet(vtxNewA, vtxNewB, vtxNewC)
                    areaText = "Area Of Triangle [" + str(areaCount) + "]: " + str(round(facetArea,2))
                    cmds.textScrollList('uiPointList', edit=True, append=[areaText])
                    areaCount += 1
                    
                    facetArea2 = (areaOfFacet(vtxNewA, vtxNewB, vtxNewC)) * 2
                    areaText2 = "Area Of Quad [" + str(areaQuadCount) + "]: " + str(round(facetArea2,2))
                    cmds.textScrollList('uiPointList', edit=True, append=[areaText2])
                    areaQuadCount += 1

                    angle = angleToGrid(PointA, PointB)
                    angleText = "Angle [" + str(angleCount) + "] " + str(round(angle,2))
                    cmds.textScrollList('uiPointList', edit=True, append=[angleText])
                    angleCount += 1
                    
                    ptText = "Intersection: [" + str(intersectCount) + "] " + str(round(PointI[0],2)) + ", " + str(round(PointI[1],2)) + ", " + str(round(PointI[2],2))
                    cmds.textScrollList('uiPointList', edit=True, append=[ptText])
                    intersectCount += 1
                    
                    normalVec = "Normal Vector: [" + str(normalCount) + "] " + str(round(normal[0],2)) + ", " + str(round(normal[1],2)) + ", " + str(round(normal[2],2))
                    cmds.textScrollList('uiPointList', edit=True, append=[normalVec])
                    normalCount += 1
                    
                    
