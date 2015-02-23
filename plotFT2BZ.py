#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
import re
from decimal import *
from os import listdir
from os.path import isfile, join
from scipy.optimize import curve_fit
from pylab import *
import scipy.interpolate
import matplotlib as mpl

PI = 3.14159265358979323846;
FTTOUPDATES = 4000;
CORRTOUPDATES = 20;
BZ = '1'

if BZ == '1':
    BVECTOR_1X = (4 * PI / sqrt(3.0)) * sqrt(3.0) / 2.0;
    BVECTOR_2X = 0;
    BVECTOR_1Y = (4 * PI / sqrt(3.0)) * (-1) / 2.0;
    BVECTOR_2Y = (4 * PI / sqrt(3.0)) * 1;

if BZ == '2':
    BVECTOR_1X = (4 * PI * sqrt(3.0)) * (1.0 / sqrt(3.0));
    BVECTOR_1Y = 0;
    BVECTOR_2X = (4 * PI * sqrt(3.0)) * -1 / (2.0 * sqrt(3.0) );
    BVECTOR_2Y = (4 * PI * sqrt(3.0)) * (1.0 / 2.0);

NN_ky = [
BVECTOR_1Y,
(-1)*BVECTOR_2Y,
(-1)*(BVECTOR_2Y + BVECTOR_1Y),
(-1)*BVECTOR_1Y,
BVECTOR_2Y,
BVECTOR_1Y + BVECTOR_2Y
]


NN_kx = [
BVECTOR_1X,
(-1)*BVECTOR_2X,
(-1)*(BVECTOR_2X + BVECTOR_1X),
(-1)*BVECTOR_1X,
BVECTOR_2X,
BVECTOR_1X + BVECTOR_2X
]

if BZ == '1':
    slopePerp_NN = [
    (-1.0)/(NN_ky[0]/NN_kx[0]),
    0,
    (-1.0)/(NN_ky[2]/NN_kx[2]),
    (-1.0)/(NN_ky[3]/NN_kx[3]),
    0,
    (-1.0)/(NN_ky[5]/NN_kx[5])
    ]

if BZ == '2':
    slopePerp_NN = [
    (float("inf")),
    (-1.0)/(NN_ky[1]/NN_kx[1]),
    (-1.0)/(NN_ky[2]/NN_kx[2]),
    (float("-inf")),
    (-1.0)/(NN_ky[4]/NN_kx[4]),
    (-1.0)/(NN_ky[5]/NN_kx[5])
    ]
    
if BZ == '1':
    bcoor_NN = [
    (NN_ky[0]/2.0) - slopePerp_NN[0]*(NN_kx[0]/2.0),
    (NN_ky[1]/2.0) - slopePerp_NN[1]*(NN_kx[1]/2.0),
    (NN_ky[2]/2.0) - slopePerp_NN[2]*(NN_kx[2]/2.0),
    (NN_ky[3]/2.0) - slopePerp_NN[3]*(NN_kx[3]/2.0),
    (NN_ky[4]/2.0) - slopePerp_NN[4]*(NN_kx[4]/2.0),
    (NN_ky[5]/2.0) - slopePerp_NN[5]*(NN_kx[5]/2.0)
    ]


if BZ == '2':
    bcoor_NN = [
    (float("inf")),
    (NN_ky[1]/2.0) - slopePerp_NN[1]*(NN_kx[1]/2.0),
    (NN_ky[2]/2.0) - slopePerp_NN[2]*(NN_kx[2]/2.0),
    (float("-inf")),
    (NN_ky[4]/2.0) - slopePerp_NN[4]*(NN_kx[4]/2.0),
    (NN_ky[5]/2.0) - slopePerp_NN[5]*(NN_kx[5]/2.0)
    ]

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class allPointsInBZ:
    def __init__(self):
        self.x = [[0 for i in range(0,3)] for j in range(0,3)]
        self.y = [[0 for i in range(0,3)] for j in range(0,3)]
        self.magSk = [[0 for i in range(0,3)] for j in range(0,3)]
        self.isExists = [[False for i in range(0,3)] for j in range(0,3)]

class singlePointInBZ:
    def __init__(self):
        self.x = 0
        self.y = 0
        self.magSk = 0


def isInside1BZ2DHoneycomb(px_kx, py_ky):
    #cout<<"py_ky "<<py_ky<<
    #" NN1 "<<(slopePerp_NN_1*px_kx + bcoor_NN_1)/PI<<
    #" NN2 "<<(slopePerp_NN_2*px_kx + bcoor_NN_2)/PI<<
    #" NN3 "<<(slopePerp_NN_3*px_kx + bcoor_NN_3)/PI<<
    #" NN4 "<<(slopePerp_NN_4*px_kx + bcoor_NN_4)/PI<<
    #" NN5 "<<(slopePerp_NN_5*px_kx + bcoor_NN_5)/PI<<
    #" NN6 "<<(slopePerp_NN_6*px_kx + bcoor_NN_6)/PI<<endl;
    isAllTrue = True
    #print px_kx, py_ky
    for i in range (0,3):
        if py_ky > ((slopePerp_NN[i]*px_kx + bcoor_NN[i])-1e-7):
            pass #inside BZ
        else:
            isAllTrue = False
    for i in range (3,6):
        if py_ky < ((slopePerp_NN[i]*px_kx + bcoor_NN[i])+1e-7):
            pass #inside BZ
        else:
            isAllTrue = False
    #print isAllTrue
    return isAllTrue





def pointIn1BZ( numPointsInRecSpaceX,  linSizeX,
    numPointsInRecSpaceY,  linSizeY,
    numPointsInRecSpaceSublattice,  numSubLattices,
    magSk):
    p = allPointsInBZ()
    for i in range (0,3):
        for j in range (0,3):
            p.x[i][j] = 0
            p.y[i][j] = 0
            p.magSk[i][j] = 0
            p.isExists[i][j] = False
    px_kx = ((numPointsInRecSpaceX/(linSizeX*1.0)) * BVECTOR_1X)
    px_kx += ((numPointsInRecSpaceY/(linSizeY*1.0)) * BVECTOR_2X)
    px_kx += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_1X)
    px_kx += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_2X)
    py_ky =  ((numPointsInRecSpaceX/(linSizeX*1.0)) * BVECTOR_1Y)
    py_ky += ((numPointsInRecSpaceY/(linSizeY*1.0)) * BVECTOR_2Y)
    py_ky += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_1Y)
    py_ky += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_2Y)
    for i_k in range (-1,2,1):
        for j_k in range (-1,2,1):
            p.x[i_k+1][j_k+1] = px_kx + (i_k * BVECTOR_1X + j_k * BVECTOR_2X)
            p.y[i_k+1][j_k+1] = py_ky + (i_k * BVECTOR_1Y + j_k * BVECTOR_2Y)
            #cout<<"numPointsInRecSpaceX "<<numPointsInRecSpaceX<<" numPointsInRecSpaceY "<<numPointsInRecSpaceY<<" p.x "<<p.x/PI<<" p.y "<<p.y/PI<<endl;
            if isInside1BZ2DHoneycomb(p.x[i_k+1][j_k+1], p.y[i_k+1][j_k+1]):
                p.isExists[i_k+1][j_k+1]=True
                p.magSk[i_k+1][j_k+1]=magSk
    for i in range (0,3,1):
        for j in range (0,3):
            if p.isExists[i][j]:
                return p
    print "can't find BZ for", numPointsInRecSpaceX, numPointsInRecSpaceY, numPointsInRecSpaceSublattice, "setting to zero"
    #p.isExists[0][0] = 0;
    #p.x[0][0] = 0;
    #p.y[0][0] = 0;
    #p.magSk[0][0] = 0;
    exit(1)
    return p





def isInside2BZ2DHoneycomb(px_kx, py_ky):
    #print "px_kx ", px_kx, "py_ky ", py_ky
    '''print "slopePerp_NN[0]", slopePerp_NN[0], "slopePerp_NN[1]", slopePerp_NN[1]
    print "slopePerp_NN[2]", slopePerp_NN[2], "slopePerp_NN[3]", slopePerp_NN[3]
    print "slopePerp_NN[4]", slopePerp_NN[4], "slopePerp_NN[5]", slopePerp_NN[5]
    print "bcoor_NN[0]", bcoor_NN[0], "bcoor_NN[1]", bcoor_NN[1]
    print "bcoor_NN[2]", bcoor_NN[2], "bcoor_NN[3]", bcoor_NN[3]
    print "bcoor_NN[4]", bcoor_NN[4], "bcoor_NN[5]", bcoor_NN[5]
    print " NN1 ",(slopePerp_NN[0]*px_kx + bcoor_NN[0])/math.pi
    print " NN2 ",(slopePerp_NN[1]*px_kx + bcoor_NN[1])/math.pi
    print " NN3 ",(slopePerp_NN[2]*px_kx + bcoor_NN[2])/math.pi
    print " NN4 ",(slopePerp_NN[3]*px_kx + bcoor_NN[3])/math.pi
    print " NN5 ",(slopePerp_NN[4]*px_kx + bcoor_NN[4])/math.pi
    print " NN6 ",(slopePerp_NN[5]*px_kx + bcoor_NN[5])/math.pi'''
    isAllTrue = True
    #print px_kx, py_ky
    #print "((BVECTOR_1X/2.0)-1e-7)", ((BVECTOR_1X/2.0)-1e-7)
    if px_kx < ((BVECTOR_1X/2.0)+1e-7):
        pass #inside BZ
    else:
        isAllTrue = False
    for i in range (1,3):
        #print "((slopePerp_NN[i]*px_kx + bcoor_NN[i])-1e-7)", ((slopePerp_NN[i]*px_kx + bcoor_NN[i])-1e-7)
        if py_ky > ((slopePerp_NN[i]*px_kx + bcoor_NN[i])-1e-7):
            pass #inside BZ
        else:
            isAllTrue = False
    for i in range (4,6):
        #print "((slopePerp_NN[i]*px_kx + bcoor_NN[i])+1e-7)", ((slopePerp_NN[i]*px_kx + bcoor_NN[i])+1e-7)
        if py_ky < ((slopePerp_NN[i]*px_kx + bcoor_NN[i])+1e-7):
            pass #inside BZ
        else:
            isAllTrue = False
    #print "((BVECTOR_1X/2.0)-1e-7)", ((-1.0*BVECTOR_1X/2.0)+1e-7)
    if px_kx > ((-1.0*BVECTOR_1X/2.0)-1e-7):
        pass #inside BZ
    else:
        isAllTrue = False
    #print isAllTrue
    return isAllTrue



def pointIn2BZ( numPointsInRecSpaceX,  linSizeX,
    numPointsInRecSpaceY,  linSizeY,
    numPointsInRecSpaceSublattice,  numSubLattices,
    magSk):
    p = allPointsInBZ()
    for i in range (0,3):
        for j in range (0,3):
            p.x[i][j] = 0
            p.y[i][j] = 0
            p.magSk[i][j] = 0
            p.isExists[i][j] = False
    px_kx = ((numPointsInRecSpaceX/(linSizeX*1.0)) * BVECTOR_1X)
    px_kx += ((numPointsInRecSpaceY/(linSizeY*1.0)) * BVECTOR_2X)
    #px_kx += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_1X)
    #px_kx += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_2X)
    py_ky =  ((numPointsInRecSpaceX/(linSizeX*1.0)) * BVECTOR_1Y)
    py_ky += ((numPointsInRecSpaceY/(linSizeY*1.0)) * BVECTOR_2Y)
    #py_ky += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_1Y)
    #py_ky += ((numPointsInRecSpaceSublattice/(numSubLattices*1.0)) * BVECTOR_2Y)
    for i_k in range (-1,2,1):
        for j_k in range (-1,2,1):
            p.x[i_k+1][j_k+1] = px_kx + (i_k * BVECTOR_1X + j_k * BVECTOR_2X)
            p.y[i_k+1][j_k+1] = py_ky + (i_k * BVECTOR_1Y + j_k * BVECTOR_2Y)
            #cout<<"numPointsInRecSpaceX "<<numPointsInRecSpaceX<<" numPointsInRecSpaceY "<<numPointsInRecSpaceY<<" p.x "<<p.x/PI<<" p.y "<<p.y/PI<<endl;
            if isInside2BZ2DHoneycomb(p.x[i_k+1][j_k+1], p.y[i_k+1][j_k+1]):
                p.isExists[i_k+1][j_k+1]=True
                p.magSk[i_k+1][j_k+1]=magSk
    for i in range (0,3,1):
        for j in range (0,3):
            if p.isExists[i][j]:
                return p
    print "can't find BZ for", numPointsInRecSpaceX, numPointsInRecSpaceY, 0, "setting to zero"
    #p.isExists[0][0] = 0;
    #p.x[0][0] = 0;
    #p.y[0][0] = 0;
    #p.magSk[0][0] = 0;
    exit(1)
    return p



#Start Main Function
#Get list of all sim files
#pathToSimFiles ='/Volumes/Toshi/dataFromSimulations/duplicateFiles'
pathToSimFiles ='/Users/Work/Documents/perkins_research/honV9K1K2-2BZ'
simFiles = [ f for f in listdir(pathToSimFiles) if f.find('.txt') > 0 ]

font = {'family' : 'monospace'}
mpl.rc('font', **font)

#For each Sim File extract and make a plot of the FT
for f in simFiles:
    print f
    opName = f[f.find('OP_')+3:-4]
    pieces = f.split('_')
    graphDescription = ''
    
    #Extract sim parameters
    numSpaces = 7
    graphTitle = ''
    phi = ''
    for p in range(1, len(pieces)-1, 2):
        tempStr = pieces[p+1]
        tempStr = tempStr.rjust(numSpaces - len(pieces[p])) + '\n'
        tempStr = tempStr.replace('.txt','',1)
        graphDescription = graphDescription + pieces[p] + ' = ' + tempStr
    '''
    if graphDescription.find(' 0.93') != -1 and graphDescription.find(' 0.368') != -1: graphTitle += 'Phi = 0.06';continue
    if graphDescription.find(' 0.771') != -1 and graphDescription.find(' 0.637') != -1: graphTitle += 'Phi = 0.11';continue
    if graphDescription.find(' 0.536') != -1 and graphDescription.find(' 0.844') != -1: graphTitle += 'Phi = 0.16';continue
    if graphDescription.find(' 0.249') != -1 and graphDescription.find(' 0.969') != -1: graphTitle += 'Phi = 0.21';continue
    if graphDescription.find(' -0.063') != -1 and graphDescription.find(' 0.998') != -1: graphTitle += 'Phi = 0.26';continue
    if graphDescription.find(' -0.368') != -1 and graphDescription.find(' 0.93') != -1: graphTitle += 'Phi = 0.31';continue
    if graphDescription.find(' -0.637') != -1 and graphDescription.find(' 0.771') != -1: graphTitle += 'Phi = 0.36';continue
    if graphDescription.find(' -0.844') != -1 and graphDescription.find(' 0.536') != -1: graphTitle += 'Phi = 0.41';continue
    if graphDescription.find(' -0.969') != -1 and graphDescription.find(' 0.249') != -1: graphTitle += 'Phi = 0.46';continue
    if graphDescription.find(' -0.998') != -1 and graphDescription.find(' -0.063') != -1: graphTitle += 'Phi = 0.51';continue
    if graphDescription.find(' -0.93') != -1 and graphDescription.find(' -0.368') != -1: graphTitle += 'Phi = 0.56';continue
    if graphDescription.find(' -0.771') != -1 and graphDescription.find(' -0.637') != -1: graphTitle += 'Phi = 0.61';continue
    if graphDescription.find(' -0.536') != -1 and graphDescription.find(' -0.844') != -1: graphTitle += 'Phi = 0.66';continue
    if graphDescription.find(' -0.249') != -1 and graphDescription.find(' -0.969') != -1: graphTitle += 'Phi = 0.71';continue
    if graphDescription.find(' 0.063') != -1 and graphDescription.find(' -0.998') != -1: graphTitle += 'Phi = 0.76';continue
    if graphDescription.find(' 0.368') != -1 and graphDescription.find(' -0.93') != -1: graphTitle += 'Phi = 0.81';continue
    if graphDescription.find(' 0.637') != -1 and graphDescription.find(' -0.771') != -1: graphTitle += 'Phi = 0.86';continue
    if graphDescription.find(' 0.844') != -1 and graphDescription.find(' -0.536') != -1: graphTitle += 'Phi = 0.91';continue
    if graphDescription.find(' 0.969') != -1 and graphDescription.find(' -0.249') != -1: graphTitle += 'Phi = 0.96';continue
    '''
    
    if f.find('K1_1') != -1 and f.find('K2_0') != -1: graphTitle += 'Phi = 0'
    if f.find('K1_0.707') != -1 and f.find('K2_0.707') != -1: graphTitle += 'Phi = 0.125'
    if f.find('K1_0') != -1 and f.find('K2_1') != -1: graphTitle += 'Phi = 0.25'
    if f.find('K1_-0.707') != -1 and f.find('K2_0.707') != -1: graphTitle += 'Phi = 0.375'
    if f.find('K1_-1') != -1 and f.find('K2_0') != -1: graphTitle += 'Phi = 0.5'
    if f.find('K1_-0.707') != -1 and f.find('K2_-0.707') != -1: graphTitle += 'Phi = 0.625'
    if f.find('K1_0') != -1 and f.find('K2_-1') != -1: graphTitle += 'Phi = 0.75'
    if f.find('K1_0.707') != -1 and f.find('K2_-0.707') != -1: graphTitle += 'Phi = 0.875'

    if graphDescription.find('KbT') != -1: graphTitle = graphTitle + ', T' + graphDescription[graphDescription.find('KbT')+3:]
    #print graphDescription

    print "Open File to get basic characteristics"
    file = open(str(pathToSimFiles + '/' + f),'r')
    numConfigsDone = 0
    cellsA = 0
    cellsB = 0
    cellsC = 0
    KbT = 0
    for line in file:
        if line.find('Cells-A: ') == 0:
            cellsA = int(line[len('Cells-A: '):])
        if line.find('KbT: ') == 0:
            KbT = float(line[len('KbT: '):])
        if line.find('Cells-B: ') == 0:
            cellsB = int(line[len('Cells-B: '):])
        if line.find('Cells-C: ') == 0:
            cellsC = int(line[len('Cells-C: '):])
        if line.find('numConfigsDone: ') == 0:
            numConfigsDone = int(line[len('numConfigsDone: '):])
    file.close()
    if cellsA < 80:
        pass
        #continue
    if KbT > 0.01:
        pass
        #continue
    numFTSumsDone = numConfigsDone // FTTOUPDATES #integer division
    recipLatt = np.zeros((cellsA,cellsB,2,cellsC))
    #print str(numConfigsDone) + " " + str(cellsA) + " " + str(cellsB) + " " + str(cellsC) + " " + str(numFTSumsDone)


    print "Extract FT data"
    file = open(str(pathToSimFiles + '/' + f),'r')
    for line in file:
        if line.find('Reciprocal Lattice, S(q): ') < 0:
            continue
        component = line[len('Reciprocal Lattice, S(q): '):]
        #print line[:-1]
        xPos = int(component[:component.find(","):])
        component = component[1 + component.find(","):]
        yPos = int(component[:component.find(",")])
        component = component[1 + component.find(","):]
        zPos = int(component[:component.find(",")])
        component = component[1 + component.find(","):]
        sPos = int(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        ReSq_x = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        ReSq_y = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        ReSq_z = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        ImSq_x = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        ImSq_y = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        ImSq_z = float(component[:component.find(" ")])
        component = component[1 + component.find(" "):]
        magnitude = float(component[:component.find(" ")])
        recipLatt[xPos][yPos][sPos][zPos] = magnitude
        #print xPos, yPos, " ", zPos, " ", sPos, " ", ReSq_x, " ", ReSq_y, " ", ReSq_z, " ", ImSq_x, " ", ImSq_y, " ", ImSq_z, " ", magnitude
    
    
    numPointsInBZ = cellsA * cellsB * 2

    print "initialize"
    allPointsInBZArray = [None for j in range(0,numPointsInBZ)]
    
    print "build array of Kvector points"

    indexOfPointsInBZ = -1
    for n4 in range(0, cellsA):
        for m4 in range(0, cellsB):
            k4 = 0
            #for k4 in range(0, 2-1):
            #print n4,cellsA,m4,cellsB,1,1,recipLatt[n4][m4][0][0]
            indexOfPointsInBZ = indexOfPointsInBZ + 1
            if BZ == '2':
                allPointsInBZArray[indexOfPointsInBZ] = pointIn2BZ(n4,cellsA,m4,cellsB,1,1,recipLatt[n4][m4][0][0] )#/ numFTSumsDone)
            if BZ == '1':
                allPointsInBZArray[indexOfPointsInBZ] = pointIn1BZ(n4,cellsA,m4,cellsB,k4,2,recipLatt[n4][m4][k4][0] / numFTSumsDone)


#cout<<"df1 "<<n<< " "<<m<<endl;
#Add Duplicate Points Together Before Graphing

    print "Flatten points array to a bunch of single points"
    singlePointInBZArr = []

    for i5 in range(0, indexOfPointsInBZ + 1):
        for i6 in range(0, 3):
            for j6 in range(0, 3):
                if allPointsInBZArray[i5].isExists[i6][j6]:
                    p = singlePointInBZ()
                    p.x = allPointsInBZArray[i5].x[i6][j6]
                    p.y = allPointsInBZArray[i5].y[i6][j6]
                    p.magSk = allPointsInBZArray[i5].magSk[i6][j6]
                    singlePointInBZArr.append(p)



#plt.plot(xaxis, yaxis, 'o')
#plt.axes()
#plt.xlabel('bx')
#plt.ylabel('by')
#cout<<"df2 "<<n<< " "<<m<<endl;
    '''
    print "Adding points together"
    while(len(singlePointInBZArr)>0):
        count = count + 1
        p = singlePointInBZArr[-1]
    
        for j5 in range(0, len(singlePointInBZArr)):
            if (fabs(p.x - singlePointInBZArr[j5].x) <= 1e-7) and (fabs(p.y - singlePointInBZArr[j5].y) <= 1e-7):
                p.magSk = p.magSk + singlePointInBZArr[j5].magSk
                del singlePointInBZArr[j5]
                j5 = j5 - 1

#cout<<"df5 "<<i5<< " "<<m<<endl;
#cout<<count<<" "<<p.x<< " "<<p.y<<" " <<p.magSk<<endl;
        x.append(p.x)
        y.append(p.y)
        z.append(p.magSk)
___



    while(len(singlePointInBZArr)>0):
        count = count + 1
        p = singlePointInBZArr[-1]
    
        for j5 in singlePointInBZArr:
            if (fabs(p.x - j5.x) <= 1e-7) and (fabs(p.y - j5.y) <= 1e-7):
                print p.x, p.y, j5.x, j5.y, j5.magSk
                p.magSk = p.magSk + j5.magSk
                del j5

    for i in singlePointInBZArr:
        for j in singlePointInBZArr:
            if (fabs(i.x - j.x) <= 1e-7) and (fabs(i.y - j.y) <= 1e-7) and (not j.isAdded):
                print i.x, i.y, i.magSk, j.x, j.y, j.magSk
                i.magSk = i.magSk + j.magSk
                j.isAdded = True
    '''

    x = []
    y = []
    z = []
    count = -1;
    print "Adding points together"
    for p in singlePointInBZArr:
        x.append(p.x)
        y.append(p.y)
        z.append(p.magSk)
        
    print "Set up a regular grid of interpolation points"
    xi, yi = np.linspace(min(x), max(x), 100), np.linspace(min(y), max(y), 100)
    xi, yi = np.meshgrid(xi, yi)
    
    print "Plot"
    #rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    #zi = rbf(xi, yi)
    #plt.imshow(zi, vmin=min(z), vmax=max(z), origin='lower',extent=[min(x), max(x), min(y), max(y)])

    w = 8
    h = 6    
    plt.figure(figsize=(w,h))
    plt.scatter(x, y, c=z,edgecolors='none',vmin=min(z), vmax=max(z))
    plt.scatter(x, y, c=z,edgecolors='none',vmin=min(z), vmax=max(z),s=20)#90)
    #plt.scatter(x, y, c=z,cmap=plt.cm.hot_r,edgecolors='none',vmin=min(z), vmax=max(z))
    
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)    
    plt.imshow(zi, vmin=min(z)*1.1, vmax=max(z)*1.1, origin='lower',extent=[min(x)*1.1, max(x)*1.1, min(y)*1.1, max(y)*1.1])#,cmap=plt.cm.winter)    
    
    
    plt.xlim(min(x)*1.05, max(x)*1.05)
    plt.ylim(min(y)*1.05, max(y)*1.05)
    plt.subplots_adjust(left=0.125,right = float(h/float(w)))
    ax = plt.gca()
    ax.set_aspect(1./ax.get_data_ratio())
    cbar = plt.colorbar()
    try:
        cbar.set_ticks( (int(0),int(max(z))) )
    except:
        pass
        
    plt.figtext(0.81, 0.15, graphDescription)
    plt.title(graphTitle)
   

    if BZ == '1':
        x_point = max(x)*(1.025)
        x_flat = (2*math.pi/3.0)*(1.025)
        y_top = max(y)*(1.025)
        y_bottom = min(y)*(1.025)
        line1=[(-1*x_flat, y_bottom),(-1*x_point,0)]
        line2=[(-1*x_point,0),       (-1*x_flat,y_top)]
        line3=[(-1*x_flat,y_top),    (x_flat,y_top)]
        line4=[(x_flat,y_top),       (x_point,0)]
        line5=[(x_point,0),          (x_flat,y_bottom)]
        line6=[(x_flat,y_bottom),    (-1*x_flat, y_bottom)]
    
    if BZ == '2':
        side_length_by_2 = max(y)/2.0
        x_point = max(x)*(1.025)
        x_flat = (2*math.pi/3.0)*(1.025)
        y_top = max(y)*(1.025)
        y_bottom = min(y)*(1.025)
        line1=[(max(x)*(1.025),-1*side_length_by_2*(1.025)),(0*(1.025),min(y)*(1.025))]
        line2=[(0*(1.025),min(y)*(1.025)),       (min(x)*(1.025),-1*side_length_by_2*(1.025))]
        line3=[(min(x)*(1.025),-1*side_length_by_2*(1.025)),    (min(x)*(1.025),1*side_length_by_2*(1.025))]
        line4=[(min(x)*(1.025),1*side_length_by_2*(1.025)),       (0*(1.025),max(y)*(1.025))]
        line5=[(0*(1.025),max(y)*(1.025)),          (max(x)*(1.025),1*side_length_by_2*(1.025))]
        line6=[(max(x)*(1.025),1*side_length_by_2*(1.025)),    (max(x)*(1.025),-1*side_length_by_2*(1.025))]
        
    lc = mpl.collections.LineCollection([line1,line2,line3,line4,line5,line6],colors='black', linewidths=1)
    ax.add_collection(lc)    
        
    
    plt.savefig('figuresFT/'+ f[:-4] + '.png',dpi=300)
    #plt.savefig('figuresFT/'+ f[:-4] + '.eps',dpi=100)
    plt.clf()

print "\a"
