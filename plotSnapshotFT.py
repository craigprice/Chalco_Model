#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
from decimal import *
from os import listdir
import numpy
from mpl_toolkits.mplot3d import axes3d

import numpy as np
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

PI = 3.14159265358979323846;


BZ = '2'

if BZ == '1':
    BVECTOR_1X = (4 * PI / math.sqrt(3.0)) * (math.sqrt(3.0) / 2.0);
    BVECTOR_1Y = (4 * PI / math.sqrt(3.0)) * ((-1.0) / 2.0);
    BVECTOR_2X = (4 * PI / math.sqrt(3.0)) * 0;
    BVECTOR_2Y = (4 * PI / math.sqrt(3.0)) * 1.0;

if BZ == '2':
    BVECTOR_1X = (4 * PI * math.sqrt(3.0)) * (1.0 / math.sqrt(3.0));
    BVECTOR_1Y = (4 * PI * math.sqrt(3.0)) * 0;
    BVECTOR_2X = (4 * PI * math.sqrt(3.0)) * ((-1.0) / (2.0 * math.sqrt(3.0) ));
    BVECTOR_2Y = (4 * PI * math.sqrt(3.0)) * (1.0 / 2.0);

if BZ == '1':
    NN_kx = [
    BVECTOR_1X,
    (-1)*BVECTOR_2X,
    (-1)*(BVECTOR_2X + BVECTOR_1X),
    (-1)*BVECTOR_1X,
    BVECTOR_2X,
    BVECTOR_1X + BVECTOR_2X
    ]
    
    NN_ky = [
    BVECTOR_1Y,
    (-1)*BVECTOR_2Y,
    (-1)*(BVECTOR_2Y + BVECTOR_1Y),
    (-1)*BVECTOR_1Y,
    BVECTOR_2Y,
    BVECTOR_1Y + BVECTOR_2Y
    ]
    
if BZ == '2':
    NN_kx = [
    BVECTOR_1X,
    (-1)*BVECTOR_2X,
    (-1)*(BVECTOR_2X + BVECTOR_1X),
    (-1)*BVECTOR_1X,
    BVECTOR_2X,
    BVECTOR_1X + BVECTOR_2X
    ]
    
    NN_ky = [
    BVECTOR_1Y,
    (-1)*BVECTOR_2Y,
    (-1)*(BVECTOR_2Y + BVECTOR_1Y),
    (-1)*BVECTOR_1Y,
    BVECTOR_2Y,
    BVECTOR_1Y + BVECTOR_2Y
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
    (float("-inf")),
    (NN_ky[1]/2.0) - slopePerp_NN[1]*(NN_kx[1]/2.0),
    (NN_ky[2]/2.0) - slopePerp_NN[2]*(NN_kx[2]/2.0),
    (float("-inf")),
    (NN_ky[4]/2.0) - slopePerp_NN[4]*(NN_kx[4]/2.0),
    (NN_ky[5]/2.0) - slopePerp_NN[5]*(NN_kx[5]/2.0)
    ]

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
            if isInside1BZ2DHoneycomb(p.x[i_k+1][j_k+1], p.y[i_k+1][j_k+1]):
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
#pathToSimFiles ='/Volumes/Toshi/dataFromSimulations/honV9/hep2/dataFiles'
#pathToSimFiles ='/Volumes/Toshi/dataFromSimulations/duplicateFiles'
pathToSimFiles ='/Users/Work/Documents/perkins_research/honV9K1K2-2BZ/dataFiles_2015_01_06'
#pathToSimFiles ='/Users/Work/Documents/perkins_research/honV9K1K2-2BZ'
simFiles = [ f for f in listdir(pathToSimFiles) if f.find('.txt') > 0 ]

font = {'family' : 'monospace'}
mpl.rc('font', **font)

#For each Sim File extract and make a plot of the FT
f_count = 0
count = 0
for f in simFiles:
    print f
    f_count = f_count + 1
    #outfile = open (f,'w')
    #outfile.write(f)
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
    if f.find('K1_0.707') != -1 and f.find('K2_0.707') != -1: graphTitle += 'Phi = Pi/4'
    if f.find('K1_0') != -1 and f.find('K2_1') != -1: graphTitle += 'Phi = Pi/2'
    if f.find('K1_-0.707') != -1 and f.find('K2_0.707') != -1: graphTitle += 'Phi = 3Pi/4'
    if f.find('K1_-1') != -1 and f.find('K2_0') != -1: graphTitle += 'Phi = Pi'
    if f.find('K1_-0.707') != -1 and f.find('K2_-0.707') != -1: graphTitle += 'Phi = 5Pi/4'
    if f.find('K1_0') != -1 and f.find('K2_-1') != -1: graphTitle += 'Phi = 3Pi/2'
    if f.find('K1_0.707') != -1 and f.find('K2_-0.707') != -1: graphTitle += 'Phi = 7*Pi/4'
    
    if not(f.find('K1_0.707') != -1 and f.find('K1_0.707') != -1 and f.find('KbT_0.11') != -1):
        pass
        #continue

    if graphDescription.find('KbT') != -1: graphTitle = graphTitle + ', T' + graphDescription[graphDescription.find('KbT')+3:]
    #print graphDescription

    print "Open File to get basic characteristics"
    sim_file = open(str(pathToSimFiles + '/' + f),'r')
    numConfigsDone = 0
    cellsA = 0
    cellsB = 0
    cellsC = 0
    KbT = 0
    for line in sim_file:
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
    sim_file.close()
    if cellsA < 80:
        pass
        #continue
    if KbT > 0.048:
        pass
        #continue
    
    if KbT < 0.046:
        pass
        #continue
    
    #numFTSumsDone = numConfigsDone // FTTOUPDATES #integer division
    recipLatt = np.zeros((cellsA,cellsB,2,cellsC))
    #print str(numConfigsDone) + " " + str(cellsA) + " " + str(cellsB) + " " + str(cellsC) + " " + str(numFTSumsDone)
    
    
    #cellsA = cellsA/3#edit
    #cellsB = cellsB/3#edit
    avector_1x = 1*1;
    avector_2x = 1*0.5;
    avector_1y = 1*0;
    avector_2y = 1*math.sqrt(3.0) / 2.0; 
    delta_x = 0
    delta_y = 1.0/math.sqrt(3.0)
    
    xSpinDirection_x = ((-1)*avector_2x) + delta_x
    xSpinDirection_y = ((-1)*avector_2y) + delta_y
    ySpinDirection_x = ((-1)*avector_2x) + avector_1x + delta_x
    ySpinDirection_y = ((-1)*avector_2y) + avector_1y + delta_y
    zSpinDirection_x = delta_x
    zSpinDirection_y = delta_y
    
    rvector_x = []
    rvector_y = []
    spin_x_u = []
    spin_x_v = []
    spin_y_u = []
    spin_y_v = []
    spin_z_u = []
    spin_z_v = []
    x2 = []
    y2 = []
    u2 = []
    v2 = []
    #print str(numConfigsDone) + " " + str(cellsA) + " " + str(cellsB) + " " + str(cellsC) + " " + str(numFTSumsDone)

    arrow_origin_x = []
    arrow_origin_y = []
    arrow_origin_z = []
    arrow_pointing_x = []
    arrow_pointing_y = []
    arrow_pointing_z = []
    arrow_size = []
    max_size = 0
    transparent_arrow_origin_x = []
    transparent_arrow_origin_y = []
    transparent_arrow_origin_z = []
    transparent_arrow_pointing_x = []
    transparent_arrow_pointing_y = []
    transparent_arrow_pointing_z = []
    transparent_arrow_size = []
    
    '''
    #w = 16
    #h = 10
    w = 10
    h = 10
    fig = plt.figure(figsize=(w,h))
    #plt.axis((0,w/10.0,0,h/10.0))
    
    #plt.figure()
    #plt.axis()
    plt.title(graphTitle)
    #plt.figtext(0.81, 0.15, graphDescription)
    #plt.axis((1.2,1.6,0.6,1.0))#edit
    
    ax = fig.gca(projection='3d')  
    '''
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    #ax.view_init(30,-60)
    ax.view_init(30,-30)
    
    max_size = 0
    min_size = 1e10
    max_x = 0
    max_y = 0
    min_x = 1e10
    min_y = 1e10
    print "Extract FT data"
    sim_file = open(str(pathToSimFiles + '/' + f),'r')
    for line in sim_file:
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
        if xPos%1 != 0:
            continue
        if yPos%1 != 0:
            continue
        if abs(magnitude) < 1e-9: continue
        spin_x = math.sqrt(ReSq_x * ReSq_x + ImSq_x * ImSq_x)
        spin_y = math.sqrt(ReSq_y * ReSq_y + ImSq_y * ImSq_y)
        spin_z = math.sqrt(ReSq_z * ReSq_z + ImSq_z * ImSq_z)
        if (spin_x*spin_x + spin_y*spin_y + spin_z*spin_z) > max_size:
            max_size = (spin_x*spin_x + spin_y*spin_y + spin_z*spin_z)
        if (spin_x*spin_x + spin_y*spin_y + spin_z*spin_z) < min_size:
            min_size = (spin_x*spin_x + spin_y*spin_y + spin_z*spin_z)
            
        p = pointIn1BZ(xPos,cellsA,yPos,cellsB,1,1, 0)
        r_x = 0
        r_y = 0
        for i in range (0,3,1):
            for j in range (0,3):
                if p.isExists[i][j]:
                    if p.x[i][j] > max_x:
                        max_x = p.x[i][j]
                    if p.x[i][j] < min_x:
                        min_x = p.x[i][j]
                    if p.y[i][j] > max_y:
                        max_y = p.y[i][j]
                    if p.y[i][j] < min_y:
                        min_y = p.y[i][j]
    sim_file.close()

    max_size = math.sqrt(max_size)
    print "Extract FT data"
    sim_file = open(str(pathToSimFiles + '/' + f),'r')
    for line in sim_file:
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
        if xPos%1 != 0:
            continue
        if yPos%1 != 0:
            continue
        if abs(magnitude) < 1e-9: continue
        spin_x = math.sqrt(ReSq_x * ReSq_x + ImSq_x * ImSq_x)
        spin_y = math.sqrt(ReSq_y * ReSq_y + ImSq_y * ImSq_y)
        spin_z = math.sqrt(ReSq_z * ReSq_z + ImSq_z * ImSq_z)

        #if ((xPos * avector_1x + yPos * avector_2x + sPos * delta_x) / cellsA) < 1.2:#edit
        #    continue;
        #if ((xPos * avector_1y + yPos * avector_2y + sPos * delta_y) / cellsB) < 0.6:#edit
        #    continue;
        '''
        if ((xPos * avector_1x + yPos * avector_2x + sPos * delta_x) / cellsA) > 0.7:#edit
            continue;
        if ((xPos * avector_1x + yPos * avector_2x + sPos * delta_x) / cellsA) < 0.3:#edit
            continue;
        if ((xPos * avector_1y + yPos * avector_2y + sPos * delta_y) / cellsB) > 0.7:#edit
            continue;
        if ((xPos * avector_1y + yPos * avector_2y + sPos * delta_y) / cellsB) < 0.3:#edit
            continue;
        '''
        '''
        s_xu = spin_xu * xSpinDirection_x
        s_xv = spin_xu * xSpinDirection_y
        s_yu = spin_yu * xSpinDirection_x
        s_yv = spin_yu * xSpinDirection_y
        s_zu = spin_zu * xSpinDirection_x
        s_zv = spin_zu * xSpinDirection_y
        size = sqrt(s_x*s_x + s_y*s_y + s_z*s_z)
        '''
        #max_size = 1.0
        size = 1
        #print xPos, yPos, zPos, sPos, spin_x, spin_y, spin_z
        p = pointIn1BZ(xPos,cellsA,yPos,cellsB,1,1, 0)
        r_x = 0
        r_y = 0
        count = count + 1
        if count > 1000 *1000*1000: break
        for i in range (0,3,1):
            for j in range (0,3):
                if p.isExists[i][j]:
                    r_x = p.x[i][j] #(p.x[i][j] * BVECTOR_1X + p.y[i][j] * BVECTOR_2X )# / cellsA
                    r_y = p.y[i][j] #(p.x[i][j] * BVECTOR_1Y + p.y[i][j] * BVECTOR_2Y) #/ cellsB
                    
                    spin_x_u.append(spin_x * xSpinDirection_x / size)
                    spin_x_v.append(spin_x * xSpinDirection_y / size)
                    spin_y_u.append(spin_y * ySpinDirection_x / size)
                    spin_y_v.append(spin_y * ySpinDirection_y / size)
                    spin_z_u.append(spin_z * zSpinDirection_x / size)
                    spin_z_v.append(spin_z * zSpinDirection_y / size)
                    rvector_x.append(r_x)
                    rvector_y.append(r_y)
                    
                    spin_x = (1)*2*PI*(spin_x)/max_size
                    spin_y = (1)*2*PI*(spin_y)/max_size
                    spin_z = (1)*2*PI*spin_z/max_size
                    size = math.sqrt(spin_x*spin_x + spin_y*spin_y + spin_z*spin_z)
                    #alpha_arrow = (-1)/math.log(size/max_size)
                    alpha_arrow = 10*size/max_size
                    #print alpha_arrow
                    
                    spin_x_u_ = spin_x * xSpinDirection_x
                    spin_x_v_ = spin_x * xSpinDirection_y
                    spin_y_u_ = spin_y * ySpinDirection_x
                    spin_y_v_ = spin_y * ySpinDirection_y
                    spin_z_u_ = spin_z * zSpinDirection_x
                    spin_z_v_ = spin_z * zSpinDirection_y
                    
                    spin_x_size = math.sqrt(spin_x_u_*spin_x_u_ + spin_x_v_*spin_x_v_)
                    spin_y_size = math.sqrt(spin_y_u_*spin_y_u_ + spin_y_v_*spin_y_v_)
                    spin_z_size = math.sqrt(spin_z_u_*spin_z_u_ + spin_z_v_*spin_z_v_)
                    
                    spin_color = ''
                    if spin_x_size > spin_y_size and spin_x_size > spin_z_size:
                        spin_color = 'red'
                    if spin_y_size > spin_x_size and spin_y_size > spin_z_size:
                        spin_color = 'blue'
                    if spin_z_size > spin_x_size and spin_z_size > spin_y_size:
                        spin_color = 'green'
                    
                    ar = mpl.patches.ArrowStyle.CurveFilledB(head_length=.1, head_width=0.03)
                    if alpha_arrow > 1:
                        alpha_arrow = 1
                    if alpha_arrow < 3e-2:
                        alpha_arrow = 3e-2
                        #ar = mpl.patches.ArrowStyle.CurveFilledB(head_length=.1, head_width=.03)

                    spin_x = 0     
                    spin_y = 0   
                    spin_z = size                      
                        
                    a = Arrow3D([r_x, r_x+spin_x], [r_y, r_y+spin_y], [0, spin_z], 
                                mutation_scale=200*alpha_arrow,
                                lw=2,
                                arrowstyle=ar,
                                alpha = alpha_arrow,
                                color=spin_color)
                    ax.add_artist(a)
                        
                    '''
                    ax.quiver(r_x+(1)*spin_x,
                              r_y+(1)*spin_y,
                              0+(1)*spin_z,
                              r_x+spin_x,
                              r_y+spin_y,
                              0+spin_z,
                              alpha=1, length=size)
                    '''
                    
                    #print r_x, r_y, spin_x, spin_y, spin_z
                    #outfile.write(str(r_x)+" "+ str(r_y)+" "+ str(spin_x)+" "+ str(spin_y)+" "+ str(spin_z)+"\n")

    #cout<<"df1 "<<n<< " "<<m<<endl;
    #Add Duplicate Points Together Before Graphing

    sim_file.close()
    print "Plot"
    

    
    '''
    Q = ax.quiver(transparent_arrow_origin_x,
                  transparent_arrow_origin_y,
                  transparent_arrow_origin_z,
                  transparent_arrow_pointing_x,
                  transparent_arrow_pointing_y,
                  transparent_arrow_pointing_z, alpha=0.05, norm=transparent_arrow_size)
    Q = ax.quiver(arrow_origin_x,
                  arrow_origin_y,
                  arrow_origin_z,
                  arrow_pointing_x,
                  arrow_pointing_y,
                  arrow_pointing_z, alpha=1, norm=arrow_size)
    '''
    
    ax.auto_scale_xyz([min_x, max_x], [min_y, max_y], [0, 6])
    
    '''
    Q = plt.quiver(rvector_x,rvector_y,spin_x_u,spin_x_v,color='red',pivot='middle')#,headaxislength=1.1,headlength=1.1,scale=60,width=w/10000.0,headwidth=2)
    Q = plt.quiver(rvector_x,rvector_y,spin_y_u,spin_y_v,color='blue',pivot='middle')#,headaxislength=1.1,headlength=1.1,scale=60,width=w/10000.0,headwidth=2)
    Q = plt.quiver(rvector_x,rvector_y,spin_z_u,spin_z_v,color='green',pivot='middle')#,headaxislength=1.1,headlength=1.1,scale=60,width=w/10000.0,headwidth=2)
    #Q = plt.quiver(rvector_x,rvector_y,spin_x_u,spin_x_v,color='red',pivot='middle',headaxislength=1.1,headlength=1.1,scale=8,width=w/1000.0,headwidth=2)#edit
    #Q = plt.quiver(rvector_x,rvector_y,spin_y_u,spin_y_v,color='blue',pivot='middle',headaxislength=1.1,headlength=1.1,scale=8,width=w/1000.0,headwidth=2)#edit
    #Q = plt.quiver(rvector_x,rvector_y,spin_z_u,spin_z_v,color='green',pivot='middle',headaxislength=1.1,headlength=1.1,scale=8,width=w/1000.0,headwidth=2)#edit
    '''
    '''
    plt.tick_params(\
    axis='both',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    right='off',         # ticks along the top edge are off
    left='off',         # ticks along the top edge are off
    labelbottom='off',
    labelleft='off')#edit  
    '''
    #BZ='1'
    ex = (4 * PI / math.sqrt(3.0)) * 0.001
    side_by_2 = ((4 * PI / math.sqrt(3.0)) * 1.0 / 2.0 ) * math.tan(30*PI/180.0) 
    max_x = ((4 * PI / math.sqrt(3.0)) * 1.0 / 2.0 ) / math.cos(30*PI/180.0)
    min_x = -1*max_x
    max_y = ((4 * PI / math.sqrt(3.0)) * 1.0 / 2.0 )
    min_y = -1*max_y
    line1_BZ1=[(max_x, 0), (side_by_2,min_y)]
    line2_BZ1=[(line1_BZ1[1][0], line1_BZ1[1][1]), (-1*side_by_2,min_y)]
    line3_BZ1=[(line2_BZ1[1][0], line2_BZ1[1][1]), (min_x,0)]
    line4_BZ1=[(line3_BZ1[1][0], line3_BZ1[1][1]), (-1*side_by_2,max_y)]
    line5_BZ1=[(line4_BZ1[1][0], line4_BZ1[1][1]), (side_by_2,max_y)]
    line6_BZ1=[(line5_BZ1[1][0], line5_BZ1[1][1]), (max_x, 0)]
    
    a = Arrow3D([line1_BZ1[0][0], line1_BZ1[1][0]], [line1_BZ1[0][1], line1_BZ1[1][1]], [0, 0], 
                lw=2,
                arrowstyle=mpl.patches.ArrowStyle.Curve(),
                color='black')
    ax.add_artist(a)
    a = Arrow3D([line2_BZ1[0][0], line2_BZ1[1][0]], [line2_BZ1[0][1], line2_BZ1[1][1]], [0, 0], 
                lw=2,
                arrowstyle=mpl.patches.ArrowStyle.Curve(),
                color='black')
    ax.add_artist(a)
    a = Arrow3D([line3_BZ1[0][0], line3_BZ1[1][0]], [line3_BZ1[0][1], line3_BZ1[1][1]], [0, 0], 
                lw=2,
                arrowstyle=mpl.patches.ArrowStyle.Curve(),
                color='black')
    ax.add_artist(a)
    a = Arrow3D([line4_BZ1[0][0], line4_BZ1[1][0]], [line4_BZ1[0][1], line4_BZ1[1][1]], [0, 0], 
                lw=2,
                arrowstyle=mpl.patches.ArrowStyle.Curve(),
                color='black')
    ax.add_artist(a)
    a = Arrow3D([line5_BZ1[0][0], line5_BZ1[1][0]], [line5_BZ1[0][1], line5_BZ1[1][1]], [0, 0], 
                lw=2,
                arrowstyle=mpl.patches.ArrowStyle.Curve(),
                color='black')
    ax.add_artist(a)
    a = Arrow3D([line6_BZ1[0][0], line6_BZ1[1][0]], [line6_BZ1[0][1], line6_BZ1[1][1]], [0, 0], 
                lw=2,
                arrowstyle=mpl.patches.ArrowStyle.Curve(),
                color='black')
    ax.add_artist(a)
    
    
    if BZ == '2':
        oversize = 1.025
        max_x = 6.28318530718 * oversize
        max_y = 7.25519745694 * oversize
        min_x = -6.28318530718 * oversize
        min_y = -7.25519745694 * oversize
        side_by_2 = ((4 * PI * math.sqrt(3.0)) * (1.0 / math.sqrt(3.0)) / 2.0 ) * math.tan(30*PI/180.0)
        side_by_2 = side_by_2 * oversize
        line1_BZ2=[(max_x,side_by_2), (max_x,-1*side_by_2)]
        line2_BZ2=[(line1_BZ2[1][0], line1_BZ2[1][1]), (0,min_y)]
        line3_BZ2=[(line2_BZ2[1][0], line2_BZ2[1][1]), (min_x,-1*side_by_2)]
        line4_BZ2=[(line3_BZ2[1][0], line3_BZ2[1][1]), (min_x,side_by_2)]
        line5_BZ2=[(line4_BZ2[1][0], line4_BZ2[1][1]), (0,max_y)]
        line6_BZ2=[(line5_BZ2[1][0], line5_BZ2[1][1]), (max_x,side_by_2)]
        a = Arrow3D([line1_BZ2[0][0], line1_BZ2[1][0]], [line1_BZ2[0][1], line1_BZ2[1][1]], [0, 0], 
                lw=2,
                arrowstyle=mpl.patches.ArrowStyle.Curve(),
                color='black')
        ax.add_artist(a)
        a = Arrow3D([line2_BZ2[0][0], line2_BZ2[1][0]], [line2_BZ2[0][1], line2_BZ2[1][1]], [0, 0], 
                    lw=2,
                    arrowstyle=mpl.patches.ArrowStyle.Curve(),
                    color='black')
        ax.add_artist(a)
        a = Arrow3D([line3_BZ2[0][0], line3_BZ2[1][0]], [line3_BZ2[0][1], line3_BZ2[1][1]], [0, 0], 
                    lw=2,
                    arrowstyle=mpl.patches.ArrowStyle.Curve(),
                    color='black')
        ax.add_artist(a)
        a = Arrow3D([line4_BZ2[0][0], line4_BZ2[1][0]], [line4_BZ2[0][1], line4_BZ2[1][1]], [0, 0], 
                    lw=2,
                    arrowstyle=mpl.patches.ArrowStyle.Curve(),
                    color='black')
        ax.add_artist(a)
        a = Arrow3D([line5_BZ2[0][0], line5_BZ2[1][0]], [line5_BZ2[0][1], line5_BZ2[1][1]], [0, 0], 
                    lw=2,
                    arrowstyle=mpl.patches.ArrowStyle.Curve(),
                    color='black')
        ax.add_artist(a)
        a = Arrow3D([line6_BZ2[0][0], line6_BZ2[1][0]], [line6_BZ2[0][1], line6_BZ2[1][1]], [0, 0], 
                    lw=2,
                    arrowstyle=mpl.patches.ArrowStyle.Curve(),
                    color='black')
        ax.add_artist(a)
    #lc = mpl.collections.LineCollection([line1,line2,line3,line4,line5,line6],colors='black', linewidths=1)
    #ax.add_collection(lc) 
    plt.title(graphTitle)
    #plt.figtext(0.88, 0.4, graphDescription + str('\nS(k) is Oriented:\nRed = X\nBlue = Y\nGreen = Z'))
    #plt.savefig('figuresFT/'+ f[:-4] + '.eps',dpi=150)
    plt.savefig('figuresFT/'+ graphTitle.replace('/','_').strip() + '.png',dpi=100)
    plt.show()
    #plt.close('all')
    #outfile.close()
    if f_count > 1:
        pass
        #break
    
print "\a"
