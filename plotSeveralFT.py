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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *
import scipy.interpolate
import matplotlib as mpl

PI = 3.14159265358979323846;
FTTOUPDATES = 4000;
CORRTOUPDATES = 20;

BVECTOR_1X = (1/1.0)*(4 * PI / sqrt(3.0)) * sqrt(3.0) / 2.0;
BVECTOR_2X = (1/1.0)*0;
BVECTOR_1Y = (1/1.0)*(4 * PI / sqrt(3.0)) * (-1) / 2.0;
BVECTOR_2Y = (1/1.0)*(4 * PI / sqrt(3.0)) * 1;

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

slopePerp_NN = [
(-1.0)/(NN_ky[0]/NN_kx[0]),
0,
(-1.0)/(NN_ky[2]/NN_kx[2]),
(-1.0)/(NN_ky[3]/NN_kx[3]),
0,
(-1.0)/(NN_ky[5]/NN_kx[5])
]

bcoor_NN = [
(NN_ky[0]/2.0) - slopePerp_NN[0]*(NN_kx[0]/2.0),
(NN_ky[1]/2.0) - slopePerp_NN[1]*(NN_kx[1]/2.0),
(NN_ky[2]/2.0) - slopePerp_NN[2]*(NN_kx[2]/2.0),
(NN_ky[3]/2.0) - slopePerp_NN[3]*(NN_kx[3]/2.0),
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



#Start Main Function
#Get list of all sim files
pathToSimFiles ='/Users/Work/Documents/perkinsResearch/honV9/figuresFT/sampleFT/'
simFiles = [ f for f in listdir(pathToSimFiles) if f.find('.txt') > 0 ]

#font = {'family' : 'monospace'}
#mpl.rc('font', **font)

w = 8*6
h = 0.8*w   
fig = plt.figure(figsize=(w,h))
count = 1

#For each Sim File extract and make a plot of the FT
for f in simFiles:
    if count >= 3:
        pass
        #continue
    print f
    opName = f[f.find('OP_')+3:-4]
    pieces = f.split('_')
    graphDescription = ''
    
    #Extract sim parameters
    numSpaces = 7
    for p in range(1, len(pieces)-3, 2):
        tempStr = pieces[p+1]
        tempStr = tempStr.rjust(numSpaces - len(pieces[p])) + '\n'
        graphDescription = graphDescription + pieces[p] + ' = ' + tempStr
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
        continue
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
    

    indexOfPointsInBZ = -1
    for n4 in range(0, cellsA):
        for m4 in range(0, cellsB):
            for k4 in range(0, 2-1):
                #cout<<n4<< " "<<m4<<" "<<k4<<endl;
                indexOfPointsInBZ = indexOfPointsInBZ + 1
                allPointsInBZArray[indexOfPointsInBZ] = pointIn1BZ(n4,cellsA,m4,cellsB,k4,2,recipLatt[n4][m4][k4][0] / numFTSumsDone)
        

#cout<<"df1 "<<n<< " "<<m<<endl;
#Add Duplicate Points Together Before Graphing

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
    for p in singlePointInBZArr:
        x.append(p.x)
        y.append(p.y)
        z.append(p.magSk)
        
    scale1 = 1.06
    x_point = max(x)*scale1
    x_flat = (2*math.pi/3.0)*scale1
    y_top = max(y)*scale1
    y_bottom = min(y)*scale1

    line1=[(-1*x_flat, y_bottom),(-1*x_point,0)]
    line2=[(-1*x_point,0),       (-1*x_flat,y_top)]
    line3=[(-1*x_flat,y_top),    (x_flat,y_top)]
    line4=[(x_flat,y_top),       (x_point,0)]
    line5=[(x_point,0),          (x_flat,y_bottom)]
    line6=[(x_flat,y_bottom),    (-1*x_flat, y_bottom)]

    line1Fit = np.polyfit( [line1[0][0],line1[1][0]], [line1[0][1],line1[1][1]], 1)
    line1x = arange(line1[1][0],line1[0][0],0.1)
    line1y = line1Fit[1] + line1Fit[0] * line1x        
    
    line2Fit = np.polyfit( [line2[0][0],line2[1][0]], [line2[0][1],line2[1][1]], 1)
    line2x = arange(line2[0][0],line2[1][0],0.1)
    line2y = line2Fit[1] + line2Fit[0] * line2x  
    
    line3Fit = np.polyfit( [line3[0][0],line3[1][0]], [line3[0][1],line3[1][1]], 1)
    line3x = arange(line3[0][0],line3[1][0],0.1)
    line3y = line3Fit[1] + line3Fit[0] * line3x  

    line4Fit = np.polyfit( [line4[0][0],line4[1][0]], [line4[0][1],line4[1][1]], 1)
    line4x = arange(line4[0][0],line4[1][0],0.01)
    line4y = line4Fit[1] + line4Fit[0] * line4x  

    line5Fit = np.polyfit( [line5[0][0],line5[1][0]], [line5[0][1],line5[1][1]], 1)
    line5x = arange(line5[1][0],line5[0][0],0.1)
    line5y = line5Fit[1] + line5Fit[0] * line5x  

    line6Fit = np.polyfit( [line6[0][0],line6[1][0]], [line6[0][1],line6[1][1]], 1)
    line6x = arange(line6[1][0],line6[0][0],0.1)
    line6y = line6Fit[1] + line6Fit[0] * line6x      
    
    x.extend(line1x)
    y.extend(line1y)
    zz = [1e-6 for i in line1x]
    z.extend(zz)
    
    x.extend(line2x)
    y.extend(line2y)
    zz = [1e-6 for i in line2x]
    z.extend(zz)
    
    x.extend(line3x)
    y.extend(line3y)
    zz = [1e-6 for i in line3x]
    z.extend(zz)
    
    x.extend(line4x)
    y.extend(line4y)
    zz = [1e-6 for i in line4x]
    z.extend(zz)
    
    x.extend(line5x)
    y.extend(line5y)
    zz = [1e-6 for i in line5x]
    z.extend(zz)
    
    x.extend(line6x)
    y.extend(line6y)  
    zz = [1e-6 for i in line6x]
    z.extend(zz)
    
    xi, yi = np.linspace(min(x)*1.1, max(x)*1.1, 100), np.linspace(min(y)*1.1, max(y)*1.1, 100)
    xi, yi = np.meshgrid(xi, yi)
    
    print "Plot"

    ax = fig.add_subplot(1,6,count)
    subplots_adjust(wspace=0)
    ax.set_aspect(1./ax.get_data_ratio())
    plt.scatter(x, y, c=z,edgecolors='none',vmin=min(z), vmax=max(z),s=20)#90)


    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)
    plt.imshow(zi, vmin=min(z)*1.1, vmax=max(z)*1.1, origin='lower',extent=[min(x)*1.1, max(x)*1.1, min(y)*1.1, max(y)*1.1])#,cmap=plt.cm.winter)    
    
    plt.xlim(min(x)*1.05, max(x)*1.05)
    plt.ylim(min(y)*1.05, max(y)*1.05)
    if count == 6:
        pass
        #divider = make_axes_locatable(ax)
        #cax = divider.append_axes("right", size="5%", pad=0.05)
        #cbar = plt.colorbar(cax=cax)
        #cbar.set_ticks( [] )
    #ax = plt.gca()
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)

    #plt.figtext(0.81, 0.15, graphDescription)
   

    
    
    lc = mpl.collections.LineCollection([line1,line2,line3,line4,line5,line6],colors='white', linewidths=4)
    ax.add_collection(lc)  
    count = count + 1
        
textheight = 0.405
textsize = 6*10
plt.figtext(0.126,0.548,'(c)', fontsize=textsize, color='white')
plt.figtext(0.176, textheight, 'FM',fontsize=textsize)
plt.figtext(0.277, textheight, '3Q-Spiral',fontsize=textsize)
plt.figtext(0.424, textheight, 'Stripy',fontsize=textsize)
plt.figtext(0.521, textheight, 'Intermediate',fontsize=textsize)
plt.figtext(0.657, textheight, '120$^\circ$ State',fontsize=textsize)
plt.figtext(0.805, textheight, 'Zigzag',fontsize=textsize)
#plt.savefig('figuresFT/fourier.png',dpi=100)
plt.savefig('figuresFT/fourier.eps',dpi=40)
plt.clf()

print "\a"
