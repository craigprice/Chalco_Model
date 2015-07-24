#!/usr/bin/python


import numpy
import matplotlib.pyplot as plt
import math
import sys
import re
from decimal import *
import os
import matplotlib as mpl
from scipy.optimize import curve_fit

#font = {'family' : 'monospace'}
#mpl.rc('font', **font)
mpl.rc('font',**{'serif':['Helvetica']})
#mpl.rc('text', usetex=True)
mpl.rc('text.latex', preamble=r'\usepackage{cmbright}')

pathToGraphFiles ='readyToPlot/'
#graphFiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(pathToGraphFiles)) for f in fn  if f.find('.txt') > 0]
graphFiles = [ f for f in os.listdir(pathToGraphFiles) if f.find('.txt') > 0 ]
#print graphFiles

for f in graphFiles:
    opName = f[f.find('OP_')+3:-4]
    pieces = f.split('_')
    graphDescription = ''
    
    for p in range(1, len(pieces)-3, 2):
        tempStr = pieces[p+1]
        tempStr = tempStr.rjust(7 - len(pieces[p])) + '\n'
        graphDescription = graphDescription + pieces[p] + ' = ' + tempStr
    if pieces[len(pieces)-1] != "M.txt":
        pass
        #continue
    #myfile = open(f,'r')
    myfile = open('readyToPlot/'+f,'r')
    graphTitle = pieces[len(pieces)-1][:-4]
    data = {}#{kbt,{cellsA, value}}
    cells = []
    print f
    for line in myfile:
        line = line.strip()
        pieces = line.split()
        data.update({pieces[0]: {}})
        for p in range(1, len(pieces), 2):
            data[pieces[0]].update({pieces[p]: pieces[p+1]})
            if pieces[p] not in cells:
                cells.append(pieces[p])
    myfile.close()
    fig = plt.figure()
    curve = []
    xaxis = []
    #print data
    #print cells
    color=iter(plt.cm.rainbow(numpy.linspace(0,1,len(cells))))
    for c in range(len(cells)):
        for kbt in data.keys():
            try:
                #print kbt
                #print data[kbt]
                #print data[kbt][cells[c]]
                curve.append(float(data[kbt][cells[c]]))
                xaxis.append(float(kbt))
            except KeyError:
                print "missing cells:", cells[c], "kbt:", kbt
                
        c=next(color)
        plt.plot(xaxis, curve, color=c, marker='o', linestyle='None')
        xaxis = []
        curve = []
    ax = plt.axes()
    #ax.set_aspect(1./ax.get_data_ratio())
    #plt.xlabel('B')
    font_size = 24
    plt.xlabel('T',fontsize=font_size)
    ax.yaxis.labelpad = 15
    if opName == "SpecificHeat":
        pass
        plt.ylabel('C',rotation='horizontal', fontsize=font_size)
    else:
        pass
        #plt.ylabel('E', rotation='horizontal', fontsize=font_size)
            
        
    #plt.ylabel(opName)
    plt.minorticks_on
    #plt.grid(True, which='both')
    plt.title(opName,y=1.03,fontsize=font_size)
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    #plt.subplots_adjust(right = 0.8)
    #plt.ylim(min(yaxis)-0.1,max(yaxis)+0.1)
    if opName[:6] == "Binder":
        #plt.ylim(0,0.67)
        pass
    if opName[:13] == "Magnetization":
        #plt.ylim(0,1)
        pass
    plt.figtext(0.81, 0.15, graphDescription)
    #plt.savefig('figures/'+ f[f.find('/'):-4] + '.png')
    #plt.savefig('figures/'+ f[:-4] + '.png')
    #plt.savefig('figures/'+ graphTitle.replace('/','_') + '.png', dpi = 100)
    #plt.savefig('figures/eps/'+ graphTitle.replace('/','_') + '.eps', dpi = 150)
    plt.savefig('figures/'+ f[:-4].strip().replace('\\','') + '.pdf',dpi=100)
    plt.close()
    #exit()

print "\a"
