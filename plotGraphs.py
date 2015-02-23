#!/usr/bin/python


import numpy as np
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
    graphTitle = ''
    if f.find('K1_1') != -1         and f.find('K2_0') != -1:       graphTitle += r'$\varphi = 0$'
    if f.find('K1_0.707') != -1     and f.find('K2_0.707') != -1:   graphTitle += r'$\varphi = \frac{\Pi}{4}$'
    if f.find('K1_0') != -1         and f.find('K2_1') != -1:       graphTitle += r'$\varphi = \frac{\Pi}{2}$'
    if f.find('K1_-0.707') != -1    and f.find('K2_0.707') != -1:   graphTitle += r'$\varphi = \frac{3\Pi}{4}$'
    if f.find('K1_-1') != -1        and f.find('K2_0') != -1:       graphTitle += r'$\varphi = \Pi$'
    if f.find('K1_-0.707') != -1    and f.find('K2_-0.707') != -1:  graphTitle += r'$\varphi = \frac{5\Pi}{4}$'
    if f.find('K1_0') != -1         and f.find('K2_-1') != -1:      graphTitle += r'$\varphi = \frac{3\Pi}{2}$'
    if f.find('K1_0.707') != -1     and f.find('K2_-0.707') != -1:  graphTitle += r'$\varphi = \frac{7\Pi}{4}$'
    if f.find('SpecificHeat') != -1: graphTitle += ', C'
    if f.find('Energy') != -1: graphTitle += ', E'
    xaxis = []
    yaxis = []
    pair = []
    print f
    for line in myfile:
        xaxis.append(float(line[:line.find(' ')]))
        yaxis.append(float(line[line.find(' '):]))
        pair.append([float(line[:line.find(' ')]),float(line[line.find(' '):])])
    fig = plt.figure()
    def getKey(item):
        return item[0]
    pair = sorted(pair, key=getKey)
    xaxis = []
    yaxis = []
    for p in pair:
        xaxis.append(p[0])
        yaxis.append(p[1])
#ave = (yaxis[0] + yaxis[1] + yaxis[2])/3.0
 #   if yaxis[0] > 0.5:
#      print yaxis[0]
#continue
    plt.plot(xaxis, yaxis, 'r-o')
    ax = plt.axes()
    #ax.set_aspect(1./ax.get_data_ratio())
 #   plt.xlabel('B')
    font_size = 24
    plt.xlabel('T',fontsize=font_size)
    ax.yaxis.labelpad = 15
    if opName == "SpecificHeat":
        plt.ylabel('C',rotation='horizontal', fontsize=font_size)
    else:
        plt.ylabel('E', rotation='horizontal', fontsize=font_size)
            
        
    #plt.ylabel(opName)
    plt.minorticks_on
    #plt.grid(True, which='both')
    plt.title(graphTitle[0:-3],y=1.03,fontsize=font_size)
    ax.tick_params(axis='x', labelsize=font_size)
    ax.tick_params(axis='y', labelsize=font_size)
    #plt.subplots_adjust(right = 0.8)
    plt.ylim(min(yaxis)-0.1,max(yaxis)+0.1)
    if opName[:6] == "Binder":
        plt.ylim(0,0.67)
    if opName[:13] == "Magnetization":
        plt.ylim(0,1)
#popt = []
#pcov = []
#if opName[:6] == "Binder":
#guess = (0.4,0.01,2.0/3.0)#mu, t, n
#plt.plot(xaxis,fermiDirac(np.asarray(xaxis), guess[0], guess[1], guess[2]), 'r-',ls='-')
#popt,pcov = curve_fit(fermiDirac, xaxis, yaxis, guess)
#for p in popt:
#fit = fermiDirac(xaxis, *popt)
#plt.plot(xaxis, fit, 'g-',ls='-')
    #plt.figtext(0.81, 0.25, graphDescription)
    #plt.savefig('figures/'+ f[f.find('/'):-4] + '.png')
    #plt.savefig('figures/'+ f[:-4] + '.png')
    #plt.savefig('figures/'+ graphTitle.replace('/','_') + '.png', dpi = 100)
    #plt.savefig('figures/eps/'+ graphTitle.replace('/','_') + '.eps', dpi = 150)
    plt.savefig('figures/'+ graphTitle.strip().replace('\\','') + '.pdf',dpi=100)

print "\a"
