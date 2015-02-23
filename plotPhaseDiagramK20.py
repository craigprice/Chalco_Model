#!/usr/bin/python

import numpy as np
import pylab as plt
from matplotlib import collections  as mc

fig, ax = plt.subplots()


zColor = 'lightgreen'

lineZ1 = [(-0.74, 1), (-0.2, 0.32)]
plt.fill_between((lineZ1[0][0],lineZ1[1][0]),(lineZ1[0][1],lineZ1[1][1]), (1,1), color=zColor)
lineZ2 = [lineZ1[1], (0.16, 0.38)]
plt.fill_between((lineZ2[0][0],lineZ2[1][0]),(lineZ2[0][1],lineZ2[1][1]), (1,1), color=zColor)
lineZ3 = [lineZ2[1], (1, 0.66)]
plt.fill_between((lineZ3[0][0],lineZ3[1][0]),(lineZ3[0][1],lineZ3[1][1]), (1,1), color=zColor)

mColor = 'lightblue'
lineM1 = [(-1, 0.25), (-0.62, 0)]
plt.fill_between((lineM1[0][0],lineM1[1][0]),(lineM1[0][1],lineM1[1][1]), color=mColor)

sColor = 'pink'
lineS1 = [(-0.38, 0), lineZ1[1]]
plt.fill_between((lineS1[0][0],lineS1[1][0]),(lineS1[0][1],lineS1[1][1]), (0,0), color=sColor)
lineS2 = [lineS1[1], (0.16, 0.26)]
plt.fill_between((lineS2[0][0],lineS2[1][0]),(lineS2[0][1],lineS2[1][1]), (0,0), color=sColor)
lineS3 = [lineS2[1], (0.86, 0)]
plt.fill_between((lineS3[0][0],lineS3[1][0]),(lineS3[0][1],lineS3[1][1]), (0,0), color=sColor)

qColor = 'cyan'
plt.fill_between((lineZ2[0][0],lineZ2[1][0]), (lineS2[0][1],lineS2[1][1]), (lineZ2[0][1],lineZ2[1][1]), color=qColor)
plt.fill_between((lineZ3[0][0],lineZ3[1][0]), (lineS3[0][1],-0.05), (lineZ3[0][1],lineZ3[1][1]), color=qColor)

cColor = 'Turquoise'
lineC1 = [(0.16, 0.38), (1,0.14)]
plt.fill_between((lineS3[0][0],lineS3[1][0]),(lineS3[0][1],lineS3[1][1]), (lineC1[0][1],0.18), color=cColor)
plt.fill_between((lineS3[1][0],1),(0,0), (0.18,lineC1[1][1]), color=cColor)
plt.fill_between((lineS2[0][0],lineS2[1][0]),(lineS2[0][1],lineS2[1][1]), (lineZ2[0][1],lineZ2[1][1]), color=cColor)

linesShade = [[(lineZ3[0][0]-0.01,lineZ3[0][1]-0.01),lineZ3[1]]]

lines = [lineZ1, lineZ2, lineZ3, lineM1, lineS1, lineS2, lineS3, lineC1, lineS2, lineS3, lineZ2]
c =     [zColor, zColor, zColor, mColor, sColor, sColor, sColor, cColor, sColor, sColor, zColor]

lc = mc.LineCollection(lines, colors=c, linewidths=1)
lc2 = mc.LineCollection(linesShade, colors='Turquoise',linewidths=4)

ax.add_collection(lc)
ax.add_collection(lc2)
ax.autoscale()
ax.margins(0.1)

fontsizeLabel = 18
bigger = 1.4
plt.title('K$_2$=0', y=1.02, fontsize=bigger*fontsizeLabel)
plt.text(-0.95,0.9,'(a)', fontsize=bigger*fontsizeLabel)
plt.text(-0.95,0.05,'FM', fontsize=fontsizeLabel)
plt.text(-0.12,0.7,'Zigzag', fontsize=fontsizeLabel)
plt.text(-0.8,0.4,'3Q-Spiral', fontsize=fontsizeLabel)
plt.text(-0.12,0.15,'Stripy', fontsize=fontsizeLabel)
plt.text(0.5,0.4,'120$^\circ$ State', fontsize=fontsizeLabel)
plt.text(-0.12,0.3,'Intermediate', fontsize=fontsizeLabel)
plt.xlabel('J$_2$/J$_1$', fontsize = bigger*fontsizeLabel)
plt.ylabel('J$_3$/J$_1$', rotation=0, fontsize=bigger*fontsizeLabel)
plt.tick_params(axis = 'both',labelsize = bigger*fontsizeLabel)
ax.xaxis.set_label_coords(0.5, -0.08)
ax.yaxis.set_label_coords(-0.02, 0.5)


plt.ylim(0,1)
plt.xlim(-1,1)

plt.tight_layout()
plt.savefig('J2J3_K2=0_Phase_Diagram.png',dpi=600)
plt.savefig('J2J3_K2=0_Phase_Diagram.eps',dpi=600)

