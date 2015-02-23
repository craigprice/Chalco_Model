#!/usr/bin/python

import numpy as np
import pylab as plt
from matplotlib import collections  as mc

fig, ax = plt.subplots()
#plt.figure(figsize=(6,3))

zColor = 'lightgreen'
lineZ1 = [(-0.34, 0), (0.12, 0.48)]
plt.fill_between((lineZ1[0][0],lineZ1[1][0]),(lineZ1[0][1],lineZ1[1][1]), (1,1), color=zColor)
plt.fill_between((-1,lineZ1[0][0]),(0,0), (1,1), color=zColor)
lineZ2 = [lineZ1[1], (0.5, 1)]
plt.fill_between((lineZ2[0][0],lineZ2[1][0]),(lineZ2[0][1],lineZ2[1][1]), (1,1), color=zColor)

sColor = 'pink'
lineS1 = [(-0.18, 0), (0.5, 1)]
plt.fill_between((lineS1[0][0],lineS1[1][0]),(lineS1[0][1],lineS1[1][1]), (0,0), color=sColor)
plt.fill_between((lineS1[1][0],1),(0,0), (1,1), color=sColor)

qColor = 'cyan'
plt.fill_between((-0.34,lineS1[0][0]),         (0,0), (0,0.16), color=qColor)
plt.fill_between((lineS1[0][0],lineZ1[1][0]), (lineS1[0][1],0.45), (0.16,lineZ1[1][1]), color=qColor)
plt.fill_between((lineZ1[1][0],0.36),         (0.44,0.8), (0.47,0.8), color=qColor)

cColor = 'Turquoise'

lines = [lineZ1, lineZ2, lineS1]

c = [cColor, cColor,cColor]

lc = mc.LineCollection(lines, colors=c, linewidths=2)
ax.add_collection(lc)
ax.autoscale()
ax.margins(0.15)

fontsizeLabel = 18
bigger = 1.4
plt.title('K$_2$=-2J$_2$', y=1.02, fontsize=bigger*fontsizeLabel)
plt.text(-0.95,0.9,'(b)', fontsize=bigger*fontsizeLabel)
plt.text(-0.7,0.7,'Zigzag', fontsize=fontsizeLabel)
plt.text(0.5,0.3,'Stripy', fontsize=fontsizeLabel)
plt.text(-0.28,0.25,'120$^\circ$ State', fontsize=fontsizeLabel,rotation=57)
plt.xlabel('J$_2$/J$_1$', fontsize = bigger*fontsizeLabel)
plt.ylabel('J$_3$/J$_1$', rotation=0, fontsize=bigger*fontsizeLabel)
plt.tick_params(axis = 'both',labelsize = bigger*fontsizeLabel)
ax.xaxis.set_label_coords(0.5, -0.08)
ax.yaxis.set_label_coords(-0.02, 0.5)

plt.ylim(0,1)
plt.xlim(-1,1)

xs = -0.75
ys = 0.2
plt.plot(xs, ys, 'r*', markersize=20)

plt.tight_layout()
plt.savefig('J2J3_K2=-2J2_Phase_Diagram.png',dpi=600)
plt.savefig('J2J3_K2=-2J2_Phase_Diagram.eps',dpi=600)

