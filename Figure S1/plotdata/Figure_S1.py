#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 12 22:33:15 2025

@author: robert
"""


import numpy as np

import matplotlib.pyplot as plt

import matplotlib.pylab as pylab
# import make_figure_data as mfd
# import matplotlib.image as img
import matplotlib as mpl
import matplotlib.image as img

savefig=0

params = {'legend.fontsize': 30,
          'figure.figsize': (12, 9),
          'axes.labelsize': 'x-large',
          'axes.titlesize': 20,
          'xtick.labelsize':20,
          'ytick.labelsize':20,
          'legend.title_fontsize':15}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath}')




plt.figure()


ax1=plt.axes([0.03,0.55,0.62,0.5])
motifnone=img.imread('panel_A.png')
# img[80:280, 150:330]
plt.imshow(motifnone[10:1869, 10:3986])
# ax1.axis([0,1,0,1])
plt.axis('off')

ax2=plt.axes([0.03,0.115,0.62,0.55])
plt.axis('off')
motifbranch=img.imread('panel_B.png')
plt.imshow(motifbranch[10:2034, 10:4014])

ax3=plt.axes([0.03,-0.25,0.62,0.5])
plt.axis('off')
motiffull=img.imread('panel_C.png')
plt.imshow(motiffull[10:1869, 10:4390])

ax4=plt.axes( [0.65,-0.18,0.65,1.18])
plt.axis('off')
motifeq=img.imread('panel_D.png')
plt.imshow(motifeq)

print(motifnone.shape)
print(motifbranch.shape)
print(motiffull.shape)
print(motifeq.shape)



if savefig==1:
    plt.savefig('Figure_S1.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S1.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S1.svg',bbox_inches='tight',format='svg')
    plt.savefig('Figure_S1.pdf',bbox_inches='tight',format='pdf')
    
plt.show()
