#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 12 11:48:19 2025

@author: robert
"""



import numpy as np

import time
from datetime import timedelta
import multiprocessing as multi
import sys

#from pathlib import Path
# import os
# import time
from pathlib import Path

# import itertools

# sys.path.append('../data/megadata')
# sys.path.append('../data/heatmapdata')
sys.path.append('../simulation_functions')
sys.path.append('../data')


# import make_heatmaps_compartment as mh_cm
# import ode_funcs_compartment as ode_cm
# import matplotlib.image as img

import matplotlib.pyplot as plt
# import numpy.random as rn
# import datetime as dt

import matplotlib.pylab as pylab
import matplotlib as mpl

savefig=0            #change to 1 to save generated figure
# numkin=101
# numkin=41

params = {'legend.fontsize': 12,
          'legend.title_fontsize': 12,
          'axes.labelsize': 15,
          'axes.titlesize': 15,
          'xtick.labelsize':15,
          'ytick.labelsize':15}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

allo19_desiredafteroxidase=[11,70,78,85,115,116,129,151,158,202,240,254,255,
                            328,421,445,453,454,487,502,555,560,564,565,570,
                            571,602,643,674,681,697,718,739,785,791,874,978,
                            995]
binedges=np.linspace(0,1,21)
#202 421 445 454 487 570 874
indplot_switchretained=202 #570

indplot_noncosub=102
plt.figure()
# ax1=plt.axes([-0.05,-0.05,0.9,0.6])       # motif
# motifnone=img.imread('motif_compartment_orig.png')
# plt.imshow(motifnone)
# plt.axis('off')

# ax2=plt.axes([1,-0.6,0.9,1.15])    # scatter
# ax3=plt.axes([-0.04,-0.65,0.4,0.48])    # fluxfraction
# ax4=plt.axes([0.44,-0.65,0.4,0.48])    # expt

# # ax3=plt.axes([2,-0.6,0.9,1.1])    # fluxfraction
# # ax4=plt.axes([3,-0.6,0.9,1.1])    # expt

# # ax5=plt.axes([0.05,-0.6,0.8,0.5])   # box plots

# # ax5.text(0.2,0.45,'Experimental data',fontsize=30)



# ax5=plt.axes([2.0,0.02,0.2,0.48])   # box plots
# ax6=plt.axes([2.325,0.02,0.2,0.48])
# ax7=plt.axes([2.65,0.02,0.2,0.48])

# ax10=plt.axes([2.0,-0.67,0.2,0.48])   # box plots
# ax11=plt.axes([2.325,-0.67,0.2,0.48])
# ax12=plt.axes([2.65,-0.67,0.2,0.48])
# ax8=plt.axes([2.975,-0.67,0.2,0.48])
# ax9=plt.axes([3.3,-0.67,0.2,0.48])


#ax14=plt.axes([1.05,-1.35,0.85,0.5])
#ax15=plt.axes([1.05,-2,0.85,0.5])

ax14=plt.axes([0,0,0.9,1.3])
ax15=plt.axes([1,0,0.9,1.3])

ax17=plt.axes([0.17,0.75,0.4,0.5])
ax16=plt.axes([1.45,0.75,0.4,0.5])
# ax13=plt.axes([1.4,0.5,0.45,0.6])    # scatter 2
# ax16=plt.axes([0.15,0.5,0.45,0.6])    # scatter 2


ax15.text(-1.37,7.9,'$\\rm A$',fontsize=40)
ax15.text(-0.15,8.5,'$\\rm B$',fontsize=40)
# ax16=plt.axes([2.05,-1.35,0.85,0.5])
# ax17=plt.axes([2.05,-2,0.85,0.5])

# get the data for new params
allo=[1.9,1.,1.,1.]
n=1000

seed=1000
np.random.seed(seed)
if seed!=False:
    seedstr=str(seed)
else:
    seedstr='none'

delwant_max=1
delwant_min=0.5
n0deltawant_max=1.
n0deltawant_min=0.5
delwant_noswitch=0.4

# args_to_run=mh_cm.make_params_newtest_3_nonfixed(allo, n)

      
savepath='../data/allo_'+str(allo)

filename='nonfixed_apr_final_seed'+seedstr+'_n_'+str(n)
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewant)

# kinvals=args_to_run[0][1]

dat=np.load(filewant)
ff_delta=dat['fluxfrac_delta']
gradmaxind=dat['gradmaxind']
sensemaxind=dat['sensemaxind']
fluxfrac=dat['fluxfrac']
ff_grad=dat['fluxfrac_grad']
sensemax=dat['sensemax']
gradmax=dat['gradmax']
ffsense_vec=dat['fluxfrac_sense']
b0delta=dat['b0delta']
c0delta=dat['c0delta']
n0delta=dat['n0delta']

ff_delta=np.abs(ff_delta)
# ax9.grid('on')


# ax16.scatter(ff_delta,n0delta,s=20,color='b',alpha=0.5,edgecolors='k',zorder=20)#,label='')
# ax16.axis([-0.05,1.05,-0.05,1.05])
# # ax2.axis([-1,1,-1,1])
# ax16.set_xlabel('$\\Delta F_f$',fontsize=20)
# ax16.set_ylabel('$\\Delta N_{\\rm H}$',fontsize=20,rotation=0)
# ax16.yaxis.set_label_coords(-0.15,0.6)
# ax16.xaxis.set_label_coords(0.7,-0.05)
# ax16.set_xticks(np.linspace(0,1,3))
# ax16.set_yticks(np.linspace(0,1,3))

ax14.hist(ff_delta,bins=binedges,alpha=1,color='cornflowerblue',edgecolor='cornflowerblue',density=True, 
            label='With CRs')
# ax14.set_xlabel('flux frac delta')
ax14.set_xlabel('$\\Delta F_f$',fontsize=30)
ax14.set_ylabel('$\\rm Density$',fontsize=20)

ax17.scatter(ff_delta,c0delta,s=10,color='cornflowerblue')
ax17.axis([-0.05,1.05,-0.05,1.05])
ax17.set_xlabel('$\\Delta F_f$',fontsize=20)
ax17.set_ylabel('$\\Delta N_{\\rm H}^{m}$',fontsize=20,rotation=0)
ax17.yaxis.set_label_coords(-0.15,0.7)
ax17.xaxis.set_label_coords(0.7,-0.1)
   

filename='nocosub_nonfixed_apr_final_seed'+seedstr+'_n_'+str(n)
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewant)

# kinvals=args_to_run[0][1]

dat=np.load(filewant)
ff_delta=dat['fluxfrac_delta']
gradmaxind=dat['gradmaxind']
sensemaxind=dat['sensemaxind']
fluxfrac=dat['fluxfrac']
ff_grad=dat['fluxfrac_grad']
sensemax=dat['sensemax']
gradmax=dat['gradmax']
ffsense_vec=dat['fluxfrac_sense']
b0delta=dat['b0delta']
c0delta=dat['c0delta']
n0delta=dat['n0delta']

ff_delta=np.abs(ff_delta)
c0delta=np.abs(c0delta)
# ax9.grid('on')



ax14.hist(ff_delta,bins=binedges,alpha=1,fill=False,edgecolor='tomato',density=True, 
            label='Without CRs', linewidth =2)
# ax14.set_xlabel('flux frac delta')
ax14.legend(loc=(0.65,0.85),title='\\underline{Parameters as in Table 1}')
ax14.axis([-0.05,1.05,0,2.75])

allo=[1.,1.,1.,1.]
savepath='../data/allo_'+str(allo)

filename='fixed_apr_final_seed'+seedstr+'_n_'+str(n)
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewant)

# kinvals=args_to_run[0][1]

dat=np.load(filewant)
ff_delta=dat['fluxfrac_delta']
ff_delta=np.abs(ff_delta)
gradmaxind=dat['gradmaxind']
sensemaxind=dat['sensemaxind']
fluxfrac=dat['fluxfrac']
ff_grad=dat['fluxfrac_grad']
sensemax=dat['sensemax']
gradmax=dat['gradmax']
ffsense_vec=dat['fluxfrac_sense']
b0delta=dat['b0delta']
c0delta=dat['c0delta']
n0delta=dat['n0delta']

# ax13.scatter(ff_delta,n0delta,s=20,color='seagreen',alpha=0.5,edgecolors='k',zorder=2)
# # ax2.axis([-1,1,-1,1])
# ax13.set_xlabel('$\\Delta F_f$',fontsize=20)
# ax13.set_ylabel('$\\Delta N_{\\rm H}$',fontsize=20,rotation=0)
# ax13.yaxis.set_label_coords(-0.15,0.3)
# ax13.xaxis.set_label_coords(0.7,-0.05)
# ax13.set_xticks(np.linspace(0,1,3))
# ax13.set_yticks(np.linspace(0,1,3))

# ax15.hist(ff_delta,bins=binedges,alpha=0.5,color='seagreen',density=True, 
#             label='Table 2, with \n co-substrates')


ax15.hist(ff_delta,bins=binedges,color='mediumseagreen',
          edgecolor='mediumseagreen',density=True, 
            label="With CRs")
ax15.set_xlabel('$\\Delta F_f$',fontsize=30)
ax15.set_ylabel('$\\rm Density$',fontsize=20)

ax16.scatter(ff_delta,c0delta,s=10,color='mediumseagreen')
ax16.axis([-0.05,1.05,-0.05,1.05])
ax16.set_xlabel('$\\Delta F_f$',fontsize=20)
ax16.set_ylabel('$\\Delta N_{\\rm H}^{m}$',fontsize=20,rotation=0)
ax16.yaxis.set_label_coords(-0.15,0.7)
ax16.xaxis.set_label_coords(0.7,-0.1)

filename='nocosub_fixed_apr_final_seed'+seedstr+'_n_'+str(n)
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewant)

# kinvals=args_to_run[0][1]

dat=np.load(filewant)
ff_delta=dat['fluxfrac_delta']
ff_delta=np.abs(ff_delta)
gradmaxind=dat['gradmaxind']
sensemaxind=dat['sensemaxind']
fluxfrac=dat['fluxfrac']
ff_grad=dat['fluxfrac_grad']
sensemax=dat['sensemax']
gradmax=dat['gradmax']
ffsense_vec=dat['fluxfrac_sense']
b0delta=dat['b0delta']
c0delta=dat['c0delta']
n0delta=dat['n0delta']



# ax15.hist(ff_delta,bins=binedges,alpha=0.5,color='r',density=True,
#             label='Table 2, without \n co-substrates')
ax15.axis([-0.05,1.05,0,9.9])


ax15.hist(ff_delta,bins=binedges,edgecolor='tomato',density=True,fill=False,
            label="Without CRs",linewidth=2)
ax15.legend(loc=(0.12,0.35),title='\\underline{Branch-point enzyme} \n \\underline{parameters set equal}')

if savefig==1:
    plt.savefig('Figure_S10.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S10.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S10.pdf',bbox_inches='tight',format='pdf')
    plt.savefig('Figure_S10.svg',bbox_inches='tight',format='svg')

plt.show()