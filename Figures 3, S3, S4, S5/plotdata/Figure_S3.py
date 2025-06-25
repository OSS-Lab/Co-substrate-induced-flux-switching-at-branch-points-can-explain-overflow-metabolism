#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 16:10:43 2025

@author: robert
"""



# import os
import sys
# import inspect

# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0, parentdir) 

sys.path.append('../simulation_functions')
sys.path.append('../data')

import numpy as np

import matplotlib.pyplot as plt

# import make_odes as mp
# import time_series_plots as ts
import matplotlib.pylab as pylab
import make_figure_data as mfd
import matplotlib.image as img
import matplotlib as mpl



savefig=0               # change to 1 to save generated figure
num_kin=201
# num_kin=51
# mfd.FIGURE2_2_ff_vs_kin_vary24()

cmapuse='RdBu'


params = {'legend.fontsize': 30,
          'figure.figsize': (12, 9),
          'axes.labelsize': 40,
          'axes.titlesize': 40,
          'xtick.labelsize':30,
          'ytick.labelsize':30}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

# Btot_vals=[10**-1,10**-2,10**1]
# VMB_vals=[10**-1,10**1.5,10**-1]
# Btot_vals=[10**-1,10**-1,10**1]
Btot_vals=[10**-2,10**-2,10**-1]
VMB_vals=[10**-1,10**1,10**-1]



Btot_vals2=[10**0,10**0,10**1]
VMB_vals2=[10**-2,10**-1,10**-2]

KeqB_val=1.
CKcat_top_val=0.1
CKcat_bottom_val=1.
cols=['k','r','b']


plt.figure()

ax1=plt.axes([0.1,0.05,0.5,0.5])
inset1=plt.axes([0.1025,0.345,0.23,0.25]) 

motifnone=img.imread('fig3_symm.png')
plt.imshow(motifnone)
ax1.text(10**-6,0.85,'$\\rm A$', fontsize=40)
plt.axis('off')


#############################


ax2=plt.axes([0.75,0.05,0.5,0.5])
inset2=plt.axes([0.7525,0.345,0.23,0.25]) 

motifnone=img.imread('fig3_symm_upstream.png')
plt.imshow(motifnone)
ax2.text(10**-6,0.85,'$\\rm B$', fontsize=40)
plt.axis('off')

#############################


ax3=plt.axes([0.1,-0.55,0.5,0.5])
# inset3=plt.axes([0.105,-0.5,0.23,0.25]) 
inset3=plt.axes([0.105,-0.26,0.23,0.25]) 

motifnone=img.imread('fig3_minimal.png')
plt.imshow(motifnone)
ax3.text(10**-6,0.85,'$\\rm C$', fontsize=40)
plt.axis('off')

#############################


ax4=plt.axes([0.75,-0.55,0.5,0.5])
inset4=plt.axes([0.755,-0.26,0.23,0.25]) 

motifnone=img.imread('fig3_minimal_upstream.png')
plt.imshow(motifnone)
ax4.text(10**-6,0.85,'$\\rm D$', fontsize=40)
plt.axis('off')



xvar='Btot'
yvar='CKcatB'

# KeqB_val=1.



reversible=False


CKcat=1.
Ks=1.
Kp=1.
Keq=1.               # params_global
KsB=1.
KpB=1.                             # params_b unchanging
Ks_top=1.
Kp_top=1.
Ks_bottom=1. 
Kp_bottom=1.                # params_branch unchanging
kout_top=1.
kout_bottom=1.

kin_min=-5
kin_max=0
# num_kin=41

Tmax=100000
tpts=10000

titlefontsize=20
rev=False
####### these dont change in the irrev case
Keq_top=1.
Keq_bottom=1.

kinvals=np.logspace(kin_min,kin_max,num_kin)

lengthparams=[2,1,1]
B_position_params=([0,0],[0,0],[1],[1])

for idb,Btot in enumerate(Btot_vals):
    VMB=VMB_vals[idb]
    col=cols[idb]
    print(Btot,VMB,col)
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                       B_position_params,
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB_val,Btot],
                          params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
        
        
    ax1.plot(kinvals,ff_vec,color=col,label='$B_T = '+str(Btot)+', ~V_M^B = '+str(VMB)+'$')#,label='$B_T = 10^{'+str(np.log10(Btot))+'}, V_M^B = 10^{'+str(np.log10(VMB))+'}$')
ax1.set_xscale('log')
ax1.set_xlabel('$v_{\\rm in}$')    
ax1.set_ylabel('$\\rm Flux~ Fraction$',rotation=90)  
ax1.set_ylim(0,1.05) 
ax1.set_xlim(10**-5.1,10**0)
ax1.legend(loc='upper right',fontsize=15)
ax1.xaxis.set_label_coords(0.85,-0.03)
ax1.yaxis.set_label_coords(-0.15,0.3)


####################################

lengthparams=[2,1,1]
B_position_params=([1,0],[0,0],[1],[1])


for idb,Btot in enumerate(Btot_vals):
    VMB=VMB_vals[idb]
    col=cols[idb]
    print(Btot,VMB,col)
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                       B_position_params,
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB_val,Btot],
                          params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
        
        
    ax2.plot(kinvals,ff_vec,color=col,label='$B_T = '+str(Btot)+', V_M^B = '+str(VMB)+'$')
ax2.set_xscale('log')
ax2.set_xlabel('$v_{\\rm in}$')    
ax2.set_ylabel('$\\rm Flux~ Fraction$',rotation=90)  
ax2.set_ylim(0,1.05) 
ax2.set_xlim(10**-5.1,10**0)
# ax2.legend(loc='lower right',fontsize=15)
ax2.xaxis.set_label_coords(0.85,-0.03)
ax2.yaxis.set_label_coords(-0.15,0.3)


# ####################################

lengthparams=[2,1,1]
B_position_params=([0,0],[0,0],[0],[1])


for idb,Btot in enumerate(Btot_vals2):
    VMB=VMB_vals2[idb]
    col=cols[idb]
    print(Btot,VMB,col)
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                       B_position_params,
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB_val,Btot],
                          params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
        
        
    ax3.plot(kinvals,ff_vec,color=col,label='$B_T = '+str(Btot)+',~ V_M^B = '+str(VMB)+'$')
ax3.set_xscale('log')
ax3.set_xlabel('$v_{\\rm in}$')    
ax3.set_ylabel('$\\rm Flux~ Fraction$',rotation=90)  
ax3.set_ylim(0,1.05) 
ax3.set_xlim(10**-5.1,10**0)
# ax3.legend(loc='upper right',fontsize=15)
ax3.legend(loc=(0.025,0.4),fontsize=15)
ax3.xaxis.set_label_coords(0.85,-0.03)
ax3.yaxis.set_label_coords(-0.15,0.3)


# ####################################

lengthparams=[2,1,1]
B_position_params=([1,0],[0,0],[0],[1])


for idb,Btot in enumerate(Btot_vals2):
    VMB=VMB_vals2[idb]
    col=cols[idb]
    print(Btot,VMB,col)
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                       B_position_params,
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB_val,Btot],
                          params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
        
        
    ax4.plot(kinvals,ff_vec,color=col,label='$B_T = '+str(Btot)+', V_M^B = '+str(VMB)+'$')
ax4.set_xscale('log')
ax4.set_xlabel('$v_{\\rm in}$')    
ax4.set_ylabel('$\\rm Flux~ Fraction$',rotation=90)  
ax4.set_ylim(0,1.05) 
ax4.set_xlim(10**-5.1,10**0)
# ax4.legend(loc='lower right',fontsize=15)
ax4.xaxis.set_label_coords(0.85,-0.03)
ax4.yaxis.set_label_coords(-0.15,0.3)


if savefig==1:
    plt.savefig('Figure_S3.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S3.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S3.svg',bbox_inches='tight',format='svg')
    plt.savefig('Figure_S3.pdf',bbox_inches='tight',format='pdf')

plt.show()
