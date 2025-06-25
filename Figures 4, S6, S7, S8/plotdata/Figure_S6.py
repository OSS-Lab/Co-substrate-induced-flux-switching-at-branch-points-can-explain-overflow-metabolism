#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 22 14:27:20 2025

@author: robert
"""


import numpy as np

import time
from datetime import timedelta
import multiprocessing as multi
import sys

#from pathlib import Path
import os
import time
from pathlib import Path

import itertools

# sys.path.append('../data/megadata')
# sys.path.append('../data/heatmapdata')
# sys.path.append('../simulation_functions')
sys.path.append('../data')


import make_heatmaps_compartment as mh_cm
import ode_funcs_compartment as ode_cm
import matplotlib.image as img

import matplotlib.pyplot as plt
import numpy.random as rn
import datetime as dt

import matplotlib.pylab as pylab
import matplotlib as mpl

savefig=0       # change to 1 to save generated figures
# numkin=101
numkin=41

params = {'legend.fontsize': 12,
          'axes.labelsize': 15,
          'axes.titlesize': 15,
          'xtick.labelsize':15,
          'ytick.labelsize':15}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

# allo19_desiredafteroxidase=[11,70,78,85,115,116,129,151,158,202,240,254,255,
#                             328,421,445,453,454,487,502,555,560,564,565,570,
#                             571,602,643,674,681,697,718,739,785,791,874,978,
#                             995]

#202 421 445 454 487 570 874
# indplot_switchretained=202 #570

# indplot_noncosub=102
plt.figure()
# ax1=plt.axes([-0.05,-0.05,0.9,0.6])       # motif
# motifnone=img.imread('motif_compartment_orig.png')
# plt.imshow(motifnone)
# plt.axis('off')

ax2=plt.axes([1,-0.6,0.9,1.15])    # scatter
# ax3=plt.axes([-0.04,-0.65,0.4,0.48])    # fluxfraction
# ax4=plt.axes([0.44,-0.65,0.4,0.48])    # expt

# ax3=plt.axes([2,-0.6,0.9,1.1])    # fluxfraction
# ax4=plt.axes([3,-0.6,0.9,1.1])    # expt

# ax5=plt.axes([0.05,-0.6,0.8,0.5])   # box plots

# ax5.text(0.2,0.45,'Experimental data',fontsize=30)



ax5=plt.axes([2.0,0.02,0.3,0.48])   # box plots
ax6=plt.axes([2.4,0.02,0.3,0.48])
ax7=plt.axes([2.8,0.02,0.3,0.48])

ax10=plt.axes([2.18,-0.67,0.3,0.48])   # box plots
ax11=plt.axes([2.62,-0.67,0.3,0.48])
# ax12=plt.axes([2.65,-0.67,0.2,0.48])
# ax8=plt.axes([2.975,-0.67,0.2,0.48])
# ax9=plt.axes([3.3,-0.67,0.2,0.48])

# ax13=plt.axes([0,-2,0.9,1.15])    # scatter 2

# ax14=plt.axes([1.05,-1.35,0.85,0.5])
# ax15=plt.axes([1.05,-2,0.85,0.5])
# ax16=plt.axes([2.05,-1.35,0.85,0.5])
# ax17=plt.axes([2.05,-2,0.85,0.5])
ax2.text(-0.2,1.05,'$\\rm A$',fontsize=40)
ax2.text(1.09,1.05,'$\\rm B$',fontsize=40)
ax2.text(1.22,0.35,'$\\rm C$',fontsize=40)

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

args_to_run=mh_cm.make_params_newtest_6_nonfixed(allo, n)

      
savepath='../data/allo_'+str(allo)

filename='nonfixed_apr_final_seed'+seedstr+'_n_'+str(n)
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewant)

kinvals=args_to_run[0][1]

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



para01_v=list()
para01_ks=list()
para01_kp=list()
para01_keq=list()


para1E_v=list()
para1E_ks=list()
para1E_kp=list()
para1E_keq=list()

paraE2_v=list()
paraE2_ks=list()
paraE2_kp=list()
paraE2_keq=list()

para23_v=list()
para23_ks=list()
# para23_kp=list()
# para23_keq=list()

para35_v=list()
para35_ks=list()
para35_kp=list()
para35_keq=list()

para_24_v=list()
para_24_ks=list()
# para_24_kp=list()
# para_24_keq=list()

para_35_v=list()
para_35_ks=list()
para_35_kp=list()
para_35_keq=list()

para_t_v=list()
para_t_ks=list()

para_diff_3_in=list()
para_diff_3_out=list()
para_diff_5_in=list()
para_diff_5_out=list()

paraB_v=list()
paraB_ks=list()
paraB_kp=list()
paraB_keq=list()

paraC_v=list()
paraC_ks=list()
paraC_kp=list()
paraC_keq=list()

btot=list()
ctot=list()

btot_p=list()
ctot_p=list()



para01_v_n=list()
para01_ks_n=list()
para01_kp_n=list()
para01_keq_n=list()

para1E_v_n=list()
para1E_ks_n=list()
para1E_kp_n=list()
para1E_keq_n=list()


paraE2_v_n=list()
paraE2_ks_n=list()
paraE2_kp_n=list()
paraE2_keq_n=list()


para23_v_n=list()
para23_ks_n=list()
# para23_kp_n=list()
# para23_keq_n=list()

para35_v_n=list()
para35_ks_n=list()
para35_kp_n=list()
para35_keq_n=list()

para_24_v_n=list()
para_24_ks_n=list()
# para_24_kp_n=list()
# para_24_keq_n=list()

para_35_v_n=list()
para_35_ks_n=list()
para_35_kp_n=list()
para_35_keq_n=list()

para_t_v_n=list()
para_t_ks_n=list()

para_diff_3_in_n=list()
para_diff_3_out_n=list()
para_diff_5_in_n=list()
para_diff_5_out_n=list()

paraB_v_n=list()
paraB_ks_n=list()
paraB_kp_n=list()
paraB_keq_n=list()

paraC_v_n=list()
paraC_ks_n=list()
paraC_kp_n=list()
paraC_keq_n=list()

btot_p_n=list()
ctot_p_n=list()

para01=list()
para1E=list()
paraE2=list()
para23=list()
para35=list()
para_24=list()
para_35=list()

para_t=list()
para_diff_3=list()
para_diff_5=list()

paraB=list()
paraC=list()

para01_n=list()
para1E_n=list()
paraE2_n=list()
para23_n=list()
para35_n=list()
para_24_n=list()
para_35_n=list()

para_t_n=list()
para_diff_3_n=list()
para_diff_5_n=list()

paraB_n=list()
paraC_n=list()

btot_n=list()
ctot_n=list()

ff_delta=np.abs(ff_delta)
c0delta=np.abs(c0delta)
b0delta=np.abs(b0delta)
n0delta=np.abs(n0delta)
indexes_noswitch=list()



para01_v_n2=list()
para01_ks_n2=list()
para01_kp_n2=list()
para01_keq_n2=list()

para1E_v_n2=list()
para1E_ks_n2=list()
para1E_kp_n2=list()
para1E_keq_n2=list()


paraE2_v_n2=list()
paraE2_ks_n2=list()
paraE2_kp_n2=list()
paraE2_keq_n2=list()


para23_v_n2=list()
para23_ks_n2=list()
# para23_kp_n2=list()
# para23_keq_n2=list()

para35_v_n2=list()
para35_ks_n2=list()
para35_kp_n2=list()
para35_keq_n2=list()

para_24_v_n2=list()
para_24_ks_n2=list()
# para_24_kp_n2=list()
# para_24_keq_n2=list()

para_35_v_n2=list()
para_35_ks_n2=list()
para_35_kp_n2=list()
para_35_keq_n2=list()

para_t_v_n2=list()
para_t_ks_n2=list()

para_diff_3_in_n2=list()
para_diff_3_out_n2=list()
para_diff_5_in_n2=list()
para_diff_5_out_n2=list()

paraB_v_n2=list()
paraB_ks_n2=list()
paraB_kp_n2=list()
paraB_keq_n2=list()

paraC_v_n2=list()
paraC_ks_n2=list()
paraC_kp_n2=list()
paraC_keq_n2=list()

btot_p_n2=list()
ctot_p_n2=list()

para01=list()
para1E=list()
paraE2=list()
para23=list()
para35=list()
para_24=list()
para_35=list()

para_t=list()
para_diff_3=list()
para_diff_5=list()

paraB=list()
paraC=list()

para01_n=list()
para1E_n=list()
paraE2_n=list()
para23_n=list()
para35_n=list()
para_24_n=list()
para_35_n=list()

para_t_n=list()
para_diff_3_n=list()
para_diff_5_n=list()

paraB_n=list()
paraC_n=list()

btot_n=list()
ctot_n=list()

para01_n2=list()
para1E_n2=list()
paraE2_n2=list()
para23_n2=list()
para35_n2=list()
para_24_n2=list()
para_35_n2=list()

para_t_n2=list()
para_diff_3_n2=list()
para_diff_5_n2=list()

paraB_n2=list()
paraC_n2=list()

btot_n2=list()
ctot_n2=list()

indexes_switch=list()
indexes_switch_retained=list()

ind_notpicked=list()
for i in range(n):

    # if ff_delta[i]<delwant_min:
    #     if ff_delta[i]>delwant_max:
    pinit=args_to_run[i][0]
    kinvals=args_to_run[i][1]
    params_out=args_to_run[i][2]
    
    params01=args_to_run[i][3]
    params1E=args_to_run[i][4]
    paramsE2=args_to_run[i][5]
    params23=args_to_run[i][6]
    params35=args_to_run[i][7]
    params_24=args_to_run[i][8]
    
    params_35=args_to_run[i][9]
    params_trans=args_to_run[i][10]
    params_diff_3=args_to_run[i][11]
    params_diff_5=args_to_run[i][12]
    
    paramsB=args_to_run[i][13]
    paramsC=args_to_run[i][14]
    allo=args_to_run[i][15]
    Tmax=args_to_run[i][16]
    
    # print(params35)
    # print(params_35)
    # if ff_delta[i]>delwant_min and ff_delta[i]<delwant_max and \
    #     c0delta[i]>n0deltawant_min and c0delta[i]<n0deltawant_max and \
    #     n0delta[i]>n0deltawant_min and n0delta[i]<n0deltawant_max:
            # print(i)
          #  count+=1
            # print(args_to_run[i])
   
        
        #tpts=args_to_run[i][16]
    
    #if ff_delta[i]>0.5 and ff_delta[i]<1.0 and \
        # c0delta[i]>0.3 and c0delta[i]<1.0: 
    if ff_delta[i]>0.5 and c0delta[i]>0.2: 
    # if i in allo19_desiredafteroxidase and n0delta[i]>0.2:        
        para01.append(params01)
        para1E.append(params1E)
        paraE2.append(paramsE2)
        para23.append(params23)
        para35.append(params35)
    
        para_24.append(params_24)
        para_35.append(params_35)
    
        para_t.append(params_trans)
        para_diff_3.append(params_diff_3)
        para_diff_5.append(params_diff_5)
        
        paraB.append(paramsB)
        paraC.append(paramsC)
    
        ctot.append(pinit[-3]+pinit[-2])
        btot.append(pinit[-5]+pinit[-4])
        
        
        para01_v.append(np.log10(1000*params01[0]))
        para01_ks.append(np.log10(params01[1]))
        para01_kp.append(np.log10(params01[2]))
        para01_keq.append(np.log10(params01[3]))
        
        para1E_v.append(np.log10(1000*params1E[0]))
        para1E_ks.append(np.log10(params1E[1]))
        para1E_kp.append(np.log10(params1E[2]))
        para1E_keq.append(np.log10(params1E[3]))
        
        paraE2_v.append(np.log10(1000*paramsE2[0]))
        paraE2_ks.append(np.log10(paramsE2[1]))
        paraE2_kp.append(np.log10(paramsE2[2]))
        paraE2_keq.append(np.log10(paramsE2[3]))
        
        para23_v.append(np.log10(1000*params23[0]))
        para23_ks.append(np.log10(params23[1]))
        # para23_kp.append(np.log10(params23[2]))
        # para23_keq.append(np.log10(params23[3]))
        
        para35_v.append(np.log10(1000*params35[0]))
        para35_ks.append(np.log10(params35[1]))
        para35_kp.append(np.log10(params35[2]))
        para35_keq.append(np.log10(params35[3]))
        
        para_24_v.append(np.log10(1000*params_24[0]))
        para_24_ks.append(np.log10(params_24[1]))
        # para_24_kp.append(np.log10(params_24[2]))
        # para_24_keq.append(np.log10(params_24[3]))
        
        para_35_v.append(np.log10(1000*params_35[0]))
        para_35_ks.append(np.log10(params_35[1]))
        para_35_kp.append(np.log10(params_35[2]))
        para_35_keq.append(np.log10(params_35[3]))
        
        para_t_v.append(np.log10(1000*params_trans[0]))
        para_t_ks.append(np.log10(params_trans[1]))
        
        para_diff_3_in.append(np.log10(params_diff_3[0]))
        para_diff_3_out.append(np.log10(params_diff_3[1]))
        para_diff_5_in.append(np.log10(params_diff_5[0]))
        para_diff_5_out.append(np.log10(params_diff_5[1]))
        
        paraB_v.append(np.log10(1000*paramsB[0]))
        paraB_ks.append(np.log10(paramsB[1]))
        paraB_kp.append(np.log10(paramsB[2]))
        paraB_keq.append(np.log10(paramsB[3]))
        
        paraC_v.append(np.log10(1000*paramsC[0]))
        paraC_ks.append(np.log10(paramsC[1]))
        paraC_kp.append(np.log10(paramsC[2]))
        paraC_keq.append(np.log10(paramsC[3]))
        
        btot_p.append(np.log10(pinit[-5]+pinit[-4]))
        ctot_p.append(np.log10(pinit[-3]+pinit[-2]))
    

        indexes_switch_retained.append(i)
    
    elif ff_delta[i]<0.4:
            
        para01_v_n.append(np.log10(1000*params01[0]))
        para01_ks_n.append(np.log10(params01[1]))
        para01_kp_n.append(np.log10(params01[2]))
        para01_keq_n.append(np.log10(params01[3]))
        
        para1E_v_n.append(np.log10(1000*params1E[0]))
        para1E_ks_n.append(np.log10(params1E[1]))
        para1E_kp_n.append(np.log10(params1E[2]))
        para1E_keq_n.append(np.log10(params1E[3]))
        
        paraE2_v_n.append(np.log10(1000*paramsE2[0]))
        paraE2_ks_n.append(np.log10(paramsE2[1]))
        paraE2_kp_n.append(np.log10(paramsE2[2]))
        paraE2_keq_n.append(np.log10(paramsE2[3]))
        
        para23_v_n.append(np.log10(1000*params23[0]))
        para23_ks_n.append(np.log10(params23[1]))
        # para23_kp_n.append(np.log10(params23[2]))
        # para23_keq_n.append(np.log10(params23[3]))
        
        para35_v_n.append(np.log10(1000*params35[0]))
        para35_ks_n.append(np.log10(params35[1]))
        para35_kp_n.append(np.log10(params35[2]))
        para35_keq_n.append(np.log10(params35[3]))
        
        para_24_v_n.append(np.log10(1000*params_24[0]))
        para_24_ks_n.append(np.log10(params_24[1]))
        # para_24_kp_n.append(np.log10(params_24[2]))
        # para_24_keq_n.append(np.log10(params_24[3]))
        
        para_35_v_n.append(np.log10(1000*params_35[0]))
        para_35_ks_n.append(np.log10(params_35[1]))
        para_35_kp_n.append(np.log10(params_35[2]))
        para_35_keq_n.append(np.log10(params_35[3]))
        
        para_t_v_n.append(np.log10(1000*params_trans[0]))
        para_t_ks_n.append(np.log10(params_trans[1]))
        
        para_diff_3_in_n.append(np.log10(params_diff_3[0]))
        para_diff_3_out_n.append(np.log10(params_diff_3[1]))
        para_diff_5_in_n.append(np.log10(params_diff_5[0]))
        para_diff_5_out_n.append(np.log10(params_diff_5[1]))
        
        paraB_v_n.append(np.log10(1000*paramsB[0]))
        paraB_ks_n.append(np.log10(paramsB[1]))
        paraB_kp_n.append(np.log10(paramsB[2]))
        paraB_keq_n.append(np.log10(paramsB[3]))
        
        paraC_v_n.append(np.log10(1000*paramsC[0]))
        paraC_ks_n.append(np.log10(paramsC[1]))
        paraC_kp_n.append(np.log10(paramsC[2]))
        paraC_keq_n.append(np.log10(paramsC[3]))
        
        btot_p_n.append(np.log10(pinit[-5]+pinit[-4]))
        ctot_p_n.append(np.log10(pinit[-3]+pinit[-2]))
        
        
        para01_n.append(params01)
        para1E_n.append(params1E)
        paraE2_n.append(paramsE2)
        para23_n.append(params23)
        para35_n.append(params35)
    
        para_24_n.append(params_24)
        para_35_n.append(params_35)
    
        para_t_n.append(params_trans)
        para_diff_3_n.append(params_diff_3)
        para_diff_5_n.append(params_diff_5)
        
        paraB_n.append(paramsB)
        paraC_n.append(paramsC)
    
        ctot_n.append(pinit[-3]+pinit[-2])
        btot_n.append(pinit[-5]+pinit[-4])
        indexes_noswitch.append(i)
       
    #elif  ff_delta[i]<delwant_max and ff_delta[i]>delwant_min and n0delta[i]<0.2 \
     #    and c0delta[i]<0.2:
    # elif ff_delta[i]>0.5 and ff_delta[i]<1.0 and \
    #      c0delta[i]<0.3: 
    elif ff_delta[i]>0.5 and c0delta[i]<0.1 and n0delta[i]<0.1: 
        
        
        para01_v_n2.append(np.log10(1000*params01[0]))
        para01_ks_n2.append(np.log10(params01[1]))
        para01_kp_n2.append(np.log10(params01[2]))
        para01_keq_n2.append(np.log10(params01[3]))
        
        para1E_v_n2.append(np.log10(1000*params1E[0]))
        para1E_ks_n2.append(np.log10(params1E[1]))
        para1E_kp_n2.append(np.log10(params1E[2]))
        para1E_keq_n2.append(np.log10(params1E[3]))
        
        paraE2_v_n2.append(np.log10(1000*paramsE2[0]))
        paraE2_ks_n2.append(np.log10(paramsE2[1]))
        paraE2_kp_n2.append(np.log10(paramsE2[2]))
        paraE2_keq_n2.append(np.log10(paramsE2[3]))
        
        para23_v_n2.append(np.log10(1000*params23[0]))
        para23_ks_n2.append(np.log10(params23[1]))
        # para23_kp_n2.append(np.log10(params23[2]))
        # para23_keq_n2.append(np.log10(params23[3]))
        
        para35_v_n2.append(np.log10(1000*params35[0]))
        para35_ks_n2.append(np.log10(params35[1]))
        para35_kp_n2.append(np.log10(params35[2]))
        para35_keq_n2.append(np.log10(params35[3]))
        
        para_24_v_n2.append(np.log10(1000*params_24[0]))
        para_24_ks_n2.append(np.log10(params_24[1]))
        # para_24_kp_n2.append(np.log10(params_24[2]))
        # para_24_keq_n2.append(np.log10(params_24[3]))
        
        para_35_v_n2.append(np.log10(1000*params_35[0]))
        para_35_ks_n2.append(np.log10(params_35[1]))
        para_35_kp_n2.append(np.log10(params_35[2]))
        para_35_keq_n2.append(np.log10(params_35[3]))
        
        para_t_v_n2.append(np.log10(1000*params_trans[0]))
        para_t_ks_n2.append(np.log10(params_trans[1]))
        
        para_diff_3_in_n2.append(np.log10(params_diff_3[0]))
        para_diff_3_out_n2.append(np.log10(params_diff_3[1]))
        para_diff_5_in_n2.append(np.log10(params_diff_5[0]))
        para_diff_5_out_n2.append(np.log10(params_diff_5[1]))
        
        paraB_v_n2.append(np.log10(1000*paramsB[0]))
        paraB_ks_n2.append(np.log10(paramsB[1]))
        paraB_kp_n2.append(np.log10(paramsB[2]))
        paraB_keq_n2.append(np.log10(paramsB[3]))
        
        paraC_v_n2.append(np.log10(1000*paramsC[0]))
        paraC_ks_n2.append(np.log10(paramsC[1]))
        paraC_kp_n2.append(np.log10(paramsC[2]))
        paraC_keq_n2.append(np.log10(paramsC[3]))
        
        btot_p_n2.append(np.log10(pinit[-5]+pinit[-4]))
        ctot_p_n2.append(np.log10(pinit[-3]+pinit[-2]))
        
        
        para01_n2.append(params01)
        para1E_n2.append(params1E)
        paraE2_n2.append(paramsE2)
        para23_n2.append(params23)
        para35_n2.append(params35)
    
        para_24_n2.append(params_24)
        para_35_n2.append(params_35)
    
        para_t_n2.append(params_trans)
        para_diff_3_n2.append(params_diff_3)
        para_diff_5_n2.append(params_diff_5)
        
        paraB_n2.append(paramsB)
        paraC_n2.append(paramsC)
    
        ctot_n2.append(pinit[-3]+pinit[-2])
        btot_n2.append(pinit[-5]+pinit[-4])
        indexes_switch.append(i)
    else:
        ind_notpicked.append(i)
#[812, 649, 416, 603,682]
  

# allo19_desiredafteroxidase=[11,70,78,85,115,116,129,151,158,202,240,254,255,
#                             328,421,445,453,454,487,502,555,560,564,565,570,
#                             571,602,643,674,681,697,718,739,785,791,874,978,
#                             995]

#choose from [85, 115, 151, 202, 255, 421, 445, 453, 454, 487, 555, 
                #560, 564, 570, 571, 602, 681, 697, 718, 739, 785, 874, 978, 995]
  #    158, 254    ,     328 502
  
# vgood [85, 115, 202, 328, 453, 454, 487, 555, 
                #560, 564, 570, 571, 602, 681, 697, 718, 739, 785, 874, 978, 995]

# indplot=[33]

#202
        
#ax2.scatter(np.abs(ff_delta),np.abs(n0delta),s=40,color='lightpink',alpha=0.5,edgecolors='k',zorder=10)#,label='')
ax2.scatter(ff_delta[indexes_switch],c0delta[indexes_switch],s=40,color='lightblue',alpha=1,edgecolors='k',zorder=10, 
            label='Enzyme parameter mediated swiching')
ax2.scatter(ff_delta[indexes_switch_retained],c0delta[indexes_switch_retained],s=40,color='limegreen',alpha=1,edgecolors='k',zorder=10, 
            label='Co-sub mediated switching')
ax2.scatter(ff_delta[ind_notpicked],c0delta[ind_notpicked],s=40,color='bisque',alpha=0.5,edgecolors='k',zorder=0, 
            label='Other')
ax2.scatter(ff_delta[indexes_noswitch],c0delta[indexes_noswitch],s=40,color='lightcoral',alpha=0.5,edgecolors='k',zorder=20, 
            label='No switching')

# ax2.scatter(ff_delta[indplot_switchretained],n0delta[indplot_switchretained],s=40,color='red',alpha=1,edgecolors='k',zorder=20)
ax2.axis([-0.05,1.05,-0.05,1.18])
ax2.legend(loc='upper left',ncols=2)
ax2.set_xlabel('$\\Delta F_f$',fontsize=30)
ax2.set_ylabel('$\\Delta N_{\\rm H}^{ m}$',fontsize=30,rotation=0)
ax2.yaxis.set_label_coords(-0.09,0.56)
# ax3.text(10**-7,0.8,str(indplot_switchretained),fontsize=20)
# ax2.set_xticks(np.logspace(-9,-2,8), ['$10^{-9}$','','$10^{-7}$','',
#                                      '$10^{-5}$','','$10^{-3}$',''])
# ind_temp=indplot_switchretained
# args_temp=args_to_run[ind_temp]
# print('\n')
# print(ind_temp)
# print(args_temp[0])
# print(args_temp[1])
# print(args_temp[2])
# print(args_temp[3])
# print('\n')
# print(args_temp[4])
# print(args_temp[5])
# print(args_temp[6])
# print(args_temp[7])
# print('\n')
# print(args_temp[8])
# print(args_temp[9])
# print(args_temp[10])
# print(args_temp[11])
# print('\n')
# print(args_temp[12])
# print(args_temp[13])
# print('\n')

# print(args_temp[4])
# print(args_temp[6])



# kinplot=np.logspace(-8,-3,numkin)

# result=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
#                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
#                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
#                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
#                              args_temp[16],args_temp[17])

    
# # flux01= ode_cm.flux_enz_forward(result[:,0], result[:,1], result[:,9], result[:,10], 
# #                          args_temp[3][0], args_temp[3][1], args_temp[3][2])            # gap -> pep

# # flux12= ode_cm.flux_enz_forward(result[:,1], result[:,2], 1.,1., 
# #                          args_temp[4][0], args_temp[4][1], args_temp[4][2])           # pep -> pyr

# flux23= ode_cm.flux_enz_forward(result[:,2]**allo[0], 1.,
#                              args_temp[6][0],  args_temp[6][1])            # pyr -> acet

# # flux35= ode_cm.flux_enz_forward(result[:,3], result[:,4], result[:,10],result[:,9], 
# #                          args_temp[6][0], args_temp[6][1], args_temp[6][2])              # acet -> etn

# # flux22=ode_cm.transp(result[:,2],args_temp[9][0],args_temp[9][1])             # pyr -> pyr_mito

# flux24= ode_cm.flux_enz_forward(result[:,5],  result[:,11], 
#                              args_temp[8][0], args_temp[8][1])

# ax3.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
# ax3.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
# ax3.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
# ax3.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')

# ax3.set_xlabel('$k_{\\rm in}$',fontsize=20)
# # ax3.set_ylabel('flux fraction')
# ax3.legend(loc=(0.02,1.015),ncols=2)
# ax3.set_xscale('log')
# ax3.set_xticks(np.logspace(-8,-3,6), ['$10^{-8}$','','$10^{-6}$','',
#                                      '$10^{-4}$',''])#,'$10^{-3}$',''])

# # args_temp[12][3]*=0.1
# args_temp[13][3]*=0.1
# args_temp[14][3]*=0.1
# # args_temp[12][0]*=10.
# # args_temp[13][0]*=10.

# result2=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
#                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
#                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
#                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
#                              args_temp[16],args_temp[17])

  
# flux23_2= ode_cm.flux_enz_forward(result2[:,2]**allo[0], 1.,
#                              args_temp[6][0],  args_temp[6][1])             # pyr -> acet

# flux24_2= ode_cm.flux_enz_forward(result2[:,5],  result2[:,11], 
#                              args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


# ax3.plot(kinplot,1-(result2[:,9]+result2[:,11])/(result2[:,9]+result2[:,10]+result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='firebrick')
# ax3.plot(kinplot,1-(result2[:,9])/(result2[:,9]+result2[:,10]),linestyle='--',linewidth=2,color='peru')
# ax3.plot(kinplot,1-(result2[:,11])/(result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
# ax3.plot(kinplot,1-flux24_2/(flux23_2+flux24_2),linestyle='--',linewidth=2,color='royalblue')


boxwidth=0.8

bplot=ax5.boxplot([paraC_v , paraC_v_n2,paraC_v_n],widths=(boxwidth,boxwidth,boxwidth),
            patch_artist=True)

colors=['limegreen','lightblue','lightcoral']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

for median in bplot['medians']:
    median.set_color('black')

# ax5.set_title('Total cosubstrate (Mito)', fontsize=20)
ax5.set_title('$V_{\\rm M}^{ N^m}$', fontsize=15,y=1.01)

ax5.set_xticklabels(labels=['C','E','N'], 
                       rotation=0)
ax5.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
ax5.grid('on')



bplot=ax6.boxplot([paraC_ks, paraC_ks_n2,paraC_ks_n], widths=(boxwidth,boxwidth,boxwidth),
            patch_artist=True)

colors=['limegreen','lightblue','lightcoral']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
for median in bplot['medians']:
    median.set_color('black')                 
# ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
ax6.set_title('$K_{\\rm S}^{N_{+}^m}$', fontsize=15,y=1.01)

ax6.set_xticklabels(labels=['C','E','N'], 
                       rotation=0)
ax6.set_yticks(ticks=[-1,0,1,2], labels=['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'])
# ax6.set_yticks(ticks=np.linspace(-4,0))

ax6.grid('on')



bplot=ax7.boxplot([paraC_kp, paraC_kp_n2,paraC_kp_n],widths=(boxwidth,boxwidth,boxwidth),
            patch_artist=True)

colors=['limegreen','lightblue','lightcoral']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
for median in bplot['medians']:
    median.set_color('black')                 
# ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
ax7.set_title('$K_{\\rm P}^{N_{\\rm H}^m}$', fontsize=15,y=1.01)

ax7.set_xticklabels(labels=['C','E','N'], 
                       rotation=0)
ax7.set_yticks(ticks=[-3,-2,-1,0], labels=['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'])

ax7.grid('on')

# # get the data for original params
# ind_temp=indplot_noncosub
# args_temp=args_to_run[ind_temp]
# print('\n')
# print(ind_temp)
# print(args_temp[0])
# print(args_temp[1])
# print(args_temp[2])
# print(args_temp[3])
# print('\n')
# print(args_temp[4])
# print(args_temp[5])
# print(args_temp[6])
# print(args_temp[7])
# print('\n')
# print(args_temp[8])
# print(args_temp[9])
# print(args_temp[10])
# print(args_temp[11])
# print('\n')
# print(args_temp[12])
# print(args_temp[13])
# print('\n')

# print(args_temp[4])
# print(args_temp[6])



# kinplot=np.logspace(-8,-3,numkin)

# result=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
#                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
#                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
#                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
#                              args_temp[16],args_temp[17])

    
# # flux01= ode_cm.flux_enz_forward(result[:,0], result[:,1], result[:,9], result[:,10], 
# #                          args_temp[3][0], args_temp[3][1], args_temp[3][2])            # gap -> pep

# # flux12= ode_cm.flux_enz_forward(result[:,1], result[:,2], 1.,1., 
# #                          args_temp[4][0], args_temp[4][1], args_temp[4][2])           # pep -> pyr

# flux23= ode_cm.flux_enz_forward(result[:,2]**allo[0], 1.,
#                              args_temp[6][0],  args_temp[6][1])            # pyr -> acet

# # flux35= ode_cm.flux_enz_forward(result[:,3], result[:,4], result[:,10],result[:,9], 
# #                          args_temp[6][0], args_temp[6][1], args_temp[6][2])              # acet -> etn

# # flux22=ode_cm.transp(result[:,2],args_temp[9][0],args_temp[9][1])             # pyr -> pyr_mito

# flux24= ode_cm.flux_enz_forward(result[:,5],  result[:,11], 
#                              args_temp[8][0], args_temp[8][1])

# ax4.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
# ax4.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
# ax4.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
# ax4.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')

# ax4.set_xlabel('$k_{\\rm in}$',fontsize=20)
# # ax3.set_ylabel('flux fraction')
# ax4.set_xscale('log')
# ax4.set_xticks(np.logspace(-8,-3,6), ['$10^{-8}$','','$10^{-6}$','',
#                                      '$10^{-4}$',''])#,'$10^{-3}$',''])

# # args_temp[12][3]*=0.1
# args_temp[13][3]*=0.1
# args_temp[14][3]*=0.1
# # args_temp[12][0]*=10.
# # args_temp[13][0]*=10.

# result2=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
#                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
#                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
#                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
#                              args_temp[16],args_temp[17])

  
# flux23_2= ode_cm.flux_enz_forward(result2[:,2]**allo[0], 1.,
#                              args_temp[6][0],  args_temp[6][1])             # pyr -> acet

# flux24_2= ode_cm.flux_enz_forward(result2[:,5],  result2[:,11], 
#                              args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


# ax4.plot(kinplot,1-(result2[:,9]+result2[:,11])/(result2[:,9]+result2[:,10]+result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='firebrick')
# ax4.plot(kinplot,1-(result2[:,9])/(result2[:,9]+result2[:,10]),linestyle='--',linewidth=2,color='peru')
# ax4.plot(kinplot,1-(result2[:,11])/(result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
# ax4.plot(kinplot,1-flux24_2/(flux23_2+flux24_2),linestyle='--',linewidth=2,color='royalblue')



print('param diffs here')
print(10**np.array(np.mean(np.log10(para01_n2),0)-np.mean(np.log10(para01_n),0)))
# print(10**np.array(np.var(np.log10(para01),0)-np.log10(params01var)))
# print(np.mean(para01,0)/params01avg)
print(10**np.array(np.mean(np.log10(para1E_n2),0)-np.mean(np.log10(para1E_n),0)))
print(10**np.array(np.mean(np.log10(paraE2_n2),0)-np.mean(np.log10(paraE2_n),0)))
# print(10**np.array(np.var(np.log10(para12),0)-np.log10(params12var)))
print(10**np.array(np.mean(np.log10(para23_n2),0)-np.mean(np.log10(para23_n),0)))
# print(10**np.array(np.var(np.log10(para23),0)-np.log10(params23var)))
print(10**np.array(np.mean(np.log10(para35_n2),0)-np.mean(np.log10(para35_n),0)))
# print(10**np.array(np.var(np.log10(para35),0)-np.log10(params35var)))
print('\n')
print(10**np.array(np.mean(np.log10(para_24_n2),0)-np.mean(np.log10(para_24_n),0)))
# print(10**np.array(np.var(np.log10(para_24),0)-np.log10(params_24var)))
print(10**np.array(np.mean(np.log10(para_35_n2),0)-np.mean(np.log10(para_35_n),0)))
# print(10**np.array(np.var(np.log10(para_35),0)-np.log10(params_35var)))
print(10**np.array(np.mean(np.log10(para_t_n2),0)-np.mean(np.log10(para_t_n),0)))
# print(10**np.array(np.var(np.log10(para_t),0)-np.log10(params_trans_var)))
print('\n')
print(10**np.array(np.mean(np.log10(paraB_n2),0)-np.mean(np.log10(paraB_n),0)))
# print(10**np.array(np.var(np.log10(paraB),0)-np.log10(paramsBvar)))
print(10**np.array(np.mean(np.log10(paraC_n2),0)-np.mean(np.log10(paraC_n),0)))
# print(10**np.array(np.var(np.log10(paraC),0)-np.log10(paramsCvar)))
print(10**np.array(np.mean(np.log10(ctot_n2),0)-np.mean(np.log10(ctot_n),0)))
# print(10**np.array(np.var(np.log10(ctot),0)-np.log10(ctotvar)))
print(10**np.array(np.mean(np.log10(btot_n2),0)-np.mean(np.log10(btot_n),0)))
# print(10**np.array(np.var(np.log10(btot),0)-np.log10(btotvar)))
# np.savez(savename, flux_cyto=ftlist, flux_mito=fblist, flux_diffn=fdlist,
#          fluxfrac=fflist, fluxfrac_grad=ffgradlist, fluxfrac_sense=ffsenselist,
#          fluxfrac_delta=ffdeltalist, gradmax=gradmaxlist, gradmaxind=gradmaxindlist,
#          sensemax=sensemaxlist, sensemaxind=sensemaxindlist)
print('param diffs end')

print()
boxwidth=0.8

# ax8.boxplot([para23_v_n2 , para23_v_n],widths=(boxwidth,boxwidth))
                  
# # ax5.set_title('Total cosubstrate (Mito)', fontsize=20)
# ax8.set_title('$V_{\\rm M}~ {\\rm Pyr. \\rightarrow Acet.}$', fontsize=15,y=1.01)

# ax8.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# ax8.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
# ax8.grid('on')

# p23vn2=np.array(para23_v_n2)
# ptvn2=np.array(para_t_v_n2)
# p23vn=np.array(para23_v_n)
# ptvn=np.array(para_t_v_n)

# vratios_n2 = p23vn2-ptvn2
# vratios_n= p23vn-ptvn

# ax8.boxplot([vratios_n2 , vratios_n],widths=(boxwidth,boxwidth))
                  
# # ax5.set_title('Total cosubstrate (Mito)', fontsize=20)
# ax8.set_title('$V_{\\rm M,23} / V_{\\rm M,T}  $', fontsize=15,y=1.01)

# ax8.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# # ax8.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
# ax8.grid('on')

# ax9.boxplot([para_t_v_n2, para_t_v_n], widths=(boxwidth,boxwidth))
                  
# # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
# ax9.set_title('$V_{\\rm M}~ {\\rm Pyr. ~Trans.}$', fontsize=15,y=1.01)

# ax9.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# ax9.set_yticks(ticks=[-2,-1], labels=['$10^{-2}$','$10^{-1}$'])
# # ax6.set_yticks(ticks=np.linspace(-4,0))

# ax9.grid('on')

# p23ksn2=np.array(para23_ks_n2)
# ptksn2=np.array(para_t_ks_n2)
# p23ksn=np.array(para23_ks_n)
# ptksn=np.array(para_t_ks_n)

# ksratios_n2 = p23ksn2-ptksn2
# ksratios_n= p23ksn-ptksn

# ax9.boxplot([ksratios_n2, ksratios_n], widths=(boxwidth,boxwidth))
                  
# # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
# ax9.set_title('$K_{\\rm S,T} / K_{\\rm S,T}  $', fontsize=15,y=1.01)

# ax9.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# # ax9.set_yticks(ticks=[-2,-1], labels=['$10^{-2}$','$10^{-1}$'])
# # ax6.set_yticks(ticks=np.linspace(-4,0))

# ax9.grid('on')

# p23vn2=np.array(para23_v_n2)
# ptvn2=np.array(para_t_v_n2)
# p23vn=np.array(para23_v_n)
# ptvn=np.array(para_t_v_n)

# vratios_n2 = p23vn2-ptvn2
# vratios_n= p23vn-ptvn

ptrans_v=np.array(para_t_v)
ptrans_v_n2=np.array(para_t_v_n2)
ptrans_v_n=np.array(para_t_v_n)

p23_v=np.array(para23_v)
p23_v_n2=np.array(para23_v_n2)
p23_v_n=np.array(para23_v_n)

vrat=p23_v-ptrans_v
vrat_n2=p23_v_n2-ptrans_v_n2
vrat_n=p23_v_n-ptrans_v_n

ptrans_ks=np.array(para_t_ks)
ptrans_ks_n2=np.array(para_t_ks_n2)
ptrans_ks_n=np.array(para_t_ks_n)

p23_ks=np.array(para23_ks)
p23_ks_n2=np.array(para23_ks_n2)
p23_ks_n=np.array(para23_ks_n)

ksrat=p23_ks-ptrans_ks
ksrat_n2=p23_ks_n2-ptrans_ks_n2
ksrat_n=p23_ks_n-ptrans_ks_n


bplot=ax10.boxplot([vrat,vrat_n2, vrat_n],widths=(boxwidth,boxwidth,boxwidth),
            patch_artist=True)

colors=['limegreen','lightblue','lightcoral']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
for median in bplot['medians']:
    median.set_color('black')
               
# ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
ax10.set_title('$V_{\\rm M}^{{\\rm Pyr}~ \\rightarrow~ {\\rm AcCHO}} / V_{\\rm M}^{\\rm Pyr. ~Trans.}$', fontsize=15,y=1.01,x=0.6)

ax10.set_xticklabels(labels=['C','E','N'], 
                       rotation=0)
ax10.set_yticks(ticks=[-0,1,2,3], labels=['$10^{0}$','$10^{1}$',\
                                          '$10^{2}$','$10^{3}$'])

ax10.grid('on')


bplot=ax11.boxplot([ksrat,ksrat_n2, ksrat_n],widths=(boxwidth,boxwidth,boxwidth),
            patch_artist=True)

colors=['limegreen','lightblue','lightcoral']
for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)
for median in bplot['medians']:
    median.set_color('black')               
# ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
ax11.set_title('$K_{\\rm S}^{{\\rm Pyr}~ \\rightarrow ~{\\rm AcCHO}} / K_{\\rm S}^{\\rm Pyr. ~Trans.}$', fontsize=15,y=1.01,x=0.6)

ax11.set_xticklabels(labels=['C','E','N'], 
                       rotation=0)
ax11.set_yticks(ticks=[0,1,2], labels=['$10^{0}$',\
                                           '$10^{1}$','$10^{2}$'])

ax11.grid('on')




# bplot=ax12.boxplot([paraB_ks_n2, paraB_ks_n],widths=(boxwidth,boxwidth),
#             patch_artist=True)

# colors=['lightblue','lightcoral']
# for patch, color in zip(bplot['boxes'], colors):
#     patch.set_facecolor(color)
# for median in bplot['medians']:
#     median.set_color('black')            
# # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
# ax12.set_title('$K_{\\rm S}^{\\rm N}~ {\\rm (Cyto)}$', fontsize=15,y=1.01)

# ax12.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# ax12.set_yticks(ticks=[-2,-1,-0,1], labels=['$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$'])

# ax12.grid('on')


# ax14.hist(ff_delta,bins=20)
# ax14.set_xlabel('flux frac delta')


# allo=[1.,1.,1.,1.]
# savepath='../data/allo_'+str(allo)

# filename='fixed_new_3_seed'+seedstr+'_n_'+str(n)
# savename=savepath+filename
# print(savepath)
# print(filename)
# p=Path(savepath)

# filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
# print(filewant)

# kinvals=args_to_run[0][1]

# dat=np.load(filewant)
# ff_delta=dat['fluxfrac_delta']
# ff_delta=np.abs(ff_delta)
# gradmaxind=dat['gradmaxind']
# sensemaxind=dat['sensemaxind']
# fluxfrac=dat['fluxfrac']
# ff_grad=dat['fluxfrac_grad']
# sensemax=dat['sensemax']
# gradmax=dat['gradmax']
# ffsense_vec=dat['fluxfrac_sense']
# b0delta=dat['b0delta']
# c0delta=dat['c0delta']
# n0delta=dat['n0delta']

# ax13.scatter(ff_delta,n0delta,s=40,color='lightpink',alpha=0.5,edgecolors='k',zorder=20)#,label='')


# ax15.hist(ff_delta,bins=20)
# ax15.set_xlabel('flux frac delta')

# print(indexes_switch_retained)

# axs2[0,0].set_xticklabels(labels=['co','No Switch'], 
#                         rotation=90)

if savefig==1:
    plt.savefig('Figure_S6.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S6.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S6.pdf',bbox_inches='tight',format='pdf')
    plt.savefig('Figure_S6.svg',bbox_inches='tight',format='svg')


    
plt.show()