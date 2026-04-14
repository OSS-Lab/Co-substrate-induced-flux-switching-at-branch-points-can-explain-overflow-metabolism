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
sys.path.append('../simulation_functions')
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
# numkin=41


diffvals=[10**-2]#np.logspace(-3,0,4)
diffscales=[2.]

VBmod=1.
VCmod=1.
KBmod=1.
KCmod=1.

params = {'legend.fontsize': 12,
          'axes.labelsize': 15,
          'axes.titlesize': 15,
          'xtick.labelsize':15,
          'ytick.labelsize':15}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)
fs=18

for diffval in diffvals:
    for diffscale in diffscales:
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
        
        # ax2=plt.axes([1,-0.6,0.9,1.15])    # scatter
        # ax3=plt.axes([-0.04,-0.65,0.4,0.48])    # fluxfraction
        # ax4=plt.axes([0.44,-0.65,0.4,0.48])    # expt
        
        # ax3=plt.axes([2,-0.6,0.9,1.1])    # fluxfraction
        # ax4=plt.axes([3,-0.6,0.9,1.1])    # expt
        
        # ax5=plt.axes([0.05,-0.6,0.8,0.5])   # box plots
        
        # ax5.text(0.2,0.45,'Experimental data',fontsize=30)
        
        
        ax16=plt.axes([2.0,-3.5,0.5,0.48])   # box plots
        ax17=plt.axes([2.6,-3.5,0.5,0.48])
        ax18=plt.axes([3.2,-3.5,0.5,0.48])
        ax19=plt.axes([3.8,-3.5,0.5,0.48])
        # ax199=plt.axes([3.8,-3.5,0.5,0.48])
        # plt.axis('off')
        # ax199.text(1.1,0.5,'${\\rm NAD}^{\\rm +}_{\\rm C} \\rightarrow {\\rm NADH}_{\\rm C} $',fontsize=30)
        ax16.set_title('$V_{\\rm max}^{{\\rm NAD}^{\\rm +}_{\\rm C} \\rightarrow {\\rm NADH}_{\\rm C}}$',fontsize=fs, y=1.05)#, pad=-14)
        ax17.set_title('$ K_{\\rm S}^{{\\rm NAD}^{\\rm +}_{\\rm C} \\rightarrow {\\rm NADH}_{\\rm C}}$',fontsize=fs,y=1.05)
        ax18.set_title('$ K_{\\rm P}^{{\\rm NAD}^{\\rm +}_{\\rm C} \\rightarrow {\\rm NADH}_{\\rm C}}$',fontsize=fs,y=1.05)
        ax19.set_title('$N_{\\rm C}^{\\rm Tot}$',fontsize=fs,y=1.05)
        
        
        ax20=plt.axes([2.0,-4.2,0.5,0.48])   # box plots
        ax21=plt.axes([2.6,-4.2,0.5,0.48])
        ax22=plt.axes([3.2,-4.2,0.5,0.48])
        ax23=plt.axes([3.8,-4.2,0.5,0.48])
        # ax233=plt.axes([3.8,-4.2,0.5,0.48])
        # plt.axis('off')
        ax20.set_title('$V_{\\rm max}^{{\\rm NAD}^{\\rm +}_{\\rm M} \\rightarrow {\\rm NADH}_{\\rm M}}$',fontsize=fs, y=1.05)#, pad=-14)
        ax21.set_title('$ K_{\\rm S}^{{\\rm NAD}^{\\rm +}_{\\rm M} \\rightarrow {\\rm NADH}_{\\rm M}}$',fontsize=fs,y=1.05)
        ax22.set_title('$ K_{\\rm P}^{{\\rm NAD}^{\\rm +}_{\\rm M} \\rightarrow {\\rm NADH}_{\\rm M}}$',fontsize=fs,y=1.05)
        ax23.set_title('$N_{\\rm M}^{\\rm Tot}$',fontsize=fs,y=1.05)
        # ax233.text(1.1,0.5,'${\\rm NAD}^{\\rm +}_{\\rm M} \\rightarrow {\\rm NADH}_{\\rm M} $',fontsize=30)
        
        ax31=plt.axes([2.0,-4.9,0.5,0.48])
        ax311=plt.axes([2.6,-4.9,0.5,0.48])
        # ax3111=plt.axes([2.6,-4.9,0.5,0.48])
        # plt.axis('off')
        ax31.set_title('$V_{\\rm max}^{{\\rm Pyr}~ \\rightarrow~ {\\rm AcCHO}} / V_{\\rm max}^{\\rm Pyr. ~Trans.} $',fontsize=fs, y=1.05)#, pad=-14)
        ax311.set_title('$ K_{\\rm M}^{{\\rm Pyr}~ \\rightarrow ~{\\rm AcCHO}} / K_{\\rm M}^{\\rm Pyr. ~Trans.}$',fontsize=fs,y=1.05)
        
        ax32=plt.axes([3.2,-4.9,0.5,0.48])
        ax322=plt.axes([3.8,-4.9,0.5,0.48])
        # ax3222=plt.axes([3.8,-4.9,0.5,0.48])
        # plt.axis('off')
        ax32.set_title('$V_{\\rm max}^{{\\rm AcCHO_C}~ \\rightarrow~ {\\rm EtOH_C}} / V_{\\rm max}^{{\\rm AcCHO_M}~ \\rightarrow~ {\\rm EtOH_M}}$',fontsize=fs,y=1.05)
        ax322.set_title('$K_{\\rm M}^{{\\rm AcCHO_C}~ \\rightarrow~ {\\rm EtOH_C}} / K_{\\rm M}^{{\\rm AcCHO_M}~ \\rightarrow~ {\\rm EtOH_M}}$',fontsize=fs,y=1.05,x=0.6)
        
        # ax5=plt.axes([2.0,0.02,0.3,0.48])   # box plots
        # ax6=plt.axes([2.4,0.02,0.3,0.48])
        # ax7=plt.axes([2.8,0.02,0.3,0.48])
        
        # ax10=plt.axes([2.18,-0.67,0.3,0.48])   # box plots
        # ax11=plt.axes([2.62,-0.67,0.3,0.48])
        # ax12=plt.axes([2.65,-0.67,0.2,0.48])
        # ax8=plt.axes([2.975,-0.67,0.2,0.48])
        # ax9=plt.axes([3.3,-0.67,0.2,0.48])
        
        # ax13=plt.axes([0,-2,0.9,1.15])    # scatter 2
        
        # ax14=plt.axes([1.05,-1.35,0.85,0.5])
        # ax15=plt.axes([1.05,-2,0.85,0.5])
        # ax16=plt.axes([2.05,-1.35,0.85,0.5])
        # ax17=plt.axes([2.05,-2,0.85,0.5])
        # ax2.text(-0.2,0.9,'$\\rm A$',fontsize=40)
        # ax2.text(1.09,0.9,'$\\rm B$',fontsize=40)
        # ax2.text(1.22,0.3,'$\\rm C$',fontsize=40)
        
        # get the data for new params
        allo=[1.9,1.,1.,1.]
        n=10000
        
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
        
        args_to_run=mh_cm.make_params_nov25_allranges_nonfixed(allo,n,VBmod=1.,VCmod=1.,KBmod=1.,KCmod=1.,diffval=1.,diffscale=2.)
        
              
        savepath='../data/allo_'+str(allo)
        
        filename='nonfixed_apr_final_seed'+seedstr+'_n_'+str(n)+\
            '_diffvalscale_'+str(np.log10(diffval))+'_'+str(diffscale)+\
                '_datetime'
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
        
        
        ff_delta=-ff_delta
        # c0delta=np.abs(c0delta)
        # b0delta=np.abs(b0delta)
        # n0delta=np.abs(n0delta)
        
        maxdelta=np.maximum(c0delta,b0delta)
        
        ####################
        para01_ns=list()
        para1E_ns=list()
        paraE2_ns=list()
        para23_ns=list()
        para35_ns=list()
        para_24_ns=list()
        para_35_ns=list()
        
        para_t_ns=list()
        para_diff_3_ns=list()
        para_diff_5_ns=list()
        
        paraB_ns=list()
        paraC_ns=list()
        
        btot_ns=list()
        ctot_ns=list()
        
        para01_v_ns=list()
        para01_ks_ns=list()
        para01_kp_ns=list()
        para01_keq_ns=list()
        
        para1E_v_ns=list()
        para1E_ks_ns=list()
        para1E_kp_ns=list()
        para1E_keq_ns=list()
        
        
        paraE2_v_ns=list()
        paraE2_ks_ns=list()
        paraE2_kp_ns=list()
        paraE2_keq_ns=list()
        
        
        para23_v_ns=list()
        para23_ks_ns=list()
        # para23_kp_ns=list()
        # para23_keq_ns=list()
        
        para35_v_ns=list()
        para35_ks_ns=list()
        para35_kp_ns=list()
        para35_keq_ns=list()
        
        para_24_v_ns=list()
        para_24_ks_ns=list()
        # para_24_kp_ns=list()
        # para_24_keq_ns=list()
        
        para_35_v_ns=list()
        para_35_ks_ns=list()
        para_35_kp_ns=list()
        para_35_keq_ns=list()
        
        para_t_v_ns=list()
        para_t_ks_ns=list()
        
        para_diff_3_in_ns=list()
        para_diff_3_out_ns=list()
        para_diff_5_in_ns=list()
        para_diff_5_out_ns=list()
        
        paraB_v_ns=list()
        paraB_ks_ns=list()
        paraB_kp_ns=list()
        paraB_keq_ns=list()
        
        paraC_v_ns=list()
        paraC_ks_ns=list()
        paraC_kp_ns=list()
        paraC_keq_ns=list()
        
        btot_p_ns=list()
        ctot_p_ns=list()
        #######################
        
        #######################
        para01_nc=list()
        para1E_nc=list()
        paraE2_nc=list()
        para23_nc=list()
        para35_nc=list()
        para_24_nc=list()
        para_35_nc=list()
        
        para_t_nc=list()
        para_diff_3_nc=list()
        para_diff_5_nc=list()
        
        paraB_nc=list()
        paraC_nc=list()
        
        btot_nc=list()
        ctot_nc=list()
        
        para01_v_nc=list()
        para01_ks_nc=list()
        para01_kp_nc=list()
        para01_keq_nc=list()
        
        para1E_v_nc=list()
        para1E_ks_nc=list()
        para1E_kp_nc=list()
        para1E_keq_nc=list()
        
        
        paraE2_v_nc=list()
        paraE2_ks_nc=list()
        paraE2_kp_nc=list()
        paraE2_keq_nc=list()
        
        
        para23_v_nc=list()
        para23_ks_nc=list()
        # para23_kp_nc=list()
        # para23_keq_nc=list()
        
        para35_v_nc=list()
        para35_ks_nc=list()
        para35_kp_nc=list()
        para35_keq_nc=list()
        
        para_24_v_nc=list()
        para_24_ks_nc=list()
        # para_24_kp_nc=list()
        # para_24_keq_nc=list()
        
        para_35_v_nc=list()
        para_35_ks_nc=list()
        para_35_kp_nc=list()
        para_35_keq_nc=list()
        
        para_t_v_nc=list()
        para_t_ks_nc=list()
        
        para_diff_3_in_nc=list()
        para_diff_3_out_nc=list()
        para_diff_5_in_nc=list()
        para_diff_5_out_nc=list()
        
        paraB_v_nc=list()
        paraB_ks_nc=list()
        paraB_kp_nc=list()
        paraB_keq_nc=list()
        
        paraC_v_nc=list()
        paraC_ks_nc=list()
        paraC_kp_nc=list()
        paraC_keq_nc=list()
        
        btot_p_nc=list()
        ctot_p_nc=list()
        ###########################
        
        #######################
        para01_b=list()
        para1E_b=list()
        paraE2_b=list()
        para23_b=list()
        para35_b=list()
        para_24_b=list()
        para_35_b=list()
        
        para_t_b=list()
        para_diff_3_b=list()
        para_diff_5_b=list()
        
        paraB_b=list()
        paraC_b=list()
        
        btot_b=list()
        ctot_b=list()
        
        para01_v_b=list()
        para01_ks_b=list()
        para01_kp_b=list()
        para01_keq_b=list()
        
        para1E_v_b=list()
        para1E_ks_b=list()
        para1E_kp_b=list()
        para1E_keq_b=list()
        
        
        paraE2_v_b=list()
        paraE2_ks_b=list()
        paraE2_kp_b=list()
        paraE2_keq_b=list()
        
        
        para23_v_b=list()
        para23_ks_b=list()
        # para23_kp_b=list()
        # para23_keq_b=list()
        
        para35_v_b=list()
        para35_ks_b=list()
        para35_kp_b=list()
        para35_keq_b=list()
        
        para_24_v_b=list()
        para_24_ks_b=list()
        # para_24_kp_b=list()
        # para_24_keq_b=list()
        
        para_35_v_b=list()
        para_35_ks_b=list()
        para_35_kp_b=list()
        para_35_keq_b=list()
        
        para_t_v_b=list()
        para_t_ks_b=list()
        
        para_diff_3_in_b=list()
        para_diff_3_out_b=list()
        para_diff_5_in_b=list()
        para_diff_5_out_b=list()
        
        paraB_v_b=list()
        paraB_ks_b=list()
        paraB_kp_b=list()
        paraB_keq_b=list()
        
        paraC_v_b=list()
        paraC_ks_b=list()
        paraC_kp_b=list()
        paraC_keq_b=list()
        
        btot_p_b=list()
        ctot_p_b=list()
        ###########################
        
        
        ###########################
        para01_c=list()
        para1E_c=list()
        paraE2_c=list()
        para23_c=list()
        para35_c=list()
        para_24_c=list()
        para_35_c=list()
        
        para_t_c=list()
        para_diff_3_c=list()
        para_diff_5_c=list()
        
        paraB_c=list()
        paraC_c=list()
        
        btot_c=list()
        ctot_c=list()
        
        para01_v_c=list()
        para01_ks_c=list()
        para01_kp_c=list()
        para01_keq_c=list()
        
        para1E_v_c=list()
        para1E_ks_c=list()
        para1E_kp_c=list()
        para1E_keq_c=list()
        
        
        paraE2_v_c=list()
        paraE2_ks_c=list()
        paraE2_kp_c=list()
        paraE2_keq_c=list()
        
        
        para23_v_c=list()
        para23_ks_c=list()
        # para23_kp_c=list()
        # para23_keq_c=list()
        
        para35_v_c=list()
        para35_ks_c=list()
        para35_kp_c=list()
        para35_keq_c=list()
        
        para_24_v_c=list()
        para_24_ks_c=list()
        # para_24_kp_c=list()
        # para_24_keq_c=list()
        
        para_35_v_c=list()
        para_35_ks_c=list()
        para_35_kp_c=list()
        para_35_keq_c=list()
        
        para_t_v_c=list()
        para_t_ks_c=list()
        
        para_diff_3_in_c=list()
        para_diff_3_out_c=list()
        para_diff_5_in_c=list()
        para_diff_5_out_c=list()
        
        paraB_v_c=list()
        paraB_ks_c=list()
        paraB_kp_c=list()
        paraB_keq_c=list()
        
        paraC_v_c=list()
        paraC_ks_c=list()
        paraC_kp_c=list()
        paraC_keq_c=list()
        
        btot_p_c=list()
        ctot_p_c=list()
        #######################
        
        #######################
        para01_m=list()
        para1E_m=list()
        paraE2_m=list()
        para23_m=list()
        para35_m=list()
        para_24_m=list()
        para_35_m=list()
        
        para_t_m=list()
        para_diff_3_m=list()
        para_diff_5_m=list()
        
        paraB_m=list()
        paraC_m=list()
        
        btot_m=list()
        ctot_m=list()
        
        para01_v_m=list()
        para01_ks_m=list()
        para01_kp_m=list()
        para01_keq_m=list()
        
        para1E_v_m=list()
        para1E_ks_m=list()
        para1E_kp_m=list()
        para1E_keq_m=list()
        
        
        paraE2_v_m=list()
        paraE2_ks_m=list()
        paraE2_kp_m=list()
        paraE2_keq_m=list()
        
        
        para23_v_m=list()
        para23_ks_m=list()
        # para23_kp_m=list()
        # para23_keq_m=list()
        
        para35_v_m=list()
        para35_ks_m=list()
        para35_kp_m=list()
        para35_keq_m=list()
        
        para_24_v_m=list()
        para_24_ks_m=list()
        # para_24_kp_m=list()
        # para_24_keq_m=list()
        
        para_35_v_m=list()
        para_35_ks_m=list()
        para_35_kp_m=list()
        para_35_keq_m=list()
        
        para_t_v_m=list()
        para_t_ks_m=list()
        
        para_diff_3_in_m=list()
        para_diff_3_out_m=list()
        para_diff_5_in_m=list()
        para_diff_5_out_m=list()
        
        paraB_v_m=list()
        paraB_ks_m=list()
        paraB_kp_m=list()
        paraB_keq_m=list()
        
        paraC_v_m=list()
        paraC_ks_m=list()
        paraC_kp_m=list()
        paraC_keq_m=list()
        
        btot_p_m=list()
        ctot_p_m=list()
        #######################
        
        
       
        
        
        
        indexes_noswitching=list()      # no switching in ff or cosub
        indexes_switch_nocosub=list()   # switching in ff but none in cosub
        indexes_switch_cyto=list()      # switching in ff and cyto
        indexes_switch_mito=list()      # switching in ff and mito
        indexes_switch_both=list()      # switching in ff and both
        indexes_uncat=list()            # others that dont fit these
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
            # print(ff_delta[i])
            if ff_delta[i]<0.4: 
            # if i in allo19_desiredafteroxidase and n0delta[i]>0.2:        
                para01_ns.append(params01)
                para1E_ns.append(params1E)
                paraE2_ns.append(paramsE2)
                para23_ns.append(params23)
                para35_ns.append(params35)
            
                para_24_ns.append(params_24)
                para_35_ns.append(params_35)
            
                para_t_ns.append(params_trans)
                para_diff_3_ns.append(params_diff_3)
                para_diff_5_ns.append(params_diff_5)
                
                paraB_ns.append(paramsB)
                paraC_ns.append(paramsC)
            
                ctot_ns.append(pinit[-3]+pinit[-2])
                btot_ns.append(pinit[-5]+pinit[-4])
                
                
                para01_v_ns.append(np.log10(1000*params01[0]))
                para01_ks_ns.append(np.log10(params01[1]))
                para01_kp_ns.append(np.log10(params01[2]))
                para01_keq_ns.append(np.log10(params01[3]))
                
                para1E_v_ns.append(np.log10(1000*params1E[0]))
                para1E_ks_ns.append(np.log10(params1E[1]))
                para1E_kp_ns.append(np.log10(params1E[2]))
                para1E_keq_ns.append(np.log10(params1E[3]))
                
                paraE2_v_ns.append(np.log10(1000*paramsE2[0]))
                paraE2_ks_ns.append(np.log10(paramsE2[1]))
                paraE2_kp_ns.append(np.log10(paramsE2[2]))
                paraE2_keq_ns.append(np.log10(paramsE2[3]))
                
                para23_v_ns.append(np.log10(1000*params23[0]))
                para23_ks_ns.append(np.log10(params23[1]))
                # para23_kp.append(np.log10(params23[2]))
                # para23_keq.append(np.log10(params23[3]))
                
                para35_v_ns.append(np.log10(1000*params35[0]))
                para35_ks_ns.append(np.log10(params35[1]))
                para35_kp_ns.append(np.log10(params35[2]))
                para35_keq_ns.append(np.log10(params35[3]))
                
                para_24_v_ns.append(np.log10(1000*params_24[0]))
                para_24_ks_ns.append(np.log10(params_24[1]))
                # para_24_kp.append(np.log10(params_24[2]))
                # para_24_keq.append(np.log10(params_24[3]))
                
                para_35_v_ns.append(np.log10(1000*params_35[0]))
                para_35_ks_ns.append(np.log10(params_35[1]))
                para_35_kp_ns.append(np.log10(params_35[2]))
                para_35_keq_ns.append(np.log10(params_35[3]))
                
                para_t_v_ns.append(np.log10(1000*params_trans[0]))
                para_t_ks_ns.append(np.log10(params_trans[1]))
                
                para_diff_3_in_ns.append(np.log10(params_diff_3[0]))
                para_diff_3_out_ns.append(np.log10(params_diff_3[1]))
                para_diff_5_in_ns.append(np.log10(params_diff_5[0]))
                para_diff_5_out_ns.append(np.log10(params_diff_5[1]))
                
                paraB_v_ns.append(np.log10(1000*paramsB[0]))
                paraB_ks_ns.append(np.log10(paramsB[1]))
                paraB_kp_ns.append(np.log10(paramsB[2]))
                paraB_keq_ns.append(np.log10(paramsB[3]))
                
                paraC_v_ns.append(np.log10(1000*paramsC[0]))
                paraC_ks_ns.append(np.log10(paramsC[1]))
                paraC_kp_ns.append(np.log10(paramsC[2]))
                paraC_keq_ns.append(np.log10(paramsC[3]))
                
                btot_p_ns.append(np.log10(pinit[-5]+pinit[-4]))
                ctot_p_ns.append(np.log10(pinit[-3]+pinit[-2]))
            
        
                indexes_noswitching.append(i)
            
            elif ff_delta[i]>=0.4 and maxdelta[i]<0.2: 
                    
                para01_nc.append(params01)
                para1E_nc.append(params1E)
                paraE2_nc.append(paramsE2)
                para23_nc.append(params23)
                para35_nc.append(params35)
            
                para_24_nc.append(params_24)
                para_35_nc.append(params_35)
            
                para_t_nc.append(params_trans)
                para_diff_3_nc.append(params_diff_3)
                para_diff_5_nc.append(params_diff_5)
                
                paraB_nc.append(paramsB)
                paraC_nc.append(paramsC)
            
                ctot_nc.append(pinit[-3]+pinit[-2])
                btot_nc.append(pinit[-5]+pinit[-4])
                
                
                para01_v_nc.append(np.log10(1000*params01[0]))
                para01_ks_nc.append(np.log10(params01[1]))
                para01_kp_nc.append(np.log10(params01[2]))
                para01_keq_nc.append(np.log10(params01[3]))
                
                para1E_v_nc.append(np.log10(1000*params1E[0]))
                para1E_ks_nc.append(np.log10(params1E[1]))
                para1E_kp_nc.append(np.log10(params1E[2]))
                para1E_keq_nc.append(np.log10(params1E[3]))
                
                paraE2_v_nc.append(np.log10(1000*paramsE2[0]))
                paraE2_ks_nc.append(np.log10(paramsE2[1]))
                paraE2_kp_nc.append(np.log10(paramsE2[2]))
                paraE2_keq_nc.append(np.log10(paramsE2[3]))
                
                para23_v_nc.append(np.log10(1000*params23[0]))
                para23_ks_nc.append(np.log10(params23[1]))
                # para23_kp.append(np.log10(params23[2]))
                # para23_keq.append(np.log10(params23[3]))
                
                para35_v_nc.append(np.log10(1000*params35[0]))
                para35_ks_nc.append(np.log10(params35[1]))
                para35_kp_nc.append(np.log10(params35[2]))
                para35_keq_nc.append(np.log10(params35[3]))
                
                para_24_v_nc.append(np.log10(1000*params_24[0]))
                para_24_ks_nc.append(np.log10(params_24[1]))
                # para_24_kp.append(np.log10(params_24[2]))
                # para_24_keq.append(np.log10(params_24[3]))
                
                para_35_v_nc.append(np.log10(1000*params_35[0]))
                para_35_ks_nc.append(np.log10(params_35[1]))
                para_35_kp_nc.append(np.log10(params_35[2]))
                para_35_keq_nc.append(np.log10(params_35[3]))
                
                para_t_v_nc.append(np.log10(1000*params_trans[0]))
                para_t_ks_nc.append(np.log10(params_trans[1]))
                
                para_diff_3_in_nc.append(np.log10(params_diff_3[0]))
                para_diff_3_out_nc.append(np.log10(params_diff_3[1]))
                para_diff_5_in_nc.append(np.log10(params_diff_5[0]))
                para_diff_5_out_nc.append(np.log10(params_diff_5[1]))
                
                paraB_v_nc.append(np.log10(1000*paramsB[0]))
                paraB_ks_nc.append(np.log10(paramsB[1]))
                paraB_kp_nc.append(np.log10(paramsB[2]))
                paraB_keq_nc.append(np.log10(paramsB[3]))
                
                paraC_v_nc.append(np.log10(1000*paramsC[0]))
                paraC_ks_nc.append(np.log10(paramsC[1]))
                paraC_kp_nc.append(np.log10(paramsC[2]))
                paraC_keq_nc.append(np.log10(paramsC[3]))
                
                btot_p_nc.append(np.log10(pinit[-5]+pinit[-4]))
                ctot_p_nc.append(np.log10(pinit[-3]+pinit[-2]))
            
        
                indexes_switch_nocosub.append(i)
                
            elif ff_delta[i]>=0.4 and b0delta[i]>=0.2 and c0delta[i]>=0.2: 
                    
                para01_b.append(params01)
                para1E_b.append(params1E)
                paraE2_b.append(paramsE2)
                para23_b.append(params23)
                para35_b.append(params35)
            
                para_24_b.append(params_24)
                para_35_b.append(params_35)
            
                para_t_b.append(params_trans)
                para_diff_3_b.append(params_diff_3)
                para_diff_5_b.append(params_diff_5)
                
                paraB_b.append(paramsB)
                paraC_b.append(paramsC)
            
                ctot_b.append(pinit[-3]+pinit[-2])
                btot_b.append(pinit[-5]+pinit[-4])
                
                
                para01_v_b.append(np.log10(1000*params01[0]))
                para01_ks_b.append(np.log10(params01[1]))
                para01_kp_b.append(np.log10(params01[2]))
                para01_keq_b.append(np.log10(params01[3]))
                
                para1E_v_b.append(np.log10(1000*params1E[0]))
                para1E_ks_b.append(np.log10(params1E[1]))
                para1E_kp_b.append(np.log10(params1E[2]))
                para1E_keq_b.append(np.log10(params1E[3]))
                
                paraE2_v_b.append(np.log10(1000*paramsE2[0]))
                paraE2_ks_b.append(np.log10(paramsE2[1]))
                paraE2_kp_b.append(np.log10(paramsE2[2]))
                paraE2_keq_b.append(np.log10(paramsE2[3]))
                
                para23_v_b.append(np.log10(1000*params23[0]))
                para23_ks_b.append(np.log10(params23[1]))
                # para23_kp.append(np.log10(params23[2]))
                # para23_keq.append(np.log10(params23[3]))
                
                para35_v_b.append(np.log10(1000*params35[0]))
                para35_ks_b.append(np.log10(params35[1]))
                para35_kp_b.append(np.log10(params35[2]))
                para35_keq_b.append(np.log10(params35[3]))
                
                para_24_v_b.append(np.log10(1000*params_24[0]))
                para_24_ks_b.append(np.log10(params_24[1]))
                # para_24_kp.append(np.log10(params_24[2]))
                # para_24_keq.append(np.log10(params_24[3]))
                
                para_35_v_b.append(np.log10(1000*params_35[0]))
                para_35_ks_b.append(np.log10(params_35[1]))
                para_35_kp_b.append(np.log10(params_35[2]))
                para_35_keq_b.append(np.log10(params_35[3]))
                
                para_t_v_b.append(np.log10(1000*params_trans[0]))
                para_t_ks_b.append(np.log10(params_trans[1]))
                
                para_diff_3_in_b.append(np.log10(params_diff_3[0]))
                para_diff_3_out_b.append(np.log10(params_diff_3[1]))
                para_diff_5_in_b.append(np.log10(params_diff_5[0]))
                para_diff_5_out_b.append(np.log10(params_diff_5[1]))
                
                paraB_v_b.append(np.log10(1000*paramsB[0]))
                paraB_ks_b.append(np.log10(paramsB[1]))
                paraB_kp_b.append(np.log10(paramsB[2]))
                paraB_keq_b.append(np.log10(paramsB[3]))
                
                paraC_v_b.append(np.log10(1000*paramsC[0]))
                paraC_ks_b.append(np.log10(paramsC[1]))
                paraC_kp_b.append(np.log10(paramsC[2]))
                paraC_keq_b.append(np.log10(paramsC[3]))
                
                btot_p_b.append(np.log10(pinit[-5]+pinit[-4]))
                ctot_p_b.append(np.log10(pinit[-3]+pinit[-2]))
            
        
                indexes_switch_both.append(i)
                
            elif ff_delta[i]>=0.4 and b0delta[i]<0.2 and c0delta[i]>=0.2:
                
                para01_m.append(params01)
                para1E_m.append(params1E)
                paraE2_m.append(paramsE2)
                para23_m.append(params23)
                para35_m.append(params35)
            
                para_24_m.append(params_24)
                para_35_m.append(params_35)
            
                para_t_m.append(params_trans)
                para_diff_3_m.append(params_diff_3)
                para_diff_5_m.append(params_diff_5)
                
                paraB_m.append(paramsB)
                paraC_m.append(paramsC)
            
                ctot_m.append(pinit[-3]+pinit[-2])
                btot_m.append(pinit[-5]+pinit[-4])
                
                
                para01_v_m.append(np.log10(1000*params01[0]))
                para01_ks_m.append(np.log10(params01[1]))
                para01_kp_m.append(np.log10(params01[2]))
                para01_keq_m.append(np.log10(params01[3]))
                
                para1E_v_m.append(np.log10(1000*params1E[0]))
                para1E_ks_m.append(np.log10(params1E[1]))
                para1E_kp_m.append(np.log10(params1E[2]))
                para1E_keq_m.append(np.log10(params1E[3]))
                
                paraE2_v_m.append(np.log10(1000*paramsE2[0]))
                paraE2_ks_m.append(np.log10(paramsE2[1]))
                paraE2_kp_m.append(np.log10(paramsE2[2]))
                paraE2_keq_m.append(np.log10(paramsE2[3]))
                
                para23_v_m.append(np.log10(1000*params23[0]))
                para23_ks_m.append(np.log10(params23[1]))
                # para23_kp.append(np.log10(params23[2]))
                # para23_keq.append(np.log10(params23[3]))
                
                para35_v_m.append(np.log10(1000*params35[0]))
                para35_ks_m.append(np.log10(params35[1]))
                para35_kp_m.append(np.log10(params35[2]))
                para35_keq_m.append(np.log10(params35[3]))
                
                para_24_v_m.append(np.log10(1000*params_24[0]))
                para_24_ks_m.append(np.log10(params_24[1]))
                # para_24_kp.append(np.log10(params_24[2]))
                # para_24_keq.append(np.log10(params_24[3]))
                
                para_35_v_m.append(np.log10(1000*params_35[0]))
                para_35_ks_m.append(np.log10(params_35[1]))
                para_35_kp_m.append(np.log10(params_35[2]))
                para_35_keq_m.append(np.log10(params_35[3]))
                
                para_t_v_m.append(np.log10(1000*params_trans[0]))
                para_t_ks_m.append(np.log10(params_trans[1]))
                
                para_diff_3_in_m.append(np.log10(params_diff_3[0]))
                para_diff_3_out_m.append(np.log10(params_diff_3[1]))
                para_diff_5_in_m.append(np.log10(params_diff_5[0]))
                para_diff_5_out_m.append(np.log10(params_diff_5[1]))
                
                paraB_v_m.append(np.log10(1000*paramsB[0]))
                paraB_ks_m.append(np.log10(paramsB[1]))
                paraB_kp_m.append(np.log10(paramsB[2]))
                paraB_keq_m.append(np.log10(paramsB[3]))
                
                paraC_v_m.append(np.log10(1000*paramsC[0]))
                paraC_ks_m.append(np.log10(paramsC[1]))
                paraC_kp_m.append(np.log10(paramsC[2]))
                paraC_keq_m.append(np.log10(paramsC[3]))
                
                btot_p_m.append(np.log10(pinit[-5]+pinit[-4]))
                ctot_p_m.append(np.log10(pinit[-3]+pinit[-2]))
            
        
                indexes_switch_mito.append(i)
               
            elif ff_delta[i]>=0.4 and b0delta[i]>=0.2 and c0delta[i]<0.2: 
                
                para01_c.append(params01)
                para1E_c.append(params1E)
                paraE2_c.append(paramsE2)
                para23_c.append(params23)
                para35_c.append(params35)
            
                para_24_c.append(params_24)
                para_35_c.append(params_35)
            
                para_t_c.append(params_trans)
                para_diff_3_c.append(params_diff_3)
                para_diff_5_c.append(params_diff_5)
                
                paraB_c.append(paramsB)
                paraC_c.append(paramsC)
            
                ctot_c.append(pinit[-3]+pinit[-2])
                btot_c.append(pinit[-5]+pinit[-4])
                
                
                para01_v_c.append(np.log10(1000*params01[0]))
                para01_ks_c.append(np.log10(params01[1]))
                para01_kp_c.append(np.log10(params01[2]))
                para01_keq_c.append(np.log10(params01[3]))
                
                para1E_v_c.append(np.log10(1000*params1E[0]))
                para1E_ks_c.append(np.log10(params1E[1]))
                para1E_kp_c.append(np.log10(params1E[2]))
                para1E_keq_c.append(np.log10(params1E[3]))
                
                paraE2_v_c.append(np.log10(1000*paramsE2[0]))
                paraE2_ks_c.append(np.log10(paramsE2[1]))
                paraE2_kp_c.append(np.log10(paramsE2[2]))
                paraE2_keq_c.append(np.log10(paramsE2[3]))
                
                para23_v_c.append(np.log10(1000*params23[0]))
                para23_ks_c.append(np.log10(params23[1]))
                # para23_kp.append(np.log10(params23[2]))
                # para23_keq.append(np.log10(params23[3]))
                
                para35_v_c.append(np.log10(1000*params35[0]))
                para35_ks_c.append(np.log10(params35[1]))
                para35_kp_c.append(np.log10(params35[2]))
                para35_keq_c.append(np.log10(params35[3]))
                
                para_24_v_c.append(np.log10(1000*params_24[0]))
                para_24_ks_c.append(np.log10(params_24[1]))
                # para_24_kp.append(np.log10(params_24[2]))
                # para_24_keq.append(np.log10(params_24[3]))
                
                para_35_v_c.append(np.log10(1000*params_35[0]))
                para_35_ks_c.append(np.log10(params_35[1]))
                para_35_kp_c.append(np.log10(params_35[2]))
                para_35_keq_c.append(np.log10(params_35[3]))
                
                para_t_v_c.append(np.log10(1000*params_trans[0]))
                para_t_ks_c.append(np.log10(params_trans[1]))
                
                para_diff_3_in_c.append(np.log10(params_diff_3[0]))
                para_diff_3_out_c.append(np.log10(params_diff_3[1]))
                para_diff_5_in_c.append(np.log10(params_diff_5[0]))
                para_diff_5_out_c.append(np.log10(params_diff_5[1]))
                
                paraB_v_c.append(np.log10(1000*paramsB[0]))
                paraB_ks_c.append(np.log10(paramsB[1]))
                paraB_kp_c.append(np.log10(paramsB[2]))
                paraB_keq_c.append(np.log10(paramsB[3]))
                
                paraC_v_c.append(np.log10(1000*paramsC[0]))
                paraC_ks_c.append(np.log10(paramsC[1]))
                paraC_kp_c.append(np.log10(paramsC[2]))
                paraC_keq_c.append(np.log10(paramsC[3]))
                
                btot_p_c.append(np.log10(pinit[-5]+pinit[-4]))
                ctot_p_c.append(np.log10(pinit[-3]+pinit[-2]))
            
        
                indexes_switch_cyto.append(i)
                
            else:
                indexes_uncat.append(i)
                
        boxwidth=0.8 
        colors=['grey','lightblue','firebrick','mediumseagreen','peru']

        bplot=ax16.violinplot([paraB_v_ns , paraB_v_nc,paraB_v_b,paraB_v_m,paraB_v_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax16.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax16.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        # ax6.set_yticks(ticks=[-1,0,1,2], labels=['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'])

        ax16.grid('on')
        
        bplot=ax17.violinplot([paraB_ks_ns , paraB_ks_nc,paraB_ks_b,paraB_ks_m,paraB_ks_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax17.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax2.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax17.set_yticks(ticks=[-1,-0, 1], labels=['$10^{-1}$','$10^{0}$','$10^{1}$'])
        ax17.grid('on')
        
        bplot=ax18.violinplot([paraB_kp_ns, paraB_kp_nc,paraB_kp_b,paraB_kp_m,paraB_kp_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax18.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax18.set_yticks(ticks=[-2,-1, 0], labels=['$10^{-2}$','$10^{-1}$','$10^{0}$'])
        ax18.grid('on')
        
        bplot=ax19.violinplot([btot_p_ns, btot_p_nc,btot_p_b,btot_p_m,btot_p_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax19.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax19.set_yticks(ticks=[-0.5,-0, 0.5], labels=['$10^{-0.5}$','$10^{0}$','$10^{0.5}$'])
        ax19.grid('on')
        
        
        ######################################################################################################
        
        bplot=ax20.violinplot([paraC_v_ns, paraC_v_nc,paraC_v_b,paraC_v_m,paraC_v_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax20.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax20.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        ax20.grid('on')
        
        bplot=ax21.violinplot([paraC_ks_ns, paraC_ks_nc,paraC_ks_b,paraC_ks_m,paraC_ks_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax21.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax2.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax21.set_yticks(ticks=[0,1,2], labels=['$10^{0}$','$10^{1}$','$10^{2}$'])
        ax21.grid('on')
        
        bplot=ax22.violinplot([paraC_kp_ns, paraC_kp_nc,paraC_kp_b,paraC_kp_m,paraC_kp_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax22.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax22.set_yticks(ticks=[-3,-2, -1], labels=['$10^{-3}$','$10^{-2}$','$10^{-1}$'])
        ax22.grid('on')
        
        bplot=ax23.violinplot([ctot_p_ns, ctot_p_nc,ctot_p_b,ctot_p_m,ctot_p_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax23.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax23.set_yticks(ticks=[-1.5,-1,-0.5], labels=['$10^{-1.5}$','$10^{-1}$','$10^{-0.5}$'])
        ax23.grid('on')
        
        # #[812, 649, 416, 603,682]
          
        
        py_tr_v_ns=np.array(para_t_v_ns)
        py_tr_v_nc=np.array(para_t_v_nc)
        py_tr_v_b=np.array(para_t_v_b)
        py_tr_v_m=np.array(para_t_v_m)
        py_tr_v_c=np.array(para_t_v_c)
        
        p23_v_ns=np.array(para23_v_ns)
        p23_v_nc=np.array(para23_v_nc)
        p23_v_b=np.array(para23_v_b)
        p23_v_m=np.array(para23_v_m)
        p23_v_c=np.array(para23_v_c)
        
        vrat_py_ns=p23_v_ns-py_tr_v_ns
        vrat_py_nc=p23_v_nc-py_tr_v_nc
        vrat_py_b=p23_v_b-py_tr_v_b
        vrat_py_m=p23_v_m-py_tr_v_m
        vrat_py_c=p23_v_c-py_tr_v_c
        
        py_tr_ks_ns=np.array(para_t_ks_ns)
        py_tr_ks_nc=np.array(para_t_ks_nc)
        py_tr_ks_b=np.array(para_t_ks_b)
        py_tr_ks_m=np.array(para_t_ks_m)
        py_tr_ks_c=np.array(para_t_ks_c)
        
        p23_ks_ns=np.array(para23_ks_ns)
        p23_ks_nc=np.array(para23_ks_nc)
        p23_ks_b=np.array(para23_ks_b)
        p23_ks_m=np.array(para23_ks_m)
        p23_ks_c=np.array(para23_ks_c)
        
        ksrat_py_ns=p23_ks_ns-py_tr_ks_ns
        ksrat_py_nc=p23_ks_nc-py_tr_ks_nc
        ksrat_py_b=p23_ks_b-py_tr_ks_b
        ksrat_py_m=p23_ks_m-py_tr_ks_m
        ksrat_py_c=p23_ks_c-py_tr_ks_c
        
        bplot=ax31.violinplot([vrat_py_ns,vrat_py_nc,vrat_py_b,vrat_py_m,vrat_py_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax31.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        ax31.set_yticks(ticks=[0,1,2,3], labels=['$10^{0}$','$10^{1}$','$10^{2}$','$10^{3}$'])
        ax31.grid('on')
        
        
        bplot=ax311.violinplot([ksrat_py_ns,ksrat_py_nc,ksrat_py_b,ksrat_py_m,ksrat_py_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax311.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        # ax1.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        ax311.grid('on')
        ax311.set_yticks(ticks=[-1,0,1,2], labels=['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'])

        
        p35_v_ns=np.array(para35_v_ns)
        p35_v_nc=np.array(para35_v_nc)
        p35_v_b=np.array(para35_v_b)
        p35_v_m=np.array(para35_v_m)
        p35_v_c=np.array(para35_v_c)
        
        p_35_v_ns=np.array(para_35_v_ns)
        p_35_v_nc=np.array(para_35_v_nc)
        p_35_v_b=np.array(para_35_v_b)
        p_35_v_m=np.array(para_35_v_m)
        p_35_v_c=np.array(para_35_v_c)
        
        v35rat_py_ns=p35_v_ns-p_35_v_ns
        v35rat_py_nc=p35_v_nc-p_35_v_nc
        v35rat_py_b=p35_v_b-p_35_v_b
        v35rat_py_m=p35_v_m-p_35_v_m
        v35rat_py_c=p35_v_c-p_35_v_c
        
        p35_ks_ns=np.array(para35_ks_ns)
        p35_ks_nc=np.array(para35_ks_nc)
        p35_ks_b=np.array(para35_ks_b)
        p35_ks_m=np.array(para35_ks_m)
        p35_ks_c=np.array(para35_ks_c)
        
        p_35_ks_ns=np.array(para_35_ks_ns)
        p_35_ks_nc=np.array(para_35_ks_nc)
        p_35_ks_b=np.array(para_35_ks_b)
        p_35_ks_m=np.array(para_35_ks_m)
        p_35_ks_c=np.array(para_35_ks_c)
        
        ks35rat_py_ns=p35_ks_ns-p_35_ks_ns
        ks35rat_py_nc=p35_ks_nc-p_35_ks_nc
        ks35rat_py_b=p35_ks_b-p_35_ks_b
        ks35rat_py_m=p35_ks_m-p_35_ks_m
        ks35rat_py_c=p35_ks_c-p_35_ks_c
        
        bplot=ax32.violinplot([v35rat_py_ns,v35rat_py_nc,v35rat_py_b,v35rat_py_m,v35rat_py_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax32.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        # ax1.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        ax32.set_yticks(ticks=[-1,0,1,2], labels=['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'])

        ax32.grid('on')
        
        
        bplot=ax322.violinplot([ks35rat_py_ns,ks35rat_py_nc,ks35rat_py_b,ks35rat_py_m,ks35rat_py_c], widths=(boxwidth,boxwidth,boxwidth,boxwidth,boxwidth),
                    showmeans=False,
                          showmedians=True)#,
        # Set the color of the violin patches
        for pc, color in zip(bplot['bodies'], colors):
            pc.set_facecolor(color)
            pc.set_alpha(1)
        # Set the color of the median lines
        bplot['cmedians'].set_colors('k')
        bplot['cbars'].set_colors('k')
        bplot['cmins'].set_colors('k')
        bplot['cmaxes'].set_colors('k')
        # Set the labels
        ax322.set_xticks([1, 2, 3, 4, 5], labels=['NS','NC','B','M','C'])#patch_artist=True)
        # ax1.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        # ax1.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        ax322.set_yticks(ticks=[-2,0,2,4], labels=['$10^{-2}$','$10^{0}$','$10^{2}$','$10^{4}$'])

        ax322.grid('on')
        # # allo19_desiredafteroxidase=[11,70,78,85,115,116,129,151,158,202,240,254,255,
        # #                             328,421,445,453,454,487,502,555,560,564,565,570,
        # #                             571,602,643,674,681,697,718,739,785,791,874,978,
        # #                             995]
        
        # #choose from [85, 115, 151, 202, 255, 421, 445, 453, 454, 487, 555, 
        #                 #560, 564, 570, 571, 602, 681, 697, 718, 739, 785, 874, 978, 995]
        #   #    158, 254    ,     328 502
          
        # # vgood [85, 115, 202, 328, 453, 454, 487, 555, 
        #                 #560, 564, 570, 571, 602, 681, 697, 718, 739, 785, 874, 978, 995]
        
        # # indplot=[33]
        
        # #202
                
        # # #ax2.scatter(np.abs(ff_delta),np.abs(n0delta),s=40,color='lightpink',alpha=0.5,edgecolors='k',zorder=10)#,label='')
        # # ax2.scatter(ff_delta[indexes_switch],c0delta[indexes_switch],s=40,color='lightblue',alpha=1,edgecolors='k',zorder=10, 
        # #             label='Enzyme parameter mediated swiching')
        # # ax2.scatter(ff_delta[indexes_switch_retained],c0delta[indexes_switch_retained],s=40,color='limegreen',alpha=1,edgecolors='k',zorder=10, 
        # #             label='Co-sub mediated switching')
        # # ax2.scatter(ff_delta[ind_notpicked],c0delta[ind_notpicked],s=40,color='bisque',alpha=1,edgecolors='k',zorder=0, 
        # #             label='Other')
        # # ax2.scatter(ff_delta[indexes_noswitch],c0delta[indexes_noswitch],s=40,color='lightcoral',alpha=1,edgecolors='k',zorder=20, 
        # #             label='No switching')
        
        # # # ax2.scatter(ff_delta[indplot_switchretained],n0delta[indplot_switchretained],s=40,color='red',alpha=1,edgecolors='k',zorder=20)
        # # ax2.axis([-0.05,1.05,-0.05,1.05])
        # # ax2.legend(loc=(0.15,1.05),ncols=4,fontsize=15,markerscale=3)
        # # ax2.set_xlabel('$\\Delta F_f$',fontsize=30)
        # # ax2.set_ylabel('$\\Delta N_{\\rm H}^{ m}$',fontsize=30,rotation=0)
        # # ax2.yaxis.set_label_coords(-0.09,0.46)
        # # ax3.text(10**-7,0.8,str(indplot_switchretained),fontsize=20)
        # # ax2.set_xticks(np.logspace(-9,-2,8), ['$10^{-9}$','','$10^{-7}$','',
        # #                                      '$10^{-5}$','','$10^{-3}$',''])
        # # ind_temp=indplot_switchretained
        # # args_temp=args_to_run[ind_temp]
        # # print('\n')
        # # print(ind_temp)
        # # print(args_temp[0])
        # # print(args_temp[1])
        # # print(args_temp[2])
        # # print(args_temp[3])
        # # print('\n')
        # # print(args_temp[4])
        # # print(args_temp[5])
        # # print(args_temp[6])
        # # print(args_temp[7])
        # # print('\n')
        # # print(args_temp[8])
        # # print(args_temp[9])
        # # print(args_temp[10])
        # # print(args_temp[11])
        # # print('\n')
        # # print(args_temp[12])
        # # print(args_temp[13])
        # # print('\n')
        
        # # print(args_temp[4])
        # # print(args_temp[6])
        
        
        
        # # kinplot=np.logspace(-8,-3,numkin)
        
        # # result=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
        # #                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
        # #                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
        # #                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
        # #                              args_temp[16],args_temp[17])
        
            
        # # # flux01= ode_cm.flux_enz_forward(result[:,0], result[:,1], result[:,9], result[:,10], 
        # # #                          args_temp[3][0], args_temp[3][1], args_temp[3][2])            # gap -> pep
        
        # # # flux12= ode_cm.flux_enz_forward(result[:,1], result[:,2], 1.,1., 
        # # #                          args_temp[4][0], args_temp[4][1], args_temp[4][2])           # pep -> pyr
        
        # # flux23= ode_cm.flux_enz_forward(result[:,2]**allo[0], 1.,
        # #                              args_temp[6][0],  args_temp[6][1])            # pyr -> acet
        
        # # # flux35= ode_cm.flux_enz_forward(result[:,3], result[:,4], result[:,10],result[:,9], 
        # # #                          args_temp[6][0], args_temp[6][1], args_temp[6][2])              # acet -> etn
        
        # # # flux22=ode_cm.transp(result[:,2],args_temp[9][0],args_temp[9][1])             # pyr -> pyr_mito
        
        # # flux24= ode_cm.flux_enz_forward(result[:,5],  result[:,11], 
        # #                              args_temp[8][0], args_temp[8][1])
        
        # # ax3.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
        # # ax3.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
        # # ax3.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
        # # ax3.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')
        
        # # ax3.set_xlabel('$k_{\\rm in}$',fontsize=20)
        # # # ax3.set_ylabel('flux fraction')
        # # ax3.legend(loc=(0.02,1.015),ncols=2)
        # # ax3.set_xscale('log')
        # # ax3.set_xticks(np.logspace(-8,-3,6), ['$10^{-8}$','','$10^{-6}$','',
        # #                                      '$10^{-4}$',''])#,'$10^{-3}$',''])
        
        # # # args_temp[12][3]*=0.1
        # # args_temp[13][3]*=0.1
        # # args_temp[14][3]*=0.1
        # # # args_temp[12][0]*=10.
        # # # args_temp[13][0]*=10.
        
        # # result2=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
        # #                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
        # #                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
        # #                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
        # #                              args_temp[16],args_temp[17])
        
          
        # # flux23_2= ode_cm.flux_enz_forward(result2[:,2]**allo[0], 1.,
        # #                              args_temp[6][0],  args_temp[6][1])             # pyr -> acet
        
        # # flux24_2= ode_cm.flux_enz_forward(result2[:,5],  result2[:,11], 
        # #                              args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA
        
        
        # # ax3.plot(kinplot,1-(result2[:,9]+result2[:,11])/(result2[:,9]+result2[:,10]+result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='firebrick')
        # # ax3.plot(kinplot,1-(result2[:,9])/(result2[:,9]+result2[:,10]),linestyle='--',linewidth=2,color='peru')
        # # ax3.plot(kinplot,1-(result2[:,11])/(result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
        # # ax3.plot(kinplot,1-flux24_2/(flux23_2+flux24_2),linestyle='--',linewidth=2,color='royalblue')
        
        
        # boxwidth=0.8
        
        # bplot=ax5.violinplot([paraC_v , paraC_v_n2,paraC_v_n],widths=(boxwidth,boxwidth,boxwidth),
        #             showmeans=False,
        #                   showmedians=True)#,
        #             # patch_artist=True)
        
        # colors=['limegreen','lightblue','lightcoral']
        # # Set the color of the violin patches
        # for pc, color in zip(bplot['bodies'], colors):
        #     pc.set_facecolor(color)
        #     pc.set_alpha(1)
        
        # # Set the color of the median lines
        # bplot['cmedians'].set_colors('k')
        # bplot['cbars'].set_colors('k')
        # bplot['cmins'].set_colors('k')
        # bplot['cmaxes'].set_colors('k')
        
        # # Set the labels
        # ax5.set_xticks([1, 2, 3], labels=['C', 'E', 'N'])#patch_artist=True)
        
        # # colors=['limegreen','lightblue','lightcoral']
        # # for patch, color in zip(bplot['boxes'], colors):
        # #     patch.set_facecolor(color)
        
        # # for median in bplot['medians']:
        # #     median.set_color('black')
        
        # # ax5.set_title('Total cosubstrate (Mito)', fontsize=20)
        # ax5.set_title('$V_{\\rm max}^{ N^m}$', fontsize=15,y=1.01)
        
        # # ax5.set_xticklabels(labels=['C','E','N'], 
        # #                        rotation=0)
        # ax5.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        # ax5.grid('on')
        
        
        
        # bplot=ax6.violinplot([paraC_ks, paraC_ks_n2,paraC_ks_n], widths=(boxwidth,boxwidth,boxwidth),
        #             showmeans=False,
        #                   showmedians=True)#,
        #             # patch_artist=True)
        
        # colors=['limegreen','lightblue','lightcoral']
        # # Set the color of the violin patches
        # for pc, color in zip(bplot['bodies'], colors):
        #     pc.set_facecolor(color)
        #     pc.set_alpha(1)
        
        # # Set the color of the median lines
        # bplot['cmedians'].set_colors('k')
        # bplot['cbars'].set_colors('k')
        # bplot['cmins'].set_colors('k')
        # bplot['cmaxes'].set_colors('k')
        
        # # Set the labels
        # ax6.set_xticks([1, 2, 3], labels=['C', 'E', 'N'])#patch_artist=True)
        
        # # colors=['limegreen','lightblue','lightcoral']
        # # for patch, color in zip(bplot['boxes'], colors):
        # #     patch.set_facecolor(color)
        # # for median in bplot['medians']:
        # #     median.set_color('black')                 
        # # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        # ax6.set_title('$K_{\\rm S}^{N_{+}^m}$', fontsize=15,y=1.01)
        
        # # ax6.set_xticklabels(labels=['C','E','N'], 
        # #                        rotation=0)
        # ax6.set_yticks(ticks=[-1,0,1,2], labels=['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'])
        # # ax6.set_yticks(ticks=np.linspace(-4,0))
        
        # ax6.grid('on')
        
        
        
        # bplot=ax7.violinplot([paraC_kp, paraC_kp_n2,paraC_kp_n],widths=(boxwidth,boxwidth,boxwidth),
        #             showmeans=False,
        #                   showmedians=True)#,
        #             # patch_artist=True)
        
        # colors=['limegreen','lightblue','lightcoral']
        # # Set the color of the violin patches
        # for pc, color in zip(bplot['bodies'], colors):
        #     pc.set_facecolor(color)
        #     pc.set_alpha(1)
        
        # # Set the color of the median lines
        # bplot['cmedians'].set_colors('k')
        # bplot['cbars'].set_colors('k')
        # bplot['cmins'].set_colors('k')
        # bplot['cmaxes'].set_colors('k')
        
        # # Set the labels
        # ax7.set_xticks([1, 2, 3], labels=['C', 'E', 'N'])#patch_artist=True)
        
        # # colors=['limegreen','lightblue','lightcoral']
        # # for patch, color in zip(bplot['boxes'], colors):
        # #     patch.set_facecolor(color)
        # # for median in bplot['medians']:
        # #     median.set_color('black')                 
        # # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        # ax7.set_title('$K_{\\rm P}^{N_{\\rm H}^m}$', fontsize=15,y=1.01)
        
        # # ax7.set_xticklabels(labels=['C','E','N'], 
        # #                        rotation=0)
        # ax7.set_yticks(ticks=[-3,-2,-1,0], labels=['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'])
        
        # ax7.grid('on')
        
        # # # get the data for original params
        # # ind_temp=indplot_noncosub
        # # args_temp=args_to_run[ind_temp]
        # # print('\n')
        # # print(ind_temp)
        # # print(args_temp[0])
        # # print(args_temp[1])
        # # print(args_temp[2])
        # # print(args_temp[3])
        # # print('\n')
        # # print(args_temp[4])
        # # print(args_temp[5])
        # # print(args_temp[6])
        # # print(args_temp[7])
        # # print('\n')
        # # print(args_temp[8])
        # # print(args_temp[9])
        # # print(args_temp[10])
        # # print(args_temp[11])
        # # print('\n')
        # # print(args_temp[12])
        # # print(args_temp[13])
        # # print('\n')
        
        # # print(args_temp[4])
        # # print(args_temp[6])
        
        
        
        # # kinplot=np.logspace(-8,-3,numkin)
        
        # # result=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
        # #                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
        # #                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
        # #                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
        # #                              args_temp[16],args_temp[17])
        
            
        # # # flux01= ode_cm.flux_enz_forward(result[:,0], result[:,1], result[:,9], result[:,10], 
        # # #                          args_temp[3][0], args_temp[3][1], args_temp[3][2])            # gap -> pep
        
        # # # flux12= ode_cm.flux_enz_forward(result[:,1], result[:,2], 1.,1., 
        # # #                          args_temp[4][0], args_temp[4][1], args_temp[4][2])           # pep -> pyr
        
        # # flux23= ode_cm.flux_enz_forward(result[:,2]**allo[0], 1.,
        # #                              args_temp[6][0],  args_temp[6][1])            # pyr -> acet
        
        # # # flux35= ode_cm.flux_enz_forward(result[:,3], result[:,4], result[:,10],result[:,9], 
        # # #                          args_temp[6][0], args_temp[6][1], args_temp[6][2])              # acet -> etn
        
        # # # flux22=ode_cm.transp(result[:,2],args_temp[9][0],args_temp[9][1])             # pyr -> pyr_mito
        
        # # flux24= ode_cm.flux_enz_forward(result[:,5],  result[:,11], 
        # #                              args_temp[8][0], args_temp[8][1])
        
        # # ax4.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
        # # ax4.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
        # # ax4.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
        # # ax4.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')
        
        # # ax4.set_xlabel('$k_{\\rm in}$',fontsize=20)
        # # # ax3.set_ylabel('flux fraction')
        # # ax4.set_xscale('log')
        # # ax4.set_xticks(np.logspace(-8,-3,6), ['$10^{-8}$','','$10^{-6}$','',
        # #                                      '$10^{-4}$',''])#,'$10^{-3}$',''])
        
        # # # args_temp[12][3]*=0.1
        # # args_temp[13][3]*=0.1
        # # args_temp[14][3]*=0.1
        # # # args_temp[12][0]*=10.
        # # # args_temp[13][0]*=10.
        
        # # result2=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
        # #                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
        # #                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
        # #                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
        # #                              args_temp[16],args_temp[17])
        
          
        # # flux23_2= ode_cm.flux_enz_forward(result2[:,2]**allo[0], 1.,
        # #                              args_temp[6][0],  args_temp[6][1])             # pyr -> acet
        
        # # flux24_2= ode_cm.flux_enz_forward(result2[:,5],  result2[:,11], 
        # #                              args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA
        
        
        # # ax4.plot(kinplot,1-(result2[:,9]+result2[:,11])/(result2[:,9]+result2[:,10]+result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='firebrick')
        # # ax4.plot(kinplot,1-(result2[:,9])/(result2[:,9]+result2[:,10]),linestyle='--',linewidth=2,color='peru')
        # # ax4.plot(kinplot,1-(result2[:,11])/(result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
        # # ax4.plot(kinplot,1-flux24_2/(flux23_2+flux24_2),linestyle='--',linewidth=2,color='royalblue')
        
        
        
        # print('param diffs here')
        # print(10**np.array(np.mean(np.log10(para01_n2),0)-np.mean(np.log10(para01_n),0)))
        # # print(10**np.array(np.var(np.log10(para01),0)-np.log10(params01var)))
        # # print(np.mean(para01,0)/params01avg)
        # print(10**np.array(np.mean(np.log10(para1E_n2),0)-np.mean(np.log10(para1E_n),0)))
        # print(10**np.array(np.mean(np.log10(paraE2_n2),0)-np.mean(np.log10(paraE2_n),0)))
        # # print(10**np.array(np.var(np.log10(para12),0)-np.log10(params12var)))
        # print(10**np.array(np.mean(np.log10(para23_n2),0)-np.mean(np.log10(para23_n),0)))
        # # print(10**np.array(np.var(np.log10(para23),0)-np.log10(params23var)))
        # print(10**np.array(np.mean(np.log10(para35_n2),0)-np.mean(np.log10(para35_n),0)))
        # # print(10**np.array(np.var(np.log10(para35),0)-np.log10(params35var)))
        # print('\n')
        # print(10**np.array(np.mean(np.log10(para_24_n2),0)-np.mean(np.log10(para_24_n),0)))
        # # print(10**np.array(np.var(np.log10(para_24),0)-np.log10(params_24var)))
        # print(10**np.array(np.mean(np.log10(para_35_n2),0)-np.mean(np.log10(para_35_n),0)))
        # # print(10**np.array(np.var(np.log10(para_35),0)-np.log10(params_35var)))
        # print(10**np.array(np.mean(np.log10(para_t_n2),0)-np.mean(np.log10(para_t_n),0)))
        # # print(10**np.array(np.var(np.log10(para_t),0)-np.log10(params_trans_var)))
        # print('\n')
        # print(10**np.array(np.mean(np.log10(paraB_n2),0)-np.mean(np.log10(paraB_n),0)))
        # # print(10**np.array(np.var(np.log10(paraB),0)-np.log10(paramsBvar)))
        # print(10**np.array(np.mean(np.log10(paraC_n2),0)-np.mean(np.log10(paraC_n),0)))
        # # print(10**np.array(np.var(np.log10(paraC),0)-np.log10(paramsCvar)))
        # print(10**np.array(np.mean(np.log10(ctot_n2),0)-np.mean(np.log10(ctot_n),0)))
        # # print(10**np.array(np.var(np.log10(ctot),0)-np.log10(ctotvar)))
        # print(10**np.array(np.mean(np.log10(btot_n2),0)-np.mean(np.log10(btot_n),0)))
        # # print(10**np.array(np.var(np.log10(btot),0)-np.log10(btotvar)))
        # # np.savez(savename, flux_cyto=ftlist, flux_mito=fblist, flux_diffn=fdlist,
        # #          fluxfrac=fflist, fluxfrac_grad=ffgradlist, fluxfrac_sense=ffsenselist,
        # #          fluxfrac_delta=ffdeltalist, gradmax=gradmaxlist, gradmaxind=gradmaxindlist,
        # #          sensemax=sensemaxlist, sensemaxind=sensemaxindlist)
        # print('param diffs end')
        
        # print()
        # boxwidth=0.8
        
        # # ax8.boxplot([para23_v_n2 , para23_v_n],widths=(boxwidth,boxwidth))
                          
        # # # ax5.set_title('Total cosubstrate (Mito)', fontsize=20)
        # # ax8.set_title('$V_{\\rm M}~ {\\rm Pyr. \\rightarrow Acet.}$', fontsize=15,y=1.01)
        
        # # ax8.set_xticklabels(labels=['S','N'], 
        # #                        rotation=0)
        # # ax8.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        # # ax8.grid('on')
        
        # # p23vn2=np.array(para23_v_n2)
        # # ptvn2=np.array(para_t_v_n2)
        # # p23vn=np.array(para23_v_n)
        # # ptvn=np.array(para_t_v_n)
        
        # # vratios_n2 = p23vn2-ptvn2
        # # vratios_n= p23vn-ptvn
        
        # # ax8.boxplot([vratios_n2 , vratios_n],widths=(boxwidth,boxwidth))
                          
        # # # ax5.set_title('Total cosubstrate (Mito)', fontsize=20)
        # # ax8.set_title('$V_{\\rm M,23} / V_{\\rm M,T}  $', fontsize=15,y=1.01)
        
        # # ax8.set_xticklabels(labels=['S','N'], 
        # #                        rotation=0)
        # # # ax8.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
        # # ax8.grid('on')
        
        # # ax9.boxplot([para_t_v_n2, para_t_v_n], widths=(boxwidth,boxwidth))
                          
        # # # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        # # ax9.set_title('$V_{\\rm M}~ {\\rm Pyr. ~Trans.}$', fontsize=15,y=1.01)
        
        # # ax9.set_xticklabels(labels=['S','N'], 
        # #                        rotation=0)
        # # ax9.set_yticks(ticks=[-2,-1], labels=['$10^{-2}$','$10^{-1}$'])
        # # # ax6.set_yticks(ticks=np.linspace(-4,0))
        
        # # ax9.grid('on')
        
        # # p23ksn2=np.array(para23_ks_n2)
        # # ptksn2=np.array(para_t_ks_n2)
        # # p23ksn=np.array(para23_ks_n)
        # # ptksn=np.array(para_t_ks_n)
        
        # # ksratios_n2 = p23ksn2-ptksn2
        # # ksratios_n= p23ksn-ptksn
        
        # # ax9.boxplot([ksratios_n2, ksratios_n], widths=(boxwidth,boxwidth))
                          
        # # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        # # ax9.set_title('$K_{\\rm S,T} / K_{\\rm S,T}  $', fontsize=15,y=1.01)
        
        # # ax9.set_xticklabels(labels=['S','N'], 
        # #                        rotation=0)
        # # # ax9.set_yticks(ticks=[-2,-1], labels=['$10^{-2}$','$10^{-1}$'])
        # # # ax6.set_yticks(ticks=np.linspace(-4,0))
        
        # # ax9.grid('on')
        
        # # p23vn2=np.array(para23_v_n2)
        # # ptvn2=np.array(para_t_v_n2)
        # # p23vn=np.array(para23_v_n)
        # # ptvn=np.array(para_t_v_n)
        
        # # vratios_n2 = p23vn2-ptvn2
        # # vratios_n= p23vn-ptvn
        
        # ptrans_v=np.array(para_t_v)
        # ptrans_v_n2=np.array(para_t_v_n2)
        # ptrans_v_n=np.array(para_t_v_n)
        
        # p23_v=np.array(para23_v)
        # p23_v_n2=np.array(para23_v_n2)
        # p23_v_n=np.array(para23_v_n)
        
        # vrat=p23_v-ptrans_v
        # vrat_n2=p23_v_n2-ptrans_v_n2
        # vrat_n=p23_v_n-ptrans_v_n
        
        # ptrans_ks=np.array(para_t_ks)
        # ptrans_ks_n2=np.array(para_t_ks_n2)
        # ptrans_ks_n=np.array(para_t_ks_n)
        
        # p23_ks=np.array(para23_ks)
        # p23_ks_n2=np.array(para23_ks_n2)
        # p23_ks_n=np.array(para23_ks_n)
        
        # ksrat=p23_ks-ptrans_ks
        # ksrat_n2=p23_ks_n2-ptrans_ks_n2
        # ksrat_n=p23_ks_n-ptrans_ks_n
        
        
        # bplot=ax10.violinplot([vrat,vrat_n2, vrat_n],widths=(boxwidth,boxwidth,boxwidth),
        #                       showmeans=False,
        #                   showmedians=True)#,
        #             # patch_artist=True)
        
        # colors=['limegreen','lightblue','lightcoral']
        # # Set the color of the violin patches
        # for pc, color in zip(bplot['bodies'], colors):
        #     pc.set_facecolor(color)
        #     pc.set_alpha(1)
        
        # # Set the color of the median lines
        # bplot['cmedians'].set_colors('k')
        # bplot['cbars'].set_colors('k')
        # bplot['cmins'].set_colors('k')
        # bplot['cmaxes'].set_colors('k')
        
        # # Set the labels
        # ax10.set_xticks([1, 2, 3], labels=['C', 'E', 'N'])
        # # for patch, color in zip(bplot['boxes'], colors):
        # #     patch.set_facecolor(color)
        # # for median in bplot['medians']:
        # #     median.set_color('black')
                       
        # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        # ax10.set_title('$V_{\\rm max}^{{\\rm Pyr}~ \\rightarrow~ {\\rm AcCHO}} / V_{\\rm max}^{\\rm Pyr. ~Trans.}$', fontsize=15,y=1.01,x=0.6)
        
        # # ax10.set_xticklabels(labels=['C','E','N'], 
        #                        # rotation=0)
        # ax10.set_yticks(ticks=[-0,1,2,3], labels=['$10^{0}$','$10^{1}$',\
        #                                           '$10^{2}$','$10^{3}$'])
        
        # ax10.grid('on')
        
        
        # bplot=ax11.violinplot([ksrat,ksrat_n2, ksrat_n],widths=(boxwidth,boxwidth,boxwidth),
        #                   showmeans=False,
        #               showmedians=True)#,
        #             # patch_artist=True)
        
        
        # colors=['limegreen','lightblue','lightcoral']
        # # Set the color of the violin patches
        # for pc, color in zip(bplot['bodies'], colors):
        #     pc.set_facecolor(color)
        #     pc.set_alpha(1)
        
        # # Set the color of the median lines
        # bplot['cmedians'].set_colors('k')
        # bplot['cbars'].set_colors('k')
        # bplot['cmins'].set_colors('k')
        # bplot['cmaxes'].set_colors('k')
        
        # # Set the labels
        # ax11.set_xticks([1, 2, 3], labels=['C', 'E', 'N'])
        # # for patch, color in zip(bplot['boxes'], colors):
        # #     patch.set_facecolor(color)
        # # for median in bplot['medians']:
        # #     median.set_color('black')
                       
        # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        
        
        # # colors=['limegreen','lightblue','lightcoral']
        # # for patch, color in zip(bplot['boxes'], colors):
        # #     patch.set_facecolor(color)
        # # for median in bplot['medians']:
        # #     median.set_color('black')               
        # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        # ax11.set_title('$K_{\\rm M}^{{\\rm Pyr}~ \\rightarrow ~{\\rm AcCHO}} / K_{\\rm M}^{\\rm Pyr. ~Trans.}$', fontsize=15,y=1.01,x=0.6)
        
        # # ax11.set_xticklabels(labels=['C','E','N'], 
        # #                        rotation=0)
        # ax11.set_yticks(ticks=[0,1,2], labels=['$10^{0}$',\
        #                                            '$10^{1}$','$10^{2}$'])
        
        # ax11.grid('on')
        
        
        
        
        # # bplot=ax12.boxplot([paraB_ks_n2, paraB_ks_n],widths=(boxwidth,boxwidth),
        # #             patch_artist=True)
        
        # # colors=['lightblue','lightcoral']
        # # for patch, color in zip(bplot['boxes'], colors):
        # #     patch.set_facecolor(color)
        # # for median in bplot['medians']:
        # #     median.set_color('black')            
        # # # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
        # # ax12.set_title('$K_{\\rm S}^{\\rm N}~ {\\rm (Cyto)}$', fontsize=15,y=1.01)
        
        # # ax12.set_xticklabels(labels=['S','N'], 
        # #                        rotation=0)
        # # ax12.set_yticks(ticks=[-2,-1,-0,1], labels=['$10^{-2}$','$10^{-1}$','$10^{0}$','$10^{1}$'])
        
        # # ax12.grid('on')
        
        
        # # ax14.hist(ff_delta,bins=20)
        # # ax14.set_xlabel('flux frac delta')
        
        
        # # allo=[1.,1.,1.,1.]
        # # savepath='../data/allo_'+str(allo)
        
        # # filename='fixed_new_3_seed'+seedstr+'_n_'+str(n)
        # # savename=savepath+filename
        # # print(savepath)
        # # print(filename)
        # # p=Path(savepath)
        
        # # filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
        # # print(filewant)
        
        # # kinvals=args_to_run[0][1]
        
        # # dat=np.load(filewant)
        # # ff_delta=dat['fluxfrac_delta']
        # # ff_delta=np.abs(ff_delta)
        # # gradmaxind=dat['gradmaxind']
        # # sensemaxind=dat['sensemaxind']
        # # fluxfrac=dat['fluxfrac']
        # # ff_grad=dat['fluxfrac_grad']
        # # sensemax=dat['sensemax']
        # # gradmax=dat['gradmax']
        # # ffsense_vec=dat['fluxfrac_sense']
        # # b0delta=dat['b0delta']
        # # c0delta=dat['c0delta']
        # # n0delta=dat['n0delta']
        
        # # ax13.scatter(ff_delta,n0delta,s=40,color='lightpink',alpha=0.5,edgecolors='k',zorder=20)#,label='')
        
        
        # # ax15.hist(ff_delta,bins=20)
        # # ax15.set_xlabel('flux frac delta')
        
        # # print(indexes_switch_retained)
        
        # # axs2[0,0].set_xticklabels(labels=['co','No Switch'], 
        # #                         rotation=90)
        
        # # ax2.text(0.4,1.3,'log10 diffn rate ='+str(np.round(np.log10(diffval)))+', diffn scale ='+str(np.round(diffscale)), fontsize=30)



        # # print(indexes_switch_retained)
        if savefig==1:
            plt.savefig('Figure_S6.png',bbox_inches='tight',format='png')
            plt.savefig('Figure_S6.eps',bbox_inches='tight',format='eps')
            plt.savefig('Figure_S6.pdf',bbox_inches='tight',format='pdf')
            plt.savefig('Figure_S6.svg',bbox_inches='tight',format='svg')
        

    
plt.show()