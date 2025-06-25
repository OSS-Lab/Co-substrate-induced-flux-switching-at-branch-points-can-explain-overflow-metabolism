#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 12:49:32 2024

@author: robert
"""

import numpy as np

import matplotlib.pyplot as plt

import make_odes as mp
import time_series_plots as ts
import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 9),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)


def get_fluxfractions(lengthparams=[3,2,2],B_position_params=([0,0,0],[0,0,0],[1,0],[1,0]),
            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,10.],allo=[1.,1.,1.,1.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.1,1.,1.],
            reversible=False,Tmax=1000,tpts=100):
    
    result,time_points,odefunc,params = mp.make_time_series(lengthparams,B_position_params,
            params_global,params_B,allo,
            params_branch,
            params_inout,
            reversible,Tmax,tpts)
    
    isbuildup=mp.check_buildup(result,odefunc,params,reversible)
    
    flux_top,flux_bottom,bfrac,brat=mp.get_branch_fluxes(result,odefunc,params,lengthparams,B_position_params,reversible)
    fluxfrac=flux_bottom/(flux_top+flux_bottom)
    flux_rat=flux_bottom/flux_top
    
    return isbuildup, fluxfrac, bfrac, flux_rat,brat

def get_fluxfractions_pinit(lengthparams=[3,2,2],B_position_params=([0,0,0],[0,0,0],[1,0],[1,0]),
            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,10.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.1,1.,1.],
            reversible=False,Tmax=1000,tpts=100,pinit=[],usepinit=False):
    
    if usepinit==False:
        result,time_points,odefunc,params = mp.make_time_series(lengthparams,B_position_params,
                params_global,params_B,
                params_branch,
                params_inout,
                reversible,Tmax,tpts)
    else:
        result,time_points,odefunc,params = mp.make_time_series_pinit(lengthparams,B_position_params,
                params_global,params_B,
                params_branch,
                params_inout,
                reversible,Tmax,tpts,pinit,True)
    
    isbuildup=mp.check_buildup(result,odefunc,params,reversible)
    
    flux_top,flux_bottom,bfrac,brat=mp.get_branch_fluxes(result,odefunc,params,lengthparams,B_position_params,reversible)
    fluxfrac=flux_bottom/(flux_top+flux_bottom)
    flux_rat=flux_bottom/flux_top
    
    return isbuildup, fluxfrac, bfrac, flux_rat,brat


def FIGURE1_ff_vs_kin_vary24():
    
    fig, axs = plt.subplots(2, 2)
    V24_vals=np.array([0.1,0.5,1,2,10])
    
    cols=['blue','lightblue','black','pink','red']
    
    kinvals=np.logspace(-4,1,51)
    # no cosubs
    for j in range(len(V24_vals)):
        V24=V24_vals[j]
        print(V24)
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac=get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[0],[0]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[0.1,1.,1.,1.,0.01],
                              params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=False,Tmax=1000000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=fluxfrac
            
            
        axs[1,0].plot(kinvals,ff_vec,color=cols[j])
        axs[1,0].set_xscale('log')
        
    # cosubs on branches
    for V24 in V24_vals:
        print(V24)
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac=get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[1,1.,1.,1.,0.01],
                              params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=False,Tmax=1000000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=fluxfrac
            
            
        axs[1,1].plot(kinvals,ff_vec)
        axs[1,1].set_xscale('log')

#FIGURE1_ff_vs_kin_vary24()

def FIGURE2_ff_vs_kin_vary24():
    
    fig, axs = plt.subplots(2, 2)
    Vmaxb_vals=np.logspace(-3,1,5)
    Vmaxb_vals=np.logspace(-7,1,5)

  #  Vmaxb_vals=[0]
    cols=['blue','lightblue','black','pink','red']
    
    kinvals=np.logspace(-8,1,21)
    # no cosubs
   
        
    # cosubs on branches
    for j in range(len(Vmaxb_vals)):
        Vmaxb=Vmaxb_vals[j]
        print(Vmaxb)
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac=get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[Vmaxb,1.,1.,1,1],
                              params_branch=[1.,1.,1.,1.,10.,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=False,Tmax=1000000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=fluxfrac
            
            
        axs[1,1].plot(kinvals,ff_vec,color=cols[j])
        axs[1,1].set_xscale('log')
        
    for j in range(len(Vmaxb_vals)):
        Vmaxb=Vmaxb_vals[j]
        print(Vmaxb)
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac=get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[Vmaxb,1.,1.,1,1],
                              params_branch=[1.,1.,1.,1.,0.01,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=False,Tmax=1000000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=fluxfrac
            
            
        axs[1,0].plot(kinvals,ff_vec,color=cols[j])
        axs[1,0].set_xscale('log')

# FIGURE2_ff_vs_kin_vary24()

def FIGURE2_2_ff_vs_kin_vary24():
    
    fig, axs = plt.subplots(2, 2)
    # Btotvals=np.logspace(-3,1,5)
    Btotvals=np.logspace(-2,2,5)

    cols=['blue','lightblue','black','pink','red']
    
    kinvals=np.logspace(-8,1,21)
    # no cosubs
   
        
    # cosubs on branches
    for j in range(len(Btotvals)):
        Btot=Btotvals[j]
        print(Btot)
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        
        b_vec=np.zeros(np.size(kinvals))
        b_vec[:]=np.nan
        
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac=get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[0.1,1.,1.,1,Btot],
                              params_branch=[1.,1.,1.,1.,10.,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=False,Tmax=1000000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=fluxfrac
                b_vec[i]=bfrac
            
            
        axs[1,1].plot(kinvals,ff_vec,color=cols[j])
        axs[1,1].plot(kinvals,b_vec,color=cols[j],linestyle=':')
        axs[1,1].set_xscale('log')
        
    for j in range(len(Btotvals)):
        Btot=Btotvals[j]
        print(Btot)
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        b_vec=np.zeros(np.size(kinvals))
        b_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac=get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[0.1,1.,1.,1,Btot],
                              params_branch=[1.,1.,1.,1.,10.,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=False,Tmax=1000000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=fluxfrac
                b_vec[i]=bfrac

            
        axs[1,0].plot(kinvals,ff_vec,color=cols[j])
        axs[1,0].plot(kinvals,b_vec,color=cols[j],linestyle=':')
        axs[1,0].set_xscale('log')





# keqbvals=np.logspace(-3,3,7)
# for keqb in keqbvals:
#     print(keqb)
#     isbuildup, flux_frac=get_fluxfractions(lengthparams=[3,2,2],B_position_params=([1,0,0],[0,0,0],[1,0],[1,0]),
#                 params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,keqb,10.],
#                 params_branch=[1.,1.,1.,1.,
#                                1.,1.,1.,1.],
#                 params_inout=[0.1,1.,1.],
#                 reversible=False,Tmax=1000000,tpts=1000)
    
#     print(isbuildup)
#     print(flux_frac)
