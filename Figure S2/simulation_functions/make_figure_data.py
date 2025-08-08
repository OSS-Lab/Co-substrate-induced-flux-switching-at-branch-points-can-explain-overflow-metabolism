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
            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,10.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.1,1.,1.],
            bsynthdeg=[10**-2,10**-2,10**-2],
            reversible=False,Tmax=1000,tpts=100):
    
    result,time_points,odefunc,params = mp.make_time_series(lengthparams,B_position_params,
            params_global,params_B,
            params_branch,
            params_inout,
            bsynthdeg,
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



