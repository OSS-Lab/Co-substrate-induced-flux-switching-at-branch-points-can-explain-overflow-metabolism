#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 16:14:47 2024

@author: robert
"""
import numpy as np

import matplotlib.pyplot as plt

import make_odes as mp

import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 9),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)

def plot_timeseries(lengthparams=[3,2,2],B_position_params=([0,0,0],[0,0,0],[1,0],[1,0]),
            params_global=[1.,1.,1.,1.],params_B=[1.,1.,1.,1.,10.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.1,1.,1.],
            reversible=False,Tmax=100,tpts=100,linest='-'):
    
    result,time_points,odefunc,params = mp.make_time_series(lengthparams,B_position_params,
            params_global,params_B,
            params_branch,
            params_inout,
            reversible,Tmax,tpts)

    # Plot the concentrations of all metabolites as a time series
    fig, axs = plt.subplots(2, 2)
   # print(result)
    # Plot upstream metabs
  #  print(result[-1,:])
  #  print(result[-1,:-2])
    for i in range(lengthparams[0]+1):
        axs[0,0].plot(time_points, result[:, i], label=f'M{i}',linestyle=linest)
        axs[0,0].legend()
        axs[0,0].set_title('Upstream metabolites')
   
    # Plot topbranch
    for i in range(lengthparams[0]+1, lengthparams[0]+lengthparams[1]+1):
   #     print(result[:,i])
        axs[0,1].plot(time_points, result[:, i], label=f'M{i}',linestyle=linest)
        axs[0,1].legend()
        axs[0,1].set_title('Top branch metabolites')
    
    # Plot bottombranch
    for i in range(lengthparams[0]+lengthparams[1]+1, lengthparams[0]+lengthparams[1]+lengthparams[2]+1):
     #   print(result[:,i])
        axs[1,0].plot(time_points, result[:, i], label=f'M{i}',linestyle=linest)
        axs[1,0].legend()
        axs[1,0].set_title('Bottom branch metabolites')
    
    # Plot B0 and B1
   # print(result[:,-2])
    axs[1,1].plot(time_points,result[:,-2],label='B0',linestyle=linest)
    axs[1,1].plot(time_points,result[:,-1],label='B1',linestyle=linest)
    axs[1,1].legend()
    axs[1,1].set_title('Cosubs')
    
    for ax in axs.flat:
        ax.set(xlabel='Time', ylabel='Concentration',yscale='log')
        
        
    plt.subplots_adjust(wspace=0.2, hspace=0.4)
    
    print(result[-1,-2]+result[-1,-1])
    
def plot_timeseries_nov(lengthparams=[5,2,2],B_position_params=([1,0,0,0,0],[0,0,0,0,0],[0,1],[1,0]),
            params_global=[1.,1.,1.,1.],params_B=[1.,1.,1.,1.,10.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.1,1.,1.],
            reversible=False,Tmax=100,tpts=100,linest='-'):
    
    result,time_points,odefunc,params = mp.make_time_series(lengthparams,B_position_params,
            params_global,params_B,
            params_branch,
            params_inout,
            reversible,Tmax,tpts)

    # Plot the concentrations of all metabolites as a time series
    fig, axs = plt.subplots(2, 2)
   # print(result)
    # Plot upstream metabs
  #  print(result[-1,:])
  #  print(result[-1,:-2])
    for i in range(lengthparams[0]+1):
        axs[0,0].plot(time_points, result[:, i], label=f'M{i}',linestyle=linest)
        axs[0,0].legend()
        axs[0,0].set_title('Upstream metabolites')
   
    # Plot topbranch
    for i in range(lengthparams[0]+1, lengthparams[0]+lengthparams[1]+1):
   #     print(result[:,i])
        axs[0,1].plot(time_points, result[:, i], label=f'M{i}',linestyle=linest)
        axs[0,1].legend()
        axs[0,1].set_title('Top branch metabolites')
    
    # Plot bottombranch
    for i in range(lengthparams[0]+lengthparams[1]+1, lengthparams[0]+lengthparams[1]+lengthparams[2]+1):
     #   print(result[:,i])
        axs[1,0].plot(time_points, result[:, i], label=f'M{i}',linestyle=linest)
        axs[1,0].legend()
        axs[1,0].set_title('Bottom branch metabolites')
    
    # Plot B0 and B1
   # print(result[:,-2])
    axs[1,1].plot(time_points,result[:,-2],label='B0',linestyle=linest)
    axs[1,1].plot(time_points,result[:,-1],label='B1',linestyle=linest)
    axs[1,1].legend()
    axs[1,1].set_title('Cosubs')
    
    for ax in axs.flat:
        ax.set(xlabel='Time', ylabel='Concentration',yscale='log')
        
        
    plt.subplots_adjust(wspace=0.2, hspace=0.4)
    
    print(result[-1,-2]+result[-1,-1])
    plt.show()
    
#plot_timeseries_nov()
# vmbvals=np.logspace(-3,0,7)
# for vmb in vmbvals:				
#     print(vmb)
#     plot_timeseries(lengthparams=[3,1,1],B_position_params=([1,1,0],[0,0,0],[1],[1]),
#                   params_global=[1.,1.,1.,1.],params_B=[vmb,1.,1.,1.,10.],
#                   params_branch=[1.,1.,1.,1.,
#                                   1.,1.,1.,1.],
#                   params_inout=[0.0001,1.,1.],
#                   reversible=False,Tmax=100000,tpts=10000,linest='-')
    
#     plot_timeseries(lengthparams=[3,1,1],B_position_params=([1,1,0],[0,0,0],[1],[1]),
#                   params_global=[1.,1.,1.,1.],params_B=[vmb,1.,1.,1.,10.],
#                   params_branch=[1.,1.,1.,1.,
#                                   1.,1.,1.,1.],
#                   params_inout=[0.001,1.,1.],
#                   reversible=False,Tmax=100000,tpts=10000,linest='--')

# plot_timeseries(lengthparams=[3,2,2],B_position_params=([1,0,0],[0,0,0],[0,1],[0,1]),
#               params_global=[1.,1.,1.,1.],params_B=[1.,1.,1.,1.,0.01],
#               params_branch=[1.,1.,1.,1.,
#                               1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#               reversible=False,Tmax=1000)


# plot_timeseries(lengthparams=[3,2,2],B_position_params=([1,0,0],[0,0,1],[0,1],[0,1]),
#               params_global=[1.,1.,1.,1.],params_B=[1.,1.,1.,1.,0.01],
#               params_branch=[1.,1.,1.,1.,
#                               1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#               reversible=False,Tmax=1000)

# keqbvals=np.logspace(-3,3,7)
# for keqb in keqbvals:
#     plot_timeseries(lengthparams=[3,2,2],B_position_params=([1,0,0],[0,0,0],[1,0],[1,0]),
#                 params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,keqb,10.],
#                 params_branch=[1.,1.,1.,1.,
#                                 1.,1.,1.,1.],
#                 params_inout=[0.1,1.,1.],
#                 reversible=False,Tmax=1000000,tpts=1000)

# keqbvals=np.logspace(-2,2,5)
# for keqb in keqbvals:
#     plot_timeseries(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
#                   params_global=[1.,1.,1.,1.],params_B=[1.,1.,1.,keqb,100],
#                   params_branch=[1.,1.,1.,1.,
#                                   1.,1.,1.,1.],
#                   params_inout=[0.01,1.,1.],
#                   reversible=False,Tmax=1000)


# plot_timeseries(lengthparams=[3,2,2],B_position_params=([1,0,0],[0,0,0],[1,0],[0,1]),
#               params_global=[1.,1.,1.,1.],params_B=[1.,1.,1.,10.,0.1],
#               params_branch=[1.,1.,1.,1.,
#                               1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#               reversible=False,Tmax=1000)

# plot_timeseries(lengthparams=[3,2,2],B_position_params=([1,0,0],[0,0,1],[0,1],[0,1]),
#               params_global=[1.,1.,1.,1.],params_B=[1.,1.,1.,1.,100.],
#               params_branch=[100.,1.,1.,1.,
#                               100.,1.,1.,1.],
#               params_inout=[0.01,0.1,0.1],
#               reversible=False,Tmax=1000)
