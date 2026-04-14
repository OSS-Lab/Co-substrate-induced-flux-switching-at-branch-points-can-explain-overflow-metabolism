#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  1 11:28:40 2024

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

import itertools

# sys.path.append('../data/megadata')
# sys.path.append('../data/heatmapdata')
sys.path.append('../simulation_functions')


import make_heatmaps as mh





####### these parameters pretty much never change #######
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

kin_min=-8
kin_max=0
num_kin=41

Tmax=100000
tpts=10000


reversible=False


vals_b_in=np.logspace(-5,2,29)#np.logspace(-5,2,8)
vals_b_out=np.logspace(-5,2,29)#np.logspace(-5,2,8)


vals_CKcB=[0.01]#np.logspace(-1,1,3)
vals_Btot=[10.]#np.logspace(-2,0,3)
vals_CKct=[1.]
vals_CKcb=[1.]

### ONLY need Keqb/t when reversible =true
vals_keqt=[1.]
vals_keqb=[1.]
Keq_top_vals=vals_keqt
Keq_bottom_vals=vals_keqb
Keqb_vals=[1.]#np.logspace(-2,0,3)



xvars=['b_in']
yvars=['b_out']

allvars=['Btot','CKcatB','KeqB','CKcat_bottom','CKcat_top','Keq_top','Keq_bottom','b_in','b_out']




lengthparams_set=([2,1,1],
                  [2,1,1],
                  [2,1,1])

B_position_params_set=(([0,0],[0,0],[1],[1]),
					   ([1,0],[0,0],[1],[1]),
                       ([0,0],[0,0],[0],[1]))

for idstruct, struct in enumerate(lengthparams_set):
    #print(struct)
    #print(B_position_params_set[idstruct])
    
    lengthparams=struct
    B_position_params=B_position_params_set[idstruct]

    print(lengthparams)
    print(B_position_params)    
    print('here1')
    for idxvar, xvar in enumerate(xvars):
        for idyvar, yvar in enumerate(yvars):
            #othervars_list=allvars
            print(xvar,yvar)

            # if idyvar < idxvar:
            print('here2')
            print(xvar,yvar)
            othervars_list=allvars.copy()
            othervars_tuple=()
            print(othervars_list)

            if xvar=='b_in':
                xvarvals=vals_b_in
            
            if yvar=='b_out':
                yvarvals=vals_b_out
            
            othervars_list.remove(xvar)
            othervars_list.remove(yvar)
            print(othervars_list)
            
            for id_other, other in enumerate(othervars_list):
                if other == 'KeqB':
                    add_vals=(Keqb_vals,)
                elif other == 'Keq_top':
                    add_vals=(Keq_top_vals,)
                elif other == 'Keq_bottom':
                    add_vals=(Keq_bottom_vals,)
                elif other == 'CKcatB':
                    add_vals=(vals_CKcB,)
                elif other == 'Btot':
                    add_vals=(vals_Btot,)
                elif other == 'CKcat_bottom':
                    add_vals=(vals_CKcb,)
                elif other == 'CKcat_top':
                    add_vals=(vals_CKct,)
                    
                print(other,add_vals)
                othervars_tuple+=add_vals
            print(othervars_tuple,'\n')
            
            variable_params=othervars_list
            variable_param_values=othervars_tuple
            
            mh.make_heatmaps(xvar,xvarvals,yvar,yvarvals,
                  lengthparams,B_position_params,
                  CKcat,Ks,Kp,Keq,               # params_global
                  KsB,KpB,                             # params_b unchanging
                  Ks_top, Kp_top,
                  Ks_bottom, Kp_bottom,                # params_branch unchanging
                  kout_top,kout_bottom,
                  variable_params,variable_param_values,
                  kin_min,kin_max,num_kin,
                  reversible,Tmax,tpts)  


# make_heatmaps(xvar='Btot',xvarvals=np.logspace(-3,3,25),
#                   yvar='CKcatB',yvarvals=np.logspace(-3,3,25),
#                   lengthparams=[2,2,2],B_position_params=([1,0],[0,0],[0,1],[0,1]),
#                   CKcat=1.,Ks=1.,Kp=1.,Keq=1.,               # params_global
#                   KsB=1.,KpB=1.,                             # params_b unchanging
#                   Ks_top=1., Kp_top=1.,
#                   Ks_bottom=1., Kp_bottom=1.,                # params_branch unchanging
#                   kout_top=1.,kout_bottom=1.,
#                   variable_params=['KeqB',
#                                    'CKcat_top',
#                                    'Keq_top',
#                                    'CKcat_bottom',
#                                    'Keq_bottom'],
#                   variable_param_values=(np.logspace(-2,0,3),
#                                          [1.],
#                                          [1.],
#                                          np.logspace(-1,1,3),
#                                          [1.]),
#                   kin_min=-8,kin_max=0,num_kin=41,
#                   reversible=False,Tmax=1000000,tpts=100000):