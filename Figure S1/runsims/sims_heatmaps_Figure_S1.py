#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 16:08:57 2024

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


#lengthparams_set=([2,1,1],
#                  [2,2,2],
#                  [2,1,1],
#                  [2,2,2])
                 # [2,2,2])

#B_position_params_set=(([1,0],[0,0],[1],[1]),
#                       ([1,0],[0,0],[0,1],[1,0]),
#                       ([0,0],[0,0],[1],[1]),
#                       ([0,0],[0,0],[0,1],[1,0]))

# lengthparams_set=([3,1,1],
# 				  [3,1,1],
# 				  [2,2,2],
# 				  [2,2,2],
# 				  [2,2,2],
# 				  [3,2,2],
# 				  [2,3,2])

# B_position_params_set=(([1,0,0],[0,0,0],[1],[1]),
# 					   ([1,1,0],[0,0,0],[1],[1]),
# 					   ([1,0],[0,0],[0,1],[0,1]),
# 					   ([0,1],[0,0],[0,1],[1,0]),
# 					   ([1,1],[0,0],[0,1],[1,0]),
# 					   ([1,0,1],[0,0,0],[0,1],[1,0]),
# 					   ([1,0],[0,0],[0,0,1],[1,0]))
lengthparams_set=([2,1,1],
                  [2,1,1],
                  [2,1,1])

B_position_params_set=(([0,0],[0,0],[1],[1]),
					   ([1,0],[0,0],[1],[1]),
                       ([0,0],[0,0],[0],[1]))

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
#vals_CKcB_hm=np.logspace(-3,3,25)
#vals_Bt_hm=np.logspace(-3,3,25)
# vals_CKct_hm=np.logspace(-3,3,25)
# vals_CKcb_hm=np.logspace(-3,3,25)

#vals_b_in=np.logspace(-5,2,8)#np.logspace(-5,2,29)
#vals_b_out=np.logspace(-5,2,8)#np.logspace(-5,2,29)
vals_b_in=np.logspace(-5,2,29)
vals_b_out=np.logspace(-5,2,29)

#vals_CKcB=np.logspace(-1,1,3)
#vals_Bt=np.logspace(-1,1,3)
#vals_CKct=np.logspace(-1,1,3)
vals_CKcB=[0.01]#np.logspace(-1,1,3)
vals_Btot=[10.]#np.logspace(-2,0,3)
vals_CKct=[1.]
vals_CKcb=[1.]

#vals_CKcB_hm_all=(vals_CKcB_hm,vals_CKcB)
#vals_Bt_hm_all=(vals_Bt_hm,vals_Bt)
#vals_CKct_hm_all=(vals_CKct_hm,vals_CKct)
#vals_CKcb_hm_all=(vals_CKcb_hm,vals_CKcb)

### ONLY need Keqb/t when reversible =true
vals_keqt=[1.]
vals_keqb=[1.]
Keq_top_vals=vals_keqt
Keq_bottom_vals=vals_keqb
Keqb_vals=[1.]#np.logspace(-2,0,3)


perms = itertools.permutations([1, 1])

# Convert the permutations to a set to remove duplicates, and then convert back to a list
unique_perms = list(set(perms))

print(unique_perms)







t0=time.time()
t0f=time.time()

for idstruct, struct in enumerate(lengthparams_set):
   # print(struct)
   # print(B_position_params_set[idstruct])

    lengthparams=struct
    B_position_params=B_position_params_set[idstruct]

    print(lengthparams)
    print(B_position_params)

    for perm in unique_perms:

        print(perm)
        CKcatB_vals=vals_CKcB#[perm[0]]
        Btot_vals=vals_Btot#[perm[1]]
        CKcat_top_vals=vals_CKct#_hm_all[perm[2]]
        CKcat_bottom_vals=vals_CKcb#_hm_all[perm[3]]



       # print(len(CKcatB_vals),len(Keqb_vals),len(Btot_vals),
           #   len(CKcat_top_vals),len(CKcat_bottom_vals))

        t0=time.time()
        params_to_add=mh.Gen_all_data(lengthparams,B_position_params,
                                      CKcat,Ks,Kp,Keq,               # params_global
                                      KsB,KpB,                             # params_b unchanging
                                      Ks_top, Kp_top,
                                      Ks_bottom, Kp_bottom,                # params_branch unchanging
                                      kout_top,kout_bottom,                # params_inout unchanging
                                      CKcatB_vals,Keqb_vals,Btot_vals,             # params_b variable
                                      CKcat_top_vals,Keq_top_vals,
                                      CKcat_bottom_vals,Keq_bottom_vals,       # params_branch variable
                                      kin_min,kin_max,num_kin,           # params_inout variable
                                      vals_b_in,vals_b_out,
                                      reversible,Tmax,tpts)


      #  print(params_to_add)
        t1=time.time()

        timetaken=timedelta(seconds=t1-t0)
        print(f'making list took {timetaken} seconds')


        t0=time.time()
        if len(params_to_add)>0:
            mh.make_data_parallel(mh.get_ff_delta,params_to_add)
        t1=time.time()

        timetaken=timedelta(seconds=t1-t0)
        print(f'simulating data took {timetaken} seconds')

t1f=time.time()
timetakenf=timedelta(seconds=t1f-t0f)

print(f'simulating all data took {timetakenf} seconds')
# mh.Gen_all_data(lengthparams,B_position_params,
#                   CKcat,Ks,Kp,Keq,               # params_global
#                   KsB,KpB,                             # params_b unchanging
#                   Ks_top, Kp_top,
#                   Ks_bottom, Kp_bottom,                # params_branch unchanging
#                   kout_top,kout_bottom,                # params_inout unchanging
#                   CKcatB_vals,Keqb_vals,Btot_vals,             # params_b variable
#                   CKcat_top_vals,Keq_top_vals,
#                   CKcat_bottom_vals,Keq_bottom_vals,       # params_branch variable
#                   kin_min,kin_max,num_kin,           # params_inout variable
#                   reversible,Tmax,tpts)



