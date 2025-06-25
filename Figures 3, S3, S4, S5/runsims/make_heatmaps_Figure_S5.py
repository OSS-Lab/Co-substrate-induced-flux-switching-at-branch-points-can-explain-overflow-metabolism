#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 14:47:00 2025

@author: robert
"""


import numpy as np

import time
import sys



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

kin_min=-10
kin_max=0
num_kin=201

Tmax=1000000
tpts=100000


reversible=False
vals_CKcB_hm=np.logspace(-3,3,25)
vals_Bt_hm=np.logspace(-3,3,25)
# vals_CKct_hm=np.logspace(-3,3,25)
# vals_CKcb_hm=np.logspace(-3,3,25)


# vals_CKcB=np.logspace(-1,1,3)
vals_CKct=[0.1]#np.logspace(-1,1,3)
vals_CKcb=[1.]#np.logspace(-1,1,3)

# vals_CKcB_hm_all=(vals_CKcB_hm,vals_CKcB)
# vals_Bt_hm_all=(vals_Bt_hm,vals_Bt)
# vals_CKct_hm_all=(vals_CKct_hm,vals_CKct)
# vals_CKcb_hm_all=(vals_CKcb_hm,vals_CKcb)

### ONLY need Keqb/t when reversible =true
vals_keqt=[1.]#np.logspace(0,4,4)#(0,2,3)
vals_keqb=[1.]#np.logspace(0,4,4)#(0,2,3)
Keq_top_vals=vals_keqt
Keq_bottom_vals=vals_keqb
Keqb_vals=[1.]#np.logspace(-2,0,3)



xvars=['Btot','CKcatB',]
yvars=['Btot','CKcatB',]

allvars=['Btot','CKcatB','KeqB','CKcat_bottom','CKcat_top','Keq_top','Keq_bottom']

lengthparams_set=([2,2,2],
                  [2,2,2],
				  [3,1,1],
                  [3,1,1])

B_position_params_set=(([0,0],[0,0],[0,1],[1,0]),
                       ([1,0],[0,0],[0,1],[1,0]),
					   ([0,0,0],[0,0,0],[1],[1]),
                       ([0,1,0],[0,0,0],[1],[1]))

for idstruct, struct in enumerate(lengthparams_set):
    #print(struct)
    #print(B_position_params_set[idstruct])
    
    lengthparams=struct
    B_position_params=B_position_params_set[idstruct]

    print(lengthparams)
    print(B_position_params)    
    for idxvar, xvar in enumerate(xvars):
        for idyvar, yvar in enumerate(yvars):
            
            #othervars_list=allvars
            print(xvar,yvar)
            if idyvar > idxvar:
                print(xvar,yvar)
                othervars_list=allvars.copy()
                othervars_tuple=()
                print(othervars_list)

             #   print(allvars)
                if xvar=='KeqB':
                    xvarvals=np.logspace(-3,0,25)
                else:
                    xvarvals=np.logspace(-3,3,25)
                    
                if yvar=='KeqB':
                    yvarvals=np.logspace(-3,0,25)
                else:
                    yvarvals=np.logspace(-3,3,25)
                
           #     print(xvar, xvarvals,'\n')
            #    print(yvar, yvarvals,'\n')
                
             #   print(othervars_list)
                othervars_list.remove(xvar)
                othervars_list.remove(yvar)
                print(othervars_list)
                
                for id_other, other in enumerate(othervars_list):
                    if other == 'KeqB':
                        # add_vals=(np.logspace(-2,0,3),)
                        add_vals=(Keqb_vals,)
                    elif other == 'Keq_top':
                        # add_vals=([1.],)
                        add_vals=(Keq_top_vals,)
                    elif other == 'Keq_bottom':
                        # add_vals=([1.],)
                        add_vals=(Keq_bottom_vals,)
                    elif other =='CKcat_bottom':
                        add_vals=(vals_CKcb,)
                    elif other =='CKcat_top':
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