#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 13:11:43 2024

@author: robert
"""

import numpy as np
import numba 
from scipy.integrate import odeint
import matplotlib.pyplot as plt

import checkinputs as check
import reaction_funcs as reac


def make_ode_vectors(lengthparams=[3,2,2],B_position_params=([0,1,0],[0,0,0],[1,0],[1,0]),
              params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,0.1],
              params_branch=[1.,1.,1.,1.,
                             1.,1.,1.,1.],
              params_inout=[0.01,1.,1.],
              reversible=False):
    """
    vectors to pass to the odefunction for making the pdes
    
    ############################################################################
                                    INPUTS
    ############################################################################
    
    lengthparams=[upstream_length,topbranch_length,bottombranch_length]
    
    length params in in terms of reactions, so there are n+1 metabolites
    
    upstream_length > 2
    topbranch_length > 1
    bottombranch_length > 1
    
    
    ############################################################################
    
    B_position_params=[upstream_b0_positions,upstream_b1_positions,
                       topbranch_b1_position,bottombranch_b0_position]
    
    position of b0 or b1 mediated reactions, given as a true/false vector
    for the upstream ones, cannot have 1s in the same index
    length must be the same relevant length param
    must have at least one B reaction in each branch
    
    
    ############################################################################
    
    params_global= [CKcat,Ks,Kp,Keq]
    
    sets the global kinetic parameters for metabolites not on the branch or 
    background conversion. 
    
    if reversible = False then only first two are used
    
    set to 1 arbitrarily
    
    
    ############################################################################
    
    params_B=[CKcatB,KsB,KpB,KeqB,Btot]
    
    sets the kinetic parameters for background conversion of B (B0->B1), as 
    well as the total amount of B. this is always a reversible reaction. 
    Typically KsB and KpB are 1, and the others are changed
    
    
    ############################################################################
    
    params_branch=[CKcat_top,   Ks_top,   Kp_top,   Keq_top,
                   CKcat_bottom,Ks_bottom,Kp_bottom,Keq_bottom]
    
    sets the kinetic parameters for the reactions at the branch point, typically
    Ks and Kp are set to 1, and the others are changed.
    
    
    ############################################################################
    
    params_inout=[kin, kout_top, kout_bottom]
    
    define the influx and outflux parameters
    
    
    ############################################################################
    
    reversible=true/false
    
    define the metabolite reactions as reversible or irreversible
    
    
    ############################################################################
    
    ############################################################################
    ############################################################################
    ############################################################################
    
    ############################################################################
                                    OUTPUTS
    ############################################################################
    
    Pinitial   -   initial concentrations 
    Pinitial=[upstream_concs, topbranch_concs, bottombranch_concs,Bconcs]
    
    
    ############################################################################
    
    PDE_func     -   function for simulating the pdes
    
    
    ############################################################################
    
    ############################################################################
    ############################################################################
    
    """
   
    # First check for errors in inputs
    failed_lengthparams=check.check_length_params(lengthparams)
    if failed_lengthparams!=0:
        print('check lengthparams')
        return
    
    # Check B positions
    failed_B_position_params=check.check_B_position_params(lengthparams,B_position_params)
    if failed_B_position_params!=0:
        print('check B_position_params')
        return
    
    # Check parameters
    failed_params=check.check_params(params_global,params_B,params_branch)
    if failed_params!=0:
        print('check parameters')
        return
    
    # failed_reac_type=check.check_reaction_type(reversible)
    # if failed_reac_type!=0:
    #     print('reversible reactions are not yet implemented')
    #     return
    #make the initial concs vector
    Pinitial=reac.make_initial_concs(lengthparams,params_B)
   # print(Pinitial)
    
    # make the kin vector
    kinvec=reac.make_kin_vec(lengthparams)
    # print(kinvec)
    
    #make the kout vector
    kouttopvec, koutbottomvec=reac.make_kout_vec(lengthparams)
    # print(kouttopvec)
    # print(koutbottomvec)
    
    MM_mat_noB=reac.make_MM_mat_noB(lengthparams,B_position_params)
    # print(MM_mat_noB)
    
    MM_mat_B01=reac.make_MM_mat_B01(lengthparams,B_position_params)
    # print(MM_mat_B01)
    
    MM_mat_B10=reac.make_MM_mat_B10(lengthparams,B_position_params)
    # print(MM_mat_B10)

    
    return Pinitial,kinvec,kouttopvec,koutbottomvec,MM_mat_noB,MM_mat_B01,MM_mat_B10

def make_ode_func(lengthparams=[3,2,2],B_position_params=([0,1,0],[0,0,0],[1,0],[1,0]),
            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,0.1],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.01,1.,1.],
            reversible=False):  
    
    
    Pinitial,kinvec,kouttopvec,koutbottomvec,MM_mat_noB,MM_mat_B01,MM_mat_B10=make_ode_vectors(lengthparams,B_position_params,
              params_global,params_B,
              params_branch,
              params_inout,
              reversible)
    
    
    params=np.concatenate((params_inout, params_global, params_branch,params_B))
    # kin =             params[0]
    # kouttop =         params[1]
    # koutbottom =      params[2]
    # CKcat_gobal =     params[3]
    # Ks_gobal =        params[4]
    # Kp_gobal =        params[5]
    # Keq_gobal =       params[6]
    # CKcat_top =       params[7]
    # Ks_top =          params[8]
    # Kp_top =          params[9]
    # Keq_top =         params[10]
    # CKcat_bottom =    params[11]
    # Ks_bottom =       params[12]
    # Kp_bottom =       params[13]
    # Keq_bottom =      params[14]
    # CKcatB =          params[15]
    # KsB =             params[16]
    # KpB =             params[17]
    # KeqB =            params[18]
    # Btot =            params[19]
    # 
    
    def odefunc(P,t,params):
        
        rate_kin=reac.influx(kinvec,params)
        # print(len(rate_kin))
        
        rate_kouttop=reac.kouttoprate(P,params,kouttopvec)
        # print(len(rate_kouttop))
        
        rate_koutbottom=reac.koutbottomrate(P,params,koutbottomvec)
        # print(len(rate_koutbottom))
        
        rate_noB=reac.noBrate(P,params,lengthparams,MM_mat_noB,reversible)
        # print((rate_noB))
        
        rate_B01=reac.B01rate(P,params,lengthparams,MM_mat_B01,reversible)
        # print(len(rate_B01))
        
        rate_B10=reac.B10rate(P,params,lengthparams,MM_mat_B10,reversible)
        # print(len(rate_B10))
        
        rate_Bback=reac.Bback(P,params,lengthparams)
        # print(len(rate_Bback))

        rate_total = rate_kin + rate_kouttop + rate_koutbottom +rate_noB + rate_B01 + rate_B10 + rate_Bback
        return rate_total

    return Pinitial, odefunc, params
    # if reversible==False:
    #     rate_kinvec=
    # else:
    #     print('irreversible is not yet inplemented')
    #     return
    #return odefunc,Pinitial

#@numba.njit  
def make_time_series(lengthparams=[3,2,2],B_position_params=([0,0,0],[0,0,0],[1,0],[1,0]),
            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,10.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.1,1.,1.],
            reversible=False,Tmax=1000,tpts=100):
    
    # print(params_global,'\n')
    # print(params_B,'\n')
    # print(params_branch,'\n')
    
    
    pinit,odefunc,params=make_ode_func(lengthparams,B_position_params,
              params_global,params_B,
              params_branch,
              params_inout,
               reversible)
    

    # Set time points for simulation
    time_points = np.linspace(0, Tmax,tpts)# 10*Tmax)

    # params.type
    result = odeint(odefunc, pinit, time_points,args=(params,))
    #check_buildup(result,odefunc,params,reversible)
    return result,time_points,odefunc,params

def make_time_series_pinit(lengthparams=[3,2,2],B_position_params=([0,0,0],[0,0,0],[1,0],[1,0]),
            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,10.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            params_inout=[0.1,1.,1.],
            reversible=False,Tmax=1000,tpts=100,p_init=np.zeros(5),usepinit=False):
    
    # print(params_global,'\n')
    # print(params_B,'\n')
    # print(params_branch,'\n')
    
    
    pinit,odefunc,params=make_ode_func(lengthparams,B_position_params,
              params_global,params_B,
              params_branch,
              params_inout,
               reversible)
    
    if usepinit==True:
        pinit=p_init
    # Set time points for simulation
    time_points = np.linspace(0, Tmax,tpts)# 10*Tmax)

    # params.type
    result = odeint(odefunc, pinit, time_points,args=(params,))
    #check_buildup(result,odefunc,params,reversible)
    return result,time_points,odefunc,params

def check_buildup(result,odefunc,params,reversible):
    
    concs_final=result[-1,:]            # get final concs
    # print(concs_final)
    grad_final=odefunc(concs_final,0,params)        # use them to calc final gradient
    # print(grad_final)
    grad_final=grad_final[:-2]     
    # print(grad_final)                     # remove the entries for B0 and B1
    
    isbuildup=np.array(grad_final)>10**(-5)
    
#    print(grad_final)
 #   print(grad_final/concs_final[:-2])
  #  print(isbuildup)
    return isbuildup

def get_branch_fluxes(result,odefunc,params,lengthparams,B_position_params,reversible):
    concs_final=result[-1,:]
    
    Mbranch_final=concs_final[lengthparams[0]]
    Mtop_final=concs_final[lengthparams[0]+1]
    Mbottom_final=concs_final[lengthparams[0]+lengthparams[1]+1]
    B0_final=concs_final[-2]
    B1_final=concs_final[-1]
    
    bfrac=B0_final/(B0_final+B1_final)
    brat=B0_final/B1_final
    topbranch10=B_position_params[2]
    bottombranch01=B_position_params[3]
    
    CKcat_top =       params[7]
    Ks_top =          params[8]
    Kp_top =          params[9]
    Keq_top =         params[10]
    CKcat_bottom =    params[11]
    Ks_bottom =       params[12]
    Kp_bottom =       params[13]
    Keq_bottom =      params[14]
    
    if reversible == False:
        if topbranch10[0]==0:
            flux_top=reac.mu(Mbranch_final, Mtop_final, 1., 1., CKcat_top, Ks_top, Kp_top, Keq_top)
        else:
            flux_top=reac.mu(Mbranch_final, Mtop_final,B1_final , B0_final, CKcat_top, Ks_top, Kp_top, Keq_top)
            
        if bottombranch01[0]==0:
            flux_bottom=reac.mu(Mbranch_final, Mbottom_final, 1., 1., CKcat_bottom, Ks_bottom, Kp_bottom, Keq_bottom)
        else:
            flux_bottom=reac.mu(Mbranch_final, Mbottom_final, B0_final, B1_final, CKcat_bottom, Ks_bottom, Kp_bottom, Keq_bottom)
    else:
        if topbranch10[0]==0:
            flux_top=reac.forwardflux_rev(Mbranch_final, Mtop_final, 1., 1., CKcat_top, Ks_top, Kp_top, Keq_top)
        else:
            flux_top=reac.forwardflux_rev(Mbranch_final, Mtop_final,B1_final , B0_final, CKcat_top, Ks_top, Kp_top, Keq_top)
            
        if bottombranch01[0]==0:
            flux_bottom=reac.forwardflux_rev(Mbranch_final, Mbottom_final, 1., 1., CKcat_bottom, Ks_bottom, Kp_bottom, Keq_bottom)
        else:
            flux_bottom=reac.forwardflux_rev(Mbranch_final, Mbottom_final, B0_final, B1_final, CKcat_bottom, Ks_bottom, Kp_bottom, Keq_bottom)
     
    return flux_top,flux_bottom,bfrac,brat
# print('\n######################################################\n')
# print('3, 2, 2')

# pinit,odefunc,params=make_ode_func(lengthparams=[3,2,2],B_position_params=([0,1,0],[0,0,0],[1,0],[1,0]),
#               params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,0.1],
#               params_branch=[1.,1.,1.,1.,
#                              1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#                reversible=False)

# print(pinit)
# print(odefunc)


# print('\n######################################################\n')
# print('3, 3, 3')

# make_ode_func(lengthparams=[3,3,3],B_position_params=([0,1,0],[0,0,0],[1,0,0],[1,0,0]),
#               params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,0.1],
#               params_branch=[1.,1.,1.,1.,
#                              1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#               reversible=False)

# print('\n######################################################\n')
# print('3, 2, 2')

# make_ode_func(lengthparams=[3,2,2],B_position_params=([1,0,0],[0,0,0],[0,1],[0,1]),
#               params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,0.1],
#               params_branch=[1.,1.,1.,1.,
#                              1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#               reversible=False)

# print('\n######################################################\n')
# print('2, 2, 2')


# make_ode_func(lengthparams=[2,2,2],B_position_params=([1,0],[0,0],[0,1],[0,1]),
#               params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,0.1],
#               params_branch=[1.,1.,1.,1.,
#                              1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#               reversible=False)

# print('\n######################################################\n')
# print('2,1,1')


# make_ode_func(lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
#               params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,0.1],
#               params_branch=[1.,1.,1.,1.,
#                              1.,1.,1.,1.],
#               params_inout=[0.01,1.,1.],
#               reversible=False)
   
   
   
   
   