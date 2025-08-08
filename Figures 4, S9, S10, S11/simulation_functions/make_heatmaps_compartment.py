#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 14:03:13 2024

@author: robert
"""


import numpy as np

import time
from datetime import timedelta
import multiprocessing as multi
import ode_funcs_compartment as ode_cm
import sys
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#from pathlib import Path
import os
import time
import numpy.random as rn
import numpy as np


sys.path.append('../data_compartment/megadata')
sys.path.append('../data_compartment/heatmapdata')


def _nanargmin(arr):
    """
    Returns the min of an array ignoring nan values
    """
    try:
        return np.nanargmin(arr)
    except ValueError:
       return np.nan

def _nanargmax(arr):
    """
    Returns the max of an array ignoring nan values
    """
    try:
        return np.nanargmax(arr)
    except ValueError:
       return np.nan

def rd_str(x):
    """
    Returns a string from an array of numbers rounded to the nearest 0.001
    """
    return str(np.round(np.log10(x),6))

def make_data_parallel(fn,inp_args):#,totleng):
    #print('this is doing now')
    """
    runs a function in parrallel with arguments from a list using 6 less
    threads the total
    # """
    with multi.Pool(16) as pool:
    # with multi.Pool(26) as pool:
        result=pool.starmap(fn,inp_args)
    # pool.close() 
    # pool.join()
    return result

def findvarindict(key, d):
    return key in d
        
def print_vars_from_dict(d):
    for key, value in d.items():
        if findvarindict(key, d):
            globals()[key] = value
            
def make_time_series(pinit,params_inout, params01,params1E,paramsE2,params23,params35,
         params_24,params_35,
         params_trans,params_diff_3,params_diff_5,
         paramsB,paramsC,allo,Tmax,tpts):
    
    # print(params_global,'\n')
    # print(params_B,'\n')
    # print(params_branch,'\n')
    

    
    odefunc=ode_cm.dPdt

    # Set time points for simulation
    time_points = np.linspace(0, Tmax,tpts)# 10*Tmax)
 #   print(time_points)
   # params.type
    result = odeint(odefunc, pinit, time_points,args=(params_inout, params01,
         params1E,paramsE2,params23,params35,
         params_24,params_35,
         params_trans,params_diff_3,params_diff_5,
         paramsB,paramsC,allo))
    
    # isbuildup=check_buildup(result,odefunc,params_inout, params01,
    #           params1E,paramsE2,params23,params35,
    #           params_24,params_35,
    #           params_trans,params_diff_3,params_diff_5,
    #           paramsB,paramsC,allo)
    return result,time_points,odefunc


def check_buildup(result,odefunc,params_inout, params01,
          params1E,paramsE2,params23,params35,
          params_24,params_35,
          params_trans,params_diff_3,params_diff_5,
          paramsB,paramsC,allo):
    
    concs_final=result[-1,:]            # get final concs
    # print(concs_final)
    
    
   
    
    
    
    grad_final=odefunc(concs_final,0,params_inout, params01,params1E,paramsE2,params23,params35,
         params_24,params_35,
         params_trans,params_diff_3,params_diff_5,
         paramsB,paramsC,allo)        # use them to calc final gradient
    # print(grad_final)
    # print('\n ######## HERE ####### \n')
    # print(grad_final)
    # print(grad_final[:-5])
    # print(grad_final[-1])
    # print('\n ######## HERE ####### \n')
    grad_final_1=grad_final[:-5]     
    grad_final_2=grad_final[-1]
    # print(grad_final_1)
    # print(grad_final_2)

    # print(grad_final)                     # remove the entries for B0 and B1
    
    isbuildup_1=np.array(grad_final_1)>10**(-5)
    isbuildup_2=np.array(grad_final_2)>10**(-5)
    # print(isbuildup_1)
    # print(isbuildup_2)
    
    isbuildup=np.append(isbuildup_1,isbuildup_2)
#    print(grad_final)
 #   print(grad_final/concs_final[:-2])
  #  print(isbuildup)
    return isbuildup

def find_sw_pt(fluxfrac,kinplot,choosepoint=False):
    
    max_orig=max(fluxfrac)
    min_orig=min(fluxfrac)
    
    if choosepoint==False:
        sw_orig=(max_orig+min_orig)/2
    else:
        sw_orig=choosepoint
    
    if max(fluxfrac)>sw_orig:
    
        arg_sw=np.argmax(fluxfrac>sw_orig)
        
        in_min=kinplot[arg_sw-1]
        in_max=kinplot[arg_sw]
        kin_value=(in_min+in_max)/2
        
        sw_min=fluxfrac[arg_sw-1]
        sw_max=fluxfrac[arg_sw]
        
        sw_orig=(sw_min+sw_max)/2
        if max_orig-min_orig>0.1:
            return kin_value,sw_orig
        else: 
            return np.nan, np.nan
        
    else: 
        return np.nan,np.nan
def get_fluxfractions(pinit,params_inout, params01,
          params1E,paramsE2,params23,params35,
          params_24,params_35,
          params_trans,params_diff_3,params_diff_5,
          paramsB,paramsC,allo,Tmax,tpts):
    """
    Return the flux fraction, b fraction and whether there is buildup for a given
    set of parameters
    """
    
 #   print(Tmax,tpts)
    result,timepts,odefunc= make_time_series(np.array(pinit),np.array(params_inout),
             np.array(params01),np.array(params1E),np.array(paramsE2),np.array(params23),
             np.array(params35),
             np.array(params_24),np.array(params_35),
             np.array(params_trans),np.array(params_diff_3),np.array(params_diff_5),
             np.array(paramsB),np.array(paramsC),np.array(allo),
             Tmax,tpts)  
    
    final_concs=result[-1,:]
    
    m2=final_concs[2]
    m3=final_concs[3]
    
    n2=final_concs[5]
    n3=final_concs[6]
    n4=final_concs[8]
    
    c0=final_concs[11]
    c1=final_concs[12]
    
    cosub_concs=[final_concs[9],final_concs[10],c0,c1]
    
    al_p_23=allo[0]
    al_b_23=allo[1]
    al_p_24=allo[2]
    al_b_24=allo[3]
    
    if np.any(final_concs):

        isbuildup=check_buildup(result,odefunc,params_inout, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo)

        flux_top=ode_cm.flux_enz_forward(m2**al_p_23, 1., params23[0], params23[1])
        flux_bottom=ode_cm.flux_enz_forward(n2**al_p_24, c0**al_b_24, params_24[0], params_24[1])
        flux_diff=params_diff_3[1]*n3 
        
    else:
        isbuildup=np.array([np.nan])
        flux_top=np.nan
        flux_bottom=np.nan
        flux_diff=np.nan
        
    return isbuildup, flux_top, flux_bottom, flux_diff, cosub_concs



def get_flux_vecs(pinit,kinvals,params_out, params01,
          params1E,paramsE2,params23,params35,
          params_24,params_35,
          params_trans,params_diff_3,params_diff_5,
          paramsB,paramsC,allo,Tmax,tpts,i):
   
    
    ft_vec=np.zeros(np.size(kinvals))
    ft_vec[:]=np.nan
    
    fb_vec=np.zeros(np.size(kinvals))
    fb_vec[:]=np.nan
    
    fd_vec=np.zeros(np.size(kinvals))
    fd_vec[:]=np.nan
    
    b0frac_vec=np.zeros(np.size(kinvals))
    b0frac_vec[:]=np.nan
    
    c0frac_vec=np.zeros(np.size(kinvals))
    c0frac_vec[:]=np.nan
    
    n0frac_vec=np.zeros(np.size(kinvals))
    n0frac_vec[:]=np.nan
    
    for idk, kin in enumerate(kinvals):
        # print(idk)
        params_inout=[kin,params_out[0],params_out[1]]
        isbuildup, flux_top, flux_bottom, flux_diff, cosub_concs=get_fluxfractions(pinit,params_inout, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts)
        if np.any(isbuildup)==False:
            ft_vec[idk]=flux_top
            fb_vec[idk]=flux_bottom
            fd_vec[idk]=flux_diff
     #   print(cosub_concs)
            b0frac_vec[idk]=cosub_concs[0]/(cosub_concs[0]+cosub_concs[1])
            c0frac_vec[idk]=cosub_concs[2]/(cosub_concs[2]+cosub_concs[3])
            n0frac_vec[idk]=(cosub_concs[0]+cosub_concs[2])/(cosub_concs[0]+cosub_concs[1]+cosub_concs[2]+cosub_concs[3])
    
    ff_vec=fb_vec/(fb_vec+ft_vec)
    
    ff_delta_unsigned=np.nanmax(ff_vec)-np.nanmin(ff_vec)
    arg_min_flux=_nanargmin(ff_vec)
    arg_max_flux=_nanargmax(ff_vec)
    ffsign=np.sign(arg_max_flux- arg_min_flux)
    ff_delta=ffsign*ff_delta_unsigned
    
    b0_delta_unsigned=np.nanmax(b0frac_vec)-np.nanmin(b0frac_vec)
    arg_min_b0=_nanargmin(b0frac_vec)
    arg_max_b0=_nanargmax(b0frac_vec)
    b0sign=np.sign(arg_min_b0-arg_max_b0)
    b0_delta=b0sign*b0_delta_unsigned
    
    c0_delta_unsigned=np.nanmax(c0frac_vec)-np.nanmin(c0frac_vec)
    arg_min_c0=_nanargmin(c0frac_vec)
    arg_max_c0=_nanargmax(c0frac_vec)
    c0sign=np.sign(arg_min_c0-arg_max_c0)
    c0_delta=c0sign*c0_delta_unsigned
    
    n0_delta_unsigned=np.nanmax(n0frac_vec)-np.nanmin(n0frac_vec)
    arg_min_n0=_nanargmin(n0frac_vec)
    arg_max_n0=_nanargmax(n0frac_vec)
    n0sign=np.sign(arg_min_n0-arg_max_n0)
    n0_delta=n0sign*n0_delta_unsigned
    
    
  #  print(ff_delta)
    
    ff_grad=np.gradient(ff_vec,kinvals)    
    
    ff_sens=np.abs(ff_grad/(ff_vec/kinvals))
    
    #if i>60:
  
    # print(~np.isnan(ff_grad))
    # print(np.any(~np.isnan(ff_grad)))
    if np.any(~np.isnan(ff_grad)):
        gradmax = ff_grad[_nanargmax(np.abs(ff_grad))]
        gradmaxind=np.log10(kinvals[_nanargmax(np.abs(ff_grad))])

        sensemax=ff_sens[_nanargmax(np.abs(ff_sens))]
        sensemaxind=np.log10(kinvals[_nanargmax(np.abs(ff_sens))])
        
        swpt,swptval=find_sw_pt(1-ff_vec,kinvals)
    
    else:
        gradmax = np.nan
        gradmaxind=np.nan

        sensemax=np.nan
        sensemaxind=np.nan
        
        swpt=np.nan
        swptval=np.nan
    print(i)
   # if savepath != '':
    #    np.savez(savepath,deltadat=ff_delta,fluxtop=ft_vec,fluxbottom=fb_vec,fluxdiff=fd_vec)
    return ft_vec,fb_vec,fd_vec,ff_vec, ff_grad, ff_sens,ff_delta,gradmax,\
            gradmaxind,sensemax,sensemaxind,b0_delta,c0_delta,n0_delta,i,swpt,swptval
            

def get_fluxfractions_2(pinit,params_inout, params01,
          params1E,paramsE2,params23,params35,
          params_24,params_35,
          params_trans,params_diff_3,params_diff_5,
          paramsB,paramsC,allo,Tmax,tpts):
    """
    Return the flux fraction, b fraction and whether there is buildup for a given
    set of parameters
    """
    
 #   print(Tmax,tpts)
    result,timepts,odefunc= make_time_series(np.array(pinit),np.array(params_inout),
             np.array(params01),np.array(params1E),np.array(paramsE2),np.array(params23),
             np.array(params35),
             np.array(params_24),np.array(params_35),
             np.array(params_trans),np.array(params_diff_3),np.array(params_diff_5),
             np.array(paramsB),np.array(paramsC),np.array(allo),
             Tmax,tpts)  
    
    final_concs=result[-1,:]
    
    m2=final_concs[2]
    m3=final_concs[3]
    
    n2=final_concs[5]
    n3=final_concs[6]
    n4=final_concs[8]
    
    c0=final_concs[11]
    c1=final_concs[12]
    
    cosub_concs=[final_concs[9],final_concs[10],c0,c1]
    
    al_p_23=allo[0]
    al_b_23=allo[1]
    al_p_24=allo[2]
    al_b_24=allo[3]
    
    if np.any(final_concs):

        isbuildup=check_buildup(result,odefunc,params_inout, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo)

        flux_top=ode_cm.flux_enz_forward(m2**al_p_23,  1., params23[0], params23[1])
        flux_bottom=ode_cm.flux_enz_forward(n2**al_p_24,  c0**al_b_24,  params_24[0], params_24[1])
        flux_diff=params_diff_3[1]*n3 
        
    else:
        isbuildup=np.array([np.nan])
        flux_top=np.nan
        flux_bottom=np.nan
        flux_diff=np.nan
        
    return isbuildup, flux_top, flux_bottom, flux_diff, cosub_concs, final_concs



def get_flux_vecs_2(pinit,kinvals,params_out, params01,
          params1E,paramsE2,params23,params35,
          params_24,params_35,
          params_trans,params_diff_3,params_diff_5,
          paramsB,paramsC,allo,Tmax,tpts):
   
    
    ft_vec=np.zeros(np.size(kinvals))
    ft_vec[:]=np.nan
    
    fb_vec=np.zeros(np.size(kinvals))
    fb_vec[:]=np.nan
    
    fd_vec=np.zeros(np.size(kinvals))
    fd_vec[:]=np.nan
    
    b0frac_vec=np.zeros(np.size(kinvals))
    b0frac_vec[:]=np.nan
    
    c0frac_vec=np.zeros(np.size(kinvals))
    c0frac_vec[:]=np.nan
    
    n0frac_vec=np.zeros(np.size(kinvals))
    n0frac_vec[:]=np.nan
    
    concs_vec=np.zeros((len(kinvals),14))
    
    for idk, kin in enumerate(kinvals):
        # print(idk)
        params_inout=[kin,params_out[0],params_out[1]]
        isbuildup, flux_top, flux_bottom, flux_diff, cosub_concs,final_concs=get_fluxfractions_2(pinit,params_inout, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts)
        if np.any(isbuildup)==False:
            ft_vec[idk]=flux_top
            fb_vec[idk]=flux_bottom
            fd_vec[idk]=flux_diff
            b0frac_vec[idk]=cosub_concs[0]/(cosub_concs[0]+cosub_concs[1])
            c0frac_vec[idk]=cosub_concs[2]/(cosub_concs[2]+cosub_concs[3])
            n0frac_vec[idk]=(cosub_concs[0]+cosub_concs[2])/(cosub_concs[0]+cosub_concs[1]+cosub_concs[2]+cosub_concs[3])
            concs_vec[idk,:]=final_concs
    ff_vec=fb_vec/(fb_vec+ft_vec)
    
    ff_delta_unsigned=np.nanmax(ff_vec)-np.nanmin(ff_vec)
    arg_min_flux=_nanargmin(ff_vec)
    arg_max_flux=_nanargmax(ff_vec)
    ffsign=np.sign(arg_max_flux- arg_min_flux)
    ff_delta=ffsign*ff_delta_unsigned
    
    b0_delta_unsigned=np.nanmax(b0frac_vec)-np.nanmin(b0frac_vec)
    arg_min_b0=_nanargmin(b0frac_vec)
    arg_max_b0=_nanargmax(b0frac_vec)
    b0sign=np.sign(arg_min_b0-arg_max_b0)
    b0_delta=b0sign*b0_delta_unsigned
    
    c0_delta_unsigned=np.nanmax(c0frac_vec)-np.nanmin(c0frac_vec)
    arg_min_c0=_nanargmin(c0frac_vec)
    arg_max_c0=_nanargmax(c0frac_vec)
    c0sign=np.sign(arg_min_c0-arg_max_c0)
    c0_delta=c0sign*c0_delta_unsigned
    
    n0_delta_unsigned=np.nanmax(n0frac_vec)-np.nanmin(n0frac_vec)
    arg_min_n0=_nanargmin(n0frac_vec)
    arg_max_n0=_nanargmax(n0frac_vec)
    n0sign=np.sign(arg_min_n0-arg_max_n0)
    n0_delta=n0sign*n0_delta_unsigned
    
    
  #  print(ff_delta)
    
    ff_grad=np.gradient(ff_vec,kinvals)    
    
    ff_sens=np.abs(ff_grad/(ff_vec/kinvals))
    
    #if i>60:
  
    # print(~np.isnan(ff_grad))
    # print(np.any(~np.isnan(ff_grad)))
    if np.any(~np.isnan(ff_grad)):
        gradmax = ff_grad[_nanargmax(np.abs(ff_grad))]
        gradmaxind=np.log10(kinvals[_nanargmax(np.abs(ff_grad))])

        sensemax=ff_sens[_nanargmax(np.abs(ff_sens))]
        sensemaxind=np.log10(kinvals[_nanargmax(np.abs(ff_sens))])
    
    else:
        gradmax = np.nan
        gradmaxind=np.nan

        sensemax=np.nan
        sensemaxind=np.nan
    # print(i)
   # if savepath != '':
    #    np.savez(savepath,deltadat=ff_delta,fluxtop=ft_vec,fluxbottom=fb_vec,fluxdiff=fd_vec)
    return  concs_vec

def make_params_newtest(allo,n,fixed=False,newtots=False):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        print(i)
        numkin=101
        kinvals=np.logspace(-12,0,numkin)
        params_out=[0.01,0.01]



        ########  do between 0.1, 1 rather than 0.1 and 10
        Btot=10**rn.uniform(-1,1)
        Ctot=10**rn.uniform(-1,1)

        # print(Btot,Ctot)
        
        params01=[4*10**rn.uniform(-1,1),10**rn.uniform(-2,0),10**rn.uniform(-3,0),0.0056]   #GAP->BPG
        params1E=[10**rn.uniform(0,1),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]     #BPG->PEP
        paramsE2=[4*10**rn.uniform(-1,1),10**rn.uniform(-2,np.log10(50)),10**rn.uniform(-2,np.log10(50)),6.5*10**3]     #PEP->PYR
        
        params35=[3*10**rn.uniform(-1,1),10**rn.uniform(-2,np.log10(50)),10**rn.uniform(-2,np.log10(50)),1.45*10**4]     #ACET->ETH cyto
        params_35=[3*10**rn.uniform(-1,1),10**rn.uniform(-2,np.log10(50)),10**rn.uniform(-2,np.log10(50)),1.45*10**4]    #ACET->ETH mito
        
        paramsB=[10**rn.uniform(-1,1),10**rn.uniform(-2,np.log10(50)),10**rn.uniform(-2,np.log10(50)),0.1]      #N+->NH cyto
        paramsC=[10**rn.uniform(-1,1),10**rn.uniform(-2,np.log10(50)),10**rn.uniform(-2,np.log10(50)),0.01]      #N+->NH mito
        
        params23=[6*10**rn.uniform(-2,0),43*10**rn.uniform(-2,0)]     #PYR->ACET 
        params_trans=[2*10**rn.uniform(-3,-1),6*10**rn.uniform(-2,0)] #PYR->PYR cyto->mito
        params_24=[10**rn.uniform(-2,0),6*10**rn.uniform(-2,0)]    #PYR->AcCoA
        
        
     
        param_diff=10**rn.uniform(0,2)
        params_diff_3=[param_diff*0.4,param_diff*1.15]
        params_diff_5=[param_diff*0.6,param_diff*1.27]
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

        if newtots==True:
            Btot*=10.
            Ctot*=0.1
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        if fixed==True:
            params_trans=params23
            params_24=params23
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args

def make_params_newtest_2(allo,n,gapkeq=False,newtots=False):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        print(i)
        numkin=101
        kinvals=np.logspace(-12,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals
        # print(kpvals)

        Etotvals=rn.uniform(-5,-1,9)
        Etotvals=10**Etotvals
        # print(Etotvals)

        kcatvals=rn.uniform(1,7,9)
        kcatvals=10**kcatvals
        # print(kcatvals)

        VMvals=Etotvals*kcatvals
        # print(VMvals)
        
        if newtots==True:
            Btot=10**rn.uniform(-0.5,0.5)*3
            tt=10**rn.uniform(-1,1)*3
            Ctot=Btot*0.005
        else:
            Btot=10**rn.uniform(-1,1)*3
            Ctot=10**rn.uniform(-1,1)/5
        
#         if newtots==True:
# 			Btot=10**rn.uniform(-0.5,0.5)*3
# 			tt=10**rn.uniform(-1,1)*3
# 			Ctot=Btot*0.005
# 		else:
# 			Btot=10**rn.uniform(-1,1)*3
# 			Ctot=10**rn.uniform(-1,1)/5

        ########  do between 0.1, 1 rather than 0.1 and 10
        
        # Btot=10**rn.uniform(-1,1)  
        # Ctot=10**rn.uniform(-1,1)

        # print(Btot,Ctot)
        
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-2,-1)/2.,10**rn.uniform(-2,-1)/10.,10**rn.uniform(-1,3)*2.5]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-2,-1)/2.,10**rn.uniform(-2,-1)/10.,10**rn.uniform(-1,3)*2.5]    #BPG->PEP
        paramsE2=[4.05*10**rn.uniform(-4,-2)*1.7,0.014,0.53,2.4*10**4]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[2*10**rn.uniform(-1,1),17.,0.17,1.5*10**3]

        params_24=params23
        params_35=[params35[0],17.,0.17,1.5*10**3]
        
        paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]
        
        params_trans=[params23[0], params23[1]] # vt kt
        param_diff=10**rn.uniform(0,2)
        params_diff_3=[param_diff*0.4,param_diff*1.15]
        params_diff_5=[param_diff*0.6,param_diff*1.27]
        
        
     
        if gapkeq==True:
            params01=[4*10**rn.uniform(-1,1),10**rn.uniform(-2,0),10**rn.uniform(-3,0),0.0056] 

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

        
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
      
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args

def make_params_newtest_2_orig(allo,n,gapkeq=False):#,largercosubtot=False):
    
    ####################################################################
    #### Btot, Ctot and transport are the same as the fixed case #######
    #### paramsB and paramsC are the same as the fixed case ############
    ####################################################################
    
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        print(i)
        numkin=101
        kinvals=np.logspace(-12,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals
        # print(kpvals)

        Etotvals=rn.uniform(-5,-1,9)
        Etotvals=10**Etotvals
        # print(Etotvals)

        kcatvals=rn.uniform(1,7,9)
        kcatvals=10**kcatvals
        # print(kcatvals)

        VMvals=Etotvals*kcatvals
        # print(VMvals)
        
        Btot=10**rn.uniform(-1,1)*3
        Ctot=10**rn.uniform(-1,1)/5

        ########  do between 0.1, 1 rather than 0.1 and 10
        
        # Btot=10**rn.uniform(-1,1)  
        # Ctot=10**rn.uniform(-1,1)

        # print(Btot,Ctot)
        
        params01=[10**rn.uniform(-3,-2),ksvals[0],kpvals[0],10**rn.uniform(-1,3)] #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),ksvals[0],kpvals[0],10**rn.uniform(-1,3)]   #BPG->PEP
        paramsE2=[4.05*10**rn.uniform(-4,-2),0.014,0.53,2.4*10**4]  #PEP->PYR
        
        params23=[0.65*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-1,1),17.,0.17,1.5*10**3]

        params_24=[VMvals[8],0.5]#,0.5,10**6.]
        params_35=[params35[0],17.,0.17,1.5*10**3]
        
        paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]
        
        params_trans=10**np.array([rn.uniform(-5,-2),rn.uniform(-4,-1)])# vt kt
        
        param_diff=10**rn.uniform(0,2)
        params_diff_3=[param_diff*0.4,param_diff*1.15]
        params_diff_5=[param_diff*0.6,param_diff*1.27]
        
        if gapkeq==True:
            params01=[4*10**rn.uniform(-1,1),10**rn.uniform(-2,0),10**rn.uniform(-3,0),0.0056] 
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
      
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args

def make_params_newtest_3(allo,n):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals

        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=params23
        params_35=params35
        
        paramsB=[10**rn.uniform(-5,-3)*6 ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1]
        paramsC=[10**rn.uniform(-5,-3)*4 ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01]
        
        params_trans=[params23[0], params23[1]] # vt kt
        param_diff=10**rn.uniform(0,2)
        params_diff_3=[param_diff*0.4,param_diff*1.15]
        params_diff_5=[param_diff*0.6,param_diff*1.27]
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args

def make_params_newtest_3_nonfixed(allo,n):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        # print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals

        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=[10**rn.uniform(-5,-3),6*10**rn.uniform(-2,-0)] 
        params_35=params35
        
        paramsB=[10**rn.uniform(-5,-3)*6 ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1]
        paramsC=[10**rn.uniform(-5,-3)*4 ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01]
        
        params_trans=[2*10**rn.uniform(-6,-4),6*10**rn.uniform(-2,-0)]# vt kt
        param_diff=10**rn.uniform(0,2)
        params_diff_3=[param_diff*0.4,param_diff*1.15]
        params_diff_5=[param_diff*0.6,param_diff*1.27]
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args


def make_params_newtest_4(allo,n):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals

        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=params23
        params_35=params35
        
        paramsB=[10**rn.uniform(-5,-3)*6 ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1]
        paramsC=[10**rn.uniform(-5,-3)*4 ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01]
        
        params_trans=[params23[0], params23[1]] # vt kt
        param_diff=10**rn.uniform(0,2)
        params_diff_5=[param_diff*0.4,param_diff*1.15]
        params_diff_3=[param_diff*0.6,param_diff*1.27]
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args


def make_params_newtest_4_nonfixed(allo,n):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        # print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals

        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=[10**rn.uniform(-5,-3),6*10**rn.uniform(-2,-0)] 
        params_35=params35
        
        paramsB=[10**rn.uniform(-5,-3)*6 ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1]
        paramsC=[10**rn.uniform(-5,-3)*4 ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01]
        
        params_trans=[2*10**rn.uniform(-6,-4),6*10**rn.uniform(-2,-0)]# vt kt
        param_diff=10**rn.uniform(0,2)
        params_diff_5=[param_diff*0.4,param_diff*1.15]
        params_diff_3=[param_diff*0.6,param_diff*1.27]
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args


def make_params_newtest_5(allo,n):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals

        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=params23
        params_35=params35
        
        paramsB=[10**rn.uniform(-5,-3)*6 ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1]
        paramsC=[10**rn.uniform(-5,-3)*4 ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01]
        
        params_trans=[params23[0], params23[1]] # vt kt
        param_diff=10**rn.uniform(-2,0)
        params_diff_5=[param_diff*0.4,param_diff*1.15]
        params_diff_3=[param_diff*0.6,param_diff*1.27]
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args


def make_params_newtest_5_nonfixed(allo,n):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        # print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-3,1,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-3,1,8)
        kpvals=10**kpvals

        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=[10**rn.uniform(-5,-3),6*10**rn.uniform(-2,-0)] 
        params_35=params35
        
        paramsB=[10**rn.uniform(-5,-3)*6 ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1]
        paramsC=[10**rn.uniform(-5,-3)*4 ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01]
        
        params_trans=[2*10**rn.uniform(-6,-4),6*10**rn.uniform(-2,-0)]# vt kt
        param_diff=10**rn.uniform(-2,0)
        params_diff_5=[param_diff*0.4,param_diff*1.15]
        params_diff_3=[param_diff*0.6,param_diff*1.27]
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args


def make_params_newtest_6(allo,n):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-4,0,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-4,0,8)
        kpvals=10**kpvals
        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=params23
        params_35=params35
        
        paramsB=[10**rn.uniform(-5,-3)*6 ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1]
        paramsC=[10**rn.uniform(-5,-3)*4 ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01]
        
        params_trans=[params23[0], params23[1]] # vt kt
        dummy=10**rn.uniform(-4,0)
        # param_diff2=10**rn.uniform(-4,0)
       
        param_diff1=ksvals[0]
        param_diff2=kpvals[0]
        params_diff_5=[param_diff1,param_diff1*2.]
        params_diff_3=[param_diff2,param_diff2*2.]
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args

def make_params_newtest_6_nonfixed(allo,n,VBmod=1.,VCmod=1.,KBmod=1.,KCmod=1.):#,largercosubtot=False):
    
    list_of_args=list()
   # print(fixed,fasttrans,fastdiff,smallercosubs,fasterbackground,asymdiff)
    for i in range(n):
        # print(i)
        numkin=81
        kinvals=np.logspace(-8,0,numkin)
        params_out=[0.01,0.01]

        ksvals=rn.uniform(-4,0,8)
        ksvals=10**ksvals
      #  print(ksvals)

        kpvals=rn.uniform(-4,0,8)
        kpvals=10**kpvals

        ########  do between 0.1, 1 rather than 0.1 and 10
    
        Btot=10**rn.uniform(-0.5,0.5)
        # tt=10**rn.uniform(-1,1)*3
        Ctot=0.1*10**rn.uniform(-0.5,0.5)
   

        # print(Btot,Ctot)
        params01=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),0.0056]  #GAP->BPG
        params1E=[10**rn.uniform(-3,-2),10**rn.uniform(-3,0),10**rn.uniform(-3,0),10**rn.uniform(-1,1)]    #BPG->PEP
        paramsE2=[4.0*10**rn.uniform(-4,-2),0.014,0.53,6.5*10**3]   #PEP->PYR
        
        params23=[3*10**rn.uniform(-5,-3),5.]#,5.,2.4*10**2]
        params35=[3*10**rn.uniform(-5,-3),17.,0.17,1.45*10**4]
 
        params_24=[10**rn.uniform(-5,-3),6*10**rn.uniform(-2,-0)] 
        params_35=params35
        
        paramsB=[(10**rn.uniform(-5,-3)*6)*VBmod ,10**rn.uniform(-2,1),10**rn.uniform(-2,1),0.1*KBmod]
        paramsC=[(10**rn.uniform(-5,-3)*4)*VCmod ,10**rn.uniform(-1,2),10**rn.uniform(-3,-0),0.01*KCmod]
        
        params_trans=[2*10**rn.uniform(-6,-4),6*10**rn.uniform(-2,-0)]# vt kt
        
        dummy=10**rn.uniform(-4,0)
        # param_diff2=10**rn.uniform(-4,0)
        
        param_diff1=ksvals[0]
        param_diff2=kpvals[0]
        params_diff_5=[param_diff1,param_diff1*2.]
        params_diff_3=[param_diff2,param_diff2*2.]
     
        

        # paramsB=[10**rn.uniform(-4,-2)*0.6 ,ksvals[6]*0.6,kpvals[6]*0.9,0.1]
        # paramsC=[10**rn.uniform(-4,-2)*0.4 ,ksvals[7]*44,kpvals[7]*0.06,0.01]

       

        Tmax=100000.
        tpts=int(10*Tmax)

        # make all 0.001        


        pinit=0.001*np.ones(14)    
        pinit[0]*=10

       
          
        pinit[-2]=Ctot/2.
        pinit[-3]=Ctot/2.
        pinit[-4]=Btot/2.
        pinit[-5]=Btot/2.
        
        
        
       
            
       
            
       

        
        list_of_args.append(tuple((pinit,kinvals,params_out, params01,
                  params1E,paramsE2,params23,params35,
                  params_24,params_35,
                  params_trans,params_diff_3,params_diff_5,
                  paramsB,paramsC,allo,Tmax,tpts,i)))
        
    return list_of_args