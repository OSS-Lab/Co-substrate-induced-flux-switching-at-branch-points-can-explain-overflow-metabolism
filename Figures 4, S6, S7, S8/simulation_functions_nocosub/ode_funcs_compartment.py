#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 12:01:46 2024

@author: robert
"""

import numba
import numpy as np


@numba.njit()
def mu(sub, cosub, CKcat, Ks):
    
    s=sub*cosub/Ks
    
    return CKcat  * s/(1.  + s)

@numba.njit()
def mu_rev(sub, prod, cosub, coprod, CKcat, Ks, Kp, Keq):
    
    s=sub*cosub/Ks
    p=prod*coprod/Kp
    
    
    return (CKcat/Ks)*(sub*cosub-prod*coprod/Keq)/(1.+p+s)
    #rat=prod*coprod/(sub*cosub*Keq)
   
    #return CKcat * (1. - rat) * s/(1. + p + s)

@numba.njit()
def transp(sub,v,k):
    return v*sub/(k+sub)

@numba.njit()
def diffn(sub,prod,ds,df):
    return ds*sub - df*prod

@numba.njit()
def flux_enz_forward(sub, cosub, CKcat, Ks):
    s=sub*cosub/Ks
    return (CKcat*s)/(1. + s)

# P=[m0,  #0
#    m1,  #1
#    m2,  #2
#    m3,  #3
#    m5,  #4
#    n2,  #5
#    n3,  #6
#    n5,  #7
#    n4,  #8
#    n6,  #9
#    b0,  #10
#    b1,  #11
#    c0,  #12
#    c1]  #13

@numba.njit()
def dm0dt(P,params_inout,params01,t):
    m0=P[0]
    m1=P[1]
    # b0=P[9]
    # b1=P[10]
    
    
    kin=params_inout[0]
    VM01,Ks01,Kp01,Keq01=params01
    
    # print(m0,b0,Keq01)
    # print(m0*b0*Keq01)
    
   # print(1+m0*b0/Ks01+m1*b1/Kp01)
    return kin-mu_rev(m0,m1,1.,1.,VM01,Ks01,Kp01,Keq01)

@numba.njit()
def dm1dt(P,params01,params1E,t):
    m0=P[0]
    m1=P[1]
    mE=P[13]
    # b0=P[9]
    # b1=P[10]
    
    VM01,Ks01,Kp01,Keq01=params01
    VM1E,Ks1E,Kp1E,Keq1E=params1E

    return mu_rev(m0,m1,1.,1.,VM01,Ks01,Kp01,Keq01)-\
           mu_rev(m1,mE,1.,1.,VM1E,Ks1E,Kp1E,Keq1E)
           
@numba.njit()
def dmEdt(P,params1E,paramsE2,t):
    m1=P[1]
    m2=P[2]
    mE=P[13]
    # b0=P[9]
    # b1=P[10]
    
    VM1E,Ks1E,Kp1E,Keq1E=params1E
    VME2,KsE2,KpE2,KeqE2=paramsE2

    return mu_rev(m1,mE,1.,1.,VM1E,Ks1E,Kp1E,Keq1E)-\
           mu_rev(mE,m2,1.,1.,VME2,KsE2,KpE2,KeqE2)

@numba.njit()           
def dm2dt(P,paramsE2,params23,params_trans,allo,t):
    mE=P[13]
    m2=P[2]
    m3=P[3]
    
    VME2,KsE2,KpE2,KeqE2=paramsE2
    VM23,Ks23=params23
    vt,kt=params_trans
    al_p_23=allo[0]
    
    # return mu_rev(mE,m2,1.,1.,VME2,KsE2,KpE2,KeqE2)-\
    #        mu_rev(m2**al_p_23,m3,1.,1.,VM23,Ks23,Kp23,Keq23)-\
    #        transp(m2,vt,kt)
    return mu_rev(mE,m2,1.,1.,VME2,KsE2,KpE2,KeqE2)-\
           mu(m2**al_p_23, 1., VM23, Ks23)-\
           transp(m2,vt,kt)

@numba.njit()    
def dm3dt(P,params23,params35,params_diff_3,allo,t):
    m2=P[2]
    m3=P[3]
    m5=P[4]
    n3=P[6]
    # b0=P[9]
    # b1=P[10]
    al_p_23=allo[0]
    
    VM23,Ks23=params23
    VM35,Ks35,Kp35,Keq35=params35
    ds3,df3=params_diff_3
    
    return mu(m2**al_p_23, 1., VM23, Ks23) -\
           mu_rev(m3,m5,1.,1.,VM35,Ks35,Kp35,Keq35)-\
           diffn(m3,n3,ds3,df3)

@numba.njit()           
def dm5dt(P,params_inout,params35,params_diff_5,t):
    m3=P[3]
    m5=P[4]
    n5=P[7]
    # b0=P[9]
    # b1=P[10]
    
    VM35,Ks35,Kp35,Keq35=params35
    ds5,df5=params_diff_5
    kout5=params_inout[1]

    return mu_rev(m3,m5,1.,1.,VM35,Ks35,Kp35,Keq35) -\
           diffn(m5,n5,ds5,df5)- kout5*m5

@numba.njit()    
def dn2dt(P,params_24,params_trans,allo,t):
    m2=P[2]
    n2=P[5]
    n3=P[6]
    n4=P[8]
    # c0=P[11]
    # c1=P[12]
    
    vt,kt=params_trans
    VM_24,Ks_24=params_24
    al_p_24=allo[2]
    al_b_24=allo[3]
    
    return transp(m2,vt,kt) -\
        mu(n2**al_p_24,1.,VM_24,Ks_24)

@numba.njit()        
def dn3dt(P,params_35,params_diff_3,t):
    m3=P[3]
    n2=P[5]
    n3=P[6]
    n5=P[7]
    # c0=P[11]
    # c1=P[12]
    
    VM_35,Ks_35,Kp_35,Keq_35=params_35
    ds3,df3=params_diff_3
    
    return mu_rev(n3,n5,1.,1.,VM_35,Ks_35,Kp_35,Keq_35)+\
        diffn(m3,n3,ds3,df3)

@numba.njit()        
def dn5dt(P,params_35,params_diff_5,t):
    m5=P[4]
    n3=P[6]
    n5=P[7]
    # c0=P[11]
    # c1=P[12]

    VM_35,Ks_35,Kp_35,Keq_35=params_35
    ds5,df5=params_diff_5
    
    
    return mu_rev(n3,n5,1.,1.,VM_35,Ks_35,Kp_35,Keq_35) +  diffn(m5,n5,ds5,df5)

@numba.njit()
def dn4dt(P,params_24,params_inout,allo,t):
    n2=P[5]
    n4=P[8]
    n6=P[9]
    # c0=P[11]
    # c1=P[12]
    al_p_24=allo[2]
    al_b_24=allo[3]
    deg=params_inout[2]
    
    VM_24,Ks_24=params_24
    
    return mu(n2**al_p_24,1.,VM_24,Ks_24)-\
        deg*n4#mu_rev(n4,n6,1.,1.,VM_46,Ks_46,Kp_46,Keq_46)

# def dn6dt(P,params_46,params_inout):
#     n4=P[8]
#     n6=P[9]
    
#     VM_46,Ks_46,Kp_46,Keq_46=params_46
#     deg=params_inout[2]
    
#     return mu_rev(n4,n6,1.,1.,VM_46,Ks_46,Kp_46,Keq_46) - deg*n6\

@numba.njit()
def db0dt(P,params01,params35,paramsB,t):
    # m0=P[0]
    # m1=P[1]
    # m3=P[3]
    # m5=P[4]
    # # b0=P[9]
    # # b1=P[10]
    
    # VM01,Ks01,Kp01,Keq01=params01
    # VM35,Ks35,Kp35,Keq35=params35
    # VMB,KsB,KpB,KeqB=paramsB
    
    return 0

@numba.njit()        
def db1dt(P,params01,params35,paramsB,t):
    return 0

@numba.njit()
def dc0dt(P,params_35,params_24,paramsC,allo,t):
    # n2=P[5]
    # n3=P[6]
    # n5=P[7]
    # n4=P[8]
    # c0=P[11]
    # c1=P[12]
    
    
    # VM_24,Ks_24=params_24
    # VM_35,Ks_35,Kp_35,Keq_35=params_35
    # VMC,KsC,KpC,KeqC=paramsC
    
    # al_p_24=allo[2]
    # al_b_24=allo[3]
    
    return 0

@numba.njit()
def dc1dt(P,params_35,params_24,paramsC,allo,t):
    return 0
    
    
    
@numba.njit()    
def dPdt(P,t,params_inout, params01,params1E,paramsE2,params23,params35,
         params_24,params_35,
         params_trans,params_diff_3,params_diff_5,
         paramsB,paramsC,allo):
    
   # print(P,params_24,params_trans,allo,t)
    return [dm0dt(P,params_inout,params01,t),
            dm1dt(P,params01,params1E,t),
            dm2dt(P,paramsE2,params23,params_trans,allo,t),
            dm3dt(P,params23,params35,params_diff_3,allo,t),
            dm5dt(P,params_inout,params35,params_diff_5,t),
            dn2dt(P,params_24,params_trans,allo,t),
            dn3dt(P,params_35,params_diff_3,t),
            dn5dt(P,params_35,params_diff_5,t),
            dn4dt(P,params_24,params_inout,allo,t),
            db0dt(P,params01,params35,paramsB,t),
            db1dt(P,params01,params35,paramsB,t),
            dc0dt(P,params_35,params_24,paramsC,allo,t),
            dc1dt(P,params_35,params_24,paramsC,allo,t),
            dmEdt(P,params1E,paramsE2,t)]  



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    