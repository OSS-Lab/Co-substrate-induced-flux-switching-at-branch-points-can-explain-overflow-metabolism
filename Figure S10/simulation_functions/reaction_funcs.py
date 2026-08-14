#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:01:01 2024

@author: robert
"""

import numba
import numpy as np

@numba.njit
def mu(sub, prod, cosub, coprod, CKcat, Ks, Kp, Keq):
    
    s=sub*cosub/Ks
    
    return CKcat  * s/(1.  + s)

@numba.njit
def mu_rev(sub, prod, cosub, coprod, CKcat, Ks, Kp, Keq):
    
    s=sub*cosub/Ks
    p=prod*coprod/Kp
    
    rat=prod*coprod/(sub*cosub*Keq)
    return CKcat * (1. - rat) * s/(1. + p + s)

@numba.njit
def forwardflux_rev(sub, prod, cosub, coprod, CKcat, Ks, Kp, Keq):
    
    s=sub*cosub/Ks
    p=prod*coprod/Kp
    
    rat=prod*coprod/(sub*cosub*Keq)
    return CKcat * s/(1. + p + s)


#@numba.njit
def make_initial_concs(lengthparams,params_B):
   # print(np.sum(lengthparams))
    veclength=np.sum(lengthparams)+3 # there is one more metabolite than 
    # number of reactions and two co substrates
    Btot=float(params_B[-1])
    
    # make initial conc vector with concs of 0.01
    P0=0.01*np.ones(veclength)
    # make the first two represent B0 and B1
    P0[[-2,-1]]=[Btot/2,Btot/2]
    # first metabolite has higher concs
    P0[0]=0.1
    
    return P0


#@numba.njit
def make_kin_vec(lengthparams):
  #  kin=params_inout[0]
    veclength=sum(lengthparams)+3 # there is one more metabolite than 
    # number of reactions and two co substrates
    kinvec=np.zeros(veclength)
    kinvec[0]=1
    
    return kinvec
    
def make_kout_vec(lengthparams):
   # kout_top=params_inout[1]
   # kout_bottom=params_inout[2]
    veclength=sum(lengthparams)+3 # there is one more metabolite than 
    # number of reactions and two co substrates
    
    upstream_length=lengthparams[0]
    topbranch_length=lengthparams[1]
    bottombranch_length=lengthparams[2]
    
    kouttopvec=np.zeros(veclength)
    koutbottomvec=np.zeros(veclength)
    kouttopvec[upstream_length+topbranch_length]=1
    koutbottomvec[upstream_length+topbranch_length+bottombranch_length]=1
    return kouttopvec, koutbottomvec
    
def make_MM_mat_noB(lengthparams,B_position_params):
    veclength=sum(lengthparams)+3
    
    upstream_length=lengthparams[0]
    topbranch_length=lengthparams[1]
    bottombranch_length=lengthparams[2]
    
    upstream01=B_position_params[0]
    upstream10=B_position_params[1]
    topbranch10=B_position_params[2]
    bottombranch01=B_position_params[3]
    
    MM_mat_noB=np.zeros((veclength,veclength))
    
    
    
    # first do the upstream chain up to but not including the branch point
    for i in range(upstream_length):
        if upstream01[i]==0 and upstream10[i]==0:
            MM_mat_noB[i+1,i]=1
    
    # now do branch points
    if topbranch10[0]==0:
        MM_mat_noB[upstream_length+1,upstream_length]=1
    if bottombranch01[0]==0:
        MM_mat_noB[upstream_length+topbranch_length+1,upstream_length]=1
   
    # now do top branch after branch point
    for i in range(topbranch_length):
        if i>0:
            if topbranch10[i]==0:
                MM_mat_noB[upstream_length+i+1,upstream_length+i]=1
    
    # now do bottom branch after branch point
    for i in range(bottombranch_length):
        if i>0:
            if bottombranch01[i]==0:
                MM_mat_noB[upstream_length+topbranch_length+i+1,upstream_length+topbranch_length+i]=1
    return MM_mat_noB

def make_MM_mat_B01(lengthparams,B_position_params):
    veclength=sum(lengthparams)+3
    
    upstream_length=lengthparams[0]
    topbranch_length=lengthparams[1]
    bottombranch_length=lengthparams[2]
    
    upstream01=B_position_params[0]
    bottombranch01=B_position_params[3]
    
    MM_mat_B01=np.zeros((veclength,veclength))
    
    # first do the upstream chain up to but not including the branch point
    for i in range(upstream_length):
        if upstream01[i]==1:
            MM_mat_B01[i+1,i]=1
            
    # now do bottom branch point
    if bottombranch01[0]==1:
        MM_mat_B01[upstream_length+topbranch_length+1,upstream_length]=1
        
    # now do bottom branch after branch point
    for i in range(bottombranch_length):
        if i>0:
            if bottombranch01[i]==1:
                MM_mat_B01[upstream_length+topbranch_length+i+1,upstream_length+topbranch_length+i]=1
    
    return MM_mat_B01

def make_MM_mat_B10(lengthparams,B_position_params):
    veclength=sum(lengthparams)+3
    
    upstream_length=lengthparams[0]
    topbranch_length=lengthparams[1]
 #   bottombranch_length=lengthparams[2]
    
 #   upstream01=B_position_params[0]
    upstream10=B_position_params[1]
    topbranch10=B_position_params[2]
#    bottombranch01=B_position_params[3]
    
    MM_mat_B10=np.zeros((veclength,veclength))
    
    # first do the upstream chain up to but not including the branch point
   # print(upstream10)

    for i in range(upstream_length):
      #  print(i)
        if upstream10[i]==1:
          #  print('yes')
          #  print(i+1,i)
            MM_mat_B10[i+1,i]=1
    
    # now do top branch point
   # print(topbranch10[0])
    if topbranch10[0]==1:
        #print('yes')
       # print(upstream_length+1,upstream_length)
        MM_mat_B10[upstream_length+1,upstream_length]=1
    
   # print(topbranch10)
    # now do top branch after branch point
    for i in range(topbranch_length):
        if i>0:
            if topbranch10[i]==1:
               # print('yes')
               # print(upstream_length+i+1,upstream_length+i)
                MM_mat_B10[upstream_length+i+1,upstream_length+i]=1
    
    
    
    return MM_mat_B10


@numba.njit
def influx(kinvec,params):
    kin = params[0]
    kin_rate=kin*kinvec
    return kin_rate

@numba.njit
def kouttoprate(P,params,kouttopvec):
    kouttop = params[1]
    rate_kouttop=-kouttop*kouttopvec*P
#    print(rate_kouttop,P)
    return rate_kouttop

@numba.njit    
def koutbottomrate(P,params,koutbottomvec):
    koutbottom = params[2]
    rate_koutbottom=-koutbottom*koutbottomvec*P
#    print(rate_koutbottom)
    return rate_koutbottom

#@numba.njit
def noBrate(P,params,lengthparams,MM_mat_noB,reversible):
    # print(P)
    # print(params)
    # print(lengthparams)
    # print(MM_mat_noB)
    # make the vector, last 2 entries are for B0 and B1 and are therefore unused
    rate_noB=np.zeros(np.sum(lengthparams)+3)
    
    CKcat_gobal =     params[3]
    Ks_gobal =        params[4]
    Kp_gobal =        params[5]
    Keq_gobal =       params[6]
    
    CKcat_top =       params[7]
    Ks_top =          params[8]
    Kp_top =          params[9]
    Keq_top =         params[10]
    
    CKcat_bottom =    params[11]
    Ks_bottom =       params[12]
    Kp_bottom =       params[13]
    Keq_bottom =      params[14]
    
   # print(lengthparams)
    if reversible == True:
        reacf=mu_rev
    elif reversible ==False:
        reacf=mu
        
    # for each metabolite that is not a cosub
   # print(np.sum(lengthparams)+1)
    for i in range(np.sum(lengthparams)+1):
        # print(i)
        reac_out_indx=MM_mat_noB[:,i]

        # do the branch point separately
        if i != lengthparams[0]:
            # print('in this branch')
            # 
            for j in range(len(reac_out_indx)):
                if reac_out_indx[j]==1:
                    rate_noB[i] += -reacf(P[i], P[j], 1., 1., CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    rate_noB[j] += reacf(P[i], P[j], 1., 1., CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
        
        else:
           # print('in that branch')
           for j in range(len(reac_out_indx)):
               if reac_out_indx[j]==1:
                   if j==i+1:
                        ckcat=CKcat_top
                        ks=Ks_top
                        kp=Kp_top
                        keq=Keq_top
                   else:
                        ckcat=CKcat_bottom
                        ks=Ks_bottom
                        kp=Kp_bottom
                        keq=Keq_bottom
                   
                   rate_noB[i] += -reacf(P[i], P[j], 1., 1., ckcat, ks, kp, keq)
                   rate_noB[j] += reacf(P[i], P[j], 1., 1., ckcat, ks, kp, keq)
        
                # if reac_in_indx[j]==1:
                #    rate += mu(P[j], P[i], 1., 1., CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    
    return rate_noB

def B01rate(P,params,lengthparams,MM_mat_B01,reversible):
    
    
    rate_B01=np.zeros(np.sum(lengthparams)+3)
    
    CKcat_gobal =     params[3]
    Ks_gobal =        params[4]
    Kp_gobal =        params[5]
    Keq_gobal =       params[6]
    
    CKcat_top =       params[7]
    Ks_top =          params[8]
    Kp_top =          params[9]
    Keq_top =         params[10]
    
    CKcat_bottom =    params[11]
    Ks_bottom =       params[12]
    Kp_bottom =       params[13]
    Keq_bottom =      params[14]
    
    if reversible == True:
        reacf=mu_rev
    elif reversible ==False:
        reacf=mu
    
    # for each metabolite that is not a cosub
    for i in range(np.sum(lengthparams)+1):
        reac_out_indx=MM_mat_B01[:,i]
        # do the branch point separately
        if i != lengthparams[0]:
            
            
            # 
            for j in range(len(reac_out_indx)):
                if reac_out_indx[j]==1:
                    rate_B01[i] += -reacf(P[i], P[j], P[-2], P[-1], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    rate_B01[j] += reacf(P[i], P[j], P[-2], P[-1], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    
                    rate_B01[-2] += -reacf(P[i], P[j], P[-2], P[-1], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    rate_B01[-1] += reacf(P[i], P[j], P[-2], P[-1], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    
        else:
        #   print('on branch point' +str(i))

           for j in range(len(reac_out_indx)):
               if reac_out_indx[j]==1:
                #   print(str(j))
                   if j==i+1:
                    #    print('top')
                        ckcat=CKcat_top
                        ks=Ks_top
                        kp=Kp_top
                        keq=Keq_top
                   else:
                      #  print('bottom')
                        ckcat=CKcat_bottom
                        ks=Ks_bottom
                        kp=Kp_bottom
                        keq=Keq_bottom
                   
                   
                   rate_B01[i] += -reacf(P[i], P[j], P[-2], P[-1], ckcat, ks, kp, keq)
                   rate_B01[j] += reacf(P[i], P[j], P[-2], P[-1], ckcat, ks, kp, keq)
                   
                   rate_B01[-2] += -reacf(P[i], P[j], P[-2], P[-1], ckcat, ks, kp, keq)
                   rate_B01[-1] += reacf(P[i], P[j], P[-2], P[-1], ckcat, ks, kp, keq)
    
    return rate_B01


def B10rate(P,params,lengthparams,MM_mat_B10,reversible):
    
    rate_B10=np.zeros(np.sum(lengthparams)+3)
    
    CKcat_gobal =     params[3]
    Ks_gobal =        params[4]
    Kp_gobal =        params[5]
    Keq_gobal =       params[6]
    
    CKcat_top =       params[7]
    Ks_top =          params[8]
    Kp_top =          params[9]
    Keq_top =         params[10]
    
    CKcat_bottom =    params[11]
    Ks_bottom =       params[12]
    Kp_bottom =       params[13]
    Keq_bottom =      params[14]
    
    if reversible == True:
        reacf=mu_rev
    elif reversible ==False:
        reacf=mu
    
    # for each metabolite that is not a cosub
    for i in range(np.sum(lengthparams)+1):
        reac_out_indx=MM_mat_B10[:,i]
        # do the branch point separately
        if i != lengthparams[0]:
            
            
            # 
            for j in range(len(reac_out_indx)):
                if reac_out_indx[j]==1:
                    rate_B10[i] += -reacf(P[i], P[j], P[-1], P[-2], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    rate_B10[j] += reacf(P[i], P[j], P[-1], P[-2], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    
                    rate_B10[-2] += reacf(P[i], P[j], P[-1], P[-2], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    rate_B10[-1] += -reacf(P[i], P[j], P[-1], P[-2], CKcat_gobal, Ks_gobal, Kp_gobal, Keq_gobal)
                    
        else:
          # print('on branch point' +str(i))
           for j in range(len(reac_out_indx)):
               
               if reac_out_indx[j]==1:
                #   print(str(j))
                   if j==i+1:
                       # print('top')
                        ckcat=CKcat_top
                        ks=Ks_top
                        kp=Kp_top
                        keq=Keq_top
                   else:
                     #   print('bottom')
                        ckcat=CKcat_bottom
                        ks=Ks_bottom
                        kp=Kp_bottom
                        keq=Keq_bottom
                   
                   rate_B10[i] += -reacf(P[i], P[j], P[-1], P[-2], ckcat, ks, kp, keq)
                   rate_B10[j] += reacf(P[i], P[j], P[-1], P[-2], ckcat, ks, kp, keq)
                   
                   rate_B10[-2] += reacf(P[i], P[j], P[-1], P[-2], ckcat, ks, kp, keq)
                   rate_B10[-1] += -reacf(P[i], P[j], P[-1], P[-2], ckcat, ks, kp, keq)
    return rate_B10

def Bback(P,params,lengthparams):
    # add background B0 -> B1 reaction
    rate_Bback=np.zeros(np.sum(lengthparams)+3)
    
    CKcatB =          params[15]
    KsB =             params[16]
    KpB =             params[17]
    KeqB =            params[18]
    
    rate_Bback[-2]=-mu_rev(P[-2], P[-1], 1., 1., CKcatB, KsB, KpB, KeqB)
    rate_Bback[-1]=mu_rev(P[-2], P[-1], 1., 1., CKcatB, KsB, KpB, KeqB)
    
    return rate_Bback

def Binout(P,params,lengthparams):
    
    rate_Binout=np.zeros(np.sum(lengthparams)+3)
    b0in = params[20]
    b0out = params[21]
    b1out = params[22]
    
    
    # rate_Binout[-2]=b_in*(P[-2]/(P[-2]+P[-1])) - bdeg*P[-2]
    # rate_Binout[-1]=b_in*(P[-1]/(P[-2]+P[-1]))-bdeg*P[-1]
    
    rate_Binout[-2]=b0in - b0out*P[-2]
    rate_Binout[-1]= - b1out*P[-1]
    
    return rate_Binout

# lengthparams=[2,1,1]
# B_position_params=([1,0],[0,0],[1],[1])
# params_global=[1.,1.,1.,1.]
# params_B=[0.1,1.,1.,1.,0.1]
# params_branch=[[1.,1.,1.,1.],[1.,1.,1.,1.]]
# params_inout=[0.01,1.,1.]
# reversible=False
# MM_mat_noB=make_MM_mat_noB(lengthparams,B_position_params)
# print(MM_mat_noB)
# print(MM_mat_noB[:,1])
# print(MM_mat_noB[2,:])
# noBrate(0,0,np.array([[0,1,0],[0,0,0],[0,0,1]]))

# make_initial_concs(lengthparams,params_B)