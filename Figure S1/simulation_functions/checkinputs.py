#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 13:51:47 2024

@author: robert
"""
import numpy as np
import numba

def check_length_params(lengthparams):
    upstream_length=lengthparams[0]
    topbranch_length=lengthparams[1]
    bottombranch_length=lengthparams[2]
    
    failed = 0
    
    if upstream_length < 2:
        print('must have at least 2 metabolite reactions upstream')
        failed = 1
        
    if topbranch_length < 1:
        print('must have at least 2 metabolite reactions on top branch')
        failed = 1

    if bottombranch_length < 1:
        print('must have at least 2 metabolite reactions on bottom branch')
        failed = 1

    return failed

def check_B_position_params(lengthparams,B_position_params):
    upstream_length=lengthparams[0]
    topbranch_length=lengthparams[1]
    bottombranch_length=lengthparams[2]
    
    upstream_b0_positions=B_position_params[0]
    upstream_b1_positions=B_position_params[1]                  
    topbranch_b1_position=B_position_params[2]
    bottombranch_b0_position=B_position_params[3]
    
    failed  = 0
    
    if len(upstream_b0_positions)!=upstream_length:
        print('upstream B0 positions do not match number of reactions')
        failed = 1
    if len(upstream_b1_positions)!=upstream_length:
        print('upstream B1 positions do not match number of reactions')
        failed = 1
    if len(topbranch_b1_position)!=topbranch_length:
        print('topbranch B1 positions do not match number of reactions')
        failed = 1
    if len(bottombranch_b0_position)!=bottombranch_length:
        print('bottombranch B0 positions do not match number of reactions')
        failed = 1
        
    if failed == 1:
        return failed
    else:
        checkupstream=np.array(upstream_b0_positions)+np.array(upstream_b1_positions)
        if max(checkupstream)>1:
            print('cannot have B0 -> B1 and B1 -> B0 on the same reaction \n'+ 
                  'check upstream B positions'+ str(checkupstream))
            failed = 1
     #   if max(checkupstream)<1:
      #      print('need at least one B reaction upstream')
       #     failed = 1
#        if max(np.array(topbranch_b1_position)) < 1:
 #           print('need at least one B1 -> B0 reaction on top branch')
  #          failed = 1
   #     if max(np.array(bottombranch_b0_position)) < 1:
    #        print('need at least one B0 -> B1 reaction on bottom branch')
            
    return failed
        
def check_params(params_global,params_B,params_branch):
    failed = 0
    
    if len(params_global)!=4:
        print('params_global needs 4 parameters')
        failed = 1
        
    if len(params_B)!=5:
        print('params_B needs 5 parameters')
        failed = 1
        
    if len(params_branch)!=8:
        print('params_branch needs 8 parameters')
        failed = 1
        
    return failed
        
# def check_reaction_type(reversible):
#     failed = 0
#     if reversible != False:
#         failed = 1
#     return failed
    