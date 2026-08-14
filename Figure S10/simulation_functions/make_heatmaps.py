#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 12:57:55 2024

@author: robert
"""

import numpy as np

import time
from datetime import timedelta
import multiprocessing as multi
import make_odes as mp
import sys

#from pathlib import Path
import os
import time



sys.path.append('../data/megadata')
sys.path.append('../data/heatmapdata')


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
    return str(np.round(np.log10(x),3))

def make_data_parallel(fn,inp_args):#,totleng):
    #print('this is doing now')
    """
    runs a function in parrallel with arguments from a list using 2 less
    threads the total
    # """
    with multi.Pool(os.cpu_count()-4) as pool:
        pool.starmap(fn,inp_args)
    # with multi.Pool(1) as pool:
    #     pool.starmap(fn,inp_args)

# def findvarindict(lab,d):
#     if lab in d:
#         return True
#     else:
#         return False

def findvarindict(key, d):
    return key in d

def print_vars_from_dict(d):
    for key, value in d.items():
        if findvarindict(key, d):
            globals()[key] = value


def get_fluxfractions(lengthparams,B_position_params,
            params_global,params_B,
            params_branch,
            params_inout,
            bsynthdeg,
            reversible,Tmax,tpts):
    """
    Return the flux fraction, b fraction and whether there is buildup for a given
    set of parameters
    """
    result,time_points,odefunc,params = mp.make_time_series(lengthparams,B_position_params,
            params_global,params_B,
            params_branch,
            params_inout,
            bsynthdeg,
            reversible,Tmax,tpts)
    #print(result)
  #  print('this is the final result',result[-1,:])
  #  print(np.any(result[-1,:]))
    if np.any(result[-1,:]):
    # if result[-1,:] == not np.any(a):
    #     print('all zeros')
        isbuildup=mp.check_buildup(result,odefunc,params,reversible)

    #    print('this checks for buildup', isbuildup)
        flux_top,flux_bottom,bfrac,brat=mp.get_branch_fluxes(result,odefunc,params,lengthparams,B_position_params,reversible)
        fluxfrac=flux_bottom/(flux_top+flux_bottom)
    else:
        isbuildup=np.array([np.nan])
        fluxfrac=np.array([np.nan])
        bfrac=np.array([np.nan])
        flux_top=np.array([np.nan])
        flux_bottom=np.array([np.nan])
    return isbuildup, fluxfrac, bfrac, flux_top, flux_bottom






def get_ff_delta(lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,10.],
            params_branch=[1.,1.,1.,1.,
                           1.,1.,1.,1.],
            kouts=[1.,1.],#params_inout=[0.1,1.,1.],
            kin_min=-8,kin_max=0,num_kin=41,
            binout=[1.,1.,1.],
            reversible=False,Tmax=1000000,tpts=100000,savepath='../data/test.npy'):
    """


    Parameters
    ----------
    lengthparams :
        gives the number of upstream, topbranch and bottom branch reactions.
        The default is [2,1,1].

    B_position_params :
        reactions where there is a b0->b1 reaction upstream,
        b1->b0 reaction upstream,
        b1->b0 on the top branch,
        b0->b1 on the bottom branch.
        The default is ([1,0],[0,0],[1],[1]).

    params_global :
        CKcat,Ks,Kp,Keq for all reactions other than the branch reactions and
        the off pathway b reaction
        . The default is [1.,1.,1.,1.].

    params_B :
        CKcatB,KsB,KpB,KeqB,Btot.
        The default is [0.1,1.,1.,1.,10.].

    params_branch :
        [CKcat_top,   Ks_top,   Kp_top,   Keq_top,
        CKcat_bottom,Ks_bottom,Kp_bottom,Keq_bottom]
        The default is [1.,1.,1.,1.,
                        1.,1.,1.,1.].

    kouts :
        kout_top,kout_bottom.
        The default is [1.,1.].

    kin_min :
        lowest power of 10 for kin.
        The default is -8.

    kin_max :
        highest power of 10 for kin.
        The default is 0.

    num_kin :
        number of values for kin.
        The default is 41.

    reversible :
        use reversible or irreversible reactions. Affects all reactions apart
        from the in/out fluxes (always irreversible) and the off pathway B
        reactions (always reversible).
        The default is False.

    Tmax :
        maximum time for the ode simulations.
        The default is 1000000.

    tpts :
        number of time points for the ode simulations, should be Tmax/10.
        The default is 100000.

    savepath :
        path to save the data. The default is '../data/test.npy'.


    Returns
    -------
    ff_delta :
        signed change in flux fraction in to bottom branch as kin varies.

    """
    kinvals=np.logspace(kin_min,kin_max,num_kin)
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    b_vec=np.zeros(np.size(kinvals))
    b_vec[:]=np.nan
    f_top=np.zeros(np.size(kinvals))
    f_top[:]=np.nan
    f_bottom=np.zeros(np.size(kinvals))
    f_bottom[:]=np.nan
  #  print(binout)
    for id_kin, kin in enumerate(kinvals):
        params_inout=np.concatenate(([kin],kouts))
   #     print(np.log10(kin))


        isbuildup, fluxfrac, bfrac, flux_top, flux_bottom =get_fluxfractions(lengthparams,B_position_params,
            params_global,params_B,params_branch,
            params_inout,binout,
            reversible,Tmax,tpts)
   #     print(isbuildup)
        if sum(isbuildup)==0:
         #   print('here')
            ff_vec[id_kin]=fluxfrac
            b_vec[id_kin]=bfrac
            f_top[id_kin]=flux_top
            f_bottom[id_kin]=flux_bottom

    ff_delta_unsigned=np.nanmax(ff_vec)-np.nanmin(ff_vec)
    arg_min_flux=_nanargmin(ff_vec)
    arg_max_flux=_nanargmax(ff_vec)
    ffsign=np.sign(arg_max_flux- arg_min_flux)
    ff_delta=ffsign*ff_delta_unsigned

    bf_delta_unsigned=np.nanmax(b_vec)-np.nanmin(b_vec)
    barg_min_flux=_nanargmin(b_vec)
    barg_max_flux=_nanargmax(b_vec)
    bfsign=np.sign(barg_max_flux- barg_min_flux)
    bf_delta=bfsign*bf_delta_unsigned

    ff_grad=np.gradient(ff_vec,kinvals)
    if np.any(~np.isnan(ff_grad)) == True:
        gradmax = ff_grad[_nanargmax(np.abs(ff_grad))]
    else:
        gradmax=np.nan


    ff_sens=ff_grad/(ff_vec/kinvals)
    if np.any(~np.isnan(ff_sens))==True:
        sensmax= ff_sens[_nanargmax(np.abs(ff_sens))]
    else:
        sensmax=np.nan

    print([ff_delta,bf_delta,gradmax,sensmax])
    np.savez(savepath,deltadat=[ff_delta,bf_delta,gradmax,sensmax],fluxtop=f_top,fluxbottom=f_bottom)
    return ff_delta


def Gen_all_data(lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
                 CKcat=1.,Ks=1.,Kp=1.,Keq=1.,               # params_global
                 KsB=1.,KpB=1.,                             # params_b unchanging
                 Ks_top=1., Kp_top=1.,
                 Ks_bottom=1., Kp_bottom=1.,                # params_branch unchanging
                 kout_top=1.,kout_bottom=1.,                # params_inout unchanging
                 CKcatB_vals=np.logspace(-2,0,2),
                 Keqb_vals=np.logspace(-2,0,2),
                 Btot_vals=np.logspace(-2,0,2),             # params_b variable
                 CKcat_top_vals=np.logspace(-2,0,2),
                 Keq_top_vals=np.logspace(-2,0,2),
                 CKcat_bottom_vals=np.logspace(-2,0,2),
                 Keq_bottom_vals=np.logspace(-2,0,2),       # params_branch variable
                 kin_min=-8,kin_max=0,num_kin=41,           # params_inout variable
                 vals_b_in=np.logspace(-5,2,2),
                 vals_b_out=np.logspace(-5,2,2),
                 reversible=False,Tmax=1000000,tpts=100000):
    """


    Parameters
    ----------
    lengthparams :
        gives the number of upstream, topbranch and bottom branch reactions.
        The default is [2,1,1].

    B_position_params :
        reactions where there is a b0->b1 reaction upstream,
        b1->b0 reaction upstream,
        b1->b0 on the top branch,
        b0->b1 on the bottom branch.
        The default is ([1,0],[0,0],[1],[1]).

    CKcat :
        parameter for the non-focused reactions.
        The default is 1..

    Ks :
        parameter for the non-focused reactions.
        . The default is 1..

    Kp :
        parameter for the non-focused reactions.
        . The default is 1..

    Keq :
        parameter for the non-focused reactions.
        . The default is 1..

    KsB :
        parameter for the off-pathway B reaction.
        The default is 1..

    KpB :
        parameter for the off-pathway B reaction.
        The default is 1..

    Ks_top :
        parameter for the reaction into the top branch.
        The default is 1..

    Kp_top :
        parameter for the reaction into the top branch.
        The default is 1..

    Ks_bottom :
        parameter for the reaction into the bottom branch. The default is 1..

    Kp_bottom :
        parameter for the reaction into the bottom branch. The default is 1..

    kout_top:
        parameter for flux out of top branch.
        The default is 1..

    kout_bottom :
        parameter for flux out of bottom branch.
        The default is 1..

    CKcatB_vals :
        list of CKcatB vals.
        The default is np.logspace(-2,0,2).

    Keqb_vals :
        list of KeqB vals.
        The default is np.logspace(-2,0,2).

    Btot_vals :
        list of Btot vals.
        The default is np.logspace(-2,0,2).

    CKcat_top_vals :
        list of CKcat_top vals.
        The default is np.logspace(-2,0,2).

    Keq_top_vals :
        list of Keq_top vals.
        The default is np.logspace(-2,0,2).

    CKcat_bottom_vals :
        list of CKcat_bottom vals.
        The default is np.logspace(-2,0,2).

    Keq_bottom_vals :
        list of Keq_bottom vals.
        The default is np.logspace(-2,0,2).

    kin_min :
        lowest power of 10 for kin.
        The default is -8.

    kin_max :
        highest power of 10 for kin.
        The default is 0.

    num_kin :
        number of values for kin.
        The default is 41.

    reversible :
        use reversible or irreversible reactions. Affects all reactions apart
        from the in/out fluxes (always irreversible) and the off pathway B
        reactions (always reversible).
        The default is False.

    Tmax :
        maximum time for the ode simulations.
        The default is 1000000.

    tpts :
        number of time points for the ode simulations, should be Tmax/10.
        The default is 100000.

    Returns
    -------
    list of sets of parameters for which simulations are needed.

    """

    #num_kin=int(num_kin[0])
    path='../data/megadata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
        '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
        str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
        '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
        '_'+str([kout_top,kout_bottom])

    # print(path)



    if not os.path.exists(path):
        print('Directory does not exist')
        os.makedirs(path)
        print("Directory created successfully")
    else:
        print("Directory already exists")


    list_of_args=list()

    params_global= [CKcat,Ks,Kp,Keq]
    kouts=[kout_top, kout_bottom]
    kinvals=np.logspace(kin_min,kin_max)
    print(len(CKcatB_vals))
    print(len(Keqb_vals))
    print(len(Btot_vals))
    print(len(CKcat_top_vals))
    print(len(Keq_top_vals))
    print(len(CKcat_bottom_vals))
    print(len(Keq_bottom_vals))
    print(len(vals_b_in))
    print(len(vals_b_out))

    for CKcatB in CKcatB_vals:
        for KeqB in Keqb_vals:
            for Btot in Btot_vals:
                for CKcat_top in CKcat_top_vals:
                    for Keq_top in Keq_top_vals:
                        for CKcat_bottom in CKcat_bottom_vals:
                            for Keq_bottom in Keq_bottom_vals:
                                for b_in in vals_b_in:
                                    for b_out in vals_b_out:


                                        filename='log_variable_params_'+\
                                            rd_str([CKcatB,KeqB,Btot,CKcat_top,Keq_top,CKcat_bottom,Keq_bottom])+\
                                            '__kinvals_'+str([kin_min,kin_max,num_kin])+\
                                            '__binout_'+rd_str([b_in,b_out,b_out])+\
                                            '__timevals_'+str([Tmax,tpts])+'.npy.npz'

                                      #  print(filename)
                                        params_B=[CKcatB,KsB,KpB,KeqB,Btot]

                                        params_branch=[CKcat_top,   Ks_top,   Kp_top,   Keq_top,
                                                       CKcat_bottom,Ks_bottom,Kp_bottom,Keq_bottom]

                                        params_binout=[b_in,b_out,b_out]
                                    #    print(num_kin)
                                        savepath=path+'/'+filename
        #                                print(path+filename)
                                        # print('\n\n')
                                        # print(path)
                                        # print('\n')
                                        # print(filename)
                                        # print('\n\n')
                                        if not os.path.isfile(savepath):
                                          #  print(savepath+'\n')
                                         #   print('file does not exist, adding to list')
                                            list_of_args.append(tuple((lengthparams,
                                                B_position_params,params_global,params_B,
                                                params_branch,kouts,kin_min,kin_max,num_kin,
                                                params_binout,
                                                reversible,Tmax,tpts,savepath)))
                                        #else:
                                           # print('file exists')
                                        # if file doesnt exits add it to list
            # print(list_of_args)

    howmany=len(list_of_args)
    if howmany==0:
        print('all data points exist')
        return list()
    else:
        print(f'there are {howmany} data points to add to the list')
        return list_of_args




def make_heatmaps(xvar='b_in',xvarvals=np.logspace(-5,2,8),
                  yvar='b_out',yvarvals=np.logspace(-5,2,8),
                  lengthparams=[2,2,2],B_position_params=([1,0],[0,0],[0,1],[0,1]),
                  CKcat=1.,Ks=1.,Kp=1.,Keq=1.,               # params_global
                  KsB=1.,KpB=1.,                             # params_b unchanging
                  Ks_top=1., Kp_top=1.,
                  Ks_bottom=1., Kp_bottom=1.,                # params_branch unchanging
                  kout_top=1.,kout_bottom=1.,
                  variable_params=['Btot',
                                   'CKcatB'
                                   'KeqB',
                                   'CKcat_top',
                                   'Keq_top',
                                   'CKcat_bottom',
                                   'Keq_bottom'],
                  variable_param_values=(np.logspace(-2,0,3),
                                         np.logspace(-2,0,3),
                                         np.logspace(-2,0,3),
                                         [1.],
                                         [1.],
                                         np.logspace(-1,1,3),
                                         [1.]),
                  kin_min=-8,kin_max=0,num_kin=41,
                  reversible=False,Tmax=1000000,tpts=100000):


    path='../data/megadata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
        '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
        str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
        '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
        '_'+str([kout_top,kout_bottom])



    xmin=np.log10(np.min(xvarvals))
    xmax=np.log10(np.max(xvarvals))
    xnum=len(xvarvals)

    ymin=np.log10(np.min(yvarvals))
    ymax=np.log10(np.max(yvarvals))
    ynum=len(yvarvals)

    # print(xmin,xmax,xnum)
    # print(ymin,ymax,ynum)


    all_variable_values=['CKcatB',
                         'KeqB',
                         'Btot',
                         'CKcat_top',
                         'Keq_top',
                         'CKcat_bottom',
                         'Keq_bottom']

    # Creating a dictionary to store values
    variable_dict = {}

    xvar_dict = {}

    yvar_dict = {}



    for idvar0, var_vals0 in enumerate(variable_param_values[0]):
        variable_name0 = variable_params[0]
        variable_value0 = var_vals0
        variable_dict[variable_name0] = variable_value0

        # # Accessing the value using the converted string
        # retrieved_value = variable_dict[variable_name]
        # print(f"{variable_name}: {retrieved_value}")

        for idvar1, var_vals1 in enumerate(variable_param_values[1]):

            # Converting a string into a variable name and assigning a value
            variable_name1 = variable_params[1]
            variable_value1 = var_vals1
            variable_dict[variable_name1] = variable_value1


            for idvar2, var_vals2 in enumerate(variable_param_values[2]):
                variable_name2 = variable_params[2]
                variable_value2 = var_vals2
                variable_dict[variable_name2] = variable_value2


                for idvar3, var_vals3 in enumerate(variable_param_values[3]):
                    variable_name3 = variable_params[3]
                    variable_value3 = var_vals3
                    variable_dict[variable_name3] = variable_value3

                    for idvar4, var_vals4 in enumerate(variable_param_values[4]):
                        variable_name4 = variable_params[4]
                        variable_value4 = var_vals4
                        variable_dict[variable_name4] = variable_value4

                        for idvar5, var_vals5 in enumerate(variable_param_values[5]):
                            variable_name5 = variable_params[5]
                            variable_value5 = var_vals5
                            variable_dict[variable_name5] = variable_value5

                            for idvar6, var_vals6 in enumerate(variable_param_values[6]):
                                variable_name6 = variable_params[6]
                                variable_value6 = var_vals6
                                variable_dict[variable_name6] = variable_value6


                                print_vars_from_dict(variable_dict)

                                savestring_ff='ff_'
                                savestring_bf='bf_'
                                savestring_fg='fg_'
                                savestring_fs='fs_'

                                for key in sorted(variable_dict):
                                    savestring_ff += str(key) + '_' + str(variable_dict[key]) + '__'
                                    savestring_bf += str(key) + '_' + str(variable_dict[key]) + '__'
                                    savestring_fg += str(key) + '_' + str(variable_dict[key]) + '__'
                                    savestring_fs += str(key) + '_' + str(variable_dict[key]) + '__'

                            #print(savestring)

                                savepath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
                                    '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
                                    str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
                                    '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
                                    '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
                                    '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
                                    yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)

                               # print(savepath)
                                if not os.path.exists(savepath):
                                    print('Directory does not exist')
                                    os.makedirs(savepath)
                                    print("Directory created successfully")
                                else:
                                    print("Directory already exists")
                             #   heatmap_filename='heatmap__x
                                savename_ff=savepath+'/'+savestring_ff+'.npy'
                                print(savename_ff)
                                savename_bf=savepath+'/'+savestring_bf+'.npy'
                                print(savename_bf)
                                savename_fg=savepath+'/'+savestring_fg+'.npy'
                                print(savename_fg)
                                savename_fs=savepath+'/'+savestring_fs+'.npy'
                                print(savename_fs)


                                ffdelta_mat=np.zeros((ynum,xnum))
                                bfdelta_mat=np.zeros((ynum,xnum))
                                ffgrad_mat=np.zeros((ynum,xnum))
                                ffsense_mat=np.zeros((ynum,xnum))

                                for idx, varx in enumerate(xvarvals):
                                    variable_namex = xvar
                                    variable_valuex = varx
                                    xvar_dict[xvar] = variable_valuex

                                    print_vars_from_dict(xvar_dict)

                                    for idy, vary in enumerate(yvarvals):
                                        variable_namey = yvar
                                        variable_valuey = vary
                                        yvar_dict[yvar] = variable_valuey

                                        print_vars_from_dict(yvar_dict)

                                       # print(CKcatB,KeqB,Btot,CKcat_top,Keq_top,
                                     #         CKcat_bottom,Keq_bottom,'\n')

                                        filename='log_variable_params_'+\
                                            rd_str([CKcatB,KeqB,Btot,CKcat_top,Keq_top,CKcat_bottom,Keq_bottom])+\
                                            '__kinvals_'+str([kin_min,kin_max,num_kin])+\
                                            '__binout_'+rd_str([b_in, b_out, b_out])+\
                                            '__timevals_'+str([Tmax,tpts])+'.npy.npz'


                                        loadpath=path+'/'+filename
                                     #   print(loadpath,'\n')

                                        dat=np.load(loadpath)
                                        ff_dat=dat['deltadat']
                                        ffdelta=ff_dat[0]
                                        bfdelta=ff_dat[1]
                                        ffgrad=ff_dat[2]
                                        ffsense=ff_dat[3]

                                        ffdelta_mat[idy,idx]=ffdelta
                                        bfdelta_mat[idy,idx]=bfdelta
                                        ffgrad_mat[idy,idx]=ffgrad
                                        ffsense_mat[idy,idy]=ffsense
                                    #    print(ffdelta,'\n')


                         #       print(ffdelta_mat)
                                np.save(savename_ff,ffdelta_mat)
                                np.save(savename_bf,bfdelta_mat)
                                np.save(savename_fg,ffgrad_mat)
                                np.save(savename_fs,ffsense_mat)

# def make_heatmaps_ratio(xvar='Btot',xvarvals=np.logspace(-3,3,25),
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


#     path='../data/megadata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
#         '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
#         str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
#         '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
#         '_'+str([kout_top,kout_bottom])



#     xmin=np.log10(np.min(xvarvals))
#     xmax=np.log10(np.max(xvarvals))
#     xnum=len(xvarvals)

#     ymin=np.log10(np.min(yvarvals))
#     ymax=np.log10(np.max(yvarvals))
#     ynum=len(yvarvals)

#     # print(xmin,xmax,xnum)
#     # print(ymin,ymax,ynum)


#     all_variable_values=['CKcatB',
#                          'KeqB',
#                          'Btot',
#                          'CKcat_top',
#                          'Keq_top',
#                          'CKcat_bottom',
#                          'Keq_bottom']

#     # Creating a dictionary to store values
#     variable_dict = {}

#     xvar_dict = {}

#     yvar_dict = {}



#     for idvar0, var_vals0 in enumerate(variable_param_values[0]):
#         variable_name0 = variable_params[0]
#         variable_value0 = var_vals0
#         variable_dict[variable_name0] = variable_value0

#         # # Accessing the value using the converted string
#         # retrieved_value = variable_dict[variable_name]
#         # print(f"{variable_name}: {retrieved_value}")

#         for idvar1, var_vals1 in enumerate(variable_param_values[1]):

#             # Converting a string into a variable name and assigning a value
#             variable_name1 = variable_params[1]
#             variable_value1 = var_vals1
#             variable_dict[variable_name1] = variable_value1


#             for idvar2, var_vals2 in enumerate(variable_param_values[2]):
#                 variable_name2 = variable_params[2]
#                 variable_value2 = var_vals2
#                 variable_dict[variable_name2] = variable_value2


#                 for idvar3, var_vals3 in enumerate(variable_param_values[3]):
#                     variable_name3 = variable_params[3]
#                     variable_value3 = var_vals3
#                     variable_dict[variable_name3] = variable_value3

#                     for idvar4, var_vals4 in enumerate(variable_param_values[4]):
#                         variable_name4 = variable_params[4]
#                         variable_value4 = var_vals4
#                         variable_dict[variable_name4] = variable_value4


#                         print_vars_from_dict(variable_dict)

#                         savestring_fr='fr_'
#                         savestring_br='br_'
#                         savestring_frg='frg_'
#                         savestring_frs='frs_'

#                         for key in sorted(variable_dict):
#                             savestring_fr += str(key) + '_' + str(variable_dict[key]) + '__'
#                             savestring_br += str(key) + '_' + str(variable_dict[key]) + '__'
#                             savestring_frg += str(key) + '_' + str(variable_dict[key]) + '__'
#                             savestring_frs += str(key) + '_' + str(variable_dict[key]) + '__'

#                     #print(savestring)

#                         savepath='../data/heatmapdata_ratio/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
#                             '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
#                             str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
#                             '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
#                             '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
#                             '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
#                             yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)

#                        # print(savepath)
#                         if not os.path.exists(savepath):
#                             print('Directory does not exist')
#                             os.makedirs(savepath)
#                             print("Directory created successfully")
#                         else:
#                             print("Directory already exists")
#                      #   heatmap_filename='heatmap__x
#                         savename_fr=savepath+'/'+savestring_fr+'.npy'
#                         print(savename_fr)
#                         savename_br=savepath+'/'+savestring_br+'.npy'
#                         print(savename_br)
#                         savename_frg=savepath+'/'+savestring_frg+'.npy'
#                         print(savename_frg)
#                         savename_frs=savepath+'/'+savestring_frs+'.npy'
#                         print(savename_frs)


#                         frdelta_mat=np.zeros((ynum,xnum))
#                         frdelta_scaledmin_mat=np.zeros((ynum,xnum))
#                         frdelta_scaledmax_mat=np.zeros((ynum,xnum))

#                         frdelta_useful=np.zeros((ynum,xnum))

#                         brdelta_mat=np.zeros((ynum,xnum))
#                         frgrad_mat=np.zeros((ynum,xnum))
#                         frsense_mat=np.zeros((ynum,xnum))

#             #            ffdelta_mat=np.zeros((ynum,xnum))


#                         for idx, varx in enumerate(xvarvals):
#                             variable_namex = xvar
#                             variable_valuex = varx
#                             xvar_dict[xvar] = variable_valuex

#                             print_vars_from_dict(xvar_dict)

#                             for idy, vary in enumerate(yvarvals):
#                                 variable_namey = yvar
#                                 variable_valuey = vary
#                                 yvar_dict[yvar] = variable_valuey

#                                 print_vars_from_dict(yvar_dict)

#                                # print(CKcatB,KeqB,Btot,CKcat_top,Keq_top,
#                              #         CKcat_bottom,Keq_bottom,'\n')

#                                 filename='log_variable_params_'+\
#                                     rd_str([CKcatB,KeqB,Btot,CKcat_top,Keq_top,CKcat_bottom,Keq_bottom])+\
#                                     '__kinvals_'+str([kin_min,kin_max,num_kin])+\
#                                     '__timevals_'+str([Tmax,tpts])+'.npy.npz'


#                                 loadpath=path+'/'+filename
#                              #   print(loadpath,'\n')

#                                 dat=np.load(loadpath)
#                                 ftop=dat['fluxtop']
#                                 fbot=dat['fluxbottom']

#                               #  ff_dat=dat['deltadat']
#                              #   ffdelta=ff_dat[0]

#                                 frat=fbot/ftop

#                                 # print('\n#######')
#                                 # print(frat)
#                                 # print(np.nanmin(frat))
#                                 # print(np.nanmax(frat))
#                                 # print('#######\n')
#                                 frat_delta_unsigned=np.nanmax(frat)-np.nanmin(frat)

#                                 fratlog_delta_unsigned=np.abs(np.log10(np.nanmax(frat))-np.log10(np.nanmin(frat)))

#                                 # print('\n#######')
#                                 # print(frat_delta_unsigned)
#                                 # print('#######\n')
#                               #  print(frat_delta_unsigned)
#                               #  print(np.nanmin(frat))

#                                 arg_min_rat=_nanargmin(frat)
#                                 arg_max_rat=_nanargmax(frat)
#                                 ratsign=np.sign(arg_max_rat-arg_min_rat)
#                                 ratdelta=ratsign*frat_delta_unsigned

#                                 ratdelta_scaledmin=ratdelta/np.nanmin(frat)
#                                 ratdelta_scaledmax=ratdelta/np.nanmax(frat)

#                                 logratdelta=fratlog_delta_unsigned
#                               #  logratmax=np.log10(np.nanmax(frat))

#                               #  logratmax=np.nanmax(np.abs(np.log10(frat)))
#                                 # print('\n#########')
#                                 # print(logratdelta)
#                                 fratidpick=np.nanmin([arg_min_rat,arg_max_rat])
#                                 # print(arg_min_rat)
#                                 # print(arg_max_rat)
#                                 # print(fratidpick)
#                                 # print('#########\n')
#                                 if np.isnan(fratidpick)==False:
#                                  #   print(np.sign(np.log10([np.nanmax(frat),np.nanmin(frat)])))
#                                  #   print(np.unique(np.sign(np.log10([np.nanmax(frat),np.nanmin(frat)]))))
#                                  #   print(len(np.unique(np.sign(np.log10([np.nanmax(frat),np.nanmin(frat)])))))
#                                     if len(np.unique(np.sign(np.log10([np.nanmax(frat),np.nanmin(frat)]))))==1:
#                                         logratmax=np.nanmax(np.abs(np.log10(frat)))
#                                     else:
#                                         logratmax=np.nanmax(np.abs(np.log10(frat)))
#                                         #logratmax=np.abs(np.log10(frat[fratidpick]))
#                                 else:
#                                     logratmax=np.nan
#                                 logratshow=logratdelta*ratsign/logratmax

#                                 # print('\n#######')
#                                 # print(ratdelta_scaled)
#                                 # print('#######\n')
#                                 # if np.abs(ffdelta)>0.5:
#                                 #     print('\n#########')
#                                 #     print('fratio is', frat)
#                                 #     print('fratio delta is', ratdelta)
#                                 #     print('fratio delta scaled by maximum is', ratdelta/np.nanmax(frat))
#                                 #     print('fratio delta scaled by minimum is', ratdelta/np.nanmin(frat))
#                                 #     print('fratio where the min/max is first attained is', frat[np.min([arg_max_rat,arg_min_rat])])
#                                 #     print('fratio scaled by this is', ratdelta/frat[np.min([arg_max_rat,arg_min_rat])])
#                                 #     print('fdelta is', ffdelta)
#                                 #     print('#########\n')

#                                 #     print('\n#########')
#                                 #     print('log10 fratio is', np.log10(frat))
#                                 #     print('log10 maxfratio is', np.log10(np.nanmax(frat)) )
#                                 #     print('log10 minfratio is', np.log10(np.nanmin(frat)) )
#                                 #     print('log10 maxfratio - log10 minfratio is',fratlog_delta_unsigned )
#                                 #     print('test', fratlog_delta_unsigned/np.log10(np.nanmax(frat)) )
#                                 #     print('sign is', ratsign)
#                                 #     print('logratshow is', logratshow)
#                                 #     print('#########\n')

#                                 # frdelta=ff_dat[0]
#                                 # # brdelta=ff_dat[1]
#                                 # frgrad=ff_dat[2]
#                                 # frsense=ff_dat[3]

#                                 frdelta_mat[idy,idx]=ratdelta
#                                 frdelta_scaledmin_mat[idy,idx]=ratdelta_scaledmin
#                                 frdelta_scaledmax_mat[idy,idx]=ratdelta_scaledmax
#                                 frdelta_useful[idy,idx]=logratshow

#               #                  ffdelta_mat[idy,idx]=ffdelta
#                             #     # brdelta_mat[idy,idx]=brdelta
#                             #     frgrad_mat[idy,idx]=frgrad
#                             #     frsense_mat[idy,idy]=frsense
#                             # #    print(ffdelta,'\n')

#                         # print('\n########')
#                         # print(frdelta_scaledmin_mat)
#                         # print(frdelta_scaledmax_mat)
#                         # print('########\n')

#                         # print('\n########')
#                         # print(np.nanmin(frdelta_mat[:]),np.nanmax(frdelta_mat[:]))
#                         # print(np.nanmin(frdelta_scaledmin_mat[:]),np.nanmax(frdelta_scaledmin_mat[:]))
#                         # print(np.nanmin(frdelta_scaledmax_mat[:]),np.nanmax(frdelta_scaledmax_mat[:]))
#                         # print('########\n')

#                         print('\n########')
#                    #     print(frdelta_useful)
#                         print('########\n')

#                         np.savez(savename_fr,frd=frdelta_mat, frdmin=frdelta_scaledmin_mat,
#                                  frdmax=frdelta_scaledmax_mat,
#                                  frd_use=frdelta_useful)
#                         # np.save(savename_fr,frdelta_scaledmin_mat)
#                         # np.save(savename_fr,frdelta_scaledmax_mat)
#                         # np.save(,frdelta_useful)


#                      # #   np.save(savename_br,brdelta_mat)
#                      #    np.save(savename_frg,frgrad_mat)
#                      #    np.save(savename_frs,frsense_mat)

# make_heatmaps_ratio(xvar='Btot',xvarvals=np.logspace(-3,3,25),
#                   yvar='CKcatB',yvarvals=np.logspace(-3,3,25),
#                   lengthparams=[2,2,2],B_position_params=([1,0],[0,0],[0,1],[1,0]),
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
#                   variable_param_values=([0.01],
#                                          [0.1],
#                                          [1.],
#                                          [10.],
#                                          [1.]),
#                   kin_min=-8,kin_max=0,num_kin=41,
#                   reversible=False,Tmax=100000,tpts=10000)
   # print(xvar_dict)
   # print(yvar_dict)                        #variable_name0=var_vals0
                       # print(variable_name0)
  #  print(variable_dict)

 #   print(KeqB)
    # for idlab, lab in enumerate(variable_params):
    #     print(idlab,lab)

    #     if lab == 'CKcatB':
    #         CKcatB_vals=variable_param_values[idlab]
    #         print(CKcatB_vals)

    #     elif lab == 'KeqB':
    #         KeqB_vals=variable_param_values[idlab]
    #         print(KeqB_vals)

    #     elif lab == 'Btot':
    #         Btot_vals=variable_param_values[idlab]
    #         print(Btot_vals)

    #     elif lab == 'CKcat_top':
    #         CKcat_top_vals=variable_param_values[idlab]
    #         print(CKcat_top_vals)

    #     elif lab == 'Keq_top':
    #         Keq_top_vals=variable_param_values[idlab]
    #         print(Keq_top_vals)

    #     elif lab == 'CKcat_bottom':
    #         CKcat_bottom_vals=variable_param_values[idlab]
    #         print(CKcat_bottom_vals)

    #     elif lab == 'Keq_bottom':
    #         Keq_bottom_vals=variable_param_values[idlab]
    #         print(Keq_bottom_vals)

    #     else:
    #         print('cannot find label')
    #     all_variable_values.remove(lab)
    # print(all_variable_values)





    # filename='log_variable_params_'+\
    #     rd_str([CKcatB,KeqB,Btot,CKcat_top,Keq_top,CKcat_bottom,Keq_bottom])+\
    #     '__kinvals_'+str([kin_min,kin_max,num_kin])+\
    #     '__reversible_'+str(reversible)+\
    #     '__timevals_'+str([Tmax,tpts])+'.npy'

    # if 'CKcatB_vals' in locals():
    #     print('yes')
    #     for CKcatB in CKcatB_vals:


  #  for id0, var0 in enumerate(var0_vals):
  #      print(id0,var0_vals[id0],var0,lab0)

#            params_global=[1.,1.,1.,1.],params_B=[0.1,1.,1.,1.,10.],
#            params_branch=[1.,1.,1.,1.,
#                           1.,1.,1.,1.],
#            kouts=[1.,1.],#params_inout=[0.1,1.,1.],
#            kinvals=np.logspace(-8,0,41),
#            reversible=False,Tmax=1000000,tpts=100000):



# def make_heatmap_data(xvar='Btot',yvar='KeqB'




#         kinvals=np.logspace(-8,0,41),
#                  xvar='Btot',yvar='KeqB',
#                  xvar_vals=np.logspace(-2,2,11),yvar_vals=np.logspace(-2,2,11),
#                  defaultparams=['CKcat','Ks','Kp','Keq',
#                                 'KsB','KpB',  #CKcatB,KsB,KpB,KeqB,Btot
#                                 'Ks_top','Kp_top',  #CKcat_top,   Ks_top,   Kp_top,   Keq_top,
#                                 'Ks_bottom','Kp_bottom', #CKcat_bottom,Ks_bottom,Kp_bottom,Keq_bottom]
#                                 'kout_top','kout_bottom'],
#                  changedparams=['CKcatB','CKcat_top','Keq_top','CKcat_bottom','Keq_bottom'],
#                  changedparam_vals=[0.1,10.,1.,0.1,1.],
#                  lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
#                  reversible=False,Tmax=100000,tpts=10000):

#     params_global=np.zeros(4)
#     params_B=np.zeros(5)
#     params_branch=np.zeros(8)
#     kouts=np.zeros(2)

#     params_global[:]=np.nan
#     params_B[:]=np.nan
#     params_branch[:]=np.nan
#     kouts[:]=np.nan


