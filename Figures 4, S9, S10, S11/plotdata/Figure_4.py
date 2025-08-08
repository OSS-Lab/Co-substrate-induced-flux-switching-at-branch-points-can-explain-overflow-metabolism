#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 11:09:05 2025

@author: robert
"""




import numpy as np
# import time
# from datetime import timedelta
# import multiprocessing as multi
import sys

#from pathlib import Path
# import os
# import time
from pathlib import Path

# import itertools

# sys.path.append('../data/megadata')
# sys.path.append('../data/heatmapdata')
sys.path.append('../simulation_functions')
sys.path.append('../data')


import make_heatmaps_compartment as mh_cm
import ode_funcs_compartment as ode_cm
import matplotlib.image as img

import matplotlib.pyplot as plt
# import numpy.random as rn
# import datetime as dt

import matplotlib.pylab as pylab
import matplotlib as mpl

savefig=0       # change to 1 to save generated figures
# numkin=201
numkin=81
# 
params = {'legend.fontsize': 12,
          'legend.title_fontsize': 16,
          'axes.labelsize': 20,
          'axes.titlesize': 20,
          'xtick.labelsize':20,
          'ytick.labelsize':20}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)



#202 421 445 454 487 570 874
indplot_switchretained=697 #570

indplot_noncosub=933
plt.figure()
# ax1=plt.axes([-0.15,-0.1,1,0.6])       # motif
ax1=plt.axes([1,0.5,1.5,1])       # motif

# motifnone=img.imread('motif_compartment_orig.png')
#motifnone=img.imread('motif_compartment_sonal.png')
motifnone=img.imread('motif_fig4.png')
# motifnone=img.imread('yeast_model_s.png')
plt.imshow(motifnone)
plt.axis('off')

ax2=plt.axes([1,-0.7,0.95,1.2])    # scatter
# ax3=plt.axes([-0.04,-0.65,0.4,0.48])    # fluxfraction
# ax4=plt.axes([0.44,-0.65,0.4,0.48])    # expt
ax4=plt.axes([2.1,-0.1,0.9,0.6])    # fluxfraction
ax11=plt.axes([2.1,-0.7,0.9,0.6])    # fluxfraction

ax3=plt.axes([3.05,-0.1,0.9,0.6])    # expt
ax10=plt.axes([3.05,-0.7,0.9,0.6])    # fluxfraction
ax13=plt.axes([3.05,-0.7,0.9,0.6])    # fluxfraction
ax13.set_xlim((0,1))
ax13.set_ylim((0,1))
ax13.axis('off')
# ax3=plt.axes([2,-0.6,0.9,1.1])    # fluxfraction
# ax4=plt.axes([3,-0.6,0.9,1.1])    # expt
#ax5=plt.axes([-0.1,-0.6,0.65,0.5])
ax5=plt.axes([2.65,1.05,0.5,0.35])
ax12=plt.axes([2.65,0.7,0.5,0.35])
# ax5.text(0.5,0.5,'Expt.',fontsize=30)
# ax5=plt.axes([0.05,-0.6,0.8,0.5])   # box plots
# ax6=plt.axes([0.65,-0.6,0.25,0.25])
# ax7=plt.axes([0.65,-0.35,0.25,0.25])

ax6=plt.axes([3.3,0.7,0.35,0.35])
ax7=plt.axes([3.3,1.05,0.35,0.35])

ax8=plt.axes([3.3,0.7,0.35,0.35])
ax9=plt.axes([3.3,1.05,0.35,0.35])
ax8.axis('off')
ax9.axis('off')

# ax8=plt.axes([3.85,0.7,0.35,0.35])
# ax9=plt.axes([3.85,1.05,0.35,0.35])

# ax5.text(0.2,0.45,'Experimental data',fontsize=30)

# ax2.text(-1.55,0.85,'$\\rm A$',fontsize=40)
ax2.text(-0.15,1.8,'$\\rm A$',fontsize=40)
ax2.text(1.725,1.8,'$\\rm B$',fontsize=40)
ax2.text(2.5,1.8,'$\\rm C$',fontsize=40)
ax2.text(-0.15,0.85,'$\\rm D$',fontsize=40)
ax3.text(10**-8,0.8,'$\\rm F$',fontsize=40,color='darkviolet')
ax4.text(10**-8,0.8,'$\\rm E$',fontsize=40,color='fuchsia')
# ax3.text(10**(-10.4),0.75,'$\\rm C$',fontsize=30)
# ax4.text(10**(-10.5),0.75,'$\\rm D$',fontsize=30)


# ax5=plt.axes([2.0,0.02,0.2,0.48])   # box plots
# ax6=plt.axes([2.325,0.02,0.2,0.48])
# ax7=plt.axes([2.65,0.02,0.2,0.48])

# ax8=plt.axes([2.0,-0.67,0.2,0.48])   # box plots
# ax9=plt.axes([2.325,-0.67,0.2,0.48])
# ax10=plt.axes([2.65,-0.67,0.2,0.48])
# ax11=plt.axes([2.975,-0.67,0.2,0.48])
# ax12=plt.axes([3.3,-0.67,0.2,0.48])

# get the data for new params
allo=[1.9,1.,1.,1.]
n=1000

seed=1000
np.random.seed(seed)
if seed!=False:
    seedstr=str(seed)
else:
    seedstr='none'

delwant_max=1
delwant_min=0.5
n0deltawant_max=1.
n0deltawant_min=0.5
delwant_noswitch=0.4

args_to_run=mh_cm.make_params_newtest_6_nonfixed(allo, n)

      
savepath='../data/allo_'+str(allo)

#filename='nonfixed_new_3_seed'+seedstr+'_n_'+str(n)
filename='nonfixed_apr_final_seed'+seedstr+'_n_'+str(n)+'_datetime'
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewant=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewant)

kinvals=args_to_run[0][1]

dat=np.load(filewant)
ff_delta=dat['fluxfrac_delta']
gradmaxind=dat['gradmaxind']
sensemaxind=dat['sensemaxind']
fluxfrac=dat['fluxfrac']
ff_grad=dat['fluxfrac_grad']
sensemax=dat['sensemax']
gradmax=dat['gradmax']
ffsense_vec=dat['fluxfrac_sense']
b0delta=dat['b0delta']
c0delta=dat['c0delta']
n0delta=dat['n0delta']
n0delta=dat['n0delta']
swpt=dat['swpt']
swpt_orig=dat['swpt']
swptvals=dat['swptvals']
swptvals_orig=dat['swptvals']

para01_v=list()
para01_ks=list()
para01_kp=list()
para01_keq=list()


para1E_v=list()
para1E_ks=list()
para1E_kp=list()
para1E_keq=list()

paraE2_v=list()
paraE2_ks=list()
paraE2_kp=list()
paraE2_keq=list()

para23_v=list()
para23_ks=list()
# para23_kp=list()
# para23_keq=list()

para35_v=list()
para35_ks=list()
para35_kp=list()
para35_keq=list()

para_24_v=list()
para_24_ks=list()
# para_24_kp=list()
# para_24_keq=list()

para_35_v=list()
para_35_ks=list()
para_35_kp=list()
para_35_keq=list()

para_t_v=list()
para_t_ks=list()

para_diff_3_in=list()
para_diff_3_out=list()
para_diff_5_in=list()
para_diff_5_out=list()

paraB_v=list()
paraB_ks=list()
paraB_kp=list()
paraB_keq=list()

paraC_v=list()
paraC_ks=list()
paraC_kp=list()
paraC_keq=list()

btot=list()
ctot=list()

btot_p=list()
ctot_p=list()


para01=list()
para1E=list()
paraE2=list()
para23=list()
para35=list()
para_24=list()
para_35=list()

para_t=list()
para_diff_3=list()
para_diff_5=list()

paraB=list()
paraC=list()

ff_delta=np.abs(ff_delta)
c0delta=np.abs(c0delta)
b0delta=np.abs(b0delta)
n0delta=np.abs(n0delta)


indexes_switch_retained=list()

for i in range(n):

    # if ff_delta[i]<delwant_min:
    #     if ff_delta[i]>delwant_max:
    pinit=args_to_run[i][0]
    kinvals=args_to_run[i][1]
    params_out=args_to_run[i][2]
    
    params01=args_to_run[i][3]
    params1E=args_to_run[i][4]
    paramsE2=args_to_run[i][5]
    params23=args_to_run[i][6]
    params35=args_to_run[i][7]
    params_24=args_to_run[i][8]
    
    params_35=args_to_run[i][9]
    params_trans=args_to_run[i][10]
    params_diff_3=args_to_run[i][11]
    params_diff_5=args_to_run[i][12]
    
    paramsB=args_to_run[i][13]
    paramsC=args_to_run[i][14]
    allo=args_to_run[i][15]
    Tmax=args_to_run[i][16]
    
    # print(params35)
    # print(params_35)
    # if ff_delta[i]>delwant_min and ff_delta[i]<delwant_max and \
    #     c0delta[i]>n0deltawant_min and c0delta[i]<n0deltawant_max and \
    #     n0delta[i]>n0deltawant_min and n0delta[i]<n0deltawant_max:
            # print(i)
          #  count+=1
            # print(args_to_run[i])
   
        
        #tpts=args_to_run[i][16]
    
    #if ff_delta[i]>0.5 and ff_delta[i]<1.0 and \
        # c0delta[i]>0.3 and c0delta[i]<1.0: 
    if ff_delta[i]>0.5 and c0delta[i]>0.2: 
    # if i in allo19_desiredafteroxidase and n0delta[i]>0.2:        
        para01.append(params01)
        para1E.append(params1E)
        paraE2.append(paramsE2)
        para23.append(params23)
        para35.append(params35)
    
        para_24.append(params_24)
        para_35.append(params_35)
    
        para_t.append(params_trans)
        para_diff_3.append(params_diff_3)
        para_diff_5.append(params_diff_5)
        
        paraB.append(paramsB)
        paraC.append(paramsC)
    
        ctot.append(pinit[-3]+pinit[-2])
        btot.append(pinit[-5]+pinit[-4])
        
        
        para01_v.append(np.log10(1000*params01[0]))
        para01_ks.append(np.log10(params01[1]))
        para01_kp.append(np.log10(params01[2]))
        para01_keq.append(np.log10(params01[3]))
        
        para1E_v.append(np.log10(1000*params1E[0]))
        para1E_ks.append(np.log10(params1E[1]))
        para1E_kp.append(np.log10(params1E[2]))
        para1E_keq.append(np.log10(params1E[3]))
        
        paraE2_v.append(np.log10(1000*paramsE2[0]))
        paraE2_ks.append(np.log10(paramsE2[1]))
        paraE2_kp.append(np.log10(paramsE2[2]))
        paraE2_keq.append(np.log10(paramsE2[3]))
        
        para23_v.append(np.log10(1000*params23[0]))
        para23_ks.append(np.log10(params23[1]))
        # para23_kp.append(np.log10(params23[2]))
        # para23_keq.append(np.log10(params23[3]))
        
        para35_v.append(np.log10(1000*params35[0]))
        para35_ks.append(np.log10(params35[1]))
        para35_kp.append(np.log10(params35[2]))
        para35_keq.append(np.log10(params35[3]))
        
        para_24_v.append(np.log10(1000*params_24[0]))
        para_24_ks.append(np.log10(params_24[1]))
        # para_24_kp.append(np.log10(params_24[2]))
        # para_24_keq.append(np.log10(params_24[3]))
        
        para_35_v.append(np.log10(1000*params_35[0]))
        para_35_ks.append(np.log10(params_35[1]))
        para_35_kp.append(np.log10(params_35[2]))
        para_35_keq.append(np.log10(params_35[3]))
        
        para_t_v.append(np.log10(1000*params_trans[0]))
        para_t_ks.append(np.log10(params_trans[1]))
        
        para_diff_3_in.append(np.log10(params_diff_3[0]))
        para_diff_3_out.append(np.log10(params_diff_3[1]))
        para_diff_5_in.append(np.log10(params_diff_5[0]))
        para_diff_5_out.append(np.log10(params_diff_5[1]))
        
        paraB_v.append(np.log10(1000*paramsB[0]))
        paraB_ks.append(np.log10(paramsB[1]))
        paraB_kp.append(np.log10(paramsB[2]))
        paraB_keq.append(np.log10(paramsB[3]))
        
        paraC_v.append(np.log10(1000*paramsC[0]))
        paraC_ks.append(np.log10(paramsC[1]))
        paraC_kp.append(np.log10(paramsC[2]))
        paraC_keq.append(np.log10(paramsC[3]))
        
        btot_p.append(np.log10(pinit[-5]+pinit[-4]))
        ctot_p.append(np.log10(pinit[-3]+pinit[-2]))
    

        indexes_switch_retained.append(i)

        
#ax2.scatter(np.abs(ff_delta),np.abs(n0delta),s=40,color='lightpink',alpha=0.5,edgecolors='k',zorder=10)#,label='')
ax2.scatter(ff_delta,c0delta,s=45,color='bisque',edgecolors='k',zorder=0)#,label=''),alpha=0.5
ax2.scatter(ff_delta[indexes_switch_retained],c0delta[indexes_switch_retained],s=45,color='teal',edgecolors='teal',zorder=0,label='${\\rm co}$-${\\rm sub.~mediated~switching}$')#,label=''),alpha=0.5
ax2.scatter(ff_delta[indplot_switchretained],c0delta[indplot_switchretained],s=45,color='darkviolet',alpha=1,edgecolors='darkviolet',zorder=20)#,label='')
ax2.scatter(ff_delta[indplot_noncosub],c0delta[indplot_noncosub],s=45,color='fuchsia',alpha=1,edgecolors='fuchsia',zorder=20)#,label='')
leg=ax2.legend(loc=(0.2,0.935))
leg.get_frame().set_alpha(None)
leg.get_frame().set_facecolor((1, 1, 1, 0))

# ax2.scatter(ff_delta[indexes_noswitch],n0delta[indexes_noswitch],s=40,color='bisque',alpha=0.5,edgecolors='k',zorder=20)#,label='')
# ax2.scatter(ff_delta[indexes_switch],n0delta[indexes_switch],s=40,color='lightblue',alpha=1,edgecolors='k',zorder=10)#,label='')
# ax2.scatter(ff_delta[indexes_switch_retained],n0delta[indexes_switch_retained],s=40,color='limegreen',alpha=1,edgecolors='k',zorder=10)#,label='')
# ax2.scatter(ff_delta[ind_notpicked],n0delta[ind_notpicked],s=40,color='bisque',alpha=0.5,edgecolors='k',zorder=0)#,label='')
# ax2.scatter(ff_delta[indplot_switchretained],n0delta[indplot_switchretained],s=40,color='red',alpha=1,edgecolors='k',zorder=20)
ax2.axis([-0.05,1.05,-0.05,1.05])
# ax2.axis([-1,1,-1,1])
ax2.set_xlabel('${\\rm Change ~in ~flux ~fraction,~} \\Delta F_f$',fontsize=30)
# ax2.set_ylabel('${\\rm Change ~in ~mitochondrial}$ \n ${\\rm NADH~ fraction,~} \\Delta N_{\\rm H}^{m}$',fontsize=30,rotation=90)
ax2.set_ylabel('${\\rm Change ~in ~mitochondrial~ NADH~ fraction,~} \\Delta N_{\\rm H}^{m}$',fontsize=20,rotation=90)
ax2.yaxis.set_label_coords(-0.09,0.5)

ax2.text(0.95,c0delta[indplot_switchretained],'$\\rm F$',fontsize=30,color='darkviolet')
ax2.text(1,c0delta[indplot_noncosub],'$\\rm E$',fontsize=30,color='fuchsia')# ax3.text(10**-7,0.8,str(indplot_switchretained),fontsize=20)
# ax2.set_xticks(np.logspace(-9,-2,8), ['$10^{-9}$','','$10^{-7}$','',
#                                      '$10^{-5}$','','$10^{-3}$',''])
ind_temp=indplot_switchretained
args_temp=args_to_run[ind_temp]
print('\n')
print(ind_temp)
print(args_temp[0])
print(args_temp[1])
print(args_temp[2])
print(args_temp[3])
print('\n')
print(args_temp[4])
print(args_temp[5])
print(args_temp[6])
print(args_temp[7])
print('\n')
print(args_temp[8])
print(args_temp[9])
print(args_temp[10])
print(args_temp[11])
print('\n')
print(args_temp[12])
print(args_temp[13])
print('\n')

print(args_temp[4])
print(args_temp[6])



kinplot=np.logspace(-8,0,numkin)

result=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],args_temp[13],args_temp[14],args_temp[15],
                             args_temp[16],args_temp[17])

    
# flux01= ode_cm.flux_enz_forward(result[:,0], result[:,1], result[:,9], result[:,10], 
#                          args_temp[3][0], args_temp[3][1], args_temp[3][2])            # gap -> pep

# flux12= ode_cm.flux_enz_forward(result[:,1], result[:,2], 1.,1., 
#                          args_temp[4][0], args_temp[4][1], args_temp[4][2])           # pep -> pyr

flux23= ode_cm.flux_enz_forward(result[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])            # pyr -> acet

# flux35= ode_cm.flux_enz_forward(result[:,3], result[:,4], result[:,10],result[:,9], 
#                          args_temp[6][0], args_temp[6][1], args_temp[6][2])              # acet -> etn

# flux22=ode_cm.transp(result[:,2],args_temp[9][0],args_temp[9][1])             # pyr -> pyr_mito

flux24= ode_cm.flux_enz_forward(result[:,5],  result[:,11], 
                             args_temp[8][0], args_temp[8][1])

ax10.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Whole cell NADH Fraction',linewidth=2,color='firebrick')
ax10.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
ax10.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
ax3.plot(kinplot,1-flux24/(flux23+flux24),linewidth=2,color='royalblue')#,label='Flux Fraction'
ax3.text(10**-7,0.85,'$\\rm \\underline{Model}$',fontsize=20)

# ax10.text(10**-7.5,0.8,'mito.', fontsize=25, color='mediumseagreen')
# ax10.text(10**-4.5,0.5,'cyto.', fontsize=25, color='peru')
# ax10.text(10**-6,0.5,'total.', fontsize=25, color='firebrick')
# ax3.text(10**-6,0.6,'$F_f$', fontsize=20, color='royalblue')


ax10.set_xlabel('$v_{\\rm in}$',fontsize=40)
# ax3.set_ylabel('$ F_f$',fontsize=30)
ax4.set_ylabel('$ {\\rm Flux ~Fraction}, F_f$',fontsize=20)
ax11.set_ylabel('$ {\\rm NADH ~Fraction}, N_{\\rm H}$',fontsize=20)
# ax3.legend(loc='upper left',ncols=1)
ax3.set_xscale('log')
ax10.set_xscale('log')

ax3.set_xticks(np.logspace(-8,-4,5),labels=['','','','',''])#, ['$10^{-9}$','','$10^{-7}$','',
                                    # '$10^{-5}$','','$10^{-3}$',''])

ax10.set_xticks(np.logspace(-8,-4,5))#, ['$10^{-9}$','','$10^{-7}$','',
ax3.set_yticks(np.linspace(0,1,3),labels=['','',''])
ax3.set_ylim((-0.05,1.05))
ax10.set_ylim((-0.05,1.05))
ax10.set_yticks(np.linspace(0,1,3),labels=['','',''])
ax10.xaxis.set_label_coords(0.7,-0.1)

ax3.set_xlim((10**-(8*1.025),10**-(4*0.975)))
ax10.set_xlim((10**-(8*1.025),10**-(4*0.975)))

ax3.yaxis.set_label_coords(-0.05,0.5)
ax10.yaxis.set_label_coords(-0.05,0.5)
ax10.xaxis.set_label_coords(0.7,-0.1)

# print(swpt[ind_temp])
# print(fluxfrac[ind_temp])
# ax3.scatter(swpt[ind_temp],fluxfrac[ind_temp],'o',filled=True,size=20)

# ff_new=1-flux24_4/(flux23_4+flux24_4)
swpt_4,sw_orig=mh_cm.find_sw_pt(1-flux24/(flux23+flux24),kinplot)
ax3.scatter(swpt_4,sw_orig,60, color='royalblue',edgecolor='royalblue',linewidth=2,
            label='$\\rm Without~ NADH~Oxidase$')



  

# args_temp[12][3]*=0.1
# args_temp[13][3]*=0.1
args_temp[14][3]*=0.1
# args_temp[12][0]*=10.
# args_temp[13][0]*=10.
print('########\n check')
print(args_temp[13])
print(args_temp[14])

result2=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],args_temp[13],args_temp[14],args_temp[15],
                             args_temp[16],args_temp[17])

  
flux23_2= ode_cm.flux_enz_forward(result2[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])             # pyr -> acet

flux24_2= ode_cm.flux_enz_forward(result2[:,5],  result2[:,11], 
                             args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


ax10.plot(kinplot,1-(result2[:,9]+result2[:,11])/(result2[:,9]+result2[:,10]+result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='firebrick')
ax10.plot(kinplot,1-(result2[:,9])/(result2[:,9]+result2[:,10]),linestyle='--',linewidth=2,color='peru')
ax10.plot(kinplot,1-(result2[:,11])/(result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
ax3.plot(kinplot,1-flux24_2/(flux23_2+flux24_2),linestyle='--',linewidth=2,color='royalblue')

ax10.text(10**-6.025,0.3,'${\\rm NADH}$ \n ${\\rm oxidase}$',color='k',fontsize=15)

sw_o=sw_orig
swpt_kc,sw_orig_n=mh_cm.find_sw_pt(1-flux24_2/(flux23_2+flux24_2),kinplot,sw_orig)
ax3.scatter(swpt_kc,sw_orig,60, color='none',edgecolor='royalblue',linewidth=2,
            label='$\\rm With~NADH~oxidase$')
 
ax3.legend(title='$\\rm \\underline{Switch~ Points}$')
# swpt_kc,sw_orig=mh_cm.find_sw_pt(1-flux24_2/(flux23_2+flux24_2),kinplot)
# ax3.scatter(swpt_kc,sw_orig,60, color='none',edgecolor='royalblue',linewidth=2)

ax13.arrow(0.45, 0.5, 0.225,0,length_includes_head=True,
          width=0.005,head_width=0.05,head_length=0.03,
          color='k')#, **kwargs)

swpt_kc_newtest,swpt_kc_vals_newtest=mh_cm.find_sw_pt(1-flux24_2/(flux23_2+flux24_2), kinplot,swptvals_orig[ind_temp])
print(swpt_kc_newtest,swpt_kc_vals_newtest)
# ax3.scatter(swpt_kc_newtest,swpt_kc_vals_newtest,60,'k',)
# # args_temp[12][3]*=0.1
# # args_temp[13][3]*=0.1
# args_temp[14][3]*=100
# # args_temp[12][0]*=10.
# # args_temp[13][0]*=10.
# print('########\n check')
# print(args_temp[13])
# print(args_temp[14])

# result2=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
#                              args_temp[4],args_temp[5],args_temp[6],args_temp[7],
#                              args_temp[8],args_temp[9],args_temp[10],args_temp[11],
#                              args_temp[12],args_temp[13],args_temp[14],args_temp[15],
#                              args_temp[16],args_temp[17])

  
# flux23_2= ode_cm.flux_enz_forward(result2[:,2]**allo[0], 1.,
#                              args_temp[6][0],  args_temp[6][1])             # pyr -> acet

# flux24_2= ode_cm.flux_enz_forward(result2[:,5],  result2[:,11], 
#                              args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


# ax3.plot(kinplot,1-(result2[:,9]+result2[:,11])/(result2[:,9]+result2[:,10]+result2[:,11]+result2[:,12]),linestyle=':',linewidth=2,color='firebrick')
# ax3.plot(kinplot,1-(result2[:,9])/(result2[:,9]+result2[:,10]),linestyle=':',linewidth=2,color='peru')
# ax3.plot(kinplot,1-(result2[:,11])/(result2[:,11]+result2[:,12]),linestyle=':',linewidth=2,color='mediumseagreen')
# ax3.plot(kinplot,1-flux24_2/(flux23_2+flux24_2),linestyle=':',linewidth=2,color='royalblue')

# boxwidth=0.8

# ax5.boxplot([paraC_v , paraC_v_n],widths=(boxwidth,boxwidth))
                  
# # ax5.set_title('Total cosubstrate (Mito)', fontsize=20)
# ax5.set_title('$V_{\\rm M}^{\\rm N}~ {\\rm (Mito)}$', fontsize=15,y=1.01)

# ax5.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# ax5.set_yticks(ticks=[-1,-0], labels=['$10^{-1}$','$10^{0}$'])
# ax5.grid('on')



# ax6.boxplot([paraC_ks, paraC_ks_n], widths=(boxwidth,boxwidth))
                  
# # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
# ax6.set_title('$K_{\\rm S}^{N_{+}}~{\\rm (Mito)}$', fontsize=15,y=1.01)

# ax6.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# ax6.set_yticks(ticks=[-1,0,1,2], labels=['$10^{-1}$','$10^{0}$','$10^{1}$','$10^{2}$'])
# # ax6.set_yticks(ticks=np.linspace(-4,0))

# ax6.grid('on')



# ax7.boxplot([paraC_kp, paraC_kp_n],widths=(boxwidth,boxwidth))
                  
# # ax6.set_title('$K_{\\rm S}$, NAD+ $\\leftrightarrow$ NADH (Mito)', fontsize=20)
# ax7.set_title('$K_{\\rm P}^{N_H}~{\\rm (Mito)}$', fontsize=15,y=1.01)

# ax7.set_xticklabels(labels=['S','N'], 
#                        rotation=0)
# ax7.set_yticks(ticks=[-3,-2,-1,0], labels=['$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'])

# ax7.grid('on')

# get the data for original params
ind_temp=indplot_noncosub
args_temp=args_to_run[ind_temp]
print('\n')
print(ind_temp)
print(args_temp[0])
print(args_temp[1])
print(args_temp[2])
print(args_temp[3])
print('\n')
print(args_temp[4])
print(args_temp[5])
print(args_temp[6])
print(args_temp[7])
print('\n')
print(args_temp[8])
print(args_temp[9])
print(args_temp[10])
print(args_temp[11])
print('\n')
print(args_temp[12])
print(args_temp[13])
print('\n')

print(args_temp[4])
print(args_temp[6])



kinplot_2=np.logspace(-8,-3,numkin)

result=mh_cm.get_flux_vecs_2(args_temp[0],kinplot_2,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],args_temp[13],args_temp[14],args_temp[15],
                             args_temp[16],args_temp[17])

    
# flux01= ode_cm.flux_enz_forward(result[:,0], result[:,1], result[:,9], result[:,10], 
#                          args_temp[3][0], args_temp[3][1], args_temp[3][2])            # gap -> pep

# flux12= ode_cm.flux_enz_forward(result[:,1], result[:,2], 1.,1., 
#                          args_temp[4][0], args_temp[4][1], args_temp[4][2])           # pep -> pyr

flux23= ode_cm.flux_enz_forward(result[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])            # pyr -> acet

# flux35= ode_cm.flux_enz_forward(result[:,3], result[:,4], result[:,10],result[:,9], 
#                          args_temp[6][0], args_temp[6][1], args_temp[6][2])              # acet -> etn

# flux22=ode_cm.transp(result[:,2],args_temp[9][0],args_temp[9][1])             # pyr -> pyr_mito

flux24= ode_cm.flux_enz_forward(result[:,5],  result[:,11], 
                             args_temp[8][0], args_temp[8][1])

ax11.plot(kinplot_2,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='$\\rm Whole ~Cell$',linewidth=2,color='firebrick')
ax11.plot(kinplot_2,1-(result[:,9])/(result[:,9]+result[:,10]),label='$\\rm Cytosolic$',linewidth=2,color='peru')
ax11.plot(kinplot_2,1-(result[:,11])/(result[:,11]+result[:,12]),label='$\\rm Mitochondrial$',linewidth=2,color='mediumseagreen')
ax4.plot(kinplot_2,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')
ax11.legend(loc='upper left',ncols=1,fontsize=15,title='$\\frac{\\rm [NADH]}{\\rm [NADH]+[NAD^+]}$',title_fontsize=20)

ax11.set_xlabel('$v_{\\rm in}$',fontsize=40)
# ax4.set_ylabel('$ \\rm Flux ~Fraction $',fontsize=30)
ax4.set_xscale('log')
ax11.set_xscale('log')
# ax4.legend(loc='upper left',ncols=1)
ax4.set_xticks(np.logspace(-8,-4,5),labels=['','','','',''])#, ['$10^{-9}$','','$10^{-7}$','',
                                    # '$10^{-5}$','','$10^{-3}$',''])

ax11.set_xticks(np.logspace(-8,-4,5))#, ['$10^{-9}$','','$10^{-7}$','',
ax4.set_yticks(np.linspace(0,1,3),labels=['$0$','','$1$'])
ax4.set_ylim((-0.05,1.05))
ax11.set_ylim((-0.05,1.05))
ax4.set_xlim((10**-(8*1.025),10**-(4*0.975)))
ax11.set_xlim((10**-(8*1.025),10**-(4*0.975)))

ax11.set_yticks(np.linspace(0,1,3),labels=['$0$','','$1$'])
ax4.set_xlim((10**-(8*1.025),10**-(4*0.975)))
ax11.set_xlim((10**-(8*1.025),10**-(4*0.975)))

ax4.set_ylabel('$ {\\rm Flux ~Fraction}, F_f$',fontsize=20)
ax11.set_ylabel('$ {\\rm NADH ~Fraction}, N_{\\rm H}$',fontsize=20)

ax4.text(10**-7,0.85,'$\\rm \\underline{Model}$',fontsize=20)

# args_temp[12][3]*=0.1
args_temp[14][3]*=0.1
# args_temp[14][3]*=0.1
# args_temp[12][0]*=10.
# args_temp[13][0]*=10.
print('########\n check')
print(args_temp[13])
print(args_temp[14])
result2=mh_cm.get_flux_vecs_2(args_temp[0],kinplot_2,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],args_temp[13],args_temp[14],args_temp[15],
                             args_temp[16],args_temp[17])

  
flux23_2= ode_cm.flux_enz_forward(result2[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])             # pyr -> acet

flux24_2= ode_cm.flux_enz_forward(result2[:,5],  result2[:,11], 
                             args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


ax11.plot(kinplot_2,1-(result2[:,9]+result2[:,11])/(result2[:,9]+result2[:,10]+result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='firebrick')
ax11.plot(kinplot_2,1-(result2[:,9])/(result2[:,9]+result2[:,10]),linestyle='--',linewidth=2,color='peru')
ax11.plot(kinplot_2,1-(result2[:,11])/(result2[:,11]+result2[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
ax4.plot(kinplot_2,1-flux24_2/(flux23_2+flux24_2),linestyle='--',linewidth=2,color='royalblue')






xvals=[0.06,0.06,0.06,0.06,0.11,0.11,0.11,0.11,0.16,0.16,0.16,0.16,
       0.31,0.31,0.31,0.31,0.24,0.24,0.24,0.24]

yvals_NADplus=[167512.24540,132341.72260,18742.81813,60240.29535,108574.66520,
               72207.64215,40378.42461,48398.24931,89811.10767,92739.67721,
               36268.65257,78182.84713,53368.24383,56771.21420,30889.54444,
               17953.57795,202405.80600,95370.78405,30688.82054,49933.18369]

yvals_NADH=[1080.220631,1034.355194,300.000000,300.000000,3102.862200,
            1134.519565,300.000000,428.974770,3320.129937,692.435146,
            454.754043,300.000000,7121.705525,3499.396985,300.000000,
            300.000000,1576.795062,496.231033,300.000000,300.000000]

xvals_dilution_rate=[0.06,0.11,0.16,0.23,0.31]
yvals_eth_prod_rate=[0,0,0,0.021,3.708]  #mol/litre cell volume/h

# ax5.plot(xvals_dilution_rate,yvals_eth_prod_rate,'-',linewidth=2)
lin1=ax5.plot(xvals_dilution_rate,yvals_eth_prod_rate,'-x',color='royalblue',
         ms = 8, mec = 'k', mfc = 'k',
         label='$\\rm Ethanol ~Production ~Rate$')# \n $\\rm (Experiment)$')
ax12.set_xlabel('${\\rm Dilution ~Rate}~(h^{-1})$',fontsize=20)
# ax5.set_ylabel('${\\rm Ethanol}$'+'\n'+'${\\rm Production}$'+'\n'+'${\\rm Rate~(mol/litre~cell~colume)/}h$',ha="left")
# ax5.set_ylabel('${\\rm mol/litre~cell~volume}/h$',fontsize=20)
ax5.set_yticks(np.linspace(0,4,3))
ax5.set_xticks(np.linspace(0,0.3,4),labels=['','','',''])
# ax5.text(0.05,1.5,'$\\rm \\underline{Experiment}$',fontsize=16)

ax12.set_xticks(np.linspace(0,0.3,4))

# ax12.set_yticks(np.linspace(0,1,3))
ax5.yaxis.set_label_coords(-0.1,0.5)
# ax5.legend()

dilRate=np.array([	0.05	,0.11	,0.16	,0.22	,0.3])
NAD=np.array([	0.001171528	,0.001096782	, 0.001090527,	0.001043792,	0.000807479])
NADH	=np.array([3.93345E-06	,4.55817E-06	,5.18715E-06	,8.31691E-06	,1.97706E-05])

# ax12.plot(dilRate,NADH/(NAD+NADH),'-',color='firebrick',linewidth=2)
lin2=ax12.plot(dilRate,NADH/(NAD+NADH),'-x',color='firebrick',linewidth=2,
          ms = 8, mec = 'k', mfc = 'k',
          label='$\\rm Whole ~Cell~ NADH ~Fraction$')# \n $\\rm (Experiment)$')

lns = lin1+lin2
labs = [l.get_label() for l in lns]
ax12.legend(lns, labs, loc='upper left', title='$\\rm \\underline{Experimental~ Data:}$')

# ax12.legend(loc='upper left')

ax12.set_xlim((0,0.33))
ax5.set_xlim((0.0,0.33))
ax12.set_ylim((-0.005,0.105))
ax5.set_ylim((-0.5,4.1))

mods=[[1.,1.,0.1,1.],
      [1.,1.,1.,0.1]]


filename='nonfixed_apr_final_KCmod_seed'+seedstr+'_n_'+str(n)
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewantkcmod=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewantkcmod)

dat=np.load(filewantkcmod)
ffvals_kc=dat['fluxfrac']

swpt_diffs_c=list()
print(swpt_diffs_c)

# swpt_diffs_c_a=np.zeros(len(indexes_switch_retained))
# print(swpt_diffs_c_a)
for i in range(len(indexes_switch_retained)):
    ind_temp=indexes_switch_retained[i]

# for i in range(n):
#     ind_temp=i
    
    # swpt_orig=swpt[ind_temp]
    # swpt_val_orig=swptvals[ind_temp]
    ffvals_kc_temp=1-ffvals_kc[ind_temp]
    ffvals_orig=1-fluxfrac[ind_temp]
    # print(ffvals_kc_temp)
    # print(ffvals_orig)
    
    swpt_orig,swpt_val_orig=mh_cm.find_sw_pt(ffvals_orig, kinplot)
    # print(swpt_orig,swpt_val_orig)
    
    swpt_kc,swpt_kc_val=mh_cm.find_sw_pt(ffvals_kc_temp, kinplot,swpt_val_orig)
    
    swptdiff_c=100*(swpt_kc-swpt_orig)/swpt_orig
    # print(swpt_diffs_c_a)
    if np.isnan(swptdiff_c)==False:
        swpt_diffs_c.append(swptdiff_c)
    # swpt_diffs_c_a[i]=swptdiff_c
    print('\n #######')
    print(swpt_orig,swpt_val_orig)
    print(swpt_kc,swpt_kc_val)
    print(swpt_kc,swpt_orig)
    print('####### \n')
    # print(swpt,swpt_kc) 
    # print('####### \n')
    # if ind_temp==indplot_switchretained:
    #     print('\n #######')
    #     print('\n #######')
    #     print(swpt_kc,swpt_orig)
    #     print(swptdiff_c)
    #     print(swpt_kc, swpt_kc_val)
    #     ax3.scatter(swpt_orig,swpt_val_orig,60,'k')
    #     ax3.scatter(swpt_kc,swpt_kc_val,60,'k')
    #     ax3.plot(kinplot,ffvals_orig,'-k')
    #     ax3.plot(kinplot,ffvals_kc_temp,'--k')
    #     print('####### \n')
    #     print('####### \n')
# print(swpt_diffs_c)
# print(swpt_diffs_c_a)
# print(np.max(swpt_diffs_c_a))
print(swpt_diffs_c)
mean_c=np.mean(swpt_diffs_c)
std_c=np.std(swpt_diffs_c)
err_c=std_c/np.sqrt(len(swpt_diffs_c))

print(mean_c)
print(std_c)
print(err_c)




filename='nonfixed_apr_final_KBmod_seed'+seedstr+'_n_'+str(n)
savename=savepath+filename
print(savepath)
print(filename)
p=Path(savepath)

filewantkcmod=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
print(filewantkcmod)

dat=np.load(filewantkcmod)
ffvals_kb=dat['fluxfrac']

swpt_diffs_b=list()
print(swpt_diffs_b)

# swpt_diffs_c_a=np.zeros(len(indexes_switch_retained))
# print(swpt_diffs_c_a)
for i in range(len(indexes_switch_retained)):
    ind_temp=indexes_switch_retained[i]

# for i in range(n):
#     ind_temp=i

    # swpt_orig=swpt[ind_temp]
    # swpt_val_orig=swptvals[ind_temp]
    ffvals_kb_temp=1-ffvals_kb[ind_temp]
    ffvals_orig=1-fluxfrac[ind_temp]
    # print(ffvals_kc_temp)
    # print(ffvals_orig)
    
    swpt_orig,swpt_val_orig=mh_cm.find_sw_pt(ffvals_orig, kinplot)
    # print(swpt_orig,swpt_val_orig)
    
    swpt_kb,swpt_kb_val=mh_cm.find_sw_pt(ffvals_kb_temp, kinplot,swpt_val_orig)
    
    swptdiff_b=100*(swpt_kb-swpt_orig)/swpt_orig
    # print(swpt_diffs_c_a)
    if np.isnan(swptdiff_b)==False:
        swpt_diffs_b.append(swptdiff_b)
    # swpt_diffs_c_a[i]=swptdiff_c
    # print('\n #######')
    # print(swpt_orig,swpt_val_orig)
    # print(swpt_kb,swpt_kb_val)
    # print(swpt_kb,swpt_orig)
    # print('####### \n')
    
    if swptdiff_b>0:
        print(swptdiff_b)
    # print(swpt,swpt_kc) 
    # print('####### \n')
    # if ind_temp==indplot_switchretained:
    #     print('\n #######')
    #     print('\n #######')
    #     print(swpt_kc,swpt_orig)
    #     print(swptdiff_c)
    #     print(swpt_kc, swpt_kc_val)
    #     ax3.scatter(swpt_orig,swpt_val_orig,60,'k')
    #     ax3.scatter(swpt_kc,swpt_kc_val,60,'k')
    #     ax3.plot(kinplot,ffvals_orig,'-k')
    #     ax3.plot(kinplot,ffvals_kc_temp,'--k')
    #     print('####### \n')
    #     print('####### \n')
# print(swpt_diffs_c)
# print(swpt_diffs_c_a)
# print(np.max(swpt_diffs_c_a))
# print(swpt_diffs_b)
mean_b=np.mean(swpt_diffs_b)
std_b=np.std(swpt_diffs_b)
err_b=std_b/np.sqrt(len(swpt_diffs_b))

print(mean_b)
print(std_b)
print(err_b)

print(mean_c)
print(std_c)
print(err_c)
# print(swpt_kcmod[indexes_switch_retained])
# print(swpt[indexes_switch_retained])

# swpt_diff_c=100*(swpt_kcmod[indexes_switch_retained]-swpt[indexes_switch_retained])/swpt[indexes_switch_retained]
# print(swpt_diff_c)

# mean_c=np.mean(swpt_diff_c)
# std_c=np.std(swpt_diff_c)
# err_c=std_c/np.sqrt(len(swpt_diff_c))


# filename='nonfixed_apr_final_KBmod_seed'+seedstr+'_n_'+str(n)
# savename=savepath+filename
# print(savepath)
# print(filename)
# p=Path(savepath)

# filewantkcmod=max([fn for fn in p.glob(filename+'*')], key=lambda f: f.stat().st_mtime)
# print(filewantkcmod)

# dat=np.load(filewantkcmod)
# swpt_kbmod=dat['swpt']
# swpt_kb_orig=dat['swpt']
# print(swpt_kbmod[indexes_switch_retained])
# print(swpt[indexes_switch_retained])

# swpt_diff_b=100*(swpt_kbmod[indexes_switch_retained]-swpt[indexes_switch_retained])/swpt[indexes_switch_retained]
# print(swpt_diff_b)

# mean_b=np.mean(swpt_diff_b)
# std_b=np.std(swpt_diff_b)
# err_b=std_b/np.sqrt(len(swpt_diff_b))

# print(mean_c)
# print(std_c)
# print(err_c)

# print(mean_b)
# print(std_b)
# print(err_b)

# if mean_b==0:
#     mean_b=0.5
    
means = np.array([mean_b,mean_c])

errors = np.array([err_b,err_c])

compartments = ('$\\rm Cytosol$', '$\\rm Mitchondria$')
y_pos = np.arange(len(compartments))
means = np.array([mean_b,mean_c])
errors = np.array([err_b,err_c])
clr = ('darkorange','teal')

h = ax6.barh(y_pos, means, xerr=errors, align='center',  label=compartments, color=clr)
ax6.set_yticks(y_pos, compartments)
ax6.set_xticks(np.linspace (-20,30,6),labels=['-$20$','-$10$','$0$','$10$','$20$','$30$'])
ax6.set_xlim(-20,35)
ax6.get_yaxis().set_visible(False) 
ax6.invert_yaxis() 
# ax6.legend(h,compartments,loc=(1.05,0.31),fontsize=15,title='$\\rm \\underline{Model}$')
ax6.set_xlabel('$\\rm \% ~change ~in ~switch ~point$')
ax6.text(10,0.2,'$\\rm \\underline{Model}$',fontsize=16)


means_th = np.array([-6.90,10.34])
errors_th = np.abs(np.array([5.87,-10.35]))


compartments = ('$\\rm Cytosol$', '$\\rm Mitchondria$')
y_pos = np.arange(len(compartments))
means_th = np.array([-6.90,10.34])
errors_th = np.abs(np.array([-10.35,5.87]))
# clr = ('y','c')

h = ax7.barh(y_pos, means_th, xerr=errors_th, align='center',  label=compartments, color=clr)
ax7.set_yticks(y_pos, compartments)
ax7.set_xticks(np.linspace (-20,30,6),labels=['','','','','',''])
ax7.set_xlim(-20,35)
ax7.get_yaxis().set_visible(False) 
ax7.invert_yaxis() 
ax7.text(5,0.3,'$\\rm \\underline{Experiment}$',fontsize=16)



# ax7.set_title('How fast do you want to go today?')

ax7.legend(h,compartments,loc=(1.05,-0.3),fontsize=15,title='$\\rm \\underline{Expressing~NADH}$ \n $\\rm \\underline{Oxidase~ in:}$')

# swpt_b=100*(swpt_kb_orig-swpt_orig)/swpt_orig
# print(swpt_orig)
# print(swpt_kb_orig)
# print(swpt_b)

# print(indexes_switch_retained)

# bsw_indexes=list()
# for i in range(n):#indexes_switch_retained:
#     if np.abs(swpt_b[i])<0:
#         print(swpt_b[i],i)
#         bsw_indexes.append(i)
        
# print(bsw_indexes)

# ind=indplot_switchretained
# print(ind)
# print(swpt[ind])
ax8.plot([0,0],[0.55,0.95],linewidth=3,color='darkorange')
ax8.plot([0,0],[0,1],'--k')
ax9.plot([0,0],[0,1],'--k')

ax8.set_xlim(-20,35)
ax9.set_xlim(-20,35)
ax8.set_ylim(0,1)
ax9.set_ylim(0,1)
        
if savefig==1:
    # plt.savefig('Figure4_supp.png',bbox_inches='tight',format='png')
    # plt.savefig('Figure4_supp.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_4.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_4.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_4.pdf',bbox_inches='tight',format='pdf')
    plt.savefig('Figure_4.svg',bbox_inches='tight',format='svg')

plt.show( )
