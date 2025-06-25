#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 27 14:34:15 2025

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
# from pathlib import Path

# import itertools

# sys.path.append('../data/megadata')
# sys.path.append('../data/heatmapdata')
sys.path.append('../simulation_functions')
sys.path.append('../data')


import make_heatmaps_compartment as mh_cm
import ode_funcs_compartment as ode_cm
# import matplotlib.image as img

import matplotlib.pyplot as plt
# import numpy.random as rn
# import datetime as dt

import matplotlib.pylab as pylab
import matplotlib as mpl

savefig=0           # change to 1 to save generated figure
numkin=81
# numkin=21

params = {'legend.fontsize': 12,
          'axes.labelsize': 15,
          'axes.titlesize': 15,
          'xtick.labelsize':15,
          'ytick.labelsize':15}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)


indplot_switchretained=697 #570

plt.figure()

ax3=plt.axes([2.1,-0.6,0.9,1.1])    # fluxfraction
ax4=plt.axes([3.2,-0.6,0.9,1.1])    # expt
ax5=plt.axes([2.1,-2,0.9,1.1])    # fluxfraction
ax6=plt.axes([3.2,-2,0.9,1.1])    # expt

ax3.text(10**-8.75,0.85,'$\\rm A$',fontsize=40)
ax4.text(10**-8.75,0.85,'$\\rm B$',fontsize=40)
ax5.text(10**-8.75,0.85,'$\\rm C$',fontsize=40)
ax6.text(10**-8.75,0.85,'$\\rm D$',fontsize=40)

ax3.text(10**-8, 0.7, '$K_{\\rm eq}^c \\downarrow$',fontsize=50)
ax4.text(10**-8, 0.7, '$K_{\\rm eq}^c \\uparrow$',fontsize=50)
ax5.text(10**-8, 0.7, '$K_{\\rm eq}^m \\downarrow$',fontsize=50)
ax6.text(10**-8, 0.7, '$K_{\\rm eq}^m \\uparrow$',fontsize=50)

allo=[1.9,1.,1.,1.]
n=1000

seed=1000
np.random.seed(seed)
if seed!=False:
    seedstr=str(seed)
else:
    seedstr='none'


args_to_run=mh_cm.make_params_newtest_6_nonfixed(allo, n)

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

ax3.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
ax3.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
ax3.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
ax3.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')

ax3.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax3.set_ylabel('$ \\rm Flux ~Fraction$',fontsize=30)
# ax3.legend(loc='upper left',ncols=1)
ax3.set_xscale('log')

ax3.set_xticks(np.logspace(-8,-5,4))#, ['$10^{-9}$','','$10^{-7}$','',
                                    # '$10^{-5}$','','$10^{-3}$',''])
                                    

ax4.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
ax4.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
ax4.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
ax4.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')

ax4.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax4.set_ylabel('$ \\rm Flux ~Fraction$',fontsize=30)
# ax4.legend(loc='upper left',ncols=1)
ax4.set_xscale('log')

ax4.set_xticks(np.logspace(-8,-5,4))#, ['$10^{-9}$','','$10^{-7}$','',
                                    # '$10^{-5}$','','$10^{-3}$',''])


ax5.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
ax5.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
ax5.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
ax5.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')

ax5.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax5.set_ylabel('$ \\rm Flux ~Fraction$',fontsize=30)
# ax5.legend(loc='upper left',ncols=1)
ax5.set_xscale('log')

ax5.set_xticks(np.logspace(-8,-5,4))#, ['$10^{-9}$','','$10^{-7}$','',
                                    # '$10^{-5}$','','$10^{-3}$',''])


ax6.plot(kinplot,1-(result[:,9]+result[:,11])/(result[:,9]+result[:,10]+result[:,11]+result[:,12]),label='Total NADH Fraction',linewidth=2,color='firebrick')
ax6.plot(kinplot,1-(result[:,9])/(result[:,9]+result[:,10]),label='Cytosolic NADH Fraction',linewidth=2,color='peru')
ax6.plot(kinplot,1-(result[:,11])/(result[:,11]+result[:,12]),label='Mitochondrial NADH Fraction',linewidth=2,color='mediumseagreen')
ax6.plot(kinplot,1-flux24/(flux23+flux24),label='Flux Fraction',linewidth=2,color='royalblue')

ax6.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax6.set_ylabel('$ \\rm Flux ~Fraction$',fontsize=30)
# ax6.legend(loc='upper left',ncols=1)
ax6.set_xscale('log')

ax6.set_xticks(np.logspace(-8,-5,4))#, ['$10^{-9}$','','$10^{-7}$','',
                                    # '$10^{-5}$','','$10^{-3}$',''])

print('###########\n orig')
print(args_temp[13])
print(args_temp[14])

params_B=np.array(args_temp[13])
params_C=np.array(args_temp[14])
print(params_B)
print(params_C)

params_B_inhib=params_B.copy()
params_B_prom=params_B.copy()
params_C_inhib=params_C.copy()
params_C_prom=params_C.copy()
print('##########\n orig')
print(params_B_inhib)
print(params_B_prom)
print(params_C_inhib)
print(params_C_prom)

params_B_inhib[3]=0.01
params_B_prom[3]=1.
params_C_inhib[3]=0.001
params_C_prom[3]=0.1
print('##########\n changed')
print(params_B_inhib)
print(params_B_prom)
print(params_C_inhib)
print(params_C_prom)
print(args_temp[13])
print(args_temp[14])


result3=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],params_B_inhib,args_temp[14],args_temp[15],
                             args_temp[16],args_temp[17])

  
flux23_3= ode_cm.flux_enz_forward(result3[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])             # pyr -> acet

flux24_3= ode_cm.flux_enz_forward(result3[:,5],  result3[:,11], 
                             args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


ax3.plot(kinplot,1-(result3[:,9]+result3[:,11])/(result3[:,9]+result3[:,10]+result3[:,11]+result3[:,12]),linestyle='--',linewidth=2,color='firebrick')
ax3.plot(kinplot,1-(result3[:,9])/(result3[:,9]+result3[:,10]),linestyle='--',linewidth=2,color='peru')
ax3.plot(kinplot,1-(result3[:,11])/(result3[:,11]+result3[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
ax3.plot(kinplot,1-flux24_3/(flux23_3+flux24_3),linestyle='--',linewidth=2,color='royalblue')


result4=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],params_B_prom,args_temp[14],args_temp[15],
                             args_temp[16],args_temp[17])

  
flux23_4= ode_cm.flux_enz_forward(result4[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])             # pyr -> acet

flux24_4= ode_cm.flux_enz_forward(result4[:,5],  result4[:,11], 
                             args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


ax4.plot(kinplot,1-(result4[:,9]+result4[:,11])/(result4[:,9]+result4[:,10]+result4[:,11]+result4[:,12]),linestyle='--',linewidth=2,color='firebrick')
ax4.plot(kinplot,1-(result4[:,9])/(result4[:,9]+result4[:,10]),linestyle='--',linewidth=2,color='peru')
ax4.plot(kinplot,1-(result4[:,11])/(result4[:,11]+result4[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
ax4.plot(kinplot,1-flux24_4/(flux23_4+flux24_4),linestyle='--',linewidth=2,color='royalblue')



result5=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],args_temp[13],params_C_inhib,args_temp[15],
                             args_temp[16],args_temp[17])

  
flux23_5= ode_cm.flux_enz_forward(result5[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])             # pyr -> acet

flux24_5= ode_cm.flux_enz_forward(result5[:,5],  result5[:,11], 
                             args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


ax5.plot(kinplot,1-(result5[:,9]+result5[:,11])/(result5[:,9]+result5[:,10]+result5[:,11]+result5[:,12]),linestyle='--',linewidth=2,color='firebrick')
ax5.plot(kinplot,1-(result5[:,9])/(result5[:,9]+result5[:,10]),linestyle='--',linewidth=2,color='peru')
ax5.plot(kinplot,1-(result5[:,11])/(result5[:,11]+result5[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
ax5.plot(kinplot,1-flux24_5/(flux23_5+flux24_5),linestyle='--',linewidth=2,color='royalblue')



result6=mh_cm.get_flux_vecs_2(args_temp[0],kinplot,args_temp[2],args_temp[3],
                             args_temp[4],args_temp[5],args_temp[6],args_temp[7],
                             args_temp[8],args_temp[9],args_temp[10],args_temp[11],
                             args_temp[12],args_temp[13] ,params_C_prom ,args_temp[15],
                             args_temp[16],args_temp[17])

  
flux23_6= ode_cm.flux_enz_forward(result6[:,2]**allo[0], 1.,
                             args_temp[6][0],  args_temp[6][1])             # pyr -> acet

flux24_6= ode_cm.flux_enz_forward(result6[:,5],  result6[:,11], 
                             args_temp[8][0], args_temp[8][1])            # pyr_mito -> AcCoA


ax6.plot(kinplot,1-(result6[:,9]+result6[:,11])/(result6[:,9]+result6[:,10]+result6[:,11]+result6[:,12]),linestyle='--',linewidth=2,color='firebrick')
ax6.plot(kinplot,1-(result6[:,9])/(result6[:,9]+result6[:,10]),linestyle='--',linewidth=2,color='peru')
ax6.plot(kinplot,1-(result6[:,11])/(result6[:,11]+result6[:,12]),linestyle='--',linewidth=2,color='mediumseagreen')
ax6.plot(kinplot,1-flux24_6/(flux23_6+flux24_6),linestyle='--',linewidth=2,color='royalblue')

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



if savefig==1:
    # plt.savefig('Figure4_supp.png',bbox_inches='tight',format='png')
    # plt.savefig('Figure4_supp.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S8.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S8.eps',bbox_inches='tight',format='eps')
    # plt.savefig('Figure_compartments_main_v2.pdf',bbox_inches='tight',format='pdf')
    # plt.savefig('Figure_compartments_main_v2.svg',bbox_inches='tight',format='svg')

plt.show( )