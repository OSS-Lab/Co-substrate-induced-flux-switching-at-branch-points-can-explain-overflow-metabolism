#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:35:33 2025

@author: robert
"""
import sys

sys.path.append('../simulation_functions')
sys.path.append('../data')
import numpy as np

import matplotlib.pyplot as plt

# import make_odes as mp
# import time_series_plots as ts
import matplotlib.pylab as pylab
import make_figure_data as mfd
import matplotlib.image as img
import matplotlib as mpl




savefig=0
numkin=21
# numkin=51 

kinvals=np.logspace(-4,1,numkin)


xvar='Btot'
yvar='CKcatB'

# KeqB_val=1.

KeqB_val=0.1
CKcat_top_val=1.
CKcat_bottom_val=0.1

reversible=False



##### THIS IS FIGURE S1

cmapuse='RdBu_r'

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 9),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':20,
         'ytick.labelsize':20}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

plt.figure()
ax1=plt.axes([0.13,0.48,0.5,0.5])
plt.axis('off')
motifnone=img.imread('figS4_mot_nochange.png')
plt.imshow(motifnone)
ax2=plt.axes([0.7,0.6,0.3,0.3])
ax2.text(10**-19,10**2.1,'$\\rm A$', fontsize=35)
ax2.text(10**-4.5,10**1.7,'$\\rm B$', fontsize=35)

ax5=plt.axes([1.15,0.6,0.3,0.3])

# ax3=plt.axes([0.1,0.12,0.5,0.5])
# plt.axis('off')
# motifnone=img.imread('motif_branch_downstream_tr.png')
# plt.imshow(motifnone)
# ax4=plt.axes([0.7,0.22,0.3,0.3])
# ax4.text(10**-18.5,10**1.7,'$\\rm C$', fontsize=35)
# ax4.text(10**-4.5,10**1.7,'$\\rm D$', fontsize=35)

# ax5=plt.axes([0.1,-0.26,0.5,0.5])
# plt.axis('off')
# motifnone=img.imread('222noup_asym.png')
# plt.imshow(motifnone)
# ax6=plt.axes([0.7,-0.16,0.3,0.3])
# ax6.text(10**-18.5,10**1.7,'$\\rm E$', fontsize=35)
# ax6.text(10**-4.5,10**1.7,'$\\rm F$', fontsize=35)#,color='w')

# ax7=plt.axes([0.1,-0.64,0.5,0.5])
# plt.axis('off')
# motifnone=img.imread('222noup_1branch.png')
# plt.imshow(motifnone)
# ax8=plt.axes([0.7,-0.54,0.3,0.3])

# ax9=plt.axes([0.1,-1.02,0.5,0.5])
# plt.axis('off')
# motifnone=img.imread('asym_twoupstreamonbranch.png')
# plt.imshow(motifnone)
# ax10=plt.axes([0.7,-0.92,0.3,0.3])

# ax11=plt.axes([0.1,-1.4,0.5,0.5])
# plt.axis('off')
# motifnone=img.imread('asym_twoupstreamonbranch_gap.png')
# plt.imshow(motifnone)
# ax12=plt.axes([0.7,-1.3,0.3,0.3])

# ax13=plt.axes([0.1,-1.78,0.5,0.5])
# plt.axis('off')
# motifnone=img.imread('asym_gap_upstream.png')
# plt.imshow(motifnone)
# ax14=plt.axes([0.7,-1.68,0.3,0.3])

lengthparams=[2,2,2]
				  
B_position_params=([0,0],[0,0],[0,1],[0,1])


 
        
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

titlefontsize=20

####### these dont change in the irrev case
Keq_top=1.
Keq_bottom=1.

allvars=['Btot','CKcatB','KeqB','CKcat_bottom','CKcat_top','Keq_top','Keq_bottom']


xlabstr='${\\rm B_T}$'
ylabstr='${\\rm V_{M}^B}$'


xvarvals=np.logspace(-3,3,25)
yvarvals=np.logspace(-3,3,25)



xmin=np.log10(np.min(xvarvals))
xmax=np.log10(np.max(xvarvals))
xnum=len(xvarvals)

ymin=np.log10(np.min(yvarvals))
ymax=np.log10(np.max(yvarvals))
ynum=len(yvarvals)

othervars_list=sorted(allvars.copy())
print(othervars_list)
othervars_list.remove(xvar)
othervars_list.remove(yvar)
print(othervars_list)


othervars_tuple=()
for id_other, other in enumerate(othervars_list):
    print(other)
    if other=='KeqB':
        add_vals=(KeqB_val,)
    elif other=='CKcat_bottom':
        add_vals=(CKcat_bottom_val,)
    elif other=='CKcat_top':
        add_vals=(CKcat_top_val,)
    elif other=='Keq_bottom':
        add_vals=(1.,)
    elif other=='Keq_top':
        add_vals=(1.,)
        
    othervars_tuple+=add_vals


# for idstruc,struc in enumerate(lengthparams_set):
    

#     lengthparams=struc
#     B_position_params=B_position_params_set[idstruc]

#     print(lengthparams)
#     print(B_position_params)  
    
loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
        '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
        str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
        '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
        '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
        '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
        yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
            
    #     loadstring_ff='ff_'
        
    #     for idop, op in enumerate(othervars_list):
    #         loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'
        
    #     loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      
    
    
    #     ffdelta_mat=np.load(loadname_ff)
    #     print(loadpath)
        
    # else:
# loadpath='../../../../general branching model change cosubs_new_aug24/data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
#         '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
#         str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
#         '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
#         '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
#         '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
#         yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
        
loadstring_ff='ff_'

for idop, op in enumerate(othervars_list):
    loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'

loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      


ffdelta_mat=np.load(loadname_ff)
print(loadpath)
   
        
im=ax2.pcolormesh(np.logspace(xmin,xmax,xnum),
                        np.logspace(ymin,ymax,ynum),
                        -ffdelta_mat,edgecolors='face',
                        vmin=-1,vmax=1,cmap=cmapuse)
        #edgecolors{'none', None, 'face', color, color sequence}, optional
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(xlabstr,fontsize=25)    
ax2.set_ylabel(ylabstr,fontsize=25,rotation=0)   
ax2.xaxis.set_label_coords(0.8,-0.08)
ax2.yaxis.set_label_coords(-0.15,0.45)
ax2.set_xticks(np.logspace(-3,3,4))
ax2.set_yticks(np.logspace(-3,3,4))
print('checking max for both downstream, '+\
      str(np.max(ffdelta_mat[:]))+'\n')
print('checking min for both downstream, '+\
      str(np.min(ffdelta_mat[:]))+'\n')
 
  

cbar = plt.colorbar(im)
#cbar.set_label('$\\Delta F_f$',fontsize=30,rotation=0,labelpad=-40, y=1.15)
cbar.set_label('$\\Delta F_f$',fontsize=30,rotation=0,labelpad=0, y=0.9)

V24_vals=np.array([0.1,1,10])
    
cols=['blue','black','red']
inset1=plt.axes([1.335,0.787,0.1,0.1])

Btot=10.
VMB=10**-2.
KeqB=1.
rev=False
for j in range(len(V24_vals)):
    print(j)
    V24=V24_vals[j]
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    b_vec=np.zeros(np.size(kinvals))
    b_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,2,2],B_position_params=([0,0],[0,0],[0,1],[0,1]),
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB,Btot],
                          params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
            b_vec[i]=1.-bfrac
        
    ax5.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$')
    inset1.plot(kinvals,b_vec,color=cols[j])#,linestyle=':') # add as inset
    
ax5.set_xscale('log')
# ax5.set_yscale('log')
ax5.legend(loc=(0.6,0.05),fontsize='large')
ax5.set_ylim(0,1)
ax5.set_xticks(np.logspace(-4,0,5),['$10^{-4}$','','$10^{-2}$','','$10^0$'])
ax5.set_xlabel('$v_{\\rm in}$',fontsize=25)
ax5.xaxis.set_label_coords(0.7,-0.05)
ax5.set_ylabel('$\\rm Flux~ Fraction$',fontsize=25)
ax5.yaxis.set_label_coords(-0.2,0.4)
ax5.text(10**-5,0.8,'$\\rm C$', fontsize=35)


inset1.set_xscale('log')
# inset1.set_yscale('log')
inset1.set_ylim(0,1)
inset1.set_xticks(np.logspace(-4,0,3),['$10^{-4}$','','$10^0$'])
# inset1.set_xticklabels(['$10^{-4}$','$$','$10^0$'])
inset1.set_xlabel('$v_{\\rm in}$',fontsize=20)
inset1.xaxis.set_label_coords(0.7,0)
inset1.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=20)
inset1.yaxis.set_label_coords(-0.2,0.2)

if savefig==1:
    plt.savefig('Figure_S4.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S4.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S4.svg',bbox_inches='tight',format='svg')
    
plt.show()
print(np.min(ffdelta_mat[:]))