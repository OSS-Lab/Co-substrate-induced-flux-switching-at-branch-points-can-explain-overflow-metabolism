#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 10:46:08 2024

@author: robert
"""



import sys


sys.path.append('../simulation_functions/')
sys.path.append('../data/')

import numpy as np

import matplotlib.pyplot as plt

import matplotlib.pylab as pylab
import matplotlib.image as img
import matplotlib as mpl



savefig=0           # change to 1 to save generated figure


cmapuse='RdBu_r'


params = {'legend.fontsize': 30,
          'figure.figsize': (12, 9),
          'axes.labelsize': 40,
          'axes.titlesize': 40,
          'xtick.labelsize':30,
          'ytick.labelsize':30}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

plt.figure()

ax1=plt.axes([0.1,0.05,0.6,0.5])
inset1=plt.axes([0.3,0.3,0.28,0.28]) 

motifnone=img.imread('fig3_symm.png')
plt.imshow(motifnone)
ax1.text(10**-4,10**2,'$\\rm A$', fontsize=40)
plt.axis('off')


#############################


ax2=plt.axes([0.85,0.05,0.6,0.5])
inset2=plt.axes([1.05,0.3,0.28,0.28]) 

motifnone=img.imread('fig3_symm_upstream.png')
plt.imshow(motifnone)
ax2.text(10**-4,10**2,'$\\rm B$', fontsize=40)
plt.axis('off')

#############################


ax3=plt.axes([0.1,-0.55,0.6,0.5])
inset3=plt.axes([0.3,-0.3,0.28,0.28]) 

motifnone=img.imread('fig3_minimal.png')
plt.imshow(motifnone)
ax3.text(10**-4,10**2,'$\\rm D$', fontsize=40)
plt.axis('off')

#############################


ax4=plt.axes([0.85,-0.55,0.6,0.5])
inset4=plt.axes([1.05,-0.3,0.28,0.28]) 

motifnone=img.imread('fig3_minimal_upstream.png')
plt.imshow(motifnone)
ax4.text(10**-4,10**2,'$\\rm E$', fontsize=40)
plt.axis('off')

#################################

ax5=plt.axes([1.6,-0.55,0.6,0.5])
inset5=plt.axes([1.8,-0.3,0.28,0.28]) 

motifnone=img.imread('fig3_beforebranch.png')
plt.imshow(motifnone)
ax5.text(10**-4,10**2,'$\\rm F$', fontsize=40)
plt.axis('off')

##################################

ax6=plt.axes([1.6,0.05,0.6,0.5])
inset6=plt.axes([1.8,0.3,0.28,0.28]) 

motifnone=img.imread('fig3_downstream.png')
plt.imshow(motifnone)
ax6.text(10**-4,10**2,'$\\rm C$', fontsize=40)
plt.axis('off')

xvar='Btot'
yvar='CKcatB'


KeqB_val=1.
CKcat_top_val=0.1
CKcat_bottom_val=1.

reversible=False


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

Tmax=100000
tpts=10000

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

lengthparams=[2,1,1]
B_position_params=([0,0],[0,0],[1],[1])

loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
            '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
            str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
            '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
            '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
            '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
            yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
            
loadstring_ff='ff_'

for idop, op in enumerate(othervars_list):
    loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'

loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      

ffdelta_mat=np.load(loadname_ff)

im=ax1.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(xlabstr)    
ax1.set_ylabel(ylabstr,rotation=0)   
ax1.xaxis.set_label_coords(0.85,-0.04)
ax1.yaxis.set_label_coords(-0.12,0.45)
ax1.set_xticks(np.logspace(-3,3,4))
ax1.set_yticks(np.logspace(-3,3,4))

cbar = plt.colorbar(im)
cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-20, y=0.9)


####################################

lengthparams=[2,1,1]
B_position_params=([1,0],[0,0],[1],[1])

loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
            '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
            str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
            '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
            '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
            '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
            yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
            
loadstring_ff='ff_'

for idop, op in enumerate(othervars_list):
    loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'

loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      

ffdelta_mat=np.load(loadname_ff)

im=ax2.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(xlabstr)    
ax2.set_ylabel(ylabstr,rotation=0)   
ax2.xaxis.set_label_coords(0.85,-0.04)
ax2.yaxis.set_label_coords(-0.12,0.45)
ax2.set_xticks(np.logspace(-3,3,4))
ax2.set_yticks(np.logspace(-3,3,4))

cbar = plt.colorbar(im)
cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-20, y=0.9)

####################################

lengthparams=[2,1,1]
B_position_params=([0,0],[0,0],[0],[1])

loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
            '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
            str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
            '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
            '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
            '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
            yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
            
loadstring_ff='ff_'

for idop, op in enumerate(othervars_list):
    loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'

loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      

ffdelta_mat=np.load(loadname_ff)

im=ax3.pcolormesh(np.logspace(xmin,xmax,xnum),
                    np.logspace(ymin,ymax,ynum),
                    -ffdelta_mat,edgecolors='face',
                    vmin=-1,vmax=1,cmap=cmapuse)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel(xlabstr)    
ax3.set_ylabel(ylabstr,rotation=0) 
ax3.xaxis.set_label_coords(0.85,-0.04)
ax3.yaxis.set_label_coords(-0.12,0.45)
ax3.set_xticks(np.logspace(-3,3,4))
ax3.set_yticks(np.logspace(-3,3,4))

cbar = plt.colorbar(im)
cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-20, y=0.9)


####################################

lengthparams=[2,1,1]
B_position_params=([1,0],[0,0],[0],[1])

loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
            '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
            str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
            '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
            '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
            '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
            yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
            
loadstring_ff='ff_'

for idop, op in enumerate(othervars_list):
    loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'

loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      

ffdelta_mat=np.load(loadname_ff)

im=ax4.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlabel(xlabstr)    
ax4.set_ylabel(ylabstr,rotation=0)  
ax4.xaxis.set_label_coords(0.85,-0.04)
ax4.yaxis.set_label_coords(-0.12,0.45)
ax4.set_xticks(np.logspace(-3,3,4))
ax4.set_yticks(np.logspace(-3,3,4))

cbar = plt.colorbar(im)
cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-20, y=0.9)


####################################

lengthparams=[2,1,1]
B_position_params=([0,1],[0,0],[0],[0])

loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
            '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
            str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
            '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
            '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
            '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
            yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
            
loadstring_ff='ff_'

for idop, op in enumerate(othervars_list):
    loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'

loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      

ffdelta_mat=np.load(loadname_ff)

im=ax5.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax5.set_xscale('log')
ax5.set_yscale('log')
ax5.set_xlabel(xlabstr)    
ax5.set_ylabel(ylabstr,rotation=0)  
ax5.xaxis.set_label_coords(0.85,-0.04)
ax5.yaxis.set_label_coords(-0.12,0.45)
ax5.set_xticks(np.logspace(-3,3,4))
ax5.set_yticks(np.logspace(-3,3,4))

cbar = plt.colorbar(im)
cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-20, y=0.9)

####################################

lengthparams=[2,1,2]
B_position_params=([0,0],[0,0],[0],[0,1])

loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
            '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
            str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
            '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
            '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
            '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
            yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)
            
loadstring_ff='ff_'

for idop, op in enumerate(othervars_list):
    loadstring_ff += str(op) + '_' + str(othervars_tuple[idop]) + '__'

loadname_ff=loadpath+'/'+loadstring_ff+'.npy'      

ffdelta_mat=np.load(loadname_ff)

im=ax6.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax6.set_xscale('log')
ax6.set_yscale('log')
ax6.set_xlabel(xlabstr)    
ax6.set_ylabel(ylabstr,rotation=0)  
ax6.xaxis.set_label_coords(0.85,-0.04)
ax6.yaxis.set_label_coords(-0.12,0.45)
ax6.set_xticks(np.logspace(-3,3,4))
ax6.set_yticks(np.logspace(-3,3,4))

cbar = plt.colorbar(im)
cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-20, y=0.9)

if savefig==1:
    plt.savefig('Figure_3.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_3.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_3.svg',bbox_inches='tight',format='svg')

plt.show()
