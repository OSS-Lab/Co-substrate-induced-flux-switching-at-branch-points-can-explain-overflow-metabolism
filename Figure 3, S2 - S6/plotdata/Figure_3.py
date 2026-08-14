#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 15:02:59 2025

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

from matplotlib.patches import Rectangle

savefig=0         # change to 1 to save generated figure


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

ax1=plt.axes([0.1,0.05,0.48,0.5])
inset1=plt.axes([0.3,0.322,0.28,0.28]) 

motifnone=img.imread('fig3_minimal.png')
plt.imshow(motifnone)
# ax1.text(10**-4,10**2,'$\\rm A$', fontsize=40)
plt.axis('off')


#############################


ax2=plt.axes([0.61,0.05,0.48,0.5])
inset2=plt.axes([0.81,0.322,0.28,0.28]) 

motifnone=img.imread('fig3_minimal_upstream.png')
plt.imshow(motifnone)
# ax2.text(10**-4,10**2,'$\\rm B$', fontsize=40)
plt.axis('off')

#############################


ax3=plt.axes([1.12,-0.48,0.48,0.5])
inset3=plt.axes([1.32,-0.208,0.28,0.28]) 

motifnone=img.imread('fig3_beforebranch.png')
plt.imshow(motifnone)
# ax3.text(10**-4,10**2,'$\\rm D$', fontsize=40)
plt.axis('off')

# ax6=plt.axes([1.6,0.05,0.6,0.5])
inset6=plt.axes([1.32,-0.38,0.28,0.28]) 

motifnone=img.imread('fig3_downstream.png')
plt.imshow(motifnone)
# ax6.text(10**-4,10**2,'$\\rm C$', fontsize=40)
plt.axis('off')
#############################


ax4=plt.axes([0.1,-0.48,0.48,0.5])
inset4=plt.axes([0.3,-0.208,0.28,0.28]) 

motifnone=img.imread('fig3_symm.png')
plt.imshow(motifnone)
# ax4.text(10**-4,10**2,'$\\rm E$', fontsize=40)
plt.axis('off')

#################################

ax5=plt.axes([0.61,-0.48,0.48,0.5])
inset5=plt.axes([0.81,-0.208,0.28,0.28]) 

motifnone=img.imread('fig3_symm_upstream.png')
plt.imshow(motifnone)
# ax5.text(10**-4,10**2,'$\\rm F$', fontsize=40)
plt.axis('off')

##################################

ax6=plt.axes([1.12,0.05,0.48,0.5])
inset7=plt.axes([1.32,0.322,0.28,0.28]) 

motifnone=img.imread('fig3_branchdownstream.png')
plt.imshow(motifnone)
# motifnone=img.imread('fig3_symm_upstream.png')
# plt.imshow(motifnone)
# ax5.text(10**-4,10**2,'$\\rm F$', fontsize=40)
plt.axis('off')

##################################

ax1.text(10**-3,10**-3,'$\\rm A$', fontsize=40)
ax2.text(10**-3,10**-3,'$\\rm B$', fontsize=40)
ax3.text(10**-3,10**-3,'$\\rm F$', fontsize=40)
ax4.text(10**-3,10**-3,'$\\rm D$', fontsize=40)
ax5.text(10**-3,10**-3,'$\\rm E$', fontsize=40,color='w')
ax6.text(10**-3,10**-3,'$\\rm C$', fontsize=40)






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
num_kin=101

Tmax=100000
tpts=10000

titlefontsize=20

####### these dont change in the irrev case
Keq_top=1.
Keq_bottom=1.

allvars=['Btot','CKcatB','KeqB','CKcat_bottom','CKcat_top','Keq_top','Keq_bottom']


xlabstr='${\\rm Total ~co}$-${\\rm substrate ~pool ~size,~} {B_{\\rm T}}$'
ylabstr='${\\rm Maximal~ rate~ of ~co}$-${\\rm substrate ~conversion,~} V_{\\rm max}^{B_0 \\leftrightarrow B_1}$'
# ylabstr='${\\rm V_{M}^B}$'


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

im=ax1.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(ylabstr,rotation=90,fontsize=35)   
ax1.yaxis.set_label_coords(-0.15,-0.05)
ax1.set_xticks(np.logspace(-3,3,4))
ax1.set_xticklabels('')

ax1.set_yticks(np.logspace(-3,3,4))


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

im=ax2.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xticks(np.logspace(-3,3,4))

ax2.set_xticklabels('')
ax2.set_yticklabels('')


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

im=ax3.pcolormesh(np.logspace(xmin,xmax,xnum),
                    np.logspace(ymin,ymax,ynum),
                    -ffdelta_mat,edgecolors='face',
                    vmin=-1,vmax=1,cmap=cmapuse)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xticks(np.logspace(-3,3,4))
ax3.set_yticks(np.logspace(-3,3,4),['','','',''])

cax = ax3.inset_axes([1.09, 0, 0.04, 1.85])
cbar=plt.colorbar(im, cax=cax, orientation='vertical',ticks=[-1, -0.5,0,0.5, 1])
cbar.ax.tick_params(labelsize=25)

# cbar = plt.colorbar(im,location='top',shrink=0.8)#pad=0.1)
# cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-40, y=2.9)
# ax=cbar.ax
ax3.text(10**3.3,10**9,'$\\Delta F_f$',fontsize=40,rotation=0)

####################################

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

im=ax4.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xticks(np.logspace(-3,3,4))
ax4.set_yticks(np.logspace(-3,3,4))
# ax4.set_yticklabels('')




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

im=ax5.pcolormesh(np.logspace(xmin,xmax,xnum),
                   np.logspace(ymin,ymax,ynum),
                   -ffdelta_mat,edgecolors='face',
                   vmin=-1,vmax=1,cmap=cmapuse)
ax5.set_xscale('log')
ax5.set_yscale('log')

ax5.set_xticks(np.logspace(-3,3,4))
ax5.set_yticks(np.logspace(-3,3,4))
ax5.set_yticklabels('')



ax5.set_xlabel(xlabstr,fontsize=35)    
ax5.xaxis.set_label_coords(0.5,-0.13)
####################################

lengthparams=[2,2,2]
B_position_params=([0,0],[0,0],[0,1],[1,0])

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

ax6.set_xticks(np.logspace(-3,3,4),['','','',''])
ax6.set_yticks(np.logspace(-3,3,4))
ax6.set_yticklabels('')

# ax6=plt.axes([1.12,0.15,0.48,0.4])
# # ax6.patch.set_facecolor('lightgrey')
# plt.axis('off')


# ax6.add_patch(Rectangle((0, 0), 1, 1,
#              edgecolor = 'lightgrey',
#              facecolor = 'lightgrey',
#              fill=True,
#              lw=0))
# ax6.text(0.5,0.85,'\\bf \\underline{Regulation by} \n \\bf \\underline{pathway structure}',fontsize=30,
#          horizontalalignment='center',
#          verticalalignment='center')

# ax6.text(0.05,0.65,'$\\bullet$ Driven by the asymmetric usage of a \n co-substrate immediately downstream \n of a branchpoint.',
#          horizontalalignment='left',
#          verticalalignment='top',
#          fontsize=25)

# ax6.text(0.05,0.3,'$\\bullet$ Asymmetry modulated by upstream \n co-substrate usages.',
#          horizontalalignment='left',
#          verticalalignment='top',
#          fontsize=25)

if savefig==1:
    plt.savefig('Figure_3.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_3.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_3.svg',bbox_inches='tight',format='svg')
    plt.savefig('Figure_3.pdf',bbox_inches='tight',format='pdf')

plt.show()
