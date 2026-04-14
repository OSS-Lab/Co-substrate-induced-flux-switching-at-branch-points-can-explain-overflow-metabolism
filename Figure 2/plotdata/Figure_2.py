#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 12:00:48 2025

@author: robert
"""


import sys

sys.path.append('../simulation_functions')
import numpy as np

import matplotlib.pyplot as plt

import matplotlib.pylab as pylab
import make_figure_data as mfd
import matplotlib.image as img
import matplotlib as mpl
# import matplotlib.patheffects as pe
import matplotlib.patches as patches



   
    
savefig=0          # Change to 1 in order to save generated figure
numkin=501
# numkin=21
kinvals=np.logspace(-4,1,numkin)
TMAX=100000
TPTS=TMAX

# Ks_top=10**(2.)
V23=1.
K23=1.

Btot=10.
VMB=0.01
KeqB=1.

# V24_vals=np.array([1.,0.5, 2., 0.5,2.])
# K24_vals=np.array([1.,0.5, 2.,2. ,0.5])

# V24_vals=np.array([1.,0.1,10.])
# K24_vals=np.array([1.,0.1 ,0.1])
V24_vals=np.array([1.,0.5,2.])
K24_vals=np.array([1.,0.5 ,0.5])

# V24_vals=np.array([1.,10., 0.1, 10.,0.1])
# K24_vals=np.array([1.,10.,0.1,0.1,10.])

rev=False


al_sub23_vals=np.array([1.,1.3,1.,1.3])

lengthparams=[2,1,1]
B_position_params_sets=(([0,0],[0,0],[0],[0]),
                        ([0,0],[0,0],[0],[0]),
                        ([0,0],[0,0],[0],[1]),
                        ([0,0],[0,0],[0],[1]))

params = {'legend.fontsize': 20,
          'figure.figsize': (12, 9),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':25,
          'ytick.labelsize':25}

# cols=['black','blue','red','darkgoldenrod','darkmagenta']
cols=['black','cornflowerblue','orangered']
lw=3
col_labs=['A','B','C','D']
axes_labs=['$\\rm Branch~ Point$','$\\rm With~ Allostery$','$\\rm With~ Cosubstrates$','$\\rm With~ Both$']

pylab.rcParams.update(params)
mpl.rc('text', usetex = True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{color}')

plt.figure()
ax1=plt.axes([0.091,0.65,0.62,0.5])
motifnone=img.imread('fig2_none_tr.png')
plt.imshow(motifnone)
plt.axis('off')

ax2=plt.axes([0.73,0.65,0.62,0.5])
motifnone=img.imread('branch_with_allo_2.png')
plt.imshow(motifnone)
plt.axis('off')

ax3=plt.axes([1.368,0.65,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('fig2_minimal_tr.png')
plt.imshow(motifbranch)

ax4=plt.axes([2.01,0.65,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('fig2_allo_cosub.png')
plt.imshow(motifbranch)



ax5=plt.axes([0.1,0.08,0.6,0.6])
ax6=plt.axes([0.74,0.08,0.6,0.6])
ax7=plt.axes([1.38,0.08,0.6,0.6])
ax8=plt.axes([2.02,0.08,0.6,0.6])

inset1=plt.axes([1.76,0.18,0.2,0.2])
inset2=plt.axes([2.4,0.18,0.2,0.2])


ax9=plt.axes([0.1,-0.55,0.6,0.6])
ax10=plt.axes([0.74,-0.55,0.6,0.6])
ax11=plt.axes([1.38,-0.55,0.6,0.6])
ax12=plt.axes([2.02,-0.55,0.6,0.6])

axsets=[[ax5,ax9],[ax6,ax10],[ax7,ax11],[ax8,ax12]]
ins=[inset1,inset2]

for idmot,al23sub in enumerate(al_sub23_vals):
    B_position=B_position_params_sets[idmot]    
    axes_plot=axsets[idmot]
    print(B_position,al23sub)
    
    al=np.array([al23sub,1.,1.,1.])
    for idv24, V24 in enumerate(V24_vals):
        K24=K24_vals[idv24]
        print(str(V23/V24),str(K23/K24))
        
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        
        b_vec=np.zeros(np.size(kinvals))
        b_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],
                              B_position_params=B_position,
                              params_global=[1.,1.,1.,1.],
                              params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
                              params_branch=[V23,K23,1.,1.,V24,K24,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=rev,Tmax=TMAX,tpts=TPTS)
            if sum(isbuildup)==0:
                ff_vec[i]=1.-fluxfrac
                b_vec[i]=1.-bfrac
            
    
        gradient=np.abs(np.gradient(ff_vec, kinvals))
        sense=gradient/(ff_vec/kinvals)
        
        axes_plot[0].plot(kinvals,ff_vec,color=cols[idv24],linewidth=lw,
                          label='$\\rm V_{23}/V_{24}='+str(V23/V24)+',~ K_{23}/K_{24}='+str(K23/K24)+'$')        
        axes_plot[1].plot(kinvals,sense,color=cols[idv24],linewidth=lw,
                          label='$\\rm V_{23}/V_{24}='+str(V23/V24)+',~ K_{23}/K_{24}='+str(K23/K24)+'$')    
        # axes_plot[1].plot(kinvals,gradient,color=cols[idv24],linewidth=lw,
        #                   label='$\\rm V_{23}/V_{24}='+str(V23/V24)+',~ K_{23}/K_{24}='+str(K23/K24)+'$')    
        print(np.nanmax(sense))
        if idmot==1 and idv24 ==1:
            ###### THIS IS THE DELTA F_f LABEL BIT
            min_y=np.nanmin(ff_vec)
            min_y_arg=np.nanargmin(ff_vec)
            print(min_y)
            print(min_y_arg)
            
            max_y=np.nanmax(ff_vec)
            max_y_arg=np.nanargmax(ff_vec)
            print(max_y)
            print(max_y_arg)
            ax6.plot([10**-4,10**0],[max_y,max_y],color=cols[idv24],linewidth=3,linestyle='--')

    # axes_plot[0].xaxis.set_label_coords(0.8,-0.05)
    if idmot==0:
        axes_plot[0].set_ylabel('${\\rm Fraction ~of ~upper ~branch ~flux}, F_f$',fontsize=30)
        axes_plot[0].yaxis.set_label_coords(-0.1,0.5)
        axes_plot[0].set_yticks(np.linspace(0,1,6))
        
        axes_plot[1].set_ylabel('${\\rm Sensitivity}, \\left\\lvert \\frac{d F_f}{dv_{\\rm in}} \\right\\rvert \\big/ \\frac{F_f}{v_{\\rm in}}$',fontsize=30)
        axes_plot[1].yaxis.set_label_coords(-0.08,0.5)
        axes_plot[1].set_yticks(np.linspace(0,8,5))

    else:
        # 
        axes_plot[0].set_yticks(np.linspace(0,1,6),labels=['','','','','',''])
        axes_plot[1].set_yticks(np.linspace(0,8,5),labels=['','','','',''])
        
        
    if idmot>1:
        ins_ax=ins[idmot-2]
        
        ins_ax.plot(kinvals,b_vec,color=cols[idv24],linewidth=lw)
        ins_ax.set_xscale('log')
        ins_ax.set_ylim(-0.05,1.05)
        ins_ax.set_xticks(np.logspace(-4,0,3))
        ins_ax.set_xlabel('$v_{\\rm in}$',fontsize=35)
        ins_ax.xaxis.set_label_coords(0.75,-0.15)
        ins_ax.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
        ins_ax.yaxis.set_label_coords(-0.15,0.3)
        ins_ax.set_yticks(np.linspace(0,1,3),labels=['$0$','','$1$'])

    axes_plot[0].text(10**-3.3,1.65,axes_labs[idmot],fontsize=30,bbox=dict(facecolor='none', edgecolor='black',boxstyle='round'))
    axes_plot[0].set_xscale('log')
    axes_plot[0].set_ylim(0,1.05)
    axes_plot[0].set_xlim(10**-4.2,10**0.2)
    axes_plot[0].text(10**-4.1,0.94,'$\\rm ' +str(col_labs[idmot])+'(ii)$', fontsize=40)
    axes_plot[0].text(10**-4.1,1.65,'$\\rm ' +str(col_labs[idmot])+'(i)$', fontsize=40)
    axes_plot[0].set_xticks(np.logspace(-4,0,5),labels=['','','','',''])
    axes_plot[0].text(10**-4.15, 1.47, '$v_{\\rm in}$',fontsize=40)
    axes_plot[0].set_xticks(np.logspace(-4,0,5),labels=['','','','',''])
    
    axes_plot[1].set_xscale('log')
    axes_plot[1].text(10**-4.1,8,'$\\rm ' +str(col_labs[idmot])+'(iii)$', fontsize=40)
    axes_plot[1].set_ylim(0,9)
    axes_plot[1].set_xlabel('$v_{\\rm in}$',fontsize=40)
    axes_plot[1].xaxis.set_label_coords(0.65,-0.05)
    axes_plot[1].set_xticks(np.logspace(-4,0,5))
    axes_plot[1].set_xlim(10**-4.2,10**0.2)
    
    
    
    
    
ax5.legend(ncols=1)

ax_labs=plt.axes([0.7,0.08,0.6,0.6])
ax_labs.axis('off')

ax_labs.set_ylim(0,1.05)
ax_labs.set_xlim(0,1)


start1=0.4
start2=0.3

ax_labs.arrow(0.0675, start1, 0, max_y-start1,length_includes_head=True,
          width=0.005,head_width=0.03,head_length=0.05,
          color=cols[1])#, **kwargs)

ax_labs.arrow(0.0675, start2, 0, min_y-start2,length_includes_head=True,
          width=0.005,head_width=0.03,head_length=0.05,
          color=cols[1])#, **kwargs)


ax_labs.plot([0.0475,0.0875],[min_y,min_y], color=cols[1],linewidth=3)
ax_labs.plot([0.0475,0.0875],[max_y,max_y], color=cols[1],linewidth=3)

ax6.text(10**-1,0.7,'\\boldmath $F_f(\\max(v_{\\rm in}))$',fontsize=25,backgroundcolor='deepskyblue',color='w',ha='right')
ax6.text(10**-3.5,0.22,'\\boldmath $F_f(\\min(v_{\\rm in}))$',fontsize=25,backgroundcolor='teal',color='w',ha='left')
# ax6.text(10**-1,0.87,'\\boldmath $\\max(F_f)$',fontsize=25,backgroundcolor='deepskyblue',color='w',ha='right')
# ax6.text(10**-3.5,0.22,'\\boldmath $\\min(F_f)$',fontsize=25,backgroundcolor='teal',color='w',ha='left')


style = "Simple, tail_width=0.5, head_width=10, head_length=12"
kw1 = dict(arrowstyle=style, color="teal")
kw2 = dict(arrowstyle=style, color="deepskyblue")
a3 = patches.FancyArrowPatch((10**-3.5, 0.22), (kinvals[min_y_arg],min_y),
                             connectionstyle="arc3,rad=.5", **kw1)
a2 = patches.FancyArrowPatch((10**-1, 0.7), (kinvals[max_y_arg],max_y),
                             connectionstyle="arc3,rad=-.5", **kw2)

ax6.add_patch(a3)
ax6.add_patch(a2)

ax_labs.text(0.065,0.35,'\\boldmath $\\Delta F$',color=cols[1],fontsize=35,
         ha='center',va='center',backgroundcolor='w')
ax_labs.text(0.115,0.32,'\\boldmath $f$',color=cols[1],fontsize=25,
         ha='center',va='center')
ax_labs.text(0.16,0.36,'$=$',fontsize=35,ha='center',va='center')
ax_labs.text(0.275,0.38,'$-$',fontsize=35,ha='center',va='center')

r1=patches.Rectangle((0.3,0.32), 0.05, 0.07,color='teal')
r2=patches.Rectangle((0.19,0.32), 0.05, 0.07,color='deepskyblue')
ax_labs.add_patch(r1)
ax_labs.add_patch(r2)


if savefig==1:
    plt.savefig('Figure_2.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_2.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_2.svg',bbox_inches='tight',format='svg')
    plt.savefig('Figure_2.pdf',bbox_inches='tight',format='pdf')
    
plt.show()
    
    
    
    
    
    