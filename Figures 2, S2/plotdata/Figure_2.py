#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 12:19:16 2025

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



savefig=0           # Change to 1 in order to save generated figure
numkin=201 
wantorig=1
kinvals=np.logspace(-4,1,numkin)


Ks_top=10**(2.)
Ks_bottom=1.

Btot=10.
VMB=0.01
KeqB=1.

V24_vals=np.array([0.1,1,10])


rev=False
cmapuse='RdBu'


params = {'legend.fontsize': 30,
          'figure.figsize': (12, 9),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':20,
          'ytick.labelsize':20}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath}')

plt.figure()
ax1=plt.axes([0.03,0.55,0.62,0.5])
motifnone=img.imread('fig2_none_tr.png')
plt.imshow(motifnone)
plt.axis('off')

ax2=plt.axes([0.68,0.55,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('fig2_minimal_tr.png')
plt.imshow(motifbranch)

ax3=plt.axes([1.33,0.55,0.62,0.5])
plt.axis('off')

motifbranch=img.imread('fig2_branch_upstream_tr.png')
plt.imshow(motifbranch)


ax4=plt.axes([0.1,0.05,0.5,0.5])
ax8=plt.axes([0.75,0.05,0.5,0.5])
ax6=plt.axes([1.4,0.05,0.5,0.5])

inset3=plt.axes([1.085,0.12,0.15,0.15])
inset2=plt.axes([1.73,0.15,0.15,0.15])


ax9=plt.axes([0.75,0.05,0.5,0.5])
plt.axis('off')


ax10=plt.axes([0.1,-0.55,0.5,0.5])
ax11=plt.axes([0.75,-0.55,0.5,0.5])
ax13=plt.axes([1.4,-0.55,0.5,0.5])


    
cols=['blue','black','red']



al=[1., 1., 1., 1.]  # pyv_23, b1_23, pyv_24, b0_24

for j in range(len(V24_vals)):
    V24=V24_vals[j]
    print(V24)
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[0],[0]),
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
                          params_branch=[1.,Ks_top,1.,1.,V24,Ks_bottom,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
        
    gradient=np.abs(np.gradient(ff_vec, kinvals))
    sense=gradient/(ff_vec/kinvals)
    ax10.plot(kinvals,gradient,color=cols[j])
    
    ax4.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$')
ax4.set_xscale('log')
ax4.legend(loc=(0.05,0.15),fontsize='large')
ax4.set_ylim(0,1.05)
ax4.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax4.xaxis.set_label_coords(0.8,-0.05)
ax4.set_ylabel('$\\rm Flux~ Fraction$',fontsize=40)
ax4.yaxis.set_label_coords(-0.1,0.4)
ax4.text(10**-4.6,0.86,'$\\rm D$', fontsize=40)
ax4.text(10**-4.5,1.75,'$\\rm A$', fontsize=40)

ax10.set_xscale('log')
ax10.set_ylim(0,200)
ax10.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax10.xaxis.set_label_coords(0.8,-0.05)
ax10.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax10.yaxis.set_label_coords(-0.15,0.25)
ax10.text(10**-4.6,160,'$\\rm G$', fontsize=40)
ax10.set_yticks(np.linspace(0,200,5))



for j in range(len(V24_vals)):
    print(j)
    V24=V24_vals[j]
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    b_vec=np.zeros(np.size(kinvals))
    b_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
                          params_branch=[1.,Ks_top,1.,1.,V24,Ks_bottom,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
            b_vec[i]=1.-bfrac
    
    gradient=np.abs(np.gradient(ff_vec, kinvals))
    sense=gradient/(ff_vec/kinvals)
    ax13.plot(kinvals,gradient,color=cols[j])
    
    ax6.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$')
    inset2.plot(kinvals,b_vec,color=cols[j])
    
ax6.set_xscale('log')
ax6.legend(loc=(0.7,0.6),fontsize='large')
ax6.set_ylim(0,1.05)
ax6.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax6.xaxis.set_label_coords(0.8,-0.05)
ax6.text(10**-4.6,0.86,'$\\rm F$', fontsize=40)
ax6.set_ylabel('$\\rm Flux~ Fraction$',fontsize=40)
ax6.yaxis.set_label_coords(-0.1,0.4)
ax6.text(10**-4.5,1.75,'$\\rm C$', fontsize=40)


ax13.set_xscale('log')
ax13.set_ylim(0,200)
ax13.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax13.xaxis.set_label_coords(0.8,-0.05)
ax13.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax13.yaxis.set_label_coords(-0.15,0.25)
ax13.text(10**-4.6,160,'$\\rm I$', fontsize=40)
ax13.set_yticks(np.linspace(0,200,5))




inset2.set_xscale('log')
inset2.set_ylim(0,1)
inset2.set_xticks(np.logspace(-4,0,3))
inset2.set_xlabel('$v_{\\rm in}$',fontsize=25)
inset2.xaxis.set_label_coords(0.75,-0.15)
inset2.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=25)
inset2.yaxis.set_label_coords(-0.35,0.5)



for j in range(len(V24_vals)):
    print(j)
    V24=V24_vals[j]
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    b_vec=np.zeros(np.size(kinvals))
    b_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[0],[1]),
                          params_global=[1.,1.,1.,1.],
                          params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
                          params_branch=[1.,Ks_top,1.,1.,V24,Ks_bottom,1.,1.],
                          params_inout=[kin,1.,1.],
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
            b_vec[i]=1.-bfrac
        if cols[j]=='black':
            ffvec_save=ff_vec
    
    gradient=np.abs(np.gradient(ff_vec, kinvals))
    sense=gradient/(ff_vec/kinvals)
    ax11.plot(kinvals,gradient,color=cols[j])
    
    ax8.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$')
    inset3.plot(kinvals,b_vec,color=cols[j])
  
ax8.set_xscale('log')
ax8.legend(loc=(0.7,0.6),fontsize='large')
ax8.set_ylim(0,1.05)
ax8.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax8.xaxis.set_label_coords(0.8,-0.05)
ax8.text(10**-4.6,0.86,'$\\rm E$', fontsize=40)
ax8.set_ylabel('$\\rm Flux~ Fraction$',fontsize=40)
ax8.yaxis.set_label_coords(-0.1,0.4)
ax8.text(10**-4.5,1.75,'$\\rm B$', fontsize=40)


ax11.set_xscale('log')
ax11.set_ylim(0,200)
ax11.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax11.xaxis.set_label_coords(0.8,-0.05)
ax11.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax11.yaxis.set_label_coords(-0.15,0.25)
ax11.text(10**-4.6,160,'$\\rm H$', fontsize=40)
ax11.set_yticks(np.linspace(0,200,5))


min_y=np.nanmin(ffvec_save)
print(min_y)
ax8.plot([10**-4,10**0],[min_y,min_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')

max_y=np.nanmax(ffvec_save)
print(max_y)
ax8.plot([10**-4,10**0],[max_y,max_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')




ax9.set_ylim(0,1.05)
ax9.set_xlim(0,1)

start1=0.57
start2=0.4

ax9.arrow(0.3, start1, 0, max_y-start1,length_includes_head=True,
          width=0.01,head_width=0.05,head_length=0.07,
          color=[0.5,0.5,0.5])#, **kwargs)

ax9.arrow(0.3, start2, 0, min_y-start2,length_includes_head=True,
          width=0.01,head_width=0.05,head_length=0.07,
          color=[0.5,0.5,0.5])#, **kwargs)

ax9.text(0.3,0.47,'$\\Delta F_f$',color=[0.5,0.5,0.5],fontsize=50,
         ha='center',va='center')

inset3.set_xscale('log')
inset3.set_ylim(0,1)
inset3.set_xticks(np.logspace(-4,0,3))
inset3.set_xlabel('$v_{\\rm in}$',fontsize=25)
inset3.xaxis.set_label_coords(0.75,-0.15)
inset3.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=25)
inset3.yaxis.set_label_coords(-0.35,0.5)


if wantorig==1:
    al=[1., 1., 1., 1.] 
        
    linst=':'
    for j in range(len(V24_vals)):
        V24=V24_vals[j]
        print(V24)
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[0],[0]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
                              params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=rev,Tmax=100000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=1.-fluxfrac
            
        gradient=np.abs(np.gradient(ff_vec, kinvals))
        sense=gradient/(ff_vec/kinvals)
        ax10.plot(kinvals,gradient,color=cols[j],linestyle=linst)
        
        ax4.plot(kinvals,ff_vec,color=cols[j],linestyle=linst)
        
    
   
    
    
    
    
    for j in range(len(V24_vals)):
        print(j)
        V24=V24_vals[j]
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        b_vec=np.zeros(np.size(kinvals))
        b_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
                              params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=rev,Tmax=100000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=1.-fluxfrac
                b_vec[i]=1.-bfrac
        
        gradient=np.abs(np.gradient(ff_vec, kinvals))
        sense=gradient/(ff_vec/kinvals)
        ax13.plot(kinvals,gradient,color=cols[j],linestyle=linst)
        
        ax6.plot(kinvals,ff_vec,color=cols[j],linestyle=linst)
        inset2.plot(kinvals,b_vec,color=cols[j],linestyle=':')

   
    
    for j in range(len(V24_vals)):
        print(j)
        V24=V24_vals[j]
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        b_vec=np.zeros(np.size(kinvals))
        b_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[0],[1]),
                              params_global=[1.,1.,1.,1.],
                              params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
                              params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                              params_inout=[kin,1.,1.],
                              reversible=rev,Tmax=100000,tpts=10000)
            if sum(isbuildup)==0:
                ff_vec[i]=1.-fluxfrac
                b_vec[i]=1.-bfrac
            if cols[j]=='black':
                ffvec_save=ff_vec
        
        gradient=np.abs(np.gradient(ff_vec, kinvals))
        sense=gradient/(ff_vec/kinvals)
        ax11.plot(kinvals,gradient,color=cols[j],linestyle=linst)
        
        ax8.plot(kinvals,ff_vec,color=cols[j],linestyle=linst)
        inset3.plot(kinvals,b_vec,color=cols[j],linestyle=':')

  
    
          
      
    
if savefig==1:
    plt.savefig('Figure_2.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_2.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_2.svg',bbox_inches='tight',format='svg')
    
plt.show()
