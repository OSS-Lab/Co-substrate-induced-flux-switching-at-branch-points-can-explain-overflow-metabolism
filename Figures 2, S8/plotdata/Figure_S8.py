#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 11:53:45 2025

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
import matplotlib.patches as patches



savefig=0       # Change to 1 in order to save generated figure
# numkin=11
numkin=201 
wantorig=1
kinvals=np.logspace(-4,1,numkin)

V24_vals=np.array([0.1,1,10])
    

Btot=10.
VMB=0.01
KeqB=1.

rev=False



al=[2., 1., 1., 1.]  # pyv_23, b1_23, pyv_24, b0_24
cmapuse='RdBu'


params = {'legend.fontsize': 30,
          'figure.figsize': (12, 9),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':25,
          'ytick.labelsize':25}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath}')

plt.figure()
ax1=plt.axes([0.03,0.55,0.62,0.5])
motifnone=img.imread('fig2_none_tr.png')
plt.imshow(motifnone)
plt.axis('off')

ax2=plt.axes([0.65,0.55,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('fig2_minimal_tr.png')
plt.imshow(motifbranch)

ax3=plt.axes([1.27,0.55,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('figS8_branch_tr.png')
plt.imshow(motifbranch)

ax7=plt.axes([1.89,0.55,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('fig2_branch_upstream_tr.png')
plt.imshow(motifbranch)


ax4=plt.axes([0.04,0,0.6,0.6])
ax5=plt.axes([1.28,0,0.6,0.6])
ax6=plt.axes([1.9,0,0.6,0.6])
ax8=plt.axes([0.66,0,0.6,0.6])

inset3=plt.axes([1.04,0.07,0.2,0.2])
inset1=plt.axes([1.67,0.385,0.2,0.2])
inset2=plt.axes([2.285,0.07,0.2,0.2])



# ax9=plt.axes([0.66,0,0.6,0.6])
# plt.axis('off')


ax10=plt.axes([0.04,-0.62,0.6,0.6])
ax11=plt.axes([0.66,-0.62,0.6,0.6])
ax12=plt.axes([1.28,-0.62,0.6,0.6])
ax13=plt.axes([1.9,-0.62,0.6,0.6])

cols=['blue','olive','red']

lw=3
legfont=25
# no cosubs
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
    ax10.plot(kinvals,gradient,color=cols[j],linewidth=lw)
    
    ax4.plot(kinvals,ff_vec,color=cols[j],label='$V_{\\rm max}^{24} = '+str(V24)+'$',linewidth=lw)
ax4.set_xscale('log')
# ax4.legend(loc=(0.05,0.52),fontsize=legfont)
ax4.set_ylim(0,1.05)
# ax4.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax4.xaxis.set_label_coords(0.8,-0.05)
# ax4.set_ylabel('$\\rm Flux~ Fraction$',fontsize=40)
ax4.set_ylabel('$F_f$',fontsize=40,rotation=0)
ax4.yaxis.set_label_coords(-0.1,0.42)
ax4.text(10**-4.1,0.95,'$\\rm A(ii)$', fontsize=40)
ax4.text(10**-4.1,1.6,'$\\rm A(i)$', fontsize=40)
ax4.text(10**-4.1, 1.43, '$v_{\\rm in}$',fontsize=30)

ax4.text(10**-4,0.3,'$\\frac{}{V_{\\rm max}^{24}} =$',fontsize=40)
ax4.text(10**-4,0.36,'$V_{\\rm max}^{23} $',fontsize=29)
ax4.text(10**-2.95,0.365,'$10$',fontsize=25,horizontalalignment='right',color='b')
ax4.text(10**-2.95,0.305,'$1$',fontsize=25,horizontalalignment='right',color='olive')
ax4.text(10**-2.95,0.245,'$0.1$',fontsize=25,horizontalalignment='right',color='r')

ax4.text(10**-3.3,0.77,'${\\rm without~allostery}$',fontsize=25)
ax4.text(10**-3.5,0.55,'${\\rm with~allostery}$',fontsize=25)

style = "Simple, tail_width=0.5, head_width=10, head_length=12"
kw3=dict(arrowstyle=style, color="black")

a4 = patches.FancyArrowPatch((10**-2.2, 0.81), (10**-2.2,0.9),
                             **kw3)
a5 = patches.FancyArrowPatch((10**-2.35, 0.57), (10**-1.7,0.57),
                             **kw3)

ax4.add_patch(a4)
ax4.add_patch(a5)


ax10.set_xscale('log')
ax10.set_ylim(0,180)
ax10.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax10.xaxis.set_label_coords(0.8,-0.05)
# ax10.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax10.set_ylabel('$ \\left\\lvert \\frac{d  F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax10.yaxis.set_label_coords(-0.1,0.35)
ax10.text(10**-4.1,160,'$\\rm A(iii)$', fontsize=40)



for j in range(len(V24_vals)):
    print(j)
    V24=V24_vals[j]
    ff_vec=np.zeros(np.size(kinvals))
    ff_vec[:]=np.nan
    b_vec=np.zeros(np.size(kinvals))
    b_vec[:]=np.nan
    for i in range(len(kinvals)):
        kin=kinvals[i]
        
        isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
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
    ax12.plot(kinvals,gradient,color=cols[j],linewidth=lw)
    
    ax5.plot(kinvals,ff_vec,color=cols[j],label='$V_{\\rm max}^{24} = '+str(V24)+'$',linewidth=lw)
    inset1.plot(kinvals,b_vec,color=cols[j],linewidth=lw)
    
ax5.set_xscale('log')
# ax5.legend(loc=(0.05,0.5),fontsize=legfont)
ax5.set_ylim(0,1.05)
# ax5.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax5.xaxis.set_label_coords(0.8,-0.05)
ax5.text(10**-4.1,0.95,'$\\rm C(ii)$', fontsize=40)
ax5.text(10**-4.1,1.6,'$\\rm C(i)$', fontsize=40)
ax5.text(10**-4.1, 1.43, '$v_{\\rm in}$',fontsize=30)


ax12.set_xscale('log')
ax12.set_ylim(0,180)
ax12.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax12.xaxis.set_label_coords(0.8,-0.05)
# ax12.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax12.yaxis.set_label_coords(-0.15,0.3)
ax12.text(10**-4.1,160,'$\\rm C(iii)$', fontsize=40)



inset1.set_xscale('log')
inset1.set_ylim(-0.05,1)
inset1.set_xlabel('$v_{\\rm in}$',fontsize=30)
inset1.xaxis.set_label_coords(0.8,-0.1)
inset1.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
inset1.yaxis.set_label_coords(-0.15,0.4)





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
    ax13.plot(kinvals,gradient,color=cols[j],linewidth=lw)
    
    ax6.plot(kinvals,ff_vec,color=cols[j],label='$V_{\\rm max}^{24} = '+str(V24)+'$',linewidth=lw)
    inset2.plot(kinvals,b_vec,color=cols[j],linewidth=lw)
    
ax6.set_xscale('log')
# ax6.legend(loc=(0.7,0.6),fontsize=legfont)
ax6.set_ylim(0,1.05)
# ax6.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax6.xaxis.set_label_coords(0.8,-0.05)
ax6.text(10**-4.1,0.95,'$\\rm D(ii)$', fontsize=40)
ax6.text(10**-4.1,1.6,'$\\rm D(i)$', fontsize=40)
ax6.text(10**-4.1, 1.43, '$v_{\\rm in}$',fontsize=30)


ax13.set_xscale('log')
ax13.set_ylim(0,180)
ax13.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax13.xaxis.set_label_coords(0.8,-0.05)
# ax13.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax13.yaxis.set_label_coords(-0.15,0.3)
ax13.text(10**-4.1,160,'$\\rm D(iii)$', fontsize=40)




inset2.set_xscale('log')
inset2.set_ylim(-0.05,1)
inset2.set_xlabel('$v_{\\rm in}$',fontsize=30)
inset2.xaxis.set_label_coords(0.8,-0.1)
inset2.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
inset2.yaxis.set_label_coords(-0.15,0.4)



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
        if cols[j]=='olive':
            ffvec_save=ff_vec
    
    gradient=np.abs(np.gradient(ff_vec, kinvals))
    sense=gradient/(ff_vec/kinvals)
    ax11.plot(kinvals,gradient,color=cols[j],linewidth=lw)
    
    ax8.plot(kinvals,ff_vec,color=cols[j],label='$V_{\\rm max}^{24} = '+str(V24)+'$',linewidth=lw)
    inset3.plot(kinvals,b_vec,color=cols[j],linewidth=lw)
  
ax8.set_xscale('log')
# ax8.legend(loc=(0.7,0.6),fontsize=legfont)
ax8.set_ylim(0,1.05)
# ax8.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax8.xaxis.set_label_coords(0.8,-0.05)
ax8.text(10**-4.1,0.95,'$\\rm B(ii)$', fontsize=40)
ax8.text(10**-4.1,1.6,'$\\rm B(i)$', fontsize=40)
ax8.text(10**-4.1, 1.43, '$v_{\\rm in}$',fontsize=30)


ax11.set_xscale('log')
ax11.set_ylim(0,180)
ax11.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax11.xaxis.set_label_coords(0.8,-0.05)
# ax11.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax11.yaxis.set_label_coords(-0.15,0.3)
ax11.text(10**-4.1,160,'$\\rm B(iii)$', fontsize=40)


# min_y=np.nanmin(ffvec_save)
# print(min_y)
# ax8.plot([10**-3,10**0],[min_y,min_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')

# max_y=np.nanmax(ffvec_save)
# print(max_y)
# ax8.plot([10**-3,10**0],[max_y,max_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')




# ax9.set_ylim(0,1.05)
# ax9.set_xlim(0,1)

# start1=0.57
# start2=0.4

# ax9.arrow(0.3, start1, 0, max_y-start1,length_includes_head=True,
#           width=0.01,head_width=0.05,head_length=0.07,
#           color=[0.5,0.5,0.5])#, **kwargs)

# ax9.arrow(0.3, start2, 0, min_y-start2,length_includes_head=True,
#           width=0.01,head_width=0.05,head_length=0.07,
#           color=[0.5,0.5,0.5])#, **kwargs)

# ax9.text(0.3,0.47,'$\\Delta F_f$',color=[0.5,0.5,0.5],fontsize=50,
#          ha='center',va='center')

inset3.set_xscale('log')
inset3.set_ylim(-0.05,1)
inset3.set_xlabel('$v_{\\rm in}$',fontsize=30)
inset3.xaxis.set_label_coords(0.8,-0.1)
inset3.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
inset3.yaxis.set_label_coords(-0.15,0.4)

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
        ax10.plot(kinvals,gradient,color=cols[j],linestyle=linst,linewidth=lw)
        
        ax4.plot(kinvals,ff_vec,color=cols[j],linestyle=linst,linewidth=lw)
    
    
    
    for j in range(len(V24_vals)):
        print(j)
        V24=V24_vals[j]
        ff_vec=np.zeros(np.size(kinvals))
        ff_vec[:]=np.nan
        b_vec=np.zeros(np.size(kinvals))
        b_vec[:]=np.nan
        for i in range(len(kinvals)):
            kin=kinvals[i]
            
            isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[1],[1]),
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
        ax12.plot(kinvals,gradient,color=cols[j],linestyle=linst,linewidth=lw)
        
        ax5.plot(kinvals,ff_vec,color=cols[j],linestyle=linst,linewidth=lw)
        inset1.plot(kinvals,b_vec,color=cols[j],linestyle=':',linewidth=lw)
   
    
    
    
    
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
        ax13.plot(kinvals,gradient,color=cols[j],linestyle=linst,linewidth=lw)
        
        ax6.plot(kinvals,ff_vec,color=cols[j],linestyle=linst,linewidth=lw)
        inset2.plot(kinvals,b_vec,color=cols[j],linestyle=':',linewidth=lw)
   
    
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
            if cols[j]=='olive':
                ffvec_save=ff_vec
        
        gradient=np.abs(np.gradient(ff_vec, kinvals))
        sense=gradient/(ff_vec/kinvals)
        ax11.plot(kinvals,gradient,color=cols[j],linestyle=linst,linewidth=lw)
        
        ax8.plot(kinvals,ff_vec,color=cols[j],linestyle=linst,linewidth=lw)
        inset3.plot(kinvals,b_vec,color=cols[j],linestyle=':',linewidth=lw)
  
    
ax4.set_yticks(np.linspace(0,1,6),labels=None)
ax4.set_xticks(np.logspace(-4,0,5),labels=['','','','',''])

ax5.set_yticks(np.linspace(0,1,6),labels=['','','','','',''])
ax5.set_xticks(np.logspace(-4,0,5),labels=['','','','',''])

ax6.set_yticks(np.linspace(0,1,6),labels=['','','','','',''])
ax6.set_xticks(np.logspace(-4,0,5),labels=['','','','',''])

ax8.set_yticks(np.linspace(0,1,6),labels=['','','','','',''])
ax8.set_xticks(np.logspace(-4,0,5),labels=['','','','',''])

ax10.set_yticks(np.linspace(0,150,4))
ax10.set_xticks(np.logspace(-4,0,5))

ax11.set_yticks(np.linspace(0,150,4),labels=['','','',''])
ax11.set_xticks(np.logspace(-4,0,5))

ax12.set_yticks(np.linspace(0,150,4),labels=['','','',''])
ax12.set_xticks(np.logspace(-4,0,5))

ax13.set_yticks(np.linspace(0,150,4),labels=['','','',''])
ax13.set_xticks(np.logspace(-4,0,5))

inset1.set_yticks(np.linspace(0,1,5),labels=['$0$','','','','$1$'])
inset2.set_yticks(np.linspace(0,1,5),labels=['$0$','','','','$1$'])
inset3.set_yticks(np.linspace(0,1,5),labels=['$0$','','','','$1$'])

inset1.set_xticks(np.logspace(-4,0,3))
inset2.set_xticks(np.logspace(-4,0,3))
inset3.set_xticks(np.logspace(-4,0,3))
    
if savefig==1:
    plt.savefig('Figure_S8.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S8.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S8.svg',bbox_inches='tight',format='svg')
    
plt.show()
