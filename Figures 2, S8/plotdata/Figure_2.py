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
# import matplotlib.patheffects as pe
import matplotlib.patches as patches

savefig=0          # Change to 1 in order to save generated figure
numkin=201
# numkin=21
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
          'xtick.labelsize':25,
          'ytick.labelsize':25}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath} \usepackage{color}')

plt.figure()
ax1=plt.axes([0.091,0.65,0.62,0.5])
motifnone=img.imread('fig2_none_tr.png')
plt.imshow(motifnone)
plt.axis('off')

ax2=plt.axes([0.73,0.65,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('fig2_minimal_tr.png')
plt.imshow(motifbranch)

# ax3=plt.axes([1.368,0.65,0.62,0.5])
# plt.axis('off')

# motifbranch=img.imread('fig2_branch_upstream_tr.png')
# plt.imshow(motifbranch)


ax4=plt.axes([0.1,0.08,0.6,0.6])
ax8=plt.axes([0.74,0.08,0.6,0.6])
# ax6=plt.axes([1.38,0.08,0.6,0.6])

inset3=plt.axes([1.12,0.18,0.2,0.2])
# inset2=plt.axes([1.75,0.18,0.2,0.2])




ax10=plt.axes([0.1,-0.55,0.6,0.6])
ax11=plt.axes([0.74,-0.55,0.6,0.6])
# ax13=plt.axes([1.38,-0.55,0.6,0.6])

ax9=plt.axes([0.1,0.08,0.8,0.6])
ax9.axis('off')

    
cols=['blue','olive','red']
lw=3


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
        
    if cols[j]=='blue':
         ffvec_save=ff_vec   
    gradient=np.abs(np.gradient(ff_vec, kinvals))
    sense=gradient/(ff_vec/kinvals)
    ax10.plot(kinvals,gradient,color=cols[j],linewidth=lw)
    
    ax4.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$',linewidth=lw)
ax4.set_xscale('log')
# ax4.legend(loc=(0.05,0.15),fontsize='large')
ax4.set_ylim(0,1.05)
# ax4.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax4.xaxis.set_label_coords(0.8,-0.05)
ax4.set_ylabel('${\\rm Fraction ~of ~upper ~branch ~flux}, F_f$',fontsize=30)
ax4.yaxis.set_label_coords(-0.1,0.5)
ax4.text(10**-4.1,0.94,'$\\rm A(ii)$', fontsize=50)
ax4.text(10**-4.1,1.62,'$\\rm A(i)$', fontsize=50)
ax4.set_xticks(np.logspace(-4,0,5),labels=['','','','',''])

ax10.set_xscale('log')
ax10.set_ylim(0,110)
ax10.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax10.xaxis.set_label_coords(0.7,-0.05)
# ax10.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax10.set_ylabel('${\\rm Degree ~of ~non}$-${\\rm linearity}, \\left\\lvert \\frac{d F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=30)
ax10.yaxis.set_label_coords(-0.08,0.5)
ax10.text(10**-4.1,95,'$\\rm A(iii)$', fontsize=50)
ax10.set_yticks(np.linspace(0,100,5))

###### THIS IS THE DELTA F_f LABEL BIT
min_y=np.nanmin(ffvec_save)
print(min_y)
# ax4.plot([10**-4,10**0],[min_y,min_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')

max_y=np.nanmax(ffvec_save)
print(max_y)
# ax4.plot([10**-4,10**0],[max_y,max_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')


ax9.set_ylim(0,1.05)
ax9.set_xlim(0,1)

start1=0.3
start2=0.2

ax9.arrow(0.75, start1, 0, max_y-start1,length_includes_head=True,
          width=0.005,head_width=0.03,head_length=0.05,
          color='b')#, **kwargs)

ax9.arrow(0.75, start2, 0, min_y-start2,length_includes_head=True,
          width=0.005,head_width=0.03,head_length=0.05,
          color='b')#, **kwargs)



ax9.plot([0.73,0.77],[min_y,min_y],'b',linewidth=3)
ax9.plot([0.73,0.77],[max_y,max_y],'b',linewidth=3)

ax4.text(10**-1.8,0.97,'\\boldmath $F_f(\\max(v_{\\rm in}))$',fontsize=25,backgroundcolor='deepskyblue',color='w')
         # path_effects=[pe.withStroke(linewidth=1,foreground="k")])

ax4.text(10**-4,0.25,'\\boldmath $F_f(\\min(v_{\\rm in}))$',fontsize=25,backgroundcolor='teal',color='w')
         # path_effects=[pe.withStroke(linewidth=1,foreground="k")])

ax9.text(0.74,0.25,'\\boldmath $\\Delta F$',color='b',fontsize=35,
         ha='center',va='center',backgroundcolor='w')
ax9.text(0.78,0.22,'\\boldmath $f$',color='b',fontsize=25,
         ha='center',va='center')

style = "Simple, tail_width=0.5, head_width=10, head_length=12"
kw1 = dict(arrowstyle=style, color="teal")
kw2 = dict(arrowstyle=style, color="deepskyblue")

# a1 = patches.FancyArrowPatch((-0.4, -0.6), (0, 0.6), **kw)
# a2 = patches.FancyArrowPatch((0, 0.6), (0.4, -0.6), **kw)
a3 = patches.FancyArrowPatch((10**-4, 0.25), (10**-3.9,min_y),
                             connectionstyle="arc3,rad=.5", **kw1)
a2 = patches.FancyArrowPatch((10**-0.5, 0.97), (10**0,max_y),
                             connectionstyle="arc3,rad=-.5", **kw2)

ax4.add_patch(a3)
ax4.add_patch(a2)

r1=patches.Rectangle((0.631,0.22), 0.035, 0.07,color='teal')
r2=patches.Rectangle((0.562,0.22), 0.035, 0.07,color='deepskyblue')
ax9.add_patch(r1)
ax9.add_patch(r2)


# ax9.text(0.67,0.24,'\\bf{=}',fontsize=25)
# ax9.text(0.602,0.25,'\\boldmath $-$',fontsize=25)
ax9.plot([0.6035,0.6235],[0.255,0.255],'-k',linewidth=2)

ax9.plot([0.672,0.692],[0.2475,0.2475],'-k',linewidth=2)
ax9.plot([0.672,0.692],[0.2625,0.2625],'-k',linewidth=2)

# ax4.text(10**-4,0.7,'$\\frac{V_{\\rm max}^{23}}{V_{\\rm max}^{24}} =$',fontsize=40)
ax4.text(10**-4,0.7,'$\\frac{}{V_{\\rm max}^{24}} =$',fontsize=40)
ax4.text(10**-4,0.76,'$V_{\\rm max}^{23} $',fontsize=29)
ax4.text(10**-2.9,0.785,'$10$',fontsize=25,horizontalalignment='right',color='b')
ax4.text(10**-2.9,0.71,'$1$',fontsize=25,horizontalalignment='right',color='olive')
ax4.text(10**-2.9,0.635,'$0.1$',fontsize=25,horizontalalignment='right',color='r')

ax4.text(10**-2.21,0.79,'$K_{\\rm M}^{23} = K_{\\rm M}^{24}$',fontsize=25)
ax4.text(10**-2.6,0.55,'$K_{\\rm M}^{23} = 100 \\cdot K_{\\rm M}^{ 24}$',fontsize=25)

kw3=dict(arrowstyle=style, color="black")
style = "Simple, tail_width=0.5, head_width=8, head_length=10"

a4 = patches.FancyArrowPatch((10**-1.75, 0.82), (10**-1.75,0.9),
                             **kw3)
a5 = patches.FancyArrowPatch((10**-1.2, 0.57), (10**-0.7,0.57),
                             **kw3)

ax4.add_patch(a4)
ax4.add_patch(a5)

ax4.text(10**-4.2, 1.47, '$v_{\\rm in}$',fontsize=40)
ax8.text(10**-4.15, 1.47, '$v_{\\rm in}$',fontsize=40)
# ax6.text(10**-4.15, 1.47, '$v_{\\rm in}$',fontsize=40)
         #'$\\textcolor{blue} 10$ \n $\\textcolor{black} 1$ \n $\\textcolor{red} 0.1$',fontsize=25,horizontalalignment='right')


# ax9.annotate("", xytext=(0.75, start2), xy=(0.75, min_y),
#             arrowprops=dict(arrowstyle="simple,head_width=0.8, head_length=0.5",
#                             color='b'))
             # arrowprops=dict(arrowstyle="-|>, head_width=0.8, head_length=2",
            #                 linewidth=5,
            #                 color='b'))

# for j in range(len(V24_vals)):
#     print(j)
#     V24=V24_vals[j]
#     ff_vec=np.zeros(np.size(kinvals))
#     ff_vec[:]=np.nan
#     b_vec=np.zeros(np.size(kinvals))
#     b_vec[:]=np.nan
#     for i in range(len(kinvals)):
#         kin=kinvals[i]
        
#         isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
#                           params_global=[1.,1.,1.,1.],
#                           params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
#                           params_branch=[1.,Ks_top,1.,1.,V24,Ks_bottom,1.,1.],
#                           params_inout=[kin,1.,1.],
#                           reversible=rev,Tmax=100000,tpts=10000)
#         if sum(isbuildup)==0:
#             ff_vec[i]=1.-fluxfrac
#             b_vec[i]=1.-bfrac
    
#     gradient=np.abs(np.gradient(ff_vec, kinvals))
#     sense=gradient/(ff_vec/kinvals)
#     ax13.plot(kinvals,gradient,color=cols[j],linewidth=lw)
    
#     ax6.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$',linewidth=lw)
#     inset2.plot(kinvals,b_vec,color=cols[j],linewidth=lw)
    
# ax6.set_xscale('log')
# # ax6.legend(loc=(0.7,0.6),fontsize='large')
# ax6.set_ylim(0,1.05)
# # ax6.set_xlabel('$v_{\\rm in}$',fontsize=30)
# ax6.xaxis.set_label_coords(0.8,-0.05)
# ax6.text(10**-4.1,0.94,'$\\rm C(ii)$', fontsize=50)
# # ax6.set_ylabel('$\\rm Flux~ Fraction$',fontsize=40)
# ax6.yaxis.set_label_coords(-0.1,0.4)
# ax6.text(10**-4.1,1.62,'$\\rm C(i)$', fontsize=50)
# ax6.set_xticks(np.logspace(-4,0,5),labels=['','','','',''])
# ax6.set_yticks(np.linspace(0,1,6),labels=['','','','','',''])



# ax13.set_xscale('log')
# ax13.set_ylim(0,205)
# ax13.set_xlabel('$v_{\\rm in}$',fontsize=40)
# ax13.xaxis.set_label_coords(0.7,-0.05)
# # ax13.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
# ax13.yaxis.set_label_coords(-0.15,0.25)
# ax13.text(10**-4.1,180,'$\\rm C(iii)$', fontsize=50)
# ax13.set_yticks(np.linspace(0,200,5),labels=['','','','',''])




# inset2.set_xscale('log')
# inset2.set_ylim(0,1)
# inset2.set_xticks(np.logspace(-4,0,3))
# inset2.set_yticks(np.linspace(0,1,3),labels=['$0$','','$1$'])
# inset2.set_xlabel('$v_{\\rm in}$',fontsize=35)
# inset2.xaxis.set_label_coords(0.75,-0.15)
# inset2.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
# inset2.yaxis.set_label_coords(-0.15,0.3)



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
        if cols[j]=='olive':
            ffvec_save=ff_vec
    
    gradient=np.abs(np.gradient(ff_vec, kinvals))
    sense=gradient/(ff_vec/kinvals)
    ax11.plot(kinvals,gradient,color=cols[j],linewidth=lw)
    
    ax8.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$',linewidth=lw)
    inset3.plot(kinvals,b_vec,color=cols[j],linewidth=lw)
  
ax8.set_xscale('log')
# ax8.legend(loc=(0.7,0.6),fontsize='large')
ax8.set_ylim(0,1.05)
# ax8.set_xlabel('$v_{\\rm in}$',fontsize=30)
ax8.xaxis.set_label_coords(0.8,-0.05)
ax8.text(10**-4.1,0.94,'$\\rm B(ii)$', fontsize=50)
# ax8.set_ylabel('${\\rm fraction ~of ~upper ~branch ~flux}, F_f$',fontsize=40)
ax8.yaxis.set_label_coords(-0.1,0.4)
ax8.text(10**-4.1,1.62,'$\\rm B(i)$', fontsize=50)
ax8.set_xticks(np.logspace(-4,0,5),labels=['','','','',''])
ax8.set_yticks(np.linspace(0,1,6),labels=['','','','','',''])


ax11.set_xscale('log')
ax11.set_ylim(0,110)
ax11.set_xlabel('$v_{\\rm in}$',fontsize=40)
ax11.xaxis.set_label_coords(0.7,-0.05)
# ax11.set_ylabel('$ \\left\\lvert \\frac{d \\Delta F_f}{dv_{\\rm in}} \\right\\rvert$',fontsize=40,rotation=0)
ax11.yaxis.set_label_coords(-0.15,0.25)
ax11.text(10**-4.1,95,'$\\rm B(iii)$', fontsize=50)
ax11.set_yticks(np.linspace(0,100,5),labels=['','','','',''])





inset3.set_xscale('log')
inset3.set_ylim(0,1)
inset3.set_xticks(np.logspace(-4,0,3))
inset3.set_xlabel('$v_{\\rm in}$',fontsize=35)
inset3.xaxis.set_label_coords(0.75,-0.15)
inset3.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
inset3.yaxis.set_label_coords(-0.15,0.3)
inset3.set_yticks(np.linspace(0,1,3),labels=['$0$','','$1$'])


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
        
    
   
    
    
    
    
    # for j in range(len(V24_vals)):
    #     print(j)
    #     V24=V24_vals[j]
    #     ff_vec=np.zeros(np.size(kinvals))
    #     ff_vec[:]=np.nan
    #     b_vec=np.zeros(np.size(kinvals))
    #     b_vec[:]=np.nan
    #     for i in range(len(kinvals)):
    #         kin=kinvals[i]
            
    #         isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([1,0],[0,0],[1],[1]),
    #                           params_global=[1.,1.,1.,1.],
    #                           params_B=[VMB,1.,1.,KeqB,Btot],allo=al,
    #                           params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
    #                           params_inout=[kin,1.,1.],
    #                           reversible=rev,Tmax=100000,tpts=10000)
    #         if sum(isbuildup)==0:
    #             ff_vec[i]=1.-fluxfrac
    #             b_vec[i]=1.-bfrac
        
    #     gradient=np.abs(np.gradient(ff_vec, kinvals))
    #     sense=gradient/(ff_vec/kinvals)
    #     ax13.plot(kinvals,gradient,color=cols[j],linestyle=linst,linewidth=lw)
        
    #     ax6.plot(kinvals,ff_vec,color=cols[j],linestyle=linst,linewidth=lw)
    #     inset2.plot(kinvals,b_vec,color=cols[j],linestyle=':',linewidth=lw)

   
    
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

  
    
          
      
    
if savefig==1:
    plt.savefig('Figure_2.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_2.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_2.svg',bbox_inches='tight',format='svg')
    plt.savefig('Figure_2.pdf',bbox_inches='tight',format='pdf')
    
plt.show()
