#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 19:50:05 2024

@author: robert
"""



# import os
import sys
# import inspect

# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# sys.path.insert(0, parentdir) 

sys.path.append('../simulation_functions')
sys.path.append('../data/heatmapdata')

import numpy as np

import matplotlib.pyplot as plt

# import make_odes as mp
#import time_series_plots as ts
import matplotlib.pylab as pylab
import make_figure_data as mfd
import matplotlib.image as img
import matplotlib as mpl



savefig=0          # change to 1 to save generated figure
# numkin=11
numkin=201 

kinvals=np.logspace(-4,1,numkin)

# mfd.FIGURE2_2_ff_vs_kin_vary24()

rev=False
cmapuse='RdBu_r'

lw=3

params = {'legend.fontsize': 30,
          'figure.figsize': (12, 9),
          'axes.labelsize': 30,
          'axes.titlesize':'x-large',
          'xtick.labelsize':30,
          'ytick.labelsize':30}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

# binoutvals=[10**-2,10**0,10**0] #b0in b0out b1out
binoutvals=[10**-4,10**-4,10**-4] #b0in b0out b1out
    
plt.figure()
ax1=plt.axes([0.03,0.55,0.62,0.5])
motifnone=img.imread('figS2_minimal_tr.png')
plt.imshow(motifnone)
plt.axis('off')

ax2=plt.axes([0.68,0.55,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('figS2_branch_tr.png')
plt.imshow(motifbranch)

ax3=plt.axes([1.33,0.55,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('figS2_branch_upstream_tr.png')
plt.imshow(motifbranch)

ax9=plt.axes([0.04,-0.1,0.6,0.7])
ax7=plt.axes([0.69,-0.1,0.6,0.7])
ax8=plt.axes([1.34,-0.1,0.6,0.7])

#ax7=plt.axes([1.98,0.55,0.62,0.5])
#plt.axis('off')
#motifbranch=img.imread('fig1_minimal_tr.png')
#plt.imshow(motifbranch)

ax6=plt.axes([0.04,-0.95,0.6,0.7])
ax4=plt.axes([0.69,-0.95,0.6,0.7])
ax5=plt.axes([1.34,-0.95,0.6,0.7])

inset3=plt.axes([0.42,-0.8,0.20,0.20])
inset1=plt.axes([1.07,-0.47,0.20,0.20])
inset2=plt.axes([1.7,-0.8,0.20,0.20])

# inset1=plt.axes([0.82,-0.27,0.15,0.15])
# inset2=plt.axes([1.47,-0.35,0.15,0.15])
# inset3=plt.axes([0.16,-0.265,0.15,0.15])

# ax4=plt.axes([0.1,-0.6,0.52,0.55])
# ax5=plt.axes([0.75,-0.6,0.52,0.55])
# ax6=plt.axes([1.4,-0.6,0.52,0.55])

# ax8=plt.axes([2.05,0.05,0.5,0.5])
# ax9=plt.axes([2.05,0.05,0.5,0.5])
# plt.axis('off')

V24_vals=np.array([0.1,1,10])
    
cols=['blue','olive','red']

# Btot=0.01
# VMB=0.01
# KeqB=1.


#vals_Btot=np.logspace(-2,0,3)
#vals_CKcB=np.logspace(-1,1,3)
#Keqb_vals=np.logspace(-2,0,3)

# Btot=0.01
# VMB=10**-1.
# KeqB=1.

Btot=10.# [0.1,1.,10.,100.] #   #[10]
VMB=0.01
KeqB=1.



xvar='b_in'
vals_b_in=np.logspace(-5,2,29)#np.logspace(-5,2,8)
xmin=np.min(np.log10(vals_b_in))
xmax=np.max(np.log10(vals_b_in))
xnum=len(vals_b_in)

yvar='b_out'
vals_b_out=np.logspace(-5,2,29)#np.logspace(-5,2,8)
ymin=np.min(np.log10(vals_b_out))
ymax=np.max(np.log10(vals_b_out))
ynum=len(vals_b_out)

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

kin_min=-8
kin_max=0
num_kin=41

Tmax=100000
tpts=10000
reversible=False

CKctop=1.
CKcbottom=1.

### ONLY need Keqb/t when reversible =true
vals_keqt=[1.]
vals_keqb=[1.]
Keq_top_vals=vals_keqt
Keq_bottom_vals=vals_keqb
lengthparams=[2,1,1]



B_position_params=([0,0],[0,0],[1],[1])
loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
        '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
        str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
        '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
        '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
        '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
        yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)

loadstring='ff_Btot_'+str(Btot)+'__CKcatB_'+str(VMB)+'__CKcat_bottom_'+str(CKcbottom)+\
    '__CKcat_top_'+str(CKctop)+'__KeqB_'+str(KeqB)+'__Keq_bottom_1.0__Keq_top_1.0__.npy'

ffdeltamat=np.load(loadpath+'/'+loadstring)
#print(ffdeltamat)
im=ax7.pcolormesh(vals_b_in,vals_b_out,-ffdeltamat,vmin=-1,vmax=1,cmap=cmapuse,edgecolors='face')
ax7.set_xscale('log')
ax7.set_yscale('log')





ax7.set_xlabel('${B_{\\rm 0}\\rm ~synthesis~rate}$, ${v^{B_0}_{\\rm in}}$')
ax7.xaxis.set_label_coords(0.5,-0.08)
# ax7.set_ylabel('${v^{B}_{\\rm deg}}$',fontsize=30)
# ax7.yaxis.set_label_coords(-0.1,0.6)
ax7.set_xticks(np.logspace(-4,2,4))
ax7.set_yticks(np.logspace(-4,2,4),['','','',''])
ax7.scatter(binoutvals[0],binoutvals[1],marker='o',color='olive',s=150)
ax7.text(10**-5,10**1.3,'$\\rm B(ii)$', fontsize=50)

# cbar = plt.colorbar(im, ticks=np.linspace(-1,1,5))
# #cbar.set_label('$\\Delta F_f$',fontsize=30,rotation=0,labelpad=-40, y=1.15)
# cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=0, y=0.9)

B_position_params=([1,0],[0,0],[1],[1])
loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
        '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
        str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
        '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
        '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
        '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
        yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)

loadstring='ff_Btot_'+str(Btot)+'__CKcatB_'+str(VMB)+'__CKcat_bottom_'+str(CKcbottom)+\
    '__CKcat_top_'+str(CKctop)+'__KeqB_'+str(KeqB)+'__Keq_bottom_1.0__Keq_top_1.0__.npy'

ffdeltamat=np.load(loadpath+'/'+loadstring)

im=ax8.pcolormesh(vals_b_in,vals_b_out,-ffdeltamat,vmin=-1,vmax=1,cmap=cmapuse,edgecolors='face')
ax8.set_xscale('log')
ax8.set_yscale('log')

ax8.set_xlabel('${B_{\\rm 0}\\rm ~synthesis~rate}$, ${v^{B_0}_{\\rm in}}$')
ax8.xaxis.set_label_coords(0.5,-0.08)
# ax7.set_ylabel('${v^{B}_{\\rm deg}}$',fontsize=30)
# ax7.yaxis.set_label_coords(-0.1,0.6)
ax8.set_xticks(np.logspace(-4,2,4))
ax8.set_yticks(np.logspace(-4,2,4),['','','',''])
ax8.scatter(binoutvals[0],binoutvals[1],marker='o',color='olive',s=150)
ax8.text(10**-5,10**1.3,'$\\rm C(ii)$', fontsize=50)


# cbar = plt.colorbar(im, ticks=np.linspace(-1,1,5))
# #cbar.set_label('$\\Delta F_f$',fontsize=30,rotation=0,labelpad=-40, y=1.15)
# cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=0, y=0.9)

cax = ax8.inset_axes([1.05, 0, 0.04, 0.85])
cbar=plt.colorbar(im, cax=cax, orientation='vertical',ticks=[-1, -0.5,0,0.5, 1])
cbar.ax.tick_params(labelsize=25)

# cbar = plt.colorbar(im,location='top',shrink=0.8)#pad=0.1)
# cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=-40, y=2.9)
# ax=cbar.ax
ax8.text(10**2.3,10**1.5,'$\\Delta F_f$',fontsize=40,rotation=0)

B_position_params=([0,0],[0,0],[0],[1])
loadpath='../data/heatmapdata/Reversible_'+str(reversible)+'/Structure_'+str(lengthparams)+\
        '/cosub_struc_'+str(B_position_params)+'/'+'defaultparams_'+\
        str([CKcat,Ks,Kp,Keq])+'_'+str([KsB,KpB])+\
        '_'+str([Ks_top,Kp_top,Ks_bottom,Kp_bottom])+\
        '_'+str([kout_top,kout_bottom])+'/xvar_'+xvar+\
        '_'+str(xmin)+'_'+str(xmax)+'_'+str(xnum)+'__yvar_'+\
        yvar+'_'+str(ymin)+'_'+str(ymax)+'_'+str(ynum)

loadstring='ff_Btot_'+str(Btot)+'__CKcatB_'+str(VMB)+'__CKcat_bottom_'+str(CKcbottom)+\
    '__CKcat_top_'+str(CKctop)+'__KeqB_'+str(KeqB)+'__Keq_bottom_1.0__Keq_top_1.0__.npy'

ffdeltamat=np.load(loadpath+'/'+loadstring)

im=ax9.pcolormesh(vals_b_in,vals_b_out,-ffdeltamat,vmin=-1,vmax=1,cmap=cmapuse,edgecolors='face')
ax9.set_xscale('log')
ax9.set_yscale('log')

ax9.set_xlabel('${B_{\\rm 0}\\rm ~synthesis~rate}$, ${v^{B_0}_{\\rm in}}$')
ax9.xaxis.set_label_coords(0.5,-0.08)
# ax7.set_ylabel('${v^{B}_{\\rm deg}}$',fontsize=30)
# ax7.yaxis.set_label_coords(-0.1,0.6)
ax9.set_xticks(np.logspace(-4,2,4))
ax9.set_ylabel('${\\rm Co}$-${\\rm substrate ~degredation ~rate}, {v^{B}_{\\rm deg}}$')#,fontsize=30)
ax9.yaxis.set_label_coords(-0.12,0.5)
# ax9.set_xticks(np.logspace(-4,2,4),['','','',''])
ax9.set_yticks(np.logspace(-4,2,4))
ax9.scatter(binoutvals[0],binoutvals[1],marker='o',color='olive',s=150)
ax9.text(10**-5,10**1.3,'$\\rm A(ii)$', fontsize=50)

# cbar = plt.colorbar(im, ticks=np.linspace(-1,1,5))
# #cbar.set_label('$\\Delta F_f$',fontsize=30,rotation=0,labelpad=-40, y=1.15)
# cbar.set_label('$\\Delta F_f$',fontsize=40,rotation=0,labelpad=0, y=0.9)
# # no cosubs
# for j in range(len(V24_vals)):
#     V24=V24_vals[j]
#     print(V24)
#     ff_vec=np.zeros(np.size(kinvals))
#     ff_vec[:]=np.nan
#     for i in range(len(kinvals)):
#         kin=kinvals[i]
        
#         isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams=[2,1,1],B_position_params=([0,0],[0,0],[0],[0]),
#                           params_global=[1.,1.,1.,1.],
#                           params_B=[VMB,1.,1.,KeqB,Btot],
#                           params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
#                           params_inout=[kin,1.,1.],
#                           reversible=rev,Tmax=100000,tpts=10000)
#         if sum(isbuildup)==0:
#             ff_vec[i]=fluxfrac
        
        
#     ax4.plot(kinvals,ff_vec,color=cols[j],label='$V_{24} = '+str(V24)+'$')
# ax4.set_xscale('log')
# # ax4.set_yscale('log')
# ax4.legend(loc=(0.05,0.15),fontsize='large')
# ax4.set_ylim(0,1)
# ax4.set_xlabel('$v_{\\rm in}$',fontsize=30)
# ax4.xaxis.set_label_coords(0.8,-0.05)
# ax4.set_ylabel('$\\rm Flux~ Fraction$',fontsize=40)
# ax4.yaxis.set_label_coords(-0.1,0.4)
# ax4.text(10**-4.6,0.86,'$\\rm E$', fontsize=40)
# ax4.text(10**-4.5,1.75,'$\\rm A$', fontsize=40)



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
                          params_B=[VMB,1.,1.,KeqB,Btot],
                          params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          bsynthdeg=binoutvals,
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
            b_vec[i]=1.-bfrac
        
    ax4.plot(kinvals,ff_vec,color=cols[j],label='$V_{\\rm max}^{24} = '+str(V24)+'$',linewidth=lw)
    inset1.plot(kinvals,b_vec,color=cols[j],linewidth=lw)#,linestyle=':') # add as inset

legfont=20
ax4.set_xscale('log')
# ax5.set_yscale('log')
# ax4.legend(loc=(0.6,0.05),fontsize=legfont)
ax4.set_ylim(0,1)
ax4.set_xlabel('$v_{\\rm in}$',fontsize=35)
ax4.xaxis.set_label_coords(0.85,-0.05)
ax4.text(10**-4.1,0.8,'$\\rm B(iii)$', fontsize=50)
# ax4.set_ylabel('${\\rm Fraction ~of ~upper ~branch ~flux}, F_f$')
# ax4.yaxis.set_label_coords(-0.1,0.5)
ax4.set_yticks(np.linspace(0,1,6),['','','','','',''])

ax4.text(10**-4.1,2.65,'$\\rm B(i)$', fontsize=50)

ax4.text(10**-4.1,2.52,'$v_{\\rm in}$',fontsize=35)
ax4.text(10**-4.1,2.32,'${v^{B_0}_{\\rm in}}$',fontsize=25)
ax4.text(10**-3.4,2.32,'${v^{B}_{\\rm deg}}$',fontsize=25)
ax4.set_xticks(np.logspace(-4,0,5))


inset1.set_xscale('log')
# inset1.set_yscale('log')
inset1.set_ylim(-0.05,1.05)
inset1.set_xticks(np.logspace(-4,0,3))
inset1.set_xlabel('$v_{\\rm in}$')#,fontsize=25)
inset1.xaxis.set_label_coords(0.75,-0.15)
inset1.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
inset1.yaxis.set_label_coords(-0.2,0.4)





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
                          params_B=[VMB,1.,1.,KeqB,Btot],
                          params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          bsynthdeg=binoutvals,
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
            b_vec[i]=1.-bfrac
        
    ax5.plot(kinvals,ff_vec,color=cols[j],label='$V_{\\rm max}^{24}= '+str(V24)+'$',linewidth=lw)
    inset2.plot(kinvals,b_vec,color=cols[j],linewidth=lw)#,linestyle=':') # add as inset
    
ax5.set_xscale('log')
# ax6.set_yscale('log')
# ax6.legend(loc=(0.6,0.05),fontsize='large')
# ax5.legend(loc=(0.6,0.5),fontsize=legfont)
ax5.set_ylim(0,1)
ax5.set_xlabel('$v_{\\rm in}$',fontsize=35)
ax5.xaxis.set_label_coords(0.85,-0.05)
ax5.text(10**-4.1,0.8,'$\\rm C(iii)$', fontsize=50)

# ax5.set_ylabel('${\\rm Fraction ~of ~upper ~branch ~flux}, F_f$')
# ax5.yaxis.set_label_coords(-0.1,0.5)
ax5.set_yticks(np.linspace(0,1,6),['','','','','',''])

ax5.text(10**-4.1,2.65,'$\\rm C(i)$', fontsize=50)

ax5.text(10**-4.1,2.52,'$v_{\\rm in}$',fontsize=35)
ax5.text(10**-4.1,2.32,'${v^{B_0}_{\\rm in}}$',fontsize=25)
ax5.text(10**-3.4,2.32,'${v^{B}_{\\rm deg}}$',fontsize=25)
ax5.set_xticks(np.logspace(-4,0,5))

inset2.set_xscale('log')
# inset2.set_yscale('log')
inset2.set_ylim(-0.05,1.05)
inset2.set_xticks(np.logspace(-4,0,3))
inset2.set_xlabel('$v_{\\rm in}$')#,fontsize=25)
inset2.xaxis.set_label_coords(0.75,-0.15)
inset2.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
inset2.yaxis.set_label_coords(-0.2,0.4)



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
                          params_B=[VMB,1.,1.,KeqB,Btot],
                          params_branch=[1.,1.,1.,1.,V24,1.,1.,1.],
                          params_inout=[kin,1.,1.],
                          bsynthdeg=binoutvals,
                          reversible=rev,Tmax=100000,tpts=10000)
        if sum(isbuildup)==0:
            ff_vec[i]=1.-fluxfrac
            b_vec[i]=1.-bfrac
        if cols[j]=='black':
            ffvec_save=ff_vec
            
    ax6.plot(kinvals,ff_vec,color=cols[j],label='$V_{\\rm max}^{24} = '+str(V24)+'$',linewidth=lw)
    inset3.plot(kinvals,b_vec,color=cols[j],linewidth=lw)#,linestyle=':') # add as inset
  
ax6.set_xscale('log')
# ax6.set_yscale('log')
# ax6.legend(loc=(0.6,0.05),fontsize='large')
# ax6.legend(loc=(0.6,0.5),fontsize=legfont)
ax6.set_ylim(0,1)
ax6.set_xlabel('$v_{\\rm in}$',fontsize=35)
ax6.xaxis.set_label_coords(0.85,-0.05)
ax6.text(10**-4.1,0.8,'$\\rm A(iii)$', fontsize=50)
ax6.set_ylabel('${\\rm Fraction ~of ~upper ~branch ~flux}, F_f$')
ax6.yaxis.set_label_coords(-0.1,0.5)
ax6.set_yticks(np.linspace(0,1,6))
ax6.text(10**-4.1,2.65,'$\\rm A(i)$', fontsize=50)

ax6.text(10**-4.1,2.52,'$v_{\\rm in}$',fontsize=35)
ax6.text(10**-4.1,2.32,'${v^{B_0}_{\\rm in}}$',fontsize=25)
ax6.text(10**-3.4,2.32,'${v^{B}_{\\rm deg}}$',fontsize=25)
ax6.set_xticks(np.logspace(-4,0,5))

ax6.text(10**-4,0.4,'$\\frac{}{V_{\\rm max}^{24}} =$',fontsize=40)
ax6.text(10**-4,0.46,'$V_{\\rm max}^{23} $',fontsize=29)
ax6.text(10**-2.9,0.47,'$10$',fontsize=25,horizontalalignment='right',color='b')
ax6.text(10**-2.9,0.41,'$1$',fontsize=25,horizontalalignment='right',color='olive')
ax6.text(10**-2.9,0.35,'$0.1$',fontsize=25,horizontalalignment='right',color='r')

inset3.set_xscale('log')
# inset2.set_yscale('log')
inset3.set_ylim(-0.05,1.05)
inset3.set_xticks(np.logspace(-4,0,3))
inset3.set_xlabel('$v_{\\rm in}$')#,fontsize=25/)
inset3.xaxis.set_label_coords(0.7,-0.15)
inset3.set_ylabel('$\\frac{B_1}{B_T}$',rotation=0,fontsize=35)
inset3.yaxis.set_label_coords(-0.2,0.4)

# min_y=np.nanmin(ffvec_save)
# print(min_y)
# ax8.plot([10**-4,10**0],[min_y,min_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')

# max_y=np.nanmax(ffvec_save)
# print(max_y)
# ax8.plot([10**-4,10**-2],[max_y,max_y],color=[0.5,0.5,0.5],linewidth=3,linestyle='--')




# ax9.set_ylim(0,1)
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

if savefig==1:
    plt.savefig('Figure_S2.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S2.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S2.svg',bbox_inches='tight',format='svg')
    
plt.show()
