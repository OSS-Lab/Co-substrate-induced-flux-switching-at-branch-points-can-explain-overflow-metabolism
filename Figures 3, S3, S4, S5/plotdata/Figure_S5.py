#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 22:02:47 2024

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


savefig=0               # change to 1 to save generated figure
num_kin=201
# num_kin=51

# mfd.FIGURE2_2_ff_vs_kin_vary24()


xvar='Btot'
yvar='CKcatB'

# KeqB_val=1.

KeqB_val=1.
CKcat_top_val=0.1
CKcat_bottom_val=1.

reversible=False
rev=False

kin_min=-5
kin_max=0


kinvals=np.logspace(kin_min,kin_max,num_kin)


##### THIS IS FIGURE S1

cmapuse='RdBu_r'

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (12, 9),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

plt.figure()


# BT_dots_1=[10**-1,10**-2,10,10**-1.5]
# VMB_dots_1=[0.01,10**2,0.01,10**-1]

BT_dots_1=[10**-1,10**-2,10]
VMB_dots_1=[0.01,10**2,0.01]
cols_1=['grey','coral','violet']

BT_dots_2=[0.01,10]
VMB_dots_2=[100,0.01]
cols_2=['coral','violet']

BT_dots_3=[0.01,10**-2.75,10]
VMB_dots_3=[10**-2,10**0,0.01]
cols_3=['grey','coral','violet']

BT_dots_4=[0.01,10**-2,10]
VMB_dots_4=[10**-2,10**-0.75,0.01]
cols_4=['grey','coral','violet']


ax1=plt.axes([0.07,0.5,0.55,0.5])
plt.axis('off')
motifnone=img.imread('Fig_S5_updown_line.png')
plt.imshow(motifnone)
ax2=plt.axes([0.7,0.6,0.3,0.3])
ax9=plt.axes([1.13,0.6,0.3,0.3])


ax3=plt.axes([0.07,0.12,0.55,0.5])
plt.axis('off')
motifnone=img.imread('Fig_S5_updown_stronger.png')
plt.imshow(motifnone)
ax4=plt.axes([0.7,0.22,0.3,0.3])
ax10=plt.axes([1.13,0.22,0.3,0.3])

ax5=plt.axes([0.07,-0.26,0.55,0.5])
plt.axis('off')
motifnone=img.imread('fig_S5_light_triangle.png')
plt.imshow(motifnone)
ax6=plt.axes([0.7,-0.16,0.3,0.3])
ax11=plt.axes([1.13,-0.16,0.3,0.3])

ax7=plt.axes([0.07,-0.64,0.55,0.5])
plt.axis('off')
motifnone=img.imread('fig_S5_dark_triangle.png')
plt.imshow(motifnone)
ax8=plt.axes([0.7,-0.54,0.3,0.3])
ax12=plt.axes([1.13,-0.54,0.3,0.3])


lengthparams_set=([2,2,2],
                  [2,2,2],
				  [3,1,1],
                  [3,1,1])

B_position_params_set=(([0,0],[0,0],[0,1],[1,0]),
                       ([1,0],[0,0],[0,1],[1,0]),
					   ([0,0,0],[0,0,0],[1],[1]),
                       ([0,1,0],[0,0,0],[1],[1]))


 
        
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

count=0

for idstruc,struc in enumerate(lengthparams_set):
    

    lengthparams=struc
    B_position_params=B_position_params_set[idstruc]

    print(lengthparams)
    print(B_position_params)  
    
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
    
    print(loadpath)
    if count==0:
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
        ax2.yaxis.set_label_coords(-0.12,0.45)
        ax2.set_xticks(np.logspace(-3,3,4))
        ax2.set_yticks(np.logspace(-3,3,4))
        
        for i in range(len(BT_dots_1)):
            ax2.plot(BT_dots_1[i],VMB_dots_1[i],'o',markersize=12,
                     markeredgewidth=2,markeredgecolor='k',markerfacecolor=cols_1[i])
            
            
        for idb,Btot in enumerate(BT_dots_1):
            VMB=VMB_dots_1[idb]
            col=cols_1[idb]
            #print(Btot,VMB,col)
            ff_vec=np.zeros(np.size(kinvals))
            ff_vec[:]=np.nan
            for i in range(len(kinvals)):
                kin=kinvals[i]
                
                isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                               B_position_params,
                                  params_global=[1.,1.,1.,1.],
                                  params_B=[VMB,1.,1.,KeqB_val,Btot],
                                  params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                                  params_inout=[kin,1.,1.],
                                  reversible=rev,Tmax=100000,tpts=10000)
                if sum(isbuildup)==0:
                    ff_vec[i]=1.-fluxfrac
                
                
            ax9.plot(kinvals,ff_vec,color=col,linewidth=3)
            #,label='$B_T = 10^{'+str(np.log10(Btot))+'}, V_M^B = 10^{'+str(np.log10(VMB))+'}$')

        ax9.set_xlabel('$v_{\\rm in}$',fontsize=30)    
        ax9.set_ylabel('$\\rm Flux~ Fraction$',rotation=90,fontsize=30)  
        ax9.set_ylim(0,1.05) 
        ax9.set_xlim(10**-5.1,10**0)
        ax9.set_xscale('log')
       # ax9.legend(loc='lower right',fontsize=15)
        ax9.xaxis.set_label_coords(0.9,-0.05)
        ax9.yaxis.set_label_coords(-0.15,0.40)
        
    if count==1:
        im=ax4.pcolormesh(np.logspace(xmin,xmax,xnum),
                                np.logspace(ymin,ymax,ynum),
                                -ffdelta_mat,edgecolors='face',
                                vmin=-1,vmax=1,cmap=cmapuse)
                #edgecolors{'none', None, 'face', color, color sequence}, optional
        ax4.set_xscale('log')
        ax4.set_yscale('log')
        ax4.set_xlabel(xlabstr,fontsize=25)    
        ax4.set_ylabel(ylabstr,fontsize=25,rotation=0)   
        ax4.xaxis.set_label_coords(0.8,-0.08)
        ax4.yaxis.set_label_coords(-0.12,0.45)
        ax4.set_xticks(np.logspace(-3,3,4))
        ax4.set_yticks(np.logspace(-3,3,4))
        
        for i in range(len(BT_dots_2)):
            ax4.plot(BT_dots_2[i],VMB_dots_2[i],'o',markersize=12,
                     markeredgewidth=2,markeredgecolor='k',markerfacecolor=cols_2[i])
        
        for idb,Btot in enumerate(BT_dots_2):
            VMB=VMB_dots_2[idb]
            col=cols_2[idb]
            #print(Btot,VMB,col)
            ff_vec=np.zeros(np.size(kinvals))
            ff_vec[:]=np.nan
            for i in range(len(kinvals)):
                kin=kinvals[i]
                
                isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                               B_position_params,
                                  params_global=[1.,1.,1.,1.],
                                  params_B=[VMB,1.,1.,KeqB_val,Btot],
                                  params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                                  params_inout=[kin,1.,1.],
                                  reversible=rev,Tmax=100000,tpts=10000)
                if sum(isbuildup)==0:
                    ff_vec[i]=1.-fluxfrac
                
                
            ax10.plot(kinvals,ff_vec,color=col,linewidth=3)
            #,label='$B_T = 10^{'+str(np.log10(Btot))+'}, V_M^B = 10^{'+str(np.log10(VMB))+'}$')

        ax10.set_xlabel('$v_{\\rm in}$',fontsize=30)    
        ax10.set_ylabel('$\\rm Flux~ Fraction$',rotation=90,fontsize=30)  
        ax10.set_ylim(0,1.05) 
        ax10.set_xlim(10**-5.1,10**0)
        ax10.set_xscale('log')
       # ax9.legend(loc='lower right',fontsize=15)
        ax10.xaxis.set_label_coords(0.9,-0.05)
        ax10.yaxis.set_label_coords(-0.15,0.40)
        
            
    if count==2:
        im=ax6.pcolormesh(np.logspace(xmin,xmax,xnum),
                                np.logspace(ymin,ymax,ynum),
                                -ffdelta_mat,edgecolors='face',
                                vmin=-1,vmax=1,cmap=cmapuse)
                #edgecolors{'none', None, 'face', color, color sequence}, optional
        ax6.set_xscale('log')
        ax6.set_yscale('log')
        ax6.set_xlabel(xlabstr,fontsize=25)    
        ax6.set_ylabel(ylabstr,fontsize=25,rotation=0)   
        ax6.xaxis.set_label_coords(0.8,-0.08)
        ax6.yaxis.set_label_coords(-0.12,0.45)
        ax6.set_xticks(np.logspace(-3,3,4))
        ax6.set_yticks(np.logspace(-3,3,4))
        
        for i in range(len(BT_dots_3)):
            ax6.plot(BT_dots_3[i],VMB_dots_3[i],'o',markersize=12,
                     markeredgewidth=2,markeredgecolor='k',markerfacecolor=cols_3[i])
            
        for idb,Btot in enumerate(BT_dots_3):
            VMB=VMB_dots_3[idb]
            col=cols_3[idb]
            #print(Btot,VMB,col)
            ff_vec=np.zeros(np.size(kinvals))
            ff_vec[:]=np.nan
            for i in range(len(kinvals)):
                kin=kinvals[i]
                
                isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                               B_position_params,
                                  params_global=[1.,1.,1.,1.],
                                  params_B=[VMB,1.,1.,KeqB_val,Btot],
                                  params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                                  params_inout=[kin,1.,1.],
                                  reversible=rev,Tmax=100000,tpts=10000)
                if sum(isbuildup)==0:
                    ff_vec[i]=1.-fluxfrac
                
                
            ax11.plot(kinvals,ff_vec,color=col,linewidth=3)
            #,label='$B_T = 10^{'+str(np.log10(Btot))+'}, V_M^B = 10^{'+str(np.log10(VMB))+'}$')

        ax11.set_xlabel('$v_{\\rm in}$',fontsize=30)    
        ax11.set_ylabel('$\\rm Flux~ Fraction$',rotation=90,fontsize=30)  
        ax11.set_ylim(0,1.05) 
        ax11.set_xlim(10**-5.1,10**0)
        ax11.set_xscale('log')
       # ax9.legend(loc='lower right',fontsize=15)
        ax11.xaxis.set_label_coords(0.9,-0.05)
        ax11.yaxis.set_label_coords(-0.15,0.40)
        
    if count==3:
        im=ax8.pcolormesh(np.logspace(xmin,xmax,xnum),
                                np.logspace(ymin,ymax,ynum),
                                -ffdelta_mat,edgecolors='face',
                                vmin=-1,vmax=1,cmap=cmapuse)
                #edgecolors{'none', None, 'face', color, color sequence}, optional
        ax8.set_xscale('log')
        ax8.set_yscale('log')
        ax8.set_xlabel(xlabstr,fontsize=25)    
        ax8.set_ylabel(ylabstr,fontsize=25,rotation=0)   
        ax8.xaxis.set_label_coords(0.8,-0.08)
        ax8.yaxis.set_label_coords(-0.12,0.45)
        ax8.set_xticks(np.logspace(-3,3,4))
        ax8.set_yticks(np.logspace(-3,3,4))
        
        for i in range(len(BT_dots_4)):
            ax8.plot(BT_dots_4[i],VMB_dots_4[i],'o',markersize=12,
                     markeredgewidth=2,markeredgecolor='k',markerfacecolor=cols_4[i])
            
        for idb,Btot in enumerate(BT_dots_4):
            VMB=VMB_dots_4[idb]
            col=cols_4[idb]
            #print(Btot,VMB,col)
            ff_vec=np.zeros(np.size(kinvals))
            ff_vec[:]=np.nan
            for i in range(len(kinvals)):
                kin=kinvals[i]
                
                isbuildup, fluxfrac,bfrac, flux_rat,brat=mfd.get_fluxfractions(lengthparams,
                                                                               B_position_params,
                                  params_global=[1.,1.,1.,1.],
                                  params_B=[VMB,1.,1.,KeqB_val,Btot],
                                  params_branch=[CKcat_top_val,1.,1.,1.,CKcat_bottom_val,1.,1.,1.],
                                  params_inout=[kin,1.,1.],
                                  reversible=rev,Tmax=100000,tpts=10000)
                if sum(isbuildup)==0:
                    ff_vec[i]=1.-fluxfrac
                
                
            ax12.plot(kinvals,ff_vec,color=col,linewidth=3)
            #,label='$B_T = 10^{'+str(np.log10(Btot))+'}, V_M^B = 10^{'+str(np.log10(VMB))+'}$')

        ax12.set_xlabel('$v_{\\rm in}$',fontsize=30)    
        ax12.set_ylabel('$\\rm Flux~ Fraction$',rotation=90,fontsize=30)  
        ax12.set_ylim(0,1.05) 
        ax12.set_xlim(10**-5.1,10**0)
        ax12.set_xscale('log')
       # ax9.legend(loc='lower right',fontsize=15)
        ax12.xaxis.set_label_coords(0.9,-0.05)
        ax12.yaxis.set_label_coords(-0.15,0.40)
        


    cbar = plt.colorbar(im)
    #cbar.set_label('$\\Delta F_f$',fontsize=30,rotation=0,labelpad=-40, y=1.15)
    cbar.set_label('$\\Delta F_f$',fontsize=30,rotation=0,labelpad=0, y=0.9)
    count+=1


 


if savefig==1:
    plt.savefig('Figure_S5.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S5.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S5.svg',bbox_inches='tight',format='svg')
    
plt.show()