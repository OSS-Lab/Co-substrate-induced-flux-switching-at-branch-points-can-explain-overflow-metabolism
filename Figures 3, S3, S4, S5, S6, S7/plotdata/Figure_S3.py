#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 12:03:33 2025

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
from matplotlib.colors import LogNorm

from matplotlib.patches import Rectangle
import make_heatmaps as mh

savefig=0         # change to 1 to save generated figure
nkin=101
# nkin=21

cmapuse='Spectral_r'


params = {'legend.fontsize': 18,
          'figure.figsize': (12, 9),
          'axes.labelsize': 30,
          'axes.titlesize': 30,
          'xtick.labelsize':30,
          'ytick.labelsize':30}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)

plt.figure()

inset1=plt.axes([0.075,0.7,0.48,0.48]) 
motifnone=img.imread('figS3_none.png')
plt.imshow(motifnone)
plt.axis('off')

inset2=plt.axes([0.74,0.7,0.48,0.48]) 
motifnone=img.imread('figS3_none.png')
plt.imshow(motifnone)
plt.axis('off')

inset3=plt.axes([0.075,-0.5,0.48,0.48]) 
motifnone=img.imread('figS3_minimal.png')
plt.imshow(motifnone)
plt.axis('off')

inset4=plt.axes([0.74,-0.5,0.48,0.48]) 
motifnone=img.imread('figS3_minimal.png')
plt.imshow(motifnone)
plt.axis('off')

# ax11=plt.axes([1.5,-1,0.5,0.5])
ax10=plt.axes([-0.05,-1.3,1.35,2.4])   
ax10.set_ylim(0,1)     
ax10.set_xlim(0,1)     
plt.axis('off')

ax1=plt.axes([0.1,0.05,0.48,0.5])

Mvalues = np.logspace(-2,2,101)

K=1
V23=10
V24=1

f_simple = np.zeros(shape=(np.size(Mvalues)))
v23 = np.zeros(shape=(np.size(Mvalues)))
v24 = np.zeros(shape=(np.size(Mvalues)))

for idM,M in enumerate(Mvalues):
    v23[idM]=V23*M/(K+M)
    v24[idM]=V24*M/(K+M)
    f_simple[idM]=V23/(V23+V24)
    # print(V23/(V23+V24))

ax1.plot(Mvalues,v23,'k',linewidth=3,label='$v_{23}$')# (mmol/min)')
ax1.plot(Mvalues,v24,'r',linewidth=3,label='$v_{24}$')# (mmol/min)')
ax1.plot(Mvalues,f_simple,'b',linewidth=3,label='$F_f$')
# im = ax1.pcolor(Bvalues,Mvalues,f_simple, cmap=cmapuse,vmin=0,vmax=1)
# print(np.min(f_simple))
# print(np.max(f_simple))
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel('$\\rm M_2~ Concentration$')
ax1.legend(loc='lower right')
# ax1.set_ylabel('$\\rm M_2~ Concentration$')

# cbar = plt.colorbar(im, orientation='vertical')
# cbar.set_label('upper to lower branch flux ratio')
#############################


ax2=plt.axes([0.75,0.05,0.48,0.5])

Mvalues = np.logspace(-2,2,101)

f_simple_kdiff = np.zeros(shape=np.size(Mvalues))
v23_simple_kdiff = np.zeros(shape=np.size(Mvalues))
v24_simple_kdiff = np.zeros(shape=np.size(Mvalues))
K24=1
K23=100
V23=10
V24=1
    
for idM,M in enumerate(Mvalues):
    v23_simple_kdiff[idM]=V23*M/(K23+M)
    v24_simple_kdiff[idM]=V24*M/(K24+M)
    
    ratbit=(K23+M)/(K24+M)
    frat=V23/(V23+V24*ratbit)
    # f_simple_kdiff[idM,idB]=flux23/(flux23+flux24)
    f_simple_kdiff[idM]=frat
        
        
      
    
 

ax2.plot(Mvalues,v23_simple_kdiff,'k',linewidth=3,label='$v_{23}$')# (mmol/min)')
ax2.plot(Mvalues,v24_simple_kdiff,'r',linewidth=3,label='$v_{24}$')# (mmol/min)')
ax2.plot(Mvalues,f_simple_kdiff,'b',linewidth=3,label='$F_f$')
# im = ax1.pcolor(Bvalues,Mvalues,f_simple, cmap=cmapuse,vmin=0,vmax=1)
# print(np.min(f_simple))
# print(np.max(f_simple))
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('$\\rm M_2~ Concentration$')
ax2.legend(loc='lower right')
# lengthparams=[2,1,1]
# B_position_params=([0,0],[0,0],[0],[0])

# Btot=1
# Vmb=1
# keqB=1

# B0conc,M2conc=mh.return_m2_b0(lengthparams,B_position_params,
#             params_global=[1.,1.,1.,1.],params_B=[Vmb,1.,1.,keqB,Btot],
#             params_branch=[V23,K23,1.,1.,
#                            V24,K24,1.,1.],
#             kouts=[1.,1.],#params_inout=[0.1,1.,1.],
#             kin_min=-8,kin_max=1,num_kin=41,
#             reversible=False,Tmax=1000000,tpts=100000)

# print(B0conc)
# print(M2conc)

# ax2.plot(B0conc,M2conc,'-k',linewidth=2)
# #############################



ax3=plt.axes([0.1,-1.15,0.4,0.5])


Bvalues = np.logspace(-5,3,100)
Mvalues = np.logspace(-5,3,100)

f_simple_B = np.zeros(shape=(np.size(Mvalues), np.size(Bvalues)))
K=1
V=1
for idB,B in enumerate(Bvalues):
    
    for idM,M in enumerate(Mvalues):
        flux23=V*M/(K+M)
        flux24=V*M*B/(K+M*B)
        
        ratbit=(K+M)/(K+M*B)
        frat=1/(1+B*ratbit)
        # f_simple_B[idM,idB]=flux23/(flux23+flux24)
        f_simple_B[idM,idB]=frat
        
        
      
    
im=ax3.pcolormesh(Bvalues,Mvalues,f_simple_B, cmap=cmapuse,vmin=0,vmax=1,
                  edgecolors='face')
# print(np.min(f_simple_B))
# print(np.max(f_simple_B))  
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel('$\\rm B_0~ Concentration$')
ax3.set_ylabel('$\\rm M_2~ Concentration$')
# cbar = plt.colorbar(im, orientation='vertical',label='$F_f$')
cax = ax3.inset_axes([1.05, 0, 0.04, 1])
cbar=plt.colorbar(im, cax=cax, orientation='vertical',ticks=[-1, -0.5,0,0.5, 1])
cbar.set_label('$F_f$', labelpad=-10, y=0.8, rotation=0,fontsize=40)
# cbar.ax.tick_params(labelsize=25)

lengthparams=[2,1,1]
B_position_params=([0,0],[0,0],[0],[1])

Btot=10
Vmb=0.01
keqB=1

B0conc,M2conc=mh.return_m2_b0(lengthparams,B_position_params,
            params_global=[1.,1.,1.,1.],params_B=[Vmb,1.,1.,keqB,Btot],
            params_branch=[V,K,1.,1.,
                           V,K,1.,1.],
            kouts=[1.,1.],#params_inout=[0.1,1.,1.],
            kinvals=np.logspace(-4,0,nkin),
            reversible=False,Tmax=100000,tpts=10000)

# print(B0conc)
# print(M2conc)
bconc=[B0conc[0],B0conc[-1]]
mconc=[M2conc[0],M2conc[-1]]
ax3.plot(B0conc,M2conc,'-k',linewidth=4)
ax3.plot(bconc,mconc,'ok',markersize=10)




#################################

ax4=plt.axes([0.75,-1.15,0.4,0.5])

Bvalues = np.logspace(-5,3,100)
Mvalues = np.logspace(-5,3,100)

f_simple_B_Vdiff = np.zeros(shape=(np.size(Mvalues), np.size(Bvalues)))
K=1
V23=0.1
V24=1
for idB,B in enumerate(Bvalues):
    
    for idM,M in enumerate(Mvalues):
        flux23=V23*M/(K+M)
        flux24=V24*M*B/(K+M*B)
        
        ratbit=(K+M)/(K+M*B)
        frat=V23/(V23+V24*B*ratbit)
        # f_simple_B[idM,idB]=flux23/(flux23+flux24)
        f_simple_B_Vdiff[idM,idB]=frat
        
        
      
    
im=ax4.pcolormesh(Bvalues,Mvalues,f_simple_B_Vdiff, cmap=cmapuse,vmin=0,vmax=1,
                  edgecolors='face')
# print(np.min(f_simple_B))
# print(np.max(f_simple_B))  
ax4.set_xscale('log')
ax4.set_yscale('log')
ax4.set_xlabel('$\\rm B_0~ Concentration$')
ax4.set_ylabel('$\\rm M_2~ Concentration$')
cax = ax4.inset_axes([1.05, 0, 0.04, 1])
cbar=plt.colorbar(im, cax=cax, orientation='vertical',ticks=[-1, -0.5,0,0.5, 1])
cbar.set_label('$F_f$', labelpad=-10, y=0.8, rotation=0,fontsize=40)

lengthparams=[2,1,1]
B_position_params=([0,0],[0,0],[0],[1])

Btot=100
Vmb=10**-2
keqB=1

B0conc,M2conc=mh.return_m2_b0(lengthparams,B_position_params,
            params_global=[1.,1.,1.,1.],params_B=[Vmb,1.,1.,keqB,Btot],
            params_branch=[V23,K,1.,1.,
                           V24,K,1.,1.],
            kouts=[1.,1.],#params_inout=[0.1,1.,1.],
            kinvals=np.logspace(-2.9,-0.95,nkin),
            reversible=False,Tmax=100000,tpts=10000)


# 
ax4.plot(B0conc,M2conc,'-k',linewidth=4,label='High $B_{\\rm T}$, Low $V_{\\rm max}^B$')
bconc=[B0conc[0],B0conc[-1]]
mconc=[M2conc[0],M2conc[-1]]
ax4.plot(bconc,mconc,'ok',markersize=10)
# ax4.set_ylim(np.min(Mvalues),np.max(Mvalues))

# ax11.plot(kinv,B0conc,'r')
# ax11.plot(kinv,M2conc,'b')
# ax11.set_xscale('log')
# ax11.set_yscale('log')

lengthparams=[2,1,1]
B_position_params=([0,0],[0,0],[0],[1])

Btot=0.01
Vmb=10**2
keqB=1

B0conc,M2conc=mh.return_m2_b0(lengthparams,B_position_params,
            params_global=[1.,1.,1.,1.],params_B=[Vmb,1.,1.,keqB,Btot],
            params_branch=[V23,K,1.,1.,
                           V24,K,1.,1.],
            kouts=[1.,1.],#params_inout=[0.1,1.,1.],
            kinvals=np.logspace(-5,-0.2,nkin),
            reversible=False,Tmax=100000,tpts=10000)

# print(B0conc)
# print(M2conc)

print(np.min(B0conc))
print(np.max(B0conc))
print(np.min(M2conc))
print(np.max(M2conc))

ax4.plot(B0conc,M2conc,'--k',linewidth=4,label='Low $B_{\\rm T}$, High $V_{\\rm max}^B$')
bconc=[B0conc[0],B0conc[-1]]
mconc=[M2conc[0],M2conc[-1]]
ax4.plot(bconc,mconc,'ok',markersize=10)
ax4.legend(loc='upper right',fontsize=18)

ax1.text(10**-2.3,10**3.8,'$v_{\\rm in}$',fontsize=30)
ax1.text(10**0.32,10**4.,'$v_{23}$',fontsize=30,rotation=40)
ax1.text(10**0.35,10**3.2,'$v_{24}$',fontsize=30,rotation=-30)

ax1.text(10**-2.5,10**2.4,'$v_{23} = \\frac{V_{\\rm max}^{23}m_2}{K_{\\rm M}+m_2},$',
         fontsize=30)
ax1.text(10**0,10**2.4,'$v_{24} = \\frac{V_{\\rm max}^{24}m_2}{K_{\\rm M}+m_2}$',
         fontsize=30)

ax1.text(10**-2.5,10**1.5,'$F_f = \\frac{v_{23}}{v_{23}+v_{24}} = \\frac{V_{\\rm max}^{23}}{V_{\\rm max}^{23}+V_{\\rm max}^{24}}$',
         fontsize=35)

ax1.text(10**-3.1,10**4.2,'$\\rm A$',fontsize=60)


ax2.text(10**-2.15,10**4.15,'$v_{\\rm in}$',fontsize=30)
ax2.text(10**0.45,10**4.4,'$v_{23}$',fontsize=30,rotation=40)
ax2.text(10**0.45,10**3.4,'$v_{24}$',fontsize=30,rotation=-30)

ax2.text(10**-2.5,10**2.4,'$v_{23} = \\frac{V_{\\rm max}^{23}m_2}{K_{\\rm M}^{23}+m_2},$',
         fontsize=30)
ax2.text(10**0,10**2.4,'$v_{24} = \\frac{V_{\\rm max}^{24}m_2}{K_{\\rm M}^{24}+m_2}$',
         fontsize=30)

ax2.text(10**-2.5,10**1.6,'$F_f = \\frac{v_{23}}{v_{23}+v_{24}} = \\frac{V_{\\rm max}^{23}}{V_{\\rm max}^{23}+V_{\\rm max}^{24}\\left(\\frac{K_{\\rm M}^{23}+m_2}{K_{\\rm M}^{24}+m_2}\\right)}$',
         fontsize=35)

ax2.text(10**-3.,10**4.6,'$\\rm B$',fontsize=60)


ax3.text(10**-5.2,10**9.5,'$v_{\\rm in}$',fontsize=30)
ax3.text(10**0.6,10**10,'$v_{23}$',fontsize=30,rotation=40)
ax3.text(10**1,10**8.7,'$v_{24}$',fontsize=30,rotation=-30)

ax3.text(10**-5,10**6.2,'$v_{23} = \\frac{V_{\\rm max}m_2}{K_{\\rm M}+m_2},$',
         fontsize=30)
ax3.text(10**0,10**6.2,'$v_{24} = \\frac{V_{\\rm max}m_2b_0}{K_{\\rm M}+m_2b_0}$',
         fontsize=30)
ax3.text(10**-5,10**4.5,'$F_f = \\frac{v_{23}}{v_{23}+v_{24}} = \\frac{1}{1+b_0\\left( \\frac{K_{\\rm M}+m_2}{K_{\\rm M}+m_2b_0}\\right)}$',fontsize=35)

ax3.text(10**-7,10**10.5,'$\\rm C$',fontsize=60)



ax4.text(10**-4.9,10**9.5,'$v_{\\rm in}$',fontsize=30)
ax4.text(10**0.8,10**10,'$v_{23}$',fontsize=30,rotation=40)
ax4.text(10**1.4,10**8.7,'$v_{24}$',fontsize=30,rotation=-30)

ax4.text(10**-5,10**6.1,'$v_{23} = \\frac{V_{\\rm max}^{23}m_2}{K_{\\rm M}+m_2},$',
         fontsize=30)
ax4.text(10**0.5,10**6.1,'$v_{24} = \\frac{V_{\\rm max}^{24}m_2b_0}{K_{\\rm M}+m_2b_0}$',
         fontsize=30)
ax4.text(10**-5.5,10**4.5,'$F_f = \\frac{v_{23}}{v_{23}+v_{24}} = \\frac{V_{\\rm max}^{23}}{V_{\\rm max}^{23}+V_{\\rm max}^{24}b_0\\left( \\frac{K_{\\rm M}+m_2}{K_{\\rm M}+m_2b_0}\\right)}$',fontsize=35)

ax4.text(10**-6.9,10**10.5,'$\\rm D$',fontsize=60)

rect = Rectangle((0.02, 0.515), width=0.47, height=0.485, edgecolor='k',facecolor='none',linewidth=2)
ax10.add_patch(rect)

rect = Rectangle((0.51, 0.515), width=0.47, height=0.485, edgecolor='k',facecolor='none',linewidth=2)
ax10.add_patch(rect)

rect = Rectangle((0.02, 0.02), width=0.47, height=0.48, edgecolor='k',facecolor='none',linewidth=2)
ax10.add_patch(rect)

rect = Rectangle((0.51, 0.02), width=0.47, height=0.48, edgecolor='k',facecolor='none',linewidth=2)
ax10.add_patch(rect)



ax33=plt.axes([0.1,-1.15,0.4,0.5])
ax33.set_ylim(0,1)
ax33.set_xlim(0,1)
ax33.axis('off')


ax33.arrow(0.6, 0.5, -0.2, 0.2,length_includes_head=True,
          width=0.005,head_width=0.03,head_length=0.05,
          color='k')#, **kwargs)

ax33.text(0.5,0.65,'Increasing $v_{\\rm in}$',fontsize=30)


ax44=plt.axes([0.75,-1.15,0.4,0.5])
ax44.set_ylim(0,1)
ax44.set_xlim(0,1)
ax44.axis('off')


ax44.arrow(0.75, 0.35, -0.25, 0.25,length_includes_head=True,
          width=0.005,head_width=0.03,head_length=0.05,
          color='k')#, **kwargs)

ax44.arrow(0.25, 0.2, 0, 0.25,length_includes_head=True,
          width=0.005,head_width=0.03,head_length=0.05,
          color='k')#, **kwargs)


if savefig==1:
    plt.savefig('Figure_S3.png',bbox_inches='tight',format='png')
    plt.savefig('Figure_S3.eps',bbox_inches='tight',format='eps')
    plt.savefig('Figure_S3.svg',bbox_inches='tight',format='svg')
    plt.savefig('Figure_S3.pdf',bbox_inches='tight',format='pdf')

plt.show()
