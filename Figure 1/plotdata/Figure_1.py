#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 12:58:11 2025

@author: robert
"""



import numpy as np

import matplotlib.pyplot as plt

import matplotlib.pylab as pylab
# import make_figure_data as mfd
# import matplotlib.image as img
import matplotlib as mpl
import matplotlib.image as img

savefig=0

params = {'legend.fontsize': 30,
          'figure.figsize': (12, 9),
          'axes.labelsize': 'x-large',
          'axes.titlesize': 20,
          'xtick.labelsize':20,
          'ytick.labelsize':20,
          'legend.title_fontsize':15}
pylab.rcParams.update(params)
mpl.rc('text', usetex = True)
mpl.rc('text.latex', preamble=r'\usepackage{amsmath}')




plt.figure()

ax8=plt.axes([0.62,-0.065,0.55,1.034])
motifnone=img.imread('bg.png')
plt.imshow(motifnone)
plt.axis('off')
# rect = patches.FancyBboxPatch((0,0),
#                               boxstyle="round, pad=0, rounding_size=0.0",
#                               width=1,
#                               height=1,
#                               fill=False,linewidth=3)
# ax8.add_patch(rect)

ax1=plt.axes([0.03,0.55,0.62,0.5])
motifnone=img.imread('lower_glycolysis.png')
plt.imshow(motifnone)
plt.axis('off')

ax3=plt.axes([0.03,0.2,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('upper_glycolysis.png')
plt.imshow(motifbranch)

ax3=plt.axes([0.03,-0.15,0.62,0.5])
plt.axis('off')
motifbranch=img.imread('ammonia_assimilation.png')
plt.imshow(motifbranch)

# ax4=plt.axes([0.7,0.15,0.3,0.4])
# xlabs_y= ('NADH','ATP','NADPH')
# counts_y=np.array([60,189,66])
# width=0.6
# p=ax4.bar(xlabs_y,counts_y,width, color='tab:olive')
# ax4.set_title('Yeast',fontsize=30)
# ax4.bar_label(p,label_type='center',fontsize=20)

# ax5=plt.axes([1.,0.15,0.3,0.4])
# xlabs_e= ('NADH','ATP','NADPH')
# counts_e=np.array([60,112,51])
# p=ax5.bar(xlabs_e,counts_e,width, color='firebrick')
# ax5.set_yticks([])
# ax5.set_title('E. coli',fontsize=30)
# ax5.bar_label(p,label_type='center',fontsize=20)

def fmt(x):
    return '{:.0f}'.format(total*x/100)

def total(x):
    return np.sum(x)

ax4=plt.axes([0.62,0.62,0.3,0.3])
xlabs_cosub_y= ('NAD(H)','ATP','NADP(H)','Multiple')
# sizes_cyto=
sizes_cosub=np.array([10,83,16,29])
p, tx, autotexts=ax4.pie(sizes_cosub, autopct="",radius=0.7,
        colors=['deepskyblue', 'thistle', 'mediumseagreen','darkkhaki'],startangle=90,
        textprops={'size': 20})
for i, a in enumerate(autotexts):
    a.set_text("{}".format(sizes_cosub[i]))
    a.set_fontsize(20)
# ax4.text(-0.8,1.1,'Yeast Cytosol Branch Points',fontsize=30)
ax4.legend(p, xlabs_cosub_y,
          title='\\underline{Co-substrate Usage}',
          loc=(0.12,-0.5), fontsize=15)
#0.35
ax6=plt.axes([0.77,0.25,0.5,0.5])
xlabs_metabs_y=('Non-Branch','Branch with CR','Branch~Without CR')
# xlabs_metabs_y=('Non-Branch','Branch with \n Co-substrate Utilising \n Reactions (CR)','Branch~Without\n Co-substrate Usage')
sizes_metabs_y=np.array([612,138,130])
p, tx, autotexts=ax6.pie(sizes_metabs_y, autopct="",radius=0.7,
        colors=['lightgrey','salmon','wheat'],startangle=230,explode=(0,0.1,0),
        textprops={'size': 20})
for i, a in enumerate(autotexts):
    a.set_text("{}".format(sizes_metabs_y[i]))
    a.set_fontsize(20)
    
#-1.25,0.8
ax6.annotate("", xytext=(-0.7, 0.45), xy=(-1.3, 1.05),
            arrowprops=dict(arrowstyle="-|>, head_width=0.5, head_length=1",
                            color='k'),annotation_clip=False)
ax6.legend(p, xlabs_metabs_y,
          loc=(0.15,0.9), fontsize=15)
# ax6.legend(p, xlabs_metabs_y,
#           loc=(0.2,-0.105), fontsize=15)
ax6.text(-0.6,0.8,'$\\rm \\it{S.~cerevisiae}$',fontsize=25)
# ax6.legend(p, xlabs_metabs_y,
#           loc=(0.3,0.8), fontsize=15)


ax5=plt.axes([0.62,0.15,0.3,0.3])
xlabs_cosub_e= ['NADH','ATP','NADPH','Multiple']
sizes_cosub_e=np.array([21,101,14,42])
p, tx, autotexts=ax5.pie(sizes_cosub_e, autopct='%1.1f%%',radius=0.7,
        colors=['deepskyblue', 'thistle', 'mediumseagreen','darkkhaki'],startangle=90,
        textprops={'size': 20})
for i, a in enumerate(autotexts):
    a.set_text("{}".format(sizes_cosub_e[i]))
    a.set_fontsize(20)
# ax5.text(-0.8,1.1,'E. coli Cytosol Branch Points',fontsize=30)

ax7=plt.axes([0.77,-0.15, 0.5,0.5])
xlabs_metabs_e=('Non-Branch','Co-substrate Utilising \n Reactions (CR)','Without co-substrate')
sizes_metabs_e=np.array([767,178,128])
p, tx, autotexts=ax7.pie(sizes_metabs_e, autopct="",radius=0.7,
        colors=['lightgrey','salmon','wheat'],startangle=220,explode=(0,0.1,0),
        textprops={'size': 20})
for i, a in enumerate(autotexts):
    a.set_text("{}".format(sizes_metabs_e[i]))
    a.set_fontsize(20)
# ax7.annotate("", xytext=(-0.7, 0.45), xy=(-1.2, 0.75),
#             arrowprops=dict(arrowstyle="-|>, head_width=0.5, head_length=1",
#                             color='k'))

ax7.annotate("", xytext=(-0.7, 0.45), xy=(-1.23, 0.8),
            arrowprops=dict(arrowstyle="-|>, head_width=0.5, head_length=1",
                            color='k'))

plt.axis('off')
ax7.text(-0.3,0.8,'$\\rm \\it{E. coli}$',fontsize=25)

# counts_y={
#     'Other': np.array([13,80,13]),
#     'Cytosol':np.array([36,95,41]),
#     }

# width=0.6
# # bottom=np.zeros(3)

# cols=['wheat','salmon']
# # i=0
# for comp, counts_y in counts_y.items():
#         p=ax4.bar(xlabs_y,counts_y,width, label=comp,bottom=bottom,color=cols[i])
#         bottom+=counts_y
        
#         labels_y = [count if count > 0 else '' for count in counts_y] 
#         # print(labels)
        
#         ax4.bar_label(p, labels=labels_y,label_type='center',fontsize=20)
#         i+=1
# ax4.set_title('Yeast',fontsize=30)
# ax4.bar_label(p,label_type='center',fontsize=20)

# handles, labels = plt.gca().get_legend_handles_labels()
# print(handles)
# print(labels)
# order=[1,0]
# ax4.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
#            fontsize=14)#,loc=(0.62,0.8))



# counts_e={
#     'Other': np.array([0,37,0]),
#     'Cytosol':np.array([60,112,51]),
#     }
# bottom=np.zeros(3)
# i=0
# for comp, counts_e in counts_e.items():
#         p=ax5.bar(xlabs_e,counts_e,width, label=comp,bottom=bottom,color=cols[i])
#         bottom+=counts_e
        
#         labels_e = [count if count > 0 else '' for count in counts_e] 
#         # print(labels)
        
#         ax5.bar_label(p, labels=labels_e,label_type='center',fontsize=20)
#         # if comp=='other':
#         #     ax5.bar_label(p,labels=['','37','0'],label_type='center',fontsize=20)
#         i+=1
# ax5.set_yticks([])
# ax5.set_title('E. coli',fontsize=30)
# handles, labels = plt.gca().get_legend_handles_labels()
# print(handles)
# print(labels)
# order=[1,0]
# ax5.legend([handles[idx] for idx in order],[labels[idx] for idx in order],
#            fontsize=14)
# # ax5.bar_label(p,label_type='center',fontsize=20)

# # for comp, counts_e in counts_e.items():
# #     f = counts_e != 0
# #     print(f)
# #     print(xlabs_e[f])
# #     print(counts_e[f])
# #     p = ax5.bar(xlabs_e[f], counts_e[f], width, label=comp,  bottom=bottom)
# #     bottom += counts_e
# #     ax5.bar_label(p, label_type='center')


# # globalmax=max(max(counts_e),max(counts_y))
# ax4.set_ylim(0,200)
# ax5.set_ylim(0,200)
# -1.2
# ax4.text(-8,-1.3,'$\\rm A\\colon Lower~Glycolysis$',fontsize=30)
ax4.text(-8,1.2,'${\\rm A}$: ${\\rm Lower~Glycolysis}$',fontsize=28)
ax4.text(-8,-1.7,'${\\rm B}$: ${\\rm Upper~Glycolysis}$',fontsize=28)
ax4.text(-8,-4.6,'${\\rm C}$: ${\\rm Ammonia~Assimilation}$',fontsize=28)

ax4.text(-1,1.2,'${\\rm D}$: ${\\rm Cytosol~Branch~Points}$',fontsize=28)
# ax4.text(-1.3,-3.8,'$\\rm E$',fontsize=40)


if savefig==1:
    plt.savefig('Figure_1.png',bbox_inches=mpl.transforms.Bbox([[0.39, -0.55], [13.63, 8.68]]),format='png')
    plt.savefig('Figure_1.eps',bbox_inches=mpl.transforms.Bbox([[0.39, -0.55], [13.63, 8.68]]),format='eps')
    plt.savefig('Figure_1.svg',bbox_inches=mpl.transforms.Bbox([[0.39, -0.55], [13.63, 8.68]]),format='svg')
    plt.savefig('Figure_1.pdf',bbox_inches=mpl.transforms.Bbox([[0.39, -0.55], [13.63, 8.68]]),format='pdf')
    
plt.show()
