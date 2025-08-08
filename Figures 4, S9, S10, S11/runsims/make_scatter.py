#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 24 14:38:23 2025

@author: robert
"""

import numpy as np

import time
from datetime import timedelta
import multiprocessing as multi
import sys

#from pathlib import Path
import os
import time

import itertools

# sys.path.append('../data/megadata')
# sys.path.append('../data/heatmapdata')
sys.path.append('../simulation_functions')
sys.path.append('../data')


import make_heatmaps_compartment as mh_cm

import matplotlib.pyplot as plt
import numpy.random as rn
import datetime as dt


aldo=[1.9]
for al in aldo:
	t0=time.time()

	allo=[al,1.,1.,1.]
	n=1000

	seed=1000

	if seed!=False:
		np.random.seed(seed)
		seedstr=str(seed)
	else:
		seedstr='none'
		np.random.seed(None)
		
	args_to_run=mh_cm.make_params_newtest_6_nonfixed(allo,n)

	now=dt.datetime.now()
	print(now)

	savepath='../data/allo_'+str(allo)
	filename='/nonfixed_apr_final_seed'+seedstr+'_n_'+str(n)+'_datetime'+str(now)
	savename=savepath+filename

	print(savepath)
	if not os.path.exists(savepath):
		print('Directory does not exist')
		os.makedirs(savepath)
		print("Directory created successfully")
	else:
		print("Directory already exists")

	# print(args_to_run)
	result=mh_cm.make_data_parallel(mh_cm.get_flux_vecs,args_to_run)


	# ft_vec=result[0][0]
	# fb_vec=result[0][1]
	# fd_vec=result[0][2]
	# ff_vec=result[0][3]
	# ff_grad=result[0][4]
	# ff_sens=result[0][5]
	# ff_delta=result[0][6]
	# gradmax=result[0][7]
	# gradmaxind=result[0][8]
	# sensemax=result[0][9]
	# sensemaxind=result[0][10]

	ftlist=list()
	fblist=list()
	fdlist=list()
	fflist=list()
	ffgradlist=list()
	ffsenselist=list()
	ffdeltalist=list()
	gradmaxlist=list()
	gradmaxindlist=list()
	sensemaxlist=list()
	sensemaxindlist=list()
	b0deltalist=list()
	c0deltalist=list()
	n0deltalist=list()
	swptlist=list()
	swptvalslist=list()


	#concentrationslist=list()

	for i in range(n):
		ft_vec=result[i][0]
		ftlist.append(ft_vec)

		fb_vec=result[i][1]
		fblist.append(fb_vec)

		fd_vec=result[i][2]
		fdlist.append(fd_vec)

		ff_vec=result[i][3]
		fflist.append(ff_vec)

		ff_grad=result[i][4]
		ffgradlist.append(ff_grad)

		ff_sense=result[i][5]
		ffsenselist.append(ff_sense)

		ff_delta=result[i][6]
		ffdeltalist.append(ff_delta)

		gradmax=result[i][7]
		gradmaxlist.append(gradmax)

		gradmaxind=result[i][8]
		gradmaxindlist.append(gradmaxind)

		sensemax=result[i][9]
		sensemaxlist.append(sensemax)

		sensemaxind=result[i][10]
		sensemaxindlist.append(sensemaxind)

		b0delta=result[i][11]
		b0deltalist.append(b0delta)

		c0delta=result[i][12]
		c0deltalist.append(c0delta)

		n0delta=result[i][13]
		n0deltalist.append(n0delta)

		swpt=result[i][15]
		swptlist.append(swpt)
        
		swptval=result[i][16]
		swptvalslist.append(swptval)
        
        

	 #   concentrations=result[i][11]
	  #  concentrationslist.append(concentrations)
	# print(ffdeltalist)
	# print(ffdeltalist[0])
	# print(ffdeltalist[1])

	# print(gradmaxindlist)
	# print(gradmaxindlist[0])
	# print(gradmaxindlist[1])

	# print(sensemaxindlist)
	# print(sensemaxindlist[0])
	# print(sensemaxindlist[1])

	t1=time.time()

	timetaken=dt.timedelta(seconds=t1-t0)
	print(f'making data took {timetaken} seconds')


	np.savez(savename, flux_cyto=ftlist, flux_mito=fblist, flux_diffn=fdlist,
			 fluxfrac=fflist, fluxfrac_grad=ffgradlist, fluxfrac_sense=ffsenselist,
			 fluxfrac_delta=ffdeltalist, gradmax=gradmaxlist, gradmaxind=gradmaxindlist,
			 sensemax=sensemaxlist, sensemaxind=sensemaxindlist, b0delta=b0deltalist,
			 c0delta=c0deltalist, n0delta=n0deltalist, swpt=swptlist,swptvals=swptvalslist)


# # ff figure
# plt.figure()
# plt.scatter(ffdeltalist,sensemaxindlist)
# #plt.xscale('log')

# plt.axes([1.05, 0.55, 0.4, 0.4])
# plt.scatter(ffdeltalist,gradmaxindlist)

# #plt.xscale('log')

# # plt.legend(loc='lower right')

# plt.axes([1.05, 0.05, 0.4,0.4])
# plt.scatter(sensemaxindlist,gradmaxindlist)

# #plt.xscale('log')

# plt.figure()
# plt.plot(concentrations)
# plt.yscale('log')
# plt.show()


aldo=[1.]

for al in aldo:
	t0=time.time()

	allo=[al,1.,1.,1.]
	n=1000

	seed=1000

	if seed!=False:
		np.random.seed(seed)
		seedstr=str(seed)
	else:
		seedstr='none'
		np.random.seed(None)
		
	args_to_run=mh_cm.make_params_newtest_6(allo,n)

	now=dt.datetime.now()
	print(now)

	savepath='../data/allo_'+str(allo)
	filename='/fixed_apr_final_seed'+seedstr+'_n_'+str(n)+'_datetime'+str(now)
	savename=savepath+filename

	print(savepath)
	if not os.path.exists(savepath):
		print('Directory does not exist')
		os.makedirs(savepath)
		print("Directory created successfully")
	else:
		print("Directory already exists")

	# print(args_to_run)
	result=mh_cm.make_data_parallel(mh_cm.get_flux_vecs,args_to_run)


	# ft_vec=result[0][0]
	# fb_vec=result[0][1]
	# fd_vec=result[0][2]
	# ff_vec=result[0][3]
	# ff_grad=result[0][4]
	# ff_sens=result[0][5]
	# ff_delta=result[0][6]
	# gradmax=result[0][7]
	# gradmaxind=result[0][8]
	# sensemax=result[0][9]
	# sensemaxind=result[0][10]

	ftlist=list()
	fblist=list()
	fdlist=list()
	fflist=list()
	ffgradlist=list()
	ffsenselist=list()
	ffdeltalist=list()
	gradmaxlist=list()
	gradmaxindlist=list()
	sensemaxlist=list()
	sensemaxindlist=list()
	b0deltalist=list()
	c0deltalist=list()
	n0deltalist=list()
	swptlist=list()
	swptvalslist=list()
	#concentrationslist=list()

	for i in range(n):
		ft_vec=result[i][0]
		ftlist.append(ft_vec)

		fb_vec=result[i][1]
		fblist.append(fb_vec)

		fd_vec=result[i][2]
		fdlist.append(fd_vec)

		ff_vec=result[i][3]
		fflist.append(ff_vec)

		ff_grad=result[i][4]
		ffgradlist.append(ff_grad)

		ff_sense=result[i][5]
		ffsenselist.append(ff_sense)

		ff_delta=result[i][6]
		ffdeltalist.append(ff_delta)

		gradmax=result[i][7]
		gradmaxlist.append(gradmax)

		gradmaxind=result[i][8]
		gradmaxindlist.append(gradmaxind)

		sensemax=result[i][9]
		sensemaxlist.append(sensemax)

		sensemaxind=result[i][10]
		sensemaxindlist.append(sensemaxind)

		b0delta=result[i][11]
		b0deltalist.append(b0delta)

		c0delta=result[i][12]
		c0deltalist.append(c0delta)

		n0delta=result[i][13]
		n0deltalist.append(n0delta)
        
		swpt=result[i][15]
		swptlist.append(swpt)
 
		swptval=result[i][16]
		swptvalslist.append(swptval)
        
        

	 #   concentrations=result[i][11]
	  #  concentrationslist.append(concentrations)
	# print(ffdeltalist)
	# print(ffdeltalist[0])
	# print(ffdeltalist[1])

	# print(gradmaxindlist)
	# print(gradmaxindlist[0])
	# print(gradmaxindlist[1])

	# print(sensemaxindlist)
	# print(sensemaxindlist[0])
	# print(sensemaxindlist[1])

	t1=time.time()

	timetaken=dt.timedelta(seconds=t1-t0)
	print(f'making data took {timetaken} seconds')


	np.savez(savename, flux_cyto=ftlist, flux_mito=fblist, flux_diffn=fdlist,
			 fluxfrac=fflist, fluxfrac_grad=ffgradlist, fluxfrac_sense=ffsenselist,
			 fluxfrac_delta=ffdeltalist, gradmax=gradmaxlist, gradmaxind=gradmaxindlist,
			 sensemax=sensemaxlist, sensemaxind=sensemaxindlist, b0delta=b0deltalist,
			 c0delta=c0deltalist, n0delta=n0deltalist, swpt=swptlist, swptvals=swptvalslist)


# # ff figure
# plt.figure()
# plt.scatter(ffdeltalist,sensemaxindlist)
# #plt.xscale('log')

# plt.axes([1.05, 0.55, 0.4, 0.4])
# plt.scatter(ffdeltalist,gradmaxindlist)

# #plt.xscale('log')

# # plt.legend(loc='lower right')

# plt.axes([1.05, 0.05, 0.4,0.4])
# plt.scatter(sensemaxindlist,gradmaxindlist)

# #plt.xscale('log')

# plt.figure()
# plt.plot(concentrations)
# plt.yscale('log')
# plt.show()
