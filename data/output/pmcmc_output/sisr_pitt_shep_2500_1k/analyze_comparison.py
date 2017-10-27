#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 09:33:50 2017

@author: taylor
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from pandas.tools.plotting import autocorrelation_plot


# read in data
root = "/home/taylor/ssm/mcmc/output/sisr_pitt_shep_2500_1k/"
third = os.path.join(root, "outputfile_3.csv")
pitt_shep_samps = pd.read_csv(third, sep = ",", header = None)

# betas (first nine) (true is 1)
pitt_shep_samps.iloc[:,:9].plot()
plt.xlabel('iteration')
plt.ylabel('betas')

# phi1 (next ten) (true is .92) (started at 0)
pitt_shep_samps.iloc[:,9:19].plot()
plt.axhline(y=0.92, color='r', linestyle='-')
plt.xlabel('iteration')
plt.ylabel('phis')

# mus (next 9) (true is 1.0) (start at -1)
pitt_shep_samps.iloc[:,19:28].plot()
plt.axhline(y=1.0, color='r', linestyle='-')
plt.xlabel('iteration')
plt.ylabel('mus')

# sigmas (last ten) ( true is .2)
pitt_shep_samps.iloc[:,-10:].plot()
plt.axhline(y=0.2, color='r', linestyle='-')
plt.xlabel('iteration')
plt.ylabel('sigmas')


# get a general idea of how correlated these guys are
derp = pd.concat([pitt_shep_samps.iloc[:,1:9].mean(1),
             pitt_shep_samps.iloc[:,9:19].mean(1),
             pitt_shep_samps.iloc[:,19:28].mean(1),
             pitt_shep_samps.iloc[:,-10:].mean(1)], 
             axis=1)
derp.columns = ['betas', 'phis', 'mus', 'sigmas']
pd.tools.plotting.scatter_matrix(derp, alpha=0.2)

# this takes a while
def twiceFisher(phi):
    return np.log(1.0 + phi) - np.log(1.0 - phi)
trans_data = pitt_shep_samps.iloc[:,:9] # betas aren't transformed
trans_data = pd.concat([trans_data, 
                       pitt_shep_samps.iloc[:,9:19].apply(twiceFisher)], 
                       axis=1)
trans_data = pd.concat([trans_data, 
                       pitt_shep_samps.iloc[:,19:28]],
                       axis=1)
trans_data = pd.concat([trans_data, 
                       pitt_shep_samps.iloc[:,-10:].apply(np.log)],
                       axis=1)
pd.tools.plotting.scatter_matrix(trans_data.iloc[:,1:], alpha=.2)


#plt.subplot(1,2,1)
#sisr_samples.ix[:,0].plot()
#plt.axhline(y=0.2, color='r', linestyle='-')
#plt.subplot(1,2,2)
#rbpf_samples.ix[:,0].plot()
#plt.axhline(y=0.2, color='r', linestyle='-')
#
#autocorrelation_plot(sisr_samples.ix[:,0])
#autocorrelation_plot(rbpf_samples.ix[:,0])