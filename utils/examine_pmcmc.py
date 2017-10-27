#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
True parameters:
    betas  = 1.0;
    mus    = 1.0
    phis   = .92
    sigmas = .1118 = sqrt(.0125);

Priors:
    see tv-covs document
        
Starting points: ????????????/
    beta   =  1.5 (except for the first one its fixed at 1)
    mus    = -1.0
    phis   = 0.0
    sigmas = .5
    
Created on Thu Jan 12 18:20:27 2017

@author: taylor
"""
import numpy as np
import pandas as pd
from pandas.tools.plotting import autocorrelation_plot

#myDir = '/home/taylor/SSMTools/cpp/mcmc/output/' 
#myDir = '/home/taylor/UVa/research/cluster_output/50k/3/'
myDir = '/home/taylor/UVa/research/cluster_output/100k/3/'
sampleFileName1 = 'outputfile_1.csv'
#sampleFileName2 = 'outputfile_2.csv'

pmcmcdata = pd.io.api.read_csv(myDir+sampleFileName1, header=None,sep=',')#sep ='\s+')

#pmcmcdata.plot(legend=False)

# separate out all data
betas = pmcmcdata.ix[:,:8]
phis = pmcmcdata.ix[:,9:18]
mus = pmcmcdata.ix[:,19:27]
sigmas = pmcmcdata.ix[:,28:37]

#########
# PLOTS #
#########

#betas
betas.plot(legend=False)
betas.ix[50000:,:].plot()
autocorrelation_plot(betas.ix[50000::500,1])
betas.ix[50000::500,1].hist()

#mus
mus.plot(legend=False)
mus.ix[50000:,:].plot()
autocorrelation_plot(mus.iloc[50000::500,1])
mus.iloc[50000:,1].hist()


# phis (constrained)
def twiceFisher(x): return(np.log((1+x)/(1-x)))
phis.plot(legend=False)
phis.iloc[50000:,:].plot()
phis.iloc[50000:,0].apply(twiceFisher).hist()
phis.iloc[50000:,1].hist()
phis.iloc[50000:,2].hist()
phis.iloc[50000:,3].hist()
phis.iloc[50000:,4].hist()

# sigmas
sigmas.plot(legend=False)
sigmas.ix[50000:,:].apply(np.log).plot()


# plot the market factor
#plt.plot(filtData.ix[:,0])
#
##  plot the idiosyncratic error factors
#plt.plot(filtData.ix[:,1:])
#
## plot the comparisons
#for i in range(1,11):
#    plt.subplot(10,1,i)
#    plt.plot(filtData.ix[:,i-1])
#    plt.plot(x.ix[:,i-1])
#
## plot the differences
#for i in range(1,11):
#    plt.subplot(10,1,i)
#    #plt.plot((filtData - x).ix[:,i-1])
#    plt.hist((filtData - x).ix[:,i-1])