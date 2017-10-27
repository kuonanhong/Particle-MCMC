#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 09:33:50 2017

@author: taylor
"""
import pandas as pd
import matplotlib.pyplot as plt
import os
from pandas.tools.plotting import autocorrelation_plot

root = "/home/taylor/ssm/mcmc/output/"
first = os.path.join(root, "rbpf_mod1_1kparts_10kiters.csv")
second = os.path.join(root, "sisr_mod1_1kparts_10kiters.csv")

sisr_samples = pd.io.api.read_csv(second, sep = ",", header=None)
rbpf_samples = pd.io.api.read_csv(first, sep = ",", header=None)

# betas
plt.subplot(1,2,1)
sisr_samples.ix[:,0].plot()
plt.axhline(y=0.2, color='r', linestyle='-')
plt.subplot(1,2,2)
rbpf_samples.ix[:,0].plot()
plt.axhline(y=0.2, color='r', linestyle='-')

autocorrelation_plot(sisr_samples.ix[:,0])
autocorrelation_plot(rbpf_samples.ix[:,0])


# phi1
plt.subplot(1,2,1)
sisr_samples.ix[:,1].plot()
plt.axhline(y=0.1, color='r', linestyle='-')
plt.subplot(1,2,2)
rbpf_samples.ix[:,1].plot()
plt.axhline(y=0.1, color='r', linestyle='-')

autocorrelation_plot(sisr_samples.ix[:,1])
autocorrelation_plot(rbpf_samples.ix[:,1])



# phi2
plt.subplot(1,2,1)
sisr_samples.ix[:,2].plot()
plt.axhline(y=0.9, color='r', linestyle='-')
plt.subplot(1,2,2)
rbpf_samples.ix[:,2].plot()
plt.axhline(y=0.9, color='r', linestyle='-')

autocorrelation_plot(sisr_samples.ix[:,2])
autocorrelation_plot(rbpf_samples.ix[:,2])



# mean var
plt.subplot(1,2,1)
sisr_samples.ix[:,3].plot()
plt.axhline(y=0.05, color='r', linestyle='-')
plt.subplot(1,2,2)
rbpf_samples.ix[:,3].plot()
plt.axhline(y=0.05, color='r', linestyle='-')

autocorrelation_plot(sisr_samples.ix[:,3])
autocorrelation_plot(rbpf_samples.ix[:,3])


# var var
plt.subplot(1,2,1)
sisr_samples.ix[:,4].plot()
plt.axhline(y=0.2, color='r', linestyle='-')
plt.subplot(1,2,2)
rbpf_samples.ix[:,4].plot()
plt.axhline(y=0.2, color='r', linestyle='-')

autocorrelation_plot(sisr_samples.ix[:,4])
autocorrelation_plot(rbpf_samples.ix[:,4])

