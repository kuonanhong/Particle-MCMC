#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 18:20:27 2017

@author: taylor
"""
import pandas as pd
import matplotlib.pyplot as plt

data_dir = '/home/taylor/ssm/data/some_csvs/' 
results_dir = '/home/taylor/ssm/examples/msvol_filtering/'
filterFile = 'known_param_filter_data_1kparts.csv'
realStateFile = 'msvol_x_data.csv'
realYFile = 'msvol_y_data.csv'


filtData = pd.read_csv(results_dir + filterFile, header=None,sep ='\s+')
x = pd.read_csv(data_dir + realStateFile, header=None)

# plot the market factor
filtData.iloc[:,0].plot()

#  plot the idiosyncratic error factors
filtData.iloc[:,1:].plot()

# plot the comparisons
for i in range(1,11):
    plt.subplot(10,1,i)
    plt.plot(filtData.iloc[:,i-1])
    plt.plot(x.iloc[:,i-1])

# plot the differences
for i in range(1,11):
    plt.subplot(10,1,i)
    #plt.plot((filtData - x).iloc[:,i-1])
    plt.hist((filtData - x).iloc[:,i-1])
