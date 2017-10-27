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
filterFile = 'known_param_jacq_filter_1kparts.csv'
realStateFile = 'jacq_x_data.csv'
realYFile = 'jacq_y_data.csv'


filtData = pd.read_csv(results_dir + filterFile, header=None,sep ='\s+')
x = pd.read_csv(data_dir + realStateFile, header=None)

# compare true hidden data and the estimated factor
pd.concat([filtData, x], axis=1).plot()

# or plot the difference
(x - filtData).plot(kind="hist")