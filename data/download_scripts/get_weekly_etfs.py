#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 29 18:14:37 2017

@author: taylor
"""
import pandas as pd
import numpy as np

from pandas_datareader import data
from datetime import datetime

# define window and symbols of interest
symbols= ['xle', 'xlu', 'xlk', 'xlb', 'xlp', 'xly', 'xli', 'xlv', 'xlf'] #'xlre', 
start = datetime(2005, 12, 23)
finish = datetime(2014, 5, 1)

# get data
df = data.DataReader('xle', data_source='google', start=start, end=finish)['Close']
for s in symbols[1:]:
    print ('getting data for symbol: ', s)
    df = pd.concat([df, 
                    data.DataReader(s, data_source='google', start=start, end=finish)['Close']],
                    axis=1)
    
# resample to weekly data
rweekly = df.resample('W')
wdf = rweekly.last()

# make into log returns
wr = wdf.apply(np.log).diff()
wr = wr.iloc[1:,:]

# scale by 100 
wr = wr*100

# write to csv
wr.to_csv("~/ssm/data/some_csvs/weekly_etf_data_200151223_201451", 
          header=False,
          index=False)

