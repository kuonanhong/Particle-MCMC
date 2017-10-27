# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#####################
# conditional means #
#####################
filterFilePath = "/home/taylor/SSMTools/cpp/examples/noisy_ar1/mean_output.csv"
statesFilePath = "/home/taylor/SSMTools/cpp/examples/noisy_ar1/noisy_ar1_states.csv"

df = pd.io.api.read_csv(filterFilePath)
states = pd.io.api.read_csv(statesFilePath, header=None)
df[['actualmean','SISRmean', 'APFmean']].plot()

plt.subplot(2,1,1)
plt.plot(df['actualmean'] - df['SISRmean'])
plt.subplot(2,1,2)
plt.plot(df['actualmean'] - df['APFmean'])


# compare filtered means with actual
kalman, = plt.plot(df[' closed form mean'], label='kalman filter')
pfilt, = plt.plot(df[' particle mean '], label = 'particle filter')
actual, = plt.plot(states[0], label='actual state')
plt.legend(handles=[kalman, pfilt, actual], loc='lower right')

###############################
# log conditional likelihoods #
###############################
lclPath = "/home/taylor/SSMTools/cpp/examples/noisy_ar1/cond_like_output.csv"
ldf = pd.io.api.read_csv(lclPath)
badLclPath = "/home/taylor/SSMTools/cpp/examples/noisy_ar1/bad_cond_like_output.csv"
bldf = pd.io.api.read_csv(badLclPath)

plt.subplot(2,1,1)
plt.hist(ldf.iloc[:,1] - ldf.iloc[:,2]) # difference between real and bootstrap
plt.subplot(2,1,2)
plt.hist(ldf.iloc[:,1] - ldf.iloc[:,3]) # difference between real and apf

plt.subplot(2,1,1)
plt.plot(ldf.iloc[:,3] - ldf.iloc[:,1])
plt.subplot(2,1,2)
plt.plot(bldf.iloc[:,3] - bldf.iloc[:,1] )

ldf.iloc[:,1:].plot()

plt.plot(ldf['APFloglike'] - ldf['actualloglike'])


############
# for pmcmc
#############
df = pd.io.api.read_csv("/home/taylor/Desktop/example.csv", header=None)

df.plot()
plt.axhline(y=1.5, color='blue')
plt.axhline(y=1.0, color="green")
plt.axhline(y=.91, color="red")

plt.subplots_adjust(hspace=.4)
plt.subplot(3,1,1)
plt.hist(df.iloc[:,0])
plt.axvline(1.5, color='b', linestyle='dashed', linewidth=2)
plt.title("observational std dev")
plt.subplot(3,1,2)
plt.hist(df.iloc[:,1])
plt.axvline(1.0, color='b', linestyle='dashed', linewidth=2)
plt.title("state transition std dev")
plt.subplot(3,1,3)
plt.hist(df.iloc[:,2])
plt.axvline(.91, color='b', linestyle='dashed', linewidth=2)
plt.title("ar1 coefficient")

pd.tools.plotting.scatter_matrix(df, alpha=.05)