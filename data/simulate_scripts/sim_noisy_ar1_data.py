# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 17:05:07 2016
Simulates from the noisy ar1 model:
    
    .. math:: X_n = \\alpha X_{n-1} + W_n
    .. math:: Y_n = X_n + V_n
with
    .. math:: W_n \overset{iid}{\sim} \mathcal{N}(0,\sigma^2)
    .. math:: V_n \overset{iid}{\sim} \mathcal{N}(0, \gamma^2)
    .. math:: X_1 \sim \mathcal{N}(0, \\frac{\sigma^2}{1-\\alpha^2})
    
Default values: 
    alpha = .91, sigma = 1., and gamma = 1.5

@author: taylor
"""

import pandas as pd
import numpy as np

def sim_noisy_ar1(T, alpha=.91, sigma = 1., gamma=1.5):
    """simulates from the noisy ar1 model:
    
    .. math:: X_n = \\alpha X_{n-1} + W_n
    .. math:: Y_n = X_n + V_n
    .. math:: W_n \overset{iid}{\sim} \mathcal{N}(0,\sigma^2)
    .. math:: V_n \overset{iid}{\sim} \mathcal{N}(0, \gamma^2)
    .. math:: X_1 \sim \mathcal{N}(0, \\frac{\sigma^2}{1-\\alpha^2})
    """
    v = gamma * np.random.randn(T) 
    w = sigma * np.random.randn(T)
    x = np.empty(T)
    y = np.empty(T)
    
    x[0] = np.random.randn(1)* (sigma/np.sqrt(1 - alpha**2))
    for t in range(1,T):
        x[t] = alpha * x[t-1] + w[t]
    y = x + v
    
    return x,y
    
    

x,y = sim_noisy_ar1(T=200, alpha = .31, sigma = 1.0, gamma = 1.5 )
y = pd.DataFrame(y)
x = pd.DataFrame(x)
y.to_csv("/home/taylor/ssm/data/some_csvs/noisy_ar1_data.csv",
         index=False,
         header=False)
x.to_csv("/home/taylor/ssm/data/some_csvs/noisy_ar1_states.csv",
         index=False,
         header=False)