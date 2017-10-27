# -*- coding: utf-8 -*-
"""
Generates fake data for SVOL model in c++.

simulates from the stochastic volatility model:
  
.. math:: X_n = \\alpha X_{n-1} + \sigma V_n
.. math:: Y_n = \\beta \exp(X_n/2)W_n
.. math:: X_1 \sim \mathcal{N}( 0, \\frac{\sigma^2}{1-\\alpha^2})

The parameters are 
alpha = .91
sigma=1.0
beta = .5

Created on Sun Sep 11 19:47:37 2016

@author: taylor
"""
import pandas as pd
import numpy as np

def sim_svol_data(T, alpha = .91, sigma = 1., beta = .5):
    """simulates from the stochastic volatility model:
    
    .. math:: X_n = \\alpha X_{n-1} + \sigma V_n
    .. math:: Y_n = \\beta \exp(X_n/2)W_n
    .. math:: X_1 \sim \mathcal{N}( 0, \\frac{\sigma^2}{1-\\alpha^2})
    """    
    # noise
    v = np.random.randn(T)
    w = np.random.randn(T)
    
    # state and obs respectively
    x = np.empty(T)
    y = np.empty(T)
    x[0] = np.random.randn(1)* sigma/np.sqrt(1 - alpha**2)
    y[0] = beta * np.exp(x[0]/2.) * w[0]
    for t in range(1,T):
        x[t] = alpha*x[t-1] + sigma*v[t]   
        y[t] = beta * np.exp(x[t]/2.) * w[t]
        
    return (x,y)
    
x,y = sim_svol_data(T=500)
x = pd.DataFrame(x)
y = pd.DataFrame(y)
x.to_csv("/home/taylor/ssm/data/some_csvs/svol_x_data.csv",
         index=False,
         header=False)
y.to_csv("/home/taylor/ssm/data/some_csvs/svol_y_data.csv",
         index=False,
         header=False)



