#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  2 17:15:21 2017

@author: taylor
"""

def gen_rbpfsvol(phi1, phi2, beta, varw1, varw2, T):
    """
    Simulates data from rbpf svol model.
    Follows notation from ~/ 
    
    Parameters
    ----------
    phi1
        a double representing the mean dynamics loadings
    phi2
        float representing log-vol dynamics
    beta 
        modal volatility
    varw1
        float representing wiggly-ness of mean
    varw2
        float representing vol-of-log-vol
    T
        Length of time series.
    """
    x_array = np.empty((T,2))
    y_array = np.empty(T)
    
    # generate first x and y
    x_array[0,0] = np.random.normal(loc=0,scale=np.sqrt(varw1/(1 - phi1**2)) )
    x_array[0,1] = np.random.normal(loc=0,scale=np.sqrt(varw2/(1 - phi2**2)) )
    y_array[0] = np.random.normal(loc=x_array[0,0],
                                  scale=beta*np.exp( x_array[0,1]/2 ))
    
    # generate the rest
    for t in range(1,T):
        x_array[t,0] = np.random.normal(loc=phi1*x_array[t-1,0],
                                        scale=np.sqrt(varw1))
        x_array[t,1] = np.random.normal(loc=phi2*x_array[t-1,1],
                                        scale=np.sqrt(varw2))
        y_array[t] = np.random.normal(loc=x_array[t,0],
                                      scale = beta*np.exp(x_array[t,1]/2) )
        
    return (x_array,y_array)
    
if __name__ == "__main__":
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    # these haven't really been decided yet, though
    beta = .2   # modal volatility
    phi1 = .1   # mean dynamics
    phi2 = .9   # log-vol dynamics
    varw1 = .05 # variance of mean
    varw2 = .2  # vol-of-vol
    T = 300
    x,y = gen_rbpfsvol(phi1, phi2, beta, varw1, varw2, T)
    y = pd.DataFrame(y)
    y.to_csv("/home/taylor/ssm/data/some_csvs/rbpf_svol_y_data.csv",
         index=False,
         header=False)
    
    plt.subplot(3,1,1)
    plt.plot(x[:,0])
    plt.subplot(3,1,2)
    plt.plot(x[:,1])
    plt.subplot(3,1,3)
    plt.plot(y)
    plt.show()
