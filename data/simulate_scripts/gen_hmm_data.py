#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 17:51:29 2017

@author: taylor
"""

import numpy as np

def sim_hmm(T = 50, p = .8, low_var = 1.0, high_var = 3.5):
    """Simulates data from a two state stochastic volatility hmm.
   
    starting distribution is the stationary distribution aka discrete uniform.
    observations are generated from normal distribution with variance 1
    (if in first state), or variance 1.5 (if in second state)
    """
    assert p >= 0.0
   
    init_distn = np.repeat(1/2, 2)
    trans_mat = np.array([[p, 1-p],[1-p, p]])

    x = np.empty((T,1))
    y = np.empty((T,1))
   
    for t in range(T):
       
        if t == 0: # first time
            x[t,0] = np.random.choice(a =[1,2], size=1, p = init_distn)[0]
            if x[t,0] == 1:
                y[t,0] = np.random.normal(scale=np.sqrt(low_var))
            else:
                y[t,0] = np.random.normal(scale=np.sqrt(high_var))
        else: # not first time
            if x[t-1,0] == 1:
                x[t,0] = np.random.choice(a = [1,2], size=1, p = trans_mat[0,:])[0]
            else:
                x[t,0] = np.random.choice(a = [1,2], size=1, p = trans_mat[1,:])[0]
           
            if x[t,0] == 1:
                y[t,0] = np.random.normal(scale=np.sqrt(low_var))
            else:
                y[t,0] = np.random.normal(scale=np.sqrt(high_var))
               
    return (x,y)

if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    import pandas as pd
    import numpy as np
    
    x,y = sim_hmm(T=150, low_var=1, high_var=4.5, p=.9)
    y = pd.DataFrame(y)
    x = pd.DataFrame(x)
    
    # save data
    y.to_csv("/home/taylor/ssm/data/some_csvs/hmm_y_data.csv",
         index=False,
         header=False)
    x.to_csv("/home/taylor/ssm/data/some_csvs/hmm_x_data.csv",
         index=False,
         header=False)
    
    # plot data 
    plt.subplot(2,1,1)
    plt.plot(y)
    plt.subplot(2,1,2)
    plt.plot(x)
    plt.show()