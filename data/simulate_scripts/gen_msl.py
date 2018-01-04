import numpy as np
from scipy.special import binom
import matplotlib.pyplot as plt 
import pandas as pd

def lexicographically_next_permutation(a):
    """
    Generates the lexicographically next permutation.
    https://pythonadventures.wordpress.com/tag/lexicographical-order/
     
    Input: a permutation, called "a". This method modifies
    "a" in place. Returns True if we could generate a next
    permutation. Returns False if it was the last permutation
    lexicographically.
    """
    i = len(a) - 2
    while not (i < 0 or a[i] < a[i+1]):
        i -= 1
    if i < 0:
        return False
    # else
    j = len(a) - 1
    while not (a[j] > a[i]):
        j -= 1
    a[i], a[j] = a[j], a[i]        # swap
    a[i+1:] = list(reversed(a[i+1:]))    # reverse elements from position i+1 till the end of the sequence
    return True

def sim_msl(beta_1, mu_array, phi_array, sigmaSqd_array, 
                R_array, lambdas, p, K, time_length):
    """
    
    """
    # throw in some assertions
    assert beta_1.ndim == 2
    assert R_array.ndim == 2
    assert R_array.shape[0] == beta_1.shape[0]
    assert R_array.shape[1] == beta_1.shape[0] # make sure we're plugging in a matrix
    assert phi_array.ndim == 1
    assert mu_array.ndim == 1
    assert lambdas.ndim == 1
    assert sigmaSqd_array.ndim == 1
    assert isinstance(p,float)
    assert isinstance(K, int)
    assert isinstance(time_length, int)

    
    # commonly used numbers for indexing
    dim_obs = beta_1.shape[0] 
    num_factors = mu_array.shape[0] 
    dim_state = num_factors + 1
    Tk = sum(int(binom(dim_obs, i)) for i in range(K+1))
    
    # allocate state samples
    x = np.empty((time_length, dim_state))
    y = np.empty((time_length, dim_obs))
    
    # generate all lexicographical orderings
    all_arrays = np.empty((Tk,dim_obs))
    cnt = 0
    for numHot in range(K+1):
        dyn_arr = np.array([0 for i in range(dim_obs)])
        
        for i in range(numHot):
            dyn_arr[i] = 1.0
        dyn_arr.sort()
        
        # generate all permutations and add them to collection
        cond = True
        while(cond):
            all_arrays[cnt,:] = dyn_arr
            cond = lexicographically_next_permutation(dyn_arr)
            cnt += 1
            
    # generate first time's states and obs
    x[0,0] = np.random.choice(range(Tk)) # discrete uniform
    x[0,1:3] = np.random.multivariate_normal(mean = mu_array,
                                           cov  = np.diag(sigmaSqd_array/(1 - phi_array**2)))
    y_cov = np.exp(x[0,1]) * np.dot(beta_1, beta_1.T)
    x11 = int(x[0,0])
    beta_2 = all_arrays[x11,:].reshape((dim_obs,1)) * beta_1 #sparsify first column
    y_cov = y_cov +  np.exp(x[0,2]) * np.dot(beta_2, beta_2.T) + R_array
    y_mean = lambdas[0]*beta_1.reshape((9)) + lambdas[1] * beta_2.reshape((9))
    y[0,:] = np.random.multivariate_normal(mean = y_mean, cov = y_cov)
    
    # generate the rest of the times' states and obs
    for i in range(1,time_length):
        # first state component
        c = (1-p)/(Tk-1)
        probs = np.repeat(c, Tk)
        x1tm1 = int(x[i-1,0])
        probs[x1tm1] = p
        x[i,0] = np.random.choice(range(Tk),p=probs)
        
        # second state component
        second_state_comp_mean = mu_array + phi_array * (x[i-1,1:3] - mu_array)
        x[i,1:3] = np.random.multivariate_normal(mean = second_state_comp_mean,
                                               cov  = np.diag(sigmaSqd_array))
        
        y_cov = np.exp(x[i,1])*np.dot(beta_1, beta_1.T)
        xt1 = int(x[i,0])
        beta_2 = all_arrays[ xt1, :].reshape((dim_obs,1)) * beta_1 # sparsify first col
        y_cov = y_cov + np.exp(x[i,2]) * np.dot(beta_2, beta_2.T) + R_array
        y_mean = lambdas[0]*beta_1.reshape((9)) + lambdas[1]*beta_2.reshape((9))
        y[i,:] = np.random.multivariate_normal(y_mean, y_cov)
    
    # return everything
    return (x,y)
  
if __name__ == "__main__":

    
    # set parameters (HAS NOT BEEN DONE CAREFULLY YET)
    # taking some numbers from the Jacquier paper 1999
    beta1 = np.array([1.00, 1.08, 1.26, 1.46, 0.80, 1.23, 1.10, 0.90, 1.37]).reshape((9,1))    
    mus = np.array([.005, -.01])
    phis = np.array([.62, .86])
    sigma_sqds = np.array([.38, .52])**2 # state vars
    R = np.diag([2.19, 1.45, 0.99, 1.24, 0.83, 0.69, 0.75, 0.89, 0.85])**2 
    lambdas = np.array([.01, -.01])
    p=.94
    K = 3
    time_obs = 500
    
    # generate data
    x,y = sim_msl(beta1, mus, phis, sigma_sqds, R, lambdas, p, K, time_obs)
    y = pd.DataFrame(y)
    x = pd.DataFrame(x)
    
    # save data
    y.to_csv("/home/taylor/ssm/data/some_csvs/msl_y_data.csv",
         index=False,
         header=False)
    x.to_csv("/home/taylor/ssm/data/some_csvs/msl_x_data.csv",
         index=False,
         header=False)
    
    # plot data 
    plt.subplot(3,1,1)
    plt.plot(y)
    plt.subplot(3,1,2)
    plt.plot(np.apply_over_axes(np.cumsum, y, 0))
    plt.subplot(3,1,3)
    plt.plot(x)
    plt.show()
