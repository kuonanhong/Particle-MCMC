
def sim_jacetal(beta_array, mu_array, phi_array, sigmaSqd_array, R_array, time_length):
    """
    Simulates data from "Stochastic Volatility: Univariate and Multivariate Extensions" by Jacquier et al. 
    
    
    Follows notation from preliminary exam document
    
    """
    # throw in some assertions
    assert beta_array.ndim == 2
    assert R_array.ndim == 2
    assert R_array.shape[0] == beta_array.shape[0]
    assert R_array.shape[1] == beta_array.shape[0] # make sure we're plugging in a matrix
    assert phi_array.ndim == 1
    assert mu_array.ndim == 1
    assert sigmaSqd_array.ndim == 1
    assert type(time_length) == int
    assert beta_array.shape[1] == phi_array.shape[0]
    assert beta_array.shape[1] == sigmaSqd_array.shape[0]

    
    # commonly used numbers for indexing
    dim_obs = beta_array.shape[0] 
    num_factors = beta_array.shape[1] 
    dim_state = num_factors
    
    # allocate state samples
    x = np.empty((time_length, dim_state))
    y = np.empty((time_length, dim_obs))
    
    # generate first time's states and obs
    x[0,:] = np.random.multivariate_normal(mean = mu_array,
                                           cov  = np.diag(sigmaSqd_array/(1 - phi_array**2)))
    y_cov = np.dot(beta_array, np.diag(np.exp(x[0,:])))
    y_cov = np.dot(y_cov, beta_array.T) + R_array
    y[0,:] = np.random.multivariate_normal(mean = np.zeros(dim_obs), cov = y_cov)
    
    # generate the rest of the times' states and obs
    for i in range(1,time_length):
        tmp_mean = mu_array + phi_array * (x[i-1,:] - mu_array)
        x[i,:] = np.random.multivariate_normal(mean = tmp_mean,
                                               cov  = np.diag(sigmaSqd_array))
        y_cov = np.dot(beta_array, np.diag(np.exp(x[i,:])))
        y_cov = np.dot(y_cov, beta_array.T) + R_array
        y[i,:] = np.random.multivariate_normal(np.zeros(dim_obs), y_cov)
    
    # return everything
    return (x,y)
  
if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    import pandas as pd
    import numpy as np
    
    # set parameters (HAS NOT BEEN DONE CAREFULLY YET)
    # taking some numbers from the Jacquier paper 1999
    beta = np.array([1.1, 1.0, .9, 1.1, 1.0, .9, 1.1, 1.0, .9]).reshape((9,1))    
    phi = np.repeat(.95, 1) # 0.894736842105263, 10)
    mu = np.array([.09])
    sigma_sqds = np.repeat(.08975, 1) # state vars
    time_obs = 1000
    R = np.diag(np.repeat(1.5,9)) # the paper does not assume this is diagonal
    
    # generate data
    x,y = sim_jacetal(beta, mu, phi, sigma_sqds, R, time_obs)
    y = pd.DataFrame(y)
    x = pd.DataFrame(x)
    
    # save data
    y.to_csv("/home/taylor/ssm/data/some_csvs/jacq_y_data.csv",
         index=False,
         header=False)
    x.to_csv("/home/taylor/ssm/data/some_csvs/jacq_x_data.csv",
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
