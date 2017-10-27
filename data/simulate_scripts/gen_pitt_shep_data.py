
def sim_msvol(beta_array, phi_array, mu_array, var_array, T):
    """
    Simulates data from time-varying covariances paper by Shephard and Pitt. 
    Follows notation from ~/UVa/research/TeX/time-varying-covs/tvcovs.pdf 
    
   Parameters
   ----------
   phi_array
       a 2-d numpy array representing the loadings
   mu1
       float

    """
    # throw in some assertions
    assert beta_array.ndim == 2
    assert phi_array.ndim == 1
    assert mu_array.ndim == 1
    assert var_array.ndim == 1
    assert type(T) == int
    assert beta_array.shape[0] + beta_array.shape[1] == phi_array.shape[0]
    assert beta_array.shape[0] + beta_array.shape[1] == var_array.shape[0]
    assert beta_array.shape[0] == mu_array.shape[0]

    
    # commonly used numbers for indexing
    dim_obs = beta_array.shape[0] #N
    num_factors = beta_array.shape[1] #K
    dim_state = dim_obs + num_factors
    
    # for readability
    u_var_array = var_array[0:num_factors]
    v_var_array = var_array[num_factors:]
    factor_phis = phi_array[0:num_factors]
    errors_phis = phi_array[num_factors:]
    
    # allocate state samples
    u = np.empty((T, num_factors))
    v = np.empty((T, dim_obs))
    x = np.empty((T,dim_state))
    y = np.empty((T, dim_obs))
    
    # generate states
    u[0,:] = np.random.multivariate_normal(np.repeat(0, num_factors),
                                       np.diag(u_var_array/(1 - factor_phis**2)))
    v[0,:] = np.random.multivariate_normal(np.repeat(0, dim_obs),
                                       np.diag(v_var_array/(1 - errors_phis**2)))
    y_cov = np.dot(beta_array, np.diag(np.exp(u[0,:])))
    y_cov = np.dot(y_cov, beta_array.T) + np.diag(np.exp(v[0,:]))
    y[0,:] = np.random.multivariate_normal(np.zeros(dim_obs), y_cov)
    for i in range(1,T):
        u[i,:] = np.random.multivariate_normal(factor_phis * u[i-1,:],
                                               np.diag(u_var_array))
        v[i,:] = np.random.multivariate_normal(mu_array+ errors_phis*(v[i-1,:] - mu_array),
                                               np.diag(v_var_array))
        y_cov = np.dot(beta_array, np.diag(np.exp(u[i,:])))
        y_cov = np.dot(y_cov, beta_array.T) + np.diag(np.exp(v[i,:]))
        y[i,:] = np.random.multivariate_normal(np.zeros(dim_obs), y_cov)
    
    x[:,0:num_factors] = u
    x[:,num_factors:] = v
    return (x,y)
  
if __name__ == "__main__":
    import matplotlib.pyplot as plt 
    import pandas as pd
    import numpy as np
    
    # prior means of time-varying-covs paper
    beta = np.ones(9).reshape((9,1))    
    phi = np.repeat(.92, 10)#0.894736842105263, 10)
    mu = np.repeat(1.0, 9)
    var = np.repeat(.2, 10)
    T = 2000
    x,y = sim_msvol(beta, phi, mu, var, T)
    y = pd.DataFrame(y)
    x = pd.DataFrame(x)
    y.to_csv("/home/taylor/ssm/data/some_csvs/msvol_y_data.csv",
         index=False,
         header=False)
    x.to_csv("/home/taylor/ssm/data/some_csvs/msvol_x_data.csv",
         index=False,
         header=False)
    
    plt.subplot(2,1,1)
    plt.plot(y)
    plt.subplot(2,1,2)
    plt.plot(np.apply_over_axes(np.cumsum, y, 0))
    plt.show()
