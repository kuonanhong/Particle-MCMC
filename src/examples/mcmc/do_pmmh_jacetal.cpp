#include "do_pmmh_jacetal.h"

#include "pmmh_jac_apf.h" // for Pmmh_jac_apf

#include <iostream>
void do_pmmh_jacetal()
{
    // ordering is: betas(9-1), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    std::vector<double> st {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                            .7, 0.0, .1,
                            .75, .75, .75, .75, .75, .75, .75, .75, .75};    
                                        
    Pmmh_jac_apf wut(100,                                               // number of particles 
                     st,                                                 // starting thetas
                     200,                                              // number of MCMC iterations
                     "/home/taylor/ssm/data/some_csvs/jacq_y_data.csv",  // where the data is
                     9);                                                 // number of columns of data
    
    // begin the sampling
    wut.commenceSampling("/home/taylor/Desktop/temp1", "/home/taylor/Desktop/temp2");
}
