#include "do_pmmh_jacetal.h"

#include "pmmh_jac_apf.h" // for Pmmh_jac_apf

void do_pmmh_jacetal()
{
    // ordering is: betas(9-1), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    std::vector<Vec> st;
    Vec startBeta(9);
    startBeta << 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    Vec startPhi(1);
    startPhi << .7;
    Vec startMu(1);
    startMu << 0.0;
    Vec startSigma(1);
    startSigma << .1;
    Vec startRStdDevs(9);
    startRStdDevs << .75, .75, .75, .75, .75, .75, .75, .75, .75;
    st.push_back(startBeta);
    st.push_back(startPhi);
    st.push_back(startMu);
    st.push_back(startSigma);
    st.push_back(startRStdDevs);
                                        
    Pmmh_jac_apf wut(100,                                               // number of particles 
                     st,                                                 // starting thetas
                     200,                                              // number of MCMC iterations
                     "/home/taylor/ssm/data/some_csvs/jacq_y_data.csv",  // where the data is
                     9,                                                 // number of columns of data
                     true);

    // begin the sampling
    wut.commenceSampling("/home/taylor/Desktop/temp1", "/home/taylor/Desktop/temp2");
}
