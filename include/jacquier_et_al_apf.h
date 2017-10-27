#ifndef JACQUIER_ET_AL_APF_H
#define JACQUIER_ET_AL_APF_H

#include "apf_filter.h"
#include "densities.h"

//! An APF implementation of the multivariate stochastic volatility model of Jacquier et al.
/**
 * @class JacEtAlAPF
 * @author taylor
 * @date 09/09/17
 * @file jacquier_et_al_apf.h
 * @brief An APF implementation of the multivariate stochastic volatility model 
 * from "Stochastic Volatility: Univariate and Multivariate Extensions"
 */
class JacEtAlAPF : public APFFilter
{
private:
    densities::EigenMultivariateNormalSampler m_timeOneSampler;
    densities::EigenMultivariateNormalSampler m_transJumpSampler;
    
    const Mat m_beta;   // factor loadings
    const Vec m_R_sigmas;      // observatonal covariance matrix in vector form (std. devs though)
    const Vec m_mus;
    const Vec m_phis;   // state transition diagonals
    const Vec m_sigmas; // state std devs
    int m_dim_obs;
    int m_num_factors;
public:

    /**
     * @brief The constructor.
     * @param numParts the number of particles you want (as an int).
     * @param beta a Mat representing the factor loadings matrix.
     * @param R_sigmas_as_vec a Vec representing the square root of the diagonals of the observation covariance matrix R.
     * @param mus a Vec representing the mean of the log-volatility process.
     * @param phis a Vec representing the AR parameters of the log-volatility processes.
     * @param sigmas a Vec representing the standard deviation of the log-volatility processes.
     * @param resampTechnique an enum class representing the type of resampling you want to do.
     * @param pathLength set to 0 if filtering. Otherwise, if entire path samples are desired, set to the time-length of your data.
     * @return 
     */
    JacEtAlAPF(int numParts, const Mat &beta, const Vec &R_sigmas_as_vec, const Vec &mus, const Vec &phis, const Vec &sigmas, 
              APFResampStyle resampTechnique, int pathLength = 0);
              
              
    /**
     * @brief The destructor.
     */
    ~JacEtAlAPF();


    /**
     * @brief Get the prediction covariance matrix after resampling has been performed.
     * @return a Mat representing the covariance matrix of interest.
     */
    Mat getPredictiveVar() const;
    

    /**
     * @brief Evaluate the log of q1. 
     * @param x1 a Vec representing the first time's state.
     * @param y1 a Vec representing the first time's data observation.
     * @return a double evaluation.
     */
    double logQ1Ev (const Vec &x1, const Vec &y1); 
    
    
    /**
     * @brief Evaluate the log of mu.
     * @param x1 a Vec representing the first time's state.
     * @return a double evalution.
     */
    double logMuEv (const Vec &x1); 
    
    
    /**
     * @brief Evaluate the log of g.
     * @param yt a Vec representing time t's data observation.
     * @param xt a Vec representing time t's state.
     * @return a double evaluation.
     */
    double logGEv (const Vec &yt, const Vec &xt); 
    
    
    /**
     * @brief Evaluate the mode of the state transition distribution.
     * @param xtm1 a Vec representing the previous time's state.
     * @return a Vec representing the mode.
     */
    Vec propMu (const Vec &xtm1);
    
    
    /**
     * @brief Samples from f.
     * @param xtm1 the previous time's state.
     * @return a Vec sample of the current time's state.
     */
    Vec  fSamp (const Vec &xtm1);
    
    
    /**
     * @brief Samples from q1.
     * @param y1 a Vec representing time 1's data observation.
     * @return a Vec sample for the first time's state.
     */
    Vec q1Samp (const Vec &y1); 

};

#endif //JACQUIER_ET_AL_APF_H