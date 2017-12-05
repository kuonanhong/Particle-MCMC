#ifndef JACQUIER_ET_AL_H
#define JACQUIER_ET_AL_H

#include "sisr_filter.h"
#include "densities.h"

//! A SISR implementation of the multivariate stochastic volatility model of Jacquier et al.
/**
 * @class JacEtAl
 * @author taylor
 * @date 08/09/17
 * @file jacquier_et_al.h
 * @brief A class that performs filtering for the multivariate stochastic 
 * volatility model of "Stochastic Volatility: Univariate and Multivariate Extensions"
 * by Jacuquier et al
 */
class JacEtAl : public SISRFilter
{
private:
    densities::MVNSampler m_timeOneSampler;
    densities::MVNSampler m_transJumpSampler;
    
    const Mat m_beta;   // factor loadings
    const Vec m_R;      // observatonal covariance matrix diagonals 
    const Vec m_mus;
    const Vec m_phis;   // assuming log-vols have mean 0
    const Vec m_sigmas; // state vars
    int m_dim_obs;
    int m_num_factors;
public:

    /**
     * @brief The constructor
     * @param numParts an integer representing the number of particles.
     * @param beta a Mat representing the factor loadings matrix.
     * @param RAsVec a Vec representing the diagonal elements of the observational covariance matrix.
     * @param mus a Vec representing the mean of the log-volatility processes.
     * @param phis a Vec representing the AR parameters of the log-volatility processes.
     * @param sigmas a Vec representing the standard deviations of the log-volatility processes.
     * @param resampTechnique an enum class representing the type of resampling to be performed.
     * @param pathLength set to the length of your time series data if you want to store entire paths. Otherwise set to 0.
     * @param percEss the percent of the maximum effective sample size threshold for resampling. Not important for several resampTechniques.
     */
    JacEtAl(int numParts, const Mat &beta, const Vec &RAsVec, const Vec &mus, const Vec &phis, const Vec &sigmas, 
              SISRResampStyle resampTechnique, int pathLength = 0, double percEss = 1.0);
              
              
    ~JacEtAl();
    
    /**
     * @brief Samples from q1
     * @param y1
     * @return a Vec sample for x1
     */
    Vec q1Samp(const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of mu.
     * @param x1 a sample for the first time's state.
     * @return a double evaluation.
     */
    double logMuEv (const Vec &x1);


    /**
     * @brief Evaluates the log of q1.
     * @param x1
     * @param y1
     * @return a double evaluation.
     */
    double logQ1Ev (const Vec &x1, const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of g.
     * @param yt
     * @param xt
     * @return a double evaluation.
     */
    double logGEv (const Vec &yt, const Vec &xt);
    
    
    /**
     * @brief Evaluates the log of q (assumed to be constant for all time).
     * @param xt
     * @param xtm1
     * @param yt
     * @return a double evaluation. 
     */
    double logQEv (const Vec &xt, const Vec &xtm1, const Vec &yt);
    
    
    /**
     * @brief Evalutes the log of f.
     * @param xt
     * @param xtm1
     * @return a double evaluation.
     */
    double logFEv (const Vec &xt, const Vec &xtm1);
    
    
    /**
     * @brief Samples from q (assumed to be constant for all time).
     * @param xtm1
     * @param yt
     * @return a Vec sample.
     */
    Vec qSamp (const Vec &xtm1, const Vec &yt);
};

#endif //JACQUIER_ET_AL_H