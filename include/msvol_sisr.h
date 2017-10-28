#ifndef MSVOL_SISR_H
#define MSVOL_SISR_H

#include <vector>

#include "sisr_filter.h"
#include "densities.h"

//! A SISR implementation of the multivariate stochastic volatility model of Pitt and Shephard 1999.
/**
 * @class MSVolSISR
 * @author taylor
 * @date 08/09/17
 * @file msvol_sisr.h
 * @brief Implements the multivariate stochastic volatility model from Pitt and Shephard 1999.
 */
class MSVolSISR : public SISRFilter
{
private:
    densities::MVNSampler m_stdNormSampler; // for sampling 
    
    const Mat m_beta;
    const Vec m_phis; //ordered, factors first before obs
    const Vec m_mus;
    const Vec m_sigmas; //ordered, factors first before obs
    int m_dim_obs;
    int m_num_factors;
public:

    /**
     * @brief The constructor.
     * @param numParts an int describing the number of particles.
     * @param beta a Mat representing the factor loadings matrix.
     * @param phis a Vec representing the AR parameters of the log-volatility processes.
     * @param mus a Vec representing the mean of the log-volatility processes.
     * @param sigmas a Vec representing the standard deviation of the log-volatility processes.
     * @param resampTechnique an enum class that determines the type of resampling to be performed.
     * @param pathLength an int. Set to 0 if jus filtering. To store entire path samples, set to time-length of data.
     * @param percEss the percentage effective sample size resampling threshold.
     */
    MSVolSISR(int numParts, const Mat &beta, const Vec &phis, const Vec &mus, const Vec &sigmas, 
              SISRResampStyle resampTechnique, int pathLength = 0, double percEss = 1.0);
              
    /**
     * @brief The destructor.
     */
    ~MSVolSISR();

    
    /**
     * @brief Samples from q1.
     * @param y1
     * @return a Vec sample of x1
     */
    Vec  q1Samp(const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of mu.
     * @param x1
     * @return a double evaluation.
     */
    double logMuEv  (const Vec &x1);
    
    
    /**
     * @brief Evaluates the log of q1. 
     * @param x1
     * @param y1
     * @return  a double evaluation.
     */
    double logQ1Ev  (const Vec &x1, const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of g.
     * @param yt
     * @param xt
     * @return  a double evaluation.
     */
    double logGEv   (const Vec &yt, const Vec &xt);
    
    
    /**
     * @brief Evaluates the log of q (assumed constant for all time).
     * @param xt
     * @param xtm1
     * @param yt
     * @return  a double evaluation.
     */
    double logQEv   (const Vec &xt, const Vec &xtm1, const Vec &yt);
    
    
    /**
     * @brief Evaluates the log of f.
     * @param xt
     * @param xtm1
     * @return  a double evaluation.
     */
    double logFEv   (const Vec &xt, const Vec &xtm1);
    
    
    /**
     * @brief Samples from q (assumed constant for all time).
     * @param xtm1
     * @param yt
     * @return a Vec samp for xt.
     */
    Vec qSamp (const Vec &xtm1, const Vec &yt);
};

#endif //MSVOL_SISR_H