#ifndef SVOL_BS_FILTER_H
#define SVOL_BS_FILTER_H

#include "bootstrap_filter.h"
#include "densities.h"

//! A bootstrap implementation of the classic stochastic volatility model of Taylor.
/**
 * @class SVolBSFilter
 * @author taylor
 * @date 08/09/17
 * @file svol_bs_filter.h
 * @brief Centered parametrization. A class implementing the classic stochastic volatility model of Taylor
 */
class SVolBSFilter: public BSFilter
{
private:
    densities::MVNSampler m_stdNormSampler; // for sampling 
    
    // parameters
    double m_phi;
    double m_beta;
    double m_sigma;
public:

    /**
     * @brief The constructor.
     * @param numParts the number of particles you want.
     * @param beta 
     * @param phi
     * @param sigma
     * @param pLen the length of particles. only set to non-zero number if you want to pre-allocate.
     */
    SVolBSFilter(unsigned numParts, double beta, double phi, double sigma, unsigned pLen=0);
    
    /**
     * @brief The destructor.
     */
    ~SVolBSFilter();
    

    /**
     * @brief Evaluates the log of q1.
     * @param x1 a Vec of time 1 state.
     * @param y1 time one data Vec.
     * @return a double evaluation.
     */
    double logQ1Ev(const Vec &x1, const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of mu.
     * @param x1 Vec of time 1's state.
     * @return a double evaluation.
     */
    double logMuEv(const Vec &x1);
    
    
    /**
     * @brief Evaluates the log of g.
     * @param yt a Vec of time t's data.
     * @param xt a Vec of time t's state.
     * @return a double evaluation.
     */
    double logGEv(const Vec &yt, const Vec &xt);
    
    
    /**
     * @brief Evaluates the log of q.
     * @param xt a Vec of time t's state.
     * @param xtm1 a Vec of time (t-1)'s state.
     * @param yt a Vec pf time t's datum.
     * @return a double evaluation.
     */
    double logQEv(const Vec &xt, const Vec &xtm1, const Vec &yt);
    
    
    /**
     * @brief Samples from q (assumed constant throughout time).
     * @param xtm1 a Vec of the previous time's state.
     * @param yt a Vec of the most recent data observation.
     * @return a Vec sample for time t's state.
     */
    Vec  qSamp(const Vec &xtm1, const Vec &yt);
    
    
    /**
     * @brief Samples from q1.
     * @param y1 a Vec of time 1's data observation.
     * @return a Vec sample for time 1's state.
     */
    Vec  q1Samp(const Vec &y1);
};

#endif // SVOL_BS_FILTER_H
