#ifndef SVOLAPFFILTER_H
#define SVOLAPFFILTER_H

#include "apf_filter.h"
#include "densities.h"

//! An APF implementation of the basic stochastic volatility model of Taylor.
/**
 * @class SVolAPFFilter
 * @author taylor
 * @date 09/09/17
 * @file SVolAPFFilter.h
 * @brief Implements the basic stochastic volatility model of Taylor, using an APF.
 */
class SVolAPFFilter : public APFFilter
{
private:
    densities::MVNSampler m_stdNormSampler; // for sampling 
public:

    /**
     * @brief The constructor
     * @param numParts an int describing the number of particles you would like to use.
     */
    SVolAPFFilter(int numParts);
    
    
    /**
     * @brief The desuctor.
     */
    ~SVolAPFFilter();
    

    /**
     * @brief Evaluates the log of q1. 
     * @param x1 a Vec representing time 1's state.
     * @param y1 a Vec represneting time 1's data observation.
     * @return a double evaluation.
     */
    double logQ1Ev(const Vec& x1, const Vec& y1);
    
    
    /**
     * @brief Evaluates the log of mu.
     * @param x1 a Vec representing time 1's state.
     * @return a double evaluation.
     */
    double logMuEv(const Vec &x1);
    
    
    /**
     * @brief Evaluates the log of g.
     * @param yt a Vec representing the current time's data observation.
     * @param xt a Vec representing the current tim'e state.
     * @return a double evaluation.
     */
    double logGEv(const Vec &yt, const Vec &xt);
    
    
    /**
     * @brief Evaluates the log of f.
     * @param xt a Vec representing the current state.
     * @param xtm1 a Vec representing the previous state.
     * @return a double evaluation.
     */
    double logFEv(const Vec &xt, const Vec &xtm1);
    
    
    /**
     * @brief This evaluates the mode of f.
     * @param xtm1 a Vec representing the previous time's state sample.
     * @return a Vec representing the mode.
     */
    Vec propMu(const Vec& xtm1);
    
    
    /**
     * @brief Samples from q1.
     * @param y1 a Vec representing time 1's data observation.
     * @return a Vec sample.
     */
    Vec q1Samp(const Vec &y1);
    
    
    /**
     * @brief Samples from f.
     * @param xtm1 a Vec representing the previous state sample.
     * @return a Vec sample of the current time's state.
     */
    Vec fSamp(const Vec& xtm1);
};

#endif // SVOLFILTER_H
