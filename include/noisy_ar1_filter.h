#ifndef NOISYAR1FILTER_H
#define NOISYAR1FILTER_H

#include "sisr_filter.h"
#include "densities.h"


//! A SISR implementation of an AR(1) + Noise model.
/**
 * @class NAr1Filter
 * @author taylor
 * @date 08/09/17
 * @file NoisyAr1Filter.h
 * @brief Implements a AR(1) + Noise model.
 */
class NAr1Filter : public SISRFilter
{
private:
    densities::EigenMultivariateNormalSampler m_stdNormSampler; // for sampling 
    
    //parameters
    double m_obsSd;
    double m_stateSd;
    double m_alpha;
    
public:

    /**
     * @brief The constructor.
     * @param numParts an int determining the number of particles to be used.
     * @param osd a double for the observation standard deviation.
     * @param ssd a double for the state process standard deviation.
     * @param a a double for the state ar(1) parameter.
     */
    NAr1Filter(int numParts, double osd, double ssd, double a, 
                SISRResampStyle resampTechnique, int pathLength);
                
                
    /**
     * @brief The destructor.
     */
    ~NAr1Filter();
    

    /**
     * @brief Are either of the standard deviation parameters negative?
     * @return true if either standard deviation is negative (which is bad). Otherwise false.
     */
    bool isSdNeg();
    
    
    /**
     * @brief Samples from q1.
     * @param y1 a Vec of time 1's data.
     * @return a Vec sample of x1.
     */
    Vec q1Samp(const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of mu.
     * @param x1 a Vec of time 1's state.
     * @return a double evaluation.
     */
    double logMuEv  (const Vec &x1);
    
    
    /**
     * @brief Evaluates the log of q1.
     * @param x1 a Vec of time 1's state.
     * @param y1 a Vec of time 1's data.
     * @return  a double evaluation.
     */
    double logQ1Ev (const Vec &x1, const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of g.
     * @param yt a Vec of time t's data.
     * @param xt a Vec of time t's state.
     * @return  a double evaluation.
     */
    double logGEv   (const Vec &yt, const Vec &xt);
    
    
    /**
     * @brief Evaluates the log of q (assumed constant).
     * @param xt a Vec of time t's state.
     * @param xtm1 a Vec of time (t-1)'s state.
     * @param yt a Vec of time t's data.
     * @return a double evaluation.
     */
    double logQEv   (const Vec &xt, const Vec &xtm1, const Vec &yt);
    
    
    /**
     * @brief Evaluates the log of f.
     * @param xt a Vec for time t's state.
     * @param xtm1 a Vec for time (t-1)'s state.
     * @return a double evaluation. 
     */
    double logFEv   (const Vec &xt, const Vec &xtm1);
    
    
    /**
     * @brief Samples from q (assumed constant for all time).
     * @param xtm1 a Vec for the previous time's state.
     * @param yt a Vec for the most recent data.
     * @return a Vec sample of xt.
     */
    Vec qSamp (const Vec &xtm1, const Vec &yt);
};

#endif // NOISYAR1FILTER_H