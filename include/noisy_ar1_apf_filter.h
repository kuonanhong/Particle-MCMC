#ifndef NOISYAR1APFFILTER_H
#define NOISYAR1APFFILTER_H

#include "apf_filter.h"
#include "densities.h"


//! An APF implementation of a AR(1) + Noise Model. 
/**
 * @class NAr1APFFilter
 * @author taylor
 * @date 10/09/17
 * @file NoisyAr1APFFilter.h
 * @brief An APF implementation of a AR(1) + Noise Model. 
 */
class NAr1APFFilter : public APFFilter
{
private:
    densities::MVNSampler m_stdNormSampler; 
    // parameters
    double m_alpha;
    double m_stateSd;
    double m_obsSd;
    
public:

    /**
     * @brief The constructor.
     * @param numParts the number of particles.
     * @param osd the observational standard deviation.
     * @param ssd the state transition standard deviation.
     * @param a the state autoregressive parameter.
     * @param resampTechnique the resampling style to be used.
     * @param pathLength set to 0 if filtering. Otherwise, if you want to retain entire path samples, set to time length of data.
     */
    NAr1APFFilter(int numParts, double osd, double ssd, double a, 
                    APFResampStyle resampTechnique = APFResampStyle::everytime_multinomial,
                    int pathLength = 0);
                    
                    
    /**
     * @brief The destructor.
     */
    ~NAr1APFFilter();

    
    /**
     * @brief Evaluates the log of q1. 
     * @param x1 time 1's state.
     * @param y1 time 1's data observation.
     * @return a double evaluation.
     */
    double logQ1Ev (const Vec &x1, const Vec &y1);
    
    
    /**
     * @brief Evaluates the log of mu.
     * @param x1 time 1's state.
     * @return a double evaluation.
     */
    double logMuEv  (const Vec &x1);
    
    
    /**
     * @brief Evaluats the log of G.
     * @param yt the current time's data observation.
     * @param xt the current time's state.
     * @return a double evaluation.
     */
    double logGEv (const Vec &yt, const Vec &xt);
    
    
    /**
     * @brief Evaluates the mode of the state transition density.
     * @param xtm1 the previous time's state.
     * @return the mode.
     */
    Vec propMu(const Vec &xtm1);
    
    
    /**
     * @brief Samples from the state transition density.
     * @param xtm1 the previous time's state.
     * @return a Vec sample of the current time's state.
     */
    Vec fSamp (const Vec &xtm1);
    
    
    /**
     * @brief Samples from q1.
     * @param y1 time 1's data observation.
     * @return a Vec sample for time 1's state.
     */
    Vec q1Samp(const Vec &y1);
};

#endif // NOISYAR1APFFILTER_H