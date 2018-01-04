#include "sisr_smooth.h"

SISRSmoother::SISRSmoother(int numParts, SISRResampStyle resampTechnique, unsigned int pathLength, double essPerc)
                : m_now(0), m_logLastCondLike(0.0), m_numParts(numParts), m_resampTechnique(resampTechnique), 
                  m_pathLength(pathLength), m_ESS(m_numParts), m_percentOfNumPartsThresh(essPerc)
{
    
    // make sure essPerc is a percent
    assert ((0.0 < essPerc) && (essPerc <= 1.0));
    
    // resize stuff
    m_particles.resize(m_pathLength);
    for(size_t jj = 0; jj < m_pathLength; ++jj)
        m_particles[jj].resize(m_numParts);
    m_logUnNormWeights.resize(m_numParts); 

    //set all the weights to uniform
    std::fill(m_logUnNormWeights.begin(), m_logUnNormWeights.end(), 0.0);
    
}

SISRSmoother::~SISRSmoother() {}


void SISRSmoother::smooth(const Vec &dat) //TODO: no support for ESS stuff
{

    if (m_now == 0) //time 1
    {
    int timeSelector = m_now;
       
        // initialize m_filtMean and m_dimState
        m_dimState = q1Samp(dat).rows();
       
        // only need to iterate over particles once
        //std::vector<double> currentLogWtAdjs(m_numParts);
        double sumWts(0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)
        {
            // sample particles
            m_particles[timeSelector][ii] = q1Samp(dat);
            m_logUnNormWeights[ii] = logMuEv(m_particles[timeSelector][ii]);
            m_logUnNormWeights[ii] += logGEv(dat, m_particles[timeSelector][ii]);
            m_logUnNormWeights[ii] -= logQ1Ev(m_particles[timeSelector][ii], dat);
        }
       
        // calculate log cond likelihood with log-exp-sum trick
        std::vector<double>::iterator idxOfMax = std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        double sumExp(0.0);
        for(size_t i = 0; i < m_numParts; ++i){
            sumExp += std::exp(m_logUnNormWeights[i] - (*idxOfMax));
        }
        m_logLastCondLike = -std::log(m_numParts) + (*idxOfMax) + std::log(sumExp);
   
        // resample if you should
        if (m_resampTechnique == SISRResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);
   
        // advance time step
        m_now += 1;   
    }
    else // m_now > 0
    {
        int timeSelector = m_now;
        int prevTime     = m_now - 1;
       
        double currentLogWtAdjIndiv;
        std::vector<double> oldLogUnNormWts(m_logUnNormWeights);
        double maxOldLogUnNormWts(m_logUnNormWeights[0]);
        double sumWts(0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)
        {
            // sample and get weight adjustments
            m_particles[timeSelector][ii] = qSamp(m_particles[prevTime][ii], dat);
            currentLogWtAdjIndiv = logFEv(m_particles[timeSelector][ii], m_particles[prevTime][ii]);
            currentLogWtAdjIndiv += logGEv(dat, m_particles[timeSelector][ii]);
            currentLogWtAdjIndiv -= logQEv(m_particles[timeSelector][ii], m_particles[prevTime][ii], dat);
 
            // update max of old logUnNormWts
            if (m_logUnNormWeights[ii] > maxOldLogUnNormWts)
                maxOldLogUnNormWts = m_logUnNormWeights[ii];
 
            // overwrite stuff
            m_logUnNormWeights[ii] += currentLogWtAdjIndiv;   

        }
       
        // compute estimate of log p(y_t|y_{1:t-1}) with log-exp-sum trick
        double maxNumer = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end()); //because you added log adjustments
        double sumExp1(0.0);
        double sumExp2(0.0);
        for(size_t i = 0; i < m_numParts; ++i){
            sumExp1 += std::exp(m_logUnNormWeights[i] - maxNumer);
            sumExp2 += std::exp(oldLogUnNormWts[i] - maxOldLogUnNormWts);
        }
        m_logLastCondLike = maxNumer + std::log(sumExp1) - maxOldLogUnNormWts - std::log(sumExp2);
 
        // resample
        if (m_resampTechnique == SISRResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);

        // advance time
        m_now += 1;       
    }
}


double SISRSmoother::getLogCondLike() const
{
    return m_logLastCondLike;
}


double SISRSmoother::getESS() const
{
    return m_ESS;
}


std::vector<std::vector<Vec>> SISRSmoother::getFullParts() const
{
    return m_particles;
}


std::vector<double> SISRSmoother::getLogUWeights() const
{
    return m_logUnNormWeights;
}


void SISRSmoother::multinomRsmp(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldLogWeights)
{
    m_resampler.resampLogWts(oldParts, oldLogWeights);
}
