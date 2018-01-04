#include "apf_smooth.h"

#include <iostream>
#include <cmath>
// constructors


APFSmoother::APFSmoother(int numParts, unsigned int pathLength, APFResampStyle resampTechnique) 
    : m_now(0), m_numParts(numParts), m_logLastCondLike(0.0),
      m_pathLength(pathLength), m_resampTechnique(resampTechnique)  //, m_ESS(m_numParts)
{
    // pre-allocate particle paths
    m_particles.resize(m_pathLength);
    for(size_t ii = 0; ii < m_pathLength; ++ii)
        m_particles[ii].resize(m_numParts);
    
    // resize log weights and set them all to 0.0
    m_logUnNormWeights.resize(m_numParts);
    std::fill(m_logUnNormWeights.begin(), m_logUnNormWeights.end(), 0.0);
}

// destructor
APFSmoother::~APFSmoother(){}


void APFSmoother::smooth(const Vec &data)
{

    // TODO avoid copies more!
    
    if(m_now == 0) // first time doesn't change much from SISR
    {

        int now(0);
        double max(-1.0/0.0);
        
        for(size_t ii = 0; ii < m_numParts; ++ii)
        {
            // sample particles
            m_particles[now][ii]  = q1Samp(data);
            m_logUnNormWeights[ii]  = logMuEv(m_particles[now][ii]);
            m_logUnNormWeights[ii] += logGEv(data, m_particles[now][ii]);
            m_logUnNormWeights[ii] -= logQ1Ev(m_particles[now][ii], data);
            
            // update maximum
            if( m_logUnNormWeights[ii] > max)
                max = m_logUnNormWeights[ii];

        }
        
        // calculate log-likelihood with log-exp-sum trick
        double sumExp(0.0);
        for( size_t i = 0; i < m_numParts; ++i){
            sumExp += std::exp( m_logUnNormWeights[i] - max );
        }
        m_logLastCondLike = - std::log( static_cast<double>(m_numParts) ) + max + std::log(sumExp);
        
        // resample if you should (automatically normalizes)
        if (m_resampTechnique == APFResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights); 

        // advance time step
        m_now += 1;    
    }
    else{ //m_now > 0
        
        // need this depending on whether we are filtering or (poor man) smoothing
        int now = m_now;    
        int prevTime = m_now - 1;
    
        // set up "first stage weights" to make k index sampler 
        std::vector<double> logFirstStageUnNormWeights ( m_logUnNormWeights );    
        std::vector<std::vector<Vec> > oldPartics ( m_particles );
        double m3(-1.0/0.0);
        double m2(-1.0/0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)  
        {

            // update m3
            if(m_logUnNormWeights[ii] > m3)
                m3 = m_logUnNormWeights[ii];

            
            Vec xtm1                        = oldPartics[prevTime][ii];
            logFirstStageUnNormWeights[ii] += logGEv(data, propMu(xtm1)); // build up first stage weights

            // accumulate things
            if(logFirstStageUnNormWeights[ii] > m2)
                m2 = logFirstStageUnNormWeights[ii];
            
        }
        
        // draw k (indexes) 
        std::vector<unsigned int> myKs = kGen(logFirstStageUnNormWeights); 

        // now draw xts
        double m1(-1.0/0.0);
        double first_cll_sum(0.0);
        double second_cll_sum(0.0);
        double third_cll_sum(0.0);        
        for(size_t ii = 0; ii < m_numParts; ++ii)   
        {
            
            // calclations for log p(y_t|y_{1:t-1}) (using log-sum-exp trick)
            second_cll_sum += std::exp( logFirstStageUnNormWeights[ii] - m2 );
            third_cll_sum  += std::exp( m_logUnNormWeights[ii] - m3 );   
            
            // sampling and unnormalized wieght update
            int k                   = myKs[ii];
            Vec xtm1k               = oldPartics[prevTime][k];
            m_particles[now][ii]    = fSamp(xtm1k); 
            Vec muT                 = propMu(xtm1k); 
            m_logUnNormWeights[ii] += logGEv(data, m_particles[now][ii]) - logGEv(data, muT);
            
            // update m1
            if(m_logUnNormWeights[ii] > m1)
                m1 = m_logUnNormWeights[ii];

        }
        
        // calculate estimate for log of last conditonal likelihood
        for(size_t p = 0; p < m_numParts; ++p)
             first_cll_sum += std::exp( m_logUnNormWeights[p] - m1 );
        m_logLastCondLike = m1 + std::log(first_cll_sum) + m2 + std::log(second_cll_sum) - 2*m3 - 2*std::log(third_cll_sum);

        
        // if you have to resample
        if(m_resampTechnique == APFResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);
            
        // advance time
        m_now += 1; 
    }
}


double APFSmoother::getLogCondLike() const
{
    return m_logLastCondLike;
}


std::vector<std::vector<Vec> > APFSmoother::getFullParts() const
{
    return m_particles;
}


std::vector<double> APFSmoother::getWeights() const
{
    return m_logUnNormWeights;
}


std::vector<unsigned int> APFSmoother::kGen(const std::vector<double> &logFirstStageWeights)
{
    return m_resampler.kGen(logFirstStageWeights);
}


void APFSmoother::multinomRsmp(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldLogUnNormWts)
{
    m_resampler.resampLogWts(oldParts, oldLogUnNormWts);
}
