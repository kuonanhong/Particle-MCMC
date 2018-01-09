#include "apf_filter.h"

#include <iostream>
#include <cmath>
// constructors

//double log_sum_exp(const std::vector<double> &arr, const double &someScalar = 0.0) 
//{
//    double maxVal = arr[0] + someScalar;
//    double sum (0);
//
//    for (size_t i = 1 ; i < arr.size() ; i++){
//        if (arr[i] + someScalar > maxVal){
//            maxVal = arr[i] + someScalar;
//        }
//    }
//
//    for (size_t i = 0; i < arr.size() ; i++){
//        sum += std::exp(arr[i] + someScalar - maxVal);
//    }
//    return std::log(sum) + maxVal;
//}


APFFilter::APFFilter(int numParts, APFResampStyle resampTechnique) 
    : m_now(0), m_numParts(numParts), m_logLastCondLike(0.0),
      m_resampTechnique(resampTechnique) //, m_ESS(m_numParts)
{
    // resize particles
    m_particles.resize(m_numParts);
        
    // resize log weights and set them all to 0.0
    m_logUnNormWeights.resize(m_numParts);
    std::fill(m_logUnNormWeights.begin(), m_logUnNormWeights.end(), 0.0);
}

// destructor
APFFilter::~APFFilter(){}


void APFFilter::filter(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs)
{
    
    if(m_now == 0) // first time doesn't change much from SISR
    {
        double max(-1.0/0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)
        {
            // sample particles
            m_particles[ii]  = q1Samp(data);
            m_logUnNormWeights[ii]  = logMuEv(m_particles[ii]);
            m_logUnNormWeights[ii] += logGEv(data, m_particles[ii]);
            m_logUnNormWeights[ii] -= logQ1Ev(m_particles[ii], data);
            
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
        
        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        m_dimState = q1Samp(data).rows();
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double weightNormConst;
        for(auto & h : fs){
            weightNormConst = 0.0;
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[prtcl]) * std::exp(m_logUnNormWeights[prtcl] - max);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl] - max);
            }
            m_expectations[fId] /= weightNormConst;
            fId++;
        }
        
        // resample if you should (automatically normalizes)
        if (m_resampTechnique == APFResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights); 

        // advance time step
        m_now += 1;    
    }
    else{ //m_now > 0
        
        // set up "first stage weights" to make k index sampler 
        std::vector<double> logFirstStageUnNormWeights ( m_logUnNormWeights ); // note that this is unnormalized p_u(k|...) 
        std::vector<Vec> oldPartics ( m_particles );
        double m3(-1.0/0.0);
        double m2(-1.0/0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)  
        {
            // update m3
            if(m_logUnNormWeights[ii] > m3)
                m3 = m_logUnNormWeights[ii];
            
            // sample
            Vec xtm1                        = oldPartics[ii];
            logFirstStageUnNormWeights[ii] += logGEv(data, propMu(xtm1)); // build up first stage weights
            
            // accumulate things
            if(logFirstStageUnNormWeights[ii] > m2)
                m2 = logFirstStageUnNormWeights[ii];

        }
               
        // draw ks (indexes) (handles underflow issues)
        std::vector<unsigned int> myKs(m_numParts);
        kGen(logFirstStageUnNormWeights, myKs); 
                
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
            
            // sampling and unnormalized weight update
            int k                   = myKs[ii];
            Vec xtm1k               = oldPartics[k];
            m_particles[ii]         = fSamp(xtm1k); 
            Vec muT                 = propMu(xtm1k); 
            m_logUnNormWeights[ii] += logGEv(data, m_particles[ii]) - logGEv(data, muT);
            
            // update m1
            if(m_logUnNormWeights[ii] > m1)
                m1 = m_logUnNormWeights[ii];
        }

        // calculate estimate for log of last conditonal likelihood
        for(size_t p = 0; p < m_numParts; ++p)
             first_cll_sum += std::exp( m_logUnNormWeights[p] - m1 );
        m_logLastCondLike = m1 + std::log(first_cll_sum) + m2 + std::log(second_cll_sum) - 2*m3 - 2*std::log(third_cll_sum);

        // calculate expectations before you resample
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double weightNormConst;
        for(auto & h : fs){
            weightNormConst = 0.0;
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[prtcl]) * std::exp(m_logUnNormWeights[prtcl] - m1);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl] - m1);
            }
            m_expectations[fId] /= weightNormConst;
            fId++;
        }
        
        // if you have to resample
        if(m_resampTechnique == APFResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);
            
        // advance time
        m_now += 1; 
    }
}


double APFFilter::getLogCondLike() const
{
    return m_logLastCondLike;
}


std::vector<Mat> APFFilter::getExpectations() const
{
    return m_expectations;
}


std::vector<double> APFFilter::getWeights() const
{
    return m_logUnNormWeights;
}


void APFFilter::kGen(const std::vector<double> &logFirstStageWeights, std::vector<unsigned int> &ks)
{
    m_resampler.kGen(logFirstStageWeights, ks);
}


void APFFilter::multinomRsmp(std::vector<Vec> &oldParts, std::vector<double> &oldLogUnNormWts)
{
    m_resampler.resampLogWts(oldParts, oldLogUnNormWts);
}
