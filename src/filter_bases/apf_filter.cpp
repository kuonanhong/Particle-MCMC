#include "apf_filter.h"

#include <iostream>
#include <cmath>
// constructors

double log_sum_exp(const std::vector<double> &arr, const double &someScalar = 0.0) 
{
    double maxVal = arr[0] + someScalar;
    double sum = 0;

    for (unsigned i = 1 ; i < arr.size() ; i++){
        if (arr[i] + someScalar > maxVal){
            maxVal = arr[i] + someScalar;
        }
    }

    for (unsigned i = 0; i < arr.size() ; i++){
        sum += std::exp(arr[i] + someScalar - maxVal);
    }
    return std::log(sum) + maxVal;
}


APFFilter::APFFilter(int numParts, APFResampStyle resampTechnique, unsigned int pathLength) 
    : m_now(0), m_numParts(numParts), m_logLastCondLike(0.0),
      m_resampTechnique(resampTechnique), m_pathLength(pathLength) //, m_ESS(m_numParts)
{
    
    if (m_pathLength == 0) // only filtering
    {
        m_particles.resize(1);
        m_particles[0].resize(m_numParts);
    }
    else // storing paths
    {
        m_particles.resize(m_pathLength);
        for(unsigned int ii = 0; ii < m_pathLength; ++ii)
            m_particles[ii].resize(m_numParts);
    }
    
    // resize log weights and set them all to 0.0
    m_logUnNormWeights.resize(m_numParts);
    std::vector<double>::iterator it;
    for(it = m_logUnNormWeights.begin(); it != m_logUnNormWeights.end(); ++it)
    {
        *it = 0.0;
    }
}

// destructor
APFFilter::~APFFilter(){}


void APFFilter::filterOrSmooth(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs)
{
    if (m_pathLength == 0){
        filter(data, fs);
    }else{
        smooth(data, fs);
    }
}



void APFFilter::filter(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs)
{
    
    if(m_now == 0) // first time doesn't change much from SISR
    {
        int now (0);    
        for(unsigned int ii = 0; ii < m_numParts; ++ii)
        {
            // sample particles
            m_particles[now][ii]  = q1Samp(data);
            m_logUnNormWeights[ii]  = logMuEv(m_particles[now][ii]);
            m_logUnNormWeights[ii] += logGEv(data, m_particles[now][ii]);
            m_logUnNormWeights[ii] -= logQ1Ev(m_particles[now][ii], data);
        }
        
        // calculate log-likelihood with log-exp-sum trick
        double max = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        double sumExp(0.0);
        for( unsigned int i = 0; i < m_numParts; ++i){
            sumExp += std::exp( m_logUnNormWeights[i] - max );
        }
        m_logLastCondLike = - std::log( static_cast<double>(m_numParts) ) + max + std::log(sumExp);
        
        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        m_dimState = q1Samp(data).rows();
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        for(auto & h : fs){
            double weightNormConst (0.0);
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[now][prtcl]) * std::exp(m_logUnNormWeights[prtcl]);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl]);
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
    
        int now (0);    
        int prevTime (0);
        // set up "first stage weights" to make k index sampler 
        std::vector<double> logFirstStageUnNormWeights ( m_logUnNormWeights );    
        //std::vector<double> logUnNormKProbDistn (m_logUnNormWeights);
        std::vector<std::vector<Vec> > oldPartics ( m_particles );
        double kNormConst (0.0);
        double logQNormConst (0.0);
        for(unsigned int ii = 0; ii < m_numParts; ++ii)  
        {
            
            Vec xtm1                        = oldPartics[prevTime][ii];
            logFirstStageUnNormWeights[ii] += logGEv(data, propMu(xtm1)); // build up first stage weights
            kNormConst += std::exp( m_logUnNormWeights[ii] );
            logQNormConst += std::exp( logFirstStageUnNormWeights[ii] );
            
        }
       
        // finish off qNormConst (divide by kNormConst then take the log)
        logQNormConst = std::log(logQNormConst) - std::log(kNormConst);
        
        // draw k (indexes) 
        std::vector<int> myKs = kGen(logFirstStageUnNormWeights); 
                
        // now draw xts
        //std::vector<Vec> newSamps (m_numParts);
        //std::vector<double> logCurrWeightAdj (m_numParts); 
        std::vector<double> logLikeArrayThing (m_logUnNormWeights);
        for(unsigned int ii = 0; ii < m_numParts; ++ii)   
        {
            
            int k                   = myKs[ii];
            Vec xtm1k               = oldPartics[prevTime][k];
            m_particles[now][ii]    = fSamp(xtm1k); 
            Vec muT                 = propMu(xtm1k); 
            m_logUnNormWeights[ii] += logGEv(data, m_particles[now][ii]) - logGEv(data, muT);
            
            // overwrite and adjust weights and samples
            //m_logUnNormWeights[ii] += logCurrWeightAdj[ii]; 
            //m_particles[now][ii] = newSamps[ii];
            
            // log-like thing (log old + log adjust - log normconst + logQNormConst 
            logLikeArrayThing[ii] = m_logUnNormWeights[ii] - std::log(kNormConst) + logQNormConst;

        }
        
        // finish third piece of logLikeApprox
        //double logLikePiece3 = log_sum_exp( m_logUnNormWeights );
                
        // calculate estimate for log of last conditonal likelihood
        m_logLastCondLike = log_sum_exp(logLikeArrayThing);
        
        // calculate expectations before you resample
        //m_expectations.resize(fs.size());
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        for(auto & h : fs){
            double weightNormConst (0.0);
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[now][prtcl]) * std::exp(m_logUnNormWeights[prtcl]);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl]);
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


void APFFilter::smooth(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs)
{

    
    if(m_now == 0) // first time doesn't change much from SISR
    {
        int now(0);
        for(unsigned int ii = 0; ii < m_numParts; ++ii)
        {
            // sample particles
            m_particles[now][ii]  = q1Samp(data);
            m_logUnNormWeights[ii]  = logMuEv(m_particles[now][ii]);
            m_logUnNormWeights[ii] += logGEv(data, m_particles[now][ii]);
            m_logUnNormWeights[ii] -= logQ1Ev(m_particles[now][ii], data);
        }
        
        // calculate log-likelihood with log-exp-sum trick
        double max = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        double sumExp(0.0);
        for( unsigned int i = 0; i < m_numParts; ++i){
            sumExp += std::exp( m_logUnNormWeights[i] - max );
        }
        m_logLastCondLike = - std::log( static_cast<double>(m_numParts) ) + max + std::log(sumExp);
        
        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        m_dimState = q1Samp(data).rows();
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // fill everything with zero vctors
        int fId(0);
        for(auto & h : fs){
            double weightNormConst (0.0);
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[now][prtcl]) * std::exp(m_logUnNormWeights[prtcl]);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl]);
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
    
        // need this depending on whether we are filtering or (poor man) smoothing
        int now = m_now;    
        int prevTime = m_now - 1;
    
        // set up "first stage weights" to make k index sampler 
        std::vector<double> logFirstStageUnNormWeights ( m_logUnNormWeights );    
        //std::vector<double> logUnNormKProbDistn (m_logUnNormWeights);
        std::vector<std::vector<Vec> > oldPartics ( m_particles );
        double kNormConst (0.0);
        double logQNormConst (0.0);
        for(unsigned int ii = 0; ii < m_numParts; ++ii)  
        {
            
            Vec xtm1                        = oldPartics[prevTime][ii];
            logFirstStageUnNormWeights[ii] += logGEv(data, propMu(xtm1)); // build up first stage weights
            kNormConst += std::exp( m_logUnNormWeights[ii] );
            logQNormConst += std::exp( logFirstStageUnNormWeights[ii] );
            
        }
       
        // finish off qNormConst (divide by kNormConst then take the log)
        logQNormConst = std::log(logQNormConst) - std::log(kNormConst);
        
        // draw k (indexes) 
        std::vector<int> myKs = kGen(logFirstStageUnNormWeights); 
                
        // now draw xts
        //std::vector<Vec> newSamps (m_numParts);
        //std::vector<double> logCurrWeightAdj (m_numParts); 
        std::vector<double> logLikeArrayThing (m_logUnNormWeights);
        for(unsigned int ii = 0; ii < m_numParts; ++ii)   
        {
            
            int k                   = myKs[ii];
            Vec xtm1k               = oldPartics[prevTime][k];
            m_particles[now][ii]    = fSamp(xtm1k); 
            Vec muT                 = propMu(xtm1k); 
            m_logUnNormWeights[ii] += logGEv(data, m_particles[now][ii]) - logGEv(data, muT);
            
            // log-like thing (log old + log adjust - log normconst + logQNormConst 
            logLikeArrayThing[ii] = m_logUnNormWeights[ii] - std::log(kNormConst) + logQNormConst;

        }
        
        // calculate estimate for log of last conditonal likelihood
        m_logLastCondLike = log_sum_exp(logLikeArrayThing);
        
        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // fill everything with zero vctors
        int fId(0);
        for(auto & h : fs){
            double weightNormConst (0.0);
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[now][prtcl]) * std::exp(m_logUnNormWeights[prtcl]);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl]);
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


//double APFFilter::getESS() const
//{
//    return m_ESS;
//}


const std::vector<std::vector<Vec> >& APFFilter::getFullParts() const
{
    return m_particles;
}


const std::vector<Mat>& APFFilter::getExpectations() const
{
    return m_expectations;
}


const std::vector<double>& APFFilter::getWeights() const
{
    return m_logUnNormWeights;
}


std::vector<int> APFFilter::kGen(const std::vector<double> &logFirstStageWeights)
{
    return m_resampler.kGen(logFirstStageWeights);
}


void APFFilter::multinomRsmp(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldLogUnNormWts)
{
    m_resampler.resampLogWts(oldParts, oldLogUnNormWts);
}
