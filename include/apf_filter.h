#ifndef APF_H
#define APF_H

#include <vector>
#include <functional>
#include <Eigen/Dense> 
#include <cmath>

#include "mn_resampler.h"

/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

// TODO-actually implement ess in filterOrSmooth

/** enum class for the type of resampling to be performed */
enum class APFResampStyle {everytime_multinomial, never, ess_multinomial};

//! A base-class for Auxiliary Particle Filtering. Filtering only, no smoothing.
 /**
  * @class APFFilter
  * @author taylor
  * @date 09/09/17
  * @file apf_filter.h
  * @brief A base class for Auxiliary Particle Filtering.
  * Inherit from this if you want to use an APF for your state space model. Filtering only, no smoothing.
  * @tparam N the number of particles
  */
template<size_t N>
class APFFilter
{
private:
    std::array<Vec,N>                m_particles;
    std::array<double,N>             m_logUnNormWeights;
    unsigned int                    m_now;         // time point
    unsigned int                    m_numParts;    // num particles
    unsigned int                    m_dimState;
    double                          m_logLastCondLike; // log p(y_t|y_{1:t-1}) or log p(y1)
    APFResampStyle                  m_resampTechnique;
    MNResamp<N>                     m_resampler; // for resampling particles and sampling indices
    std::vector<Mat>                m_expectations;

    // methods
    void kGen(const std::array<double,N> &logFirstStageWeights, std::array<unsigned int, N> &ks); // for first-stage sampling
    void multinomRsmp(std::array<Vec,N> &oldParts, std::array<double, N> &oldLogUnNormWts); // for final resampling


public:

     /**
      * @brief The constructor.
      * @param resampTechnique the resampling style.
      */
    APFFilter(APFResampStyle resampTechnique = APFResampStyle::everytime_multinomial);

    /**
     * @brief The destructor. 
     */
    ~APFFilter ();
    
    
     /**
      * @brief Get the latest log conditional likelihood.
      * @return a double of the most recent conditional likelihood.
      */
    double getLogCondLike () const; 
    
    
    /**
     * @brief return all stored expectations (taken with respect to $p(x_t|y_{1:t})$
     * @return return a std::vector<Mat> of expectations. How many depends on how many callbacks you gave to 
     */
    std::vector<Mat> getExpectations () const;
    

     /**
      * @brief Use a new datapoint to update the filtering distribution (or smoothing if pathLength > 0).
      * @param data a Vec representing the data
      * @param fs a std::vector of callback functions that are used to calculate expectations with respect to the filtering distribution.
      */
    void filter(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs = std::vector<std::function<const Mat(const Vec&)> >());
    

    /**
     * @brief Evaluates the log of mu.
     * @param x1 a Vec representing time 1's state.
     * @return a double evaluation.
     */
    virtual double logMuEv (const Vec &x1 ) = 0;
    
    
    /**
     * @brief Evaluates the proposal distribution taking a Vec from the previous time's state, and returning a state for the current time.
     * @param xtm1 a Vec representing the previous time's state.
     * @return a Vec representing a likely current time state, to be used by the observation density.
     */
    virtual Vec propMu (const Vec &xtm1 ) = 0;
    
    
    /**
     * @brief Samples from q1.
     * @param y1 a Vec representing time 1's data point.
     * @return a Vec sample for time 1's state.
     */
    virtual Vec q1Samp (const Vec &y1) = 0;
    
    
    /**
     * @brief Samples from f.
     * @param xtm1 a Vec representing the previous time's state.
     * @return a Vec state sample for the current time.
     */
    virtual Vec fSamp (const Vec &xtm1) = 0;
    
    
    /**
     * @brief Evaluates the log of q1.
     * @param x1 a Vec representing time 1's state.
     * @param y1 a Vec representing time 1's data observation.
     * @return a double evaluation.
     */
    virtual double logQ1Ev (const Vec &x1, const Vec &y1) = 0;
    
    
    /**
     * @brief Evaluates the log of g.
     * @param yt a Vec representing time t's data observation.
     * @param xt a Vec representing time t's state.
     * @return a double evaluation.
     */
    virtual double logGEv (const Vec &yt, const Vec &xt) = 0;

};


template<size_t N>
APFFilter<N>::APFFilter(APFResampStyle resampTechnique) 
    : m_now(0), m_logLastCondLike(0.0),
      m_resampTechnique(resampTechnique) //, m_ESS(m_numParts)
{
    std::fill(m_logUnNormWeights.begin(), m_logUnNormWeights.end(), 0.0);
}

template<size_t N>
APFFilter<N>::~APFFilter(){}


template<size_t N>
void APFFilter<N>::filter(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs)
{
    
    if(m_now == 0) // first time doesn't change much from SISR
    {
        double max(-1.0/0.0);
        for(size_t ii = 0; ii < N; ++ii)
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
        for( size_t i = 0; i < N; ++i){
            sumExp += std::exp( m_logUnNormWeights[i] - max );
        }
        m_logLastCondLike = - std::log( static_cast<double>(N) ) + max + std::log(sumExp);
        
        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        m_dimState = q1Samp(data).rows();
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double weightNormConst;
        for(auto & h : fs){
            weightNormConst = 0.0;
            for(size_t prtcl = 0; prtcl < N; ++prtcl){ // iterate over all particles
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
        std::array<double, N> logFirstStageUnNormWeights = m_logUnNormWeights;
        std::array<Vec, N> oldPartics = m_particles;
        double m3(-1.0/0.0);
        double m2(-1.0/0.0);
        for(size_t ii = 0; ii < N; ++ii)  
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
        std::array<unsigned int, N> myKs;
        kGen(logFirstStageUnNormWeights, myKs); 
                
        // now draw xts
        double m1(-1.0/0.0);
        double first_cll_sum(0.0);
        double second_cll_sum(0.0);
        double third_cll_sum(0.0);
        for(size_t ii = 0; ii < N; ++ii)   
        {
            // calclations for log p(y_t|y_{1:t-1}) (using log-sum-exp trick)
            second_cll_sum += std::exp( logFirstStageUnNormWeights[ii] - m2 );
            third_cll_sum  += std::exp( m_logUnNormWeights[ii] - m3 );            
            
            // sampling and unnormalized weight update
            unsigned int k          = myKs[ii];
            Vec xtm1k               = oldPartics[k];
            m_particles[ii]         = fSamp(xtm1k); 
            Vec muT                 = propMu(xtm1k); 
            m_logUnNormWeights[ii] += logGEv(data, m_particles[ii]) - logGEv(data, muT);
            
            // update m1
            if(m_logUnNormWeights[ii] > m1)
                m1 = m_logUnNormWeights[ii];
        }

        // calculate estimate for log of last conditonal likelihood
        for(size_t p = 0; p < N; ++p)
             first_cll_sum += std::exp( m_logUnNormWeights[p] - m1 );
        m_logLastCondLike = m1 + std::log(first_cll_sum) + m2 + std::log(second_cll_sum) - 2*m3 - 2*std::log(third_cll_sum);

        // calculate expectations before you resample
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double weightNormConst;
        for(auto & h : fs){
            weightNormConst = 0.0;
            for(size_t prtcl = 0; prtcl < N; ++prtcl){ // iterate over all particles
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


template<size_t N>
double APFFilter<N>::getLogCondLike() const
{
    return m_logLastCondLike;
}


template<size_t N>
std::vector<Mat> APFFilter<N>::getExpectations() const
{
    return m_expectations;
}


template<size_t N>
void APFFilter<N>::kGen(const std::array<double,N> &logFirstStageWeights, std::array<unsigned int,N> &ks)
{
    m_resampler.kGen(logFirstStageWeights, ks);
}


template<size_t N>
void APFFilter<N>::multinomRsmp(std::array<Vec,N> &oldParts, std::array<double,N> &oldLogUnNormWts)
{
    m_resampler.resampLogWts(oldParts, oldLogUnNormWts);
}

#endif //APF_H