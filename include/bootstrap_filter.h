#ifndef BOOTSTRAP_FILTER_H
#define BOOTSTRAP_FILTER_H

#include <Eigen/Dense> //linear algebra stuff

#include "mn_resampler.h"

/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

/** enum class for the type of resampling to be performed */
enum class BSResampStyle {everytime_multinomial, never, ess_multinomial};

//! A base class for the boostrap particle filter.
/**
 * @class BSFilter
 * @author taylor
 * @date 26/01/18
 * @file bootstrap_filter.h
 * @brief bootstrap particle filter
 * @tparam N the number of particles
 */
template<size_t N>
class BSFilter
{
private:
    // members
    std::array<Vec, N>             m_particles;
    std::array<double, N>          m_logUnNormWeights;
    unsigned int                   m_dimState;
    unsigned int                   m_now;         // time point
    unsigned int                   m_numParts;    // num particles
    double                         m_logLastCondLike; // log p(y_t|y_{1:t-1}) or log p(y1) 
    BSResampStyle                  m_resampTechnique;
    MNResamp<N>                    m_resampler;
    std::vector<Mat>               m_expectations; // stores any sample averages the user wants
    
    // methods
    void multinomRsmp(std::array<Vec, N> &oldParts, std::array<double, N> &oldLogUnNormWts);

public:

    /**
     * @brief The constructor
     * @param resampTechnique the type of resampling you want to do
     * @param essPerc the effective sample size percent threshold, below which you want to resample.
     */
    BSFilter(BSResampStyle resampTechnique = BSResampStyle::everytime_multinomial, double essPerc = 1.0);


    /**
     * @brief The destructor.
     */
    ~BSFilter();


    /**
     * @brief Returns the most recent (log-) conditiona likelihood.
     * @return log p(y_t | y_{1:t-1})
     */
    double getLogCondLike() const; 
    
    
    /**
     * @brief updates filtering distribution on a new datapoint.
     * @param data the most recent data point
     * @param fs a vector of functions if you want to calculate expectations.
     */
    void filter(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs = std::vector<std::function<const Mat(const Vec&)> >());  //depending on if m_pathLength == 0 or > 0


    /**
     * @brief return all stored expectations (taken with respect to $p(x_t|y_{1:t})$
     * @return return a std::vector<Mat> of expectations. How many depends on how many callbacks you gave to 
     */
    std::vector<Mat> getExpectations () const;
    
    
    //! 
    /**
     * @brief  Calculate muEv or logmuEv
     * @param x1 is a const Vec& describing the state sample
     * @return the density or log-density evaluation as a double
     */
    virtual double logMuEv (const Vec &x1) = 0;

    //!
    /**
     * @brief Samples from time 1 proposal 
     * @param y1 is a const Vec& representing the first observed datum 
     * @return the sample as a Vec
     */
    virtual Vec q1Samp (const Vec &y1) = 0;    
    
    //!
    /**
     * @brief Calculate q1Ev or log q1Ev
     * @param x1 is a const Vec& describing the time 1 state sample
     * @param y1 is a const Vec& describing the time 1 datum
     * @return the density or log-density evaluation as a double
     */
    virtual double logQ1Ev (const Vec &x1, const Vec &y1 ) = 0;
    
    //!
    /**
     * @brief Calculate gEv or logGEv
     * @param yt is a const Vec& describing the time t datum
     * @param xt is a const Vec& describing the time t state
     * @return the density or log-density evaluation as a double
     */
    virtual double logGEv (const Vec &yt, const Vec &xt ) = 0;
    
    
    //!
    /**
     * @brief Sample from the state transition distribution
     * @param xtm1 is a const Vec& describing the time t-1 state
     * @return the sample as a Vec
     */
    virtual Vec fSamp (const Vec &xtm1) = 0;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////// implementations ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

template<size_t N>
BSFilter<N>::BSFilter(BSResampStyle resampTechnique, double essPerc)
                : m_now(0), m_logLastCondLike(0.0), m_resampTechnique(resampTechnique) 
                  
{
    m_numParts = m_logUnNormWeights.size();
    std::fill(m_logUnNormWeights.begin(), m_logUnNormWeights.end(), 0.0);
}


template<size_t N>
BSFilter<N>::~BSFilter() {}


template<size_t N>
void BSFilter<N>::filter(const Vec &dat, const std::vector<std::function<const Mat(const Vec&)> >& fs) //TODO: no support for ESS stuff
{

    if (m_now == 0) //time 1
    {
       
        // initialize m_filtMean and m_dimState
        m_dimState = q1Samp(dat).rows();
       
        // only need to iterate over particles once
        double sumWts(0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)
        {
            // sample particles
            m_particles[ii] = q1Samp(dat);
            m_logUnNormWeights[ii] = logMuEv(m_particles[ii]);
            m_logUnNormWeights[ii] += logGEv(dat, m_particles[ii]);
            m_logUnNormWeights[ii] -= logQ1Ev(m_particles[ii], dat);
        }
       
        // calculate log cond likelihood with log-exp-sum trick
        double max = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        double sumExp(0.0);
        for(size_t i = 0; i < m_numParts; ++i){
            sumExp += std::exp(m_logUnNormWeights[i] - max);
        }
        m_logLastCondLike = -std::log(m_numParts) + (max) + std::log(sumExp);
   
        // calculate expectations before you resample
        // paying mind to underflow
        m_expectations.resize(fs.size());
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        for(auto & h : fs){
            double weightNormConst (0.0);
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[prtcl]) * std::exp( m_logUnNormWeights[prtcl] - (max) );
                weightNormConst += std::exp( m_logUnNormWeights[prtcl] - (max) );
            }
            m_expectations[fId] /= weightNormConst;
            fId++;
        }
   
        // resample if you should
        if (m_resampTechnique == BSResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);
   
        // advance time step
        m_now += 1;   
    }
    else // m_now > 0
    {
       
        // try to iterate over particles all at once
        std::array<Vec, N> newSamps;
        std::array<double, N> oldLogUnNormWts;
        double currentLogWtAdjIndiv;       
        double maxOldLogUnNormWts(m_logUnNormWeights[0]);
        double sumWts(0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)
        {
            // sample and get weight adjustments
            newSamps[ii] = fSamp(m_particles[ii]);
            currentLogWtAdjIndiv = logGEv(dat, newSamps[ii]);
 
            // update max of old logUnNormWts
            if (m_logUnNormWeights[ii] > maxOldLogUnNormWts)
                maxOldLogUnNormWts = m_logUnNormWeights[ii];
 
            // overwrite stuff
            m_logUnNormWeights[ii] += currentLogWtAdjIndiv;
            m_particles[ii] = newSamps[ii];
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

        // calculate expectations before you resample
        // paying mind to underflow
        m_expectations.resize(fs.size());
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double weightNormConst;
        for(auto & h : fs){ // iterate over all functions
            weightNormConst = 0.0;
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[prtcl]) * std::exp(m_logUnNormWeights[prtcl] - maxNumer);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl] - maxNumer);
            }
            m_expectations[fId] /= weightNormConst;
            fId++;
        }

        // resample
        if (m_resampTechnique == BSResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);

        // advance time
        m_now += 1;       
    }
}


template<size_t N>
double BSFilter<N>::getLogCondLike() const
{
    return m_logLastCondLike;
}


template <size_t N>
std::vector<Mat> BSFilter<N>::getExpectations() const
{
    return m_expectations;
}


template<size_t N>
void BSFilter<N>::multinomRsmp(std::array<Vec, N> &oldParts, std::array<double, N> &oldLogUnNormWts)
{
    m_resampler.resampLogWts(oldParts, oldLogUnNormWts);
}



#endif // BOOTSTRAP_FILTER_H
