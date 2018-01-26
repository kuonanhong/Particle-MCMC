#ifndef SISR_FILTER_H
#define SISR_FILTER_H

#include <Eigen/Dense> //linear algebra stuff

#include "mn_resampler.h"

/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

/** enum class for the type of resampling to be performed */
enum class SISRResampStyle {everytime_multinomial, never, ess_multinomial};

//! A base-class for Sequential Importance Sampling with Resampling. Filtering only; no smoothing.
/**
 * @class SISRFilter
 * @author taylor
 * @date 07/09/17
 * @file sisr_filter.h
 * @brief A base-class for Sequential Importance Sampling with Resampling. 
 * Inherit from this if you want to use a SISR for your state space model. Filtering only; no smoothing.
 * @tparam N the number of particles
 */
template <size_t N>
class SISRFilter
{
private:
    // members
    std::array<Vec, N>             m_particles;
    std::array<double, N>          m_logUnNormWeights;
    unsigned int                   m_dimState;
    unsigned int                   m_now;         // time point
    unsigned int                   m_numParts;    // num particles
    double                         m_logLastCondLike; // log p(y_t|y_{1:t-1}) or log p(y1) 
    SISRResampStyle                m_resampTechnique;
    double                         m_ESS; // effective sample size
    double                         m_percentOfNumPartsThresh; // what percent of particles is the lower bound for ESS resampling?
    MNResamp<N>                    m_resampler;
    std::vector<Mat>               m_expectations; // stores any sample averages the user wants
    
    // methods
    void multinomRsmp(std::array<Vec,N> &oldParts, std::array<double,N> &oldLogUnNormWts);

public:

    //!
    /**
     * @brief Constructor 
     * @param resampTechnique  which resampling strategy?
     * @param essPerc ignored unless SISRResampStyle is "ess." What percent of ESS is the threshold for resampling.
     */
    SISRFilter(SISRResampStyle resampTechnique = SISRResampStyle::everytime_multinomial, double essPerc = 1.0);
    
    //! Destructor.
    ~SISRFilter();


    //!
    /**
     * @brief get log p(y_t|y_{1:t-1})$ or log p(y_1).
     * @return The estimate of the most recent log conditional likelihood. 
     */
    double getLogCondLike() const; 
    
    
    //!
    /**
     * @brief get all stored expectations. With respect to $p(x_t|y_{1:t})$
     * @return returns a std::vector<Mat> of all of the approximated expectations.
     */
     std::vector<Mat> getExpectations() const;
    
    
    //!
    /**
     * @brief Get log-un-normalized weights.
     * @return The most recent std::vector of log-un-normalized weights.
     */
    std::array<double,N> getLogUWeights () const; 
    
    
    //!
    /**
     * @brief If storing whole paths, performs smoothing. If not storing whole paths, performs filtering.
     * @param data is a const Vec& representing the current observed value of the time series.
     * @param fs is a vector of functions that operate on each particle Vec. They are used to store empirical expectations (taken with respect to the filtering distribution).
     */
    void filter (const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs = std::vector<std::function<const Mat(const Vec&)> >());

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
     * @brief Calculate fEv or logFEv
     * @param xt is a const Vec& describing the time t state
     * @param xtm1 is a const Vec& describing the time t-1 state
     * @return the density or log-denity evaluation as a double
     */
    virtual double logFEv (const Vec &xt, const Vec &xtm1 ) = 0;

    //!
    /**
     * @brief Sample from the proposal distribution
     * @param xtm1 is a const Vec& describing the time t-1 state
     * @param yt is a const Vec& describing the time t datum
     * @return the sample as a Vec
     */
    virtual Vec qSamp (const Vec &xtm1, const Vec &yt ) = 0;
    
    //!
    /**
     * @brief Calculate qEv or logQEv
     * @param xt is a const Vec& describing the time t state
     * @param xtm1 is a const Vec& describing the time t-1 state
     * @param yt is a const Vec& describing the time t datum
     * @return the density or log-density evaluation as a double
     */
    virtual double logQEv (const Vec &xt, const Vec &xtm1, const Vec &yt ) = 0;    
};


template <size_t N>
SISRFilter<N>::SISRFilter(SISRResampStyle resampTechnique, double essPerc)
                : m_now(0), m_logLastCondLike(0.0), m_resampTechnique(resampTechnique), 
                  m_ESS(m_numParts), m_percentOfNumPartsThresh(essPerc)
{
    m_numParts = m_particles.size();
    std::fill(m_logUnNormWeights.begin(), m_logUnNormWeights.end(), 0.0);
}


template <size_t N>
SISRFilter<N>::~SISRFilter() {}



template <size_t N>
void SISRFilter<N>::filter(const Vec &dat, const std::vector<std::function<const Mat(const Vec&)> >& fs) 
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
        m_logLastCondLike = -std::log(m_numParts) + max + std::log(sumExp);
   
        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double weightNormConst;
        for(auto & h : fs){
            weightNormConst = 0.0;
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[prtcl]) * std::exp(m_logUnNormWeights[prtcl]);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl]);
            }
            m_expectations[fId] /= weightNormConst;
            fId++;
        }
   
        // resample if you should
        if (m_resampTechnique == SISRResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);
   
        // advance time step
        m_now += 1;   
    }
    else // m_now > 0
    {

        // try to iterate over particles all at once
        std::array<Vec, N> newSamps;
        std::array<double,N> oldLogUnNormWts;
        double currentLogWtAdjIndiv;       
        double maxOldLogUnNormWts(m_logUnNormWeights[0]);
        double sumWts(0.0);
        for(size_t ii = 0; ii < m_numParts; ++ii)
        {
            // sample and get weight adjustments
            newSamps[ii] = qSamp(m_particles[ii], dat);
            currentLogWtAdjIndiv = logFEv(newSamps[ii], m_particles[ii]);
            currentLogWtAdjIndiv += logGEv(dat, newSamps[ii]);
            currentLogWtAdjIndiv -= logQEv(newSamps[ii], m_particles[ii], dat);
 
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
        m_expectations.resize(fs.size());
        std::fill(m_expectations.begin(), m_expectations.end(), Vec::Zero(m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double weightNormConst;
        for(auto & h : fs){ // iterate over all functions
            weightNormConst = 0.0;
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ // iterate over all particles
                m_expectations[fId] += h(m_particles[prtcl]) * std::exp(m_logUnNormWeights[prtcl]);
                weightNormConst += std::exp(m_logUnNormWeights[prtcl]);
            }
            m_expectations[fId] /= weightNormConst;
            fId++;
        }
 
        // resample
        if (m_resampTechnique == SISRResampStyle::everytime_multinomial)
            multinomRsmp(m_particles, m_logUnNormWeights);

        // advance time
        m_now += 1;       
    }
}


template <size_t N>
double SISRFilter<N>::getLogCondLike() const
{
    return m_logLastCondLike;
}


template <size_t N>
std::vector<Mat> SISRFilter<N>::getExpectations() const
{
    return m_expectations;
}


template <size_t N>
std::array<double, N> SISRFilter<N>::getLogUWeights() const
{
    return m_logUnNormWeights;
}


template <size_t N>
void SISRFilter<N>::multinomRsmp(std::array<Vec,N> &oldParts, std::array<double,N> &oldLogWeights)
{
    m_resampler.resampLogWts(oldParts, oldLogWeights);
}


#endif // SISR_FILTER_H
