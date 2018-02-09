#ifndef KALMAN_RBPF_BS_H
#define KALMAN_RBPF_BS_H

#include <vector>
#include <Eigen/Dense>

#include "lgssm.h"
#include "multinomial_resampler.h" 

/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

/** enum class for the type of resampling to be performed */
enum class RBPFResampStyle {everytime_multinomial, never, ess_multinomial};

//TODO: IMPLEMENT LOG WEIGHTS


//! A base-class for Kalman Rao-Blackwellized Particle Filtering.
/*!
 * @class Kalman_RBPF_BS
 * @author taylor
 * @file kalman_rbpf.h
 * @brief Rao-Blackwellized Particle Filter with Linear Gaussian submodel
 * @tparam np the number of particles
 */
template <size_t np>
class Kalman_RBPF_BS{
private:
    //only hangs on to filtering distribution, not smoothing
    
    // members
    std::array<Lgssm,np> m_p_innerMods;
    std::array<Vec, np> m_p_samps;
    std::vector<double, np> m_logUnNormWeights;
    unsigned int m_now;
    double m_logLastCondLike; // log p(y_t|y_{1:t-1}) or log p(y1)
    RBPFResampStyle m_resampTechnique;
    double m_ESS;      // effective sample size (not functional yet)
    double m_percentOfNumPartsThresh; // threshold for resampling ess style (not functional yet)
    MNResamp<np> m_resampler; // need for resampling method

    // methods
    void ressampMultinomKRBPF(std::array<Lgssm, np> &oldMods, 
                              std::array<Vec,np> &oldSamps, 
                              std::array<double,np> &oldLogUnNormWts);


public:
    //! The constructor.
    /**
     \param resampTechnique the type of resampling you want to do.
     */
    Kalman_RBPF(RBPFResampStyle resampTechnique = RBPFResampStyle::everytime_multinomial);
    
    //! The desuctor.
    ~Kalman_RBPF();
    
    //! Filter! 
    /**
     * \brief The workhorse function
     * \param data the most recent observable portion of the time series.
     */
    void filter(const Vec &data);

    //! Get the latest log conditional likelihood.
    /**
     * \return the latest log conditional likelihood.
     */
    double getLogCondLike() const; 
    
    
    //! Sample from the first sampler.
    /**
     * @brief samples the second component of the state at time 1.
     * @param y1 most recent datum.
     * @return a Vec sample for x21.
     */
    virtual Vec muSamp(const Vec &y1) = 0;
    
    
    //! Provides the initial mean vector for each Kalman filter object.
    /**
     * @brief provides the initial mean vector for each Kalman filter object.
     * @param x21 the second state componenent at time 1.
     * @return a Vec representing the unconditional mean.
     */
    virtual Vec initKalmanMean(const Vec &x21) = 0;
    
    
    //! Provides the initial covariance matrix for each Kalman filter object.
    /**
     * @brief provides the initial covariance matrix for each Kalman filter object.
     * @param x21 the second state component at time 1.
     * @return a covariance matrix. 
     */
    virtual Mat initKalmanVar(const Vec &x21) = 0;
    
    
    //! Samples the time t second component. 
    /**
     * @brief Samples the time t second component.
     * @param x2tm1 the previous time's second state component.
     * @param yt the current observation.
     * @return a Vec sample of the second state component at the current time.
     */
    virtual Vec fSamp(const Vec &x2tm1, const Vec &yt) = 0;
    
    
    //! How to update your inner Kalman filter object at each time.
    /**
     * @brief How to update your inner Kalman filter object at each time.
     * @param aModel a Kalman filter object describing the conditional closed-form model.
     * @param yt the current time series observation.
     * @param x2t the current second state component.
     */
    virtual void updateKalman(Lgssm &aModel, const Vec &yt, const Vec &x2t) = 0;
};


////////////////////////////////////////////////////////////////////
////////////////////////////////// implementations /////////////////
////////////////////////////////////////////////////////////////////

template <size_t np>
Kalman_RBPF_SISR<np>::Kalman_RBPF_SISR(RBPFResampStyle resampTechnique)
    : m_lastCondLike(1.0), m_resampTechnique(resampTechnique), m_now(0)
{
    std::fill(m_unNormWeights.begin(), m_unNormWeights.end(), 0.0);
}


template <size_t np>
Kalman_RBPF_SISR<np>::~Kalman_RBPF_SISR() 
{
}


template <size_t np>
double Kalman_RBPF_SISR<np>::getLogCondLike() const
{
    return m_lastLogCondLike;
}


template <size_t np>
void Kalman_RBPF_SISR<np>::resampMultinomKRBPF(std::array<Lgssm,np> &oldMods, 
                                               std::array<Vec,np> &oldSamps, 
                                               std::array<double,np> &oldLogUnNormWts)
{
    m_resampler.resampKRBPF(oldMods, oldSamps, oldWts);
}


template <size_t np>
void Kalman_RBPF_SISR<np>::filter(const Vec &data)
{

    if( m_now == 0){ // first data point coming
    
        // initialize and update the closed-form mods        
        Vec tmpMean;
        Mat tmpVar;
        double m1(-1.0/0.0); // maximum log weight adjustment
        for(size_t ii = 0; ii < np; ++ii){
            
            m_p_samps[ii] = q1Samp(data); 
            tmpMean = initKalmanMean(m_p_samps[ii]);
            tmpVar  = initKalmanVar(m_p_samps[ii]);
            m_p_innerMods.emplace_back(tmpMean, tmpVar); // time 1 prior
            updateKalman(m_p_innerMods[ii], data, m_p_samps[ii]);
            m_logUnNormWeights[ii] = m_p_innerMods[ii].getLogCondLike(); 
            
            // maximum to be used in likelihood calc
            if(m_logUnNormWeights[ii] > m1)
                m1 = m_logUnNormWeights[ii];
        }
        
        // calculate log p(y_1)
        double sumexp(0.0);
        for(size_t p = 0; p < np; ++p){
            sumexp += std::exp(m_logUnNormWeights[p] - m1);
        }
        m_logLastCondLike = m1 + std::log(sumexp) - std::log(static_cast<double>(np));
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == RBPFResampStyle::everytime_multinomial)
            resampMultinomKRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
        else if ( (m_resampTechnique == RBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
            resampMultinomKRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
            
        // advance time step
        m_now ++;
    }
    else { //m_now > 0
        
        // update
        Vec newX2Samp;
        double m1(-1.0/0.0); // for the updates weights
        double m2 = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        double sumexpdenom(0.0);
        for(size_t ii = 0; ii < np; ++ii){
            newX2Samp = qSamp(m_p_samps[ii], data);
            updateKalman(m_p_innerMods[ii], data, newX2Samp);
            
            // before you update the weights
            sumexpdenom += std::exp(m_logUnNormWeights[ii] - m2);
                
            m_logUnNormWeights[ii] += m_p_innerMods[ii].getLogCondLike();
                                    
            if ( m_logUnNormWeights[ii] > m1)
                m1 = m_logUnNormWeights[ii];
                        
            m_p_samps[ii] = newX2Samp;
        }
        
        // calculate log p(y_t | y_{1:t-1})
        double sumexpnumer(0.0);
        for(size_t p = 0; p < np; ++p){
            sumexpnumer += std::exp(m_logUnNormWeights[p] - m1);
        }
        m_logLastCondLike = m1 + std::log(sumexpnumer) - m2 - std::log(sumexpdenom);
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == RBPFResampStyle::everytime_multinomial)
            resampMultinomKRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
        else if ( (m_resampTechnique == RBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
            resampMultinomKRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
        
        // update time step
        m_now ++;
    }
    
}
 


#endif //KALMAN_RBPF_BS_H