#ifndef KALMAN_RBPF_H
#define KALMAN_RBPF_H

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
 * \class Kalman_RBPF
 * \author taylor
 * \file kalman_rbpf.h
 * \brief Rao-Blackwellized Particle Filter with Linear Gaussian submodel
 */
class Kalman_RBPF{
private:
    //only hangs on to filtering distribution, not smoothing
    
    // members
    std::vector<Lgssm>      m_p_innerMods;
    std::vector<Vec>        m_p_samps;
    std::vector<double>     m_unNormWeights;
    unsigned int            m_numParts;    // num particles
    unsigned int            m_now;
    double                  m_lastCondLike; // p(y_t|y_{1:t-1}) or p(y1)
    RBPFResampStyle         m_resampTechnique;
    double                  m_ESS;      // effective sample size (not functional yet)
    double                  m_percentOfNumPartsThresh; // threshold for resampling ess style (not functional yet)
    MultinomResamp          m_resampler; // need for resampling method

    // methods
    void ressampMultinomKRBPF(std::vector<Lgssm> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldWts);


public:
    //! The constructor.
    /**
     \param numParts the number of particles. 
     \param resampTechnique the type of resampling you want to do.
     */
    Kalman_RBPF(int numParts, RBPFResampStyle resampTechnique = RBPFResampStyle::everytime_multinomial);
    
    //! The desuctor.
    ~Kalman_RBPF();
    
    //! Filter! 
    /**
     * \brief The workhorse function
     * \param data the most recent observable portion of the time series.
     */
    void filter(const Vec &data);

    //! Get the latest conditional likelihood.
    /**
     * \return the latest conditional likelihood.
     */
    double getCondLike() const; 
    
    
    //! Evaluates the first time state density.
    /**
     * @brief evaluates mu.
     * @param x21 component two at time 1
     * @return a double evaluation
     */
    virtual double muEv(const Vec &x21) = 0;
    
    
    //! Sample from the first sampler.
    /**
     * @brief samples the second component of the state at time 1.
     * @param y1 most recent datum.
     * @return a Vec sample for x21.
     */
    virtual Vec q1Samp(const Vec &y1) = 0;
    
    
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
    virtual Vec qSamp(const Vec &x2tm1, const Vec &yt) = 0;
    
    
    //! Evaluates the proposal density of the second state component at time 1.
    /**
     * @brief Evaluates the proposal density of the second state component at time 1.
     * @param x21 the second state component at time 1 you sampled. 
     * @param y1 time 1 observation.
     * @return a double evaluation of the density.
     */
    virtual double q1Ev(const Vec &x21, const Vec &y1) = 0;
    
    
    //! Evaluates the state transition density for the second state component.
    /**
     * @brief Evaluates the state transition density for the second state component.
     * @param x2t the current second state component.
     * @param x2tm1 the previous second state component.
     * @return a double evaluation.
     */
    virtual double fEv(const Vec &x2t, const Vec &x2tm1) = 0;
    
    
    //! Evaluates the proposal density at time t > 1.
    /**
     * @brief Evaluates the proposal density at time t > 1. 
     * @param x2t the current second state component.
     * @param x2tm1 the previous second state component.
     * @param yt the current time series observation.
     * @return a double evaluation.
     */
    virtual double qEv(const Vec &x2t, const Vec &x2tm1, const Vec &yt) = 0;
    
    
    //! How to update your inner Kalman filter object at each time.
    /**
     * @brief How to update your inner Kalman filter object at each time.
     * @param aModel a Kalman filter object describing the conditional closed-form model.
     * @param yt the current time series observation.
     * @param x2t the current second state component.
     */
    virtual void updateKalman(Lgssm &aModel, const Vec &yt, const Vec &x2t) = 0;
};

#endif //KALMAN_RBPF_H