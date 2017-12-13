#ifndef HMM_RBPF_H
#define HMM_RBPF_H

#include <Eigen/Dense>
#include <vector>
#include <functional> // std::function


#include "fshmm.h"
#include "multinomial_resampler.h" 
#include "densities.h" // for typedefs

/** enum class for the type of resampling to be performed */
enum class HMMRBPFResampStyle {everytime_multinomial, never, ess_multinomial};

//TODO: IMPLEMENT LOG WEIGHTS
// TODO: implement expectations

//! A base-class for HMM Rao-Blackwellized Particle Filtering. 
/*!
 * @class Hmm_Rbpf
 * @author taylor
 * @file hmm_rbpf.h
 * @brief Rao-Blackwellized particle filtering where the innder models are discrete state hmms.
 */
// TODO we might only need one class for this and the kalman one
class Hmm_Rbpf{
private:
    unsigned            m_numParts;
    unsigned            m_now;
    unsigned            m_dimState;
    double              m_lastLogCondLike;
    //double              m_ess;//todo
    //double              m_percentNumPartsThresh; //todo
    std::vector<FSHMM>  m_p_innerMods;
    std::vector<Vec>    m_p_samps;
    std::vector<double> m_logUnNormWeights;
    HMMRBPFResampStyle  m_resampTechnique;
    MultinomResamp      m_resampler;
    std::vector<Mat>    m_expectations;
    
    void resampMultinomHRBPF(std::vector<FSHMM> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldLogWts);
public:

    //! The constructor.
    /**
     * @brief constructor.
     * @param numParts the number of particles you want.
     * @param rt the type of resampling you want to do.
     */
    Hmm_Rbpf(unsigned numParts, HMMRBPFResampStyle rt = HMMRBPFResampStyle::everytime_multinomial);
    
    
    //! Destructor.
    /**
     * @brief destructor.
     */
    ~Hmm_Rbpf();


    //! Filter.
    /**
     * @brief filters everything based on a new data point.
     * @param data the most recent time series observation.
     * @param fs a vector of functions computing logE[h(x_1t, x_2t^i)| x_2t^i,y_1:t]. will access the probability vector of x_1t
     */
    void filter(const Vec &data,
                const std::vector<std::function<const Mat(const Vec &x1tProbs, const Vec &x2t)> >& fs 
                    = std::vector<std::function<const Mat(const Vec&, const Vec&)> >());//, const std::vector<std::function<const Mat(const Vec&)> >& fs);


    //! Get the latest conditional likelihood.
    /**
     * @brief Get the latest conditional likelihood.
     * @return the latest conditional likelihood.
     */
    double getLogCondLike() const;
    
    //!
    /**
     * @brief Get vector of expectations.
     * @return vector of expectations
     */
    std::vector<Mat> getExpectations() const;

    //! Evaluates the first time state density.
    /**
     * @brief evaluates mu.
     * @param x21 component two at time 1
     * @return a double evaluation
     */
    virtual double logMuEv(const Vec &x21) = 0;
    
    
    //! Sample from the first sampler.
    /**
     * @brief samples the second component of the state at time 1.
     * @param y1 most recent datum.
     * @return a Vec sample for x21.
     */
    virtual Vec q1Samp(const Vec &y1) = 0;
    
    
    //! Provides the initial mean vector for each HMM filter object.
    /**
     * @brief provides the initial probability vector for each HMM filter object.
     * @param x21 the second state componenent at time 1.
     * @return a Vec representing the probability of each state element.
     */
    virtual Vec initHMMProbVec(const Vec &x21) = 0;
    
    
    //! Provides the transition matrix for each HMM filter object.
    /**
     * @brief provides the transition matrix for each HMM filter object.
     * @param x21 the second state component at time 1. 
     * @return a transition matrix where element (ij) is the probability of transitioning from state i to state j.
     */
    virtual Mat initHMMTransMat(const Vec &x21) = 0;

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
    virtual double logQ1Ev(const Vec &x21, const Vec &y1) = 0;
    
    
    //! Evaluates the state transition density for the second state component.
    /**
     * @brief Evaluates the state transition density for the second state component.
     * @param x2t the current second state component.
     * @param x2tm1 the previous second state component.
     * @return a double evaluation.
     */
    virtual double logFEv(const Vec &x2t, const Vec &x2tm1) = 0;
    
    
    //! Evaluates the proposal density at time t > 1.
    /**
     * @brief Evaluates the proposal density at time t > 1. 
     * @param x2t the current second state component.
     * @param x2tm1 the previous second state component.
     * @param yt the current time series observation.
     * @return a double evaluation.
     */
    virtual double logQEv(const Vec &x2t, const Vec &x2tm1, const Vec &yt ) = 0;
    
    
    //! How to update your inner HMM filter object at each time.
    /**
     * @brief How to update your inner HMM filter object at each time.
     * @param aModel a HMM filter object describing the conditional closed-form model.
     * @param yt the current time series observation.
     * @param x2t the current second state component.
     */
    virtual void updateFSHMM(FSHMM &aModel, const Vec &yt, const Vec &x2t) = 0;
};

#endif //HMM_RBPF_H