#ifndef KALMAN_RBPF_H
#define KALMAN_RBPF_H

#include <vector>
#include <Eigen/Dense>

#include "lgssm.h"
#include "multinomial_resampler.h" 

typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

enum class RBPFResampStyle {everytime_multinomial, never, ess_multinomial};

//! A base-class for Kalman Rao-Blackwellized Particle Filtering.
/*!
 * \class Kalman_RBPF
 * \author taylor
 * \file Kalman_RBPF.h
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
    
    //! Get the up-to-date un-normalized weights.
    /**
     * @return the un-normalized weights of right now. 
     */
    std::vector<double> getWeights() const;
    
    // pure virtuals...need to define these
    virtual double      muEv          (const Vec                 &x21                                       ) = 0;
    virtual Vec              q1Samp        (const Vec                 &y1                                        ) = 0;
    virtual Vec              initKalmanMean(const Vec                 &x21                                       ) = 0;
    virtual Mat              initKalmanVar (const Vec                 &x21                                       ) = 0;
    virtual Vec              qSamp         (const Vec                 &x2tm1,    const Vec &yt                   ) = 0;
    virtual double      q1Ev          (const Vec                 &x1,       const Vec &y1                   ) = 0;
    virtual double      fEv           (const Vec                 &x2t,      const Vec &x2tm1                ) = 0;
    virtual double      qEv           (const Vec                 &x2t,      const Vec &x2tm1, const Vec &yt ) = 0;
    virtual void             updateKalman  (      Lgssm               &aModel,   const Vec &yt,    const Vec &x2t) = 0;
};

#endif //KALMAN_RBPF_H