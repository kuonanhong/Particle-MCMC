#ifndef BOOTSTRAP_FILTER_H
#define BOOTSTRAP_FILTER_H

#include <Eigen/Dense> //linear algebra stuff

#include "multinomial_resampler.h"

typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

enum class BSResampStyle {everytime_multinomial, never, ess_multinomial};

//! A base-class for Bootstrap Filtering (SISR with naive proposals). This class is a bit faster because it 
// can't store entire particle paths, and it doesn't give you a choice of proposal.
/**
 * @class BSFilter
 * @author taylor
 * @date 07/09/17
 * @file bootstrap_filter.h
 */
class BSFilter
{
private:
    // members
    std::vector<Vec>               m_particles;
    std::vector<double>            m_logUnNormWeights;
    unsigned int                   m_dimState;
    unsigned int                   m_now;         // time point
    unsigned int                   m_numParts;    // num particles
    double                         m_logLastCondLike; // log p(y_t|y_{1:t-1}) or log p(y1) 
    BSResampStyle                  m_resampTechnique;
    double                         m_ESS; // effective sample size
    double                         m_percentOfNumPartsThresh; // what percent of particles is the lower bound for ESS resampling?
    MultinomResamp                 m_resampler;
    std::vector<Mat>               m_expectations; // stores any sample averages the user wants
    
    // methods
    void multinomRsmp(std::vector<Vec> &oldParts, std::vector<double> &oldLogUnNormWts);

public:

    //!
    /**
     * @brief Constructor 
     * @param numParts number of particles
     * @param resampTechnique  which resampling strategy?
     * @param pathLength Set to 0 if you don't want save entire paths. Otherwise, enter time length.
     * @param essPerc ignored unless SISRResampStyle is "ess." What percent of ESS is the threshold for resampling.
     */
    BSFilter(int numParts, BSResampStyle resampTechnique = BSResampStyle::everytime_multinomial, 
                unsigned int pathLength = 0, double essPerc = 1.0);
    
    //! Destructor.
    ~BSFilter();


    //!
    /**
     * @brief get log p(y_t|y_{1:t-1})$ or log p(y_1).
     * @return The estimate of the most recent log conditional likelihood. 
     */
    double getLogCondLike() const; 
    
    
    //!
    /**
     * @brief get effective sample size.
     * @return The current estimate for the effective sample size.
     */
    double getESS() const; 
    
        
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
    std::vector<double> getLogUWeights () const; 
    
    
    //!
    /**
     * @brief If storing whole paths, performs smoothing. If not storing whole paths, performs filtering.
     * @param data is a const Vec& representing the current observed value of the time series.
     * @param fs is a vector of functions that operate on each particle Vec. They are used to store empirical expectations (taken with respect to the filtering distribution).
     */
    void filter(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs = std::vector<std::function<const Mat(const Vec&)> >());  //depending on if m_pathLength == 0 or > 0

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

#endif // SISRFILTER_H
