#ifndef APF_H
#define APF_H

#include <vector>
#include <functional>
#include <Eigen/Dense> 

#include "multinomial_resampler.h"

typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

// TODO-actually implement ess in filterOrSmooth
enum class APFResampStyle {everytime_multinomial, never, ess_multinomial};

//! A base-class for Auxiliary Particle Filtering.
 /**
  * @class APFFilter
  * @author taylor
  * @date 09/09/17
  * @file APFFilter.h
  * @brief A base class for Auxiliary Particle Filtering.
  * Inherit from this if you want to use an APF for your state space model.
  */
class APFFilter
{
private:
    std::vector<std::vector<Vec> >  m_particles;
    std::vector<double>             m_logUnNormWeights;
    unsigned int                    m_now;         // time point
    unsigned int                    m_numParts;    // num particles
    unsigned int                    m_dimState;
    double                          m_logLastCondLike; // log p(y_t|y_{1:t-1}) or log p(y1)
    APFResampStyle                  m_resampTechnique;
    unsigned int                    m_pathLength;
    MultinomResamp                  m_resampler; // for resampling particles and sampling indices
    std::vector<Mat>                m_expectations;

    // methods
    std::vector<unsigned int> kGen(const std::vector<double> &logFirstStageWeights); // for first-stage sampling
    void multinomRsmp(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldLogUnNormWts); // for final resampling
    void filter(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs);
    void smooth(const Vec &data, const std::vector<std::function<const Mat(const Vec&)> >& fs);



public:

     /**
      * @brief The constructor.
      * @param numParts the number of particles you want to use.
      * @param resampTechnique the resampling style.
      * @param pathLength set to 0 if you are filtering. Otherwise, if you wish to retain samples of the entire path, set to time length of data.
      */
    APFFilter (int numParts, APFResampStyle resampTechnique = APFResampStyle::everytime_multinomial, unsigned int pathLength = 0);

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
     * @brief Get the full set of particle paths. Only works if pathLength > 0 in constructor.
     * @return The full set of particle paths as std::vector<std::vector<Vec> >.
     */
    const std::vector<std::vector<Vec> >& getFullParts () const;
    
    
    /**
     * @brief return all stored expectations (taken with respect to $p(x_t|y_{1:t})$
     * @return return a std::vector<Mat> of expectations. How many depends on how many callbacks you gave to 
     */
    const std::vector<Mat>& getExpectations () const;


     /**
      * @brief Get the latest log un-normalized weights.
      * @return a std::vector<double> of the log un norm weights.
      */
    const std::vector<double>& getWeights () const;
    

     /**
      * @brief Use a new datapoint to update the filtering distribution (or smoothing if pathLength > 0).
      * @param dataconst a Vec representing the data
      * @param fs a std::vector of callback functions that are used to calculate expectations with respect to the filtering distribution.
      */
    void filterOrSmooth (const Vec &dataconst, const std::vector<std::function<const Mat(const Vec&)> >& fs = std::vector<std::function<const Mat(const Vec&)> >()); 
    

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

#endif //APF_H