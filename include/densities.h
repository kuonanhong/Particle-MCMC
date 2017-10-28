#ifndef DENSITIES_H
#define DENSITIES_H

#include <chrono>
#include <Eigen/Dense> //linear algebra stuff
#include <random>

typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

namespace densities
{
////////////////////////////////////////////////
/////////         Constants            /////////
////////////////////////////////////////////////

    
const double inv_sqrt_2pi(0.3989422804014327);

const double sqrt_two_over_pi(0.797884560802865);

const double log_two_pi (1.83787706640935);

const double log_two_over_pi (-0.451582705289455);


////////////////////////////////////////////////
/////////         Evaluators           /////////
////////////////////////////////////////////////

/**
 * @brief Evaluates the multivariate normal density
 * @param x 
 * @param meanVec the mean vector
 * @param covMat the positive definite, symmetric covariance matrix
 * @param log true if you want to return the log density. false otherwise
 * @return a double 
 */
double evalMultivNorm(const Vec &x, const Vec &meanVec, const Mat &covMat, bool log = false);

/**
 * @brief Evaluates the univariate Beta density
 * @param x the point
 * @param alpha parameter 1 
 * @param beta parameter 2
 * @return double evaluation.
*/       
double evalUnivBeta(const double &x, const double &alpha, const double &beta, bool log = false);

/**
 * @brief Evaluates the univariate Inverse Gamma density
 * @param x the point
 * @param alpha shape parameter  
 * @param beta rate parameter 
 * @return double evaluation.
*/       
double evalUnivInvGamma(const double &x, const double &alpha, const double &beta, bool log = false);

/**
 * @brief Evaluates the half-normal density
 * @param x the point you're evaluating at
 * @param sigmaSqd the scale parameter
 * @return double evaluation.
 */
double evalUnivHalfNorm(const double &x, const double &sigmaSqd, bool log = false);
 
////////////////////////////////////////////////
/////////           samplers           /////////
////////////////////////////////////////////////

//! A class that performs sampling from a multivariate normal distribution.
/**
* @class MVNSampler
* @author taylor
* @date 05/09/16
* @file densities.h
* @brief Can sample from a distribution with fixed mean and covariance, fixed mean only, fixed covariance only, or nothing fixed.
*/
class MVNSampler
{
private:
    std::mt19937 m_rng;    // mt engine
    std::normal_distribution<> m_z_gen;
    Mat m_scale_mat;
    Vec m_mean;
    
public:
    /**
     * @brief Default-constructor sets up for univariate standard Normal random variate generation.
     */
     MVNSampler();


     /**
      * @brief The user must supply both mean and covariance matrix.
      * @param meanVec a Vec for the mean vector of the sampling distribution.
      * @param covMat a Mat representing the covariance matrix of the samples.
      */
    MVNSampler(const Vec &meanVec, const Mat &covMat);


    /**
     * @brief sets the covariance matrix of the sampler.
     * @param covMat the desired covariance matrix.
     */
    void setCovar(const Mat &covMat);
    
    
    /**
     * @brief sets the mean vector of the sampler.
     * @param meanVec the desired mean vector.
     */
    void setMean(const Vec &meanVec);
    
        
     /**
      * @brief Draws a random vector.
      * @return a Vec random sample.
      */
    Vec sample();    
};


//! A class that performs sampling from a continuous uniform distribution.
/**
* @class UniformSampler
* @author taylor
* @date 04/05/17
* @file densities.h
* @brief 
*/
class UniformSampler
{
private:
    std::mt19937 m_rng;
    std::uniform_real_distribution<> m_unif_gen;
    public:


     /**
      * @brief The default constructor. Gives a lower bound of 0 and upper bound of 1.
      */
    UniformSampler();
    
    
     /**
      * @brief The constructor
      * @param lower the lower bound of the PRNG.
      * @param upper the upper bound of the PRNG.
      */
    UniformSampler(double lower, double upper);
    
    
     /**
      * @brief Draws a sample.
      * @return a sample of type double.
      */
    double sample();
};


} //namespace densities

#endif //DENSITIES_H
