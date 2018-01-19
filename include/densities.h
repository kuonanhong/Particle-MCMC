#ifndef DENSITIES_H
#define DENSITIES_H

#include <chrono>
#include <Eigen/Dense> //linear algebra stuff
#include <random>

/** Shorthand typedef for Eigen::VectorXd */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;

/** Shorthand typedef for Eigen::MatrixXd */
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

namespace densities
{
////////////////////////////////////////////////
/////////         Constants            /////////
////////////////////////////////////////////////

/** (2 pi)^(-1/2) */
const double inv_sqrt_2pi(0.3989422804014327);

/** (2/pi)^(1/2) */
const double sqrt_two_over_pi(0.797884560802865);

/** log(2pi) */
const double log_two_pi (1.83787706640935);

/** log(2/pi) */
const double log_two_over_pi (-0.451582705289455);


////////////////////////////////////////////////
/////////         Evaluators           /////////
////////////////////////////////////////////////

/**
 * @brief Evaluates the univariate Normal density.
 * @param x the point at which you're evaluating.
 * @param mu the mean.
 * @param sigma the standard deviation.
 * @param log true if you want the log-density. False otherwise.
 * @return a double evaluation.
 */
double evalUnivNorm(const double &x, const double &mu, const double &sigma, bool log = false);

/**
 * @brief Evaluates the multivariate Normal density
 * @param x the point you're evaluating at.
 * @param meanVec the mean vector.
 * @param covMat the positive definite, symmetric covariance matrix.
 * @param log true if you want to return the log density. False otherwise.
 * @return a double evaluation.
 */
double evalMultivNorm(const Vec &x, const Vec &meanVec, const Mat &covMat, bool log = false);


/**
 * @brief Evaluates the multivariate Normal density using the Woodbury Matrix Identity to speed up inversion. 
 * Sigma = A + UCU'. This function assumes A is diagonal and C is symmetric.
 * @param x the point you're evaluating at.
 * @param meanVec the mean vector.
 * @param A  of A + UCU' in vector form because we explicitly make it diagonal.
 * @param U of A + UCU'
 * @param C of A + UCU'
 * @param log true if you want to return the log density. False otherwise.
 * @return a double evaluation.
 */
double evalMultivNormWBDA(const Vec &x, const Vec &meanVec, const Vec &A, const Mat &U, const Mat &C, bool log = false);

/**
 * @brief Evaluates the univariate Beta density
 * @param x the point
 * @param alpha parameter 1 
 * @param beta parameter 2
 * @param log true if you want log density
 * @return double evaluation.
*/       
double evalUnivBeta(const double &x, const double &alpha, const double &beta, bool log = false);

/**
 * @brief Evaluates the univariate Inverse Gamma density
 * @param x the point
 * @param alpha shape parameter  
 * @param beta rate parameter 
 * @param log true if you want log density.
 * @return double evaluation.
*/       
double evalUnivInvGamma(const double &x, const double &alpha, const double &beta, bool log = false);

/**
 * @brief Evaluates the half-normal density
 * @param x the point you're evaluating at
 * @param sigmaSqd the scale parameter
 * @param log true if you want log density.
 * @return double evaluation.
 */
double evalUnivHalfNorm(const double &x, const double &sigmaSqd, bool log = false);


/**
 * @brief Evaluates the logit-Normal distribution (see Wiki for more info)
 * @param x in [0,1] the point you're evaluating at
 * @param mu location parameter that can take any real number
 * @param sigma scale parameter that needs to be positive
 * @param log true if you want to evalute the log-density. False otherwise.
 * @return a double evaluation
 */
double evalLogitNormal(const double &x, const double &mu, const double &sigma, bool log = false);
 

/**
 * @brief Evaluates what I call the "twice-fisher-Normal" distribution
 * https://stats.stackexchange.com/questions/321905/what-is-the-name-of-this-random-variable/321907#321907
 * @param x in [-1,1] the point you are evaluating at
 * @param mu the location parameter (all real numbers)
 * @param sigma the scale parameter (positive)
 * @param log true if you want to evaluate the log-density. False otherwise.
 * @return a double evaluation
 */
double evalTwiceFisherNormal(const double &x, const double &mu, const double &sigma, bool log = false);


/**
 * @brief Evaluates the lognormal density
 * @param x in (0,infty) the point you are evaluating at
 * @param mu the location parameter
 * @param sigma in (0, infty) the scale parameter
 * @param log true if you want to evaluate the log-density. False otherwise.
 * @return a double evaluation
 */
double evalLogNormal(const double &x, const double &mu, const double &sigma, bool log = false);

/**
 * @brief Evaluates the uniform density.
 * @param x in (lower, upper] the point you are evaluating at.
 * @param lower the lower bound of the support for x.
 * @param upper the upper bound for the support of x.
 * @param log true if you want to evaluate the log-density. False otherwise.
 * @return a double evaluation.
 */
double evalUniform(const double &x, const double &lower, const double &upper, bool log = false);

 
////////////////////////////////////////////////
/////////           samplers           /////////
////////////////////////////////////////////////


//! A class that performs sampling from a univariate Normal distribution.
/**
* @class UnivNormSampler
* @author taylor
* @file densities.h
* @brief Samples from univariate Normal distribution.
*/
class UnivNormSampler
{
private:
    std::mt19937 m_rng;    // mt engine
    std::normal_distribution<> m_z_gen;
    double m_mu;
    double m_sigma;
    
public:
    /**
     * @brief Default-constructor sets up for univariate standard Normal random variate generation.
     */
    UnivNormSampler();


     /**
      * @brief The user must supply both mean and covariance matrix.
      * @param mu a double for the mean of the sampling distribution.
      * @param sigma a double (> 0) representing the standard deviation of the samples.
      */
    UnivNormSampler(const double &mu, const double &sigma);


    /**
     * @brief sets the standard deviation of the sampler.
     * @param sigma the desired standard deviation.
     */
    void setStdDev(const double &sigma);
    
    
    /**
     * @brief sets the mean of the sampler.
     * @param mu the desired mean.
     */
    void setMean(const double &mu);
    
        
     /**
      * @brief Draws a random number.
      * @return a random sample of type double.
      */
    double sample();    
};

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
