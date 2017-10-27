#ifndef PMFS_H
#define PMFS_H

#include <random>
#include <vector>
#include <chrono>
#include <Eigen/Dense> //linear algebra stuff

typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

namespace pmfs
{

////////////////////////////////////////////////
/////////         Evaluators           /////////
////////////////////////////////////////////////

/**
 * @brief Evaluates discrete uniform pmf
 * @param x the hypothetical value of a rv 
 * @param k the size of the support i.e. (1,2,...k)
 * @return P(X=x) probability that X equals x
 */
double evalDiscreteUnif(const int &x, const int &k, bool log = false);


////////////////////////////////////////////////
/////////           samplers           /////////
////////////////////////////////////////////////


//! A class that performs sampling from a discrete uniform distribution.
/**
* @class DiscreteUnifSampler
* @author taylor
* @date 08/14/17
* @file pmfs.h
* @brief 
*/
class DiscreteUnifSampler
{
private:
    std::mt19937 m_rng;    // mt engine
    std::uniform_int_distribution<int> m_int_gen;
public:

    //! The constructor.
    /**
     \param k the size of the support i.e. (1,2,...k)
     */
    DiscreteUnifSampler(const int& k);
    
    //! Draws a random vector
    /**
     * \return a random sample of type int 
     */
    int sample();    
};


//! A class that performs sampling from a discrete custom distribution.
/**
* @class DiscreteCustomSampler
* @author taylor
* @date 08/14/17
* @file pmfs.h
* @brief 
*/
class DiscreteCustomSampler
{
private:
    std::mt19937 m_rng;    // mt engine
    std::discrete_distribution<int> m_int_gen;
public:

    //! The constructor.
    /**
     \param k the size of the support i.e. (1,2,...k)
     */
    DiscreteCustomSampler(const std::vector<double>& weights);
    
    //! Draws a random vector
    /**
     * \return a random sample of type int 
     */
    int sample();    
};





} //namespace pmfs

#endif //PMFS
