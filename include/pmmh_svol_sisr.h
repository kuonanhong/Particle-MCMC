#ifndef PMMH_SVOL_SISR_H
#define PMMH_SVOL_SISR_H

#include "pmmh.h"

//! Estimates simple stochastic volatility model with a particle marginal Metropolis-Hastings algorithm.
/**
 * @class Pmmh_svol_sisr
 * @author taylor
 * @file pmmh_svol_sisr.h
 * @brief Estimates simple stochastic volatility model with a particle marginal Metropolis-Hastings algorithm.
 */
class Pmmh_svol_sisr : public Pmmh
{
private:
    // tuning parameters
    unsigned int        m_numParts; 
    unsigned int        m_d;
    std::vector<Mat>    m_qSigma;
    
public:

    /**
     * @brief The constructor.
     * @param numParts the number of particles you want.
     * @param startTheta the starting point parameters.
     * @param qSigma covariance matrix of the random walk proposal distribution.
     * @param numMCMCIters the number of MCMC iterations you want.
     * @param dataFile the location of the observed time series data.
     * @param numCols the dimension of your observable time series data.
     * @param mc true if you want to evaluate likelihood functions in parallel.
     */
    Pmmh_svol_sisr(unsigned numParts, 
                    const std::vector<Vec> &startTheta, 
                    const std::vector<Mat> &qSigma, 
                    unsigned numMCMCIters, 
                    const std::string &dataFile, 
                    unsigned numCols,
                    bool mc);
    
    
    /**
     * @brief The destructor.
     */
    ~Pmmh_svol_sisr();


    /**
     * @brief Samples new proposal parameters from old parameters.
     * @param oldParams
     * @param newParams
     */
    void qSample(const std::vector<Vec> &oldParams, std::vector<Vec> &newParams);

                         
    /**
     * @brief Evaluates the logarithm of the parameter proposal distribution.
     * @param oldParams 
     * @param newParams
     * @return the logarithm of the density evaluation.
     */
    double logQEvaluate(const std::vector<Vec> &oldParams, const std::vector<Vec> &newParams);


    /**
     * @brief Evaluates the log of the prior.
     * @param theta the parameter inputs.
     * @return the logarithm of the prior evaluation.
     */
    double logPriorEvaluate(const std::vector<Vec> &theta);


    /**
     * @brief Approximates the log-likelihood with a particle filter. 
     * @param theta the parameters.
     * @param data the observable data.
     * @param cancelled function will terminate when this is set to true by the class. Must implement mc is true.
     * @return the evaluation (as a double) of the log likelihood approximation.
     */
    double logLikeEvaluate (const std::vector<Vec> &theta, const std::vector<Vec> &data, std::atomic_bool& cancelled);

};

#endif // PMMH_SVOL_SISR_H