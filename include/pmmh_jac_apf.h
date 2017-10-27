#ifndef PMMH_JAC_APF_H
#define PMMH_JAC_APF_H

#include "pmmh.h"

//! Estimates the model of Jacquier et al. with a particle marginal Metropolis-Hastings algorithm.
/**
 * @class Pmmh_jac_apf
 * @author taylor
 * @date 14/10/17
 * @file Pmmh_jac_apf.h
 * @brief Estimates the model of Jacquier et al. with a particle marginal Metropolis-Hastings algorithm.
 */
class Pmmh_jac_apf : public Pmmh
{
private:
    // tuning parameters
    unsigned int        m_numParts; 
    unsigned int        m_d;
    double              m_betaVar;
    double              m_fisherPhiVar;
    double              m_muVar; 
    double              m_logSigmaVar; 
    double              m_logRStdDevVar; 
    std::vector<double> m_qSigmaVec;
    
public:

    /**
     * @brief The constructor.
     * @param numParts the number of particles you want.
     * @param startTheta the starting point parameters.
     * @param numMCMCIters the number of MCMC iterations you want.
     * @param dataFile the location of the observed time series data.
     * @param numCols the dimension of your observable time series data.
     */
    Pmmh_jac_apf(unsigned int numParts, std::vector<double> startTheta, unsigned int numMCMCIters, const std::string& dataFile, unsigned int numCols);
    
    
    /**
     * @brief The destructor.
     */
    ~Pmmh_jac_apf();


    /**
     * @brief Convenience function that maps (-1,1) to the reals.
     * @param phi the parameter that needs to stay between -1 and 1
     * @return the evaluation of this function.
     */
    double twiceFisher(const double &phi);
    
    
    /**
     * @brief Convenience function that maps the reals to (-1,1)
     * @param psi an unbounded real number
     * @return phi which is between -1 and 1
     */
    double invTwiceFisher(const double &psi);
    
    
    /**
     * @brief Convenience function that 
     * @param flatOnes the parameters in a std::vector<double> container.
     * @param beta
     * @param phis 
     * @param mus
     * @param sigmas
     * @param RstdDevs
     */
    void flattenParams(std::vector<double> &flatOnes, const Mat &beta, const Vec &phis, const Vec &mus, const Vec &sigmas, const Vec &RstdDevs);
    
    
    /**
     * @brief 
     * @param beta
     * @param phi
     * @param mu
     * @param sigma
     * @param RstdDevs
     * @param flatOnes
     */
    void unFlattenParams(Mat &beta, Vec &phi, Vec &mu, Vec &sigma, Vec &RstdDevs, const std::vector<double> &flatOnes);


    /**
     * @brief Samples new proposal parameters from old parameters.
     * @param oldParams
     * @param newParams
     */
    void qSample(const std::vector<double> &oldParams, std::vector<double> &newParams);

                         
    /**
     * @brief Evaluates the logarithm of the parameter proposal distribution.
     * @param oldParams 
     * @param newParams
     * @return the logarithm of the density evaluation.
     */
    double logQEvaluate(const std::vector<double> &oldParams, 
                                  const std::vector<double> &newParams);


    /**
     * @brief Evaluates the log of the prior.
     * @param theta the parameter inputs.
     * @return the logarithm of the prior evaluation.
     */
    double logPriorEvaluate(const std::vector<double> &theta);


    /**
     * @brief Approximates the log-likelihood with a particle filter. 
     * @param theta the parameters.
     * @param data the observable data.
     * @param out the double to write out to.
     * @return the evaluation (as a double) of the log likelihood approximation.
     */
    double logLikeEvaluate (const std::vector<double>& theta, const std::vector<Vec>& data, std::atomic_bool& cancelled);

};

#endif // PMMH_JAC_APF_H