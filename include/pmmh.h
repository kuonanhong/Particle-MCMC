#ifndef PMMH_H
#define PMMH_H

#include <fstream>
#include <vector>
#include <mutex>
#include <future>
#include <Eigen/Dense> // linear algebra stuff

typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

//! A base-class for particle marginal Metropolis-Hastings.
 /**
  * @class Pmmh
  * @author taylor
  * @date 10/14/17
  * @file Pmmh.h
  * @brief A base class for particle marginal Metropolis-Hastings.
  * Inherit from this if you want to use pmmh to estimate your SSM.
  * The benefits are that doing so will force you to implement certain functions, 
  * and it will abstract away the threaded MH implementation.
  */
class Pmmh{
    
private:    
    std::vector<Vec>    m_data;
    std::vector<double> m_currentTheta; 
    std::ofstream       m_samplesFileStream;
    std::ofstream       m_messageStream;
    unsigned int        m_dimTheta;
    unsigned int        m_numMCMCIters;
    unsigned int        m_numExtraThreads;
    std::mutex          m_outFileMutex;
    bool                m_multicore;

    
    void logParams(const std::vector<double> &thetas, std::ofstream &ofs);
    void readInData(const std::string& fileLoc, unsigned int numCols);
    void commence_sampling_mc(std::string samplesFile, std::string messagesFile);
    void commence_sampling_sc(std::string samplesFile, std::string messagesFile);

    
public:

    /**
     * @brief The constructor
     * @param startTheta the initial parameters you want to start sampling from.
     * @param numMCMCIters the number of MCMC iterations you want to do.
     * @param dataFile the location of the observed time series data.
     * @param numCols the dimension of your observable data.
     * @param mc stands for multicore. true or false if you want to use extra cores.
     */
    Pmmh(std::vector<double> startTheta, unsigned int numMCMCIters, const std::string& dataFile, unsigned int numCols, bool mc);


    /**
     * @brief The destructor.
     */
    ~Pmmh();


    /**
     * @brief Call this to begin sampling.
     * @param samplesFileFile where to store the csv of parameter samples.
     * @param messagesFile where to store the logged messages.
     */
    void commenceSampling(std::string samplesFile, std::string messagesFile);

    
    /**
     * @brief The function that proposes new parameters. 
     * Tip: Make sure to declare your static RNG inside the implementation of this
     * @param oldParams 
     * @param newParams
     */
    virtual void qSample(const std::vector<double> &oldParams,  std::vector<double> &newParams) = 0;


    /**
     * @brief Evaluates the logarithm of the proposal density.
     * @param oldParams the old parameters.
     * @param newParams the new parameters. 
     * @return the log of the proposal density.
     */
    virtual double logQEvaluate (const std::vector<double> &oldParams, const std::vector<double> &newParams) = 0;


    /**
     * @brief Evaluates the log of the model's prior distribution.
     * @param theta the parameters argument.
     * @return the log of the prior density.
     */
    virtual double logPriorEvaluate(const std::vector<double> &theta) = 0;


    /**
     * @brief Evaluates (approximates) the log-likelihood with a particle filter.
     * @param theta the parameters with which to run the particle filter.
     * @param data the observed data with which to run the particle filter.
     * @param cancelled is a token you need to provide if doing multithreaded likelihood evals. This allows the function to terminate prematurely.
     * @return the evaluation (as a double) of the log likelihood approximation.
     */
    virtual double logLikeEvaluate (const std::vector<double>& theta, const std::vector<Vec>& data, std::atomic_bool& cancelled) = 0;
  
};

#endif // PMMH_H