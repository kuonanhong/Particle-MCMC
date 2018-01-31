#ifndef PMMH_H
#define PMMH_H

#include <fstream>
#include <vector>
#include <mutex>
#include <future>
#include <Eigen/Dense> // linear algebra stuff
#include <iostream>
#include "densities.h"
#include "convenience_funcs.h"


/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

//! A base-class for particle marginal Metropolis-Hastings.
 /**
  * @class Pmmh
  * @author taylor
  * @date 10/14/17
  * @file pmmh.h
  * @brief A base class for particle marginal Metropolis-Hastings.
  * Inherit from this if you want to use pmmh to estimate your SSM.
  * The benefits are that doing so will force you to implement certain functions, 
  * and it will abstract away the threaded MH implementation.
  * @tparam np the number of particles
  * @tparam the number of chunks (Vecs) of parameters that are stored.
  */
template <size_t np, size_t pchunks>
class Pmmh{
    
private:    
    std::vector<Vec>    m_data;
    std::array<Vec, pchunks> m_currentTheta;  
    std::ofstream       m_samplesFileStream;
    std::ofstream       m_messageStream;
    unsigned int        m_dimTheta;
    unsigned int        m_numMCMCIters;
    unsigned int        m_numExtraThreads;
    std::mutex          m_outFileMutex;
    bool                m_multicore;
    double              m_numAcceptances;

    
public:

    /**
     * @brief The constructor
     * @param startTheta the initial parameters you want to start sampling from.
     * @param numMCMCIters the number of MCMC iterations you want to do.
     * @param dataFile the location of the observed time series data.
     * @param numCols the dimension of your observable data.
     * @param mc stands for multicore. true or false if you want to use extra cores.
     */
    Pmmh(const std::array<Vec,pchunks> &startTheta, unsigned numMCMCIters, const std::string &dataFile, unsigned numCols, bool mc);


    /**
     * @brief The destructor.
     */
    ~Pmmh();


    /**
     * @brief Call this to begin sampling.
     * @param samplesFile where to store the csv of parameter samples.
     * @param messagesFile where to store the logged messages.
     */
    void commenceSampling(std::string samplesFile, std::string messagesFile);

    
    /**
     * @brief The function that proposes new parameters. 
     * Tip: Make sure to declare your static RNG inside the implementation of this
     * @param oldParams 
     * @param newParams
     */
    virtual void qSample(const std::array<Vec, pchunks> &oldParams, std::array<Vec, pchunks> &newParams) = 0;


    /**
     * @brief Evaluates the logarithm of the proposal density.
     * @param oldParams the old parameters.
     * @param newParams the new parameters. 
     * @return the log of the proposal density.
     */
    virtual double logQEvaluate (const std::array<Vec, pchunks> &oldParams, const std::array<Vec, pchunks> &newParams) = 0;


    /**
     * @brief Evaluates the log of the model's prior distribution.
     * @param theta the parameters argument.
     * @return the log of the prior density.
     */
    virtual double logPriorEvaluate(const std::array<Vec, pchunks> &theta) = 0;


    /**
     * @brief Evaluates (approximates) the log-likelihood with a particle filter.
     * @param theta the parameters with which to run the particle filter.
     * @param data the observed data with which to run the particle filter.
     * @param cancelled is a token you need to provide if doing multithreaded likelihood evals. This allows the function to terminate prematurely.
     * @return the evaluation (as a double) of the log likelihood approximation.
     */
    virtual double logLikeEvaluate (const std::array<Vec, pchunks> &theta, const std::vector<Vec>& data, std::atomic_bool& cancelled) = 0;
  
};


////////////////////////////////////////////////////////////////
//////////////// implementations ///////////////////////////////
////////////////////////////////////////////////////////////////

template <size_t np, size_t pchunks>
Pmmh<np, pchunks>::Pmmh(const std::array<Vec,pchunks> &startTheta, unsigned numMCMCIters, const std::string &dataFile, unsigned numCols, bool mc) : 
        m_currentTheta(startTheta), m_dimTheta(0), m_numMCMCIters(numMCMCIters), m_multicore(mc), m_numAcceptances(0.0)
{
    m_data = convenience_funcs::readInData(dataFile);
    m_numExtraThreads = std::thread::hardware_concurrency() - 1;
    for(size_t i = 0; i < startTheta.size(); ++i)
        m_dimTheta += startTheta[i].rows();
}


template <size_t np, size_t pchunks>
Pmmh<np, pchunks>::~Pmmh()
{
}


//template <size_t np, size_t pchunks>
//void Pmmh<np, pchunks>::commenceSampling(std::string samplesFile, std::string messagesFile){
//    if ( m_multicore ){
//        commence_sampling_mc(samplesFile, messagesFile);
//    }else{
//        commence_sampling_sc(samplesFile, messagesFile);
//    }
//}


//template <size_t np, size_t pchunks>
//void Pmmh<np, pchunks>::commence_sampling_mc(std::string samplesFile, std::string messagesFile)
//{
//    // random number stuff to decide on whether to accept or reject
//    densities::UniformSampler runif; 
//    
//    // these are where we write our results
//    // todo: check for race conditions
//    m_samplesFileStream.open(samplesFile);
//    m_messageStream.open(messagesFile);
//    
//    double oldLogLike (0.0);
//    double oldLogPrior(0.0);
//    unsigned int iter(0);
//    while(iter < m_numMCMCIters) // every iteration
//    {        
//
//        
//        // first iteration no acceptance probability
//        if (iter == 0) { 
//            
//            m_messageStream << "***Iter number: " << 1 << " out of " << m_numMCMCIters << "\n";
//            std::cout << "***Iter number: " << 1 << " out of " << m_numMCMCIters << "\n";        
//        
//            // write accepted parameters to file (initial guesses are always "accepted")
//            convenience_funcs::logParams(m_currentTheta, m_samplesFileStream);
//            
//            // get logLike
//            std::atomic_bool cancel_token(false);
//            oldLogLike = logLikeEvaluate(m_currentTheta, m_data, cancel_token);
//
//            // store prior for next round
//            oldLogPrior = logPriorEvaluate(m_currentTheta);
//            
//            // increase the iteration counter
//            iter++;
//
//        } else { // not the first iteration      
//        
//            // store a few proposed logLikes and logPriors
//            std::vector<std::future<double> > newLogLikes;
//            std::atomic_bool cancel_token(false);
//            std::vector<double> newLogPriors (m_numExtraThreads, 0.0);
//        
//            // propose several new thetas 
//            std::vector<std::array<Vec, pchunks> > proposedThetas(m_numExtraThreads);
//            for(size_t i = 0; i < m_numExtraThreads; ++i)
//                qSample(m_currentTheta, proposedThetas[i]);
//            
//            // store newPrior evaluations and transition kernel evaluations 
//            std::vector<double> logQOldToNews(m_numExtraThreads, 0.0);
//            std::vector<double> logQNewToOlds(m_numExtraThreads, 0.0);
//            for(size_t i = 0; i < m_numExtraThreads; ++i){
//                newLogPriors[i] = logPriorEvaluate(proposedThetas[i]);
//                logQOldToNews[i] = logQEvaluate(m_currentTheta, proposedThetas[i]);
//                logQNewToOlds[i] = logQEvaluate(proposedThetas[i], m_currentTheta);
//                newLogLikes.push_back(std::async(std::launch::async,
//                                                 &Pmmh::logLikeEvaluate,
//                                                 this,
//                                                 std::cref(proposedThetas[i]), 
//                                                 std::cref(m_data),
//                                                 std::ref(cancel_token)));               
//            }
//
//            // accept or reject proposal
//            std::vector<double> logARs(m_numExtraThreads, 0.0);
//            for(size_t i = 0; i < m_numExtraThreads; ++i){
//                
//                // blocks until it is available
//                std::cout << "trying to get() data from core " << i+1 << " out of " << m_numExtraThreads << "\n";
//                m_messageStream << "trying to get() data from core " << i+1 << " out of " << m_numExtraThreads << "\n";
//                double newLL = newLogLikes[i].get();
//                
//                // get acceptance ratio
//                logARs[i] = newLogPriors[i] 
//                            + logQNewToOlds[i] 
//                            + newLL 
//                            - oldLogPrior 
//                            - logQOldToNews[i] 
//                            - oldLogLike;                
//                
//                // output some stuff
//                m_messageStream << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";
//                std::cout << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";        
//                
//                m_messageStream << "acceptance rate: " << m_numAcceptances / iter << " \n";
//                std::cout << "acceptance rate: " << m_numAcceptances / iter << " \n";
//                
//                m_messageStream << "Using core number " << i+1 << " out of " << m_numExtraThreads << "\n";
//                std::cout << "Using core number " << i+1 << " out of " << m_numExtraThreads << "\n";
//                
//                m_messageStream << "AR: " << std::exp(logARs[i]) << "\n";
//                std::cout << "AR: " << std::exp(logARs[i]) << "\n";
//                
//                m_messageStream << "PriorRatio: " << std::exp(newLogPriors[i] - oldLogPrior) << "\n";
//                std::cout << "PriorRatio: " << std::exp(newLogPriors[i] - oldLogPrior) << "\n";
//                
//                m_messageStream << "oldLogLike: " << oldLogLike << "\n";
//                std::cout << "oldLogLike: " << oldLogLike << "\n";
//                
//                m_messageStream << "newLogLike: " << newLL << "\n";
//                std::cout << "newLogLike: " << newLL << "\n";
//                
//                m_messageStream << "LikeRatio: " << std::exp(newLL - oldLogLike) << "\n";
//                std::cout << "LikeRatio: " << std::exp(newLL - oldLogLike) << "\n";
//                        
//                // decide whether to accept or reject
//                double draw = runif.sample();
//                if ( std::isinf(-logARs[i])){
//                    // 0 acceptance rate
//                    std::cout << "rejecting!\n";
//                    // do not change the parameters
//                    // oldPrior stays the same 
//                    // oldLogLike stays the same
//                    iter++; // increase number of iters
//                    m_messageStream << "rejected 100 percent\n";
//                }else if (logARs[i] >= 0.0){
//                    // 100 percent accept 
//                    std::cout << "accepting!\n";
//                    iter++; // increase number of iters
//                    m_currentTheta = proposedThetas[i];
//                    oldLogPrior = newLogPriors[i];
//                    oldLogLike = newLL;
//                    m_messageStream << "accepted 100 percent\n";
//                    cancel_token = true; // cancel remaining threads
//                    m_numAcceptances += 1.0;
//                    break; // stop iterating over threads because the following loglikes will need the new parameters
//                }else if ( std::log(draw) <= logARs[i] ) {
//                    // probabilistic accept
//                    std::cout << "accepting!\n";
//                    iter++; // increase number of iters
//                    m_currentTheta = proposedThetas[i];
//                    oldLogPrior = newLogPriors[i];
//                    oldLogLike = newLL;
//                    m_messageStream << "accepted probabilistically\n";
//                    cancel_token = true; // cancel remaining threads
//                    m_numAcceptances += 1.0;
//                    break; // stop iterating over threads
//                } else if ( std::log(draw) > logARs[i] ) {
//                    std::cout << "rejecting!\n";
//                    // probabilistically reject
//                    // parameters do not change
//                    // oldPrior stays the same 
//                    // oldLogLike stays the same     
//                    iter++; // increase number of iters           
//                    m_messageStream << "rejected probabilistically\n";
//                }else if (std::isnan(logARs[i]) ){ 
//                    // this is unexpected behavior
//                    iter++; // increase number of iters
//                    std::cerr << "there was a NaN. Not accepting proposal. \n";
//                    //std::cerr << "newLogLike: " << newLogLike << "\n";
//                    //std::cerr << "oldLogLikeL " << oldLogLike << "\n";
//                    //std::cerr << "newLogPrior: " << newLogPrior << "\n";
//                    //std::cerr << "oldLogPrior: " << oldLogPrior << "\n";
//                    //std::cerr << "logQNewToOld: "<< logQNewToOld << "\n";
//                    //std::cerr << "logQOldToNew: " << logQOldToNew << "\n";
//                    // does not terminate!
//                    // parameters don't change
//                    // oldPrior stays the same 
//                    // oldLogLike stays the same  
//                    cancel_token = true;
//                    break;
//                } else {
//                    // this case should never be triggered
//                    std::cerr << "you coded your MCMC incorrectly\n";
//                    std::cerr << "stopping...";
//                    iter++; // increase number of iters
//                    cancel_token = true; // cancel remaining threads
//                    break;
//                }
//                
//                // log the theta which may have changedor not
//                m_outFileMutex.lock();
//                convenience_funcs::logParams(m_currentTheta, m_samplesFileStream);
//                m_outFileMutex.unlock();
//                
//            }  // end iteration over threads
//        }
//    }
//    
//    // stop writing thetas and messages
//    m_samplesFileStream.close();
//    m_messageStream.close();
//}


template <size_t np, size_t pchunks>
void Pmmh<np, pchunks>::commenceSampling(std::string samplesFile, std::string messagesFile)
{

    // random number stuff to decide on whether to accept or reject
    densities::UniformSampler runif; 
    
    // these are where we write our results
    // todo: check for race conditions
    m_samplesFileStream.open(samplesFile);
    m_messageStream.open(messagesFile);
    
    double oldLogLike (0.0);
    double oldLogPrior(0.0);
    unsigned int iter(0);
    while(iter < m_numMCMCIters) // every iteration
    {        

        
        // first iteration no acceptance probability
        if (iter == 0) { 
            
            m_messageStream << "***Iter number: " << 1 << " out of " << m_numMCMCIters << "\n";
            std::cout << "***Iter number: " << 1 << " out of " << m_numMCMCIters << "\n";        
        
            // write accepted parameters to file (initial guesses are always "accepted")
            convenience_funcs::logParams(m_currentTheta, m_samplesFileStream);
            
            // get logLike (we use cancel token but it never change
            std::atomic_bool cancel_token(false);
            oldLogLike = logLikeEvaluate(m_currentTheta, m_data, cancel_token);

            // store prior for next round
            oldLogPrior = logPriorEvaluate(m_currentTheta);
            
            // increase the iteration counter
            iter++;
            

        } else { // not the first iteration      
        
                
            // propose several new thetas 
            std::array<Vec, pchunks> proposedTheta;
            qSample(m_currentTheta, proposedTheta);
            
            // store the proposed logPrior, and transition densities
            double newLogPrior = logPriorEvaluate(proposedTheta);
            double logQOldToNew = logQEvaluate(m_currentTheta, proposedTheta);
            double logQNewToOld = logQEvaluate(proposedTheta, m_currentTheta);
            
            // get the likelihood
            std::atomic_bool cancel_token(false);
            double newLL(0.0);
            if (!m_multicore){
                newLL = logLikeEvaluate(proposedTheta, m_data, cancel_token);
            }else{
                std::vector<std::future<double> > newLogLikes;
                for(size_t i = 0; i < m_numExtraThreads; ++i){
                    newLogLikes.push_back(std::async(std::launch::async,
                                                     &Pmmh::logLikeEvaluate,
                                                     this,
                                                     std::cref(proposedTheta), 
                                                     std::cref(m_data),
                                                     std::ref(cancel_token)));               
                }
                for(size_t i = 0; i < m_numExtraThreads; ++i){
                    newLL += newLogLikes[i].get();
                }
                newLL /= m_numExtraThreads;
            }

            // accept or reject proposal
            double logAR = newLogPrior + logQNewToOld + newLL - oldLogPrior - logQOldToNew - oldLogLike;                
                
            // output some stuff
            m_messageStream << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";
            std::cout << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";        
            
            m_messageStream << "acceptance rate: " << m_numAcceptances / iter << " \n";
            std::cout << "acceptance rate: " << m_numAcceptances / iter << " \n";

            m_messageStream << "AR: " << std::exp(logAR) << "\n";
            std::cout << "AR: " << std::exp(logAR) << "\n";

            m_messageStream << "PriorRatio: " << std::exp(newLogPrior - oldLogPrior) << "\n";
            std::cout << "PriorRatio: " << std::exp(newLogPrior - oldLogPrior) << "\n";
            
            m_messageStream << "oldLogLike: " << oldLogLike << "\n";
            std::cout << "oldLogLike: " << oldLogLike << "\n";
            
            m_messageStream << "newLogLike: " << newLL << "\n";
            std::cout << "newLogLike: " << newLL << "\n";
            
            m_messageStream << "LikeRatio: " << std::exp(newLL - oldLogLike) << "\n";
            std::cout << "LikeRatio: " << std::exp(newLL - oldLogLike) << "\n";
            

            // decide whether to accept or reject
            double draw = runif.sample();
            if ( std::isinf(-logAR)){
                // 0 acceptance rate
                // do not change the parameters
                // oldPrior stays the same 
                // oldLogLike stays the same
                iter++; // increase number of iters
                std::cout << "rejected 100 percent\n";
                m_messageStream << "rejected 100 percent\n";
            }else if (logAR >= 0.0){
                // 100 percent accept 
                iter++; // increase number of iters
                m_currentTheta = proposedTheta;
                oldLogPrior = newLogPrior;
                oldLogLike = newLL;
                m_numAcceptances += 1.0;
                std::cout << "accepted 100 percent\n";
                m_messageStream << "accepted 100 percent\n";
            }else if ( std::log(draw) <= logAR ) {
                // probabilistic accept
                iter++; // increase number of iters
                m_currentTheta = proposedTheta;
                oldLogPrior = newLogPrior;
                oldLogLike = newLL;
                m_numAcceptances += 1.0;
                std::cout << "accepted probabilistically\n";
                m_messageStream << "accepted probabilistically\n";
            } else if ( std::log(draw) > logAR ) {
                // probabilistically reject
                // parameters do not change
                // oldPrior stays the same 
                // oldLogLike stays the same     
                iter++; // increase number of iters  
                std::cout << "rejected probabilistically\n";
                m_messageStream << "rejected probabilistically\n";
            }else if (std::isnan(logAR) ){ 
                // this is unexpected behavior
                iter++; // increase number of iters
                std::cerr << "there was a NaN. Not accepting proposal. \n";
                //std::cerr << "newLogLike: " << newLogLike << "\n";
                //std::cerr << "oldLogLikeL " << oldLogLike << "\n";
                //std::cerr << "newLogPrior: " << newLogPrior << "\n";
                //std::cerr << "oldLogPrior: " << oldLogPrior << "\n";
                //std::cerr << "logQNewToOld: "<< logQNewToOld << "\n";
                //std::cerr << "logQOldToNew: " << logQOldToNew << "\n";
                // does not terminate!
                // parameters don't change
                // oldPrior stays the same 
                // oldLogLike stays the same  
            } else {
                // this case should never be triggered
                std::cerr << "you coded your MCMC incorrectly\n";
                std::cerr << "stopping...";
                iter++; // increase number of iters
            }
                
            // log the theta which may have changedor not
            convenience_funcs::logParams(m_currentTheta, m_samplesFileStream);
                
        } // else{

    } //while(iter < m_numMCMCIters)
    
    // stop writing thetas and messages
    m_samplesFileStream.close();
    m_messageStream.close();

    
}// commence_sampling_sc






#endif // PMMH_H