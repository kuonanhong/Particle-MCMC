#include "pmmh.h"

#include <iostream>
#include <future>

#include "densities.h"
#include "convenience_funcs.h"


Pmmh::Pmmh(std::vector<Vec> startTheta, unsigned numMCMCIters, const std::string& dataFile, unsigned int numCols, bool mc) : 
        m_numMCMCIters(numMCMCIters), m_multicore(mc), m_dimTheta(0), m_currentTheta(startTheta)
{
    m_data = convenience_funcs::readInData(dataFile, numCols);
    m_numExtraThreads = std::thread::hardware_concurrency() - 1;
    for(int i = 0; i < startTheta.size(); ++i)
        m_dimTheta += startTheta[i].rows();
}


Pmmh::~Pmmh()
{
}


void Pmmh::commenceSampling(std::string samplesFile, std::string messagesFile){
    if ( m_multicore ){
        commence_sampling_mc(samplesFile, messagesFile);
    }else{
        commence_sampling_sc(samplesFile, messagesFile);
    }
}


void Pmmh::commence_sampling_mc(std::string samplesFile, std::string messagesFile)
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
            
            // get logLike
            std::atomic_bool cancel_token(false);
            oldLogLike = logLikeEvaluate(m_currentTheta, m_data, cancel_token);

            // store prior for next round
            oldLogPrior = logPriorEvaluate(m_currentTheta);
            
            // increase the iteration counter
            iter++;

        } else { // not the first iteration      
        
            // store a few proposed logLikes and logPriors
            std::vector<std::future<double> > newLogLikes;
            std::atomic_bool cancel_token(false);
            std::vector<double> newLogPriors (m_numExtraThreads, 0.0);
        
            // propose several new thetas 
            std::vector<std::vector<Vec> > proposedThetas (m_numExtraThreads, std::vector<Vec>(m_dimTheta));
            for(unsigned int i = 0; i < m_numExtraThreads; ++i)
                qSample(m_currentTheta, proposedThetas[i]);
            
            // store newPrior evaluations and transition kernel evaluations 
            std::vector<double> logQOldToNews(m_numExtraThreads, 0.0);
            std::vector<double> logQNewToOlds(m_numExtraThreads, 0.0);
            for(unsigned int i = 0; i < m_numExtraThreads; ++i){
                newLogPriors[i] = logPriorEvaluate(proposedThetas[i]);
                logQOldToNews[i] = logQEvaluate(m_currentTheta, proposedThetas[i]);
                logQNewToOlds[i] = logQEvaluate(proposedThetas[i], m_currentTheta);
                newLogLikes.push_back(std::async(std::launch::async,
                                                 &Pmmh::logLikeEvaluate,
                                                 this,
                                                 std::cref(proposedThetas[i]), 
                                                 std::cref(m_data),
                                                 std::ref(cancel_token)));               
            }

            // accept or reject proposal
            std::vector<double> logARs(m_numExtraThreads, 0.0);
            for(unsigned int i = 0; i < m_numExtraThreads; ++i){
                
                // blocks until it is available
                std::cout << "trying to get() data from core " << i+1 << " out of " << m_numExtraThreads << "\n";
                double newLL = newLogLikes[i].get();
                
                // get acceptance ratio
                logARs[i] = newLogPriors[i] 
                            + logQNewToOlds[i] 
                            + newLL 
                            - oldLogPrior 
                            - logQOldToNews[i] 
                            - oldLogLike;                
                
                // output some stuff
                m_messageStream << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";
                std::cout << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";        
                m_messageStream << "Using core number " << i+1 << " out of " << m_numExtraThreads << "\n";
                m_messageStream << "AR: " << std::exp(logARs[i]) << "\n";
                m_messageStream << "PriorRatio: " << std::exp(newLogPriors[i] - oldLogPrior) << "\n";
                m_messageStream << "oldLogLike: " << oldLogLike << "\n";
                m_messageStream << "newLogLike: " << newLL << "\n";
                m_messageStream << "LikeRatio: " << std::exp(newLL - oldLogLike) << "\n";
                std::cout << "AR: " << std::exp(logARs[i]) << "\n";
                std::cout << "PriorRatio: " << std::exp(newLogPriors[i] - oldLogPrior) << "\n";
                std::cout << "oldLogLike: " << oldLogLike << "\n";
                std::cout << "newLogLike: " << newLL << "\n";
                std::cout << "LikeRatio: " << std::exp(newLL - oldLogLike) << "\n";
                        
                // decide whether to accept or reject
                double draw = runif.sample();
                if ( std::isinf(-logARs[i])){
                    // 0 acceptance rate
                    std::cout << "rejecting!\n";
                    // do not change the parameters
                    // oldPrior stays the same 
                    // oldLogLike stays the same
                    iter++; // increase number of iters
                    m_messageStream << "rejected 100 percent\n";
                }else if (logARs[i] >= 0.0){
                    // 100 percent accept 
                    std::cout << "accepting!\n";
                    iter++; // increase number of iters
                    m_currentTheta = proposedThetas[i];
                    oldLogPrior = newLogPriors[i];
                    oldLogLike = newLL;
                    m_messageStream << "accepted 100 percent\n";
                    cancel_token = true; // cancel remaining threads
                    break; // stop iterating over threads because the following loglikes will need the new parameters
                }else if ( std::log(draw) <= logARs[i] ) {
                    // probabilistic accept
                    std::cout << "accepting!\n";
                    iter++; // increase number of iters
                    m_currentTheta = proposedThetas[i];
                    oldLogPrior = newLogPriors[i];
                    oldLogLike = newLL;
                    m_messageStream << "accepted probabilistically\n";
                    cancel_token = true; // cancel remaining threads
                    break; // stop iterating over threads
                } else if ( std::log(draw) > logARs[i] ) {
                    std::cout << "rejecting!\n";
                    // probabilistically reject
                    // parameters do not change
                    // oldPrior stays the same 
                    // oldLogLike stays the same     
                    iter++; // increase number of iters           
                    m_messageStream << "rejected probabilistically\n";
                }else if (std::isnan(logARs[i]) ){ 
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
                    cancel_token = true;
                    break;
                } else {
                    // this case should never be triggered
                    std::cerr << "you coded your MCMC incorrectly\n";
                    std::cerr << "stopping...";
                    iter++; // increase number of iters
                    cancel_token = true; // cancel remaining threads
                    break;
                }
                
                // log the theta which may have changedor not
                m_outFileMutex.lock();
                convenience_funcs::logParams(m_currentTheta, m_samplesFileStream);
                m_outFileMutex.unlock();
                
            }  // end iteration over threads
        }
    }
    
    // stop writing thetas and messages
    m_samplesFileStream.close();
    m_messageStream.close();
}

void Pmmh::commence_sampling_sc(std::string samplesFile, std::string messagesFile)
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
            std::vector<Vec> proposedTheta(m_dimTheta);
            qSample(m_currentTheta, proposedTheta);
            
            // store the proposed logLike and logPrior
            std::atomic_bool cancel_token(false);
            double newLL = logLikeEvaluate(proposedTheta, m_data, cancel_token);
            double newLogPrior = logPriorEvaluate(proposedTheta);

            // store newPrior evaluations and transition kernel evaluations 
            double logQOldToNew = logQEvaluate(m_currentTheta, proposedTheta);
            double logQNewToOld = logQEvaluate(proposedTheta, m_currentTheta);

            // accept or reject proposal
            double logAR = newLogPrior + logQNewToOld + newLL - oldLogPrior - logQOldToNew - oldLogLike;                
                
            // output some stuff
            m_messageStream << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";
            std::cout << "***Iter number: " << iter+1 << " out of " << m_numMCMCIters << "\n";        
            
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
                std::cout << "rejecting!\n";
                // do not change the parameters
                // oldPrior stays the same 
                // oldLogLike stays the same
                iter++; // increase number of iters
                m_messageStream << "rejected 100 percent\n";
            }else if (logAR >= 0.0){
                // 100 percent accept 
                std::cout << "accepting!\n";
                iter++; // increase number of iters
                m_currentTheta = proposedTheta;
                oldLogPrior = newLogPrior;
                oldLogLike = newLL;
                m_messageStream << "accepted 100 percent\n";
            }else if ( std::log(draw) <= logAR ) {
                // probabilistic accept
                std::cout << "accepting!\n";
                iter++; // increase number of iters
                m_currentTheta = proposedTheta;
                oldLogPrior = newLogPrior;
                oldLogLike = newLL;
                m_messageStream << "accepted probabilistically\n";
            } else if ( std::log(draw) > logAR ) {
                std::cout << "rejecting!\n";
                // probabilistically reject
                // parameters do not change
                // oldPrior stays the same 
                // oldLogLike stays the same     
                iter++; // increase number of iters           
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



