#ifndef MN_RESAMPLER_H
#define MN_RESAMPLER_H

#include <Eigen/Dense>
#include <array>
#include <random>
#include <chrono>

#include "lgssm.h"
#include "fshmm.h"

/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1> Vec;


//! A class template for Multinomial Resampling.
 /**
  * @class MNResamp
  * @author taylor
  * @date 07/09/17
  * @file mn_resampler.h
  * @brief A class used for all the multinomial resampling in any model.
  * Keep in mind this is usually implemented by the base classes themselves, and so
  * this would mean they are not touched by the typical end-user.
  * @tparam N the size of the array aka the number of particles
  */
template<size_t N>
class MNResamp
{
private:
    std::mt19937 m_gen;

public:
    
    /**
     * @brief The default constructor. This is the only option available. Sets the seed. 
     */
    MNResamp();
    
    
    /**
     * @brief The destructor. 
     */
    ~MNResamp();
    
    
    /**
     * @brief Function to resample from log unnormalized weights
     * @param oldParts
     * @param oldLogUnNormWts
     */
    void resampLogWts(std::array<Vec, N> &oldParts, std::array<double, N> &oldLogUnNormWts);


    /**
     * @brief Function to sample random indices. Used inside APFs.
     * @param logFirstStageWeights first-stage weights for first stage of proposal distribution in APF
     * @param ks the "returned" indices of type std::array<unsigned int, N>
     */
    void kGen(const std::array<double,N> &logFirstStageWeights, std::array<unsigned int, N> &ks);
    
    /**
      * @brief Function to resample RBPFs--recall you have to resample the samples AND models.
      * @param oldMods embedded kalman filter objects.
      * @param oldSamps samples of each particle.
      * @param oldLogWts log un-normalized weights of each particle.
      */
    void resampKRBPF(std::array<Lgssm,N> &oldMods, std::array<Vec,N> &oldSamps, std::array<double,N> &oldLogWts);


    /**
      * @brief Function to resample for RBPFs--you have to resample the samples AND models.
      * @param oldMods embedded HMM filter objects.
      * @param oldSamps samples of each particle.
      * @param oldLogUnNormWts log un-normalized weights of each particle.
      */
    void resampHRBPF(std::array<FSHMM,N> &oldMods, std::array<Vec,N> &oldSamps, std::array<double,N> &oldLogUnNormWts);
};




///////////////////////////////////////////////////////////////
/////////////////////// implementations ///////////////////////
///////////////////////////////////////////////////////////////


template<size_t N>
MNResamp<N>::MNResamp() 
        : m_gen{static_cast<std::uint32_t>(
                    std::chrono::high_resolution_clock::now().time_since_epoch().count()
                                           )}
{
}


template<size_t N>
MNResamp<N>::~MNResamp()
{
}


template<size_t N>
void MNResamp<N>::resampLogWts(std::array<Vec, N> &oldParts, std::array<double, N> &oldLogUnNormWts)
{
    // these log weights may be very negative. If that's the case, exponentiating them may cause underflow
    // so we use the "log-exp-sum" trick
    // actually not quite...we just shift the log-weights because after they're exponentiated
    // they have the same normalized probabilities
    
    // get dimensions
    unsigned int numParticles = oldParts.size();
    
    // Create the distribution with exponentiated log-weights
    std::vector<double> w;
    w.resize(oldLogUnNormWts.size());
    double m = *std::max_element(oldLogUnNormWts.begin(), oldLogUnNormWts.end());
    std::transform(oldLogUnNormWts.begin(), oldLogUnNormWts.end(), w.begin(), 
                    [&m](double& d) -> double { return std::exp( d - m ); } );
    std::discrete_distribution<> idxSampler(w.begin(), w.end());
    
    // create temporary particle vector and weight vector
    std::array<Vec, N> tmpPartics = oldParts; // TODO: check this copies!
    
    // sample from the original parts and store in tmpParts
    unsigned int whichPart;
    for(size_t part = 0; part < numParticles; ++part)
    {
        whichPart = idxSampler(m_gen);
        tmpPartics[part] = oldParts[whichPart];
    }
        
    //overwrite olds with news
    oldParts = std::move(tmpPartics);
    std::fill(oldLogUnNormWts.begin(), oldLogUnNormWts.end(), 0.0); // change back    
}


template <size_t N>
void MNResamp<N>::kGen(const std::array<double,N> &logFirstStageWeights, std::array<unsigned int, N> &ks)
{
    // these log weights may be very negative. If that's the case, exponentiating them may cause underflow
    // so we use the "log-exp-sum" trick
    // actually not quite...we just shift the log-weights because after they're exponentiated
    // they have the same normalized probabilities
    
   // Create the distribution with exponentiated log-weights
   // subtract the max first to prevent underflow
   // normalization is taken care of by std::discrete_distribution
    std::array<double, N> w;
    double m = *std::max_element(logFirstStageWeights.begin(), logFirstStageWeights.end());
    std::transform(logFirstStageWeights.begin(), 
                   logFirstStageWeights.end(), 
                   w.begin(), 
                   [&m](const double& d) -> double { return std::exp(d-m); } );
    std::discrete_distribution<> kGenerator(w.begin(), w.end());
    
    // sample ks
    for(size_t i = 0; i < N; ++i){
        ks[i] = kGenerator(m_gen);
    }
}


template <size_t N>
void MNResamp<N>::resampKRBPF(std::array<Lgssm,N> &oldMods, std::array<Vec,N> &oldSamps, std::array<double,N> &oldLogWts)
{
    // Create the distribution with exponentiated log-weights
    std::array<double, N> w;
    double m = *std::max_element(oldLogWts.begin(), oldLogWts.end());
    std::transform(oldLogWts.begin(), oldLogWts.end(), w.begin(), 
                    [&m](double& d) -> double { return std::exp( d - m ); } );
    std::discrete_distribution<> idxSampler(w.begin(), w.end());
    
    // create temporary vectors for samps and mods
    std::array<Vec, N> tmpSamps;
    std::array<Lgssm, N> tmpMods;
    
    // sample from the original parts and store in temporary
    unsigned int whichPart;
    for(size_t part = 0; part < N; ++part)
    {
        whichPart = idxSampler(m_gen);
        tmpSamps[part] = oldSamps[whichPart];
        tmpMods[part] = oldMods[whichPart];
    }
    
    //overwrite olds with news
    oldMods = std::move(tmpMods);
    oldSamps = std::move(tmpSamps);
    std::fill(oldLogWts.begin(), oldLogWts.end(), 0.0);
}


template <size_t N>
void MNResamp<N>::resampHRBPF(std::array<FSHMM,N> &oldMods, std::array<Vec,N> &oldSamps, std::array<double,N> &oldLogUnNormWts)
{
    // Create the distribution with those weights
    std::array<double, N> w;
    double m = *std::max_element(oldLogUnNormWts.begin(), oldLogUnNormWts.end());
    std::transform(oldLogUnNormWts.begin(),
                    oldLogUnNormWts.end(),
                    w.begin(),
                    [&m](double &d)->double{return std::exp(d-m);}
                    );
    std::discrete_distribution<> idxSampler(w.begin(), w.end());
    
    // create temporary vectors for samps and mods
    std::array<Vec,N>   tmpSamps;
    std::array<FSHMM,N> tmpMods;
    
    // sample from the original parts and store in temporary
    int whichPart;
    for(size_t part = 0; part < N; ++part)
    {
        whichPart = idxSampler(m_gen);
        tmpSamps[part] = oldSamps[whichPart];
        tmpMods[part] = oldMods[whichPart];
    }
    
    //overwrite olds with news
    oldMods = std::move(tmpMods);
    oldSamps = std::move(tmpSamps);
    std::fill(oldLogUnNormWts.begin(), oldLogUnNormWts.end(), 0.0);
}



#endif //MN_RESAMPLER_H