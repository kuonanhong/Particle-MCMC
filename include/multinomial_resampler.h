#ifndef MULTINOMIALRESAMPLER_H
#define MULTINOMIALRESAMPLER_H

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>

#include "lgssm.h"
#include "fshmm.h"

typedef Eigen::Matrix< double, Eigen::Dynamic, 1> Vec;


//! A class for Multinomial Resampling.
 /**
  * @class MultinomResamp
  * @author taylor
  * @date 07/09/17
  * @file multinomialResampler.h
  * @brief A class used for all the multinomial resampling in any model.
  * Keep in mind this is usually implemented by the base classes themselves, and so
  * this would mean they are not touched by the typical end-user.
  */
class MultinomResamp
{
private:
    std::mt19937 m_gen;
public:
    /**
     * @brief The default constructor. This is the only option available. Sets the seed. 
     */
    MultinomResamp();

    
    /**
     * @brief The default destructor.
     */
    ~MultinomResamp();
    
    
    /**
     * @brief Function to resample in place from weights *not* in log-space.
     * @param oldParts the old particles in form std::vector<std::vector<Eigen::VectorXd> >.
     * @param oldWeights the soon-to-be outdated std::vector of un-normalized weights.
     */
    void resamp(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldWeights);
    

    /**
     * @brief Function to resample in place from SISRs or APFS with weights in log-space.
     * @param oldParts the old particles in form std::vector<std::vector<Eigen::VectorXd> >.
     * @param oldLogUnNormWts the std::vector<double> of log-un-normalized weights.
     */
    void resampLogWts(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldLogUnNormWts);

    /**
      * @brief Function to resample RBPFs--recall you have to resample the samples AND models.
      * @param oldMods embedded kalman filter objects.
      * @param oldSamps samples of each particle.
      * @param oldWts un-normalized weights of each particle.
      */
    void ressampKRBPF(std::vector<Lgssm> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldWts);
    
    
    /**
      * @brief Function to resample for RBPFs--you have to resample the samples AND models.
      * @param oldMods embedded HMM filter objects.
      * @param oldSamps samples of each particle.
      * @param oldWts un-normalized weights of each particle.
      */
    void ressampHRBPF(std::vector<FSHMM> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldWts);


    /**
     * @brief Function to sample random indices. Used inside APFs.
     * @param logFirstStageWeights first-stage weights for first stage of proposal distribution in APF
     * @return A std::vector of integer indices.
     */
    std::vector<int> kGen(const std::vector<double> &logFirstStageWeights);
};

#endif //MULTINOMIALRESAMPLER_H