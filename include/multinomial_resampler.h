#ifndef MULTINOMIALRESAMPLER_H
#define MULTINOMIALRESAMPLER_H

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <chrono>

#include "lgssm.h"
#include "fshmm.h"

/** typedef for linear algebra stuff */
typedef Eigen::Matrix< double, Eigen::Dynamic, 1> Vec;


//! A class for Multinomial Resampling.
 /**
  * @class MultinomResamp
  * @author taylor
  * @date 07/09/17
  * @file multinomial_resampler.h
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
     * @brief Function to resample in place from SISRs or APFS with weights in log-space.
     * @param oldParts the old particles in form std::vector<std::vector<Eigen::VectorXd> >.
     * @param oldLogUnNormWts the std::vector<double> of log-un-normalized weights.
     */
    void resampLogWts(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldLogUnNormWts);
    

    /**
     * @brief Function to resample in place from SISRs or APFS with weights in log-space.
     * @param oldParts the old particles in form std::vector<Eigen::VectorXd>.
     * @param oldLogUnNormWts the std::vector<double> of log-un-normalized weights.
     */
    void resampLogWts(std::vector<Vec> &oldParts, std::vector<double> &oldLogUnNormWts);

    /**
      * @brief Function to resample RBPFs--recall you have to resample the samples AND models.
      * @param oldMods embedded kalman filter objects.
      * @param oldSamps samples of each particle.
      * @param oldLogWts log un-normalized weights of each particle.
      */
    void resampKRBPF(std::vector<Lgssm> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldLogWts);
    
    
    /**
      * @brief Function to resample for RBPFs--you have to resample the samples AND models.
      * @param oldMods embedded HMM filter objects.
      * @param oldSamps samples of each particle.
      * @param oldLogUnNormWts log un-normalized weights of each particle.
      */
    void resampHRBPF(std::vector<FSHMM> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldLogUnNormWts);


    /**
     * @brief Function to sample random indices. Used inside APFs.
     * @param logFirstStageWeights first-stage weights for first stage of proposal distribution in APF
     * @param ks the "returned" indices of type std::vector<unsigned int>
     */
    void kGen(const std::vector<double> &logFirstStageWeights, std::vector<unsigned int> &ks);
};

#endif //MULTINOMIALRESAMPLER_H