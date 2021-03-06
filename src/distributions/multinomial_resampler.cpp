#include "multinomial_resampler.h"

#include <assert.h>
#include <iostream>

MultinomResamp::MultinomResamp() 
        : m_gen{static_cast<std::uint32_t>(
                    std::chrono::high_resolution_clock::now().time_since_epoch().count()
                                           )}
{
}


MultinomResamp::~MultinomResamp()
{
}


void MultinomResamp::resampLogWts(std::vector<std::vector<Vec> > &oldParts, std::vector<double> &oldLogUnNormWts)
{
    // these log weights may be very negative. If that's the case, exponentiating them may cause underflow
    // so we use the "log-exp-sum" trick
    // actually not quite...we just shift the log-weights because after they're exponentiated
    // they have the same normalized probabilities
    
    // get dimensions
    unsigned int timeLength = oldParts.size();
    unsigned int numParticles = oldParts[0].size();
    
    // Create the distribution with exponentiated log-weights
    std::vector<double> w;
    w.resize(oldLogUnNormWts.size());
    double m = *std::max_element(oldLogUnNormWts.begin(), oldLogUnNormWts.end());
    std::transform(oldLogUnNormWts.begin(), oldLogUnNormWts.end(), w.begin(), 
                    [&m](double& d) -> double { return std::exp( d - m ); } );
    std::discrete_distribution<> idxSampler(w.begin(), w.end());
    
    // create temporary particle vector and weight vector
    std::vector<std::vector<Vec> > tmpPartics(timeLength, std::vector<Vec>(numParticles));
    
    // sample from the original parts and store in tmpParts
    int whichPart;
    for(size_t part = 0; part < numParticles; ++part)
    {
        whichPart = idxSampler(m_gen);
        for(size_t time = 0; time < timeLength; ++time)
        {
            tmpPartics[time][part] = oldParts[time][whichPart];
        }
    }
        
    //overwrite olds with news
    oldParts = std::move(tmpPartics);
    std::fill (oldLogUnNormWts.begin(), oldLogUnNormWts.end(), 0.0); // change back    
}


void MultinomResamp::resampLogWts(std::vector<Vec> &oldParts, std::vector<double> &oldLogUnNormWts)
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
    std::vector<Vec> tmpPartics = oldParts;
    
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


void MultinomResamp::resampKRBPF(std::vector<Lgssm> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldLogWts)
{
    // get dimensions
    int numParticles = oldLogWts.size();
    
    // Create the distribution with exponentiated log-weights
    std::vector<double> w;
    w.resize(oldLogWts.size());
    double m = *std::max_element(oldLogWts.begin(), oldLogWts.end());
    std::transform(oldLogWts.begin(), oldLogWts.end(), w.begin(), 
                    [&m](double& d) -> double { return std::exp( d - m ); } );
    std::discrete_distribution<> idxSampler(w.begin(), w.end());
    
    // create temporary vectors for samps and mods
    std::vector<Vec>   tmpSamps(numParticles);
    std::vector<Lgssm> tmpMods(numParticles);
    
    // sample from the original parts and store in temporary
    unsigned int whichPart;
    for(size_t part = 0; part < numParticles; ++part)
    {
        whichPart = idxSampler(m_gen);
        tmpSamps[part] = oldSamps[whichPart];
        tmpMods[part] = oldMods[whichPart];
    }
    
    //overwrite olds with news
    oldMods = std::move(tmpMods);
    oldSamps = std::move(tmpSamps);
    std::fill (oldLogWts.begin(), oldLogWts.end(), 0.0);
}


void MultinomResamp::resampHRBPF(std::vector<FSHMM> &oldMods, std::vector<Vec> &oldSamps, std::vector<double> &oldLogUnNormWts)
{
    // get dimensions
    unsigned numParticles = oldLogUnNormWts.size();
    
    // Create the distribution with those weights
    std::vector<double> w;
    w.resize(numParticles);
    double m = *std::max_element(oldLogUnNormWts.begin(), oldLogUnNormWts.end());
    std::transform(oldLogUnNormWts.begin(),
                    oldLogUnNormWts.end(),
                    w.begin(),
                    [&m](double &d)->double{return std::exp(d-m);}
                    );
    std::discrete_distribution<> idxSampler(w.begin(), w.end());
    
    // create temporary vectors for samps and mods
    std::vector<Vec>   tmpSamps(numParticles);
    std::vector<FSHMM> tmpMods(numParticles);
    
    // sample from the original parts and store in temporary
    int whichPart;
    for(size_t part = 0; part < numParticles; ++part)
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


void MultinomResamp::kGen(const std::vector<double> &logFirstStageWeights, std::vector<unsigned int> &ks)
{
    // these log weights may be very negative. If that's the case, exponentiating them may cause underflow
    // so we use the "log-exp-sum" trick
    // actually not quite...we just shift the log-weights because after they're exponentiated
    // they have the same normalized probabilities
    
   // Create the distribution with exponentiated log-weights
   // subtract the max first to prevent underflow
   // normalization is taken care of by std::discrete_distribution
    std::vector<double> w;
    unsigned int dim = logFirstStageWeights.size();
    w.resize(dim);
    double m = *std::max_element(logFirstStageWeights.begin(), logFirstStageWeights.end());
    std::transform(logFirstStageWeights.begin(), 
                   logFirstStageWeights.end(), 
                   w.begin(), 
                   [&m](const double& d) -> double { return std::exp(d-m); } );
    std::discrete_distribution<> kGenerator(w.begin(), w.end());
    
    // sample ks
    for(size_t i = 0; i < dim; ++i){
        ks[i] = kGenerator(m_gen);
    }
}
