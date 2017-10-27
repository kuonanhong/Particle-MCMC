#ifndef PMCMC_TRIAL_H
#define PMCMC_TRIAL_H

#include <vector>
#include "densities.h"

namespace pmcmc_trial{

    
void qSample(const std::vector<double> &oldParams, 
                   std::vector<double> &newParams, 
                   densities::EigenMultivariateNormalSampler &s);


double logQEvaluate(const std::vector<double> &oldParams, const std::vector<double> &newParams);


double getLogLike(const std::string& fp, unsigned int np, const double& osd, const double&ssd, const double&a);


void testpmcmc(unsigned int np, unsigned int numIters, const std::string& in, const std::string& out, bool allMessages = false);


} //namespace pmcmc_trial


#endif //PMCMC_TEST_H