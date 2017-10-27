#ifndef PMCMC_MSVOL_PITT_SHEPHARD_H
#define PMCMC_MSVOL_PITT_SHEPHARD_H

#include <vector>

#include "densities.h"
#include "msvol_sisr.h" // for model

namespace PittShep{

double twiceFisher(const double &phi);


double invTwiceFisher(const double &psi);


void qSample(const std::vector<double> &oldParams, 
                   std::vector<double> &newParams,
                   densities::EigenMultivariateNormalSampler &s);
                   
                   
double logQEvaluate(const std::vector<double> &oldParams, const std::vector<double> &newParams);

double logPriorEvaluate(const std::vector<double> &theta);


void flattenParams(std::vector<double> &flatOnes, const Mat &beta, const Vec &phis, const Vec &mus, const Vec &sigmas);


void unFlattenParams(Mat &beta, Vec &phis, Vec &mus, Vec &sigmas, const std::vector<double> &flatOnes);


void commenceSampling(std::string outFile, int numParts, int numMCMCIters);


} // PittShep

#endif // PMCMC_MSVOL_PITT_SHEPHARD_H
