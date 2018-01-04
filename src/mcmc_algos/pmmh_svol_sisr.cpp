#include "pmmh_svol_sisr.h"

#include <iostream>
#include <cmath>

#include "densities.h"
#include "transformations.h"
#include "svol_filter.h"


Pmmh_svol_sisr::Pmmh_svol_sisr(unsigned numParts, 
                                const std::vector<Vec> &startTheta,
                                const std::vector<Mat> &qSigma, 
                                unsigned numMCMCIters, 
                                const std::string &dataFile, 
                                unsigned numCols,
                                bool mc) : 
    m_numParts(numParts),
    m_d(3),
    m_qSigma(qSigma),
    Pmmh(startTheta, numMCMCIters, dataFile, numCols, mc)
{
}


Pmmh_svol_sisr::~Pmmh_svol_sisr(){}

    
void Pmmh_svol_sisr::qSample(const std::vector<Vec> &oldParams, std::vector<Vec> &newParams)
{

    // ordering is: beta, phi, sigma
    static densities::MVNSampler s;
    
    // make sure we can use push_back
    if(newParams.size() == 0){
        throw std::out_of_range("newParams in qSample() should be greater than length 0.");        
    }
    
    //beta
    s.setCovar(m_qSigma[0]);
    s.setMean(oldParams[0]);
    newParams[0] = s.sample();
   
    // fisher(phi)
    Vec newPhi(1);
    s.setCovar(m_qSigma[1]);
    s.setMean(oldParams[1].unaryExpr(&transformations::twiceFisher));
    newPhi = s.sample().unaryExpr(&transformations::invTwiceFisher);
    newParams[1] = newPhi;
    
    //sigma
    Vec newSigma(1);
    s.setCovar(m_qSigma[2]);
    s.setMean(oldParams[2].array().log().matrix());
    newSigma = s.sample().array().exp().matrix();
    newParams[2] = newSigma;

}

                         
double Pmmh_svol_sisr::logQEvaluate(const std::vector<Vec> &oldParams, const std::vector<Vec> &newParams)
{
    
    // ordering is: beta, phi, sigma
    double ans(0.0);
    
    // beta
    ans += densities::evalMultivNorm(newParams[0], oldParams[0], m_qSigma[0], true);
    
    // one phi
    ans += densities::evalMultivNorm(newParams[1].unaryExpr(&transformations::twiceFisher), 
                                     oldParams[1].unaryExpr(&transformations::twiceFisher), 
                                     m_qSigma[1], true);

    // sigma
    ans += densities::evalMultivNorm(newParams[2].array().log().matrix(),
                                     oldParams[2].array().log().matrix(), 
                                     m_qSigma[2], true);
    
    return ans;    
}


double Pmmh_svol_sisr::logPriorEvaluate(const std::vector<Vec> &theta)
{
    // real ones on sim data are alpha = .91, sigma=1.0, beta = .5
    // parameters are hard-coded
    // ordering is: beta, phi, sigma
    double returnThis(0.0);
    
    //beta
    Vec betaMean(1);
    betaMean(0) = 1.0;
    Mat betaCovMat = Mat::Identity(1,1) * 10;
    returnThis += densities::evalMultivNorm(theta[0],
                                            betaMean, 
                                            betaCovMat,
                                            true);
                                            
    // phi (mind the jacobian)
    returnThis += densities::evalUnivBeta( (theta[1](0)+1.0)/2.0, 1.0, 0.2, true) - std::log(2.0);

    //sigma (mode is 1 which lines up with truth of simulated data)
    returnThis += densities::evalUnivInvGamma(theta[2].array().square()(0), .001, .001, true);//.5, 1.5, true);  
    
    return returnThis;    
}


double Pmmh_svol_sisr::logLikeEvaluate(const std::vector<Vec>& theta, const std::vector<Vec>& data, std::atomic_bool& cancelled)
{
    // we're goign to return this
    double logLike(0.0);
    
    // instantiate model:
    double beta = theta[0](0);
    double phi = theta[1](0);
    double sigma = theta[2](0);
    SVolFilter mod(m_numParts, beta, phi, sigma);

    // iterate over time
    if(data.empty())
        throw std::length_error("Can't read data in logLikeEvaluate()\n");
        
    unsigned int row (0);
    while(row < data.size() && !cancelled){
        mod.filter(data[row]);  
        logLike += mod.getLogCondLike();
        row++;
    }

    return logLike;    
}

