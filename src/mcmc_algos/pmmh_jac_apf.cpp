#include "pmmh_jac_apf.h"

#include <iostream>

#include "densities.h"
#include "transformations.h"
#include "transformations.h"
#include "jacquier_et_al_apf.h"



Pmmh_jac_apf::Pmmh_jac_apf(unsigned int numParts, 
                            std::vector<Vec> startTheta, 
                            unsigned numMCMCIters, 
                            const std::string& dataFile, 
                            unsigned numCols,
                            bool mc) : 
    m_numParts(numParts),
    m_d(21),//m_d(21),
    m_betaVar(.005),
    m_fisherPhiVar(.04),
    m_muVar(.0005),
    m_logSigmaVar(.037),
    m_logRStdDevVar(.006),
    Pmmh(startTheta, numMCMCIters, dataFile, numCols, mc)
{
    Mat betaVar = Mat::Identity(9,9) * m_betaVar* 2.4* 2.4 / m_d;
    Mat fisherPhiCov = Mat::Identity(1,1) * m_fisherPhiVar * 2.4 * 2.4 / m_d;
    Mat muCov = Mat::Identity(1,1) * m_muVar * 2.4 * 2.4 / m_d;
    Mat logSigmaCov = Mat::Identity(1,1) * m_logSigmaVar * 2.4 * 2.4 / m_d;
    Mat logRStdDevCov = Mat::Identity(9,9) * 2.4 * 2.4 / m_d;
    std::vector<Vec> m_qSigma;
    m_qSigma.push_back(betaVar);
    m_qSigma.push_back(fisherPhiCov);
    m_qSigma.push_back(muCov);
    m_qSigma.push_back(logSigmaCov);
    m_qSigma.push_back(logRStdDevCov);
}


Pmmh_jac_apf::~Pmmh_jac_apf(){}
    
    
void Pmmh_jac_apf::qSample(const std::vector<Vec> &oldParams, std::vector<Vec> &newParams)
{

    // ordering is: betas(9-1), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    static densities::MVNSampler s;
    
    // make sure we can use push_back
    if(newParams.size() == 0){
        throw std::out_of_range("newParams in qSample() should be length  greater than 0.");        
    }
    
    //betas
    Vec newBeta(9);
    Vec oldBeta = oldParams[0];
    for(int i = 0; i < 9; ++i){
        if ( i == 0)
            newBeta(i) = oldBeta(i);  // keep the first fixed to handle identifiability
        else
            newBeta(i) = oldBeta(i) + m_qSigma[0](i) * s.sample()[0];
    }
    newParams[0] = newBeta;
    
    // only one phi to handle in this case
    Vec newPhi(1);
    s.setCovar(m_qSigma[1]);
    s.setMean(oldParams[1].unaryExpr(&transformations::twiceFisher));
    newPhi = s.sample().unaryExpr(&transformations::twiceFisher);
    newParams[1] = newPhi;
    
    // only one mu to handle
    Vec newMu(1);
    s.setCovar(m_qSigma[2]);
    s.setMean(oldParams[2]);
    newMu = s.sample();
    newParams[2] = newMu;

    // only one sigma
    Vec newSigma(1);
    s.setCovar(m_qSigma[3]);
    s.setMean(oldParams[3].array().log().matrix());
    newSigma = s.sample().array().exp().matrix();
    newParams[3] = newSigma;

    // now r std. devs on diagonal
    Vec newRSds(9);
    s.setCovar(m_qSigma[4]);
    s.setMean(oldParams[4].array().log().matrix());
    newRSds = s.sample().array().exp().matrix();
    newParams[4] = newRSds;
}

                         
double Pmmh_jac_apf::logQEvaluate(const std::vector<Vec> &oldParams, const std::vector<Vec> &newParams)
{
    
    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    double ans(0.0);
    
    // betas (do not include first beta in construction because it isn't random)
    Vec betaX (9-1);
    Vec betaMean(9-1);
    Mat betaCov = Mat::Identity(9-1, 9-1);
    for (int i = 0; i < 9-1; ++i){
        betaX(i) = newParams[0](i+1);
        betaMean(i) = oldParams[0](i+1);
        betaCov(i,i) = m_qSigma[0](i+1);
    }
    ans += densities::evalMultivNorm(betaX, betaMean, betaCov, true);
    
    // only one phi
    ans += densities::evalMultivNorm(newParams[1].unaryExpr(&transformations::twiceFisher), 
                                     oldParams[1].unaryExpr(&transformations::twiceFisher), 
                                     m_qSigma[1], true);

    // only one mu
    ans += densities::evalMultivNorm(newParams[2], oldParams[2], m_qSigma[2], true);
    
    // one sigma
    ans += densities::evalMultivNorm(newParams[3].array().log().matrix(), 
                                     oldParams[3].array().log().matrix(), 
                                     m_qSigma[3], true);
    
    // observation covariance matrix
    ans += densities::evalMultivNorm(newParams[4].array().log().matrix(), 
                                     oldParams[4].array().log().matrix(), 
                                     m_qSigma[4], true);
    return ans;    
}


double Pmmh_jac_apf::logPriorEvaluate(const std::vector<Vec> &theta)
{
    // parameters are hard-coded
    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    double returnThis(0.0);
    
    // betas(9-1) 
    int c = 1; // skipping first beta
    Vec betaVec(9-1);
    Vec betaMean(9-1);
    Mat betaCovMat = Mat::Identity(9-1,9-1);
    for (int row = 0; row < 9-1; ++row){
        betaVec(row) = theta[0][c];
        betaMean(row) = 1.0;
        betaCovMat(row,row) = .125;
        c++;
    }
    returnThis += densities::evalMultivNorm(betaVec, betaMean, betaCovMat, true);
    
    // these parameters were set to make the mean and variacne of phi .95 and .0025 respectively
    returnThis += densities::evalUnivBeta( (theta[1](0)+1.0)/2.0, 5.114, 0.131, true) - std::log(2.0);
    
    // one mu
    Vec muMean(1);
    Mat muCov(1,1);
    muMean(0) = 0.0;
    muCov(0,0) = .001;
    returnThis += densities::evalMultivNorm(theta[2], muMean, muCov, true); 

    // one sigma (stored as sigmas)
    returnThis += densities::evalUnivInvGamma(theta[3].array().square()(0), 
                                                .5, .5, true); 
    
    // r variances now (9 of them)
    // each r^2 is assumed to be independent and follow an inv gamma distn
    // the true value for the variances is .5 remember
    // infinite variance again
    for(int i = 0; i < 9; ++i){
        //returnThis += log( densities::evalUnivInvGamma( pow(theta[c],2.0), 2.0, .5) );
        returnThis += densities::evalUnivInvGamma( std::pow(theta[4](i),2.0), 1.0, 1.0, true);
    }

    return returnThis;    
}


double Pmmh_jac_apf::logLikeEvaluate(const std::vector<Vec>& theta, const std::vector<Vec>& data, std::atomic_bool& cancelled)
{
    // we're goign to return this
    double logLike(0.0);
    
    // instantiate model:
    JacEtAlAPF mod(m_numParts, theta[0], theta[4], theta[2], theta[1], theta[3], APFResampStyle::everytime_multinomial); 

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

