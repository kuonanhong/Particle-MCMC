#include "pmmh_jac_apf.h"

#include <iostream>

#include "densities.h"
#include "transformations.h"
#include "jacquier_et_al_apf.h"



Pmmh_jac_apf::Pmmh_jac_apf(unsigned int numParts, 
                            std::vector<double> startTheta, 
                            unsigned int numMCMCIters, 
                            const std::string& dataFile, 
                            unsigned int numCols,
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
    m_qSigmaVec = { 999.99, // doesn't matter
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_betaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_fisherPhiVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_muVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logSigmaVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d), 
                    std::sqrt(m_logRStdDevVar*pow(2.4,2)/m_d)};
}


Pmmh_jac_apf::~Pmmh_jac_apf(){}
    
    
void Pmmh_jac_apf::flattenParams(std::vector<double> &flatOnes, const Mat &beta, const Vec &phis, const Vec &mus, const Vec &sigmas, const Vec &RstdDevs)
{
    
    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    int c = 0;
    
    // betas
    for(int row = 0; row < beta.rows(); ++row){
        for(int col = 0; col < beta.cols(); ++col){
            flatOnes[c] = beta(row,col);
            c++;
        }
    }
    
    // phis (only one)
    flatOnes[c] = phis(0);
    c++;
    
    // mus (only one)
    flatOnes[c] = mus(0);
    c++;
    
    // sigmas (only one)
    flatOnes[c] = sigmas(0);
    c++;
    
    // r std_devs (9)
    for(int row  = 0; row < RstdDevs.rows(); ++row){
        flatOnes[c] = RstdDevs(row);
        c++;
    }
}


void Pmmh_jac_apf::unFlattenParams(Mat &beta, Vec &phi, Vec &mu, Vec &sigma, Vec &RstdDevs, const std::vector<double> &flatOnes)
{
    
    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    int c = 0;
    
    // betas 
    for(int row = 0; row < beta.rows(); ++row){
        for(int col = 0; col < beta.cols(); ++col){
            beta(row,col) = flatOnes[c];
            c++;
        }
    }
    
    // phis (only one)
    phi(0) = flatOnes[c];
    c++;
    
    // mus (only one)
    mu(0) = flatOnes[c];
    c++;
        
    // sigmas (only one)
    sigma(0) = flatOnes[c];
    c++;
    
    // R std devs 
    for(int row = 0; row < RstdDevs.rows(); ++row){
        RstdDevs(row) = flatOnes[c];
        c++;
    }    
}

    
void Pmmh_jac_apf::qSample(const std::vector<double> &oldParams, std::vector<double> &newParams)
{

    // ordering is: betas(9-1), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    static densities::MVNSampler s;
    
    // start counter
    int c = 0; 

    //betas
    for(int i = 0; i < 9; ++i){
        if ( i == 0)
            newParams[c] = oldParams[c];  // keep the first fixed to handle identifiability
        else
            newParams[c] = oldParams[c] + m_qSigmaVec[c] * s.sample()[0];
        c++;
    }
    
    // only one phi to handle in this case
    newParams[c] = transformations::invTwiceFisher( transformations::twiceFisher(oldParams[c]) + m_qSigmaVec[c] * s.sample()[0] );
    c++;

    // only one mu to handle
    newParams[c] = oldParams[c] + m_qSigmaVec[c] * s.sample()[0];
    c++;

    // only one sigma
    newParams[c] = std::exp( std::log(oldParams[c]) + m_qSigmaVec[c] * s.sample()[0]);
    c++;

    // now r std. devs on diagonal
    for(int i = 0; i < 9; ++i){
        newParams[c] = std::exp( std::log(oldParams[c]) + m_qSigmaVec[c] * s.sample()[0]);
        c++;
    }    
}

                         
double Pmmh_jac_apf::logQEvaluate  (const std::vector<double> &oldParams, const std::vector<double> &newParams)
{
    
    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    double ans(0.0);
    int c (1); 
    
    // betas (do not include first beta in construction because it isn't random)
    Vec betaX (9-1);
    Vec betaMean(9-1);
    Mat betaCov = Mat::Identity(9-1, 9-1);
    for (int i = 0; i < 9-1; ++i){
        betaX(i) = newParams[c];
        betaMean(i) = oldParams[c];
        betaCov(i,i) = std::pow(m_qSigmaVec[c], 2);
        c++;
    }
    ans += densities::evalMultivNorm(betaX, betaMean, betaCov, true);
    
    // only one phi
    Vec phisX(1);
    Vec phisMean(1);
    Mat phisCov = Mat::Identity(1, 1);
    phisX(0) = transformations::twiceFisher( newParams[c] );
    phisMean(0) = transformations::twiceFisher( oldParams[c] );
    phisCov(0,0) = pow(m_qSigmaVec[c],2);
    c++;
    ans += densities::evalMultivNorm(phisX, phisMean, phisCov, true);    

    // only one mu
    Vec musX (1);
    Vec musMean(1);
    Mat muCov(1,1);
    musX(0) = newParams[c];
    musMean(0) = oldParams[c];
    muCov(0,0) = pow(m_qSigmaVec[c],2);
    c++;        
    ans += densities::evalMultivNorm(musX, musMean, muCov, true); 
    
    // one sigma
    Vec sigmasX(1);
    Vec sigmasMean(1);
    Mat sigmaCov(1,1);
    sigmasX(0) = std::log( newParams[c] );
    sigmasMean(0) = std::log( oldParams[c] );
    sigmaCov(0,0) = std::pow(m_qSigmaVec[c], 2);
    c++;        
    ans += densities::evalMultivNorm(sigmasX, sigmasMean, sigmaCov, true);  

    // observation covariance matrix
    Vec sigmasR(9);
    Vec sigmasRMean(9);
    Mat sigmasRCov = Mat::Identity(9,9);
    for(int i = 0; i < 9; ++i){
        sigmasR(i) = std::log( newParams[c]);
        sigmasRMean(i) = std::log( oldParams[c]);
        sigmasRCov(i,i) = std::pow(m_qSigmaVec[c], 2);
        c++;
    }
    ans += densities::evalMultivNorm(sigmasR, sigmasRMean, sigmasRCov, true);
    
    return ans;    
}


double Pmmh_jac_apf::logPriorEvaluate(const std::vector<double> &theta)
{
    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    int c = 1; // skipping first beta
    double returnThis(0.0);
    
    // betas(9-1) 
    Vec betaVec(9-1);
    Vec betaMean(9-1);
    Mat betaCovMat = Mat::Identity(9-1,9-1);
    for (int row = 0; row < 9-1; ++row){
        betaVec(row) = theta[c];
        betaMean(row) = 1.0;
        betaCovMat(row,row) = .125;
        c++;
    }
    returnThis += densities::evalMultivNorm(betaVec, betaMean, betaCovMat, true);
    
    // these parameters were set to make the mean and variacne of phi .95 and .0025 respectively
    returnThis += densities::evalUnivBeta( (theta[c]+1.0)/2.0, 5.114, 0.131, true) - std::log(2.0);
    c++;
    
    // one mu
    Vec muArg(1);
    Vec muMean(1);
    Mat muCov(1,1);
    muArg(0) = theta[c];
    muMean(0) = 0.0;
    muCov(0,0) = .001;
    c++;
    returnThis += densities::evalMultivNorm(muArg, muMean, muCov, true); 

    // one sigma (stored as sigmas)
    returnThis += densities::evalUnivInvGamma(std::pow(theta[c], 2.0), .5, .5, true); 
    c++;
    
    // r variances now (9 of them)
    // each r^2 is assumed to be independent and follow an inv gamma distn
    // the true value for the variances is .5 remember
    // infinite variance again
    for(int i = 0; i < 9; ++i){
        //returnThis += log( densities::evalUnivInvGamma( pow(theta[c],2.0), 2.0, .5) );
        returnThis += densities::evalUnivInvGamma( std::pow(theta[c],2.0), 1.0, 1.0, true);
        c++;
    }

    return returnThis;    
}


double Pmmh_jac_apf::logLikeEvaluate(const std::vector<double>& theta, const std::vector<Vec>& data, std::atomic_bool& cancelled)
{
    // we're goign to return this
    double logLike(0.0);
    
    // flatten parameters
    Mat beta(9, 1); 
    Vec phi(1);
    Vec mu(1);
    Vec sigma(1);
    Vec Rsigmas(9);
    unFlattenParams(beta, phi, mu, sigma, Rsigmas, theta);

    // instantiate model:
    JacEtAlAPF mod(m_numParts, beta, Rsigmas, mu, phi, sigma, APFResampStyle::everytime_multinomial); 

    // iterate over time
    // TODO: make sure that data.size() > 0
    unsigned int row (0);
    while(row < data.size() && !cancelled){
        mod.filterOrSmooth(data[row]);  
        logLike += mod.getLogCondLike();
        row++;
    }
    
    return logLike;    
}

