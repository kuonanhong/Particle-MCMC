#include <string> 
#include <iostream>
#include <fstream> //ifstream, ofstream
#include <cmath> //log, pow, isnan
#include <Eigen/Dense>

#include "pmcmc_msvol_pitt_shephard.h"
#include "msvol_sisr.h"

// betas(9), phis(10),  mus(9), sigma(10)
// mus and phis seem to be messing up f(), and this makes the likelihood evaluate to 0.0
// right now these give us 0.07807808 mixing percent
static std::vector<double> qSigmaVec = {.06, .06, .06, .06, .06, .06, .06, .06, .06, // used to be .05
                                      .05, .05, .05, .05, .05, .05, .05, .05, .05, .05,
                                      .05, .05, .05, .05, .05, .05, .05, .05, .05,
                                      .1, .1, .1, .1, .1, .1, .1, .1, .1, .1};
                                      
double PittShep::twiceFisher(const double &phi)
{
    if ( (phi <= -1.0) || (phi >= 1.0) )
        throw std::invalid_argument( "received negative value" );
    else
        return log(1.0 + phi) - log(1.0 - phi);
}

double PittShep::invTwiceFisher(const double &psi)
{
    double ans = (1.0 - exp(psi)) / ( -1.0 - exp(psi) );
    if ( (ans <= -1.0) || (ans >= 1.0) )
        std::cerr << "ERROR\n";
    return ans;
}

void PittShep::qSample(const std::vector<double>                       &oldParams, 
                   std::vector<double>                       &newParams, 
                   densities::EigenMultivariateNormalSampler &s)
{
    int c = 0; 

    //betas
    for(int i = 0; i < 9; ++i){
        if ( i == 0)
            newParams[c] = oldParams[c];
        else
            newParams[c] = oldParams[c] + qSigmaVec[c] * s.sample()[0];
        c++;
    }
    
    // phis 
    for(int i = 0; i < 10; ++i){
        newParams[c] = invTwiceFisher( twiceFisher(oldParams[c]) + qSigmaVec[c] * s.sample()[0] );
        c++;
    }
    
    // mus
    for( int i = 0; i < 9; ++i){
        newParams[c] = oldParams[c] + qSigmaVec[c] * s.sample()[0];
        c++;
    }
    
    // sigmas
    for(int i = 0; i < 10; ++i){
        newParams[c] = exp( log(oldParams[c]) + qSigmaVec[c] * s.sample()[0] );
        c++;
    }
}


double PittShep::logQEvaluate(const std::vector<double> &oldParams, 
                 const std::vector<double> &newParams)
{
    double ans(0.0);
    int c = 1; // skipping the first beta
    
    // betas
    Vec betaX (9-1);
    Vec betaMean(9-1);
    Mat betaCov = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Identity(9-1, 9-1);
    for (int i = 0; i < 9-1; ++i){
        betaX(i) = newParams[c];
        betaMean(i) = oldParams[c];
        betaCov(i,i) = pow(qSigmaVec[c], 2);
        c++;
    }
    ans += log( densities::evalMultivNorm(betaX, betaMean, betaCov) );
    
    if (std::isnan(ans)){
        std::cerr << "betas generated nans in qEval\n";
        for(int i = 0; i < 9; ++i){
            std::cerr << betaX(i);
            std::cerr << betaMean(i);
        }
        std::cerr << "\n";
    }
    
    // phis
    Vec phisX (10);
    Vec phisMean(10);
    Mat phisCov = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Identity(10, 10);
    for (int i = 0; i < 10; ++i){
        phisX(i) = twiceFisher( newParams[c] );
        phisMean(i) = twiceFisher( oldParams[c] );
        phisCov(i,i) = pow(qSigmaVec[c],2);
        c++;
    }
    ans += log( densities::evalMultivNorm(phisX, phisMean, phisCov) );    
    
    if (std::isnan(ans)){
        std::cerr << "phis generated nans in qEval\n";
        for(int i = 0; i < 9; ++i){
            std::cerr << phisX(i) <<",";
            std::cerr << phisMean(i)<<",";
        }
        std::cerr << "\n";
    }
    
    // mus
    Vec musX (9);
    Vec musMean(9);
    Mat muCov = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Identity(9,9);
    for (int i = 0; i < 9; ++i){
        musX(i) = newParams[c];
        musMean(i) = oldParams[c];
        muCov(i,i) = pow(qSigmaVec[c],2);
        c++;        
    }
    ans += log( densities::evalMultivNorm(musX, musMean, muCov) ); 
    
    if (std::isnan(ans)){
        std::cerr << "mus generated nans in qEval\n";
        for(int i = 0; i < 9; ++i){
            std::cerr << musX(i);
            std::cerr << musMean(i);
        }
        std::cerr << "\n";
    }    
    
    // sigmas
    Vec sigmasX(10);
    Vec sigmasMean(10);
    Mat sigmaCov = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Identity(10,10);
    for (int i = 0; i < 10; ++i){
        sigmasX(i) = log( newParams[c] );
        sigmasMean(i) = log( oldParams[c] );
        sigmaCov(i,i) = pow(qSigmaVec[c], 2);
        c++;        
    }
    ans += log( densities::evalMultivNorm(sigmasX, sigmasMean, sigmaCov) );  
    if (std::isnan(ans)){
        std::cerr << "sigmas generated nans in qEval\n";
        for(int i = 0; i < 9; ++i){
            std::cerr << sigmasX(i);
            std::cerr << sigmasMean(i);
        }
        std::cerr << "\n";
    }    
    
    return ans;
}

double PittShep::logPriorEvaluate(const std::vector<double> &theta)
{
    int c = 1; // skipping first beta
    double returnThis(0.0);
    
    // betas (same prior as mus in time-varying-covs paper)
    Vec betaVec(9-1);
    Vec betaMean(9-1);
    Mat betaCovMat = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Identity(9-1,9-1);
    for (int row = 0; row < 9-1; ++row){
        betaVec(row) = theta[c];
        betaMean(row) = 1.0;
        betaCovMat(row,row) = 25.0;
        c++;
    }
    returnThis += log( densities::evalMultivNorm(betaVec, betaMean, betaCovMat) );
    
    // phis
    for (int i = 0; i < 9 + 1; ++i){
        returnThis += log( densities::evalUnivBeta( (theta[c]+1.0)/2.0, 18.0, 1.0) );
        c++;
    }
    
    // mus
    Vec muArg(9);
    Vec muMean(9);
    for (int i = 0; i < 9; ++i) { 
        muArg(i) = theta[c];
        muMean(i) = 1.0;
        c++;
    }
    returnThis += log( densities::evalMultivNorm(muArg, muMean, betaCovMat) ); // same mean and covariance in time-varying covs paper
    
    // sigma squareds (stored as sigmas)
    for (int i = 0; i < 9 + 1; ++i){
        //returnThis *= densities::evalUnivInvGamma( pow(theta[c],2) , 5.0, .05);  NOOO
        returnThis += log( densities::evalUnivInvGamma( pow(theta[c],2.0), 5.0, .05) );
        c++;
    }

    return returnThis;
}


void PittShep::flattenParams(std::vector<double> &flatOnes, const Mat &beta, const Vec &phis, const Vec &mus, const Vec &sigmas)
{
    int c = 0;
    
    for(int row = 0; row < beta.rows(); ++row){
        for(int col = 0; col < beta.cols(); ++col){
            flatOnes[c] = beta(row,col);
            c++;
        }
    }
    
    for(int row = 0; row < phis.rows(); ++row){
        flatOnes[c] = phis(row);
        c++;
    }
    
    for(int row = 0; row < mus.rows(); ++row){
        flatOnes[c] = mus(row);
        c++;
    }
    
    for(int row = 0; row < sigmas.rows(); ++row){
        flatOnes[c] = sigmas(row);
        c++;
    }
}


void PittShep::unFlattenParams(Mat &beta, Vec &phis, Vec &mus, Vec &sigmas, const std::vector<double> &flatOnes)
{   
    int c = 0;
    
    for(int row = 0; row < beta.rows(); ++row){
        for(int col = 0; col < beta.cols(); ++col){
            beta(row,col) = flatOnes[c];
            c++;
        }
    }
    
    for(int row = 0; row < phis.rows(); ++row){
        phis(row) = flatOnes[c];
        c++;
    }
    
    for(int row = 0; row < mus.rows(); ++row){
        mus(row) = flatOnes[c];
        c++;
    }
    
    for(int row = 0; row < sigmas.rows(); ++row){
        sigmas(row) = flatOnes[c];
        c++;
    } 
}

void PittShep::commenceSampling(std::string outFile, int numParts, int numMCMCIters)
{
    ////////////////////////
    // defining variables //
    ////////////////////////
    
    SISRResampStyle rt = SISRResampStyle::everytime_multinomial;
    //double percEss(.8);
    //int pathLength = 2000; // THIS IS ONLY THE LENGTH OF DATA IF YOU WANT TO STORE PATHS!
    int dimObs = 9;
    int numFactors = 1;
    int numParams = numFactors*dimObs + (numFactors + dimObs) + dimObs + (numFactors + dimObs);
    std::string dataFile = "/home/taylor/ssm/data/some_csvs/msvol_y_data.csv";
    std::ofstream outFileStream;
    outFileStream.open(outFile);
    
    // select initial parameters
    Mat init_beta(dimObs, numFactors); 
    Vec init_phis(numFactors + dimObs);
    Vec init_mus(dimObs);
    Vec init_sigmas(numFactors + dimObs);
    for (int row = 0; row < dimObs; ++row){
        if (row == 0){
            init_beta(row,0) = 1.0;
            init_mus(row) = -1.0;
        } else {
            init_beta(row,0) = 1.5;
            init_mus(row) = -1.0;            
        }
    }
    for(int row = 0; row < dimObs + numFactors; ++row){
        init_phis(row) = 0.0;
        init_sigmas(row) = .5;
    }
    std::vector<std::vector<double> > thetaChain(numMCMCIters, std::vector<double>(numParams));
    
    /////////////////////////////////////////////
    // read in data to vector<vector<double> > //
    /////////////////////////////////////////////
    std::string line;
    std::ifstream inFile(dataFile);
    std::string oneNumberOnOneLine;
    std::vector<Vec> data(2000, Vec(9));
    int time(0);
    int elemNum(0);
    while ( std::getline(inFile, line) ){ // every row of data
    
        elemNum = 0;
        try{
            std::istringstream buff(line);
            while(std::getline(buff, oneNumberOnOneLine, ',')){
                data[time](elemNum) = std::stod(oneNumberOnOneLine);
                elemNum++;
            }
        } catch (const std::invalid_argument& ia){
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }          
        time++;
    }

    //////////////
    // Iterate  //
    ////////////// 
    double oldLogLike(0.0);
    double oldLogPrior(0.0);
    double newLogLike(0.0);
    double newLogPrior(0.0);
    Vec zero(1);
    zero(0) = 0.0;
    Mat one(1, 1);
    one(0, 0) = 1.0;
    densities::EigenMultivariateNormalSampler rnorm(zero, one); //RNG
    densities::UniformSampler                 runif; //RNG
    
    for(int iter = 0; iter < numMCMCIters; ++iter) // every iteration
    {        
        std::cout << "***Iter: " << iter+1 << " out of " << numMCMCIters << " with " <<numParts<< " particles\n";
        
        if (iter == 0) { // first iteration no acceptance probability
    
            // write out the parameters (no accept-reject on first iter)
            flattenParams(thetaChain[0], init_beta, init_phis, init_mus, init_sigmas);
            for(unsigned int i = 0; i < thetaChain[0].size(); ++i){
                if( i == 0){
                    outFileStream << thetaChain[iter][i];
                } else {
                    outFileStream << "," << thetaChain[iter][i];                                        
                }
            }
            outFileStream << "\n";
            
            // instantiate model:
            MSVolSISR mod(numParts, init_beta, init_phis, init_mus, init_sigmas, rt);/*, pathLength*/ //, percEss);
                
            // iterate over rows
            for(int row = 0; row < 2000; ++row){
                mod.filterOrSmooth(data[row]);   
                std::cout << mod.getLogCondLike() << "\n";
                oldLogLike += mod.getLogCondLike();
            }

            // store prior for next round
            oldLogPrior = logPriorEvaluate(thetaChain[0]);
            
            
        } else { // not the first iteration
        
            // propose a new Theta 
            std::vector<double> proposedTheta(numParams);
            qSample(thetaChain[iter-1], proposedTheta, rnorm);
        
//            std::cout << "previous thetas were \n";
//            for(int i = 0; i < 38; ++i){
//                std::cout << thetaChain[iter-1][i] << ", ";
//            }
//            std::cout << "\n";
//            std::cout << "proposed thetas are \n";
//            for(int i = 0; i < 38; ++i){
//                std::cout << proposedTheta[i] << ", ";
//            }
//            std::cout << "\n";        
        
            // store newPrior evaluation
            newLogPrior = logPriorEvaluate(proposedTheta);
    
            // store transition kernel evaluation (DONT NEED WE'RE USING SYMMETRIC Q)
            //double logQOldToNew = logQEvaluate(thetaChain[iter-1], proposedTheta);
            //double logQNewToOld = logQEvaluate(proposedTheta, thetaChain[iter-1]);

            // instantiate a new model to approx the likelihood
            Mat propBeta(dimObs, numFactors);
            Vec propPhis(dimObs + numFactors);
            Vec propMus(dimObs);
            Vec propSigmas(dimObs + numFactors);
            unFlattenParams(propBeta, propPhis, propMus, propSigmas, proposedTheta);
            MSVolSISR mod(numParts, propBeta, propPhis, propMus, propSigmas, rt);//, /*pathLength,*/ percEss);

            // iterate through data and store the log likelihood
            newLogLike = 0.0;
            for (int row = 0; row < 2000; ++row){
                mod.filterOrSmooth(data[row]);
                newLogLike += mod.getLogCondLike();       
            }
        
            // tmp
            std::cout << "oldLoglike: " << oldLogLike << "\n";
            std::cout << "newLoglike: " << newLogLike << "\n";
            //std::cout << "old params: " << thetaChain[iter-1][0] << thetaChain[iter-1][1] << "\n";
            //std::cout << "proposed params: " << proposedTheta[0] << proposedTheta[1] << "\n";
        
            double logAR = newLogPrior + newLogLike - oldLogPrior - oldLogLike; // symmetric Q
            //double logAR = newLogPrior + logQNewToOld + newLogLike - oldLogPrior - logQOldToNew - oldLogLike;
            //double logAR = newLogLike - oldLogLike; //no priors
            double draw = runif.sample();
            std::cout << "log draw: " << log( draw ) << "\n";
            std::cout << "logAR: " << logAR<<"\n";
            
            
            if ( std::isinf(-logAR) ){ // 0 acceptance rate
                thetaChain[iter] = thetaChain[iter-1];
                // oldPrior stays the same 
                // oldLogLike stays the same
            } else if ( logAR >= 0.0 ) { // 100 percent accept 
                thetaChain[iter] = proposedTheta;
                //oldLogPrior = newLogPrior;
                oldLogLike = newLogLike;
                std::cout << "accepted 100percent\n";
            } else if ( log(draw) <= logAR  ) { // probabilistic accept
                thetaChain[iter] = proposedTheta;
                //oldLogPrior = newLogPrior;
                oldLogLike = newLogLike;
                std::cout << "accepted probabilistically\n";
            } else if ( log(draw) > logAR ) { // do not accept
                thetaChain[iter] = thetaChain[iter-1];
                // oldPrior stays the same 
                // oldLogLike stays the same                
                std::cout << "rejected probabilistically\n";
            }else if (std::isnan(logAR) ){ // log-likelihood 
                std::cerr << "there was a NaN. Not accepting proposal.\n";
                std::cerr << "newLogLike: " << newLogLike << "\n";
                std::cerr << "oldLogLikeL " << oldLogLike << "\n";
                //std::cerr << "newLogPrior: " << newLogPrior << "\n";
                //std::cerr << "oldLogPrior: " << oldLogPrior << "\n";
                //std::cerr << "logQNewToOld: "<< logQNewToOld << "\n";
                //std::cerr << "logQOldToNew: " << logQOldToNew << "\n";
                thetaChain[iter] = thetaChain[iter-1];
                // oldPrior stays the same 
                // oldLogLike stays the same                         
            } else {
                std::cerr << "you coded your MCMC incorrectly\n";
                std::cerr << "logDraw is: " << log(draw) << "\n";
                std::cerr << "logAR is: " << logAR << "\n";
                std::cerr << "old parameters were: ";
                for(int i = 0; i < proposedTheta.size(); ++i)
                    std::cerr << thetaChain[iter-1][i];
                std::cerr << "\n";
                std::cerr << "proposed parameters were: ";
                for(int i = 0; i < proposedTheta.size(); ++i)
                    std::cerr << proposedTheta[i];
                std::cerr << "\n";
                //std::cerr << "newPrior: " << newLogPrior << "\n";
                //std::cerr << "oldPrior: " << oldLogPrior << "\n";
                std::cerr << "newLogLike: " << newLogLike << "\n";
                std::cerr << "oldLogLike: " << oldLogLike << "\n";
                //std::cerr << "qNewToOld: " << logQNewToOld << "\n";
                //std::cerr << "qOldToNew: " << logQOldToNew << "\n";
                
                std::cerr << "stopping...";
                break;
            }
            
            // write out thetas
            for(unsigned int i = 0; i < thetaChain[0].size(); ++i){
                if( i == 0){
                    outFileStream << thetaChain[iter][i];
                } else {
                    outFileStream << "," << thetaChain[iter][i];                                        
                }
            }
            outFileStream << "\n";

        }
    }
    
    // stop writing thetas
    outFileStream.close();
}










