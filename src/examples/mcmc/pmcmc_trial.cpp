#include <string>
#include <vector>
#include <iostream>
#include <fstream> //ifstream, ofstream
#include <stdexcept>

#include "pmcmc_trial.h"
#include "noisy_ar1_filter.h"
#include "densities.h"


// MH random walk tuning
static double osd_sigma (.05);
static double ssd_sigma (.05);
static double a_sigma (.05);
static std::vector<double> qSigmaVec = {osd_sigma, ssd_sigma, a_sigma};


void pmcmc_trial::qSample(const std::vector<double> &oldParams, 
                                std::vector<double> &newParams, 
                                densities::EigenMultivariateNormalSampler &s)
{
    //order: osd, ssd, a
    newParams[0] = oldParams[0] + qSigmaVec[0] * s.sample()[0];
    newParams[1] = oldParams[1] + qSigmaVec[1] * s.sample()[0];
    newParams[2] = oldParams[2] + qSigmaVec[2] * s.sample()[0];
}


double pmcmc_trial::logQEvaluate(const std::vector<double> &oldParams, const std::vector<double> &newParams)
{
    int dim = qSigmaVec.size();
    
    // make x
    Vec x (dim);
    x << newParams[0], newParams[1], newParams[2];
    
    // make mean
    Vec meanVec(dim);
    meanVec << oldParams[0], oldParams[1], oldParams[2];
    
    // make covariance
    Mat covMat = Mat::Identity(dim, dim);
    for (unsigned int i = 0; i < dim; ++i){
        covMat(i,i) = qSigmaVec[i] * qSigmaVec[i];
    }
    
    // return evaluation
    return densities::evalMultivNorm(x, meanVec, covMat, true);
}


double pmcmc_trial::getLogLike(const std::string& fp, unsigned int np, const double& osd, const double&ssd, const double&a){
        
    
    // make model
    NAr1Filter mod(np, osd, ssd, a, SISRResampStyle::everytime_multinomial, 0);
    double logLike (0.0);

    // iterate through data
    std::ifstream infile(fp);
    std::string line;
    unsigned int time = 0;
    while(std::getline(infile, line)) {

        // construct yt
        Vec yt(1);
        try {
            yt(0) = std::stod(line);
        } catch(const std::invalid_argument& ia) {
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }

        // update and get log of conditional likelihood
        mod.filterOrSmooth(yt);
        logLike += mod.getLogCondLike();

        // increment time
        time += 1;
    }
    
    return logLike;
}


void pmcmc_trial::testpmcmc(unsigned int np, unsigned int numIters, const std::string& in, const std::string& out, bool allMessages){
    
    //    std::string filePath("/home/taylor/ssm/data/some_csvs/noisy_ar1_data.csv");
    //    std::string outPath("/home/taylor/ssm/examples/pmcmc/pmmh_nar1_output.csv")

    // real parameters can be found in /home/taylor/ssm/data/simulate_scripts/sim_noisy_ar1_data.py 
    // ^ assuming I saved the last version I ran 

    // initial parameter selections
    double startAlpha(.31); // real is .31
    double startStateSD(1.0); // real is 1.0
    double startObsSd(1.5); // real is 1.5
    
    // initialize sampling stuff
    Vec zero(1);
    zero(0) = 0.0;
    Mat one(1, 1);
    one(0, 0) = 1.0;
    densities::EigenMultivariateNormalSampler rnorm(zero, one); //RNG
    densities::UniformSampler                 runif;    

//    // get length of data
//    std::ifstream infile(in);
//    std::string line;
//    unsigned int time = 0;
//    while(std::getline(infile, line)) { time += 1; }
//    int pl(time); // length of data (probably should have this be automatic)
//    std::cout << "this should be 400: " << time << "\n";

    // for writing out the chain
    std::ofstream outFile;
    outFile.open(out);

    // things we store 
    std::vector<std::vector<double> > thetaChain(numIters);
    std::vector<double> logLikes(numIters);

    // get first likelihood
    logLikes[0] = getLogLike(in, np, startObsSd, startStateSD, startAlpha);
    thetaChain[0].push_back(startObsSd);
    thetaChain[0].push_back(startStateSD);
    thetaChain[0].push_back(startAlpha);
    outFile << thetaChain[0][0] << "," << thetaChain[0][1]<< "," << thetaChain[0][2] << "\n";
    
    if (allMessages)
        std::cout << "new loglike: " << logLikes[0] << "\n";

    // finally.. iterate
    for(int iter = 1; iter < numIters; ++iter)
    {
        
        // print stuff if you want
        std::cout << "iteration number: " << iter <<"\n";
                      
        // first sample a new proposal
        std::vector<double> thetaProposal(3);
        qSample(thetaChain[iter-1], thetaProposal, rnorm);
        
        // get likelihood based on this proposal (we are ignoring paths here)
        double newLogLike = getLogLike(in, np, thetaProposal[0], thetaProposal[1], thetaProposal[2]);
        if (allMessages)
            std::cout << "new loglike: " << newLogLike << "\n";
      
        // this is assuming uniform priors for ease
        double logAcceptRatio = newLogLike - logLikes[iter-1];
        logAcceptRatio += logQEvaluate(thetaProposal, thetaChain[iter-1]); 
        logAcceptRatio -= logQEvaluate(thetaChain[iter-1], thetaProposal);
        double acceptRatio = std::exp(logAcceptRatio);
        
        // construct acceptance ratio
        if ( acceptRatio >= 1.0){
            if (allMessages)
                std::cout << "100% accept\n";
            acceptRatio = 1.0;
        }else if ( acceptRatio <= 1.0 ) {
            // do nothing
        }else{
            std::cout << "something unexpected happened\n";
            std::cout << "acceptRatio turned out to be: " << acceptRatio << "\n";
            std::cout << "logAcceptRatio turned out to be: " << logAcceptRatio << "\n";
            break;
        }
        
        // probabilistically accept or reject
        double draw = runif.sample();
        if (draw < acceptRatio){  // accept
            thetaChain[iter].push_back(thetaProposal[0]);
            thetaChain[iter].push_back(thetaProposal[1]);
            thetaChain[iter].push_back(thetaProposal[2]);
            logLikes[iter] = newLogLike;
        }
        else
        {
            thetaChain[iter] = thetaChain[iter-1];
            logLikes[iter] = logLikes[iter-1];
        }
        
        // write out result 
        outFile << thetaChain[iter][0] << "," << thetaChain[iter][1]<< "," << thetaChain[iter][2] << "\n";

    }
    
    outFile.close(); // finish writing theta chain 
}