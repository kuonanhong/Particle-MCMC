#include "hmm_test.h"

#include <iostream> // for cout
#include <fstream> // for ifstream
#include <string> // for string
#include "simple_hmm.h"

// the loglikelihood was verified against the following r code on the same data set
//library(RcppHMM)
//y <- read.csv("/home/taylor/ssm/data/some_csvs/hmm_y_data.csv", header=F)
//N = c("Low", "High")
//A <- matrix(c(0.9, 0.1,
//              0.1, 0.9), ncol= length(N), byrow = TRUE)
//Mu <- matrix(c(0, 0), ncol = length(N))
//Sigma <- array(c(1, 4.5), dim = c(1,1,length(N)))
//Pi <- rep(1/length(N), length(N))
//
//newModel <- verifyModel(list( "Model"="GHMM", 
//                              "StateNames" = N,
//                              "A" = A, 
//                              "Mu" = Mu, 
//                              "Sigma" = Sigma, 
//                              "Pi" = Pi))
//yArray <- array(y[,1], dim = c(1,150,1) )
//loglikelihood(newModel, sequences = yArray)


void hmmTest()
{
    // open up the data connection
    std::string filePath("/home/taylor/ssm/data/some_csvs/hmm_y_data.csv");
    std::ifstream infile(filePath);
    
    // these are assumed to be the same as data generating process 
    // see ~/ssm/data/simulate_scripts/gen_hmm_data.py
    double lowVar(1.0);
    double highVar(4.5);
    Mat transMat(2,2);
    transMat(0,0) = transMat(1,1) = .9;
    transMat(0,1) = transMat(1,0) = .1;
    Vec initStateDistn(2);
    initStateDistn(0) = initStateDistn(1) = .5;
    
    // low var, high var, init_state, transMat
    SimpleHmm mod(lowVar, highVar, initStateDistn, transMat);

    // filter through the data
    std::string num;
    //double ll(0.0);
    while(std::getline(infile, num))
    {
        Vec yt(1);
        yt(0) = std::stod(num);
        mod.update(yt);
        ll += std::log( mod.getCondLike() );
        //std::cout << std::log(mod.getCondLike()) << "\n";
        //std::cout << yt << ", " << mod.getFilterVec()(1) << "\n";

    }
    //std::cout << ll <<"\n";
    
}