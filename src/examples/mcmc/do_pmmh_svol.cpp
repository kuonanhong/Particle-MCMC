#include "do_pmmh_svol.h"

#include "pmmh_svol_sisr.h"

void do_pmmh_svol()
{
    // variables that we canchange
    unsigned np (200);
    unsigned numMCMC(1000);
    unsigned nc(1);
    bool mc(false);

    // start parameters
    Vec startBeta(1);
    startBeta(0) = .5;
    Vec startPhi(1);
    startPhi(0) = .9;
    Vec startSigma(1);
    startSigma(0) = .6;
    std::vector<Vec> startTheta;
    startTheta.push_back(startBeta);
    startTheta.push_back(startPhi);
    startTheta.push_back(startSigma);
    
    // start qSigma
    Mat qVar1(1,1);
    Mat qVar2(1,1);
    Mat qVar3(1,1);
    qVar1(0,0) = .01;
    qVar2(0,0) =  0.001;
    qVar3(0,0) = 0.0002;
    std::vector<Mat> qSigma;
    qSigma.push_back(qVar1);
    qSigma.push_back(qVar2);
    qSigma.push_back(qVar3);
    
    // make the object
    Pmmh_svol_sisr p(np, startTheta, qSigma, numMCMC, "/home/taylor/ssm/data/some_csvs/svol_y_data.csv", nc, mc);
    
    // begin the sampling
    p.commenceSampling("/home/taylor/Desktop/tempsamples1", "/home/taylor/Desktop/tempmessages1");
}
