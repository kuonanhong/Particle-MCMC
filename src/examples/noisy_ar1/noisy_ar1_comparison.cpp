#include <cmath> //pow, sqrt, log
#include <fstream> //ifstream, ofstream
#include <iostream>

#include "lgssm.h"
#include "noisy_ar1_filter.h"
#include "noisy_ar1_apf_filter.h"

void noisyAr1Comparison()
{
    
    double phi (.31);
    double stateSd (1.0);
    double obsSd (1.5);
    
    // instantiate LGSSM object and construct some matrices for using it
    Vec initMean(1);
    initMean(0) = 0.0;
    Mat initVar(1, 1);
    initVar(0, 0) = (stateSd * stateSd) / (1.0 - (phi * phi) ); // sigma^2 / (1 -alpha^2 )
    Mat stateTrans(1, 1);
    stateTrans(0, 0) = phi;
    Mat cholStateVar(1, 1);
    cholStateVar(0, 0) = stateSd; // sqrt(1)
    Mat stateInptAffector(1, 1);
    stateInptAffector(0, 0) = 0.0;
    Vec inpt(1);
    inpt(0) = 0.0; // set to zero just in case
    Mat obsMat(1, 1);
    obsMat(0, 0) = 1.0;
    Mat obsInptAffector(1, 1);
    obsInptAffector(0, 0) = 0.0;
    Mat cholObsVar(1, 1);
    cholObsVar(0, 0) = obsSd; // gamma
    Lgssm myLgssmMod(initMean, initVar);
    
    // instantiate bootstrap filter and apf for same model
    int numParticles = 10;
    SISRResampStyle rt1 = SISRResampStyle::everytime_multinomial;
    APFResampStyle rt2 = APFResampStyle::everytime_multinomial;    
    NAr1Filter mySISRMod(numParticles, obsSd, stateSd, phi, rt1, 0);
    NAr1APFFilter myAPFMod(numParticles, obsSd, stateSd, phi, rt2, 0);
    
    // stream in data
    std::string filePath("/home/taylor/ssm/data/some_csvs/noisy_ar1_data.csv");
    std::ifstream infile(filePath);
    std::string line;
    unsigned int time = 0;
    
    // lambda function that gets filter means
    std::vector<std::function<const Mat(const Vec&)> > fs;
    auto idtyLambda = [](const Vec& x) { return Eigen::Map<const Mat>(x.data(), x.size(), 1); };
    fs.push_back(idtyLambda);

    while(std::getline(infile, line)) {
        
        // construct yt
        Vec yt(1);
        try {
            yt(0) = std::stod(line);
        } catch(const std::invalid_argument& ia) {
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }

        // update all models (maybe problem below line)
        myLgssmMod.update(yt, stateTrans, cholStateVar, stateInptAffector, inpt, obsMat, obsInptAffector, cholObsVar);
        mySISRMod.filter(yt, fs);
        myAPFMod.filter(yt, fs);


        //std::cout << myLgssmMod.getLogCondLike() << ","
        //          << mySISRMod.getLogCondLike() << ", "
        //          << myAPFMod.getLogCondLike() << "\n";



        // filter mean stuff
        std::cout << time << ","
                  << yt   << ","
                  << myLgssmMod.getFiltMean() << ","
                  << mySISRMod.getExpectations()[0] << ", "
                  << myAPFMod.getExpectations()[0] << "\n";

        // uncomment for conditional likelihood suff
//        if (time == 1){
//        std::cout << myLgssmMod.getLogCondLike() << ","
//                  << mySISRMod.getLogCondLike() << ", "
//                  << myAPFMod.getLogCondLike() << "\n";

//            std::cout << myLgssmMod.getFiltMean() << ", "
//                      << mySISRMod.getFilterMean() << ", "
//                      << myAPFMod.getFilterMean() << "\n";
//
//            break;
//        }



        time += 1;
        
    
    }
    //std::cout << totalLogLikeKalman << ", " << totalLogLikeSISR << ", "<< totalLogLikeAPF << "\n";

    // print out entire paths
//    std::vector<std::vector<Vec> > derp1 = mySISRMod.getFullParts();
//    std::vector<std::vector<Vec> > derp2 = myAPFMod.getFullParts();
//
//    // write out
//    std::ofstream outFile;
//    outFile.open("/home/taylor/Desktop/example.csv");
//    for(int time = 0; time < 200; ++time) {
//        for(int part = 0; part < numParticles; ++part) {
//            outFile << derp1[time][part] << ",";
//        }
//        outFile << "\n";
//    }
//    outFile.close();
//
//    std::ofstream outFile2;
//    outFile2.open("/home/taylor/Desktop/example2.csv");
//    for(int time = 0; time < 200; ++time)
//    {
//       for(int part = 0; part < numParticles; ++part)
//       {
//           outFile2 << derp2[time][part] << ",";
//       }
//       outFile2 << "\n";
//    }
}