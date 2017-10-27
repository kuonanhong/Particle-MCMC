#include <iostream>
#include <math.h> //pow, sqrt
#include <fstream>

#include "lgssm.h"

void noisyAr1Test()
{
    // instantiate lgssm model and associated thigns
    Vec initMean(1);
    initMean(0) = 0.0;
    Mat initVar(1,1);
    initVar(0,0) = 1.0/(1- pow(.91,2)); // sigma^2 / (1 -alpha^2 ) 
    Mat stateTrans(1,1);
    stateTrans(0,0) = 1.;
    Mat cholStateVar(1,1);
    cholStateVar(0,0) = 1.0; // sqrt(1)
    Mat stateInptAffector(1,1);
    stateInptAffector(0,0) = 0.0;
    Vec inpt(1);
    inpt(0) = 0.0; // set to zero just in case
    Mat obsMat(1,1);
    obsMat(0,0) = 1.0;
    Mat obsInptAffector(1,1);
    obsInptAffector(0,0) = 0.0;
    Mat cholObsVar(1,1);
    cholObsVar(0,0) = 1.5; // gamma
    Lgssm myMod(initMean, initVar);

    // stream in data
    std::string filePath("/home/taylor/ssm/data/some_csvs/noisy_ar1_data.csv");
    std::ifstream infile(filePath);
    std::string line;
    unsigned int time = 0;
    while(std::getline(infile, line))
    {
        // construct yt
        Vec yt(1);
        try{
            yt(0) = std::stod(line);            
        }catch(const std::invalid_argument &ia){
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }
        
        myMod.update(yt, stateTrans, cholStateVar, stateInptAffector, inpt, obsMat, obsInptAffector, cholObsVar);
        std::cout << "time: " << time << std::endl;
        std::cout << "yt : " << yt << std::endl;
        std::cout << "filt mean: " << myMod.getFiltMean() << std::endl;
        std::cout << "filt var:  " << myMod.getFiltVar() << std::endl;
        std::cout << "cond likelihood" << myMod.getLogCondLike() << std::endl;
        
                
       time += 1;
    }
}