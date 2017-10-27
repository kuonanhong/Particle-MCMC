
//#include <densities.h>
#include <iostream>
#include <string> // for string
#include <fstream> // for fstream

#include "lgssm.h"

void kFilterTest()
{
    // instantiate model
    Vec initMean(1);
    initMean(0) = -.26;
    Mat initVar(1,1);
    initVar(0,0) = .01;
    Lgssm myMod(initMean, initVar);
    
    // construct matrices to update model
    Mat stateTrans(1,1);
    stateTrans(0,0) = 1.;
    Mat cholStateVar(1,1);
    cholStateVar(0,0) = .1;
    Mat stateInptAffector(1,1);
    stateInptAffector(0,0) = .1;
    Vec inpt(1);
    inpt(0) = 1.0;
    Mat obsMat(2,1);
    obsMat(0,0) = obsMat(1,0) = 1.0;
    Mat obsInptAffector(2,1);
    obsInptAffector(0,0) = obsInptAffector(1,0) = 0.0;
    Mat cholObsVar(2,2);
    cholObsVar(0,0) = cholObsVar(1,1) = .1;
    cholObsVar(0,1) = cholObsVar(1,0) = 0.0;
        
    // stream in data
    std::string filePath("/home/taylor/ssm/data/some_csvs/gtemp.csv");
    std::ifstream infile(filePath);
    std::string line;
    std::string first, second;
    unsigned int time = 0;
    while(std::getline(infile, line))
    {
        std::stringstream linestream(line);
        std::getline(linestream, first, ',');
        std::getline(linestream, second, ',');
        
        // construct yt
        Vec yt(2);
        try{
            yt(0) = std::stod(first);
            yt(1) = std::stod(second);
            
        }catch(const std::invalid_argument &ia){
            std::cerr << "Prolly first line...Invalid Argument: " << ia.what() << "\n";
            continue;
        }
        
        myMod.update(yt, stateTrans, cholStateVar, stateInptAffector, inpt, obsMat, obsInptAffector, cholObsVar);
        std::cout << "time: " << time << std::endl;
        std::cout << "filt mean: " << myMod.getFiltMean() << std::endl;
        std::cout << "filt var:  " << myMod.getFiltVar() << std::endl;
        std::cout << "cond likelihood" << myMod.getLogCondLike() << std::endl;
        
                
       time += 1;
    }
    
}
