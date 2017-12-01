#include <iostream>
#include <string> // for string
#include <fstream> // for fstream

#include "svol_filter.h" 
#include "svol_apf_filter.h"

void sVolTest()
{
    
    // instantiate model
    int numParts(10000);
    SVolFilter myMod(numParts, .5, .91, 1.0, 10); // true parameters from simulated data
    //SVolAPFFilter myMod2(numParts);
        
    // stream in data
    std::string filePath("/home/taylor/ssm/data/some_csvs/svol_y_data.csv");
    std::ifstream infile(filePath);
    std::string line;
    unsigned int time = 0;
    
    double like1(0.0);
    double like2(0.0);
    double like3(0.0);
    
    while(std::getline(infile, line))
    {
        
        if(time == 10)
            break;
            
        
        // construct yt
        Vec yt(1);
        try{
            yt(0) = std::stod(line);           
        }catch(const std::invalid_argument &ia){
            std::cerr << "Prolly first line...Invalid Argument: " << ia.what() << "\n";
            continue;
        }
        
        // try getting the filter mean
        std::vector<std::function<const Mat(const Vec&)> > fs;
        auto meanLambda = [](const Vec& x) {return Eigen::Map<const Mat>(x.data(), x.size(), 1);};
        fs.push_back( meanLambda );
        
        
        //myMod.filterOrSmooth(yt);
        myMod.filterOrSmooth(yt, fs);
        //myMod2.filterOrSmooth(yt, fs);
        

        
        //std::cout << "time: " << time << "\n";
//        std::cout //<< "cond likes: " 
//                  << myMod.getLogCondLike() << ", "
//                  << std::log(myMod2.getCondLike()) << "\n";

        
        
        //std::cout << myMod.getExpectations()[0] << ","
        //          << myMod2.getExpectations()[0] << "\n";    
                  
        //std::cout << "ESS: " << myMod.getESS() << "\n";
        //std::cout << "mean: " << myMod.getFilterMean()<< "\n";
        
        
        time += 1;
    }    
}
