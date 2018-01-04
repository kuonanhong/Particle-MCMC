#include <iostream>
#include <string> // for string
#include <fstream> // for fstream

#include "svol_apf_filter.h"

void sVolAPFTest()
{
    
    // instantiate model
    int np (500);
    SVolAPFFilter myMod(np);
        
    // stream in data
    std::string filePath("/home/taylor/ssm/data/some_csvs/svol_y_data.csv");
    std::ifstream infile(filePath);
    std::string line;
    unsigned int time = 0;
    
    // make lambda functions that help store filter mean approx.s
    auto idtyLambda = [](const Vec& x){ return Eigen::Map<const Mat>(x.data(), x.size(), 1); };
    std::vector<std::function<const Mat(const Vec&)> > fs;
    fs.push_back(idtyLambda);
    
    while(std::getline(infile, line))
    {
        // construct yt
        Vec yt(1);
        try{
            yt(0) = std::stod(line);           
        }catch(const std::invalid_argument &ia){
            std::cerr << "Prolly first line...Invalid Argument: " << ia.what() << "\n";
            continue;
        }
        
        
        myMod.filter(yt, fs);
        
        
//        std::cout <<  time << ","
//                  <<  myMod.getLogCondLike() << std::endl;


        std::cout << myMod.getExpectations()[0] << "\n";
                
        time += 1;
    }    
}
