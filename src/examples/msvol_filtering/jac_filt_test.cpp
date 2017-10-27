//#include <string> //string
#include <iostream>
#include <fstream> //ifstream, ofstream
#include "jacquier_et_al.h" // model class

void jac_filt_test()
{
    ////////////////////////////////////
    // algorithm and example variables /
    ////////////////////////////////////
    int num_parts = 1000;
    
    ////////////////////////////////////////////////////////////
    // instantiate model object with known parameters          /
    // make sure these parameters cohere with                  /
    // /home/taylor/ssm/data/simulate_scripts/gen_jacq_data.py /
    ////////////////////////////////////////////////////////////
    int num_factors = 1;
    int dim_obs = 9;
     

    Mat beta(dim_obs, num_factors);    
    Vec R(dim_obs); 
    Vec phis(num_factors);
    Vec mus(num_factors);
    Vec sigmas(num_factors);

    // this only works for this particular case of one factor
    phis(0) = .95;//.92;
    sigmas(0) =  0.299583;//.447; //sqrt .2
    mus(0) = .09;//.74;
    beta(0,0) = 1.1;
    beta(1,0) = 1.0; 
    beta(2,0) = .9;
    beta(3,0) = 1.1;
    beta(4,0) = 1.0;
    beta(5,0) = .9;
    beta(6,0) = 1.1;
    beta(7,0) = 1.0;
    beta(8,0) = .9;
    
    for(int i = 0; i < dim_obs; ++i){
        //beta(i,0) = 1.0;
        R(i) = 1.5;//.8; // the variances not the std devs
    }
    SISRResampStyle rt = SISRResampStyle::everytime_multinomial;
    int path_length = 0; // filtering only
    JacEtAl myMod(num_parts, beta, R, mus, phis, sigmas, rt, path_length);

    ////////////////////////////////////////
    // now run through the data and filter /
    ////////////////////////////////////////
    std::string filePath("/home/taylor/ssm/data/some_csvs/jacq_y_data.csv");
    std::ifstream inFile(filePath);
    std::string dataLine;
    std::string oneNumber;
    unsigned int time = 0;    
    while(std::getline(inFile, dataLine)) {

        // construct yt
        Vec yt(dim_obs);
        int i;
        try {
            std::istringstream buff(dataLine);
            i = 0;
            while(std::getline(buff, oneNumber, ',')){
                yt(i) = std::stod(oneNumber);
                i += 1;
            }
        } catch(const std::invalid_argument& ia) {
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }

        // update
        myMod.filterOrSmooth(yt);
        //std::cout << myMod.getFilterMean() << "\n";
        std::cout << myMod.getLogCondLike() << "\n";
        
        
        // increment time
        time += 1;
    }
}