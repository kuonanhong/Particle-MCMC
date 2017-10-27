//#include <string> //string

#include <iostream>
#include <fstream> //ifstream, ofstream
#include "jacquier_et_al_apf.h" // model class


void jac_filt_test_apf()
{
    // algorithm and example variables 
    int num_parts(3000);
    APFResampStyle rt = APFResampStyle::everytime_multinomial;
    int path_length = 0; // filtering only

    // model variables
    int num_factors = 1;
    int dim_obs = 9;

    // specify parameters
    Mat beta(dim_obs, num_factors);    
    Vec R_sigmas(dim_obs);
    Vec phis(num_factors);
    Vec mus(num_factors);
    Vec sigmas(num_factors);
    phis(0) = .95;//0.778770; 
    sigmas(0) = 0.299583;//0.0002810249;//0.413279; 
    mus(0) =  .09;//.0009;//-0.002749;
    beta(0,0) = 1.0;//1.000000;
    beta(1,0) = 1.0;//0.858130;
    beta(2,0) = .9;//1.329961;
    beta(3,0) = 1.1;//1.480949;
    beta(4,0) = 1.0;//0.689042;
    beta(5,0) = .9;//1.336924;
    beta(6,0) = 1.1;//1.182215;
    beta(7,0) = 1.0;//0.847413;
    beta(8,0) = .9;//1.295376;
    R_sigmas(0) = 1.224745;// 0.7071068;//2.500885;
    R_sigmas(1) = 1.224745;// 0.7071068;//1.592154;
    R_sigmas(2) = 1.224745;// 0.7071068;//1.154554;
    R_sigmas(3) = 1.224745;// 0.7071068;//1.477701;
    R_sigmas(4) = 1.224745;// 0.7071068;//0.887661;
    R_sigmas(5) = 1.224745;// 0.7071068;//0.785638;
    R_sigmas(6) = 1.224745;// 0.7071068;//0.892934;
    R_sigmas(7) = 1.224745;// 0.7071068;//0.966721;
    R_sigmas(8) = 1.224745;// 0.7071068;//1.407837;
    
    // instantiate model
    JacEtAlAPF myMod(num_parts, beta, R_sigmas, mus, phis, sigmas, rt, path_length);
    
    // make a function that will help us store mean of filtering distributions
    std::vector<std::function<const Mat(const Vec&)> > fs;
    auto idtyLambda = [](const Vec& x){ return Eigen::Map<const Mat>(x.data(), x.size(), 1); };
    fs.push_back(idtyLambda);

    // now run through the data and filter 
    //std::string filePath("/home/taylor/ssm/data/some_csvs/weekly_etf_data_200151223_second_baby.csv");
    std::string filePath("/home/taylor/ssm/data/some_csvs/jacq_y_data.csv");
    std::ifstream inFile(filePath);
    std::string dataLine;
    std::string oneNumber;
    unsigned int time = 0;
    // for getting variance of ll estimates
    double ll(0.0);
    while(std::getline(inFile, dataLine)) {

        // construct yt
        Vec yt(dim_obs);
        try {
            std::istringstream buff(dataLine);
            int i = 0;
            while(std::getline(buff, oneNumber, ',')){
                yt(i) = std::stod(oneNumber);
                i += 1;
            }
        } catch(const std::invalid_argument& ia) {
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }

        // update
        myMod.filterOrSmooth(yt, fs);

        // print stuff if you'd like
        std::cout << myMod.getExpectations()[0].transpose() << "\n";
        //std::cout << log(myMod.getCondLike()) << "\n";

//        Mat tmpVar( myMod.getPredictiveVar() );        
//        for(int row = 0; row < dim_obs; ++row){
//            for(int col = 0; col < dim_obs; ++col){
//                if( col == 0){
//                    std::cout << tmpVar(row,col) ;
//                }else{
//                    std::cout << ", " << tmpVar(row,col);
//                }
//            }
//            std::cout << "\n";
//        }
        //ll += log(myMod.getCondLike());
        
        // increment time
        time += 1;
    }
    //std::cout << ll << "\n";
}