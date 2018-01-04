#include <iostream>
#include <string> //string
#include <fstream> //ifstream, ofstream
#include "msvol_sisr.h" // model class

void msvol_test()
{
    ////////////////////////////////////
    // algorithm and example variables /
    ////////////////////////////////////
    int num_parts = 50;
    
    /////////////////////////////////////////////////////////////////
    // instantiate model object with known parameters               /
    // make sure these parameters cohere with                       /
    // /home/taylor/ssm/data/simulate_scripts/gen_pitt_shep_data.py /
    /////////////////////////////////////////////////////////////////
    int num_factors = 1;
    int dim_obs = 9;
     

    // this was a particular configuration giving NaNs (via f, which uses the parameters phis and mus)
    Mat beta(dim_obs, num_factors);    
//    beta << 1, 2.48004, 4.18123, 1.13329, 0.986817, 1.67048, 0.708067, 0.450667, 2.63194;
    Vec phis(dim_obs + num_factors);
//    phis << 0.966727, -0.868314, -0.477499, 0.839821, 0.984374, 0.799178, 0.989602, 0.768566, 0.994537, 0.864411;
    Vec mus(dim_obs);
//    mus << 0.669864, 0.40057, 1.35178, -1.14064, 1.35712, -3.60199, 1.53658, 0.0787653, -0.833408;
    Vec sigmas(dim_obs + num_factors);
//    sigmas << 0.224666, 0.459438, 0.205295, 0.378066, 0.174906, 0.597562, 0.0265288, 0.0747781, 0.0458934, 0.409963;
    
    for(int i = 0; i < dim_obs + num_factors; ++i){
        phis(i) = .92;
        sigmas(i) = .44; // sqrt .2
        if (i >= num_factors){
            mus(i - num_factors) = 1.0; 
            beta(i - num_factors,0) = 1.0;
        }
    }
    SISRResampStyle rt = SISRResampStyle::everytime_multinomial;
    int path_length = 0; // filtering only
    MSVolSISR myMod(num_parts, beta, phis, mus, sigmas, rt, path_length);
    
    // lambdas that help get the filter mean
    std::vector<std::function<const Mat(const Vec&)> > fs;
    auto idtyLambda = [](const Vec& x) { return Eigen::Map<const Mat>(x.data(), x.size(), 1); };
    fs.push_back(idtyLambda);
    
    ////////////////////////////////////////
    // now run through the data and filter /
    ////////////////////////////////////////
    std::string filePath("/home/taylor/ssm/data/some_csvs/msvol_y_data.csv");
    std::ifstream inFile(filePath);
    std::string dataLine;
    std::string oneNumber;
    unsigned int time = 0;
    // usethis if want a total log like 
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
        myMod.filter(yt, fs);

        // print stuff if you'd like
        //std::cout << myMod.getExpectations()[0].transpose() << "\n";
        std::cout << myMod.getLogCondLike() << "\n";
        //ll += myMod.getLogCondLike();
        
        // increment time
        time += 1;
    }
    //std::cout << "loglike: " << ll << "\n";
}