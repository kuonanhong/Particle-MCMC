#include "jacetal_apf_sisr_compare.h"
#include "jacquier_et_al_apf.h" // apf model
#include "jacquier_et_al.h" //sisr model


std::vector<Vec> read_data(const std::string &dataFile, int num_obs, int size_obs)
{
    // return this
    std::vector<Vec> data(num_obs, Vec(size_obs)); // dataset is 1000 observations each 9-d

    // read in data and store in object
    std::string line;
    std::ifstream inFile(dataFile);
    std::string oneNumberOnOneLine;
    int time(0);
    int elemNum(0);
    while ( std::getline(inFile, line) ){ // every row of data
    
        elemNum = 0;
        try{
            std::istringstream buff(line);
            while(std::getline(buff, oneNumberOnOneLine, ',')){
                data[time](elemNum) = std::stod(oneNumberOnOneLine);
                elemNum++;
            }
        } catch (const std::invalid_argument& ia){
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }          
        time++;
    }
    
    // now return the data object
    return data;
}


void unFlattenParams(Mat &beta, Vec &phi, Vec &mu, Vec &sigma, Vec &RstdDevs, const std::vector<double> &flatOnes)
{   
    
    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    int c = 0;
    
    // betas 
    for(int row = 0; row < beta.rows(); ++row){
        for(int col = 0; col < beta.cols(); ++col){
            beta(row,col) = flatOnes[c];
            c++;
        }
    }
    
    // phis (only one)
    phi(0) = flatOnes[c];
    c++;
    
    // mus (only one)
    mu(0) = flatOnes[c];
    c++;
        
    // sigmas (only one)
    sigma(0) = flatOnes[c];
    c++;
    
    // R std devs 
    for(int row = 0; row < RstdDevs.rows(); ++row){
        RstdDevs(row) = flatOnes[c];
        c++;
    }
}


void jacetal_apf_sisr_compare()
{
    APFResampStyle rt_apf = APFResampStyle::everytime_multinomial;
    SISRResampStyle rt_sisr = SISRResampStyle::everytime_multinomial;
    int np = 200;
        
    // read in data 
    std::vector<Vec> data = read_data("/home/taylor/ssm/data/some_csvs/weekly_etf_data_200151223_201451", 436, 9);

    // ordering is: betas(9), phis(1),  mus(1), sigma(1), R_std_dev_vec(9)
    std::vector<double> currentThetas {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                        .7, 0.0, .1,//0.00001,
                                        .75, .75, .75, .75, .75, .75, .75, .75, .75};    
                
    // flatten parameters
    Mat beta(9, 1); 
    Vec phi(1);
    Vec mu(1);
    Vec sigma(1);
    Vec Rsigmas(9);
    unFlattenParams(beta, phi, mu, sigma, Rsigmas, currentThetas);    


    // make two models
    JacEtAlAPF modAPF(np, beta, Rsigmas, mu, phi, sigma, rt_apf); 
    JacEtAl    modSISR(np, beta, Rsigmas, mu, phi, sigma, rt_sisr); 
    
    //iterate over time
    std::cout << "APFLogLike, SISRLogLike\n"; 
    for(int row = 0; row < data.size(); ++row){
        modAPF.filterOrSmooth(data[row]);  
        modSISR.filterOrSmooth(data[row]);
        std::cout << modAPF.getLogCondLike() << ", " << modSISR.getLogCondLike() << "\n";
    }
    
}

