#include "jacetal_apf_sisr_compare.h"

#include "jacquier_et_al_apf.h" // apf model
#include "jacquier_et_al.h" //sisr model
#include "convenience_funcs.h" // readInData

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
    std::vector<Vec> data = convenience_funcs::readInData("/home/taylor/ssm/data/some_csvs/weekly_etf_data_200151223_201451", 9);

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

