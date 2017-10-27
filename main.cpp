#include <iostream>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

//#include "kfiltertest.h"
//#include "sVolTest.h"
//#include "noisyAr1Test.h"
//#include "sVolAPFTest.h"
//#include "noisyar1comparison.h"
//#include "pmcmc_trial.h"
//#include "msvol_test.h"
//#include "pmcmc_msvol_pitt_shephard.h"
//#include "jac_filt_test.h"
//#include "jac_filt_test_apf.h" 
//#include "pmcmc_jacetal_apf.h"
//#include "pmcmc_jacetal_sisr.h"
//#include "jacetal_apf_sisr_compare.h"
#include "do_pmmh_jacetal.h"
//#include "do_ada_pmmh_jacetal.h"

int main(int argc, char **argv)
{

    // remove this if you want a gmon.out
    //clock_t begin = clock();  
    auto begin = Clock::now();
    
    

    //kFilterTest();
    //sVolTest();
    //noisyAr1Test();
    //sVolAPFTest();
    //for(int i = 0; i < 300; ++i){
        //noisyAr1Comparison();
    //}
    //noisyAr1Comparison();
    //pmcmc_trial::testpmcmc(5, 
    //                       1000, 
    //                       "/home/taylor/ssm/data/some_csvs/noisy_ar1_data.csv", 
    //                       "/home/taylor/Desktop/temptestpmcmcout");  // long time!
    //msvol_test();
    //PittShep::commenceSampling("/home/taylor/Desktop/tmpOut", 1000, 1000);  //2032.78minuites or 33.87967 hours or 1.4 days
    //pmmh_rbpf_mod1::samplePmmhRbpfMod1("/home/taylor/Desktop/tmpOut", 1000, 10000); // takes about  153.934 min
    //pmmh_sisr_mod1::samplePmmhSisrMod1("/home/taylor/Desktop/tmpOut", 1000, 10000);
    //jac_filt_test();
    //jac_filt_test_apf(); 
    //Pmcmc_JacEtAl_APF::commenceSampling("/home/taylor/ssm/data/some_csvs/weekly_etf_data_200151223_baby.csv", //jacq_y_data.csv
    //                                    "/home/taylor/Desktop/tmpOut", 
    //                                    "/home/taylor/Desktop/tmpMessages", 
    //                                    1300, 
    //                                    10000);
    //Pmcmc_JacEtAl_SISR::commenceSampling("/home/taylor/ssm/data/some_csvs/weekly_etf_data_200151223_201451", //jacq_y_data.csv
    //                                    "/home/taylor/Desktop/tmpOut", 
    //                                    "/home/taylor/Desktop/tmpMessages", 
    //                                    1000, 
    //                                    5000);
    //jacetal_apf_sisr_compare();
    do_pmmh_jacetal();  
    //do_ada_pmmh_jacetal(); 
    




    

    // remove this stuff if you want a gmon.out
    //clock_t end = clock();
    auto end = Clock::now();
    std::cout << "Delta t2-t1: " 
              << std::chrono::duration_cast<std::chrono::minutes>(end - begin).count()
              << " minutes" << std::endl;
    
//    double elapsed_ms = double(end - begin) / CLOCKS_PER_SEC * 1000; //number milliseconds
//    std::cout << "elapsed time in ms: " << elapsed_ms << std::endl;
//    std::cout << "elapsed time in minutes: " << elapsed_ms / (1000*60) << std::endl;

	return 0;
}
