#include <iostream>
#include <chrono>

typedef std::chrono::high_resolution_clock Clock;

//#include "kfilter_test.h"
//#include "hmm_test.h"
#include "svol_test.h"
//#include "noisy_ar1_test.h"
//#include "svol_apf_test.h"
//#include "noisy_ar1_comparison.h"
//#include "msvol_test.h"
//#include "jac_filt_test.h"
//#include "jac_filt_test_apf.h" 
//#include "jacetal_apf_sisr_compare.h"
//#include "do_pmmh_svol.h"
//#include "do_pmmh_jacetal.h"

//#include "do_ada_pmmh_svol.h"
//#include "do_ada_pmmh_jacetal.h"
//#include "do_ada_pmmh_msl_apf.h"


int main(int argc, char **argv)
{

    // remove this if you want a gmon.out
    //clock_t begin = clock();  
    //auto begin = Clock::now();
    

    //kFilterTest();
    //hmmTest();
    sVolTest();
    //noisyAr1Test();
    //sVolAPFTest();
    //for(int i = 0; i < 300; ++i){
        //noisyAr1Comparison();
    //}
    //noisyAr1Comparison();
    //msvol_test();
    //jac_filt_test();
    //jac_filt_test_apf(); 
    //jacetal_apf_sisr_compare();
    //do_pmmh_svol();
    //do_ada_pmmh_svol();
    //do_pmmh_jacetal();  
    //do_ada_pmmh_jacetal(); 
    //do_ada_pmmh_msl_apf();

    // remove this stuff if you want a gmon.out
    //clock_t end = clock();
//    auto end = Clock::now();
//    std::cout << "Delta t2-t1: " 
//              << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count()
//              << " seconds" << std::endl;

	return 0;
}
