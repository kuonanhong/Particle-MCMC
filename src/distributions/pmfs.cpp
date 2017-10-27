#include <cmath> // tgamma, pow, exp, otherstuff
#include <exception>

#include "pmfs.h"

////////////////////////////////////////////////
/////////         Evaluators           /////////
////////////////////////////////////////////////

double pmfs::evalDiscreteUnif(const int &x, const int &k, bool log)
{
    if( (1 <= x) && (x <= k) ){
        if(log){
            return -std::log(static_cast<double>(k));
        }else{
            return 1.0 / static_cast<double>(k);
        }
    }else{ // x not in support
        if(log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}


////////////////////////////////////////////////
/////////           samplers           /////////
////////////////////////////////////////////////

//////////////// DiscreteUnifSampler
pmfs::DiscreteUnifSampler::DiscreteUnifSampler(const int& k) :
        m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())}, 
        m_int_gen(1, k)
{    
}


int pmfs::DiscreteUnifSampler::sample()
{
    return m_int_gen(m_rng);
}


/////////////// UniformSampler
pmfs::DiscreteCustomSampler::DiscreteCustomSampler(const std::vector<double>& weights) : 
        m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())},
        m_int_gen(weights.begin(), weights.end())
{    
}


int pmfs::DiscreteCustomSampler::sample()
{
    return m_int_gen(m_rng) + 1;
}
