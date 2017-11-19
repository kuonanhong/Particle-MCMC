#include "transformations.h"

#include <stdexcept>
#include <iostream>
#include <cmath>


double transformations::twiceFisher(const double &phi)
{
    if ( (phi <= -1.0) || (phi >= 1.0) )
        throw std::invalid_argument( "error: phi was not between -1 and 1" );
    else
        return std::log(1.0 + phi) - std::log(1.0 - phi);
}


double transformations::invTwiceFisher(const double &psi)
{
    double ans = (1.0 - std::exp(psi)) / ( -1.0 - std::exp(psi) );
    
    if ( (ans <= -1.0) || (ans >= 1.0) )
        std::cerr << "error: there was probably overflow for exp(psi) \n";
    
    return ans;    
}


double transformations::logit(const double &p)
{
    if ( (p <= 0.0) || (p >= 1.0))
        std::cerr << "error: p was not between 0 and 1 \n";
    
    return std::log(p) - std::log(1.0 - p);
}


double transformations::inv_logit(const double &r)
{
    double ans = 1.0/( 1.0 + std::exp(-r) );
    
    if ( (ans <= 0.0) || (ans >= 1.0))
        std::cerr << "error: there was probably underflow for exp(-r) \n";
    
    return ans;
}
