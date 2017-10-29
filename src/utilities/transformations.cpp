#include "transformations.h"

#include <stdexcept>
#include <iostream>
#include <cmath>

double transformations::twiceFisher(const double &phi)
{
    if ( (phi <= -1.0) || (phi >= 1.0) )
        throw std::invalid_argument( "inappropriate range for argument" );
    else
        return std::log(1.0 + phi) - std::log(1.0 - phi);
}

double transformations::invTwiceFisher(const double &psi)
{
    double ans = (1.0 - std::exp(psi)) / ( -1.0 - std::exp(psi) );
    
    if ( (ans <= -1.0) || (ans >= 1.0) )
        std::cerr << "ERROR\n";
    
    return ans;    
}
