#include "svol_apf_filter.h"

//constructors and destructors
// defaults to 0 pathlength
SVolAPFFilter::SVolAPFFilter(int numParts)//, APFResampStyle resampTechnique) 
        : APFFilter(numParts, APFResampStyle::everytime_multinomial) //defaults to everytime resampling
{
}


SVolAPFFilter::~SVolAPFFilter()
{
}


Vec SVolAPFFilter::q1Samp(const Vec &y1)
{
    Vec initial = m_stdNormSampler.sample();
    initial(0) *= std::sqrt(1.0/(1-.91*.91)); // variance is sigma^2/(1-alpha^2)
    return initial;
}


Vec SVolAPFFilter::fSamp(const Vec &xtm1)
{
    Vec initial = m_stdNormSampler.sample();
    initial(0) += .91 * xtm1(0); //add the mean
    return initial;
}


double SVolAPFFilter::logGEv(const Vec &yt, const Vec &xt)
{
    Vec mu(1);
    mu(0) = 0.0;
    Mat sigma(1,1);
    sigma(0,0) = .25 * std::exp(xt(0)); //beta^2 * e^{xt}
    return densities::evalMultivNorm(yt, mu, sigma, true);
}


double SVolAPFFilter::logFEv(const Vec &xt, const Vec &xtm1)
{
    Vec mu(1);
    mu(0) = .91 * xtm1(0); //alpha * xmt1
    Mat sigma(1,1);
    sigma(0,0) = 1.0; //sigma^2
    return densities::evalMultivNorm(xt, mu, sigma, true);
}


double SVolAPFFilter::logMuEv(const Vec &x1)
{
    Vec mu(1);
    mu(0) = 0.0;
    Mat sigma(1,1);
    sigma(0,0) = 1.0/(1-.91*.91); //sigma^2/(1-alpha^2)
    return densities::evalMultivNorm(x1, mu, sigma, true);
}


double SVolAPFFilter::logQ1Ev(const Vec &x1samp, const Vec &y1)
{
    return logMuEv(x1samp);
}


Vec SVolAPFFilter::propMu(const Vec& xtm1)
{
    return .91*xtm1; 
}

