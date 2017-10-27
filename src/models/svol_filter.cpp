#include "svol_filter.h"

SVolFilter::SVolFilter(int numParts) 
    : SISRFilter(numParts)
{
}

SVolFilter::~SVolFilter()
{
}

Vec SVolFilter::q1Samp(const Vec &y1)
{
    Vec x1samp = m_stdNormSampler.sample();
    return x1samp * sqrt(1./(1-.91*.91)); // variance is sigma^2/(1-alpha^2)
}

Vec SVolFilter::qSamp(const Vec &xtm1, const Vec &yt)
{
    Vec xtsamp = m_stdNormSampler.sample();
    xtsamp(0) += .91 * xtm1(0); //add the mean
    return xtsamp;
}


double SVolFilter::logGEv(const Vec &yt, const Vec &xt)
{
    Vec mu(1);
    mu(0) = 0;
    Mat sigma(1,1);
    sigma(0,0) = .25 * exp(xt(0)); //beta^2 * e^{xt}
    return densities::evalMultivNorm(yt, mu, sigma, true);
}

double SVolFilter::logFEv(const Vec &xt, const Vec &xtm1)
{
    Vec mu(1);
    mu(0) = .91 * xtm1(0); //alpha * xmt1
    Mat sigma(1,1);
    sigma(0,0) = 1.0; //sigma^2
    return densities::evalMultivNorm(xt, mu, sigma, true);
}

double SVolFilter::logMuEv(const Vec &x1)
{
    Vec mu(1);
    mu(0) = 0;
    Mat sigma(1,1);
    sigma(0,0) = 1./(1-.91*.91); //sigma^2/(1-alpha^2)
    return densities::evalMultivNorm(x1, mu, sigma, true);
}

double SVolFilter::logQ1Ev(const Vec &x1samp, const Vec &y1)
{
    return logMuEv(x1samp);
}

double SVolFilter::logQEv(const Vec &xt, const Vec &xtm1, const Vec &yt)
{
    return logFEv(xt, xtm1); // bootstrap filter
}