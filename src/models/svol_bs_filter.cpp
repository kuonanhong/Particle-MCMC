#include "svol_bs_filter.h"

#include <cmath>


SVolBSFilter::SVolBSFilter(unsigned numParts, double beta, double phi, double sigma, unsigned pLen) 
    : m_beta(beta), m_phi(phi), m_sigma(sigma), 
    BSFilter(numParts, BSResampStyle::everytime_multinomial, pLen)
{
}

SVolBSFilter::~SVolBSFilter()
{
}


Vec SVolBSFilter::q1Samp(const Vec &y1)
{
    Vec x1samp = m_stdNormSampler.sample();
    return x1samp * m_sigma / std::sqrt(1.-m_phi*m_phi); // variance is sigma^2/(1-phi^2)
}


Vec SVolBSFilter::qSamp(const Vec &xtm1, const Vec &yt)
{
    Vec xtsamp = m_stdNormSampler.sample();
    xtsamp(0) += m_phi * xtm1(0); //add the mean
    return xtsamp;
}


double SVolBSFilter::logGEv(const Vec &yt, const Vec &xt)
{
    Vec mu(1);
    mu(0) = 0;
    Mat varMat(1,1);
    varMat(0,0) = m_beta * m_beta * std::exp(xt(0)); //beta^2 * e^{xt}
    return densities::evalMultivNorm(yt, mu, varMat, true);
}


double SVolBSFilter::logMuEv(const Vec &x1)
{
    Vec mu(1);
    mu(0) = 0;
    Mat covMat(1,1);
    covMat(0,0) = m_sigma* m_sigma/(1.-m_phi*m_phi); //sigma^2/(1-phi^2)
    return densities::evalMultivNorm(x1, mu, covMat, true);
}


double SVolBSFilter::logQ1Ev(const Vec &x1samp, const Vec &y1)
{
    return logMuEv(x1samp);
}


double SVolBSFilter::logQEv(const Vec &xt, const Vec &xtm1, const Vec &yt)
{
    Vec mu(1);
    mu(0) = m_phi * xtm1(0); //phi * xmt1
    Mat covMat(1,1);
    covMat(0,0) = m_sigma * m_sigma; //sigma^2
    return densities::evalMultivNorm(xt, mu, covMat, true);
}