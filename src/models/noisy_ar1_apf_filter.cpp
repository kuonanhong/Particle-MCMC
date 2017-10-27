#include <noisy_ar1_apf_filter.h>

// xt = .91 xtm1 + wt where var(wt) = 1.0
// yt = xt + vt where var(vt) = 1.5^2

//constructors destructors
NAr1APFFilter::NAr1APFFilter(int numParts,double osd, double ssd,double a, APFResampStyle resampTechnique, int pathLength) 
        : APFFilter(numParts, resampTechnique, pathLength), m_alpha(a), m_stateSd(ssd), m_obsSd(osd)
{
}

NAr1APFFilter::~NAr1APFFilter(){}

// methods
double NAr1APFFilter::logQ1Ev  (const Vec &x1, const Vec &y1)
{
    Vec mean(1);
    mean(0) = 0.0;
    Mat cov(1,1);
    cov(0,0) = std::pow(m_stateSd, 2) / ( 1.0- std::pow(m_alpha,2) ); 
    return densities::evalMultivNorm(x1, mean, cov, true);
}

double NAr1APFFilter::logMuEv (const Vec &x1) // same as above
{
    Vec mean(1);
    mean(0) = 0.0;
    Mat cov(1,1);
    cov(0,0) = std::pow(m_stateSd, 2) / ( 1.0- std::pow(m_alpha,2) ); 
    return densities::evalMultivNorm(x1, mean, cov, true);
}

double NAr1APFFilter::logGEv (const Vec &yt, const Vec &xt)
{
    Mat cov(1,1);
    cov(0,0) = std::pow(m_obsSd, 2); 
    return densities::evalMultivNorm(yt, xt, cov, true);
}

Vec NAr1APFFilter::q1Samp(const Vec &y1)
{
    Vec samp = m_stdNormSampler.sample();
    samp *= m_stateSd;
    samp /= std::sqrt( 1.0- std::pow(m_alpha, 2) );
    return samp;
}

Vec NAr1APFFilter::propMu(const Vec& xtm1)
{
    return m_alpha*xtm1; //le mode
}

Vec NAr1APFFilter::fSamp(const Vec &xtm1)
{
        return m_alpha*xtm1 + m_stateSd * m_stdNormSampler.sample(); 
}