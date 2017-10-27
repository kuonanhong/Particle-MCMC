#include "noisy_ar1_filter.h"

// xt = .91 xtm1 + wt where var(wt) = 1.0
// yt = xt + vt where var(vt) = 1.5^2

//constructors and destructors
NAr1Filter::NAr1Filter(int numParts, double osd, double ssd, double a, 
                        SISRResampStyle resampTechnique, int pathLength) : 
        SISRFilter(numParts, resampTechnique, pathLength), m_obsSd(osd), m_stateSd(ssd), m_alpha(a) 
{
}

NAr1Filter::~NAr1Filter() {}

// methods
double NAr1Filter::logQ1Ev  (const Vec &x1, const Vec &y1)
{
    Vec mean(1);
    mean(0) = 0.0;
    Mat cov(1,1);
    cov(0,0) = std::pow(m_stateSd, 2) / (1.0- std::pow(m_alpha,2)); //sigma^2 / (1-alpha^2)
    return densities::evalMultivNorm(x1, mean, cov, true);
}

double NAr1Filter::logMuEv  (const Vec &x1) // same as above
{
    Vec mean(1);
    mean(0) = 0.0;
    Mat cov(1,1);
    cov(0,0) = std::pow(m_stateSd,2) / (1.0- std::pow(m_alpha,2)); //sigma^2 / (1-alpha^2)
    return densities::evalMultivNorm(x1, mean, cov, true);
}

double NAr1Filter::logGEv   (const Vec &yt, const Vec &xt)
{
    Mat cov(1,1);
    cov(0,0) = std::pow(m_obsSd, 2); //gamma=1.5
    return densities::evalMultivNorm(yt, xt, cov, true);
}

double NAr1Filter::logQEv   (const Vec &xt, const Vec &xtm1, const Vec &yt)
{
    Mat cov(1,1);      // bootstrap
    cov(0,0) = std::pow(m_stateSd, 2); 
    return densities::evalMultivNorm(xt, m_alpha*xtm1, cov, true);
}

double NAr1Filter::logFEv   (const Vec &xt, const Vec &xtm1) // same as above
{
    Mat cov(1,1);
    cov(0,0) = pow(m_stateSd, 2); 
    return densities::evalMultivNorm(xt, m_alpha*xtm1, cov, true);
}

Vec NAr1Filter::qSamp (const Vec &xtm1, const Vec &yt)
{
    return m_alpha*xtm1 + m_stateSd * m_stdNormSampler.sample(); 
}

Vec NAr1Filter::q1Samp(const Vec &y1)
{
    Vec x1samp = m_stdNormSampler.sample();
    x1samp *= m_stateSd;
    return x1samp / sqrt( 1-pow(m_alpha, 2) );
}

bool NAr1Filter::isSdNeg()
{
    if (m_obsSd <= 0 || m_stateSd <= 0)
        return true;
    else
        return false;
}