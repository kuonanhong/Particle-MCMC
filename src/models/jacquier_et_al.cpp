#include "jacquier_et_al.h"
              
JacEtAl::JacEtAl(int numParts, const Mat &beta, const Vec &RAsVec, const Vec &mus, const Vec &phis, const Vec &sigmas, 
                     SISRResampStyle resampTechnique, int pathLength, double percEss)
    : SISRFilter(numParts, resampTechnique, percEss), m_beta(beta), m_R(RAsVec), m_mus(mus), m_phis(phis), m_sigmas(sigmas)
{
    // num factors and dimension of obs
    m_dim_obs = m_beta.rows();
    m_num_factors = m_beta.cols();    

    // set up sampler
    m_timeOneSampler.setMean(Vec::Zero(m_num_factors));
    Mat sigma0( (m_sigmas.array().square() / (1.0 - m_phis.array().square())).matrix().asDiagonal() );
    m_timeOneSampler.setCovar(sigma0);
    
    // and set up trans jump sampler
    m_transJumpSampler.setMean(Vec::Zero(m_num_factors));
    m_transJumpSampler.setCovar(m_sigmas.array().square().matrix().asDiagonal());
}

JacEtAl::~JacEtAl(){}

Vec JacEtAl::q1Samp(const Vec &y1)
{    
    return m_mus + m_timeOneSampler.sample();
}

double JacEtAl::logMuEv(const Vec &x1) 
{
    Mat cov ( ( m_sigmas.array().square() / (1.0 - m_phis.array().square()) ).matrix().asDiagonal() );
    return densities::evalMultivNorm(x1, m_mus, cov, true);
}
    
double JacEtAl::logQ1Ev(const Vec &x1, const Vec &y1)
{
    Mat cov ( ( m_sigmas.array().square() / (1.0 - m_phis.array().square()) ).matrix().asDiagonal() );
    return densities::evalMultivNorm(x1, m_mus, cov, true); 
}

double JacEtAl::logGEv(const Vec &yt, const Vec &xt)
{
    Mat cov ( m_beta * xt.array().exp().matrix().asDiagonal() * m_beta.transpose() );
    cov += m_R.asDiagonal();
    return  densities::evalMultivNorm(yt, Vec::Zero(m_dim_obs), cov, true);
}

double JacEtAl::logQEv(const Vec &xt, const Vec &xtm1, const Vec &yt)
{
    Vec mean( m_mus + m_phis.asDiagonal() * (xtm1 - m_mus) );
    Mat cov( m_sigmas.array().square().matrix().asDiagonal() );
    return densities::evalMultivNorm(xt, mean, cov, true);
}

double JacEtAl::logFEv(const Vec &xt, const Vec &xtm1)
{
    // make the mean
    Vec mean( m_mus + m_phis.asDiagonal() * (xtm1 - m_mus) );
    Mat cov( m_sigmas.array().square().matrix().asDiagonal() );
    return densities::evalMultivNorm(xt, mean, cov, true);
}

Vec JacEtAl::qSamp(const Vec &xtm1, const Vec &yt)
{
    return m_mus + m_phis.asDiagonal() * (xtm1 - m_mus) + m_transJumpSampler.sample();
}