#include "jacquier_et_al_apf.h"
//#include <cmath>

JacEtAlAPF::JacEtAlAPF(int numParts, const Mat &beta, const Vec &R_sigmas_as_vec, const Vec &mus, const Vec &phis, const Vec &sigmas, APFResampStyle resampTechnique, int pathLength)
    : APFFilter(numParts, resampTechnique), m_beta(beta), m_R_sigmas(R_sigmas_as_vec), m_mus(mus), m_phis(phis), m_sigmas(sigmas)
{
    // num factors and dimension of obs
    m_dim_obs = m_beta.rows();
    m_num_factors = m_beta.cols();    
    
    // set up time one sampler guy
    m_timeOneSampler.setMean(Vec::Zero(m_num_factors));
    Mat sigma0( ( m_sigmas.array().square() / (1.0 - m_phis.array().square()) ).matrix().asDiagonal() );
    m_timeOneSampler.setCovar(sigma0);
    
    // and set up trans jump sampler
    m_transJumpSampler.setMean(Vec::Zero(m_num_factors));
    m_transJumpSampler.setCovar(m_sigmas.array().square().matrix().asDiagonal());
}


JacEtAlAPF::~JacEtAlAPF(){}


//Mat JacEtAlAPF::getPredictiveVar() const 
//{
//    auto getExpectExponent =  [this] ( const Vec &xt ) -> Vec { return ( m_mus + m_phis*(xt - m_mus) + m_sigmas.array().square().matrix() / 2.0 ).array().exp().matrix(); } ;
//    Vec lil_guy = getMeanFunc(getExpectExponent);
//    return m_beta * lil_guy.asDiagonal() * m_beta.transpose() + Mat(m_R_sigmas.array().square().matrix().asDiagonal());
//}


Vec JacEtAlAPF::q1Samp(const Vec &y1)
{
    return m_mus + m_timeOneSampler.sample();
}


double JacEtAlAPF::logMuEv(const Vec &x1)
{

    // construct appropriate cov mat
    Mat cov ( ( m_sigmas.array().square() / (1.0 - m_phis.array().square()) ).matrix().asDiagonal() );
    
    // return the density with the correct mean
    return densities::evalMultivNorm(x1, m_mus, cov, true);
}


double JacEtAlAPF::logQ1Ev(const Vec &x1, const Vec &y1)
{
    
    Mat cov ( ( m_sigmas.array().square() / (1.0 - m_phis.array().square()) ).matrix().asDiagonal() );
    return densities::evalMultivNorm(x1, m_mus, cov, true); 
}


double JacEtAlAPF::logGEv(const Vec &yt, const Vec &xt) 
{
    // construct covariance matrix
    Mat cov ( m_beta * xt.array().exp().matrix().asDiagonal() * m_beta.transpose() );
    cov += m_R_sigmas.array().square().matrix().asDiagonal();
    
    // return density evaluation
    return  densities::evalMultivNorm(yt, Vec::Zero(m_dim_obs), cov, true);
}


Vec JacEtAlAPF::propMu(const Vec &xtm1)
{
    return m_mus + m_phis.asDiagonal() * (xtm1 - m_mus);
}
    
    
Vec JacEtAlAPF::fSamp (const Vec &xtm1)  
{
    return m_mus + m_phis.asDiagonal() * (xtm1 - m_mus) + m_transJumpSampler.sample();
}


