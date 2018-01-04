#include "msvol_sisr.h"

MSVolSISR::MSVolSISR(int numParts, const Mat &beta, const Vec &phis, const Vec &mus, const Vec &sigmas, 
                     SISRResampStyle resampTechnique, int pathLength, double percEss)
    : SISRFilter(numParts, resampTechnique, percEss), m_beta(beta), m_phis(phis), m_mus(mus), m_sigmas(sigmas)
{
    // num factors and dimension of obs
    m_dim_obs = m_beta.rows();
    m_num_factors = m_beta.cols();     
}

MSVolSISR::~MSVolSISR(){}

Vec MSVolSISR::q1Samp(const Vec &y1)
{
    double tempStdDev;
    Vec x1samp(m_num_factors + m_dim_obs);
    for(int i = 0; i < m_num_factors + m_dim_obs; ++i)
    {
        tempStdDev = m_sigmas(i) / sqrt(1-pow(m_phis(i), 2));
        x1samp(i) = tempStdDev * m_stdNormSampler.sample()(0);
    }
    return x1samp;
}

double MSVolSISR::logMuEv(const Vec &x1)
{
    // mean 
    Vec mean = Eigen::Matrix< double, Eigen::Dynamic, 1>::Zero(m_num_factors + m_dim_obs);
    Mat cov = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Zero(m_num_factors + m_dim_obs, m_num_factors + m_dim_obs);
    for(int i = 0; i < m_num_factors + m_dim_obs; ++i){
        cov(i,i) = (pow(m_sigmas(i),2)) / (1-pow(m_phis(i), 2));
    }

    return densities::evalMultivNorm(x1, mean, cov, true);
}
    
double MSVolSISR::logQ1Ev(const Vec &x1, const Vec &y1)
{
    Vec mean = Eigen::Matrix< double, Eigen::Dynamic, 1>::Zero(m_num_factors + m_dim_obs);
    Mat cov = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Zero(m_num_factors + m_dim_obs, m_num_factors + m_dim_obs);
    for(int i = 0; i < m_num_factors + m_dim_obs; ++i){
        cov(i,i) = ( pow(m_sigmas(i),2) ) / (1-pow(m_phis(i), 2));
    }
    
    return densities::evalMultivNorm(x1, mean, cov, true); 
}

double MSVolSISR::logGEv(const Vec &yt, const Vec &xt)
{
    Vec mean = Eigen::Matrix< double, Eigen::Dynamic, 1>::Zero(m_dim_obs);
    Mat tmpMat1 = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Zero(m_num_factors, m_num_factors);
    Mat tmpMat2 = Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >::Zero(m_dim_obs, m_dim_obs);
    for(int i = 0; i < m_num_factors+m_dim_obs; ++i){
        // factors first
        if (i < m_num_factors){ 
            tmpMat1(i,i) = exp(xt(i));
        } 
        else { // obserational log variances
            tmpMat2(i-m_num_factors, i-m_num_factors) = exp(xt(i));
        }
    }
    Mat cov = m_beta * tmpMat1 * m_beta.transpose() + tmpMat2;
    return densities::evalMultivNorm(yt, mean, cov, true); 
}

double MSVolSISR::logQEv(const Vec &xt, const Vec &xtm1, const Vec &yt)
{
    Vec mean(m_num_factors + m_dim_obs);
    for (int i = 0; i < m_num_factors + m_dim_obs; ++i){
        if ( i < m_num_factors){ // factors 
            mean(i) = m_phis(i) * xtm1(i);
        }
        else { // observation log variances
            mean(i) = m_mus(i - m_num_factors) + m_phis(i) * ( xtm1(i) - m_mus(i - m_num_factors) );
        }
    }
    Mat cov = m_sigmas.array().pow(2).matrix().asDiagonal();
    return densities::evalMultivNorm(xt, mean, cov, true);
}

double MSVolSISR::logFEv(const Vec &xt, const Vec &xtm1)
{
    // same as above function qEv
    Vec mean(m_num_factors + m_dim_obs);
    for (int i = 0; i < m_num_factors + m_dim_obs; ++i){
        if ( i < m_num_factors){ // factors 
            mean(i) = m_phis(i) * xtm1(i);
        }
        else { // observation log variances
            mean(i) = m_mus(i - m_num_factors) + m_phis(i) * ( xtm1(i) - m_mus(i - m_num_factors) );
        }
    }
    Mat cov = m_sigmas.array().pow(2).matrix().asDiagonal();
    return densities::evalMultivNorm(xt, mean, cov, true);
}

Vec MSVolSISR::qSamp(const Vec &xtm1, const Vec &yt)
{
    double tempStdDev;
    double tempMean;
    Vec samp(m_num_factors + m_dim_obs);
    for(int i = 0; i < m_num_factors + m_dim_obs; ++i)
    {
        if (i < m_num_factors) // factor parts
            tempMean = m_phis(i) * xtm1(i);
        else
            tempMean = m_mus(i-m_num_factors) + m_phis(i) * ( xtm1(i) - m_mus(i-m_num_factors) );
            
        // WOWWW tempStdDev = m_sigmas(i) / sqrt(1-pow(m_phis(i), 2));
        tempStdDev = m_sigmas(i);
        samp(i) = tempMean + tempStdDev * m_stdNormSampler.sample()(0);
    }
    return samp;
}