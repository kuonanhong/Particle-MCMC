#include <cmath> // tgamma, pow, exp, otherstuff
#include <exception>

#include "densities.h"

////////////////////////////////////////////////
/////////         Evaluators           /////////
////////////////////////////////////////////////


double densities::evalMultivNorm(const Vec &x, const Vec &meanVec, const Mat &covMat, bool log)
{
    double quadform  = (x - meanVec).transpose() * covMat.inverse() * (x-meanVec);

    if (log){

        // calculate log-determinant using cholesky decomposition (assumes symmetric and positive definite)
        double ld (0.0);
        Eigen::LLT<Mat> lltM(covMat);  
        Mat L = lltM.matrixL(); // the lower diagonal L such that M = LL^T

        // add up log of diagnols of Cholesky L
        for(unsigned int i = 0; i < covMat.rows(); ++i){
            ld += std::log(L(i,i));
        }
        ld *= 2; // covMat = LL^T

        return -.5*log_two_pi * covMat.rows() - .5*ld - .5*quadform;


    }else{  // not the log density
        double normConst = pow(inv_sqrt_2pi, covMat.rows()) * pow(covMat.determinant(), -.5);
        return normConst * exp(-.5* quadform);
    }

}


double densities::evalUnivBeta(const double &x, const double &alpha, const double &beta, bool log)
{
    if( (x > 0.0) && (x < 1.0) && (alpha > 0.0) && (beta > 0.0) ){ // x in support and parameters acceptable
        if(log){
            return std::lgamma(alpha + beta) - std::lgamma(alpha) - std::lgamma(beta) + (alpha - 1.0)*std::log(x) + (beta - 1.0) * std::log(1.0 - x);
        }else{
            return pow(x, alpha-1.0) * pow(1.0-x, beta-1.0) * std::tgamma(alpha + beta) / ( std::tgamma(alpha) * std::tgamma(beta) );
        }

    }else{ //not ( x in support and parameters acceptable )
        if(log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}


double densities::evalUnivInvGamma(const double &x, const double &alpha, const double &beta, bool log)
{
    if ( (x > 0.0) && (alpha > 0.0) && (beta > 0.0) ){ // x in support and acceptable parameters
        if (log){
            return alpha * std::log(beta) - std::lgamma(alpha) - (alpha + 1.0)*std::log(x) - beta/x;
        }else{
            return pow(x, -alpha-1.0) * exp(-beta/x) * pow(beta, alpha) / std::tgamma(alpha);
        }
    }else{ // not ( x in support and acceptable parameters )
        if (log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}

double densities::evalUnivHalfNorm(const double &x, const double &sigmaSqd, bool log)
{
    if( (x >= 0.0) && (sigmaSqd > 0.0)){
        if (log){
            return .5*log_two_over_pi - .5*std::log(sigmaSqd) - pow(x, 2) / (2*sigmaSqd);
        }else{
            return exp(-pow(x,2)/(2*sigmaSqd) ) * sqrt_two_over_pi / sqrt(sigmaSqd);
        }
    }else{
        if (log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}


////////////////////////////////////////////////
/////////           samplers           /////////
////////////////////////////////////////////////


densities::MVNSampler::MVNSampler()
        : m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())}, 
          m_z_gen(0.0, 1.0)
{
    Vec zero(1);
    zero(0) = 0.0;
    setMean(zero);
    
    Mat one(1,1);
    one(0,0) = 1.0;
    setCovar(one);  
}


densities::MVNSampler::MVNSampler(const Vec &meanVec, const Mat &covMat)
        : m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())}, 
          m_z_gen(0.0, 1.0)
{
    setCovar(covMat);
    setMean(meanVec);
}


void densities::MVNSampler::setCovar(const Mat &covMat)
{
    Eigen::SelfAdjointEigenSolver<Mat> eigenSolver(covMat);
    m_scale_mat = eigenSolver.eigenvectors() * eigenSolver.eigenvalues().cwiseMax(0).cwiseSqrt().asDiagonal();
}


void densities::MVNSampler::setMean(const Vec &meanVec)
{
    m_mean = meanVec;
}


Vec densities::MVNSampler::sample()
{
    Vec Z(m_mean.rows());
    for (int jj=0; jj< m_mean.rows(); ++jj) 
    {
        Z(jj) = m_z_gen(m_rng);
    }
    return m_mean + m_scale_mat * Z;
}


/////////////// UniformSampler

densities::UniformSampler::UniformSampler() : 
        m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())},
        m_unif_gen(0.0, 1.0)
{
}


densities::UniformSampler::UniformSampler(double lower, double upper) : 
        m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())},
        m_unif_gen(lower, upper)
{    
}


double densities::UniformSampler::sample()
{
    return m_unif_gen(m_rng);
}
