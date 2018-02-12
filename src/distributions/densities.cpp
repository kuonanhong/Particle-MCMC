#include <cmath> // tgamma, pow, exp, otherstuff
#include <exception>

#include "transformations.h" // logit
#include "densities.h"

////////////////////////////////////////////////
/////////         Evaluators           /////////
////////////////////////////////////////////////


double densities::evalUnivNorm(const double &x, const double &mu, const double &sigma, bool log)
{
    double exponent = -.5*(x - mu)*(x-mu)/(sigma*sigma);
    if( sigma > 0.0){
        if(log){
            return -std::log(sigma) - .5*log_two_pi + exponent;
        }else{
            return inv_sqrt_2pi * std::exp(exponent) / sigma;
        }
    }else{
        if(log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}


double densities::evalMultivNorm(const Vec &x, const Vec &meanVec, const Mat &covMat, bool log)
{
    // from Eigen: Remember that Cholesky decompositions are not rank-revealing. 
    /// This LLT decomposition is only stable on positive definite matrices, 
    // use LDLT instead for the semidefinite case. Also, do not use a Cholesky 
    // decomposition to determine whether a system of equations has a solution.
    Eigen::LLT<Mat> lltM(covMat);
    Vec tmp = lltM.solve(x-meanVec); // tmp = Sigma^{-1}(x-meanVec)  iff Sigma tmp = (x-meanVec)
    double quadform = ((x-meanVec).transpose() * tmp)(0,0);

    if (log){

        // calculate log-determinant using cholesky decomposition too
        double ld (0.0);
        Mat L = lltM.matrixL(); // the lower diagonal L such that M = LL^T

        // add up log of diagnols of Cholesky L
        for(size_t i = 0; i < covMat.rows(); ++i){
            ld += std::log(L(i,i));
        }
        ld *= 2; // covMat = LL^T

        return -.5*log_two_pi * covMat.rows() - .5*ld - .5*quadform;


    }else{  // not the log density
        double normConst = std::pow(inv_sqrt_2pi, covMat.rows()) * std::pow(covMat.determinant(), -.5);
        return normConst * std::exp(-.5* quadform);
    }

}

double densities::evalMultivNormWBDA(const Vec &x, const Vec &meanVec, const Vec &A, const Mat &U, const Mat &C, bool log)
{
    Mat Ainv = A.asDiagonal().inverse();
    Mat Cinv = C.inverse();
    Mat invThing = (Cinv + U.transpose()*Ainv*U).inverse();
    Mat SigInv = Ainv - Ainv*U*invThing*U.transpose()*Ainv;
    double quadform = (x-meanVec).transpose() * SigInv * (x-meanVec);
    
    if (log){

        // calculate log-determinant using cholesky decomposition (assumes symmetric and positive definite)
        double halfld (0.0);
        Eigen::LLT<Mat> lltSigInv(SigInv);

  
        Mat L = lltSigInv.matrixL(); // the lower diagonal L such that SigInv = LL^T

        // add up log of diagnols of Cholesky L
        for(size_t i = 0; i < SigInv.rows(); ++i){
            halfld += std::log(L(i,i));
        }

        return -.5*log_two_pi * SigInv.rows() + halfld - .5*quadform;


    }else{  // not the log density
        double normConst = std::pow(inv_sqrt_2pi, SigInv.rows()) * std::pow(SigInv.determinant(), .5);
        return normConst * std::exp(-.5* quadform);
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
            return .5*log_two_over_pi - .5*std::log(sigmaSqd) - .5*x*x / sigmaSqd;
        }else{
            return std::exp(-.5*x*x/sigmaSqd) * sqrt_two_over_pi / std::sqrt(sigmaSqd);
        }
    }else{
        if (log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}


double densities::evalLogitNormal(const double &x, const double &mu, const double &sigma, bool log)
{
    if( (x >= 0.0) && (x <= 1.0) && (sigma > 0.0)){
        
        double exponent = -.5*(transformations::logit(x) - mu)*(transformations::logit(x) - mu) / (sigma*sigma);
        if(log){
            return -std::log(sigma) - .5*log_two_pi - std::log(x) - std::log(1.0-x) + exponent;
        }else{
            return inv_sqrt_2pi * std::exp(exponent) / (x * (1.0-x) * sigma);   
        }
    }else{
        if(log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}
 

double densities::evalTwiceFisherNormal(const double &x, const double &mu, const double &sigma, bool log)
{
    if( (x >= -1.0) && (x <= 1.0) && (sigma > 0.0)){
        
        double exponent = -.5*(std::log((1.0+x)/(1.0-x)) - mu)*(std::log((1.0+x)/(1.0-x)) - mu)/(sigma* sigma);
        if(log){
            return -std::log(sigma) - .5*log_two_pi + std::log(2.0) - std::log(1.0+x) - std::log(1.0-x) + exponent;
        }else{
            return inv_sqrt_2pi * 2.0 * std::exp(exponent)/( (1.0-x)*(1.0+x)*sigma );
        }
    }else{
        if(log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}


double densities::evalLogNormal(const double &x, const double &mu, const double &sigma, bool log)
{
    if( (x > 0.0) && (sigma > 0.0)){
        
        double exponent = -.5*(std::log(x)-mu)*(std::log(x)-mu)/(sigma*sigma);
        if(log){
            return -std::log(x) - std::log(sigma) - .5*log_two_pi + exponent;
        }else{
            return inv_sqrt_2pi*std::exp(exponent)/(sigma*x);
        }
    }else{
        if(log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
}


double densities::evalUniform(const double &x, const double &lower, const double &upper, bool log)
{

    if( (x > lower) && (x <= upper)){
        
        double width = upper-lower;
        if(log){
            return -std::log(width);
        }else{
            return 1.0/width;
        }
    }else{
        if(log){
            return -1.0/0.0;
        }else{
            return 0.0;
        }
    }
    
}



////////////////////////////////////////////////
/////////           samplers           /////////
////////////////////////////////////////////////

/////////////// Univariate Normal Sampler


densities::UnivNormSampler::UnivNormSampler()
    : m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())}, 
      m_z_gen(0.0, 1.0)
{
    setMean(0.0);
    setStdDev(1.0);
}


densities::UnivNormSampler::UnivNormSampler(const double &mu, const double &sigma)
    : m_rng{static_cast<std::uint32_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count())}, 
      m_z_gen(0.0, 1.0)
{
    setMean(mu); 
    setStdDev(sigma);
}


void densities::UnivNormSampler::setMean(const double &mu)
{
    m_mu = mu;
}


void densities::UnivNormSampler::setStdDev(const double &sigma)
{
    m_sigma = sigma;
}


double densities::UnivNormSampler::sample()
{
    return m_mu + m_sigma *m_z_gen(m_rng);
}


/////////////// Multivariate Normal Sampler


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
    for (size_t jj=0; jj< m_mean.rows(); ++jj) 
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
