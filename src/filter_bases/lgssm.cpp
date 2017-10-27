#include "lgssm.h"
#include <iostream> //for temporary testing

Lgssm::Lgssm() {}


Lgssm::Lgssm(const Vec &initStateMean, const Mat &initStateVar) 
        : m_fresh(true), m_predMean(initStateMean), m_predVar(initStateVar) {}


Lgssm::~Lgssm() {}


void Lgssm::updatePrior(const Mat &stateTransMat, 
                        const Mat &cholStateVar, 
                        const Mat &stateInptAffector, 
                        const Vec &inputData)
{
    Mat Q = cholStateVar.transpose() * cholStateVar;
    m_predMean = stateTransMat * m_filtMean + stateInptAffector * inputData;
    m_predVar  = stateTransMat * m_filtVar * stateTransMat.transpose() + Q;
}


void Lgssm::updatePosterior(const Vec &yt, 
                            const Mat &obsMat, 
                            const Mat &obsInptAffector, 
                            const Vec &inputData, 
                            const Mat &cholObsVar)
{
    Mat R = cholObsVar.transpose() * cholObsVar; //obs
    Mat sigma = obsMat * m_predVar * obsMat.transpose() + R; // pred or APA' + R 
    Mat symSigma = (sigma.transpose() + sigma )/2.0; // ensure symmetric
    Mat siginv = symSigma.inverse();
    Mat K = m_predVar * obsMat.transpose() * siginv;
    Vec obsPred = obsMat * m_predMean + obsInptAffector * inputData;
    Vec innov = yt - obsPred;
    m_filtMean = m_predMean + K*innov;
    m_filtVar  = m_predVar - K * obsMat * m_predVar;

    // conditional likelihood stuff
    Mat quadForm = innov.transpose() * siginv * innov;
    Mat cholSig ( sigma.llt().matrixL() );
    double logDet = 2.0*cholSig.diagonal().array().log().sum();
    m_lastLogCondLike = -.5*innov.rows()*log(2*pi) - .5*logDet - .5*quadForm(0,0);
}


double Lgssm::getLogCondLike() const
{
    return m_lastLogCondLike;
}


Vec Lgssm::getFiltMean() const
{
    return m_filtMean;
}


Mat Lgssm::getFiltVar() const
{
    return m_filtVar;
}

      
void Lgssm::update(const Vec &yt, 
                   const Mat &stateTrans, 
                   const Mat &cholStateVar, 
                   const Mat &stateInptAffector, 
                   const Vec &inData,
                   const Mat &obsMat,
                   const Mat &obsInptAffector, 
                   const Mat &cholObsVar)
{
    // this assumes that we have latent states x_{1:...} and y_{1:...} (NOT x_{0:...})
    // for that reason, we don't have to run updatePrior() on the first iteration
    if (m_fresh == true)
    {
        updatePosterior(yt, obsMat, obsInptAffector, inData, cholObsVar);
        m_fresh = false;
    }else 
    {
        updatePrior(stateTrans, cholStateVar, stateInptAffector, inData);
        updatePosterior(yt, obsMat, obsInptAffector, inData, cholObsVar);
    }
}