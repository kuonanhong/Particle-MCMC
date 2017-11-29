#include "simple_hmm.h"

#include "densities.h"


SimpleHmm::SimpleHmm(double lowVar, double highVar, const Vec &initState, const Mat &transMat) :
    m_lowVar(lowVar), m_highVar(highVar), FSHMM(initState, transMat)
{
}
    

Vec SimpleHmm::obsDens(const Vec &data)
{
    unsigned dim = this->dimState();
    Vec gVec(dim);
    Mat stateOneCov(1,1);
    stateOneCov(0,0) = m_lowVar;
    Mat stateTwoCov(1,1);
    stateTwoCov(0,0) = m_highVar;
    gVec(0) = densities::evalMultivNorm(data, Vec::Zero(1), stateOneCov, false);
    gVec(1) = densities::evalMultivNorm(data, Vec::Zero(1), stateTwoCov, false);
    
    return gVec;
}
