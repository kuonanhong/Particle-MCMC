#include "fshmm.h"


FSHMM::FSHMM()
{
}


FSHMM::FSHMM(const Vec &initStateDistr, const Mat &transMat) :
    m_filtVec(initStateDistr), m_transMatTranspose(transMat.transpose()), m_lastCondLike(0.0), m_fresh(false)
{
}


FSHMM::~FSHMM()
{
}


double FSHMM::getCondLike() const
{
    return m_lastCondLike;
}


Vec FSHMM::getFilterVec() const
{
    return m_filtVec;
}

unsigned FSHMM::dimState() const
{
    return m_filtVec.rows(); // this is a column vector
}


void FSHMM::update(const Vec &yt, const Vec &condDensVec)
{
    if (!m_fresh)  // hasn't seen data before and so filtVec is just time 1 state prior
    {
        m_filtVec = m_filtVec.cwiseProduct( condDensVec ); // now it's p(x_1, y_1)
        m_lastCondLike = m_filtVec.sum();
        m_filtVec /= m_lastCondLike;
        m_fresh = true;
        
    } else { // has seen data before
        m_filtVec = m_transMatTranspose * m_filtVec; // now p(x_t |y_{1:t-1})
        m_filtVec = m_filtVec.cwiseProduct( condDensVec ); // now p(y_t,x_t|y_{1:t-1})
        m_lastCondLike = m_filtVec.sum();
        m_filtVec /= m_lastCondLike; // now p(x_t|y_{1:t})
    }
}