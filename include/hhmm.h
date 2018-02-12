#ifndef HHMM_H
#define HHMM_H

#include <Eigen/Dense> //linear algebra stuff
#include <cmath>       /* log */

// shorthand names
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

//!  A base-class for a 2 level Hierarchical Hidden Markov Model
/*!
 * Inherit from this if your model is a 2 level Hierarchical Hidden Markov Model.
 * Note that this does not store the filtering distribution of one big container of probabilities 
 * (one for each configuration). Rather, it has one filtering distribution for the first autonomous
 * component p(x_t^1), then it has an array of vectors for the second p(x_t^2 | x_t^1, y_{1:t}). 
 * These are the production states' distributions conditioning on each value of the first component.
*/
template<size_t n1, size_t n2>
class HHMM
{

private:
    Eigen::Matrix<double, n1, 1> m_filtVecLevel1;
    std::array<Eigen::Matrix<double, n2, 1>, n1> m_filtVecsLevel2;
    Eigen::Matrix<double, n1, n1> m_transMatTransposeLevel1;
    std::array<Eigen::Matrix<double, n2, n2>, n1> m_transMatsTransposeLevel2;
    double m_lastCondLike; 
    bool m_fresh;

public:


    //! Default constructor (need this for resampling function)
    /**
     * @brief Default constructor (need this for resampling function)
     */
    HHMM();


    //! Constructor
    /*!
      \param initStateDistr first time state prior distribution.
      \param transMat time homogeneous transition matrix.
    */
    HHMM(const Eigen::Matrix<double,n1 ,1> &initStateDistrLevel1, 
         const std::array<Eigen::Matrix<double, n2, 1>, n1> &initStateDistrLevel2, 
         const Eigen::Matrix<double,n1 ,n1> &transMat1, 
         const std::array<Eigen::Matrix<double,n2 ,n2>, n1> transMats2);
    
    
    //! Destructor.
    ~HHMM();

    //! Get the latest conditional likelihood.
    /**
     * \return the latest conditional likelihood.
     */  
    double getCondLike() const;
    
  
    //! Perform a HMM filter update.
    /**
     * @brief Perform a HMM filter update.
     * @param yt the current datum.
     * @param condDensVec the vector (in x_t) of p(y_t| time t production states)
     */
    void update(const Vec &yt, const Eigen::Matrix<double, n2, 1> &condDensVec);

};



//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////// implementations /////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////

template<size_t n1, size_t n2>
HHMM<n1,n2>::HHMM()
{
}


template<size_t n1, size_t n2>
HHMM<n1,n2>::HHMM(const Eigen::Matrix<double,n1 ,1> &initStateDistrLevel1, 
                  const std::array<Eigen::Matrix<double, n2, 1>, n1> &initStateDistrLevel2, 
                  const Eigen::Matrix<double,n1 ,n1> &transMat1, 
                  const std::array<Eigen::Matrix<double,n2 ,n2>, n1> transMats2) :
    m_filtVecLevel1(initStateDistrLevel1),
    m_filtVecsLevel2(initStateDistrLevel2),
    m_transMatTransposeLevel1(transMat1.transpose()), 
    m_lastCondLike(0.0), 
    m_fresh(false)
{
    for(size_t i = 0; i < n1; ++i){
        m_transMatsTransposeLevel2[i] = transMats2[i].transpose();
    }
}
    

template<size_t n1, size_t n2>
HHMM<n1,n2>::~HHMM()
{
}


template<size_t n1, size_t n2>
double HHMM<n1,n2>::getCondLike() const
{
    return m_lastCondLike;
}


template<size_t n1, size_t n2>
void HHMM<n1,n2>::update(const Vec &yt, const Eigen::Matrix<double, n2, 1> &condDensVec)
{
    // see pdf for formulas
    if (!m_fresh)  
    {
        for(size_t i = 0; i < n1; ++i){
            m_filtVecsLevel2[i] = m_filtVecsLevel2[i].cwiseProduct(condDensVec); // now p(y_1, x^2 | x^1) for each x^1
            m_filtVecLevel1(i) = ( m_filtVecsLevel2[i] * m_filtVecLevel1(i) ).sum();  // now each p(y_1, x^1)
            m_filtVecsLevel2[i] = m_filtVecsLevel2[i] / m_filtVecsLevel2[i].sum(); // p(x^2 | x^1, y_1) because we divide by p(y_1 | x^1)
        }
        m_lastCondLike = m_filtVecLevel1.sum();
        m_filtVecLevel1 = m_filtVecLevel1 / m_lastCondLike;
        m_fresh = true;
        
    } else { // has seen data before
        
        m_filtVecLevel1 = m_transMatTransposeLevel1 * m_filtVecLevel1;  // now p(x_t^1 | y_{1:t-1})
        for(size_t i = 0; i < n1; ++i){
            m_filtVecsLevel2[i] = m_transMatsTransposeLevel2[i] * m_filtVecsLevel2[i]; // now p(x_t^2 | x_t^1, y_{1:t-1})
            m_filtVecsLevel2[i] = m_filtVecsLevel2[i].cwiseProduct( condDensVec ); // now p(y_t, x_t^2 | x_t^1, y_{1:t-1})
            m_filtVecLevel1(i) = ( m_filtVecsLevel2[i] * m_filtVecLevel1(i) ).sum(); // each p(y_t, x_t^1 | y_{1:t-1})
            m_filtVecsLevel2[i] = m_filtVecsLevel2[i] / m_filtVecsLevel2[i].sum(); // each p(x_t^2 | x_t^1, y_{1:t})
        }
        m_lastCondLike = m_filtVecLevel1.sum(); // p(y_t | y_{1:t-1})
        m_filtVecLevel1 = m_filtVecLevel1 / m_lastCondLike;
    }
}



#endif // HHMM_H
