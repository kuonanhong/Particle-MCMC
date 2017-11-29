#ifndef FSHMM_H
#define FSHMM_H

#include <Eigen/Dense> //linear algebra stuff
#include <cmath>       /* log */

// shorthand names
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

//!  A base-class for Finite State HMM Filtering. 
/*!
 * Inherit from this if your model is a finite state hidden Markov model. 
*/
class FSHMM
{
private:
    Vec m_filtVec;
    Mat m_transMatTranspose;
    double m_lastCondLike; 
    bool m_fresh;
public:

    //! Constructor
    /*!
      \param initStateMean first time state prior distribution's mean vector.
      \param initStateVar first time state prior distribution's covariance matrix.
    */
    FSHMM(const Vec &initStateDistr, const Mat &transMat);
    
    //! Destructor.
    ~FSHMM();

    //! Get the latest conditional likelihood.
    /**
     * \return the latest conditional likelihood.
     */  
    double getCondLike() const;
    
    //! Get the current filter vector.
    /**
     * @brief get the current filter vector.
     * @return a probability Eigen::VectorXd
     */
    Vec getFilterVec() const;
    
    //! Get the dimension of the state
    /**
     * @brief Returns the dimension of the state process.
     */
    unsigned dimState() const;
    
    //! Perform a Kalman filter predict-and-update.
    /**
     * \param yt the new data point.
     */      
    void update(const Vec &yt);
    
    //! Returns a probability p(data|x) for each x in a row vector
    /**
     * @param data the observed datum at a given time
     * @return a row vector of doubles corresponding to each state value
     */
    virtual Vec obsDens(const Vec &data) = 0;

};

#endif // FSHMM_H
