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

    //! Default constructor (need this for resampling function)
    /**
     * @brief Default constructor (need this for resampling function)
     */
    FSHMM();

    //! Constructor
    /*!
      \param initStateDistr first time state prior distribution.
      \param transMat time homogeneous transition matrix.
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
    
    //! Perform a HMM filter update.
    /**
     * @brief Perform a HMM filter update.
     * @param yt the current datum.
     * @param condDensVec the vector (in x_t) of p(y_t|x_t)
     */
    void update(const Vec &yt, const Vec &condDensVec);

};

#endif // FSHMM_H
