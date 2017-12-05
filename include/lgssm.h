#ifndef LGSSM_H
#define LGSSM_H

#include <Eigen/Dense> //linear algebra stuff
#include <math.h>       /* log */

// for p(y_t|y_{1:t-1})
const double pi = 3.14159265358979323846;

// shorthand names
typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;

//!  A base-class for Kalman Filtering. 
/*!
 * Inherit from this if your model is a linear-gaussian state space model. 
*/
class Lgssm
{
private: 
    Vec m_predMean;
    Vec m_filtMean;
    Mat m_predVar;
    Mat m_filtVar;
    double m_lastLogCondLike; 
    bool m_fresh;
    
    //TODO: handle diagonal variance matrices, and ensure symmetricness in other ways
    void updatePrior(const Mat &stateTransMat, 
                     const Mat &cholStateVar, 
                     const Mat &stateInptAffector, 
                     const Vec &inputData);
    void updatePosterior(const Vec &yt, 
                         const Mat &obsMat, 
                         const Mat &obsInptAffector, 
                         const Vec &inputData, 
                         const Mat &cholObsVar);
public:
    //! Default Constructor
    Lgssm();

    //! Constructor
    /*!
      \param initStateMean first time state prior distribution's mean vector.
      \param initStateVar first time state prior distribution's covariance matrix.
    */
    Lgssm(const Vec &initStateMean, const Mat &initStateVar);
    
    //! Destructor.
    ~Lgssm();

    //! Get the latest conditional likelihood.
    /**
     * \return the latest conditional likelihood.
     */  
    double getLogCondLike() const;

    //! Get the latest filtering mean.
    /**
     * \return the mean of the filtering distribution.
     */  
    Vec getFiltMean() const;
    
    //! Get the latest filtering covariance matrix.
    /**
     * \return the covariance matrix of the filtering distribution.
     */      
    Mat getFiltVar() const;
    
    //! Perform a Kalman filter predict-and-update.
    /**
     * \param yt the new data point.
     * \param stateTrans the transition matrix of the state
     * \param cholStateVar the Cholesky Decomposition of the state noise covariance matrix.
     * \param stateInptAffector the matrix affecting how input data affects state transition.
     * \param inputData exogenous input data
     * \param obsMat the observation/emission matrix of the observation's conditional (on the state) distn.
     * \param obsInptAffector the matrix affecting how input data affects the observational distribution.
     * \param cholObsVar the Cholesky Decomposition of the observatio noise covariance matrix.
     */      
    void update(const Vec &yt, 
                const Mat &stateTrans, 
                const Mat &cholStateVar, 
                const Mat &stateInptAffector, 
                const Vec &inputData,
                const Mat &obsMat,
                const Mat &obsInptAffector, 
                const Mat &cholObsVar);

};

#endif // LGSSM_H
