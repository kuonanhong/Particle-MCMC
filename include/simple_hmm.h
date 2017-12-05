#include "fshmm.h"



//!  A simple 2 state HMM. 
/*!
 * Just used for testing purposes. 
*/
class SimpleHmm : public FSHMM
{
private:
    double m_lowVar;
    double m_highVar;
    
public:

    //! The constructor.
    /** 
     * @brief constructs an object.
     * @param lowVar the lower variance for one of the regimes
     * @param highVar the higher variance for the other regime
     * @param initState the initial state distribution.
     * @param transMat the transition matrix. element ij is the prob of trans from i to j.
     */
    SimpleHmm(double lowVar, double highVar, const Vec &initState, const Mat &transMat);
    
    
    //! Convenience function that genrates the vector of conditional densities. 
    /**
     * @brief p(y_t|x_t)
     * @param data the most recent observed data point.
     * @return Vec of all the probabilities conditional on different xts.
     */
    Vec obsDens(const Vec &data);
  
};