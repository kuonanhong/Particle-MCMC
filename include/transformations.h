#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H


namespace transformations
{

    /**
     * @brief Maps (-1, 1) to the reals.
     * @param phi
     * @return psi
     */
    double twiceFisher(const double &phi);


    /**
     * @brief Maps a real number to the itnerval (-1,1).
     * @param psi
     * @return phi
    */
    double invTwiceFisher(const double &psi);
    
    
    /**
     * @brief Maps (0,1) to the reals.
     * @param p
     * @return logit(p)
    */
    double logit(const double &p);
    
    
    /**
     * @brief Maps the reals to (0,1)
     * @param r
     * @return p = invlogit(p)
    */
    double inv_logit(const double &r);

    
}


#endif //TRANSFORMATIONS_H