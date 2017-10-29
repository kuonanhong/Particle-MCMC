#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H


namespace transformations
{

    /**
     * @brief Maps -1 < phi < 1 to the reals.
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

    
}


#endif //TRANSFORMATIONS_H