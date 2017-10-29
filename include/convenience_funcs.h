#ifndef CONVENIENCE_FUNCS_H
#define CONVENIENCE_FUNCS_H

#include <fstream>
#include <vector>
#include <Eigen/Dense>

typedef Eigen::Matrix< double, Eigen::Dynamic, 1> Vec;

namespace convenience_funcs{

    
    /**
     * @brief writes a vec to a row of an ofstream.
     * @param vec the std::vector to be written.
     * @param ofs the target ofstream.
     */
    void logParams(const std::vector<double> &vec, std::ofstream &ofs);
    
    
    /**
     * @brief reads in data in a csv file with no header and separates by commas. 
     * @param fileLoc the string filepath of the file.
     * @param numCols the number of columns of data you're expecting.
     */
    std::vector<Vec> readInData(const std::string& fileLoc, unsigned int numCols);


    
} // namespace conv_funcs

#endif //CONVENIENCE_FUNCS_H