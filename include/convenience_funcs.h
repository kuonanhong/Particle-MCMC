#ifndef CONVENIENCE_FUNCS_H
#define CONVENIENCE_FUNCS_H

#include <fstream>
#include <vector>
#include <Eigen/Dense>

typedef Eigen::Matrix< double, Eigen::Dynamic, 1> Vec;

namespace convenience_funcs{

    /**
     * @brief writes a Eigen::VectorXd to a row of an ostream.
     * @param vec the Eigen::VectorXd.
     * @param ofs the target ofstream.
     */
    void logParams(const Vec &vec, std::ofstream &ofs);

    
    /**
     * @brief writes a std::vector<double> to a row of an ofstream.
     * @param vec the std::vector<double> to be written.
     * @param ofs the target ofstream.
     */
    void logParams(const std::vector<double> &vec, std::ofstream &ofs);
    
    /**
     * @brief writes a std::vector<Eigen::VectorXd> to a row of an ofstream.
     * @param vec the std::vector<Eigen::VectorXd> to be written.
     * @param ofs the target ofstream.
     */
    void logParams(const std::vector<Vec> &vecOfParams, std::ofstream &ofs);
    
    
    /**
     * @brief writes a std::array<Vec, N> to a row of an ofstream
     * @param arrOfParams an array of Vecs of parameters to write
     * @param ofs the target ofstream
     */
    template <size_t N>
    void logParams(const std::array<Vec, N> &arrOfParams, std::ofstream &ofs);

    
    /**
     * @brief reads in data in a csv file with no header and separates by commas. 
     * @param fileLoc the string filepath of the file.
     * @return a std::vector of your data. Each elemenet is a row in Eigen::VectorXd form.
     */
    std::vector<Vec> readInData(const std::string& fileLoc);


    
} // namespace conv_funcs

template <size_t N>
void convenience_funcs::logParams(const std::array<Vec, N> &arrOfParams, std::ofstream &ofs)
{
    for(size_t i = 0; i < arrOfParams.size(); ++i){ // every Vec 
        for(size_t j = 0; j < arrOfParams[i].rows(); ++j){  // every elmeent of Vec
            if ( i == 0 && j == 0) // first thing
                ofs << arrOfParams[i](j);
            else
                ofs << ", " << arrOfParams[i](j);
        }
    }
    ofs << "\n";
}

#endif //CONVENIENCE_FUNCS_H