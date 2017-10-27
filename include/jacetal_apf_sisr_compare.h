#ifndef JACETAL_APF_SISR_COMPARE_H
#define JACETAL_APF_SISR_COMPARE_H

#include <vector>
#include <Eigen/Dense>
#include <iostream>
#include <fstream> //ifstream, ofstream

typedef Eigen::Matrix< double, Eigen::Dynamic, 1              > Vec;
typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > Mat;


std::vector<Vec> read_data(const std::string &dataFile, int num_obs, int size_obs);


void unFlattenParams(Mat &beta, Vec &phi, Vec &mu, Vec &sigma, Vec &RstdDevs, const std::vector<double> &flatOnes);


void jacetal_apf_sisr_compare();


#endif //JACETAL_APF_SISR_COMPARE_H