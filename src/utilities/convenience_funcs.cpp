#include "convenience_funcs.h"

#include <iostream> //cerr


void convenience_funcs::logParams(const Vec &vec, std::ofstream &ofs)
{
    for(size_t i = 0; i < vec.rows(); ++i){
        if( i == 0){
            ofs << vec(i);
        } else {
            ofs << "," << vec(i);                                        
        }
    }
    ofs << "\n";
}

void convenience_funcs::logParams(const std::vector<double> &vec, std::ofstream &ofs)
{
    for(size_t i = 0; i < vec.size(); ++i){
        if( i == 0){
            ofs << vec[i];
        } else {
            ofs << "," << vec[i];                                        
        }
    }
    ofs << "\n";
}


void convenience_funcs::logParams(const std::vector<Vec> &vecOfParams, std::ofstream &ofs)
{
    for(size_t i = 0; i < vecOfParams.size(); ++i){ // every Vec 
        for(size_t j = 0; j < vecOfParams[i].rows(); ++j){  // every elmeent of Vec
            if ( i == 0 && j == 0) // first thing
                ofs << vecOfParams[i](j);
            else
                ofs << ", " << vecOfParams[i](j);
        }
    }
    ofs << "\n";
}


std::vector<Vec> convenience_funcs::readInData(const std::string& fileLoc)
{
    // return this at end
    std::vector<Vec> data;
    
    std::string line;
    std::ifstream ifs(fileLoc);
    std::string one_number;
    
    // check if we can open inFile
    if(!ifs.is_open()){
        std::cerr << "error: readInData() could not read data from " << fileLoc << std::endl;
    }
    
    while ( std::getline(ifs, line) ){     // get a whole row as a string
    
        std::vector<double> data_row;
        try{
            // get a single element on a row
            std::istringstream buff(line);
            
            // make one number between commas
            while(std::getline(buff, one_number, ',')){ 
                data_row.push_back(std::stod(one_number));
            }
            
        } catch (const std::invalid_argument& ia){
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }   
        
        // now append this Vec to your collection
        Eigen::Map<Eigen::VectorXd> drw(&data_row[0], data_row.size());
        data.push_back(drw);
    }
    
    return data;
}