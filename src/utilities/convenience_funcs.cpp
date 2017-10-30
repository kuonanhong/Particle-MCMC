#include "convenience_funcs.h"

#include <iostream> //cerr


void convenience_funcs::logParams(const Vec &vec, std::ofstream &ofs)
{
    for(unsigned int i = 0; i < vec.rows(); ++i){
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
    for(unsigned int i = 0; i < vec.size(); ++i){
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
    for(unsigned i = 0; i < vecOfParams.size(); ++i){ // every Vec 
        for(unsigned j = 0; j < vecOfParams[i].rows(); ++j){  // every elmeent of Vec
            if ( i == 0 && j == 0) // first thing
                ofs << vecOfParams[i](j);
            else
                ofs << ", " << vecOfParams[i](j);
        }
    }
    ofs << "\n";
}

std::vector<Vec> convenience_funcs::readInData(const std::string& fileLoc, unsigned int numCols)
{
    // return this at end
    std::vector<Vec> data;
    
    std::string line;
    std::ifstream inFile(fileLoc);
    std::string oneNumberOnOneLine;
    unsigned int elemNum;  // temporary variable numbering the column of data we'r ereading 
    
    // iterate over every line
    while ( std::getline(inFile, line) ){ 
    
        Vec dataRow(numCols);
        elemNum = 0;
        try{
            std::istringstream buff(line);
            while(std::getline(buff, oneNumberOnOneLine, ',')){ // read the line until you get a ','
                dataRow(elemNum) = std::stod(oneNumberOnOneLine);
                elemNum++;
            }
        } catch (const std::invalid_argument& ia){
            std::cerr << "Invalid Argument: " << ia.what() << "\n";
            continue;
        }   
        
        // now append this Vec to your collection
        data.push_back(dataRow);
    }
    
    return data;
}