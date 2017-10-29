#ifndef UTILITIES_H
#define UTILITIES_H

#include <fstream>
#include <vector>

namespace utilities{

void logParams(const std::vector<double> &vec, std::ofstream &ofs);

    
} // namespace utilities

#endif //UTILITIES_H