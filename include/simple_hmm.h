#include "fshmm.h"

class SimpleHmm : public FSHMM
{
private:
    double m_lowVar;
    double m_highVar;
    
public:

    SimpleHmm(double lowVar, double highVar, const Vec &initState, const Mat &transMat);
    
    Vec obsDens(const Vec &data);
  
};