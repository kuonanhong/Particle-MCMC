#include "kalman_rbpf.h"

// constructors
// doesn't instantiate or resize LGSSM mods because 
// their first time prior might depend on samples
Kalman_RBPF::Kalman_RBPF(int numParts, RBPFResampStyle resampTechnique)
    : m_lastCondLike(1.0), m_numParts(numParts), m_resampTechnique(resampTechnique), m_now(0)
{
    // resize everything besides the k-filter mods
    m_unNormWeights.resize(m_numParts);
    m_p_samps.resize(m_numParts);
    
    //set all the weights to uniform
    std::vector<double>::iterator it;
    for(it = m_unNormWeights.begin(); it != m_unNormWeights.end(); ++it)
    {
        *it = 1.0;
    }
}

// destructor
Kalman_RBPF::~Kalman_RBPF(){}

// workhorse function
void Kalman_RBPF::filter(const Vec &data)
{

    if( m_now == 0){ // first data point coming
    
        // initialize and update the closed-form mods        
        Vec tmpMean;
        Mat tmpVar;
        double tmpInnerCondLike; // p(y_1|x_{2,1})
        double weightAdj;
        double tmpForFirstLike(0.0);
        for(unsigned int ii = 0; ii < m_numParts; ++ii){
            m_p_samps[ii] = q1Samp(data); 
            tmpMean = initKalmanMean(m_p_samps[ii]);
            tmpVar  = initKalmanVar(m_p_samps[ii]);
            m_p_innerMods.emplace_back(tmpMean, tmpVar); // time 1 prior
            updateKalman(m_p_innerMods[ii], data, m_p_samps[ii]);
            tmpInnerCondLike = exp( m_p_innerMods[ii].getLogCondLike() );
            weightAdj = ( tmpInnerCondLike * muEv(m_p_samps[ii]) ) / q1Ev(m_p_samps[ii], data); 

            // can't divide by 0
            if( std::isnan(weightAdj)){
                throw std::runtime_error("divide by 0 error: q1Ev evaluated to 0!");
            }
     
            m_unNormWeights[ii] *= weightAdj;
            tmpForFirstLike += weightAdj;
        }
        m_lastCondLike = tmpForFirstLike / m_numParts; // store likelihood
        
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == RBPFResampStyle::everytime_multinomial)
            ressampMultinomKRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
        else if ( (m_resampTechnique == RBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
            ressampMultinomKRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
            
        // advance time step
        m_now ++;
    }
    else { //m_now > 0
        
        // update
        Vec newX2Samp;
        double tmpInnerCondLike;
        double unNormWeightUpdate;
        double tmpLikeNumer(0.0);
        double tmpLikeDenom(0.0);
        for(unsigned int ii = 0; ii < m_numParts; ++ii){
            newX2Samp = qSamp(m_p_samps[ii], data);
            updateKalman(m_p_innerMods[ii], data, newX2Samp);
            tmpInnerCondLike = exp( m_p_innerMods[ii].getLogCondLike() );
            unNormWeightUpdate = (tmpInnerCondLike * fEv(newX2Samp, m_p_samps[ii])) / qEv(newX2Samp, m_p_samps[ii], data);
            
            // can't divide by 0
            if( std::isnan(unNormWeightUpdate)){
                throw std::runtime_error("divide by 0 error: q1Ev evaluated to 0!");
            }
            
            tmpLikeDenom += m_unNormWeights[ii];
            m_unNormWeights[ii] *= unNormWeightUpdate;
            tmpLikeNumer += m_unNormWeights[ii]; 
            m_p_samps[ii] = newX2Samp;
        }
        m_lastCondLike = tmpLikeNumer / tmpLikeDenom;
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == RBPFResampStyle::everytime_multinomial)
            ressampMultinomKRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
        else if ( (m_resampTechnique == RBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
            ressampMultinomKRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
        
        // update time step
        m_now ++;
    }
    
}
 

double Kalman_RBPF::getCondLike() const
{
    return m_lastCondLike;
}

std::vector<double> Kalman_RBPF::getWeights() const
{
    return m_unNormWeights;
}

void Kalman_RBPF::ressampMultinomKRBPF(std::vector<Lgssm> &oldMods, 
                                       std::vector<Vec> &oldSamps, 
                                       std::vector<double> &oldWts)
{
    m_resampler.ressampKRBPF(oldMods, oldSamps, oldWts);
}