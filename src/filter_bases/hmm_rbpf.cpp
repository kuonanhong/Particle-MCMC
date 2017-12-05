#include "hmm_rbpf.h"


// doesn't instantiate or resize LGSSM mods because 
// their first time prior might depend on samples
Hmm_Rbpf::Hmm_Rbpf(unsigned numParts, HMMRBPFResampStyle rt)
    : m_numParts(numParts), m_now(0), m_lastCondLike(1.0), m_resampTechnique(rt)
{
    m_unNormWeights.resize(numParts);
    m_p_samps.resize(numParts);
    
    // set all the weights to uniform
    std::vector<double>::iterator it;
    for(it = m_unNormWeights.begin(); it != m_unNormWeights.end(); ++it){
        *it = 1.0;
    }
}


Hmm_Rbpf::~Hmm_Rbpf()
{
}


void Hmm_Rbpf::filter(const Vec &data)
{

    if( m_now == 0){ // first data point coming
    
        // initialize and update the closed-form mods        
        Vec tmpProbs;
        Mat tmpTransMat;
        double weightAdj;
        double tmpForFirstLike(0.0);
        for(unsigned ii = 0; ii < m_numParts; ++ii){
            
            m_p_samps[ii] = q1Samp(data); 
            tmpProbs = initHMMProbVec(m_p_samps[ii]);
            tmpTransMat = initHMMTransMat(m_p_samps[ii]);
            m_p_innerMods.emplace_back(tmpProbs, tmpTransMat); 
            updateFSHMM(m_p_innerMods[ii], data, m_p_samps[ii]);
            weightAdj = m_p_innerMods[ii].getCondLike() * muEv(m_p_samps[ii]) / q1Ev(m_p_samps[ii], data); 

            // can't divide by 0
            if( std::isnan(weightAdj)){
                throw std::runtime_error("divide by 0 error: q1Ev evaluated to 0!");
            }
     
            m_unNormWeights[ii] *= weightAdj;
            tmpForFirstLike += weightAdj;
        }
        m_lastCondLike = tmpForFirstLike / m_numParts; // store likelihood
        
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == HMMRBPFResampStyle::everytime_multinomial)
            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
//        else if ( (m_resampTechnique == HMMRBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
//            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
            
        // advance time step
        m_now ++;
    }
    else { //m_now > 0
        
        // update
        Vec newX2Samp;
        double unNormWeightUpdate;
        double tmpLikeNumer(0.0);
        double tmpLikeDenom(0.0);
        for(unsigned ii = 0; ii < m_numParts; ++ii){
            
            newX2Samp = qSamp(m_p_samps[ii], data);
            updateFSHMM(m_p_innerMods[ii], data, newX2Samp);
            unNormWeightUpdate = m_p_innerMods[ii].getCondLike() * fEv(newX2Samp, m_p_samps[ii]) / qEv(newX2Samp, m_p_samps[ii], data);
            
            // can't divide by 0
            if( std::isnan(unNormWeightUpdate)){
                throw std::runtime_error("divide by 0 error: qEv evaluated to 0!");
            }
            
            tmpLikeDenom += m_unNormWeights[ii];
            m_unNormWeights[ii] *= unNormWeightUpdate;
            tmpLikeNumer += m_unNormWeights[ii]; 
            m_p_samps[ii] = newX2Samp;
        }
        m_lastCondLike = tmpLikeNumer / tmpLikeDenom;
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == HMMRBPFResampStyle::everytime_multinomial)
            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
//        else if ( (m_resampTechnique == HMMRBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
//            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
        
        // update time step
        m_now ++;
    }
    
}


double Hmm_Rbpf::getCondLike() const
{
    return m_lastCondLike;
}


void Hmm_Rbpf::resampMultinomHRBPF(std::vector<FSHMM> &oldMods, 
                                    std::vector<Vec> &oldSamps, 
                                    std::vector<double> &oldWts)
{
    m_resampler.ressampHRBPF(oldMods, oldSamps, oldWts);
}