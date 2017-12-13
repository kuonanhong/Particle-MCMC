#include "hmm_rbpf.h"

#include <iostream>
// doesn't instantiate or resize LGSSM mods because 
// their first time prior might depend on samples
Hmm_Rbpf::Hmm_Rbpf(unsigned numParts, HMMRBPFResampStyle rt)
    : m_numParts(numParts), m_now(0), m_lastLogCondLike(0.0), m_resampTechnique(rt)
{
    m_logUnNormWeights.resize(numParts);
    m_p_samps.resize(numParts);
    
    // set all the weights to uniform
    std::vector<double>::iterator it;
    for(it = m_logUnNormWeights.begin(); it != m_logUnNormWeights.end(); ++it){
        *it = 0.0;
    }
}


Hmm_Rbpf::~Hmm_Rbpf()
{
}


void Hmm_Rbpf::filter(const Vec &data, 
                      const std::vector<std::function<const Mat(const Vec &x1tProbs, const Vec &x2t)> >& fs)
{

    if( m_now == 0){ // first data point coming
    
        // initialize and update the closed-form mods        
        Vec tmpProbs;
        Mat tmpTransMat;
        double logWeightAdj;
        double tmpForFirstLike(0.0);
        for(unsigned ii = 0; ii < m_numParts; ++ii){
            
            m_p_samps[ii] = q1Samp(data); 
            tmpProbs = initHMMProbVec(m_p_samps[ii]);
            tmpTransMat = initHMMTransMat(m_p_samps[ii]);
            m_p_innerMods.emplace_back(tmpProbs, tmpTransMat); 
            updateFSHMM(m_p_innerMods[ii], data, m_p_samps[ii]);
            logWeightAdj = std::log(m_p_innerMods[ii].getCondLike()) + logMuEv(m_p_samps[ii]) - logQ1Ev(m_p_samps[ii], data); 

            m_logUnNormWeights[ii] += logWeightAdj;
            tmpForFirstLike += std::exp(logWeightAdj);
        }
        m_lastLogCondLike = std::log(tmpForFirstLike) - std::log(m_numParts); // store likelihood

        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        m_dimState = q1Samp(data).rows();
        std::fill(m_expectations.begin(), m_expectations.end(), Mat::Zero(m_dimState, m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double m = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        for(auto & h : fs){
            
            int rows = h(m_p_innerMods[0].getFilterVec(), m_p_samps[0]).rows();
            int cols = h(m_p_innerMods[0].getFilterVec(), m_p_samps[0]).cols();
            Mat numer = Mat::Zero(rows,cols);
            Mat ones = Mat::Ones(rows,cols);
            double denom(0.0);
            Mat tmp;
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ 
                tmp = h(m_p_innerMods[prtcl].getFilterVec(), m_p_samps[prtcl]);
                tmp = tmp.array().log().matrix() + (m_logUnNormWeights[prtcl] - m)*ones;
                numer = numer + tmp.array().exp().matrix();
                denom += std::exp( m_logUnNormWeights[prtcl] - m );
            }
            m_expectations[fId] = numer/denom;
            fId++;
        }

        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == HMMRBPFResampStyle::everytime_multinomial)
            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
            
        // advance time step
        m_now ++;
    }
    else { //m_now > 0
        
        // update
        Vec newX2Samp;
        double logUnNormWeightUpdate;
        double tmpLikeNumer(0.0);
        double tmpLikeDenom(0.0);
        for(unsigned ii = 0; ii < m_numParts; ++ii){
            
            newX2Samp = qSamp(m_p_samps[ii], data);
            updateFSHMM(m_p_innerMods[ii], data, newX2Samp);
            logUnNormWeightUpdate = std::log(m_p_innerMods[ii].getCondLike())
                                    + logFEv(newX2Samp, m_p_samps[ii]) 
                                    - logQEv(newX2Samp, m_p_samps[ii], data);
            
            tmpLikeDenom += std::exp(m_logUnNormWeights[ii]);
            m_logUnNormWeights[ii] += logUnNormWeightUpdate;
            tmpLikeNumer += std::exp(m_logUnNormWeights[ii]); 
            m_p_samps[ii] = newX2Samp;
        }
        m_lastLogCondLike = std::log(tmpLikeNumer) - std::log( tmpLikeDenom );
        
        // calculate expectations before you resample
        std::fill(m_expectations.begin(), m_expectations.end(), Mat::Zero(m_dimState, m_dimState)); // TODO: should this be Mat::Zero(m_dimState, m_dimState)?
        int fId(0);
        double m = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        for(auto & h : fs){
            int rows = h(m_p_innerMods[0].getFilterVec(), m_p_samps[0]).rows();
            int cols = h(m_p_innerMods[0].getFilterVec(), m_p_samps[0]).cols();
            Mat numer = Mat::Zero(rows,cols);
            Mat ones = Mat::Ones(rows,cols);
            Mat tmp;
            double denom(0.0);
            for(size_t prtcl = 0; prtcl < m_numParts; ++prtcl){ 
                tmp = h(m_p_innerMods[prtcl].getFilterVec(), m_p_samps[prtcl]);
                tmp = tmp.array().log().matrix() + (m_logUnNormWeights[prtcl] - m)*ones;
                numer = numer + tmp.array().exp().matrix();
                denom += std::exp( m_logUnNormWeights[prtcl] - m );
            }
            m_expectations[fId] = numer/denom;
            fId++;
        }
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == HMMRBPFResampStyle::everytime_multinomial)
            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
//        else if ( (m_resampTechnique == HMMRBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
//            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
        
        // update time step
        m_now ++;
    }
    
}


double Hmm_Rbpf::getLogCondLike() const
{
    return m_lastLogCondLike;
}


std::vector<Mat> Hmm_Rbpf::getExpectations() const
{
    return m_expectations;
}


void Hmm_Rbpf::resampMultinomHRBPF(std::vector<FSHMM> &oldMods, 
                                    std::vector<Vec> &oldSamps, 
                                    std::vector<double> &oldLogWts)
{
    m_resampler.ressampHRBPF(oldMods, oldSamps, oldLogWts);
}