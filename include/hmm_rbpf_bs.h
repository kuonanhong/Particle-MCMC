#ifndef HMM_RBPF_BS_H
#define HMM_RBPF_BS_H

#include <Eigen/Dense>
#include <vector>
#include <functional> // std::function


#include "fshmm.h"
#include "mn_resampler.h" 
#include "densities.h" // for typedefs

/** enum class for the type of resampling to be performed */
enum class HMMRBPFBSResampStyle {everytime_multinomial, never, ess_multinomial};

//! A base-class for HMM Rao-Blackwellized Particle Filtering. 
/*!
 * @class Hmm_Rbpf_BS
 * @author taylor
 * @file hmm_rbpf_bs.h
 * @brief Rao-Blackwellized particle filtering where the innder models are discrete state hmms.
 * @tparam np the number of parameters
 */
template <size_t np>
class Hmm_Rbpf_BS{
private:
    std::array<FSHMM,np> m_p_innerMods;
    std::array<Vec,np> m_p_samps;
    std::array<double,np> m_logUnNormWeights;
 
    unsigned int m_now;
    unsigned int m_dimState;
    double m_logLastCondLike;
    HMMRBPFBSResampStyle  m_resampTechnique;
    MNResamp<np> m_resampler;
    std::vector<Mat> m_expectations;
    
    void resampMultinomHRBPF(std::array<FSHMM, np> &oldMods, std::array<Vec,np> &oldSamps, std::array<double,np> &oldLogWts);
public:

    //! The constructor.
    /**
     * @brief constructor.
     * @param rt the type of resampling you want to do.
     */
    Hmm_Rbpf_BS(HMMRBPFBSResampStyle rt = HMMRBPFBSResampStyle::everytime_multinomial);
    
    
    //! Destructor.
    /**
     * @brief destructor.
     */
    ~Hmm_Rbpf_BS();


    //! Filter.
    /**
     * @brief filters everything based on a new data point.
     * @param data the most recent time series observation.
     * @param fs a vector of functions computing logE[h(x_1t, x_2t^i)| x_2t^i,y_1:t]. will access the probability vector of x_1t
     */
    void filter(const Vec &data,
                const std::vector<std::function<const Mat(const Vec &x1tProbs, const Vec &x2t)> >& fs 
                    = std::vector<std::function<const Mat(const Vec&, const Vec&)> >());//, const std::vector<std::function<const Mat(const Vec&)> >& fs);


    //! Get the latest conditional likelihood.
    /**
     * @brief Get the latest conditional likelihood.
     * @return the latest conditional likelihood.
     */
    double getLogCondLike() const;
    
    //!
    /**
     * @brief Get vector of expectations.
     * @return vector of expectations
     */
    std::vector<Mat> getExpectations() const;   
    
    //! Sample from the first sampler.
    /**
     * @brief samples the second component of the state at time 1.
     * @param y1 most recent datum.
     * @return a Vec sample for x21.
     */
    virtual Vec muSamp(const Vec &y1) = 0;
    
    
    //! Provides the initial mean vector for each HMM filter object.
    /**
     * @brief provides the initial probability vector for each HMM filter object.
     * @param x21 the second state componenent at time 1.
     * @return a Vec representing the probability of each state element.
     */
    virtual Vec initHMMProbVec(const Vec &x21) = 0;
    
    
    //! Provides the transition matrix for each HMM filter object.
    /**
     * @brief provides the transition matrix for each HMM filter object.
     * @param x21 the second state component at time 1. 
     * @return a transition matrix where element (ij) is the probability of transitioning from state i to state j.
     */
    virtual Mat initHMMTransMat(const Vec &x21) = 0;

    //! Samples the time t second component. 
    /**
     * @brief Samples the time t second component.
     * @param x2tm1 the previous time's second state component.
     * @param yt the current observation.
     * @return a Vec sample of the second state component at the current time.
     */
    virtual Vec fSamp(const Vec &x2tm1, const Vec &yt) = 0;
    
    
    //! How to update your inner HMM filter object at each time.
    /**
     * @brief How to update your inner HMM filter object at each time.
     * @param aModel a HMM filter object describing the conditional closed-form model.
     * @param yt the current time series observation.
     * @param x2t the current second state component.
     */
    virtual void updateFSHMM(FSHMM &aModel, const Vec &yt, const Vec &x2t) = 0;
};

////////////////////////////////////////////////////////////////////
////////////////////////////////// implementations /////////////////
////////////////////////////////////////////////////////////////////


// doesn't instantiate or resize LGSSM mods because 
// their first time prior might depend on samples
template <size_t np>
Hmm_Rbpf_BS<np>::Hmm_Rbpf_BS(HMMRBPFBSResampStyle rt)
    : m_now(0), m_logLastCondLike(0.0), m_resampTechnique(rt)
{
    std::fill(m_logUnNormWeights.begin(), m_logUnNormWeights.end(), 0.0);
}


template <size_t np>
Hmm_Rbpf_BS<np>::~Hmm_Rbpf_BS()
{
}


template <size_t np>
double Hmm_Rbpf_BS<np>::getLogCondLike() const
{
    return m_logLastCondLike;
}


template <size_t np>
std::vector<Mat> Hmm_Rbpf_BS<np>::getExpectations() const
{
    return m_expectations;
}


template <size_t np>
void Hmm_Rbpf_BS<np>::resampMultinomHRBPF(std::array<FSHMM,np> &oldMods, 
                                       std::array<Vec,np> &oldSamps, 
                                       std::array<double,np> &oldLogWts)
{
    m_resampler.resampHRBPF(oldMods, oldSamps, oldLogWts);
}


template <size_t np>
void Hmm_Rbpf_BS<np>::filter(const Vec &data, 
                            const std::vector<std::function<const Mat(const Vec &x1tProbs, const Vec &x2t)> >& fs)
{

    if( m_now == 0){ // first data point coming
    
        // initialize and update the closed-form mods        
        Vec tmpProbs;
        Mat tmpTransMat;
        double m1(0.0);
        for(size_t ii = 0; ii < np; ++ii){
            
            m_p_samps[ii] = muSamp(data); 
            tmpProbs = initHMMProbVec(m_p_samps[ii]);
            tmpTransMat = initHMMTransMat(m_p_samps[ii]);
            m_p_innerMods[ii] = FSHMM(tmpProbs, tmpTransMat);
            //m_p_innerMods.emplace_back(tmpProbs, tmpTransMat); 
            updateFSHMM(m_p_innerMods[ii], data, m_p_samps[ii]);
            m_logUnNormWeights[ii] = std::log(m_p_innerMods[ii].getCondLike()); 

            // maximum to be used in likelihood calc
            if(m_logUnNormWeights[ii] > m1)
                m1 = m_logUnNormWeights[ii];
    
        }
        
        // calculate log p(y_1)
        double sumexp(0.0);
        for(size_t p = 0; p < np; ++p){
            sumexp += std::exp(m_logUnNormWeights[p] - m1);
        }
        m_logLastCondLike = m1 + std::log(sumexp) - std::log(static_cast<double>(np));
        
        // calculate expectations before you resample
        m_expectations.resize(fs.size());
        m_dimState = muSamp(data).rows();
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
            for(size_t prtcl = 0; prtcl < np; ++prtcl){ 
                tmp = h(m_p_innerMods[prtcl].getFilterVec(), m_p_samps[prtcl]);
                tmp = tmp.array().log().matrix() + (m_logUnNormWeights[prtcl] - m)*ones;
                numer = numer + tmp.array().exp().matrix();
                denom += std::exp( m_logUnNormWeights[prtcl] - m );
            }
            m_expectations[fId] = numer/denom;
            fId++;
        }

        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == HMMRBPFBSResampStyle::everytime_multinomial)
            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
            
        // advance time step
        m_now ++;
    }
    else { //m_now > 0
        
        // update
        Vec newX2Samp;
        double m1(-1.0/0.0); // for the updates weights
        double m2 = *std::max_element(m_logUnNormWeights.begin(), m_logUnNormWeights.end());
        double sumexpdenom(0.0);
        for(size_t ii = 0; ii < np; ++ii){
            
            newX2Samp = fSamp(m_p_samps[ii], data);
            updateFSHMM(m_p_innerMods[ii], data, newX2Samp);
            sumexpdenom += std::exp(m_logUnNormWeights[ii] - m2);
            
            m_logUnNormWeights[ii] += std::log(m_p_innerMods[ii].getCondLike());

            // update the maximum for the newest stuff/numerator
            if( m_logUnNormWeights[ii] > m1)
                m1 = m_logUnNormWeights[ii];

            m_p_samps[ii] = newX2Samp;
        }
        
        // calculate log p(y_t | y_{1:t-1})
        double sumexpnumer(0.0);
        for(size_t p = 0; p < np; ++p){
            sumexpnumer += std::exp(m_logUnNormWeights[p] - m1);
        }
        m_logLastCondLike = m1 + std::log(sumexpnumer) - m2 - std::log(sumexpdenom);
        
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
            for(size_t prtcl = 0; prtcl < np; ++prtcl){ 
                tmp = h(m_p_innerMods[prtcl].getFilterVec(), m_p_samps[prtcl]);
                tmp = tmp.array().log().matrix() + (m_logUnNormWeights[prtcl] - m)*ones;
                numer = numer + tmp.array().exp().matrix();
                denom += std::exp( m_logUnNormWeights[prtcl] - m );
            }
            m_expectations[fId] = numer/denom;
            fId++;
        }
        
        // resample (unnormalized weights ok)
        if (m_resampTechnique == HMMRBPFBSResampStyle::everytime_multinomial)
            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_logUnNormWeights);
//        else if ( (m_resampTechnique == HMMRBPFResampStyle::ess_multinomial) && (m_ESS < m_percentOfNumPartsThresh * m_numParts) )
//            resampMultinomHRBPF(m_p_innerMods, m_p_samps, m_unNormWeights);
        
        // update time step
        m_now ++;
    }
    
}
#endif //HMM_RBPF_BS_H