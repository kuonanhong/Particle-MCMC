
template<size_t N>
MNResamp<N>::MNResamp() 
        : m_gen{static_cast<std::uint32_t>(
                    std::chrono::high_resolution_clock::now().time_since_epoch().count()
                                           )}
{
}


template<size_t N>
MNResamp<N>::~MNResamp()
{
}


template<size_t N>
void MNResamp<N>::resampLogWts(std::array<Vec, N> &oldParts, std::array<double, N> &oldLogUnNormWts)
{
    // these log weights may be very negative. If that's the case, exponentiating them may cause underflow
    // so we use the "log-exp-sum" trick
    // actually not quite...we just shift the log-weights because after they're exponentiated
    // they have the same normalized probabilities
    
    // get dimensions
    unsigned int numParticles = oldParts.size();
    
    // Create the distribution with exponentiated log-weights
    std::vector<double> w;
    w.resize(oldLogUnNormWts.size());
    double m = *std::max_element(oldLogUnNormWts.begin(), oldLogUnNormWts.end());
    std::transform(oldLogUnNormWts.begin(), oldLogUnNormWts.end(), w.begin(), 
                    [&m](double& d) -> double { return std::exp( d - m ); } );
    std::discrete_distribution<> idxSampler(w.begin(), w.end());
    
    // create temporary particle vector and weight vector
    std::array<Vec, N> tmpPartics = oldParts; // TODO: check this copies!
    
    // sample from the original parts and store in tmpParts
    unsigned int whichPart;
    for(size_t part = 0; part < numParticles; ++part)
    {
        whichPart = idxSampler(m_gen);
        tmpPartics[part] = oldParts[whichPart];
    }
        
    //overwrite olds with news
    oldParts = std::move(tmpPartics);
    std::fill(oldLogUnNormWts.begin(), oldLogUnNormWts.end(), 0.0); // change back    
}


