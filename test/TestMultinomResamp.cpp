#include "UnitTest++.h"
#include "multinomial_resampler.h"
#include "mn_resampler.h"

class MRFixture
{
public:
    MultinomResamp m_mr; // only default constructors available
    MNResamp<2> m_mr2;
    
    // TestMultinomResamp_kGen
    std::vector<double> m_kw;
    std::vector<unsigned int> m_returned_ks;
    std::array<double, 2> m_kw2;
    std::array<unsigned int, 2> m_returned_ks2;
    
    // for TestMultinomResamp_resampLogWtsVecVec
    // no mnresamp analog
    std::vector<std::vector<Vec>> m_vvparts;
    std::vector<double> m_vvw;
    
    // for TestMultinomResamp_resampLogWtsVec
    std::vector<Vec> m_vparts;
    std::vector<double> m_vw;
    std::array<Vec,2> m_vparts2;
    std::array<double,2> m_vw2;
    
    // for TestMultinomResamp_ressampHRBPF
    std::vector<FSHMM> m_hmms;
    std::vector<Vec> m_hmmparts;
    std::vector<double> m_hmmwts;
    std::array<FSHMM,2> m_hmms2;
    std::array<Vec,2> m_hmmparts2;
    std::array<double,2> m_hmmwts2;
    
    // for TestMultinomResamp_ressampKRBPF
    std::vector<Lgssm> m_lgssms;
    std::vector<Vec> m_lgssmsamps;
    std::vector<double> m_lgssmwts;
    std::array<Lgssm,2> m_lgssms2;
    std::array<Vec,2> m_lgssmsamps2;
    std::array<double,2> m_lgssmwts2;


    MRFixture(){
        
        // TestMultinomResamp_kGen
        m_kw.push_back(0.0);
        m_kw.push_back(-1.0/0.0);
        m_returned_ks.resize(2);
        m_kw2[0] = 0.0;
        m_kw2[1] = -1.0/0.0;
        
        // for TestMultinomResamp_resampLogWtsVecVec
        m_vvparts.resize(1);
        m_vvparts[0].push_back(Vec::Constant(1,1,0.0));
        m_vvparts[0].push_back(Vec::Constant(1,1,1.0));
        m_vvw.push_back(0.0);
        m_vvw.push_back(-1.0/0.0);
        
        // for TestMultinomResamp_resampLogWtsVec
        m_vparts.push_back(Vec::Constant(1,1,0.0));
        m_vparts.push_back(Vec::Constant(1,1,1.0));
        m_vw.push_back(0.0);
        m_vw.push_back(-1.0/0.0);
        m_vparts2[0] = Vec::Constant(1,1,0.0);
        m_vparts2[1] = Vec::Constant(1,1,1.0);
        m_vw2[0] = 0.0;
        m_vw2[1] = -1.0/0.0;
        
        // for TestMultinomResamp_ressampHRBPF
        m_hmms.push_back(FSHMM(Vec::Constant(1,1,1.0), 
                               Mat::Constant(1,1,1.0)));
        m_hmms.push_back(FSHMM(Vec::Constant(1,1,2.0), 
                               Mat::Constant(1,1,2.0)));
        m_hmmparts.push_back(Vec::Constant(1,1,1.0));
        m_hmmparts.push_back(Vec::Constant(2,1,2.0));
        m_hmmwts.push_back(0.0);
        m_hmmwts.push_back(-1.0/0.0);
        m_hmms2[0] = FSHMM(Vec::Constant(1,1,1.0), 
                               Mat::Constant(1,1,1.0));
        m_hmms2[1] = FSHMM(Vec::Constant(1,1,2.0), 
                               Mat::Constant(1,1,2.0));
        m_hmmparts2[0] = Vec::Constant(1,1,1.0);
        m_hmmparts2[1] = Vec::Constant(2,1,2.0);
        m_hmmwts2[0] = 0.0;
        m_hmmwts2[1] = -1.0/0.0;
        
        // for TestMultinomResamp_ressampKRBPF
        m_lgssms.emplace_back(Vec::Constant(1,1,1.0),
                              Mat::Constant(1,1,1.0));
        m_lgssms.emplace_back(Vec::Constant(1,1,2.0),
                              Mat::Constant(1,1,2.0));
        m_lgssmsamps.push_back(Vec::Constant(2,1,1.0)); 
        m_lgssmsamps.push_back(Vec::Constant(2,1,2.0));
        m_lgssmwts.push_back(0.0);
        m_lgssmwts.push_back(-1.0/0.0);
        m_lgssms2[0] = Lgssm(Vec::Constant(1,1,1.0),
                              Mat::Constant(1,1,1.0));
        m_lgssms2[1] = Lgssm(Vec::Constant(1,1,2.0),
                              Mat::Constant(1,1,2.0));
        m_lgssmsamps2[0] = Vec::Constant(2,1,1.0);
        m_lgssmsamps2[1] = Vec::Constant(2,1,2.0);
        m_lgssmwts2[0] = 0.0;
        m_lgssmwts2[1] = -1.0/0.0;
    }
    
};


TEST_FIXTURE(MRFixture, TestMultinomResamp_kGen)
{
    m_mr.kGen(m_kw, m_returned_ks);
    CHECK_EQUAL(m_returned_ks[0], 0);
    CHECK_EQUAL(m_returned_ks[1], 0);

    m_mr2.kGen(m_kw2, m_returned_ks2);
    CHECK_EQUAL(m_returned_ks2[0], 0);
    CHECK_EQUAL(m_returned_ks2[1], 0);
}

TEST_FIXTURE(MRFixture, TestMultinomResamp_resampLogWtsVecVec)
{
    m_mr.resampLogWts(m_vvparts, m_vvw);
    CHECK_EQUAL(m_vvparts[0][0](0), 0.0);
    CHECK_EQUAL(m_vvparts[0][1](0), 0.0);
    CHECK_EQUAL(m_vvw[0], 0.0);
    CHECK_EQUAL(m_vvw[1], 0.0);
}

TEST_FIXTURE(MRFixture, TestMultinomResamp_resampLogWtsVec)
{
    m_mr.resampLogWts(m_vparts, m_vw);
    CHECK_EQUAL(m_vparts[0](0), 0.0);
    CHECK_EQUAL(m_vparts[1](0), 0.0);
    CHECK_EQUAL(m_vw[0], 0.0);
    CHECK_EQUAL(m_vw[1], 0.0);

    m_mr2.resampLogWts(m_vparts2, m_vw2);
    CHECK_EQUAL(m_vparts2[0](0), 0.0);
    CHECK_EQUAL(m_vparts2[1](0), 0.0);
    CHECK_EQUAL(m_vw2[0], 0.0);
    CHECK_EQUAL(m_vw2[1], 0.0);
}

TEST_FIXTURE(MRFixture, TestMultinomResamp_resampHRBPF)
{
    m_mr.resampHRBPF(m_hmms, m_hmmparts, m_hmmwts);
    // TODO: test that the MODS were resampled correctly
    CHECK_EQUAL(m_hmmparts[0](0), 1.0);
    CHECK_EQUAL(m_hmmparts[1](0), 1.0);
    CHECK_EQUAL(m_hmmwts[0], 0.0);
    CHECK_EQUAL(m_hmmwts[1], 0.0);
    
    m_mr2.resampHRBPF(m_hmms2, m_hmmparts2, m_hmmwts2);
    // TODO: test that the MODS were resampled correctly
    CHECK_EQUAL(m_hmmparts2[0](0), 1.0);
    CHECK_EQUAL(m_hmmparts2[1](0), 1.0);
    CHECK_EQUAL(m_hmmwts2[0], 0.0);
    CHECK_EQUAL(m_hmmwts2[1], 0.0);

}

TEST_FIXTURE(MRFixture, TestMultinomResamp_resampKRBPF)
{
    m_mr.resampKRBPF(m_lgssms, m_lgssmsamps, m_lgssmwts);
    // TODO: check MODS are correct
    CHECK_EQUAL(m_lgssmsamps[0](0), 1.0);
    CHECK_EQUAL(m_lgssmsamps[1](0), 1.0);
    CHECK_EQUAL(m_lgssmwts[0], 0.0);
    CHECK_EQUAL(m_lgssmwts[1], 0.0);
    
    m_mr2.resampKRBPF(m_lgssms2, m_lgssmsamps2, m_lgssmwts2);
    // TODO: check MODS are correct
    CHECK_EQUAL(m_lgssmsamps2[0](0), 1.0);
    CHECK_EQUAL(m_lgssmsamps2[1](0), 1.0);
    CHECK_EQUAL(m_lgssmwts2[0], 0.0);
    CHECK_EQUAL(m_lgssmwts2[1], 0.0);
}
