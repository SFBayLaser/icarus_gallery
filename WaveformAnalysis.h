/*************************************************************
 *
 * A basic waveform analysis object
 *
 * Adapted from an example main program from Wes Ketchum
 *
 *************************************************************/

#ifndef WaveformAnalysis_h
#define WaveformAnalysis_h

//Gallery event
#include "gallery/Event.h"

// useful art info
#include "canvas/Utilities/InputTag.h"

// Need a forward declaration here...
class TFile;
class TDirectory;

namespace Waveform
{
    
class WaveformAnalysis
{
public:
    WaveformAnalysis(const art::InputTag&, const art::InputTag&);
   ~WaveformAnalysis();
    
    void setupAnalysis(TFile*,const std::string&,const std::string&);
    
    void analyzeEvent(gallery::Event&);
private:
    
    // Define a structure to contain hits
    using HitCandidate_t = struct HitCandidate
    {
        size_t startTick;
        size_t stopTick;
        size_t maxTick;
        size_t minTick;
        float  maxDerivative;
        float  minDerivative;
        float  hitCenter;
        float  hitSigma;
        float  hitHeight;
    };
    
    using HitCandidateVec = std::vector<HitCandidate_t>;
    
    // Calculate the derivative of an input vector
    void getDerivativeVec(const std::vector<float>&, std::vector<float>&) const;
    
    // Function to search a derivative vector for candidate hits
    void findHitCandidates(std::vector<float>::const_iterator, std::vector<float>::const_iterator, size_t, HitCandidateVec&) const;
    void stdFindHitCandidates(std::vector<float>::const_iterator, std::vector<float>::const_iterator, size_t, HitCandidateVec&) const;
    
    std::string   fBranchName;
    art::InputTag fHitProducerTag;
    art::InputTag fWireProducerTag;
    std::string   fWaveName;
    std::string   fHitsName;
    TFile*        fFile;
    TDirectory*   fWaveDir;
    TDirectory*   fHitsDir;
    TTree*        fWaveTree;     // What is the right way to handle the TTree's?
};

}
#endif
