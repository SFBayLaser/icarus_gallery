/*************************************************************
 *
 * A basic waveform analysis object
 *
 * Adapted from an example main program from Wes Ketchum
 *
 *************************************************************/

#ifndef RawWaveformAnalysis_h
#define RawWaveformAnalysis_h

//Gallery event
#include "gallery/Event.h"

// useful art info
#include "canvas/Utilities/InputTag.h"

// Need a forward declaration here...
class TFile;
class TDirectory;

namespace Waveform
{
    
class RawWaveformAnalysis
{
public:
    RawWaveformAnalysis(const art::InputTag&);
   ~RawWaveformAnalysis();
    
    void setupAnalysis(TFile*,const std::string&,const std::string&);
    
    void analyzeEvent(gallery::Event&);
private:
    void  averageInputWaveform(const std::vector<float>&, std::vector<float>&) const;
    void  getErosionAndDilationVecs(const std::vector<float>&, std::vector<float>&, std::vector<float>&) const;
    float fixTheFreakingWaveform(const std::vector<float>&, std::vector<float>&) const;
    
    std::string   fBranchName;
    art::InputTag fRawDigitProducerTag;
    std::string   fWaveName;
    std::string   fHitsName;
    TFile*        fFile;
    TDirectory*   fWaveDir;
    TDirectory*   fHitsDir;
    TTree*        fWaveTree;     // What is the right way to handle the TTree's?
    
    int                fNumBinsToAve;
    std::vector<float> fWeightVec;
    float              fWeightSum;
};

}
#endif
