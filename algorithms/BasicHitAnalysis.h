/*************************************************************
 *
 * A basic hit analysis object
 *
 * Adapted from an example main program from Wes Ketchum
 *
 *************************************************************/

#ifndef BasicHitAnalysis_h
#define BasicHitAnalysis_h

//Gallery event
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"

// Also need TTree
//#include "TTree.h"

class TFile;
class TTree;
class TH1F;
class TH2F;
class TProfile;

namespace HitAnalysis
{
    
class BasicHitAnalysis
{
public:
    BasicHitAnalysis(const art::InputTag&, const art::InputTag&);
   ~BasicHitAnalysis();
    
    void setupAnalysis(TFile&, const std::string&, const std::string&);
    
    void analyzeEvent(gallery::Event&);
private:
    
    struct HitTreeObject
    {
        unsigned int run;
        unsigned int ev;
        unsigned int ch;
        unsigned int start_tick;
        unsigned int end_tick;
        unsigned int roi_start;
        unsigned int roi_end;
        unsigned int roi_size;
        float        time;
        float        time_err;
        float        rms;
        float        amp;
        float        amp_err;
        float        sumadc;
        float        integral;
        float        integral_err;
        float        gof;
        int          mult;
        int          idx;
        int          view;
        
        //function to return branch list"
        std::string BranchString(){
            std::string str("run/i:ev/i:ch/i:start_tick/i:end_tick/i:roi_start/i:roi_end/i:roi_size/i:time/F:time_err/F:rms/F:amp/F:amp_err/F:sumadc/F:integral/F:integral_err/F:gof/F:mult/I:idx/I:view/I");
            return str;
        }
    };
    
    std::string                      fBranchName;
    art::InputTag                    fHitProducerTag;
    art::InputTag                    fWireProducerTag;
    std::unique_ptr<HitTreeObject>   fHitTreeObj;
    TTree*                           fHitTree;
    
    TH1F*                            fNumHitsHist;
    std::vector<TH1F*>               fMaxMinRatHistVec;
    std::vector<TH1F*>               fDeltaMaxTickMinTickHistVec;
    std::vector<TH1F*>               fDeltaMaxValMinValHistVec;
    std::vector<TH1F*>               fPulseHeight;
    std::vector<TH1F*>               fPulseHeightSingle;
    std::vector<TH1F*>               fPulseHeightMulti;
    std::vector<TH1F*>               fChi2DOF;
    std::vector<TH1F*>               fNumDegFree;
    std::vector<TH1F*>               fChi2DOFSingle;
    std::vector<TH1F*>               fHitMult;
    std::vector<TH1F*>               fHitCharge;
    std::vector<TH1F*>               fFitWidth;
    std::vector<TH1F*>               fHitSumADC;
    std::vector<TH2F*>               fPulseHVsWidth;
    std::vector<TH2F*>               fPulseHVsCharge;
    std::vector<TProfile*>           fPulseHVsHitNo;
    std::vector<TH2F*>               fRawAdcMinMaxHistVec;
};

}
#endif
