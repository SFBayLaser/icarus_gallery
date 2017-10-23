/*************************************************************
 *
 * HitAna program
 *
 * This is a simple demonstration of reading a LArSoft file
 * and accessing recob::Hit and recob::Wire information in gallery
 *
 * Wesley Ketchum (wketchum@fnal.gov), May22, 2017
 *
 *************************************************************/

// Our class definition
#include "BasicHitAnalysis.h"

//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono>

//some ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"

//"art" includes (canvas, and gallery)
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"

//"larsoft" object includes
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"

namespace HitAnalysis
{

BasicHitAnalysis::BasicHitAnalysis(const art::InputTag& hitProducerTag, const art::InputTag& wireProducerTag) :
        fHitProducerTag(hitProducerTag),
        fWireProducerTag(wireProducerTag),
        fNumHitsHist(0)
{
    return;
}

BasicHitAnalysis::~BasicHitAnalysis()
{
    // We are presuming that root takes care of deleting all of these items
    fNumHitsHist = 0; //->Delete();
    
    for(size_t idx = 0; idx < 3; idx++)
    {
        fMaxMinRatHistVec[idx] = 0; //->Delete();
        fDeltaMaxTickMinTickHistVec[idx] = 0; //->Delete();
        fDeltaMaxValMinValHistVec[idx] = 0; //->Delete();
        fRawAdcMinMaxHistVec[idx] = 0;
        fPulseHeight[idx] = 0;
        fPulseHeightSingle[idx] = 0;
        fPulseHeightMulti[idx] = 0;
        fChi2DOF[idx] = 0;
        fNumDegFree[idx] = 0;
        fChi2DOFSingle[idx] = 0;
        fHitMult[idx] = 0;
        fHitCharge[idx] = 0;
        fFitWidth[idx]  = 0;
        fHitSumADC[idx]  = 0;
        
        fPulseHVsWidth[idx] = 0;
        fPulseHVsCharge[idx] = 0;
        fPulseHVsHitNo[idx] = 0;
    }
    
    fHitTree = 0;
    
    return;
}

void BasicHitAnalysis::setupAnalysis(TFile& file, const std::string& histName, const std::string& branchName)
{
    // Where are we now?
    TDirectory* curDir = file.GetDirectory("");
    
    // Create some useful histograms
    TDirectory* histDir = file.mkdir(histName.c_str());
    histDir->cd();
    
    fNumHitsHist = new TH1F("NumHits", "# hits", 300, 0., 30000.);
    
    fMaxMinRatHistVec.resize(3);
    fDeltaMaxTickMinTickHistVec.resize(3);
    fDeltaMaxValMinValHistVec.resize(3);
    fPulseHeight.resize(3);
    fPulseHeightSingle.resize(3);
    fPulseHeightMulti.resize(3);
    fChi2DOF.resize(3);
    fNumDegFree.resize(3);
    fChi2DOFSingle.resize(3);
    fHitMult.resize(3);
    fHitCharge.resize(3);
    fFitWidth.resize(3);
    fHitSumADC.resize(3);
    fPulseHVsWidth.resize(3);
    fPulseHVsCharge.resize(3);
    fPulseHVsHitNo.resize(3);
    fRawAdcMinMaxHistVec.resize(3);
    
    for(size_t histIdx = 0; histIdx < 3; histIdx++)
    {
        std::string baseName = "_plane" + std::to_string(histIdx);
        
        fMaxMinRatHistVec[histIdx]           = new TH1F(("MaxMinRat"           + baseName).c_str(), ";ratio",      200,  -8.,   2.);
        fDeltaMaxTickMinTickHistVec[histIdx] = new TH1F(("DeltaMaxTickMinTick" + baseName).c_str(), ";delta Tick", 100,   0., 100.);
        fDeltaMaxValMinValHistVec[histIdx]   = new TH1F(("DeltaMaxValMinVal"   + baseName).c_str(), ";delta ADC",  200, -25.,  75.);
        
        fRawAdcMinMaxHistVec[histIdx]        = new TH2F(("RawMaxVMin"          + baseName).c_str(), ";max;min",    200, -20., 180., 200, -180., 20.);
        
        // From previous analysis
        fPulseHeight[histIdx]                = new TH1F(("PulseHeight"         + baseName).c_str(), "Pulse Height; (ADC)",  300,  0.,  150.);
        fPulseHeightSingle[histIdx]          = new TH1F(("PulseHeightS"        + baseName).c_str(), "Pulse Height; (ADC)",  300,  0.,  150.);
        fPulseHeightMulti[histIdx]           = new TH1F(("PulseHeightM"        + baseName).c_str(), "Pulse Height; (ADC)",  300,  0.,  150.);
        fChi2DOF[histIdx]                    = new TH1F(("Chi2DOF"             + baseName).c_str(), "Chi2DOF",              502, -1.,  250.);
        fNumDegFree[histIdx]                 = new TH1F(("NumDegFree"          + baseName).c_str(), "NDF",                  100,  0.,  100.);
        fChi2DOFSingle[histIdx]              = new TH1F(("Chi2DOFS"            + baseName).c_str(), "Chi2DOF",              502, -1.,  250.);
        fHitMult[histIdx]                    = new TH1F(("HitMult"             + baseName).c_str(), "# hits",                15,  0.,   15.);
        fHitCharge[histIdx]                  = new TH1F(("HitCharge"           + baseName).c_str(), "Charge",              1000., 0., 2000.);
        fFitWidth[histIdx]                   = new TH1F(("FitWidth"            + baseName).c_str(), "Width",                100,  0.,   10.);
        fHitSumADC[histIdx]                  = new TH1F(("SumADC"              + baseName).c_str(), "Sum ADC",             1000,  0., 2000.);

        fPulseHVsWidth[histIdx]              = new TH2F(("PHVsWidth"           + baseName).c_str(), ";PH;Width", 100,  0.,  100., 100,  0., 10.);
        fPulseHVsCharge[histIdx]             = new TH2F(("PHVsChrg"            + baseName).c_str(), ";PH;Q",     100,  0.,  100., 100,  0., 2000.);
        fPulseHVsHitNo[histIdx]              = new TProfile(("PHVsNo"          + baseName).c_str(), ";Hit #;PH", 300,  0.,  300., 0., 25.);
    }
    
    // Create a Tree to store hit information
    curDir->cd();
    fHitTreeObj = std::make_unique<HitTreeObject>();
    fHitTree    = new TTree("tree_hit","HitAna: Hit tree");
    fHitTree->Branch("hit",fHitTreeObj.get(),fHitTreeObj->BranchString().c_str());
    
    return;
}

void BasicHitAnalysis::analyzeEvent(gallery::Event& event)
{
    //time the loop
    auto t_begin = std::chrono::high_resolution_clock::now();
    
    //to get run and event info, you use this "eventAuxillary()" object.
    std::cout << "Processing "
    << "Run " << event.eventAuxiliary().run() << ", "
    << "Event " << event.eventAuxiliary().event() << std::endl;
    
    fHitTreeObj->run = event.eventAuxiliary().run();
    fHitTreeObj->ev  = event.eventAuxiliary().event();
    
    //get handle to hits
//    auto const& hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHitProducerTag);
//    auto const& hit_vec(*hit_handle);
    
    // We'll need our ROI's too
    auto const& wire_handle = event.getValidHandle<std::vector<recob::Wire>>(fWireProducerTag);
    auto const& wireVec(*wire_handle);
    
    //get the associations to the recob::Wire objects
//    art::FindOne<recob::Wire> wire_per_hit(hit_handle,event,fHitProducerTag);
    art::FindMany<recob::Hit>   hits_per_wire(wire_handle,event,fHitProducerTag);
    art::FindOne<raw::RawDigit> rawDigits_per_wire(wire_handle,event,fWireProducerTag);
    
    int numHitsTotal(0);
    
    // Ok, the idea is to loop through the Wire data and then over each ROI
    // We'll recover the hits associated to each ROI and analyze them
    for(size_t wireIdx = 0; wireIdx < wireVec.size(); wireIdx++)
    {
        // Get our wire object
        const recob::Wire& wire = wireVec.at(wireIdx);
        
        // Recover the RawDigit assoicated to this wire
        const raw::RawDigit& rawDigit = rawDigits_per_wire.at(wireIdx).ref();
        
        std::vector<short> rawAdcVec(rawDigit.Samples());
        
        // uncompress the data
        raw::Uncompress(rawDigit.ADCs(), rawAdcVec, rawDigit.Compression());
        
        // Pedestal...
        float pedestal = rawDigit.GetPedestal();
        
        // Recover hits assoicated to this wire
        std::vector<const recob::Hit*> hitVec = hits_per_wire.at(wireIdx);
        
        std::sort(hitVec.begin(),hitVec.end(),[](const auto* left, const auto* right){return left->PeakTime() < right->PeakTime();});
        
        // For now only look at ROI's with hits
        if (hitVec.empty()) continue;
        
        // First hit index...
        int firstHit(0);
        
        // For now use view/plane interchangeable since this is MicroBooNE
        size_t planeIdx = wire.View();
        
        // Now loop over ROI's this wire
        for(const auto& roi : wire.SignalROI().get_ranges())
        {
            // Need ROI start/end times
            int startTick = roi.begin_index();
            int stopTick  = roi.end_index();
            
            // Go through and process hits for this ROI
            while(firstHit < int(hitVec.size()))
            {
                // Assume hits in time order, if this hit's peak time is after ROI stop tick then exit loop
                const recob::Hit* hit = hitVec.at(firstHit);
                
                if (hit->PeakTime() > stopTick) break;
                
                // Extract interesting hit parameters
                const geo::WireID& wireID   = hit->WireID();
                float              chi2DOF  = std::min(hit->GoodnessOfFit(),float(249.8));
                int                numDOF   = hit->DegreesOfFreedom();
                int                hitMult  = hit->Multiplicity();
                float              peakTime = hit->PeakTime();
                float              charge   = hit->Integral();
                float              sumADC   = hit->SummedADC();
                float              hitPH    = std::min(hit->PeakAmplitude(),float(249.8));
                float              hitSigma = hit->RMS();
                size_t             view     = wireID.Plane;
//                size_t             wire     = wireID.Wire;
                
                fPulseHeight[view]->Fill(hitPH, 1.);
                fChi2DOF[view]->Fill(chi2DOF, 1.);
                fNumDegFree[view]->Fill(numDOF, 1.);
                fHitMult[view]->Fill(hitMult, 1.);
                fHitCharge[view]->Fill(charge, 1.);
                fFitWidth[view]->Fill(std::min(float(9.99),hitSigma), 1.);
                fHitSumADC[view]->Fill(sumADC, 1.);
                fPulseHVsHitNo[view]->Fill(hitPH, hit->LocalIndex(), 1.);
                
                if (hitMult == 1)
                {
                    fPulseHeightSingle[view]->Fill(hitPH, 1.);
                    fChi2DOFSingle[view]->Fill(chi2DOF, 1.);
                    fPulseHVsWidth[view]->Fill(std::min(float(99.9),hitPH), std::min(float(9.99),hitSigma), 1.);
                    fPulseHVsCharge[view]->Fill(std::min(float(99.9),hitPH), std::min(float(1999.),charge), 1.);
                }
                else
                    fPulseHeightMulti[view]->Fill(hitPH, 1.);
                
                //now fill the hit info
                fHitTreeObj->start_tick   = hit->StartTick();
                fHitTreeObj->end_tick     = hit->EndTick();
                fHitTreeObj->roi_start    = startTick;
                fHitTreeObj->roi_end      = stopTick;
                fHitTreeObj->roi_size     = roi.size();
                fHitTreeObj->ch           = hit->Channel();
                fHitTreeObj->time         = peakTime;
                fHitTreeObj->time_err     = hit->SigmaPeakTime();
                fHitTreeObj->rms          = hitSigma;
                fHitTreeObj->amp          = hit->PeakAmplitude();
                fHitTreeObj->amp_err      = hit->SigmaPeakAmplitude();
                fHitTreeObj->sumadc       = sumADC;
                fHitTreeObj->integral     = charge;
                fHitTreeObj->integral_err = hit->SigmaIntegral();
                fHitTreeObj->gof          = chi2DOF;
                fHitTreeObj->mult         = hitMult;
                fHitTreeObj->idx          = hit->LocalIndex();
                fHitTreeObj->view         = view;
                
                firstHit++;
                numHitsTotal++;
                
                fHitTree->Fill();
            }
            
            // Look at the associated RawDigit in the range of this roi
            std::vector<short>::iterator maxItr = std::max_element(rawAdcVec.begin()+startTick,rawAdcVec.begin()+stopTick);
            std::vector<short>::iterator minItr = std::min_element(rawAdcVec.begin()+startTick,rawAdcVec.begin()+stopTick);
            
            float rawMaxVal = float(*maxItr) - pedestal;
            float rawMinVal = float(*minItr) - pedestal;
            
            if (planeIdx == 1 && rawMaxVal > 15. && rawMinVal > -10.) std::cout << "**> channel: " << wire.Channel() << ", rawMaxVal: " << rawMaxVal << ", min: " << rawMinVal << std::endl;
            
            fRawAdcMinMaxHistVec[planeIdx]->Fill(rawMaxVal, rawMinVal, 1.);
            
            // Test idea of detecting bad waveforms from the ROI waveform
            std::vector<float>::const_iterator roiMinItr = std::min_element(roi.data().begin(),roi.data().end());
            
            // search backwards from the maximum to look for an inflection point
            std::vector<float>::const_iterator roiMaxItr = roiMinItr;
            
            // Loop to look for the point after the minimum where the curve turns over
            while(++roiMaxItr != roi.data().end())
            {
                if (*roiMaxItr > 0. && *(roiMaxItr-2) > *(roiMaxItr-1) && *(roiMaxItr-1) > *roiMaxItr)
                {
                    roiMaxItr -= 2;
                    break;
                }
            }
            
            std::vector<float>::const_iterator roiMidItr = roiMinItr + std::distance(roiMinItr,roiMaxItr)/2;
            
            int minTick = startTick + std::distance(roi.data().begin(),roiMinItr);
//            int midTick = startTick + std::distance(roi.data().begin(),roiMidItr);
            int maxTick = startTick + std::distance(roi.data().begin(),roiMaxItr);
            
            int   deltaMaxTickMinTick = maxTick - minTick;
            float deltaMaxValMinVal   = *roiMaxItr + *roiMinItr;
            float maxValMinValRat     = std::fabs(*roiMaxItr) > 0. ? *roiMinItr / *roiMaxItr : 0.;
            
//            if (deltaMaxValMinVal > 0.) continue;
            
            if (planeIdx == 1 && *roiMinItr < -5.) std::cout << "--> channel: " << wire.Channel() << ", min/max: " << *roiMinItr << "/" << *roiMaxItr << ", first: " << std::distance(roi.data().begin(),roiMinItr) << ", roiLen: " << roi.data().size() << std::endl;
            
            fDeltaMaxTickMinTickHistVec[planeIdx]->Fill(deltaMaxTickMinTick, 1.);
            fDeltaMaxValMinValHistVec[planeIdx]->Fill(deltaMaxValMinVal, 1.);
            fMaxMinRatHistVec[planeIdx]->Fill(maxValMinValRat, 1.);
        }
    }
    
    fNumHitsHist->Fill(numHitsTotal, 1.);
    
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> time_total_ms(t_end-t_begin);
    std::cout << "\tEvent took " << time_total_ms.count() << " ms to process." << std::endl;
    
    return;
}


} // namespace
