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
#include "WaveformAnalysis.h"

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
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"

namespace Waveform
{

WaveformAnalysis::WaveformAnalysis(const art::InputTag& hitProducerTag, const art::InputTag& wireProducerTag) :
        fHitProducerTag(hitProducerTag),
        fWireProducerTag(wireProducerTag),
        fWaveName(""),
        fHitsName(""),
        fFile(0)
{
    return;
}

WaveformAnalysis::~WaveformAnalysis()
{
    return;
}

void WaveformAnalysis::setupAnalysis(TFile* tFile,const std::string& waveDir,const std::string& hitsDir)
{
    // Keep track of input info
    fFile     = tFile;
    fWaveName = waveDir;
    fHitsName = hitsDir;
    
    fWaveDir  = fFile->mkdir(fWaveName.c_str());
    fHitsDir  = fFile->mkdir(fHitsName.c_str());
    
    //Create a Tree to store event information
    fWaveDir->cd();
    fWaveTree = new TTree("waveforms","WaveformAnalysis: hists of waveforms");
    
    return;
}

void WaveformAnalysis::analyzeEvent(gallery::Event& event)
{
    //time the loop
    auto t_begin = std::chrono::high_resolution_clock::now();
    
    //get handle to hits
    auto const& hit_handle = event.getValidHandle<std::vector<recob::Hit>>(fHitProducerTag);
    auto const& hit_vec(*hit_handle);
    
    // We'll need our RawDigits too
//    auto const& wire_handle = event.getValidHandle<std::vector<recob::Wire>>(fWireProducerTag);
//    auto const& wire_vec(*wire_handle);
    
    //get the associations to the recob::Wire objects
    art::FindOne<recob::Wire> wire_per_hit(hit_handle,event,fHitProducerTag);
//    art::FindMany<recob::Hit> hits_per_wire(wire_handle,event,fHitProducerTag);
    
    // Define a map to contain mapping from ROI's to hits on that ROI
    std::map<const lar::sparse_vector<float>::datarange_t*,std::vector<const recob::Hit*>> roiToHitMap;
    std::map<raw::ChannelID_t, const recob::Wire*>                                         channelToWireMap;
    
    //loop over the hits
    for(size_t i_h = 0; i_h < hit_vec.size(); ++i_h)
    {
        const recob::Hit& hit = hit_vec[i_h];
        
        //get associated wire object
        const recob::Wire&                            wire = (wire_per_hit.at(i_h)).ref();
        const lar::sparse_vector<float>::datarange_t& roi  = wire.SignalROI().find_range((int)(hit.PeakTime()));
        
        roiToHitMap[&roi].push_back(&hit);
        channelToWireMap[hit.Channel()] = &wire;
    }//end first loop over hits.
    
    // Make a maximum size container to store waveform
    // Need to use a fixed length array here to insure continuity
    // of memory locations
    struct ROIObject
    {
        unsigned channel;
        unsigned plane;
        int      startTick;
        int      startMax;
        float    maxValue;
        int      startMin;
        float    minValue;
        int      midTick;
        float    midValue;
        int      bufferSize;
        float    buffer[6400];
    } roiObject;
    
    // We want to store candidate hits too (another variable size object!)
    struct HitObject
    {
        unsigned channel;
        unsigned plane;
        int      numHits;
        int      startTick[500];
        int      stopTick[500];
        int      maxTick[500];
        int      minTick[500];
        float    maxDerivative[500];
        float    minDerivative[500];
        float    hitCenter[500];
        float    hitSigma[500];
        float    hitHeight[500];
    } hitObject;
    
    // Now go through the map an make plots (lots of them)
    for(const auto& chanWirePair : channelToWireMap)
    {
        for(const auto& roi : chanWirePair.second->SignalROI().get_ranges())
        {
            if (roiToHitMap.find(&roi) == roiToHitMap.end()) continue;
            
            const std::vector<float>&             roiSignalVec    = roi.data();
            raw::TDCtick_t                        roiFirstBinTick = roi.begin_index();
            const std::vector<const recob::Hit*>& roiHitVector    = roiToHitMap.at(&roi);
            HitCandidateVec                       hitCandidateVec;
            std::vector<float>                    derivativeVec;
        
            // Get the derivative vector (so we only calculate once)
            getDerivativeVec(roiSignalVec,derivativeVec);
        
            if (roiHitVector.front()->Channel() > 345 && roiHitVector.front()->Channel() < 348 && roiFirstBinTick > 2580)
            {
                std::cout << "Channel " << roiHitVector.front()->Channel() << ", first tick: " << roiFirstBinTick << std::endl;
            }
        
            // Find the candidate hits
            findHitCandidates(derivativeVec.begin(), derivativeVec.end(), roiFirstBinTick, hitCandidateVec);
            //stdFindHitCandidates(roiSignalVec.begin(), roiSignalVec.end(), roiFirstBinTick, hitCandidateVec);
        
            size_t bufferSize = roiSignalVec.size();

            // Copy waveform the old fashioned way...
            for(size_t waveformIdx = 0; waveformIdx < bufferSize; waveformIdx++)
                roiObject.buffer[waveformIdx] = roiSignalVec.at(waveformIdx);
        
            // Test idea of detecting bad waveforms
            std::vector<float>::const_iterator roiMaxItr = std::max_element(roiSignalVec.begin(),roiSignalVec.end());
            //std::vector<float>::const_iterator roiMinItr = std::min_element(roiSignalVec.begin(),roiMaxItr);
        
            // search backwards from the maximum to look for an inflection point
            std::vector<float>::const_iterator roiMinItr = roiMaxItr;
        
            while(--roiMinItr != roiSignalVec.begin())
            {
                if (*roiMinItr < 0.5 * *roiMaxItr && *(roiMinItr-1) > *roiMinItr) break;
            }
        
            std::vector<float>::const_iterator roiMidItr = roiMinItr + std::distance(roiMinItr,roiMaxItr)/2;

            // Now fill out the struct before handing over to root
            roiObject.channel    = roiHitVector.front()->Channel();
            roiObject.plane      = roiHitVector.front()->View();
            roiObject.startTick  = roiFirstBinTick;
            roiObject.startMax   = std::distance(roiSignalVec.begin(),roiMaxItr);
            roiObject.maxValue   = *roiMaxItr;
            roiObject.startMin   = std::distance(roiSignalVec.begin(),roiMinItr);
            roiObject.minValue   = *roiMinItr;
            roiObject.midTick    = std::distance(roiSignalVec.begin(),roiMidItr);
            roiObject.midValue   = *roiMidItr;
            roiObject.bufferSize = bufferSize;
        
            // Set up to store in root
            std::string branchName = "SignalVector_" + std::to_string(roiHitVector.front()->Channel()) + "_" + std::to_string(roiFirstBinTick);
            std::string branchDef  = "channel/i:plane/i:startTick/i:startMax/i:maxValue/F:startMin/i:minValue/F:midTick/i:midvalue/F:bufferSize/I:buffer[bufferSize]/F";

            fWaveDir->cd();
        
            TBranch* buffBranch = fWaveTree->Branch(branchName.c_str(),&roiObject,branchDef.c_str());
        
            buffBranch->Fill();
        
            // Copy candidate hit object info to a format more friendly to root
            hitObject.channel = roiHitVector.front()->Channel();
            hitObject.plane   = roiHitVector.front()->View();
            hitObject.numHits = hitCandidateVec.size();
        
            size_t hitIdx(0);
        
            for(const auto& hitCandidate : hitCandidateVec)
            {
                hitObject.startTick[hitIdx]     = hitCandidate.startTick;
                hitObject.stopTick[hitIdx]      = hitCandidate.stopTick;
                hitObject.maxTick[hitIdx]       = hitCandidate.maxTick;
                hitObject.minTick[hitIdx]       = hitCandidate.minTick;
                hitObject.maxDerivative[hitIdx] = hitCandidate.maxDerivative;
                hitObject.minDerivative[hitIdx] = hitCandidate.minDerivative;
                hitObject.hitCenter[hitIdx]     = hitCandidate.hitCenter;
                hitObject.hitSigma[hitIdx]      = hitCandidate.hitSigma;
                hitObject.hitHeight[hitIdx]     = hitCandidate.hitHeight;
                hitIdx++;
            }
        
            fHitsDir->cd();
            TTree* hitsTree = new TTree(branchName.c_str(),"WaveformAnalysis: candidate hits");
        
            hitsTree->Branch("channel",       &hitObject.channel,       "channel/i");
            hitsTree->Branch("plane",         &hitObject.plane,         "plane/i");
            hitsTree->Branch("numHits",       &hitObject.numHits,       "numHits/i");
            hitsTree->Branch("startTick",     &hitObject.startTick,     "startTick[numHits]/i");
            hitsTree->Branch("stopTick",      &hitObject.stopTick,      "stopTick[numHits]/i");
            hitsTree->Branch("maxTick",       &hitObject.maxTick,       "maxTick[numHits]/i");
            hitsTree->Branch("minTick",       &hitObject.minTick,       "minTick[numHits]/i");
            hitsTree->Branch("maxDerivative", &hitObject.maxDerivative, "maxDerivative[numHits]/F");
            hitsTree->Branch("minDerivative", &hitObject.minDerivative, "minDerivative[numHits]/F");
            hitsTree->Branch("hitCenter",     &hitObject.hitCenter,     "hitCenter[numHits]/F");
            hitsTree->Branch("hitSigma",      &hitObject.hitSigma,      "hitSigma[numHits]/F");
            hitsTree->Branch("hitHeight",     &hitObject.hitHeight,     "hitHeight[numHits]/F");

            hitsTree->Fill();
        }
    }
    
    fWaveTree->Fill();
    
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> time_total_ms(t_end-t_begin);
    std::cout << "\tEvent took " << time_total_ms.count() << " ms to process." << std::endl;
    
    return;
}
    
void WaveformAnalysis::getDerivativeVec(const std::vector<float>& roiSignalVec, std::vector<float>& derivativeVec) const
{
    // First task is to compute the derivative of the input signal vector
    std::vector<float> tempVec;
    
    tempVec.resize(roiSignalVec.size(),0.);
    
    for(size_t idx = 1; idx < roiSignalVec.size()-1; idx++)
        tempVec[idx] = 0.5 * (roiSignalVec.at(idx+1) - roiSignalVec.at(idx-1));
    
    // Now smooth the derivative vector
    derivativeVec.resize(roiSignalVec.size(),0.);
    
    for(size_t idx = 2; idx < tempVec.size() - 2; idx++)
        derivativeVec[idx] = (tempVec.at(idx-2) + 2.*tempVec.at(idx-1) + 3.*tempVec.at(idx) + 2.*tempVec.at(idx+1) + tempVec.at(idx+2))/9.;
    
    return;
}

    
void WaveformAnalysis::findHitCandidates(std::vector<float>::const_iterator startItr,
                                         std::vector<float>::const_iterator stopItr,
                                         size_t                             roiStartTick,
                                         HitCandidateVec&                   hitCandidateVec) const
{
    // Require a minimum length
    size_t roiLength = std::distance(startItr,stopItr);
    
    if (roiLength > 2)
    {
        // Now search for candidate hits... our first is to simply find the range between the maximum and minimum bins in the
        // input derivative vector. We'll consider special cases further down.
        // So, start by looking for the peak bin, peak positive over baseline, then search from there for the largest negative bin
        std::vector<float>::const_iterator maxItr = std::max_element(startItr,stopItr);
        std::vector<float>::const_iterator minItr = std::min_element(maxItr,  stopItr);
        
        float range = *maxItr - *minItr;

        // At some point small rolling oscillations on the waveform need to be ignored...
        if (range > 0.1)
        {
            // Need to back up to find zero crossing
            std::vector<float>::const_iterator newEndItr = maxItr;
        
            while(newEndItr != startItr && *newEndItr > 0.) newEndItr--;
            
            size_t startTick = std::distance(startItr,newEndItr);
        
            // Find hits in the section of the waveform leading up to this candidate hit
            if (startTick > 2) findHitCandidates(startItr,newEndItr,roiStartTick,hitCandidateVec);
        
            // Now need to go forward to again get close to zero
            std::vector<float>::const_iterator newStartItr = minItr;
        
            while(newStartItr != stopItr && *newStartItr < 0.) newStartItr++;
            
            size_t stopTick = std::distance(startItr,newStartItr);
            
            // The range from the max derivative to the min derivative found by a simple
            // search above may contain several intervening max/min. This can occur when
            // a delta ray is ejected and you have two hits overlapping. So we need to
            // step through from max to min and look for these special cases.
            std::vector<std::vector<float>::const_iterator> maxItrVec = {maxItr};
            std::vector<std::vector<float>::const_iterator> minItrVec;
            
            std::vector<float>::const_iterator adcItr  = maxItr;
            std::vector<float>::const_iterator lastItr = maxItr;
            bool                               foundMin(false);
            
            // We are gauranteed that we start at the maximum and end at the minimum. If
            // we find an intervening minimum then there MUST be a corresponding maximum
            // to get the derivative turning down again to reach the absolute minimum
            // This loop is constructed with that assumption in mind
            while(++adcItr != minItr)
            {
                // If we found an intervening minimum then switch modes to looking for next max
                if (foundMin)
                {
                    // Peak condition is that the waveform slope is decreasing again
                    if (*adcItr < *lastItr)
                    {
                        maxItrVec.push_back(lastItr);
                        foundMin = false;
                    }
                }
                // Condition for a minimum is that the waveform slope starts increasing
                else if (*adcItr > *lastItr)
                {
                    minItrVec.push_back(lastItr);
                    foundMin = true;
                }
                lastItr = adcItr;
            }
            
            minItrVec.push_back(minItr);
            
            // Fill candidate hit vector
            for(size_t idx = 0; idx < maxItrVec.size(); idx++)
            {
                // get a new hit object and fill it
                HitCandidate_t hitCandidate;
                
                std::vector<float>::const_iterator peakItr = std::min_element(maxItrVec.at(idx),minItrVec.at(idx),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
                hitCandidate.startTick     = roiStartTick + startTick;
                hitCandidate.stopTick      = roiStartTick + stopTick;
                hitCandidate.maxTick       = roiStartTick + std::distance(startItr,maxItrVec.at(idx));
                hitCandidate.minTick       = roiStartTick + std::distance(startItr,minItrVec.at(idx));
                hitCandidate.maxDerivative = maxItrVec.at(idx) != stopItr ? *maxItrVec.at(idx) : 0.;
                hitCandidate.minDerivative = minItrVec.at(idx) != stopItr ? *minItrVec.at(idx) : 0.;
                hitCandidate.hitCenter     = roiStartTick + std::distance(startItr,peakItr); //0.5 * float(hitCandidate.maxTick + hitCandidate.minTick);
                hitCandidate.hitSigma      = 0.5 * float(hitCandidate.minTick - hitCandidate.maxTick);
                hitCandidate.hitHeight     = hitCandidate.hitSigma * (hitCandidate.maxDerivative - hitCandidate.minDerivative) / 1.2130; //0.6065;
    
                hitCandidateVec.push_back(hitCandidate);
            }
            
            // Finally, search the section of the waveform following this candidate for more hits
            findHitCandidates(newStartItr,stopItr,roiStartTick + stopTick,hitCandidateVec);
        }
    }
    
    return;
}
    
// --------------------------------------------------------------------------------------------
// Initial finding of candidate peaks
// --------------------------------------------------------------------------------------------
void WaveformAnalysis::stdFindHitCandidates(std::vector<float>::const_iterator startItr,
                                            std::vector<float>::const_iterator stopItr,
                                            size_t                             roiStartTick,
                                            HitCandidateVec&                   hitCandidateVec) const
{
    static const float roiThreshold(3.2);
    
    // Need a minimum number of ticks to do any work here
    if (std::distance(startItr,stopItr) > 4)
    {
        // Find the highest peak in the range given
        auto maxItr = std::max_element(startItr, stopItr);
        
        float maxValue = *maxItr;
        int   maxTime  = std::distance(startItr,maxItr);
        
        if (maxValue > roiThreshold)
        {
            // backwards to find first bin for this candidate hit
            auto firstItr = std::distance(startItr,maxItr) > 2 ? maxItr - 1 : startItr;
            
            while(firstItr != startItr)
            {
                // Check for pathology where waveform goes too negative
                if (*firstItr < -roiThreshold) break;
                
                // Check both sides of firstItr and look for min/inflection point
                if (*firstItr < *(firstItr+1) && *firstItr <= *(firstItr-1)) break;
                
                firstItr--;
            }
            
            int firstTime = std::distance(startItr,firstItr);
            
            // Recursive call to find all candidate hits earlier than this peak
            stdFindHitCandidates(startItr, firstItr + 1, roiStartTick, hitCandidateVec);
            
            // forwards to find last bin for this candidate hit
            auto lastItr = std::distance(maxItr,stopItr) > 2 ? maxItr + 1 : stopItr - 1;
            
            while(lastItr != stopItr - 1)
            {
                // Check for pathology where waveform goes too negative
                if (*lastItr < -roiThreshold) break;
                
                // Check both sides of firstItr and look for min/inflection point
                if (*lastItr <= *(lastItr+1) && *lastItr < *(lastItr-1)) break;
                
                lastItr++;
            }
            
            int lastTime = std::distance(startItr,lastItr);
            
            // Now save this candidate's start and max time info
            HitCandidate hitCandidate;
            hitCandidate.startTick     = roiStartTick + firstTime;
            hitCandidate.stopTick      = roiStartTick + lastTime;
            hitCandidate.maxTick       = roiStartTick + firstTime;
            hitCandidate.minTick       = roiStartTick + lastTime;
            hitCandidate.maxDerivative = *(startItr + firstTime);
            hitCandidate.minDerivative = *(startItr + lastTime);
            hitCandidate.hitCenter     = roiStartTick + maxTime;
            hitCandidate.hitSigma      = 0.5 * float(lastTime - firstTime);
            hitCandidate.hitHeight     = *(startItr + maxTime);
            
            hitCandidateVec.push_back(hitCandidate);
            
            // Recursive call to find all candidate hits later than this peak
            stdFindHitCandidates(lastItr + 1, stopItr, roiStartTick + std::distance(startItr,lastItr + 1), hitCandidateVec);
        }
    }
    
    return;
}


} // namespace
