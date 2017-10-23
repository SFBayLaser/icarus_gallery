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
#include "RawWaveformAnalysis.h"

//some standard C++ includes
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <chrono>
#include <algorithm>
#include <numeric>

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
#include "lardataobj/RawData/raw.h"

namespace Waveform
{

RawWaveformAnalysis::RawWaveformAnalysis(const art::InputTag& rawDigitProducerTag) :
        fRawDigitProducerTag(rawDigitProducerTag),
        fWaveName(""),
        fHitsName(""),
        fFile(0)
{
    return;
}

RawWaveformAnalysis::~RawWaveformAnalysis()
{
    return;
}

void RawWaveformAnalysis::setupAnalysis(TFile* tFile,const std::string& waveDir,const std::string& hitsDir)
{
    // Keep track of input info
    fFile     = tFile;
    fWaveName = waveDir;
    fHitsName = hitsDir;
    
    fWaveDir  = fFile->mkdir(fWaveName.c_str());
    fHitsDir  = fFile->mkdir(fHitsName.c_str());
    
    //Create a Tree to store event information
    fWaveDir->cd();
    fWaveTree = new TTree("waveforms","RawWaveformAnalysis: hists of waveforms");
    
    fNumBinsToAve = 10;
    
    fWeightVec.resize(fNumBinsToAve);
    
    for(int idx = 0; idx < fNumBinsToAve/2; idx++)
    {
        float weight = idx + 1;
        
        fWeightVec.at(idx)                     = weight;
        fWeightVec.at(fNumBinsToAve - idx - 1) = weight;
    }
    
    fWeightSum = std::accumulate(fWeightVec.begin(),fWeightVec.end(),0.);
    
    return;
}

void RawWaveformAnalysis::analyzeEvent(gallery::Event& event)
{
    //time the loop
    auto t_begin = std::chrono::high_resolution_clock::now();
    
    // We'll need our RawDigits
    auto const& rawDigitHandle = event.getValidHandle<std::vector<raw::RawDigit>>(fRawDigitProducerTag);
    auto const& rawDigitVec(*rawDigitHandle);
    
    // Define the tuple output structure
    struct WaveformObject
    {
        unsigned channel;
        unsigned plane;
        float    pedestal;
        float    rms;
        int      dataSize;
        float    waveform[6400];
    } waveformObject;
    
    // Define the tuple output structure
    struct DerivativeObject
    {
        unsigned channel;
        int      dataSize;
        float    derivative[6400];
    } derivativeObject;
    
//    DerivativeObject derivativeObjectAve;
    
    DerivativeObject erosionObject;
    DerivativeObject dilationObject;
    DerivativeObject erosionDilationAveObject;
    DerivativeObject erosionDilationDiffObject;
    
    //loop over the raw digits
    for(size_t idx = 0; idx < rawDigitVec.size(); idx++)
    {
        // get the reference to the current raw::RawDigit
        const raw::RawDigit& rawDigit = rawDigitVec.at(idx);
        
        raw::ChannelID_t channel  = rawDigit.Channel();
        size_t           dataSize = rawDigit.Samples();
        float            pedestal = rawDigit.GetPedestal();
        
        if (dataSize < 5)
        {
            std::cout << "--> Channel " << channel << " has " << dataSize << " elements" << std::endl;
            continue;
        }
        
        std::vector<short> adcValVec(dataSize);
        
        raw::Uncompress(rawDigit.ADCs(), adcValVec, rawDigit.Compression());
        
        // Do we need to truncate?
        if (dataSize > 6400)
        {
            std::copy(adcValVec.begin()+2400,adcValVec.begin()+8800,adcValVec.begin());
            adcValVec.resize(6400);
            dataSize = 6400;
        }
        
        // Get the pedestal subtracted data, centered in the deconvolution vector
        std::vector<float> rawAdcLessPedVec(dataSize);
        
        std::transform(adcValVec.begin(),adcValVec.end(),rawAdcLessPedVec.begin(),[pedestal](const short& adc){return std::round(float(adc) - pedestal);});
        
        // Well, believe it or not that doesn't actually fix the problem...
//        fixTheFreakingWaveform(rawAdcLessPedVec,rawAdcLessPedVec);
        std::vector<float> waveErosionVec;
        std::vector<float> waveDilationVec;
        std::vector<float> waveAveVec;
        std::vector<float> waveDiffVec;
        
        getErosionAndDilationVecs(rawAdcLessPedVec, waveErosionVec, waveDilationVec);
        
        waveAveVec.resize(rawAdcLessPedVec.size());
        
        std::transform(waveErosionVec.begin(),waveErosionVec.end(),waveDilationVec.begin(),waveAveVec.begin(),[](const auto& left,const auto& right){return 0.5 *(left + right);});
        
        waveDiffVec.resize(rawAdcLessPedVec.size());
        
        std::transform(waveErosionVec.begin(),waveErosionVec.end(),waveDilationVec.begin(),waveDiffVec.begin(),[](const auto& left,const auto& right){return right-left;});
        
//        std::transform(rawAdcLessPedVec.begin(),rawAdcLessPedVec.end(),waveAveVec.begin(),rawAdcLessPedVec.begin(),[](const auto& val, const auto& cor){return cor - val;});
        
        waveformObject.channel  = channel;
        waveformObject.plane    = 0;
        waveformObject.pedestal = pedestal;
        waveformObject.rms      = rawDigit.GetSigma();
        waveformObject.dataSize = dataSize;
        
        std::copy(rawAdcLessPedVec.begin(),rawAdcLessPedVec.end(),waveformObject.waveform);
        
        // Now handle the "derivative" vector...
        std::vector<float> derivVec(dataSize);
        
        // If we have a collection plane then take the derivative
        if (channel >= 4800)
        {
            derivVec[0]          = 0;
            derivVec[dataSize-1] = 0;
            
            for(size_t idx = 1; idx < rawAdcLessPedVec.size()-1; idx++)
                derivVec[idx] = 0.5 * (rawAdcLessPedVec.at(idx+1) - rawAdcLessPedVec.at(idx-1));
        }
        // Otherwise a straight copy
        else std::copy(rawAdcLessPedVec.begin(),rawAdcLessPedVec.end(),derivVec.begin());
        
        // Now smooth the differential waveform
        // This version for testing, we'll need a better version for production
        // We don't want the smoothing procedure to factor into the output
        // So start by making a local copy of the input vector
        std::vector<float> tempVec = derivVec;
        
        // Now run the "triangle" smoothing operation
        for(size_t idx = 2; idx < tempVec.size() - 2; idx++)
            derivVec[idx] = (tempVec.at(idx-2) + 2.*tempVec.at(idx-1) + 3.*tempVec.at(idx) + 2.*tempVec.at(idx+1) + tempVec.at(idx+2))/9.;

        // Get erosion and dilation vectors
//        std::vector<float> erosionVec;
//        std::vector<float> dilationVec;
        
//        getErosionAndDilationVecs(derivVec, erosionVec, dilationVec);
        
        erosionObject.channel = channel;
        erosionObject.dataSize = waveErosionVec.size();
        
        //std::copy(erosionVec.begin(),erosionVec.end(),erosionObject.derivative);
        std::copy(waveErosionVec.begin(),waveErosionVec.end(),erosionObject.derivative);
        
        dilationObject.channel = channel;
        dilationObject.dataSize = waveDilationVec.size();
        
        //std::copy(dilationVec.begin(),dilationVec.end(),dilationObject.derivative);
        std::copy(waveDilationVec.begin(),waveDilationVec.end(),dilationObject.derivative);
        
        // Now get the average
//        std::vector<float> erosionDilationAve(waveErosionVec.size());
        
//        std::transform(waveErosionVec.begin(),waveErosionVec.end(),waveDilationVec.begin(),erosionDilationAve.begin(),[](const auto& left,const auto& right){return 0.5 *(left + right);});
        
        // Baseline correct the derivative
//        std::transform(derivVec.begin(),derivVec.end(),erosionDilationAve.begin(),derivVec.begin(),[](const auto& val, const auto& cor){return val - cor;});
        
        erosionDilationAveObject.channel  = channel;
        erosionDilationAveObject.dataSize = waveAveVec.size();
        
        //std::copy(erosionDilationAve.begin(),erosionDilationAve.end(),erosionDilationAveObject.derivative);
        std::copy(waveAveVec.begin(),waveAveVec.end(),erosionDilationAveObject.derivative);
        
        // Finally, the difference
        erosionDilationDiffObject.channel = channel;
        erosionDilationDiffObject.dataSize = waveDiffVec.size();
        
        std::copy(waveDiffVec.begin(),waveDiffVec.end(),erosionDilationDiffObject.derivative);

        // Start saving stuff, beginning with the derivative vector
        derivativeObject.channel  = channel;
        derivativeObject.dataSize = dataSize;
        
        std::copy(derivVec.begin(),derivVec.end(),derivativeObject.derivative);
        
        // Set up to store in root
        std::string branchName = "Waveform_" + std::to_string(channel);
        std::string branchDef  = "channel/i:plane/i:pedestal/F:rms/F:dataSize/I:waveform[dataSize]/F"; //:derivative[dataSize]/F";
        
        fWaveDir->cd();
        
        TBranch* waveBranch = fWaveTree->Branch(branchName.c_str(),&waveformObject,branchDef.c_str());
        
        waveBranch->Fill();
        
        branchName += "_d";
        branchDef   = "channel/i:dataSize/I:derivative[dataSize]/F";
        
        TBranch* derivBranch = fWaveTree->Branch(branchName.c_str(),&derivativeObject,branchDef.c_str());
        
        derivBranch->Fill();
        
        std::string erosionName = "Erosion_" + std::to_string(channel);
        branchDef   = "channel/i:dataSize/I:erosion[dataSize]/F";
        
        TBranch* erosionBranch = fWaveTree->Branch(erosionName.c_str(),&erosionObject,branchDef.c_str());
        
        erosionBranch->Fill();
        
        std::string dilationName = "Dilation_" + std::to_string(channel);
        branchDef   = "channel/i:dataSize/I:dilation[dataSize]/F";
        
        TBranch* dilationBranch = fWaveTree->Branch(dilationName.c_str(),&dilationObject,branchDef.c_str());
        
        dilationBranch->Fill();
        
        std::string edaveName = "EDAve_" + std::to_string(channel);
        branchDef   = "channel/i:dataSize/I:edave[dataSize]/F";
        
        TBranch* edaveBranch = fWaveTree->Branch(edaveName.c_str(),&erosionDilationAveObject,branchDef.c_str());
        
        edaveBranch->Fill();
        
        std::string eddiffName = "EDDiff_" + std::to_string(channel);
        branchDef   = "channel/i:dataSize/I:eddiff[dataSize]/F";
        
        TBranch* eddiffBranch = fWaveTree->Branch(eddiffName.c_str(),&erosionDilationDiffObject,branchDef.c_str());
        
        eddiffBranch->Fill();
    }//end first loop over hits.
    
    fWaveTree->Fill();
    
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> time_total_ms(t_end-t_begin);
    std::cout << "\tEvent took " << time_total_ms.count() << " ms to process." << std::endl;
    
    return;
}
    
void  RawWaveformAnalysis::getErosionAndDilationVecs(const std::vector<float>& inputWaveform,
                                                     std::vector<float>&       erosionVec,
                                                     std::vector<float>&       dilationVec) const
{
    // Set the window size
    int halfWindowSize(30);
    
    // Define the "smallest" function
//    auto smaller = [](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);};
    
    // Initialize min and max elements
    std::pair<std::vector<float>::const_iterator,std::vector<float>::const_iterator> minMaxItr =
        std::minmax_element(inputWaveform.begin(),inputWaveform.begin()+halfWindowSize);
    
    std::vector<float>::const_iterator minElementItr = minMaxItr.first;
    std::vector<float>::const_iterator maxElementItr = minMaxItr.second;
    
    // Initialize the erosion and dilation vectors
    erosionVec.resize(inputWaveform.size());
    dilationVec.resize(inputWaveform.size());
    
    // Now loop through remaining elements and complete the vectors
    std::vector<float>::iterator minItr = erosionVec.begin();
    std::vector<float>::iterator maxItr = dilationVec.begin();
    
    for (std::vector<float>::const_iterator inputItr = inputWaveform.begin(); inputItr != inputWaveform.end(); inputItr++)
    {
        // There are two conditions to check:
        // 1) is the current min/max element outside the current window?
        // 2) is the new element smaller/larger than the current min/max?
        
        // Make sure we are not running off the end of the vector
        if (std::distance(inputItr,inputWaveform.end()) > halfWindowSize)
        {
            if (std::distance(minElementItr,inputItr) >= halfWindowSize)
                minElementItr = std::min_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) < *minElementItr)
                minElementItr = inputItr + halfWindowSize;
            
            if (std::distance(maxElementItr,inputItr) >= halfWindowSize)
                maxElementItr = std::max_element(inputItr - halfWindowSize + 1, inputItr + halfWindowSize + 1);
            else if (*(inputItr + halfWindowSize) > *maxElementItr)
                maxElementItr = inputItr + halfWindowSize;
        }
        
        // Update the vectors
        *minItr++ = *minElementItr;
        *maxItr++ = *maxElementItr;
    }
    
    return;
}

    
void RawWaveformAnalysis::averageInputWaveform(const std::vector<float>& inputWaveform, std::vector<float>& outputWaveform) const
{
    // Vector reduction - take the 10 bin average
    float aveSum(0.);
    
    outputWaveform.resize(inputWaveform.size()/fNumBinsToAve);
    
    for(size_t idx = 0; idx < inputWaveform.size(); idx++)
    {
        aveSum += fWeightVec.at(idx % fNumBinsToAve) * inputWaveform.at(idx);
        
        if ((idx + 1) % fNumBinsToAve == 0)
        {
            outputWaveform[idx/fNumBinsToAve] = aveSum / fWeightSum;
            
            aveSum = 0.;
        }
    }
    
    return;
}
    
float RawWaveformAnalysis::fixTheFreakingWaveform(const std::vector<float>& waveform, std::vector<float>& fixedWaveform) const
{
    // do rms calculation - the old fashioned way and over all adc values
    std::vector<float> locWaveform = waveform;
    
    // sort in ascending order so we can truncate the sume
    std::sort(locWaveform.begin(), locWaveform.end(),[](const auto& left, const auto& right){return std::fabs(left) < std::fabs(right);});
    
    float localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + locWaveform.size()/2, locWaveform.begin(), 0.);
    
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(locWaveform.size()/2)));
    
    float threshold = 6. * localRMS;
    
    std::vector<float>::iterator threshItr = std::find_if(locWaveform.begin(),locWaveform.end(),[threshold](const auto& val){return std::fabs(val) > threshold;});
    
    int minNumBins = std::max(int(0.6 * locWaveform.size()),int(std::distance(locWaveform.begin(),threshItr)));
    
    // Get the truncated sum
    localRMS = std::inner_product(locWaveform.begin(), locWaveform.begin() + minNumBins, locWaveform.begin(), 0.);
    localRMS = std::sqrt(std::max(float(0.),localRMS / float(minNumBins)));
    
    // Now get the average value
    float aveSum      = std::accumulate(locWaveform.begin(), locWaveform.begin() + minNumBins, 0.);
    float newPedestal = aveSum / minNumBins;
    
    std::transform(waveform.begin(), waveform.end(), fixedWaveform.begin(), [newPedestal](const auto& val){return val - newPedestal;});
    
    return localRMS;
}

} // namespace
