/*************************************************************
 *
 * This program is meant to drive looking at an event in a single
 * file. 
 *
 * Derived from examples by Wes Ketchum
 *
 *************************************************************/


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
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"

//our own includes!
#include "hist_utilities.h"

// Start the conversion
//#include "WaveformAnalysis.h"
#include "RawWaveformAnalysis.h"

//function declarations
void print_usage();
std::vector<std::string> create_file_list(char*);

void print_usage(){
    std::cout << "Usage:"
    << "\n\tHitAna <input_file_list> <hit_input_tag> <output_file_name>"
    << "\n\n\tHitAna will read in a list of larsoft files, and "
    << "\n\tcreate a single ROOT output file with hit and wire info."
    << "\n\t<hit_input_tag> is of format module:instance:process. You can leave off the last two."
    << std::endl;
}

//normal C++ main program
int main(int argc, char** argv)
{
    if(argc!=4){
        print_usage();
        return -1;
}
    
    //hande inputs and outputs
    std::vector<std::string> filenameVec = {argv[1]};
    art::InputTag hit_tag(argv[2]);
    TFile f_output(argv[3],"RECREATE");
    
    art::InputTag rawDigit_tag(argv[2]);
    
    //Waveform::WaveformAnalysis waveformAnalysis(hit_tag, rawDigit_tag);
    Waveform::RawWaveformAnalysis waveformAnalysis(rawDigit_tag);
    
    std::string waveform = "waveforms";
    std::string hits     = "hits";
    
    waveformAnalysis.setupAnalysis(&f_output, waveform, hits);
    
    //event loop...
    for (gallery::Event ev(filenameVec) ; !ev.atEnd(); ev.next())
    {
        waveformAnalysis.analyzeEvent(ev);
    } //end loop over events!
    
    //and ... write to file!
    f_output.Write();
    f_output.Close();
    
}
