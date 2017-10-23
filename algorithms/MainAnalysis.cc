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
#include "BasicHitAnalysis.h"

//function declarations
void print_usage();
std::vector<std::string> create_file_list(char*);

void print_usage(){
    std::cout << "Usage:"
    << "\n\tMainAnalysis <input_file_list> <hit_input_tag> <output_file_name>"
    << "\n\n\tMainAnalysis will read in a list of larsoft files, and "
    << "\n\tcreate a single ROOT output file with hit and wire info."
    << "\n\t<hit_input_tag> is of format module:instance:process. You can leave off the last two."
    << std::endl;
}

std::vector<std::string> create_file_list(char* fin_list){
    
    std::vector<std::string> filenames;
    std::string file_name;
    
    std::ifstream input_file(fin_list);
    while(getline(input_file,file_name)) filenames.push_back(file_name);
    
    return filenames;
}

//normal C++ main program
int main(int argc, char** argv) {
    
    if(argc!=4){
        print_usage();
        return -1;
    }

    //hande inputs and outputs
    std::vector<std::string> filenames = create_file_list(argv[1]);
    art::InputTag hit_tag(argv[2]);
    TFile f_output(argv[3],"RECREATE");
    
    art::InputTag cal_tag("caldata::");
    
    HitAnalysis::BasicHitAnalysis basicHitAnalysis(hit_tag,cal_tag);
    
    basicHitAnalysis.setupAnalysis(f_output, "HitHistograms", "HitTuple");
    
    //event loop...
    for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next())
    {
        basicHitAnalysis.analyzeEvent(ev);
    } //end loop over events!
    
    
    //and ... write to file!
    f_output.Write();
    f_output.Close();
    
}
