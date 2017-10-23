/*************************************************************
 * 
 * NeutrinoIDFilterPastt program
 * 
 * This is a simple check for whether NeutrinoIDFilter led to 
 * a pass or fail.
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
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"

//our own includes!
#include "hist_utilities.h"

//function declarations
void print_usage();
std::vector<std::string> create_file_list(char*);

void print_usage(){
  std::cout << "Usage:"
	    << "\n\tNeutrinoIDFilterPass <input_file_list> <input_tag>"
	    << "\n\n\tNeutrinoIDFilterPass will read in a list of larsoft files, and "
	    << "\n\tsee if each event passes or fails the filter."
	    << "\n\t<input_tag> is of format module:instance:process. You can leave off the last two."
	    << std::endl;
}

std::vector<std::string> create_file_list(char* fin_list){

  std::vector<std::string> filenames;
  std::string file_name;

  std::ifstream input_file(fin_list);
  while(getline(input_file,file_name))
    filenames.push_back(file_name);

  return filenames;
}

//normal C++ main program
int main(int argc, char** argv) {

  if(argc!=3){
    print_usage();
    return -1;
  }

  //hande inputs and outputs
  std::vector<std::string> filenames = create_file_list(argv[1]);
  art::InputTag input_tag(argv[2]);

  size_t n_events=0,n_pass=0,n_fail=0;

  bool pass;
  //event loop...
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {

    try{
      auto const& assoc_handle = ev.getValidHandle< art::Assns<recob::Vertex,recob::Track> >(input_tag);
      if(assoc_handle->size() > 0)
	pass = true;
      else pass=false;
    }
    catch(...){
      pass=false;
    }
    
    //to get run and event info, you use this "eventAuxillary()" object.
    std::cout << "Run " << ev.eventAuxiliary().run() << ", "
	      << "Event " << ev.eventAuxiliary().event()
	      << ": Pass? " << pass << std::endl;

    n_events++;
    if(pass) n_pass++;
    else n_fail++;
    
    
  } //end loop over events!

  std::cout << "Total of " << n_pass << " (" << n_fail << ") / " << n_events << " pass (fail)." << std::endl;
}
