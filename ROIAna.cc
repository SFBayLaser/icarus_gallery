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

//function declarations
void print_usage();
std::vector<std::string> create_file_list(char*);

//TTree Structs (could make ntuples, but less extendible...)
typedef struct roiTreeObj{
  unsigned int run;
  unsigned int ev;
  unsigned int ch;
  unsigned int start_tick;
  unsigned int end_tick;
  unsigned int size;

  int max_adc;
  unsigned int max_tick;
  int min_adc;
  unsigned int min_tick;

  int integral;

  int integral_post_max;
  int integral_post_max5;
  int integral_post_max10;
  
  int integral_pre_max;
  int integral_pre_max5;
  int integral_pre_max10;

  int integral_first10;
  int integral_first20;
  int integral_last10;
  int integral_last20;

  int min_premax_adc;
  unsigned int min_premax_tick;

  unsigned int nhits;
  unsigned int nhits_max10;
  
  //function to return branch list"
  std::string BranchString(){
    std::string str("run/i:ev/i:ch/i:start_tick/i:end_tick/i:size/i:max_adc/I:max_tick/i:min_adc/I:min_tick/i:integral/I:integral_post_max/I:integral_post_max5/I:integral_post_max10/I:integral_pre_max/I:integral_pre_max5/I:integral_pre_max10/I:integral_first10/I:integral_first20/I:integral_last10/I:integral_last20/I:min_premax_adc/I:min_premax_tick/i:nhits/i:nhits_max10/i");
    return str;
  }
  void Init(){
    run=999999; ev=999999; ch=999999;
    start_tick=999999; end_tick=999999; size=999999;
    max_adc=-999999; max_tick=999999; min_adc=999999; min_tick=999999;
    integral=0;
    integral_post_max=0; integral_post_max5=0; integral_post_max10=0;
    integral_pre_max=0; integral_pre_max5=0; integral_pre_max10=0;
    integral_first10=0; integral_first20=0;
    integral_last10=0; integral_last20=0;
    min_premax_adc=999999; min_premax_tick=999999;
    nhits=0; nhits_max10=0;
  }
  
} roiTreeObj_t;



void print_usage(){
  std::cout << "Usage:"
	    << "\n\tROIAna <input_file_list> <wire_input_tag> <hit_input_tag> <output_file_name>"
	    << "\n\n\tROIAna will read in a list of larsoft files, and "
	    << "\n\tcreate a single ROOT output file with roi info."
	    << "\n\t the <input_tag>s are of format module:instance:process. You can leave off the last two."
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

  if(argc!=5){
    print_usage();
    return -1;
  }

  //hande inputs and outputs
  std::vector<std::string> filenames = create_file_list(argv[1]);
  art::InputTag wire_tag(argv[2]);
  art::InputTag hit_tag(argv[3]);
  TFile f_output(argv[4],"RECREATE");

  //Create a Tree to store hit information
  roiTreeObj_t roiobj;
  TTree tree_roi("tree_roi","ROIAna: ROI tree");
  tree_roi.Branch("roi",&roiobj,roiobj.BranchString().c_str());
  
  //event loop...
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {

    //time the loop
    auto t_begin = std::chrono::high_resolution_clock::now();
    
    //to get run and event info, you use this "eventAuxillary()" object.
    std::cout << "Processing "
	 << "Run " << ev.eventAuxiliary().run() << ", "
	 << "Event " << ev.eventAuxiliary().event() << std::endl;

    //get handle to wire
    auto const& wire_handle = ev.getValidHandle<std::vector<recob::Wire>>(wire_tag);
    auto const& wire_vec(*wire_handle);

    //get the associations to the recob::Wire objects
    art::FindMany<recob::Hit> hits_per_wire(wire_handle,ev,hit_tag);
    std::vector<recob::Hit const*> hitptrs_vec; 
    
    //loop over the wires
    for(size_t i_w=0; i_w<wire_vec.size(); ++i_w){

      auto const& wire = wire_vec[i_w];
      
      hitptrs_vec.clear();
      hits_per_wire.get(i_w,hitptrs_vec); 
      
      for(auto const& roi : wire.SignalROI().get_ranges()){

	roiobj.Init();
	roiobj.run = ev.eventAuxiliary().run();
	roiobj.ev = ev.eventAuxiliary().event();
	roiobj.ch = wire.Channel();
	roiobj.start_tick = roi.begin_index();
	roiobj.end_tick = roi.end_index();
	roiobj.size = roi.size();

	auto const& wvfm = roi.data();
	for(size_t i_t=0; i_t<wvfm.size(); ++i_t){
	  roiobj.integral += wvfm[i_t];
	  
	  if(i_t<10) roiobj.integral_first10 += wvfm[i_t];
	  if(i_t<20) roiobj.integral_first20 += wvfm[i_t];
	  if(wvfm.size()>10){
	    if(wvfm.size()-i_t<10) roiobj.integral_last10 += wvfm[i_t];
	    if(wvfm.size()>20){
	      if(wvfm.size()-i_t<20) roiobj.integral_last20 += wvfm[i_t];
	    }
	  }

	  if(wvfm[i_t]>roiobj.max_adc){
	    roiobj.max_adc = wvfm[i_t];
	    roiobj.max_tick=i_t;
	  }
	  if(wvfm[i_t]<roiobj.min_adc){
	    roiobj.min_adc = wvfm[i_t];
	    roiobj.min_tick=i_t;
	  }
	  
	}

	//go back now that we have max...
	for(size_t i_t=0; i_t<wvfm.size(); ++i_t){
	  if(i_t<roiobj.max_tick) roiobj.integral_pre_max+= wvfm[i_t];
	  if(roiobj.max_tick>5) { if(i_t<roiobj.max_tick-5) roiobj.integral_pre_max5+= wvfm[i_t]; }
	  if(roiobj.max_tick>10) { if(i_t<roiobj.max_tick-10) roiobj.integral_pre_max10+= wvfm[i_t]; }
	  if(i_t>roiobj.max_tick) roiobj.integral_post_max+= wvfm[i_t];
	  if(i_t>roiobj.max_tick+5) roiobj.integral_post_max5+= wvfm[i_t];
	  if(i_t>roiobj.max_tick+10) roiobj.integral_post_max10+= wvfm[i_t];

	  if(wvfm[i_t]<roiobj.min_premax_adc){
	    roiobj.min_premax_adc=wvfm[i_t];
	    roiobj.min_premax_tick=i_t;
	  }
	}

	//readjust ticks to be for full evfm time
	roiobj.max_tick += roiobj.start_tick;
	roiobj.min_tick += roiobj.start_tick;
	roiobj.min_premax_tick += roiobj.start_tick;
	
	for(auto hitptr : hitptrs_vec){
	  if(roi.includes((int)(hitptr->PeakTime()))){
	    roiobj.nhits++;
	    if(std::abs(hitptr->PeakTime()-(float)roiobj.max_tick)<10) roiobj.nhits_max10++;
	  }//end if hit in roi
	  
	}//end loop over hits associated to wire

	tree_roi.Fill();
	
      }//end loop over rois in wire

    }//end loop over wires
      

    
    auto t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double,std::milli> time_total_ms(t_end-t_begin);
    std::cout << "\tEvent took " << time_total_ms.count() << " ms to process." << std::endl;
  } //end loop over events!


  //and ... write to file!
  f_output.Write();
  f_output.Close();

}
