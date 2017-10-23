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

//TTree Structs (could make ntuples, but less extendible...)
typedef struct evTreeObj
{
    unsigned int run;
    unsigned int ev;
    unsigned int nhits_u;
    unsigned int nhits_v;
    unsigned int nhits_w;
    
    float        hit_amp_mean_u;
    float        hit_amp_rms_u;
    float        hit_amp_mean_v;
    float        hit_amp_rms_v;
    float        hit_amp_mean_w;
    float        hit_amp_rms_w;
    
    float        hit_integral_mean_u;
    float        hit_integral_rms_u;
    float        hit_integral_mean_v;
    float        hit_integral_rms_v;
    float        hit_integral_mean_w;
    float        hit_integral_rms_w;
    
    float        hit_sumadc_mean_u;
    float        hit_sumadc_rms_u;
    float        hit_sumadc_mean_v;
    float        hit_sumadc_rms_v;
    float        hit_sumadc_mean_w;
    float        hit_sumadc_rms_w;
    
    float        hit_mult_mean_u;
    float        hit_mult_rms_u;
    float        hit_mult_mean_v;
    float        hit_mult_rms_v;
    float        hit_mult_mean_w;
    float        hit_mult_rms_w;
    
    std::string BranchString(){
        std::string str("run/i:ev/i:nhits_u/i:nhits_v/i:nhits_w/i:hit_amp_mean_u/F:hit_amp_rms_u/F:hit_amp_mean_v/F:hit_amp_rms_v/F:hit_amp_mean_w/F:hit_amp_rms_w/F:hit_integral_mean_u/F:hit_integral_rms_u/F:hit_integral_mean_v/F:hit_integral_rms_v/F:hit_integral_mean_w/F:hit_integral_rms_w/F:hit_sumadc_mean_u/F:hit_sumadc_rms_u/F:hit_sumadc_mean_v/F:hit_sumadc_rms_v/F:hit_sumadc_mean_w/F:hit_sumadc_rms_w/F:hit_mult_mean_u/F:hit_mult_rms_u/F:hit_mult_mean_v/F:hit_mult_rms_v/F:hit_mult_mean_w/F:hit_mult_rms_w/F");
        return str;
    }
    
    void Init()
    {
        run=0; ev=0;
        nhits_u=0; nhits_v=0; nhits_w=0;
        
        hit_amp_mean_u=0; hit_amp_rms_u=0;
        hit_amp_mean_v=0; hit_amp_rms_v=0;
        hit_amp_mean_w=0; hit_amp_rms_w=0;
        
        hit_integral_mean_u=0; hit_integral_rms_u=0;
        hit_integral_mean_v=0; hit_integral_rms_v=0;
        hit_integral_mean_w=0; hit_integral_rms_w=0;
        
        hit_sumadc_mean_u=0; hit_sumadc_rms_u=0;
        hit_sumadc_mean_v=0; hit_sumadc_rms_v=0;
        hit_sumadc_mean_w=0; hit_sumadc_rms_w=0;
        
        hit_mult_mean_u=0; hit_mult_rms_u=0;
        hit_mult_mean_v=0; hit_mult_rms_v=0;
        hit_mult_mean_w=0; hit_mult_rms_w=0;
    }
    
} evTreeObj_t;


typedef struct hitTreeObj
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
} hitTreeObj_t;



void print_usage(){
    std::cout << "Usage:"
    << "\n\tHitAna <input_file_list> <hit_input_tag> <output_file_name>"
    << "\n\n\tHitAna will read in a list of larsoft files, and "
    << "\n\tcreate a single ROOT output file with hit and wire info."
    << "\n\t<hit_input_tag> is of format module:instance:process. You can leave off the last two."
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
int main(int argc, char** argv)
{
    if(argc!=4)
    {
        print_usage();
        return -1;
    }
    
    //hande inputs and outputs
    std::vector<std::string> filenames = create_file_list(argv[1]);
    art::InputTag hit_tag(argv[2]);
    TFile f_output(argv[3],"RECREATE");
    
    TH1F numHitsHist("NumHits", "# hits", 300, 0., 30000.);
    
    //Create a Tree to store event information
    evTreeObj_t evobj;
    TTree tree_ev_hit("tree_ev_hit","HitAna: Hits per event tree");
    tree_ev_hit.Branch("",&evobj,evobj.BranchString().c_str());
    
    //Create a Tree to store hit information
    hitTreeObj_t hitobj;
    TTree tree_hit("tree_hit","HitAna: Hit tree");
    tree_hit.Branch("hit",&hitobj,hitobj.BranchString().c_str());
    
    //event loop...
    for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next())
    {
        //clear out event tree object
        evobj.Init();
        
        //time the loop
        auto t_begin = std::chrono::high_resolution_clock::now();
        
        //to get run and event info, you use this "eventAuxillary()" object.
        std::cout << "Processing "
        << "Run " << ev.eventAuxiliary().run() << ", "
        << "Event " << ev.eventAuxiliary().event() << std::endl;
        
        hitobj.run = ev.eventAuxiliary().run();
        hitobj.ev  = ev.eventAuxiliary().event();
        evobj.run  = ev.eventAuxiliary().run();
        evobj.ev   = ev.eventAuxiliary().event();
        
        //get handle to hits
        auto const& hit_handle = ev.getValidHandle<std::vector<recob::Hit>>(hit_tag);
        auto const& hit_vec(*hit_handle);
        
        //get the associations to the recob::Wire objects
        art::FindOne<recob::Wire> wire_per_hit(hit_handle,ev,hit_tag);
        
        //For good measure, print out the number of hits
        std::cout << "\tThere are " << hit_vec.size() << " hits in this event." << std::endl;
        
        numHitsHist.Fill(hit_vec.size(), 1.);
        
        //loop over the hits
        //for(auto const& hit : hit_vec){
        for(size_t i_h=0; i_h < hit_vec.size(); ++i_h)
        {
            auto const& hit = hit_vec[i_h];
            
            //get associated wire object
            auto const& wire = (wire_per_hit.at(i_h)).ref();
            auto const& roi  = wire.SignalROI().find_range((int)(hit.PeakTime()));
            
            /*
             //Note: The association to the wire object gives us the full wire.
             //      One (i.e. you, the reader =p )still needs to find the ROI
             //      that contains the hit. You can do that using the "find_range"
             //      utility of the sparse vector. See below!
             auto const& myrange = wire.SignalROI().find_range(hit.StartTick());
             std::cout << "\tHit with time " << hit.PeakTime()
             << " found in range ["
             << myrange.begin_index() << ","
             << myrange.end_index() << ")"
             << std::endl;
             auto const& roi_wvfm = myrange.data();
             for(size_t i_w=0; i_w<roi_wvfm.size(); ++i_w)
             std::cout << " tick " << i_w << " in roi"
	            << " or tick " i_w+myrange.begin_index() << " in full waveform"
             << " has " << roi_wvfm[i_w] << " adc." << std::endl;
             */
            
            //do the event tree info
            if(hit.View()==geo::View_t::kU)
            {
                evobj.nhits_u++;
                evobj.hit_amp_mean_u      += hit.PeakAmplitude();
                evobj.hit_integral_mean_u += hit.Integral();
                evobj.hit_sumadc_mean_u   += hit.SummedADC();
                evobj.hit_mult_mean_u     += hit.Multiplicity();
            }
            else if(hit.View()==geo::View_t::kV)
            {
                evobj.nhits_v++;
                evobj.hit_amp_mean_v      += hit.PeakAmplitude();
                evobj.hit_integral_mean_v += hit.Integral();
                evobj.hit_sumadc_mean_v   += hit.SummedADC();
                evobj.hit_mult_mean_v     += hit.Multiplicity();
            }
            else if(hit.View()==geo::View_t::kW)
            {
                evobj.nhits_w++;
                evobj.hit_amp_mean_w      += hit.PeakAmplitude();
                evobj.hit_integral_mean_w += hit.Integral();
                evobj.hit_sumadc_mean_w   += hit.SummedADC();
                evobj.hit_mult_mean_w     += hit.Multiplicity();
            }
            
            //now fill the hit info
            hitobj.start_tick   = hit.StartTick();
            hitobj.end_tick     = hit.EndTick();
            hitobj.roi_start    = roi.begin_index();
            hitobj.roi_end      = roi.end_index();
            hitobj.roi_size     = roi.size();
            hitobj.ch           = hit.Channel();
            hitobj.time         = hit.PeakTime();
            hitobj.time_err     = hit.SigmaPeakTime();
            hitobj.rms          = hit.RMS();
            hitobj.amp          = hit.PeakAmplitude();
            hitobj.amp_err      = hit.SigmaPeakAmplitude();
            hitobj.sumadc       = hit.SummedADC();
            hitobj.integral     = hit.Integral();
            hitobj.integral_err = hit.SigmaIntegral();
            hitobj.gof          = hit.GoodnessOfFit();
            hitobj.mult         = hit.Multiplicity();
            hitobj.idx          = hit.LocalIndex();
            hitobj.view         = hit.View();
            
            tree_hit.Fill();
        }//end first loop over hits.
        
        
        //mean calculations...
        evobj.hit_amp_mean_u      = evobj.hit_amp_mean_u / evobj.nhits_u;
        evobj.hit_amp_mean_v      = evobj.hit_amp_mean_v / evobj.nhits_v;
        evobj.hit_amp_mean_w      = evobj.hit_amp_mean_w / evobj.nhits_w;
        
        evobj.hit_integral_mean_u = evobj.hit_integral_mean_u / evobj.nhits_u;
        evobj.hit_integral_mean_v = evobj.hit_integral_mean_v / evobj.nhits_v;
        evobj.hit_integral_mean_w = evobj.hit_integral_mean_w / evobj.nhits_w;
        
        evobj.hit_sumadc_mean_u   = evobj.hit_sumadc_mean_u / evobj.nhits_u;
        evobj.hit_sumadc_mean_v   = evobj.hit_sumadc_mean_v / evobj.nhits_v;
        evobj.hit_sumadc_mean_w   = evobj.hit_sumadc_mean_w / evobj.nhits_w;
        
        evobj.hit_mult_mean_u     = evobj.hit_mult_mean_u / evobj.nhits_u;
        evobj.hit_mult_mean_v     = evobj.hit_mult_mean_v / evobj.nhits_v;
        evobj.hit_mult_mean_w     = evobj.hit_mult_mean_w / evobj.nhits_w;
        
        //second loop over hits to do rms calculations...
        for(auto const& hit : hit_vec)
        {
            if(hit.View()==geo::View_t::kU){
                evobj.hit_amp_rms_u      += (hit.PeakAmplitude()-evobj.hit_amp_mean_u)*(hit.PeakAmplitude()-evobj.hit_amp_mean_u);
                evobj.hit_integral_rms_u += (hit.Integral()-evobj.hit_integral_mean_u)*(hit.Integral()-evobj.hit_integral_mean_u);
                evobj.hit_sumadc_rms_u   += (hit.SummedADC()-evobj.hit_sumadc_mean_u)*(hit.SummedADC()-evobj.hit_sumadc_mean_u);
                evobj.hit_mult_rms_u     += (hit.Multiplicity()-evobj.hit_mult_mean_u)*(hit.Multiplicity()-evobj.hit_mult_mean_u);
            }
            else if(hit.View()==geo::View_t::kV){
                evobj.hit_amp_rms_v      += (hit.PeakAmplitude()-evobj.hit_amp_mean_v)*(hit.PeakAmplitude()-evobj.hit_amp_mean_v);
                evobj.hit_integral_rms_v += (hit.Integral()-evobj.hit_integral_mean_v)*(hit.Integral()-evobj.hit_integral_mean_v);
                evobj.hit_sumadc_rms_v   += (hit.SummedADC()-evobj.hit_sumadc_mean_v)*(hit.SummedADC()-evobj.hit_sumadc_mean_v);
                evobj.hit_mult_rms_v     += (hit.Multiplicity()-evobj.hit_mult_mean_v)*(hit.Multiplicity()-evobj.hit_mult_mean_v);
            }
            else if(hit.View()==geo::View_t::kW){
                evobj.hit_amp_rms_w      += (hit.PeakAmplitude()-evobj.hit_amp_mean_w)*(hit.PeakAmplitude()-evobj.hit_amp_mean_w);
                evobj.hit_integral_rms_w += (hit.Integral()-evobj.hit_integral_mean_w)*(hit.Integral()-evobj.hit_integral_mean_w);
                evobj.hit_sumadc_rms_w   += (hit.SummedADC()-evobj.hit_sumadc_mean_w)*(hit.SummedADC()-evobj.hit_sumadc_mean_w);
                evobj.hit_mult_rms_w     += (hit.Multiplicity()-evobj.hit_mult_mean_w)*(hit.Multiplicity()-evobj.hit_mult_mean_w);
            }
            
        }//end second loop over hits;
        
        //rms calculations...
        evobj.hit_amp_rms_u      = std::sqrt(evobj.hit_amp_rms_u/(evobj.nhits_u-1));
        evobj.hit_integral_rms_u = std::sqrt(evobj.hit_integral_rms_u/(evobj.nhits_u-1));
        evobj.hit_sumadc_rms_u   = std::sqrt(evobj.hit_sumadc_rms_u/(evobj.nhits_u-1));
        evobj.hit_mult_rms_u     = std::sqrt(evobj.hit_mult_rms_u/(evobj.nhits_u-1));
        
        evobj.hit_amp_rms_v      = std::sqrt(evobj.hit_amp_rms_v/(evobj.nhits_v-1));
        evobj.hit_integral_rms_v = std::sqrt(evobj.hit_integral_rms_v/(evobj.nhits_v-1));
        evobj.hit_sumadc_rms_v   = std::sqrt(evobj.hit_sumadc_rms_v/(evobj.nhits_v-1));
        evobj.hit_mult_rms_v     = std::sqrt(evobj.hit_mult_rms_v/(evobj.nhits_v-1));
        
        evobj.hit_amp_rms_w      = std::sqrt(evobj.hit_amp_rms_w/(evobj.nhits_w-1));
        evobj.hit_integral_rms_w = std::sqrt(evobj.hit_integral_rms_w/(evobj.nhits_w-1));
        evobj.hit_sumadc_rms_w   = std::sqrt(evobj.hit_sumadc_rms_w/(evobj.nhits_w-1));
        evobj.hit_mult_rms_w     = std::sqrt(evobj.hit_mult_rms_w/(evobj.nhits_w-1));
        
        tree_ev_hit.Fill();
        
        auto t_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double,std::milli> time_total_ms(t_end-t_begin);
        std::cout << "\tEvent took " << time_total_ms.count() << " ms to process." << std::endl;
    } //end loop over events!
    
    //and ... write to file!
    f_output.Write();
    f_output.Close();
    
}
