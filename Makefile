#Makefile for gallery c++ programs.
#Note, being all-incllusive here: you can cut out libraries/includes you do not need
#you can also change the flags if you want too (Werror, pedantic, etc.)

CPPFLAGS=-I $(BOOST_INC) \
         -I $(CANVAS_INC) \
         -I $(CETLIB_INC) \
         -I $(CETLIB_EXCEPT_INC) \
         -I $(FHICLCPP_INC) \
         -I $(GALLERY_INC) \
         -I $(LARCOREOBJ_INC) \
         -I $(LARDATAOBJ_INC) \
         -I $(NUSIMDATA_INC) \
         -I $(ROOT_INC)

CXXFLAGS=-std=c++14 -Wall -Werror -pedantic -g -gdwarf-2
CXX=g++
LDFLAGS=$$(root-config --libs) \
        -L $(CANVAS_LIB) -l canvas \
        -L $(CETLIB_LIB) -l cetlib \
        -L $(CETLIB_EXCEPT_LIB) -l cetlib_except \
        -L $(GALLERY_LIB) -l gallery \
        -L $(NUSIMDATA_LIB) -l nusimdata_SimulationBase \
        -L $(LARCOREOBJ_LIB) -l larcoreobj_SummaryData \
        -L $(LARDATAOBJ_LIB) -l lardataobj_RecoBase -l lardataobj_MCBase -l lardataobj_RawData -l lardataobj_OpticalDetectorData -l lardataobj_AnalysisBase

LDANALYSIS=-L . -l AnalysisObjects

# lets see if we can deal with the library first
SRCS      = RawWaveformAnalysis.cxx WaveformAnalysis.cxx BasicHitAnalysis.cxx
OBJS      = $(SRCS:.cxx=.o)
TARGETLIB = libAnalysisObjects.so

%.o: %.cxx
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -fPIC -c $^ -o $@

$(TARGETLIB): $(OBJS)
	@$(CXX) -shared $^ $(LDFLAGS) -o $@

# Below is meant to be the programs built here
HitAna: HitAna.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

MainAnalysis: MainAnalysis.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDANALYSIS) $(LDFLAGS) -o $@ $<

SingleAnalysis: SingleAnalysis.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDANALYSIS) $(LDFLAGS) -o $@ $<

ROIAna: ROIAna.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

NeutrinoIDFilterPass: NeutrinoIDFilterPass.cc
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $<

all: $(TARGETLIB) HitAna MainAnalysis SingleAnalysis NeutrinoIDFilterPass ROIAna

clean:
	rm -f $(PROGRAMS) *.o
