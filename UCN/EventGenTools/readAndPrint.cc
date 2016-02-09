#include	 <iostream>
#include	 <fstream>
#include	 <TGaxis.h>
#include	 <sstream>
#include	 <TGraph.h>
#include	 <TGraphErrors.h>
#include	 <TCanvas.h>
#include	 <TApplication.h>
#include	 <stdlib.h>
#include	 <TF1.h>
#include	 <TH1.h>
#include	 <TProfile.h>
#include	 <TObjArray.h>
#include	 <TStyle.h>
#include	 <TMarker.h>
#include	 <math.h>
#include	 <TStyle.h>
#include	 <TPaveStats.h>
#include	 <TPaveText.h>
#include	 <vector>
#include	 <string.h>
#include	 <fstream>
#include	 <TROOT.h>
#include	 <TFile.h>
#include	 <TLegend.h>
#include         <TLegendEntry.h>
#include	 <time.h>
#include	 <TH2F.h>
#include         <assert.h>
#include	 <string>
#include	 <TRandom.h>
#include	 <TChain.h>
using		 namespace std;

#define		INPUT_DATA_FILE			"n1_1000000/Evts_0.root"
#define		OUTPUT_FILE			"big_initPtclInfo.txt"

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

enum DecayType		// this lists the various integer assignments
{
  D_GAMMA = 22,
  D_ELECTRON = 11,
  D_POSITRON = -11,
  D_NEUTRINO = -12, // anti-nu_e
  D_NONEVENT = 0
};

int main()
{
  // opens the file listed in INPUT_DATA_FILE
  TFile* fRead = TFile::Open(INPUT_DATA_FILE);
  // extracts the TTree object inside. I've checked and it is called Evts
  TTree* tree = (TTree*)fRead -> Get("Evts");

  int nEntries = tree -> GetEntries();

  int eid = 0;	// these are the variable types found in EventTreeScanner
  double initE = 0;
  double p[3];
  double x[3];
  int decayFlag = 0;
  double t = 0;
  double w = 1.0;

  tree -> SetBranchAddress("num", &eid);
  tree -> SetBranchAddress("PID", &decayFlag);
  tree -> SetBranchAddress("KE", &initE);
  tree -> SetBranchAddress("vertex", &x);
  tree -> SetBranchAddress("direction", &p);
  tree -> SetBranchAddress("time", &t);
  tree -> SetBranchAddress("weight", &w);

  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);

  for(int i = 0; i < nEntries; i++)
  {
    tree -> GetEntry(i);
    outfile	<< eid << "\t"
	    	<< decayFlag << "\t"
	    	<< initE << "\t"
	    	<< x[0] << "\t" << x[1] << "\t" << x[2] << "\t"
		<< p[0] << "\t" << p[1] << "\t" << p[2] << "\t"
		<< t << "\t"
		<< w << "\n";

  }

  outfile.close();

  cout << "-------------- End of Program ---------------" << endl;
  return 0;
}
