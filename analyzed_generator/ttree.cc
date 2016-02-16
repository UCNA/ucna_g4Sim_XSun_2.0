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
#include	 <TTree.h>
using		 namespace std;

#define		N_SD	4
#define		DATA_FILE_IN	"UCNASimOutput.txt"

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

int main()
{
  Int_t eventID;
  Int_t ptclSpecies;

  Double_t primKE;
  Double_t primPos[4];	// fourth element is radius from the z axis
  Double_t primMo[3];
  Double_t primTime;
  Double_t primWeight;

  Int_t trapped;
  Double_t compTime;

  // variables related to scintillator crystal and mwpc active region
  Double_t edepScint[2];
  Double_t edepQScint[2];
  Double_t mwpcPos[2][3];	// 2 is for EAST/WEST, 3 is for 3 entries of G4ThreeVector
  Double_t mwpcPos2[2][3];
  Double_t scintPos[2][3];
  Double_t scintPos2[2][3];
  Double_t edepMWPC[2];
  Double_t scintHitTime[2];

  // variables related to all SD's. ATM it's just the above 4. But eventually 24.
  Double_t edepSD[N_SD];
  Double_t thetaInSD[N_SD];
  Double_t thetaOutSD[N_SD];
  Double_t keInSD[N_SD];
  Double_t keOutSD[N_SD];
//  Double_t edepAll;	// can't currently read all on event-by-event basis.
//  Int_t hitCountSD[N_SD];	// counts number of hits per SD. No ability to do this yet.
//  Double_t hitTimeSD[N_SD];	// earliest hit time in each SD

  TFile file("xuan_analyzed.root", "RECREATE");

  // note: for arrays you don't need to use & deallocator.
  // Probably cause the array variable is already a pointer.
  // For primitive types (including wrapper classes _t), you need &.
  TTree* anaTree = new TTree("anaTree", "tree for analysis");

/*// Note these are the variable names I would like to use.
  // however , I cannot since they need to match Michael's print out

  // initial variables read out by PrimaryGeneratorAction
  anaTree -> Branch("primaryEventID", &eventID, "eventID/I");
  anaTree -> Branch("primaryParticleSpecies", &ptclSpecies, "ptclSpecies/I");
  anaTree -> Branch("primaryKE", &primKE, "primKE/D");
  anaTree -> Branch("primaryPosition", primPos, "primPos[3]/D");
  anaTree -> Branch("primaryMomentum", primMo, "primMo[3]/D");
  anaTree -> Branch("primaryTime", &primTime, "primTime/D");
  anaTree -> Branch("primaryWeight", &primWeight, "primWeight/D");

  // global event print out from EventAction
  anaTree -> Branch("trapped", &trapped, "trapped/I");
  anaTree -> Branch("myCompTime", &compTime, "compTime/D");

  // variables stored by TrackerHit that are printed out by EventAction
  anaTree -> Branch("scintillatorEdep", edepScint, "EdepE/D:EdepW/D");
  anaTree -> Branch("scintillatorEdepQuenched", edepQScint, "EdepQE/D:EdepQW/D");
//  anaTree -> Branch("all SD edep", &edepAll, "edepAll/D");	// can't currently read
  anaTree -> Branch("mwpcEnergy", edepMWPC, "EdepMWPC_E/D:EdepMWPC_W/D");
  anaTree -> Branch("scintTimeToHit", scintHitTime, "scintTimeE/D:scintTimeW/D");

  char tmp[1024];
  sprintf(tmp, "EdepSD[%i]/D", N_SD);
  anaTree -> Branch("EdepSD", edepSD, tmp);
  sprintf(tmp, "thetaInSD[%i]/D", N_SD);
  anaTree -> Branch("thetaInSD", thetaInSD, tmp);
  sprintf(tmp, "thetaOutSD[%i]/D", N_SD);
  anaTree -> Branch("thetaOutSD", thetaOutSD, tmp);
  sprintf(tmp, "keInSD[%i]/D", N_SD);
  anaTree -> Branch("keInSD", keInSD, tmp);
  sprintf(tmp, "keOutSD[%i]/D", N_SD);
  anaTree -> Branch("keOutSD", keOutSD, tmp);

  anaTree -> Branch("MWPCPos", mwpcPos, "mwpcPosE[3]/D:mwpcPosW[3]/D");
  anaTree -> Branch("ScintPos", scintPos, "scintPosE[3]/D:scintPosW[3]/D");
  anaTree -> Branch("MWPCPosSigma", mwpcPos2, "mwpcPos2E[3]/D:mwpcPos2W[3]/D");
  anaTree -> Branch("ScintPosSigma", scintPos2, "scintPos2E[3]/D:scintPos2W[3]/D");
*/

  anaTree -> Branch("primKE", &primKE, "primKE/D");
  anaTree -> Branch("primTheta", primMo, "primMo[3]/D");
  anaTree -> Branch("primPos", primPos, "primPos[4]/D");
  anaTree -> Branch("primWeight", &primWeight, "primWeight/D");
  anaTree -> Branch("trapped", &trapped, "trapped/I");
  anaTree -> Branch("compTime", &compTime, "compTime/D");

  Int_t seed = 1;	// implementiong fake seed variable
  anaTree -> Branch("seed", &seed, "seed/I");

  anaTree -> Branch("Edep", edepScint, "EdepE/D:EdepW/D");
  anaTree -> Branch("EdepQ", edepQScint, "EdepQE/D:EdepQW/D");
  anaTree -> Branch("MWPCEnergy", edepMWPC, "MWPCEnergyE/D:MWPCEnergyW/D");
  anaTree -> Branch("time", scintHitTime, "timeE/D:timeW/D");

  Double_t trapMonTime[] = {0, 0};	// setting fake time from trap monitors
  Int_t subEvt = 0;			// setting subEvt flag as false for all events
  anaTree -> Branch("tmTime", trapMonTime, "tmTimeE/D:tmTimeW/D");
  anaTree -> Branch("subEvt", &subEvt, "subEvt/I");

  char tmp[1024];
  sprintf(tmp, "EdepSD[%i]/D", N_SD);
  anaTree -> Branch("EdepSD", edepSD, tmp);
  sprintf(tmp, "thetaInSD[%i]/D", N_SD);
  anaTree -> Branch("thetaInSD", thetaInSD, tmp);
  sprintf(tmp, "thetaOutSD[%i]/D", N_SD);
  anaTree -> Branch("thetaOutSD", thetaOutSD, tmp);
  sprintf(tmp, "keInSD[%i]/D", N_SD);
  anaTree -> Branch("keInSD", keInSD, tmp);
  sprintf(tmp, "keOutSD[%i]/D", N_SD);
  anaTree -> Branch("keOutSD", keOutSD, tmp);

  anaTree -> Branch("MWPCPos", mwpcPos, "MWPCPosE[3]/D:MWPCPosW[3]/D");
  anaTree -> Branch("ScintPos", scintPos, "ScintPosE[3]/D:ScintPosW[3]/D");
  anaTree -> Branch("MWPCPosSigma", mwpcPos2, "MWPCPosSigmaE[3]/D:MWPCPosSigmaW[3]/D");
  anaTree -> Branch("ScintPosSigma", scintPos2, "ScintPosSigmaE[3]/D:ScintPosSigmaW[3]/D");


  Int_t eid, sid;
  Double_t ke, x, y, z, px, py, pz, t, w;	// units keV, m, m, m, none, none, none, ns, none

  Int_t tFlag;
  Double_t compT;	// units s

  string se_name;
  Double_t se_e, se_eq, se_hitT;	// units keV, keV, s
  Int_t se_nTracks;
  Double_t se_ex, se_ey, se_ez, se_e2x, se_e2y, se_e2z, se_thetaIn, se_thetaOut, se_keIn, se_keOut;
	// units m*keV, m*keV, m*keV, m*m*keV... deg, deg, keV, keV

  string we_name;
  Double_t we_e, we_eq, we_hitT;
  Int_t we_nTracks;
  Double_t we_ex, we_ey, we_ez, we_e2x, we_e2y, we_e2z, we_thetaIn, we_thetaOut, we_keIn, we_keOut;

  string sw_name;
  Double_t sw_e, sw_eq, sw_hitT;        // units keV, keV, s
  Int_t sw_nTracks;
  Double_t sw_ex, sw_ey, sw_ez, sw_e2x, sw_e2y, sw_e2z, sw_thetaIn, sw_thetaOut, sw_keIn, sw_keOut;
        // units m*keV, m*keV, m*keV, m*m*keV... deg, deg, keV, keV

  string ww_name;
  Double_t ww_e, ww_eq, ww_hitT;
  Int_t ww_nTracks;
  Double_t ww_ex, ww_ey, ww_ez, ww_e2x, ww_e2y, ww_e2z, ww_thetaIn, ww_thetaOut, ww_keIn, ww_keOut;

  char fileName[] = DATA_FILE_IN;
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  int counter = 0;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);

    if(!bufstream.eof())
    {

      bufstream >> eid >> sid >> ke >> x >> y >> z >> px >> py >> pz >> t >> w
		>> tFlag >> compT
		>> se_name >> se_e >> se_eq >> se_hitT >> se_nTracks
		>> se_ex >> se_ey >> se_ez >> se_e2x >> se_e2y >> se_e2z >> se_thetaIn >> se_thetaOut >> se_keIn >> se_keOut
                >> we_name >> we_e >> we_eq >> we_hitT >> we_nTracks
                >> we_ex >> we_ey >> we_ez >> we_e2x >> we_e2y >> we_e2z >> we_thetaIn >> we_thetaOut >> we_keIn >> we_keOut
                >> sw_name >> sw_e >> sw_eq >> sw_hitT >> sw_nTracks
                >> sw_ex >> sw_ey >> sw_ez >> sw_e2x >> sw_e2y >> sw_e2z >> sw_thetaIn >> sw_thetaOut >> sw_keIn >> sw_keOut
                >> ww_name >> ww_e >> ww_eq >> ww_hitT >> ww_nTracks
                >> ww_ex >> ww_ey >> ww_ez >> ww_e2x >> ww_e2y >> ww_e2z >> ww_thetaIn >> ww_thetaOut >> ww_keIn >> ww_keOut;


    }

    // related to primaries and event info
    eventID = eid;
    ptclSpecies = sid;
    primKE = ke;
    primPos[0] = x;
    primPos[1] = y;
    primPos[2] = z;
    primPos[3] = sqrt(primPos[0]*primPos[0] + primPos[1]*primPos[1]);
    primMo[0] = px;
    primMo[1] = py;
    primMo[2] = pz;
    primTime = t;
    primWeight = w;
    trapped = tFlag;
    compTime = compT;

    // now we begin recording SD info
    edepScint[0] = se_e;
    edepScint[1] = sw_e;
    edepQScint[0] = se_eq;
    edepQScint[1] = sw_eq;
    scintHitTime[0] = se_hitT;
    scintHitTime[1] = sw_hitT;
    scintPos[0][0] = se_ex; scintPos[0][1] = se_ey; scintPos[0][2] = se_ez;
    scintPos[1][0] = sw_ex; scintPos[1][1] = sw_ey; scintPos[1][2] = sw_ez;
    scintPos2[0][0] = se_e2x; scintPos2[0][1] = se_e2y; scintPos2[0][2] = se_e2z;
    scintPos2[1][0] = sw_e2x; scintPos2[1][1] = sw_e2y; scintPos2[1][2] = sw_e2z;
    edepMWPC[0] = we_e;
    edepMWPC[1] = ww_e;
    mwpcPos[0][0] = we_ex; mwpcPos[0][1] = we_ey; mwpcPos[0][2] = we_ez;
    mwpcPos[1][0] = ww_ex; mwpcPos[1][1] = ww_ey; mwpcPos[1][2] = ww_ez;
    mwpcPos2[0][0] = we_e2x; mwpcPos2[0][1] = we_e2y; mwpcPos2[0][2] = we_e2z;
    mwpcPos2[1][0] = ww_e2x; mwpcPos2[1][1] = ww_e2y; mwpcPos2[1][2] = ww_e2z;


    anaTree -> Fill();
  }





  anaTree -> Write();








  cout << "-------------- End of Program ---------------" << endl;

  return 0;
}
