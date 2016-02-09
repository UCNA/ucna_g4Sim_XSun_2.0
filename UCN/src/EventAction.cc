#include "EventAction.hh"
#include "RunAction.hh"
#include "TrackerHit.hh"
#include "TrackerSD.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
using   namespace       std;

#define	OUTPUT_FILE	"UCNASimOutput.txt"

struct SingleSDInfo
{
  // these initalized variable values are needed or else you'll run into same errors as in TrackerHit
  SingleSDInfo() : name("no name"), energy(0), energyQuenched(0), hitTime(0),
		   NbOfTracks(0), edepWeightedPos(), edepWeightedPos2(), thetaIn(0), thetaOut(0) {}

  G4String name;	// name of SD we're storing info from.
  G4double energy;	// energy deposited in SD
  G4double energyQuenched;	// energy quenched deposited in SD

  G4double hitTime;	// time it takes to first hit inside the SD
  int NbOfTracks;	// number of tracks produced (primary only I believe)
//  int hitCount;	// number of hits in a particular SD (haven't been able to implement yet)

  G4ThreeVector edepWeightedPos;	// energy weighted by position vector
  G4ThreeVector edepWeightedPos2;	// energy weighted by position vector elements squared

  G4double thetaIn;	// angle (degrees) of electron from normal when entering each SD
  G4double thetaOut;	// angle (degrees) of electron from normal when exiting
  G4double keIn;
  G4double keOut;

};

struct EventInfo
{
  bool trapped;         // whether the particle was deemed trapped and killed
  G4double compTime;    // time it took from beginning of event generation to completion of event propagation

  SingleSDInfo Scint_crystals[2];
  SingleSDInfo MWPC_activeWireVol[2];
};

EventAction::EventAction()
: G4UserEventAction(), fStartTime(0)
{
  fMyDetectorConstruction = static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  // initialize all values to -1 since that doesn't correspond to any actual SD

  for(int i = 0; i < fNbSDs; i++)	// see comment in EventAction.hh on fNbSDs
  {
    fHitsCollectionIDs[i] = -1;
    fSDNames[i] = "initialized_in_EventAction_to_nothing";
  }

  fScintEast_index = -1;
  fScintWest_index = -1;
  fActiveWireVolEast_index = -1;
  fActiveWireVolWest_index = -1;
}


EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  fTrapped = false;
  fStartTime = 0;

  // Sets the start of the C++ 'clock' used for tracking trapped ptcl's
  fStartTime = clock();

  if((evt->GetEventID())%1000 == 0)
    G4cout << "\n -------------- Begin of event: " << evt->GetEventID() << "\n" << G4endl;

  // get hits collection IDs only once in an overall run though.
  if(fHitsCollectionIDs[0] == -1)
  {
    for(int i = 0; i < fNbSDs; i++)
    {
      fSDNames[i] = (*fMyDetectorConstruction).fSDNamesArray[i];
      fHitsCollectionIDs[i] = G4SDManager::GetSDMpointer()->GetCollectionID((*fMyDetectorConstruction).fHCNamesArray[i]);
      G4cout << "HC ID: " << fHitsCollectionIDs[i] << ", belongs to HC name: " << (*fMyDetectorConstruction).fHCNamesArray[i] << G4endl;
      // these are the special ones we want to record extra info from.
      if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_scint_EAST")
      {
        fScintEast_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_scint_WEST")
      {
        fScintWest_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_wireVol_EAST")
      {
        fActiveWireVolEast_index = i;
      }
      else if((*fMyDetectorConstruction).fHCNamesArray[i] == "HC_wireVol_WEST")
      {
        fActiveWireVolWest_index = i;
      }
    }
  }
}


void EventAction::EndOfEventAction(const G4Event* evt)
{
  const G4double kMe = 510.998910*keV;	// electron mass (keV obvi)

  // kill counter
  if(fTrapped == true)
  {
    G4cout << "Event " << evt->GetEventID() << " trapped. Incrementing kill count." << G4endl;
    // some code that M. Mendenhall came up with that allows me to access my own RunAction class
    ((RunAction*)G4RunManager::GetRunManager()->GetUserRunAction()) -> IncrementKillCount();
  }

  // manually (not GEANT4 clock) store how much time the event ran for
  clock_t timeOfEvent = clock() - fStartTime;
  double compTime = ((double)timeOfEvent)/CLOCKS_PER_SEC;

  // Get hit collections to extract info from each SD
  G4HCofThisEvent* hce = evt->GetHCofThisEvent();
  if(!hce)
  {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl;
    G4Exception("EventAction::EndOfEventAction()", "Error", JustWarning, msg);
    return;
  }

  // All tracking of energy and appropriate variables will be stored here.
  TrackerHitsCollection* SD_totalHC[fNbSDs];
  TrackerHit* SD_hits[fNbSDs];
  SingleSDInfo SD_info[fNbSDs];

  // loop over the SD's where everything is stored in arrays
  for(int i = 0; i < fNbSDs; i++)
  {
    SD_totalHC[i] = static_cast<TrackerHitsCollection*>(hce->GetHC(fHitsCollectionIDs[i]));

    SD_info[i].name = fSDNames[i];	// set name manually because I know how I stored them. Not optimal.
					// but guaranteed to get the name set.

    // deals with seg fault error when it tries to access an SD that had no track inside it
    if(SD_totalHC[i]->GetSize() == 0)
      continue;

    // loop over TrackerHits per hits collection assigned to each SD
    for(unsigned int t = 0; t < SD_totalHC[i] -> GetSize(); t++)
    {
      // assign array of TrackerHit objects to have info extracted from it
      SD_hits[i] = (*SD_totalHC[i])[t];

      if(t == 0)	// get info related to first track into SD i.e. first entry of hits collection
      {
        SD_info[i].hitTime = SD_hits[i] -> GetHitTime();
        SD_info[i].NbOfTracks = SD_totalHC[i] -> GetSize();

        G4double pIn_x = (SD_hits[i] -> GetIncidentMomentum()).x();
        G4double pIn_y = (SD_hits[i] -> GetIncidentMomentum()).y();
        G4double pIn_z = (SD_hits[i] -> GetIncidentMomentum()).z();
        if(pIn_x || pIn_y || pIn_z)
        {
          G4double magPIn2 = pIn_x*pIn_x + pIn_y*pIn_y + pIn_z*pIn_z;
	  // these values only useful/mean what they are supposed to for first track into SD
          SD_info[i].keIn = sqrt(magPIn2 + kMe*kMe) - kMe;
          SD_info[i].thetaIn = acos(fabs(pIn_z/sqrt(magPIn2)));
        }
      }

      // accumulation of variables across all tracks.
      SD_info[i].energy += SD_hits[i] -> GetEdep();
      SD_info[i].energyQuenched += SD_hits[i] -> GetEdepQuenched();

      SD_info[i].edepWeightedPos += SD_hits[i] -> GetEdepPos();
      SD_info[i].edepWeightedPos2 += SD_hits[i] -> GetEdepPos2();

      if(t == (SD_totalHC[i] -> GetSize() - 1))
      {
        G4double pOut_x = (SD_hits[i] -> GetExitMomentum()).x();
        G4double pOut_y = (SD_hits[i] -> GetExitMomentum()).y();
        G4double pOut_z = (SD_hits[i] -> GetExitMomentum()).z();
        if(pOut_x || pOut_y || pOut_z)
        {
          G4double magPOut2 = pOut_x*pOut_x + pOut_y*pOut_y + pOut_z*pOut_z;
          // Same as above, exit values only useful for last entry.
          SD_info[i].keOut = sqrt(magPOut2 + kMe*kMe) - kMe;
          SD_info[i].thetaOut = acos(fabs(pOut_z/sqrt(magPOut2)));
        }
      }
    }
  }

  EventInfo event;
  event.trapped = fTrapped;
  event.compTime = compTime*s;
  event.Scint_crystals[0] = SD_info[fScintEast_index];
  event.Scint_crystals[1] = SD_info[fScintWest_index];
  event.MWPC_activeWireVol[0] = SD_info[fActiveWireVolEast_index];
  event.MWPC_activeWireVol[1] = SD_info[fActiveWireVolWest_index];


  // print out of class member variables due to stepping action
  ofstream outfile;
  outfile.open(OUTPUT_FILE, ios::app);

  outfile << event.trapped << "\t"
	  << event.compTime/s << "\t";

  for(int j = 0; j < 2; j++)
  {
  outfile << event.Scint_crystals[j].name << "\t"
  	  << event.Scint_crystals[j].energy/keV << "\t"
  	  << event.Scint_crystals[j].energyQuenched/keV << "\t"
	  << event.Scint_crystals[j].hitTime/s << "\t"
	  << event.Scint_crystals[j].NbOfTracks << "\t"
	  << (event.Scint_crystals[j].edepWeightedPos).x()/(m*keV) << "\t"
          << (event.Scint_crystals[j].edepWeightedPos).y()/(m*keV) << "\t"
          << (event.Scint_crystals[j].edepWeightedPos).z()/(m*keV) << "\t"
	  << (event.Scint_crystals[j].edepWeightedPos2).x()/(m*m*keV) << "\t"
          << (event.Scint_crystals[j].edepWeightedPos2).y()/(m*m*keV) << "\t"
          << (event.Scint_crystals[j].edepWeightedPos2).z()/(m*m*keV) << "\t"
	  << event.Scint_crystals[j].thetaIn/deg << "\t"
	  << event.Scint_crystals[j].thetaOut/deg << "\t"
	  << event.Scint_crystals[j].keIn/keV << "\t"
	  << event.Scint_crystals[j].keOut/keV << "\t"

	  << event.MWPC_activeWireVol[j].name << "\t"
          << event.MWPC_activeWireVol[j].energy/keV << "\t"
          << event.MWPC_activeWireVol[j].energyQuenched/keV << "\t"
          << event.MWPC_activeWireVol[j].hitTime/s << "\t"
          << event.MWPC_activeWireVol[j].NbOfTracks << "\t"
          << (event.MWPC_activeWireVol[j].edepWeightedPos).x()/(m*keV) << "\t"
          << (event.MWPC_activeWireVol[j].edepWeightedPos).y()/(m*keV) << "\t"
          << (event.MWPC_activeWireVol[j].edepWeightedPos).z()/(m*keV) << "\t"
          << (event.MWPC_activeWireVol[j].edepWeightedPos2).x()/(m*m*keV) << "\t"
          << (event.MWPC_activeWireVol[j].edepWeightedPos2).y()/(m*m*keV) << "\t"
          << (event.MWPC_activeWireVol[j].edepWeightedPos2).z()/(m*m*keV) << "\t"
          << event.MWPC_activeWireVol[j].thetaIn/deg << "\t"
          << event.MWPC_activeWireVol[j].thetaOut/deg << "\t"
          << event.MWPC_activeWireVol[j].keIn/keV << "\t"
          << event.MWPC_activeWireVol[j].keOut/keV << "\t";

  }

  outfile << "\n";

  outfile.close();
}

// locFlag = 0 -> EAST
// 	     1 -> WEST
// typeFlag = 0 -> Scint
//	      1 -> MWPC

