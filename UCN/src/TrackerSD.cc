#include <cmath>
#include <cassert>

#include "TrackerSD.hh"

#include <G4SystemOfUnits.hh>
#include <G4HCofThisEvent.hh>
#include <G4Step.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4ios.hh>
#include <G4VProcess.hh>
#include <G4LossTableManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4Gamma.hh>

using namespace std;

TrackerSD::TrackerSD(G4String sdname, G4String hcname):
 G4VSensitiveDetector(sdname), kb(0.01907*cm/MeV), rho(1.032*g/cm3)
{
  SetName(sdname);

  new TrackerSDMessenger(this);
  // G4VSensitiveDetector class object maintains a "collectionName" vector
  // which is the name of the hits collection defined in teh sensitive detector object.
  // In the constructor, the name of the hits collection must be defined.
  collectionName.insert(hcname);
}

// quenching calculation... see Junhua's thesis
double TrackerSD::QuenchFactor(double E) const
{
        const G4double a = 116.7*MeV*cm*cm/g;           // dEdx fit parameter a*e^(b*E)
        const G4double b = -0.7287;                                     // dEdx fit parameter a*e^(b*E)
        const G4double dEdx = a*rho*pow(E/keV,b);       // estimated dE/dx
        return 1.0/(1+kb*dEdx);
}

// Initialize method is invoked at the beginning of each event. Here you must instantiate a hits collection object
// and set it to the G4HCofThisEvent object
void TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // make a new hits collection and register it for this event
  fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName,collectionName[0]);
  G4int hcID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection);
  hce -> AddHitsCollection(hcID, fHitsCollection);

  // make a single hit for this event. It will be remade at the conclusion of event
  // Previous statement not necessarily true. May need to reset explicitly.
  // So far, this is a debugging check to make sure the total energy summed is correct.
  // I've tested this against SteppingAction accumulation. Using it to verify Michael's tracking recording.
//  fHitsCollection->insert(new TrackerHit());


  // this is an explicit resetting of the tracker hit collection we're keeping
  // This is Michaels tracker variables.
  fTrackerHitList.clear();
  fTrackOriginEnergy.clear();
}

// Make a new TrackerHit object for each track. Store them in HitsCollection.
// Store that instances track information in the TrackerHit.
// In EventAction, loop over and sum all values.
G4bool TrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
/*
  G4double edepStep = aStep -> GetTotalEnergyDeposit();

  if(edepStep == 0.) return false;	// get out if in our SD, no energy is deposited

  // access first entry of fHitsCollection since that is the track-by-track tracker
  TrackerHit* firstHit = (*fHitsCollection)[0];

  // Accumulate values in TrackerHit objects now.
  firstHit -> Add(edepStep);
*/
  // BELOW IS THE ACTUAL CODE FOR USING TRACKER SD AND HITS FOR ACCUMULATION.
  G4Track* aTrack = aStep -> GetTrack();
  G4String creatorProcessName = "";
  const G4VProcess* creatorProcess = aTrack -> GetCreatorProcess();
  if(creatorProcess == NULL)
    creatorProcessName = "original";
  else
    creatorProcessName = creatorProcess -> GetProcessName();

  G4StepPoint* preStep = aStep -> GetPreStepPoint();
  G4StepPoint* postStep = aStep -> GetPostStepPoint();
  G4ThreeVector prePos = preStep -> GetPosition();
  G4ThreeVector postPos = postStep -> GetPosition();
  G4double Epre = preStep -> GetKineticEnergy();
  G4double Epost = postStep -> GetKineticEnergy();
  G4double avgE = 0.5*(Epre + Epost);

  // Get prior track, or initialize a new one
  G4int currentTrackID = aTrack -> GetTrackID();
  map<G4int, TrackerHit*>::iterator myTrack = fTrackerHitList.find(currentTrackID);

  if(myTrack == fTrackerHitList.end())
  {
    TrackerHit* newHit = new TrackerHit();

    newHit -> SetTrackID(currentTrackID);
    newHit -> SetPtclSpeciesID(aTrack -> GetDefinition() -> GetPDGEncoding());
    newHit -> SetProcessName(creatorProcessName);
    newHit -> SetIncidentEnergy(preStep -> GetKineticEnergy());
    newHit -> SetHitPos(postPos);
    newHit -> SetHitTime(preStep -> GetGlobalTime());
    newHit -> SetIncidentMomentum(preStep -> GetMomentum());
    // is pre-step physical volume defined? If not, set the name as "Unknown"
    G4VPhysicalVolume* preVolume = preStep -> GetPhysicalVolume();
    newHit -> SetVolumeName(preVolume? preVolume -> GetName(): "Unknown");
    newHit -> SetTrackVertex(aTrack -> GetVertexPosition());
    newHit -> SetCreatorVolumeName(aTrack -> GetLogicalVolumeAtVertex() -> GetName());
    (*newHit).fNbSecondaries = 0;

    map<const G4Track*, double>::iterator itorig = fTrackOriginEnergy.find(aTrack);
    if(itorig == fTrackOriginEnergy.end())
    {
      // fOriginEnergy = 0 for primary tracks (not a sub-track of another track in this volume)
      (*newHit).fOriginEnergy = 0;
    }
    else
    {
      // Get previously stored origin energy for secondary tracks. Remove listing entry
      (*newHit).fOriginEnergy = itorig -> second;
      fTrackOriginEnergy.erase(itorig);
    }

    int hitNb = fHitsCollection -> insert(newHit);
    fTrackerHitList.insert(pair<G4int, TrackerHit*>(currentTrackID, (TrackerHit*)fHitsCollection->GetHit(hitNb - 1)));
    myTrack = fTrackerHitList.find(currentTrackID);
  }

  // accumulate edep, edep quenched, and local position for this step
  G4double edep = aStep -> GetTotalEnergyDeposit();
  // fetch the TrackerHit fOriginEnergy. If 0, plug in avgE into call to QuenchFactor. Else, use fOriginEnergy.
  G4double edepQuenched;
  if((myTrack -> second -> fOriginEnergy) == 0)
    edepQuenched = edep*QuenchFactor(avgE);
  else
    edepQuenched = edep*QuenchFactor(myTrack -> second -> fOriginEnergy);

  G4ThreeVector localPos = preStep -> GetTouchableHandle() -> GetHistory() -> GetTopTransform().TransformPoint(prePos);
  myTrack -> second -> AddEdep(edep, localPos);
  myTrack -> second -> AddEdepQuenched(edepQuenched);
  myTrack -> second -> SetExitMomentum(postStep -> GetMomentum());

  // record origin energy for secondaries in same volume
  const G4TrackVector* secondaries = aStep -> GetSecondary();
  while(myTrack -> second -> fNbSecondaries < secondaries -> size())
  {
    // don't really understand the argument of this line.
    const G4Track* sTrack = (*secondaries)[myTrack -> second -> fNbSecondaries++];

    if(sTrack -> GetVolume() != sTrack -> GetVolume())
    {  continue; }

    const G4double eOrig = ((myTrack->second->fOriginEnergy)>0)?(myTrack->second->fOriginEnergy):(avgE);

    if(fTrackOriginEnergy.find(sTrack) != fTrackOriginEnergy.end())
    {
      // Add code that increments a kill counter or kills the entire event here.
      G4cout << "Apparently we have a duplicate secondary. \n"
             << "This event should be invalid. Make of a note of it." << G4endl;
    }
    fTrackOriginEnergy.insert(pair<const G4Track*, double>(sTrack, eOrig));
  }

  return true;

}

void TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  // This can be used for some cool stuff but currently I have no use for it.
  // Some notes from the MIT GEANT4 tutorial:
  // This method is invoked at the end of processing an event (obvi)
  // It is invoked even if the event is aborted. Could be useful.
  // It is invoked before UserEndOfEventAction.

}



// Here is some simplified code I wrote from M.Mendenhall's stuff in an attempt to debug.
// Saving it here. It's basically using only simple objects (TrackerHit, fHitsCollection)
// to record all info on primaries, with no sensible accumulation/recording of secondaries.
// Could be useful for debugging in the future but currently don't need it.
/*  G4Track* aTrack = aStep -> GetTrack();
  G4String creatorProcessName = "";
  const G4VProcess* creatorProcess = aTrack -> GetCreatorProcess();
  if(creatorProcess == NULL)
    creatorProcessName = "original";
  else
    creatorProcessName = creatorProcess -> GetProcessName();

  G4StepPoint* preStep = aStep -> GetPreStepPoint();
  G4StepPoint* postStep = aStep -> GetPostStepPoint();
  G4ThreeVector prePos = preStep -> GetPosition();
  G4ThreeVector postPos = postStep -> GetPosition();
  G4double Epre = preStep -> GetKineticEnergy();
  G4double Epost = postStep -> GetKineticEnergy();
  G4double avgE = 0.5*(Epre + Epost);

  TrackerHit* hit = new TrackerHit();

  hit -> SetTrackID(aTrack -> GetTrackID());
  hit -> SetPtclSpeciesID(aTrack -> GetDefinition() -> GetPDGEncoding());
  hit -> SetProcessName(creatorProcessName);
  hit -> SetIncidentEnergy(preStep -> GetKineticEnergy());
  hit -> SetHitPos(postPos);
  hit -> SetHitTime(preStep -> GetGlobalTime());        // there's a better method in example B5 on this
  hit -> SetIncidentMomentum(preStep -> GetMomentum());
  // is pre-step physical volume defined? If not, set the name as "Unknown"
  G4VPhysicalVolume* preVolume = preStep -> GetPhysicalVolume();
  hit -> SetVolumeName(preVolume? preVolume -> GetName(): "Unknown");
  hit -> SetTrackVertex(aTrack -> GetVertexPosition());
  hit -> SetCreatorVolumeName(aTrack -> GetLogicalVolumeAtVertex() -> GetName());
  (*hit).fNbSecondaries = 0;

  (*hit).fOriginEnergy = 0;     // explicitly set for now. Since ignoring secondaries atm

  // Copied over from Michael's code in order to calculate e quenched.
  map<const G4Track*, double>::iterator itorig = fTrackOriginEnergy.find(aTrack);
  if(itorig == fTrackOriginEnergy.end())
  {
    // fOriginEnergy = 0 for primary tracks (not a sub-track of another track in this volume)
    (*hit).fOriginEnergy = 0;
  }
  else
  {
    // Get previously stored origin energy for secondary tracks. Remove listing entry
    (*hit).fOriginEnergy = itorig -> second;
    fTrackOriginEnergy.erase(itorig);
  }

  G4double edep = aStep -> GetTotalEnergyDeposit();
  // This is always going to trigger the if part since now it goes sequentially.
  // Will need to go back and think of how to deal with this.
  G4double edepQuenched;
  if((hit -> fOriginEnergy) == 0)
    edepQuenched = edep*QuenchFactor(avgE);
  else
    edepQuenched = edep*QuenchFactor(hit -> fOriginEnergy);

  G4ThreeVector localPos = preStep -> GetTouchableHandle() -> GetHistory() -> GetTopTransform().TransformPoint(prePos);
  hit -> AddEdep(edep, localPos);
  hit -> AddEdepQuenched(edepQuenched);
  hit -> SetExitMomentum(postStep -> GetMomentum());

  // And then stuff pertaining to secondaries. Ignore for now. See if this other stuff works.

  fHitsCollection -> insert(hit);       // this gonna go at the end. Won't use complicated C++ objects.

  return true;
*/






// ----- Messenger classes. Won't use. ----- //

TrackerSDMessenger::TrackerSDMessenger(TrackerSD* T): mySD(T) {
        sdDir = new G4UIdirectory(("/SD/"+mySD->GetName()+"/").c_str());
        sdDir->SetGuidance("Sensitive detector response settings");

        kbCmd = new G4UIcmdWithADouble((sdDir->GetCommandPath()+"kb").c_str(), this);
        kbCmd->SetGuidance("Birk's Law quenching constant in cm/MeV");
        kbCmd->SetDefaultValue(0.01907);
        kbCmd->AvailableForStates(G4State_Idle);
}

TrackerSDMessenger::~TrackerSDMessenger() {
        delete kbCmd;
        delete sdDir;
}

void TrackerSDMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
        if( command == kbCmd ) {
                G4double k = kbCmd->GetNewDoubleValue(newValue);
                G4cout << "Setting Birk's Law kb = " << k << " cm/MeV for " << mySD->GetName() << G4endl;
                mySD->SetKb(k * cm/MeV);
        }
}

//----------------------------------------------------------------

