#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include <time.h>

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)
{ }


SteppingAction::~SteppingAction()
{ }


void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // kill switch in case simulation runs too long
  G4int stepNo = step -> GetTrack() -> GetCurrentStepNumber();
  clock_t timeSpentSoFar = clock() - ((EventAction*)G4EventManager::GetEventManager()->GetUserEventAction())->GetStartTime();

  double time = ((double)timeSpentSoFar)/CLOCKS_PER_SEC;

  if(stepNo >= 2000000 || time > 60)
  {
    G4cout << "----> Tracking killed by computation time limit." << G4endl;
    step -> GetTrack() -> SetTrackStatus(fStopAndKill);
    fEventAction -> SetTrappedTrue();	// sets the event action flag as true.
  }

}

