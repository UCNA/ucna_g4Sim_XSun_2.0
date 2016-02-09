#ifndef EventAction_h
#define EventAction_h 1

#include "DetectorConstruction.hh"
#include "TrackerHit.hh"

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <G4Event.hh>

#include <time.h>

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* evt);
    virtual void EndOfEventAction(const G4Event* evt);

    void SetTrappedTrue() { fTrapped = true; };

    clock_t GetStartTime() { return fStartTime; };

  private:
    const DetectorConstruction* fMyDetectorConstruction;

    bool  fTrapped;             // check if event was killed due to being trapped

    clock_t fStartTime;		// time.h uses to define 'trapped' ptcl's & kill them

    G4String fSDNames[fNbSDs];
    G4int fHitsCollectionIDs[fNbSDs]; 	// note: fNbSDs is defined in DetectorConstruction.hh
					// it's a global const variable.
					// Since #include "..." copies the code, this should (and is) fine.
    int fScintEast_index;
    int fScintWest_index;
    int fActiveWireVolEast_index;
    int fActiveWireVolWest_index;
};

#endif


