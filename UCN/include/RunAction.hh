#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class G4LogicalVolume;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void   EndOfRunAction(const G4Run*);

    int GetKillCount() { return fKillCount; };
    void IncrementKillCount() { fKillCount++; };

  private:
    int fKillCount;	// number of "trapped" events that get killed

};

#endif

