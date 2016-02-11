#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "DetectorConstruction.hh"

#include "G4VUserPrimaryGeneratorAction.hh"	// original example used these 3
#include "G4ParticleGun.hh"
#include "globals.hh"

#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4ParticleGun.hh>
#include <G4Event.hh>
#include <G4VUserEventInformation.hh>


class G4ParticleGun;
class G4Event;

// user event information for recording primary event weighting
class PrimEvtWeighting : public G4VUserEventInformation
{
  public:
    PrimEvtWeighting(double W): w(W) {}	// constructor
    void Print() const { G4cout << "Primary weighting: " << w << G4endl; }

    double w;				// event primary weight
};

using namespace std;

struct event
{
  G4int event_gen_id;
  G4double event_energy;	// turns into keV
  G4int event_speciesFlag;      // 11 means electron, 22 is gamma
  G4double event_xMo;		// momentum is unitless vector
  G4double event_yMo;
  G4double event_zMo;
  G4double event_xPos;		// turns into m
  G4double event_yPos;
  G4double event_zPos;
  G4double event_time;		// turns into s but always 0. I think supposed to be ns
  G4double event_weight;	// unitless
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction(DetectorConstruction*);
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    void GeneratePrimaries(G4Event*);

    // method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  private:
    G4ParticleGun*  fParticleGun; 	// pointer a to G4 gun class
    DetectorConstruction* fMyDetector;	// pointer to the detector geometry class

    G4double fSourceRadius;		// spread radius for source droplets
    G4ThreeVector fPosOffset;		// base positioning offset
    event fEvtsArray[10000];		// size has to be number of lines in input file

    void LoadFile(G4String fileName);
    void DiskRandom(G4double radius, G4double& x, G4double& y);
    void DisplayGunStatus();
    void SavePrimPtclInfo(int index);
    void Set_113SnSource();

};

#endif


