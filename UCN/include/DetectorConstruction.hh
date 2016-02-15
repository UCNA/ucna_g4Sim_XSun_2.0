#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "TrackerSD.hh"
#include "DetectorConstructionUtils.hh"
#include "SourceHolderConstruction.hh"
#include "DecayTrapConstruction.hh"
#include "WirechamberConstruction.hh"
#include "ScintillatorConstruction.hh"
#include "PackageDetConstruction.hh"

#include "G4VUserDetectorConstruction.hh"

#include <G4ElectroMagneticField.hh>	// Taken from WirechamberConstruction.
#include <G4MagneticField.hh>
#include <G4RotationMatrix.hh>

#include <string>

const int fNbSDs = 4;

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.
class DetectorConstruction : public G4VUserDetectorConstruction, MaterialUser
{
  public:
    DetectorConstruction();		// Constructor/destructors
    virtual ~DetectorConstruction();
    virtual G4VPhysicalVolume* Construct();

    G4LogicalVolume* experimentalHall_log;	// world volume
    G4VPhysicalVolume* experimentalHall_phys;

    SourceHolderConstruction Source;		// individual components.
    G4VPhysicalVolume* source_phys;

    DecayTrapConstruction Trap;

    PackageDetConstruction DetPackage[2];
    G4VPhysicalVolume* detPackage_phys[2];



    G4String fSDNamesArray[fNbSDs];	// needs to be public since EventAction will access all elements
    G4String fHCNamesArray[fNbSDs];

  private:
    // field constructors
    void ConstructGlobalField();
    void ConstructEastMWPCField(G4double a, G4double b, G4double c, G4double d,
				G4RotationMatrix* e, G4ThreeVector f);
    void ConstructWestMWPCField(G4double a, G4double b, G4double c, G4double d,
				G4RotationMatrix* e, G4ThreeVector f);
				// a = active region wire spacing
				// b = active region plane spacing
				// c = active region anode radius
				// d = mwpc electric potential
				// e = rotation matrix of our coordinate system
				// f = translation vector of our coordinate system

    // Register and store each SD
    TrackerSD* RegisterSD(G4String sdName, G4String hcName);

    TrackerSD* SD_scint_scintillator[2];	// all the SD objects that will be used
    TrackerSD* SD_scint_deadScint[2];
    TrackerSD* SD_scint_backing[2];
    TrackerSD* SD_mwpc_winIn[2];
    TrackerSD* SD_mwpc_winOut[2];
    TrackerSD* SD_decayTrap_windows[2];
    TrackerSD* SD_mwpc_kevStrip[2];
    TrackerSD* SD_wireVol[2];
    TrackerSD* SD_wireVol_planes[2];
    TrackerSD* SD_mwpc_container[2];
    TrackerSD* SD_source;
    TrackerSD* SD_decayTrap_innerMonitors[2];
    TrackerSD* SD_world;

    // User Interface commands from .mac files
    G4float fSourceFoilThick;		// source foil full thickness
    G4ThreeVector vSourceHolderPos;	// source holder position

    G4float fCrinkleAngle;		// crinkle angle of wiggle foils to be implemented later

    // some of my own tools to help with DetectorConstruction
    int fStorageIndex;
    bool bUseSourceHolder;
    G4double fScintStepLimit;
};

#endif

