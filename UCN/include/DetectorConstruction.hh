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
#include <G4MagneticField.hh>		// Needed since all EM fields constructed here
#include <G4RotationMatrix.hh>

#include <G4UImessenger.hh>		// Taken from DetectorConstruction.hh M.M's
#include <G4UIdirectory.hh>
#include <G4UIcommand.hh>
#include <G4UIcmdWithADoubleAndUnit.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWith3VectorAndUnit.hh>
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithAString.hh>

#include <string>

const int fNbSDs = 4;

class G4VPhysicalVolume;
class G4LogicalVolume;

/// Detector construction class to define materials and geometry.
class DetectorConstruction : public G4VUserDetectorConstruction, G4UImessenger, MaterialUser
{
  public:
    DetectorConstruction();		// Constructor/destructors
    virtual ~DetectorConstruction();
    virtual G4VPhysicalVolume* Construct();

    virtual void SetNewValue(G4UIcommand * command,G4String newValue);	// UI communicator

    G4LogicalVolume* experimentalHall_log;	// world volume
    G4VPhysicalVolume* experimentalHall_phys;

    SourceHolderConstruction Source;		// individual components.
    G4VPhysicalVolume* source_phys;

    DecayTrapConstruction Trap;

    PackageDetConstruction DetPackage[2];
    G4VPhysicalVolume* detPackage_phys[2];

    G4ThreeVector GetSourceHolderPos() { return vSourceHolderPos; };


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
    G4UIdirectory* uiDetectorDir;	// UI directory for detector-related commands

    G4UIcmdWithAString* uiDetectorGeometryCmd;	// which detector geometry to construct
    G4String sGeometry;

    G4UIcmdWith3VectorAndUnit* uiSourceHolderPosCmd;	// source holder position
    G4ThreeVector vSourceHolderPos;

    G4UIcmdWith3VectorAndUnit* uiDetOffsetCmd;	// symmetrical detector offset from center origin
    G4ThreeVector vDetOffset;

    G4UIcmdWithADouble* uiDetRotCmd;		// symmetrical detector rotation angle around Z axis (radians, hence no units)
    float fDetRot;

    G4UIcmdWithABool* uiUseFoilCmd;		// construction of Indium 10um Al source foil
    bool bUseFoil;

    G4UIcmdWithADoubleAndUnit* uiVacuumLevelCmd;	// SCS bore vacuum
    G4float fVacuumPressure;

    G4UIcmdWithADoubleAndUnit* uiScintStepLimitCmd;	// step size limiter in scintillator
    G4float fScintStepLimit;

    G4UIcmdWithADoubleAndUnit* uiSourceFoilThickCmd;	// source foil full thickness
    G4float fSourceFoilThick;

//    G4float fCrinkleAngle;		// crinkle angle of wiggle foils to be implemented later, maybe

    // some of my own tools to help with DetectorConstruction
    int fStorageIndex;
    G4float fCrinkleAngle;

};

#endif

