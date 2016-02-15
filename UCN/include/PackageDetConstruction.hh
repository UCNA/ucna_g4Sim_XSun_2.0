#ifndef PackageDetConstruction_HH
#define PackageDetConstruction_HH

#include "DetectorConstructionUtils.hh"
#include "WirechamberConstruction.hh"
#include "ScintillatorConstruction.hh"

class PackageDetConstruction: public MaterialUser
{
public:
  PackageDetConstruction();	// constructor. Always will put explicitly in .cc

  G4double GetScintFacePos() { return 0; };	// I think this is here for scintillator adjustments

  G4double dDetPackageRadius;
  G4double dMWPCEntranceThick;	// mwpc entrance tube wall thickness
  G4double dMWPCEntranceRad;	// mwpc entrance tube radius
  G4double dMWPCEntranceDepth;	// mwpc entrance tube depth
  G4double dFrontWinFrameThick;	// mwpc front window frame thickness
  G4double dBackWinFrameThick;	// mwpc exit window frame thickness

  ScintillatorConstruction Scint;	// scintillator assembly
  WirechamberConstruction Wirechamber;	// wirechamber assembly

  G4LogicalVolume* container_log;	// overall positioning container to hold the package
  G4LogicalVolume* mwpcEntrance_log;	// entrance port container
  G4LogicalVolume* entranceFront_log;	// entrance port front plate
  G4LogicalVolume* entranceMid_log;	// entrance port tube
  G4LogicalVolume* entranceBack_log;	// entrance port back plate (mwpc box cover)
  G4LogicalVolume* mwpcExit_log;	// aluminum exit window from wirechamber
  G4LogicalVolume* mwpcExitGasN2_log;	// N2 between exit window and scintillator
  G4LogicalVolume* backStuff_log;	// misc mass behind detectors

  G4double dEntranceFacePos_Z;	// entrance window port entrance relative to scint face
  G4double dEntranceWinPos_Z;	// mwpc entrance window position relative to scint face
  G4double dExitFramePos_Z;	// exit window frame position along z

  void Build(int side);		// construct only logical container volume

protected:
  G4VPhysicalVolume* scint_phys;
  G4VPhysicalVolume* mwpc_phys;

  G4VPhysicalVolume* mwpcEntrance_phys;
  G4VPhysicalVolume* entranceFront_phys;
  G4VPhysicalVolume* entranceMid_phys;
  G4VPhysicalVolume* entranceBack_phys;
  G4VPhysicalVolume* mwpcExit_phys;
  G4VPhysicalVolume* mwpcExitGasN2_phys;
  G4VPhysicalVolume* backStuff_phys;


};

#endif
