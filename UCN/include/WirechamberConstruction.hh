#ifndef WirechamberConstruction_HH
#define WirechamberConstruction_HH

#include "DetectorConstructionUtils.hh"
#include "WirechamberActiveRegionConstruction.hh"

class WirechamberConstruction: public MaterialUser
{
public:
  WirechamberConstruction();	// constructor


  G4double GetWidth() { return 2*dMWPCContainerHalf_Z; };

  G4double dWindowThick;	// Mylar window thickness
  G4double dMWPCEntranceR;	// entrance window radius
  G4double dMWPCExitR;		// exit window radius
  G4Material* mMWPCActiveRegionGas;		// MWPC fill gas
  G4double dEntranceToCathodes;	// entrance-window-to-cathode distance
  G4double dExitToCathodes;	// exit window to cathode distance

  WirechamberActiveRegionConstruction ActiveRegion;	// active gas region with wireplanes

  G4LogicalVolume* container_log;	// overall gas box
  G4LogicalVolume* kevContainer_log;	// container volume for kevlar strip array
  G4LogicalVolume* kevSeg_log;		// one segment of kevlay strip array
  G4LogicalVolume* kevStrip_log;	// kevlar strip in one segment
  G4LogicalVolume* winIn_log;		// inner window
  G4LogicalVolume* winOut_log;		// outer window

  void Build(int side);

  G4RotationMatrix* fMyRotation;	// rotation from global frame to local coordinates
  G4ThreeVector vMyTranslation;		// translation from global coordinates to center of anode plane
  G4double dE0;				// field scaling constant. Needed here to get passed to MWPCField later

protected:
  G4double dMWPCContainerHalf_Z;	// half width of wirechamber



};

#endif
