#ifndef WirechamberActiveRegionConstruction_HH
#define WirechamberActiveRegionConstruction_HH

#include "DetectorConstructionUtils.hh"

class WirechamberActiveRegionConstruction: public MaterialUser
{
public:
  WirechamberActiveRegionConstruction();	// constructor

  G4double GetWidth() { return 2*cm; };	// still can't remember why this is here

  G4Material* mMWPCGas;

  G4LogicalVolume* gas_log;		// container log volume of the mwpc active region gas volume
  G4LogicalVolume* cathSeg_log;
  G4LogicalVolume* anodeSeg_log;
  G4LogicalVolume* cathWire_log;
  G4LogicalVolume* cathPlate_log;
  G4LogicalVolume* anodeWire_log;

  // build logical volume container. Takes an int flag which designates side.
  void Build(int side);

  G4double dAnodeRadius;
  G4double dCathodeRadius;
  G4double dPlatingThick;
  G4double dWireSpacing;
  G4double iNbOfWires;
  G4double dPlaneSpacing;

};

#endif
