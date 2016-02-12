#ifndef ScintillatorConstruction_HH
#define ScintillatorConstruction_HH

#include "DetectorConstructionUtils.hh"

class ScintillatorConstruction: public MaterialUser
{
public:
  ScintillatorConstruction();	// constructor

  G4double GetScintFacePos() { return dScintFace_PosZ; };
  G4double GetWidth() { return dN2Volume_Z; };

  G4double dScintRadius;	// scintillator disc radius
  G4double dBackingRadius;	// backing veto (and overall volume) radius
  G4double dScintThick;		// scintillator disc thickness
  G4double dDeadThick;		// dead scintillator thickness (3um according to Junhua's thesis)
  G4double dBackingThick;	// backing veto thickness (M.M.'s guess)
  G4double dLightGuideThick;	// light guide thickness at scintillator edge (M.M's guess)
				// ^^ sets scintillator to backing distance

  G4LogicalVolume* container_log;	// overall container volume (filled with nitrogen)
  G4LogicalVolume* deadLayer_log;	// scintillator dead layer logical volume
  G4LogicalVolume* scintillator_log;	// actual scintillator (active region) volume
  G4LogicalVolume* lightGuide_log;	// light guide material volume
  G4LogicalVolume* backing_log;		// backing veto logical volume

  void Build(int side);		// construct only the logical container volume

protected:
  G4double dScintFace_PosZ;
  G4double dN2Volume_Z;
  G4VPhysicalVolume* deadLayer_phys;
  G4VPhysicalVolume* scintillator_phys;
  G4VPhysicalVolume* lightGuide_phys;
  G4VPhysicalVolume* backing_phys;

};

#endif
