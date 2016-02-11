#ifndef DecayTrapConstruction_HH
#define DecayTrapConstruction_HH

#include "DetectorConstructionUtils.hh"

#include <G4LogicalVolume.hh>

class DecayTrapConstruction: public MaterialUser
{
public:
  DecayTrapConstruction();	// constructor

  G4double dWindowThick;
  G4double dCoatingThick;
  G4double dInnerRadiusOfTrap;
  G4double dTubeWallThickness;
  G4double dInnerRadiusOfCollimator;
  G4Material* mTubeMat;
  G4Material* mCollimatorMat;
  G4Material* mDecayTrapWindowMat;
  G4Material* mDecayTrapCoatingMat;

  G4LogicalVolume* decayTrapTube_log; 	///< decay trap tube
  G4LogicalVolume* decayTrapWin_log[2]; 	///< trap window volume
  G4LogicalVolume* mylarWin_log[2]; 	///< mylar layer of window
  G4LogicalVolume* beWin_log[2]; 	///< berillium layer of window
  G4LogicalVolume* trapMonitor_log[2];	///< extra event monitoring region
  G4LogicalVolume* collimator_log[2]; 	///< collimator
  G4LogicalVolume* collimatorBack_log[2]; ///< bracket behind collimator

  //WiggleSheet wigglefoils[2]; ///< optional replacement crinkly endcap foils

  // construct is passed-to world volume, since we are not using an overall trap "container" vol.
  void Build(G4LogicalVolume* world, float crinkleAngle = 0);

};

#endif
