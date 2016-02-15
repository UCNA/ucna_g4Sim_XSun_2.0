#include "DecayTrapConstruction.hh"

#include <G4PVPlacement.hh>
#include <G4Tubs.hh>
#include <G4RotationMatrix.hh>


DecayTrapConstruction::DecayTrapConstruction()
: dWindowThick(0.7*um), dCoatingThick(0.3*um),
dInnerRadiusOfTrap(2.45*inch), dTubeWallThickness(2*mm), dInnerRadiusOfCollimator(2.3*inch), dCollimatorThick(0.8*inch),
mTubeMat(Cu), mCollimatorMat(Polyethylene), mDecayTrapWindowMat(Mylar), mDecayTrapCoatingMat(Be)
{
  // Can use this area to set the default class members as I have just done above.
}

void DecayTrapConstruction::Build(G4LogicalVolume* world, float crinkleAngle)
{
  G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));

  // decay tube construction
  G4double tubeOuterRadius = dInnerRadiusOfTrap + dTubeWallThickness;
  G4double tubeLength = 3.0*m;

  G4Tubs* tube = new G4Tubs("decayTrap_tube", dInnerRadiusOfTrap, tubeOuterRadius, tubeLength/2., 0., 2*M_PI);
  decayTrapTube_log = new G4LogicalVolume(tube, mTubeMat, "decayTrap_tube_log");
  decayTrapTube_log -> SetVisAttributes(new G4VisAttributes(G4Colour(1,1,0,0.5)));
  new G4PVPlacement(NULL, G4ThreeVector(), decayTrapTube_log, "decayTrap_tube", world, false, 0, true);

  // decay trap windows, collimator, monitors
  G4double totalWindowThickness = dWindowThick + dCoatingThick;
//  G4double decayTrap_collimatorThick = 0.8*inch;	// original my sim coding
//  G4double collimatorThick = 0.7*inch;	// M.Brown's 2011/2012 geometry
  G4double beWindow_PosZ = -totalWindowThickness/2. + dCoatingThick/2.;
  G4double mylarWindow_PosZ = totalWindowThickness/2. - dWindowThick/2.;
  G4double window_PosZ = (tubeLength + totalWindowThickness)/2.;
  G4double monitorThickness = 1.0*mm;
  G4double monitor_PosZ = 0.5*m;

  G4Tubs* trapWindowTube = new G4Tubs("trap_win_tube", 0., tubeOuterRadius, totalWindowThickness/2., 0, 2*M_PI);
  G4Tubs* mylarTube = new G4Tubs("mylarTube", 0., tubeOuterRadius, dWindowThick/2., 0., 2*M_PI);
  G4Tubs* beTube = new G4Tubs("beTube", 0., tubeOuterRadius, dCoatingThick/2., 0., 2*M_PI);
  G4Tubs* collimatorTube = new G4Tubs("decayTrap_collimatorTube", dInnerRadiusOfCollimator,
                                dInnerRadiusOfCollimator + dCollimatorThick, dCollimatorThick/2.,0., 2*M_PI);
  G4Tubs* collimatorBackTube = new G4Tubs("decayTrap_collimatorBackTube", tubeOuterRadius + 1.*mm,
                                dInnerRadiusOfCollimator + dCollimatorThick, dCollimatorThick/2., 0., 2*M_PI);
  G4Tubs* trapMonitorTube = new G4Tubs("trap_monitor_tube", 0., dInnerRadiusOfTrap, monitorThickness/2.0, 0., 2*M_PI);

  for(int i = 0; i <= 1; i++)
  {
    decayTrapWin_log[i] = new G4LogicalVolume(trapWindowTube, Vacuum, Append(i, "trap_win_log_"));
    decayTrapWin_log[i] -> SetVisAttributes(visWindow);
    mylarWin_log[i] = new G4LogicalVolume(mylarTube, mDecayTrapWindowMat, Append(i, "mylar_win_log_"));
    beWin_log[i] = new G4LogicalVolume(beTube, mDecayTrapCoatingMat, Append(i, "be_win_log"));
  }
  new G4PVPlacement(NULL, G4ThreeVector(0.,0.,(-1)*window_PosZ), decayTrapWin_log[0], "trap_win_EAST", world, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, window_PosZ), decayTrapWin_log[1], "trap_win_WEST", world, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., (-1)*mylarWindow_PosZ), mylarWin_log[0], "mylar_win_EAST", decayTrapWin_log[0], false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., (-1)*beWindow_PosZ), beWin_log[0], "be_win_EAST", decayTrapWin_log[0], false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., mylarWindow_PosZ), mylarWin_log[1], "mylar_win_WEST", decayTrapWin_log[1], false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0., 0., beWindow_PosZ), beWin_log[1], "be_win_WEST", decayTrapWin_log[1], false, 0);

  G4double collimator_PosZ = (tubeLength + dCollimatorThick)/2.;
  collimator_PosZ += totalWindowThickness/2.;
  G4double collimatorBack_PosZ = tubeLength/2. - dCollimatorThick;

  for(int i = 0; i <= 1; i++)
  {
    collimator_log[i] = new G4LogicalVolume(collimatorTube, mCollimatorMat, Append(i, "collimator_log_"));
    collimatorBack_log[i] = new G4LogicalVolume(collimatorBackTube, mCollimatorMat, Append(i, "collimator_back_log_"));
    trapMonitor_log[i] = new G4LogicalVolume(trapMonitorTube, Vacuum, Append(i, "trap_monitor_log_"));
  }

  // place everything at -z i.e. EAST.
  new G4PVPlacement(NULL, G4ThreeVector(0, 0, (-1)*collimator_PosZ), collimator_log[0],
                        "collimator_EAST", world, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, (-1)*collimatorBack_PosZ), collimatorBack_log[0],
                        "collimator_back_EAST", world, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, (-1)*monitor_PosZ), trapMonitor_log[0],
                        "trap_monitor_EAST", world, false, 0);
  // copy but place at +z i.e. WEST
  new G4PVPlacement(NULL, G4ThreeVector(0, 0, collimator_PosZ), collimator_log[1],
                        "collimator_WEST", world, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, collimatorBack_PosZ), collimatorBack_log[1],
                        "collimator_back_WEST", world, false, 0);
  new G4PVPlacement(NULL, G4ThreeVector(0,0, monitor_PosZ), trapMonitor_log[1],
                        "trap_monitor_WEST", world, false, 0);

}
