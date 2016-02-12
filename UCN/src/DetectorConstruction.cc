#include "DetectorConstruction.hh"
#include "GlobalField.hh"	// only need field in the .cc file
#include "MWPCField.hh"		// all the rest of the geometry goes in .hh file

#include <G4UserLimits.hh>		// stole from Michael Mendenhall's code.

#include <Randomize.hh>			// Stolen from Analysis Manager
#include <G4ios.hh>			// Pretty sure needed for TrackerSD
#include <G4Run.hh>			// Leave them here since we use registerSD in DetectorConstruction
#include <G4Event.hh>			// And the registerSD is totally not working without it
#include <G4Track.hh>
#include <G4VVisManager.hh>
#include <G4TrajectoryContainer.hh>
#include <G4Trajectory.hh>
#include <G4IonTable.hh>
#include <G4SDManager.hh>
#include <G4PrimaryVertex.hh>
#include <G4PrimaryParticle.hh>
#include <G4SDManager.hh>
#include <G4EventManager.hh>

#include <G4EqMagElectricField.hh>
#include <G4MagneticField.hh>		// Bottom half of detector construction
#include <G4FieldManager.hh>
#include <G4ChordFinder.hh>
#include <G4PropagatorInField.hh>
#include <G4TransportationManager.hh>
#include <G4UserLimits.hh>
#include <G4PVParameterised.hh>

#include <G4ClassicalRK4.hh>
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4SimpleHeum.hh"
#include "G4HelixHeum.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4HelixMixedStepper.hh"


DetectorConstruction::DetectorConstruction()
: G4VUserDetectorConstruction()
{
  // initialize some useful private class variables
  bUseSourceHolder = false;
  fScintStepLimit = 1.0*mm;
  fStorageIndex = 0;	// this loops over our TrackerHit names storage array

  fCrinkleAngle = 0;
}


DetectorConstruction::~DetectorConstruction()
{ }

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // ----- Materials gets made before here via the DetectorTools class
  SetVacuumPressure(0);	// this is the set vacuum pressure from DetectorTools

  // user step limits
  G4UserLimits* UserCoarseLimits = new G4UserLimits();
  UserCoarseLimits->SetMaxAllowedStep(10*m);
  G4UserLimits* UserGasLimits = new G4UserLimits();
  UserGasLimits->SetMaxAllowedStep(1*cm);
  G4UserLimits* UserSolidLimits = new G4UserLimits();
  UserSolidLimits->SetMaxAllowedStep(fScintStepLimit);	// default value from Messenger class.

  //----- Experimental Hall. World volume.
  G4double expHall_x = 2.0*m;
  G4double expHall_y = 2.0*m;
  G4double expHall_z = 8.0*m;
  G4Box* experimentalHall_box = new G4Box("expHall_box", expHall_x/2, expHall_y/2, expHall_z/2);
  experimentalHall_log = new G4LogicalVolume(experimentalHall_box, Vacuum, "World_log");
  experimentalHall_log -> SetVisAttributes(G4VisAttributes::Invisible);
  experimentalHall_log -> SetUserLimits(UserCoarseLimits);
  experimentalHall_phys = new G4PVPlacement(NULL, G4ThreeVector(), "World_phys", experimentalHall_log, 0, false, 0);

  //----- Source holder object. Used if it is a calibration source.
  // variables related to defaults for Source Holder. Integrate with messenger class later.
  vSourceHolderPos = G4ThreeVector(0,0,0);	// specifically, these class members should be set earlier
  fSourceFoilThick = 7.2*um;	// default value set in detector construction messenger class
  Source.dSourceWindowThickness = fSourceFoilThick/2.;
  Source.Build();
  // place entire source holder object. Comment these two lines out if don't want to use source holder object
  source_phys = new G4PVPlacement(NULL, vSourceHolderPos, Source.sourceContainer_log, "source_container_phys",
				experimentalHall_log, false, 0, false);	// explicitly don't check overlaps

  //----- Decay Trap tube (length 3m, main tube)
  // these guys are needed for the wirechamber active region
  G4double wireVol_anodeRadius = 5*um;
  G4double wireVol_cathodeRadius = 39.1*um;

  Trap.Build(experimentalHall_log, fCrinkleAngle);

/*  // Michael Brown's changes that form the 2011/2012 detector geometry
  // All the appropriate variables have been commented out everywhere else.
  // Note: some of these values take the same value as previously declared (commented out)
  G4double decayTrap_windowThick = 0.500*um;
  G4Material* decayTrap_windowMaterial = Mylar;
  G4double decayTrap_innerRadiusCollimator = 2.3*inch;
  G4double decayTrap_collimatorThick = 0.7*inch;
  G4double wireVol_anodeRadius = 5*um;
  G4double wireVol_cathodeRadius = 39.1*um;
*/  // to change back, you'll need to find EXACTLY where these are all commented out.

  //----- Scintillator construction. Used as Sensitive Volume
  for(int j = 0; j <= 1; j++)
  {
    ScintCrystal[j].Build(j);
  }

//  G4ThreeVector sideTransScintEast = G4ThreeVector(0., 0., (-1)*(2.2*m - scint_face_PosZ));
//  G4ThreeVector sideTransScintWest = G4ThreeVector(0., 0., (2.2*m - scint_face_PosZ));
  G4ThreeVector sideTransScintEast = G4ThreeVector(0., 0., (-1)*(2.2*m - ScintCrystal[0].GetScintFacePos()));
  G4ThreeVector sideTransScintWest = G4ThreeVector(0., 0., (2.2*m - ScintCrystal[1].GetScintFacePos()));
  G4RotationMatrix* EastSideRot = new G4RotationMatrix();
  EastSideRot -> rotateY(M_PI*rad);


  // Place the two scintillator containers. Rotate the EAST one so they both face towards the center.
  scint_phys[0] = new G4PVPlacement(EastSideRot, sideTransScintEast, ScintCrystal[0].container_log,
				"scint_container_phys_EAST", experimentalHall_log, false, 0, true);
  scint_phys[1] = new G4PVPlacement(NULL, sideTransScintWest, ScintCrystal[1].container_log,
				"scint_container_phys_WEST", experimentalHall_log, false, 0, true);

  //----- Begin wirechamber construction. MWPC used in front of Scintillator.
/*  G4double mwpc_windowThick = 6*um;
  G4double mwpc_entranceRadius = 7.0*cm;
  G4double mwpc_exitRadius = 7.5*cm;
  G4double mwpc_entranceToCathodes = 5.0*mm;
  G4double mwpc_exitToCathodes = 5.0*mm;
  G4double mwpc_fieldE0 = 2700*volt;		// gets passed to MWPC fields and used in SetPotential
  G4Material* mwpc_fillGas = wireVol_activeGas;	// want it to be WCPentane

  G4double mwpc_containerHalf_Z = 0.5*(mwpc_entranceToCathodes + mwpc_exitToCathodes + 2*cm);
  G4double mwpc_gasVolumeWidth = 8.0*inch;	// MWPC gas box width

  // container volume for all MWPC
  G4Box* mwpc_containerBox = new G4Box("mwpc_container_box", mwpc_gasVolumeWidth/2., mwpc_gasVolumeWidth/2., mwpc_containerHalf_Z);

  mwpc_container_log[0] = new G4LogicalVolume(mwpc_containerBox, mwpc_fillGas, "mwpc_container_log_EAST");
  mwpc_container_log[0] -> SetVisAttributes(G4VisAttributes::Invisible);
  mwpc_container_log[1] = new G4LogicalVolume(mwpc_containerBox, mwpc_fillGas, "mwpc_container_log_WEST");
  mwpc_container_log[1] -> SetVisAttributes(G4VisAttributes::Invisible);
//  mwpc_container_log[0] -> SetVisAttributes(new G4VisAttributes(G4Colour(1,0,1,1.0)));
//  mwpc_container_log[1] -> SetVisAttributes(new G4VisAttributes(G4Colour(1,0,1,1.0)));

  // MWPC active gas volume placement with  wireplane, relative to MWPC container volume
  G4ThreeVector mwpc_activeRegionTrans(0, 0, (mwpc_entranceToCathodes - mwpc_exitToCathodes)/2.);
  // Note: these lines place wireVol inside mwpc.
  new G4PVPlacement(NULL, mwpc_activeRegionTrans, wireVol_gas_log[0], "mwpc_activeReg_phys_EAST", mwpc_container_log[0], false, 0);
  new G4PVPlacement(NULL, mwpc_activeRegionTrans, wireVol_gas_log[1], "mwpc_activeReg_phys_WEST", mwpc_container_log[1], false, 0);

  // construct kevlay string. Rectangular cross section strings with equal volume to nominal 140um cylinders.
  G4double mwpc_kevRadius = 0.07*mm;
  G4double mwpc_kevSpacing = 5.0*mm;
  G4int mwpc_NbKevWires = 32;
  G4double mwpc_kevLength = 15.0*cm;
  double mwpc_kevAspectRatio = 16.0;	// aspect ratio, width:depth.
  G4double mwpc_kevArea = M_PI*mwpc_kevRadius*mwpc_kevRadius;
  G4double mwpc_kevEffWidth = sqrt(mwpc_kevArea*mwpc_kevAspectRatio);
  G4double mwpc_kevEffThick = sqrt(mwpc_kevArea/mwpc_kevAspectRatio);
  G4double mwpc_kev_PosZ = -mwpc_containerHalf_Z + mwpc_kevEffThick/2.;

  // create shapes for kevlar strings. Will be replicated into an array
  G4Box* mwpc_kevContainerBox = new G4Box("kevContainer_box", mwpc_NbKevWires*mwpc_kevSpacing/2., mwpc_kevLength/2., mwpc_kevEffThick/2.);
  G4Box* mwpc_kevSegBox = new G4Box("kevSeg_box", mwpc_kevSpacing/2., mwpc_kevLength/2., mwpc_kevEffThick/2.);
  G4Box* mwpc_kevStripBox = new G4Box("kevStrip_box", mwpc_kevEffWidth/2., mwpc_kevLength/2., mwpc_kevEffThick/2.);

  // Shapes for the Mylar windows in the MWPC
  G4Tubs* mwpc_winInnerTube = new G4Tubs("winInnerTube", 0., mwpc_entranceRadius, mwpc_windowThick/2., 0., 2*M_PI);
  G4Tubs* mwpc_winOuterTube = new G4Tubs("winOuterTube", 0., mwpc_exitRadius, mwpc_windowThick/2., 0., 2*M_PI);

  for(int i = 0; i <= 1; i++)
  {
    mwpc_kevContainer_log[i] = new G4LogicalVolume(mwpc_kevContainerBox, Vacuum, Append(i, "kevContainer_log_"));
    mwpc_kevSeg_log[i] = new G4LogicalVolume(mwpc_kevSegBox, Vacuum, Append(i, "kevSeg_log_"));
    mwpc_kevStrip_log[i] = new G4LogicalVolume(mwpc_kevStripBox, Kevlar, Append(i, "kevStrip_log_"));

    mwpc_winIn_log[i] = new G4LogicalVolume(mwpc_winInnerTube, Mylar, Append(i, "winIn_log_"));
    mwpc_winIn_log[i] -> SetVisAttributes(visWindow);
    mwpc_winOut_log[i] = new G4LogicalVolume(mwpc_winOuterTube, Mylar, Append(i, "winOut_log_"));
    mwpc_winOut_log[i] -> SetVisAttributes(visWindow);

    new G4PVPlacement(NULL, G4ThreeVector(0,0, mwpc_kev_PosZ), mwpc_kevContainer_log[i], Append(i, "kevContainer_phys_"),
			mwpc_container_log[i], false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0,0,0), mwpc_kevStrip_log[i], Append(i, "kevStrip_phys_"), mwpc_kevSeg_log[i], false, 0);
    new G4PVReplica(Append(i, "kevlar_plane_"), mwpc_kevSeg_log[i], mwpc_kevContainer_log[i], kXAxis, mwpc_NbKevWires, mwpc_kevSpacing);

    new G4PVPlacement(NULL, G4ThreeVector(0,0, -mwpc_containerHalf_Z + mwpc_kevEffThick + mwpc_windowThick/2.),
        mwpc_winIn_log[i], Append(i, "winIn_phys_"), mwpc_container_log[i], false, 0);
    new G4PVPlacement(NULL, G4ThreeVector(0,0, mwpc_containerHalf_Z - mwpc_windowThick/2.),
        mwpc_winOut_log[i], Append(i, "winOut_phys_"), mwpc_container_log[i], false, 0);
  }

  G4double frame_backWinFrameThick = 0.5*inch;	// originally placed further down but needed here
  G4double mwpc_PosZ = -mwpc_containerHalf_Z - frame_backWinFrameThick - (scint_N2Volume_Z/2. + scint_face_PosZ);
  G4ThreeVector sideTransMWPCEast = G4ThreeVector(0,0, (-1)*(2.2*m + mwpc_PosZ));
  G4ThreeVector sideTransMWPCWest = G4ThreeVector(0, 0, (2.2*m + mwpc_PosZ));

  // place the two wire chambers in experimentalHall. Rotate the East one (same as scint) so they both point towards center.
  mwpc_container_phys[0] = new G4PVPlacement(EastSideRot, sideTransMWPCEast, mwpc_container_log[0],
				"mwpc_container_phys_EAST", experimentalHall_log, false, 0, true);
  mwpc_container_phys[1] = new G4PVPlacement(NULL, sideTransMWPCWest, mwpc_container_log[1],
				"mwpc_container_phys_West", experimentalHall_log, false, 0, true);
*/

  for(int i = 0; i <= 1; i++)
  {
    Wirechamber[i].Build(i);	// creates Wirechamber objects and puts them in 2 element array to hold
  }

  mwpc_phys[0] = new G4PVPlacement(EastSideRot, G4ThreeVector(0, 0, -2*m), Wirechamber[0].container_log,
                                "mwpc_container_phys_EAST", experimentalHall_log, false, 0, true);
  mwpc_phys[1] = new G4PVPlacement(NULL, G4ThreeVector(0,0,2*m), Wirechamber[1].container_log,
                                "mwpc_container_phys_West", experimentalHall_log, false, 0, true);






  //----- Begin DetectorPackageConstruction. This is the frame that holds the scintillator and MWPC.
/*
  G4double frame_packageRadius = 6.0*inch;
  G4double frame_mwpcEntranceThick = 0.375*inch;
  G4double frame_mwpcEntranceRadius = 3.0*inch;	// not the same value as mwpc_entranceRadius
  G4double frame_mwpcEntranceDepth = 5.0*inch;
  G4double frame_frontWinFrameThick = 1.0*inch;
//  G4double frame_backWinFrameThick = 0.5*inch;	// moved a few lines up since needed for mwpc placement

  // create the shapes that will be used for geometric objects below
  // aluminum entrance collimator to detector package
  G4double frame_entranceSectionLength = frame_mwpcEntranceDepth + frame_frontWinFrameThick;
  G4Tubs* frame_mwpcEntranceTube = new G4Tubs("mwpc_entrance_tube", 0., frame_packageRadius, 0.5*frame_entranceSectionLength, 0., 2*M_PI);
  G4Tubs* frame_entranceFrontTube = new G4Tubs("entrance_front_tube", frame_mwpcEntranceRadius + frame_mwpcEntranceThick,
					frame_packageRadius, 0.5*frame_mwpcEntranceThick, 0., 2*M_PI);
  G4Tubs* frame_entranceMidTube = new G4Tubs("entrance_mid_tube", frame_mwpcEntranceRadius, frame_mwpcEntranceRadius + frame_mwpcEntranceThick,
					0.5*frame_mwpcEntranceDepth, 0., 2*M_PI);
  G4Tubs* frame_entranceBackTube = new G4Tubs("entrance_back_tube", mwpc_entranceRadius, frame_packageRadius,
					0.5*frame_frontWinFrameThick, 0., 2*M_PI);
  G4VisAttributes* visMWPCEntrance = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.8));

  // create overall detector package frame tube.
  G4double frame_detFrameHalf_Z = frame_mwpcEntranceDepth + 2*mwpc_containerHalf_Z + 1.0*inch;
  G4Tubs* frame_framePackageTube = new G4Tubs("detPackage_tube_EAST", 0, frame_packageRadius, frame_detFrameHalf_Z, 0, 2*M_PI);

  // subtract off the scintillator container volume. Not rotated (since frame will be made and then rotated).
  // But displaced by -scint_face_PosZ relative to the local coordinates of the detector package frame.
  G4SubtractionSolid* frame_frameMinusScint = new G4SubtractionSolid("DPC_frame_minus_scint_container_log", frame_framePackageTube,
					scint_N2VolTube, NULL, G4ThreeVector(0., 0., -scint_face_PosZ));
  // subtract off the mwpc container volume. Not rotated (since frame will be made then rotated).
  G4SubtractionSolid* frame_containerShape = new G4SubtractionSolid("frame_container_minus_Scint_MWPC", frame_frameMinusScint,
					mwpc_containerBox, NULL, G4ThreeVector(0., 0., mwpc_PosZ));

  // "place componenets relative to scintillator face at 0" is what old code says.
  // Already successfully placed the scint and the mwpc in experimentalHall_log. This code is remaining geom that wasn't done
  G4double frame_entrance_PosZ = mwpc_PosZ - (2*mwpc_containerHalf_Z + frame_entranceSectionLength)/2.;

  // aluminum exit window and N2 volume at back of gas box
  G4Tubs* frame_mwpcExitTube = new G4Tubs("mwpc_exit_tube", mwpc_exitRadius, frame_packageRadius, 0.5*frame_backWinFrameThick, 0., 2*M_PI);
  G4VisAttributes* visMWPCExit = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.8));

  G4double frame_exitWin_PosZ = mwpc_PosZ + (2*mwpc_containerHalf_Z + frame_backWinFrameThick)/2.;
  G4Tubs* frame_mwpcExitGasN2Tube = new G4Tubs("mwpc_exit_N2_tube", 0, mwpc_exitRadius, 0.5*frame_backWinFrameThick, 0., 2*M_PI);

  // material behind the detector. Misc stuff that can cause back scattering events.
  G4double frame_backStuffThick = 1.0*inch;
  G4Tubs* frame_backStuffTube = new G4Tubs("backstuff_tube_EAST", 0, 0.5*frame_packageRadius, frame_backStuffThick, 0., 2*M_PI);

  // create logical volumes from all the shapes defined above. As always, array of 2 for East/West (0/1)
  for(int i = 0; i <= 1; i++)
  {
    frame_mwpcEntrance_log[i] = new G4LogicalVolume(frame_mwpcEntranceTube, Vacuum, Append(i, "mwpc_entrance_log_"));
    frame_entranceFront_log[i] = new G4LogicalVolume(frame_entranceFrontTube, Al, Append(i, "entrance_front_log_"));
    frame_entranceMid_log[i] = new G4LogicalVolume(frame_entranceMidTube, Al, Append(i, "entrance_mid_log_"));
    frame_entranceBack_log[i] = new G4LogicalVolume(frame_entranceBackTube, Al, Append(i, "entrance_back_log_"));
    frame_mwpcEntrance_log[i] -> SetVisAttributes(G4VisAttributes::Invisible);
    frame_entranceFront_log[i] -> SetVisAttributes(visMWPCEntrance);
    frame_entranceMid_log[i] -> SetVisAttributes(visMWPCEntrance);
    frame_entranceBack_log[i] -> SetVisAttributes(visMWPCEntrance);

    frame_container_log[i] = new G4LogicalVolume(frame_containerShape, Vacuum, Append(i, "frame_container_log_"));
    frame_container_log[i] -> SetVisAttributes(G4VisAttributes::Invisible);
//    frame_container_log[i] -> SetVisAttributes(new G4VisAttributes(G4Color(1, 0, 1, 1)));

    frame_mwpcExit_log[i] = new G4LogicalVolume(frame_mwpcExitTube, Al, Append(i, "mwpc_exit_log_"));
    frame_mwpcExit_log[i] -> SetVisAttributes(visMWPCExit);
    frame_mwpcExitGasN2_log[i] = new G4LogicalVolume(frame_mwpcExitGasN2Tube, WCNitrogen, Append(i, "mwpc_exit_N2_log_"));
    frame_mwpcExitGasN2_log[i] -> SetVisAttributes(G4VisAttributes::Invisible);

    frame_backStuff_log[i] = new G4LogicalVolume(frame_backStuffTube, SS304, Append(i, "backStuff_log_"));
  }

  // place all the physical volumes (class members) in their respective mother volumes.
  for(int i = 0; i <= 1; i++)
  {  // These three are placed in the logical vol. mwpcEntrance.
    frame_entranceFront_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*(frame_entranceSectionLength - frame_mwpcEntranceThick)),
                                  frame_entranceFront_log[i], Append(i, "entrance_front_phys_"), frame_mwpcEntrance_log[i], false, 0);
    frame_entranceMid_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*frame_frontWinFrameThick),
                                  frame_entranceMid_log[i], Append(i, "entrance_mid_phys_"), frame_mwpcEntrance_log[i], false, 0);
    frame_entranceBack_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, 0.5*(frame_entranceSectionLength - frame_frontWinFrameThick)),
                                  frame_entranceBack_log[i], Append(i, "entrance_back_phys_"), frame_mwpcEntrance_log[i], false, 0);

    // All this stuff makes up the frame "container_log" i.e. the overall frame logical.
    frame_mwpcEntrance_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0., 0., frame_entrance_PosZ),
                                  frame_mwpcEntrance_log[i], Append(i, "frame_mwpc_entrance_"), frame_container_log[i], false, 0);
    frame_mwpcExit_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, frame_exitWin_PosZ), frame_mwpcExit_log[i],
                                  Append(i, "mwpc_exit_"), frame_container_log[i], false, 0);
    frame_mwpcExitGasN2_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0,0, frame_exitWin_PosZ), frame_mwpcExitGasN2_log[i],
                                  Append(i, "mwpc_exit_N2_phys_"), frame_container_log[i], false, 0);
    frame_backStuff_phys[i] = new G4PVPlacement(NULL, G4ThreeVector(0, 0, frame_detFrameHalf_Z - 0.5*frame_backStuffThick),
                                  frame_backStuff_log[i], Append(i, "backStuff_phys_"), frame_container_log[i], false, 0);
  }
*/
  //----- Finish up detector construction. Need to place frame_container_log in experimentalHall
/*  G4ThreeVector frameTransEast = G4ThreeVector(0., 0., (-1)*(2.2*m));	// note: scint face position is 0 in local coord.
									// Also there's no offset. So it's just -2.2m
  G4ThreeVector frameTransWest = G4ThreeVector(0., 0., 2.2*m);

  frame_container_phys[0] = new G4PVPlacement(EastSideRot, frameTransEast, frame_container_log[0],
				"Detector_Package_Frame_EAST", experimentalHall_log, false, 0, true);
  frame_container_phys[1] = new G4PVPlacement(NULL, frameTransWest, frame_container_log[1],
				"Detector_Package_Frame_WEST", experimentalHall_log, false, 0, true);
*/
  for(int i = 0; i <= 1; i++)			// set user limits in specific volumes
  {
//    decayTrap_window_log[i] -> SetUserLimits(UserSolidLimits);
/*    mwpc_container_log[i] -> SetUserLimits(UserGasLimits);
    mwpc_winIn_log[i] -> SetUserLimits(UserSolidLimits);
    mwpc_winOut_log[i] -> SetUserLimits(UserSolidLimits);
    mwpc_kevStrip_log[i] -> SetUserLimits(UserSolidLimits);
*/
  }

  // register logical volumes as sensitive detectors. Used for all info tracking during sim
  for(int i = 0; i <= 1; i++)
  {
    SD_scint_scintillator[i] = RegisterSD(Append(i, "SD_scint_"), Append(i, "HC_scint_"));
    ScintCrystal[i].scintillator_log -> SetSensitiveDetector(SD_scint_scintillator[i]);
/*    SD_scint_scintillator[i] = RegisterSD(Append(i, "SD_scint_"), Append(i, "HC_scint_"));
    scint_scintillator_log[i] -> SetSensitiveDetector(SD_scint_scintillator[i]);
*/
/*    SD_scint_deadScint[i] = RegisterSD(Append(i, "SD_deadScint_"));
    scint_deadLayer_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);
    scint_container_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);
    frame_mwpcExitGasN2_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);	// include N2 vol here
    scint_lightGuide_log[i] -> SetSensitiveDetector(SD_scint_deadScint[i]);	// and also light guides

    SD_scint_backing[i] = RegisterSD(Append(i, "SD_scint_backing_"));
    scint_backing_log[i] -> SetSensitiveDetector(SD_scint_backing[i]);

    SD_mwpc_winIn[i] = RegisterSD(Append(i, "SD_mwpc_winIn_"));
    mwpc_winIn_log[i] -> SetSensitiveDetector(SD_mwpc_winIn[i]);

    SD_mwpc_winOut[i] = RegisterSD(Append(i, "SD_mwpc_winOut_"));
    mwpc_winOut_log[i] -> SetSensitiveDetector(SD_mwpc_winOut[i]);

    SD_decayTrap_windows[i] = RegisterSD(Append(i, "SD_decayTrap_windows_"));
    decayTrap_mylarWindow_log[i] -> SetSensitiveDetector(SD_decayTrap_windows[i]);
    decayTrap_beWindow_log[i] -> SetSensitiveDetector(SD_decayTrap_windows[i]);
*/
    SD_wireVol[i] = RegisterSD(Append(i, "SD_wireVol_"), Append(i, "HC_wireVol_"));
    Wirechamber[i].ActiveRegion.gas_log -> SetSensitiveDetector(SD_wireVol[i]);
    Wirechamber[i].ActiveRegion.anodeSeg_log -> SetSensitiveDetector(SD_wireVol[i]);
    Wirechamber[i].ActiveRegion.cathSeg_log -> SetSensitiveDetector(SD_wireVol[i]);

/*    SD_wireVol[i] = RegisterSD(Append(i, "SD_wireVol_"), Append(i, "HC_wireVol_"));
    wireVol_gas_log[i] -> SetSensitiveDetector(SD_wireVol[i]);
    wireVol_anodeSeg_log[i] -> SetSensitiveDetector(SD_wireVol[i]);
    wireVol_cathSeg_log[i] -> SetSensitiveDetector(SD_wireVol[i]);
*/
/*
    SD_wireVol_planes[i] = RegisterSD(Append(i, "SD_wireVol_planes_"));
    wireVol_cathodeWire_log[i] -> SetSensitiveDetector(SD_wireVol_planes[i]);
    wireVol_cathPlate_log[i] -> SetSensitiveDetector(SD_wireVol_planes[i]);
    wireVol_anodeWire_log[i] -> SetSensitiveDetector(SD_wireVol_planes[i]);

    SD_mwpc_container[i] = RegisterSD(Append(i, "SD_mwpc_container_"));	// equivalently, dead region in mwpc
    mwpc_container_log[i] -> SetSensitiveDetector(SD_mwpc_container[i]);

    SD_mwpc_kevStrip[i] = RegisterSD(Append(i, "SD_mwpc_kevStrip_"));
    mwpc_kevStrip_log[i] -> SetSensitiveDetector(SD_mwpc_kevStrip[i]);
*/
  }

/*  SD_source = RegisterSD("SD_source");
  source_window_log -> SetSensitiveDetector(SD_source);
  for(int i = 0; i <= 1; i++)
    source_coating_log[i] -> SetSensitiveDetector(SD_source);

  for(int i = 0; i <= 1; i++)
  {
    SD_decayTrap_innerMonitors[i] = RegisterSD(Append(i, "SD_decayTrap_innerMonitors_"));
    decayTrap_innerMonitors_log[i] -> SetSensitiveDetector(SD_decayTrap_innerMonitors[i]);
  }

  // exp hall vacuum, decay tube, other inert parts. Stores all energy "lost".
  SD_world = RegisterSD("SD_world");
  experimentalHall_log -> SetSensitiveDetector(SD_world);
  decayTrap_tube_log -> SetSensitiveDetector(SD_world);
  for(int i = 0; i <= 1; i++)
  {
    frame_mwpcEntrance_log[i] -> SetSensitiveDetector(SD_world);
    frame_mwpcExit_log[i] -> SetSensitiveDetector(SD_world);	// confusing, since made of Al. But trust M.M.
    frame_container_log[i] -> SetSensitiveDetector(SD_world);
    decayTrap_collimator_log[i] -> SetSensitiveDetector(SD_world);	// and this is polyethylene?
    decayTrap_collimatorBack_log[i] -> SetSensitiveDetector(SD_world);
  }
*/


  // Create everything needed for global and local EM fields
/*  G4ThreeVector East_EMFieldLocation = mwpc_activeRegionTrans + sideTransMWPCEast;
  G4ThreeVector West_EMFieldLocation = mwpc_activeRegionTrans + sideTransMWPCWest;

  ConstructGlobalField();			// make magnetic and EM fields.
  ConstructEastMWPCField(wireVol_wireSpacing, wireVol_planeSpacing, wireVol_anodeRadius,
			mwpc_fieldE0, EastSideRot, East_EMFieldLocation);
  ConstructWestMWPCField(wireVol_wireSpacing, wireVol_planeSpacing, wireVol_anodeRadius,
			mwpc_fieldE0, NULL, West_EMFieldLocation);
*/

  return experimentalHall_phys;
}

TrackerSD* DetectorConstruction::RegisterSD(G4String sdName, G4String hcName)
{
  TrackerSD* sd = new TrackerSD(sdName, hcName);
  G4SDManager::GetSDMpointer() -> AddNewDetector(sd);

  fSDNamesArray[fStorageIndex] = sdName;
  fHCNamesArray[fStorageIndex] = hcName;
  fStorageIndex++;

  return sd;
}

void DetectorConstruction::ConstructGlobalField()
{
  G4cout << "Setting up global magnetic field. Call to global field object." << G4endl;

  GlobalField* magField = new GlobalField();
  G4FieldManager* globalFieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  globalFieldManager -> SetDetectorField(magField);
  globalFieldManager -> CreateChordFinder(magField);

  G4MagIntegratorStepper* pStepper;
  G4Mag_UsualEqRhs* equationOfMotion = new G4Mag_UsualEqRhs(magField);
  //pStepper = new G4ClassicalRK4 (fEquation); // general case for "smooth" EM fields
  //pStepper = new G4SimpleHeum( fEquation ); // for slightly less smooth EM fields
  //pStepper = new G4HelixHeum( fEquation ); // for "smooth" pure-B fields
  //pStepper = new G4HelixImplicitEuler( fEquation ); // for less smooth pure-B fields; appears ~50% faster than above
  //pStepper = new G4HelixSimpleRunge( fEquation ); // similar speed to above
  //pStepper = new G4HelixExplicitEuler( fEquation ); // about twice as fast as above
  pStepper = new G4HelixMixedStepper(equationOfMotion,6); // avoids "Stepsize underflow in Stepper" errors
  globalFieldManager -> GetChordFinder() -> GetIntegrationDriver() -> RenewStepperAndAdjust(pStepper);

  globalFieldManager -> GetChordFinder() -> SetDeltaChord(100.0*um);
  globalFieldManager -> SetMinimumEpsilonStep(1e-6);
  globalFieldManager -> SetMaximumEpsilonStep(1e-5);
  globalFieldManager -> SetDeltaOneStep(0.1*um);
  G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetMaxLoopCount(INT_MAX);

  return;
}

void DetectorConstruction::ConstructEastMWPCField(G4double a, G4double b, G4double c, G4double d, G4RotationMatrix* e, G4ThreeVector f)
{
  G4cout << "Setting up East wirechamber electromagnetic field." << G4endl;
  MWPCField* eastLocalField = new MWPCField();
  eastLocalField -> SetActiveReg_d(a);
  eastLocalField -> SetActiveReg_L(b);
  eastLocalField -> SetActiveReg_r(c);
  eastLocalField -> SetSideRot(e);
  eastLocalField -> SetSideTrans(f);
  eastLocalField -> SetPotential(d);

  G4FieldManager* eastLocalFieldManager = new G4FieldManager();
  eastLocalFieldManager -> SetDetectorField(eastLocalField);

  G4EqMagElectricField* eastlocalEquation = new G4EqMagElectricField(eastLocalField);
  G4ClassicalRK4* eastlocalStepper = new G4ClassicalRK4(eastlocalEquation,8);
  G4MagInt_Driver* eastlocalIntgrDriver = new G4MagInt_Driver(0.01*um,eastlocalStepper,eastlocalStepper->GetNumberOfVariables());
  G4ChordFinder* eastlocalChordFinder = new G4ChordFinder(eastlocalIntgrDriver);
  eastLocalFieldManager -> SetChordFinder(eastlocalChordFinder);

  eastLocalFieldManager -> GetChordFinder() -> SetDeltaChord(10*um);
  eastLocalFieldManager -> SetMinimumEpsilonStep(1e-6);
  eastLocalFieldManager -> SetMaximumEpsilonStep(1e-5);
  eastLocalFieldManager -> SetDeltaOneStep(0.1*um);

//  mwpc_container_log[0] -> SetFieldManager(eastLocalFieldManager, true);
  return;
}

void DetectorConstruction::ConstructWestMWPCField(G4double a, G4double b, G4double c, G4double d, G4RotationMatrix* e, G4ThreeVector f)
{
  G4cout << "Setting up West wirechamber electromagnetic field." << G4endl;
  MWPCField* westLocalField = new MWPCField();
  westLocalField -> SetActiveReg_d(a);
  westLocalField -> SetActiveReg_L(b);
  westLocalField -> SetActiveReg_r(c);
  westLocalField -> SetSideRot(e);
  westLocalField -> SetSideTrans(f);
  westLocalField -> SetPotential(d);

  G4FieldManager* westLocalFieldManager = new G4FieldManager();
  westLocalFieldManager -> SetDetectorField(westLocalField);

  G4EqMagElectricField* westlocalEquation = new G4EqMagElectricField(westLocalField);
  G4ClassicalRK4* westlocalStepper = new G4ClassicalRK4(westlocalEquation,8);
  G4MagInt_Driver* westlocalIntgrDriver = new G4MagInt_Driver(0.01*um,westlocalStepper,westlocalStepper->GetNumberOfVariables());
  G4ChordFinder* westlocalChordFinder = new G4ChordFinder(westlocalIntgrDriver);
  westLocalFieldManager -> SetChordFinder(westlocalChordFinder);

  westLocalFieldManager -> GetChordFinder() -> SetDeltaChord(10*um);
  westLocalFieldManager -> SetMinimumEpsilonStep(1e-6);
  westLocalFieldManager -> SetMaximumEpsilonStep(1e-5);
  westLocalFieldManager -> SetDeltaOneStep(0.1*um);

//  mwpc_container_log[1] -> SetFieldManager(westLocalFieldManager, true);
  return;
}
