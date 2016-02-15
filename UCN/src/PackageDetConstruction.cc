#include "PackageDetConstruction.hh"

PackageDetConstruction::PackageDetConstruction()
: dDetPackageRadius(6.0*inch), dMWPCEntranceThick(0.375*inch), dMWPCEntranceRad(3.0*inch),
dMWPCEntranceDepth(5.0*inch), dFrontWinFrameThick(1.0*inch), dBackWinFrameThick(0.5*inch)
{
  // again can use this initializer but class members are already set
}

void PackageDetConstruction::Build(int side)
{
  // create component objects that will be placed inside the detector package
  Scint.Build(side);
  Wirechamber.Build(side);

  // make aluminum entrance collimator to detector package
  G4double entranceSectionLength = dMWPCEntranceDepth + dFrontWinFrameThick;
  G4Tubs* mwpcEntranceTube = new G4Tubs("mwpc_entrance_tube", 0., dDetPackageRadius, 0.5*entranceSectionLength, 0., 2*M_PI);
  G4Tubs* entranceFrontTube = new G4Tubs("entrance_front_tube", dMWPCEntranceRad + dMWPCEntranceThick,
                                        dDetPackageRadius, 0.5*dMWPCEntranceThick, 0., 2*M_PI);
  G4Tubs* entranceMidTube = new G4Tubs("entrance_mid_tube", dMWPCEntranceRad, dMWPCEntranceRad + dMWPCEntranceThick,
                                        0.5*dMWPCEntranceDepth, 0., 2*M_PI);
  G4Tubs* entranceBackTube = new G4Tubs("entrance_back_tube", dMWPCEntranceRad, dDetPackageRadius,
                                        0.5*dFrontWinFrameThick, 0., 2*M_PI);
  G4VisAttributes* visMWPCEntrance = new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.8));

  mwpcEntrance_log = new G4LogicalVolume(mwpcEntranceTube, Vacuum, Append(side, "mwpc_entrance_log_"));
  mwpcEntrance_log->SetVisAttributes(G4VisAttributes::Invisible);
  entranceFront_log = new G4LogicalVolume(entranceFrontTube, Al, Append(side, "entrance_front_log_"));
  entranceMid_log = new G4LogicalVolume(entranceMidTube, Al, Append(side, "entrance_mid_log_"));
  entranceBack_log = new G4LogicalVolume(entranceBackTube, Al, Append(side, "entrance_back_log_"));
  entranceFront_log->SetVisAttributes(visMWPCEntrance);
  entranceMid_log->SetVisAttributes(visMWPCEntrance);
  entranceBack_log->SetVisAttributes(visMWPCEntrance);

  entranceFront_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*(entranceSectionLength - dMWPCEntranceThick)),
                                  entranceFront_log, Append(side, "entrance_front_phys_"), mwpcEntrance_log, false, 0);
  entranceMid_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -0.5*dFrontWinFrameThick),
                                  entranceMid_log, Append(side, "entrance_mid_phys_"), mwpcEntrance_log, false, 0);
  entranceBack_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, 0.5*(entranceSectionLength - dFrontWinFrameThick)),
                                  entranceBack_log, Append(side, "entrance_back_phys_"), mwpcEntrance_log, false, 0);

  // overall detector package, aka container log which will get placed in DetectorConstruction
  G4double detFrameHalf_Z = dMWPCEntranceDepth + Wirechamber.GetWidth() + 1.0*inch;
  G4Tubs* packageTube = new G4Tubs("detPackage_tube", 0, dDetPackageRadius, detFrameHalf_Z, 0, 2*M_PI);
  container_log = new G4LogicalVolume(packageTube, Vacuum, Append(side, "container_log_"));
  container_log -> SetVisAttributes(G4VisAttributes::Invisible);

  // place all the components of the detector package relative to scintillator face position, located at (local) z=0
  scint_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,-Scint.GetScintFacePos()),
					Scint.container_log, Append(side, "scintContainer_"), container_log, false, 0, true);

  G4double mwpcPosZ = -Wirechamber.GetWidth()/2. - dBackWinFrameThick-(Scint.GetWidth()/2. + Scint.GetScintFacePos());

  mwpc_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0.,mwpcPosZ),
					  Wirechamber.container_log, Append(side, "mwpcContainer_"), container_log, false, 0, true);

  Wirechamber.vMyTranslation[2] += mwpcPosZ;	// add the shift in Z position to the wirechamber's translation matrix

  G4double entrancePosZ = mwpcPosZ - (Wirechamber.GetWidth() + entranceSectionLength)/2.;
  dEntranceFacePos_Z = entrancePosZ - 0.5*entranceSectionLength;
  dEntranceWinPos_Z = entrancePosZ + 0.5*entranceSectionLength;
  mwpcEntrance_phys = new G4PVPlacement(NULL, G4ThreeVector(0.,0., entrancePosZ),
			mwpcEntrance_log, Append(side, "mwpc_entrance_"), container_log, false, 0);

  // construct aluminum exit window and N2 volume at back of gas box.
  G4Tubs* mwpcExitTube = new G4Tubs("mwpc_exit_tube", Wirechamber.dMWPCExitR, dDetPackageRadius, 0.5*dBackWinFrameThick,0.,2*M_PI);
  G4VisAttributes* visMWPCExit = new G4VisAttributes(G4Colour(0.3,0.3,0.3,0.8));
  mwpcExit_log = new G4LogicalVolume(mwpcExitTube, Al, Append(side, "mwpc_exit_log_"));
  mwpcExit_log->SetVisAttributes(visMWPCExit);
  dExitFramePos_Z = mwpcPosZ + (Wirechamber.GetWidth() + dBackWinFrameThick)/2.;
  mwpcExit_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0., dExitFramePos_Z),
			mwpcExit_log, Append(side, "mwpc_exit_"), container_log, false, 0);

  G4Tubs* mwpcExitN2Tube = new G4Tubs("mwpc_exit_N2_tube", 0, Wirechamber.dMWPCExitR, 0.5*dBackWinFrameThick, 0., 2*M_PI);
  mwpcExitGasN2_log = new G4LogicalVolume(mwpcExitN2Tube, WCNitrogen, Append(side, "mwpc_exit_N2_log_"));
  mwpcExitGasN2_log -> SetVisAttributes(G4VisAttributes::Invisible);
  mwpcExitGasN2_phys = new G4PVPlacement(NULL,G4ThreeVector(0.,0., dExitFramePos_Z),
			mwpcExitGasN2_log, Append(side, "mwpc_exit_"), container_log, false, 0);

  // material behind the detector. Misc stuff that can cause back scattering events.
  G4double backStuffThick = 1.0*inch;
  G4Tubs* backStuffTube = new G4Tubs("backstuff_tube", 0, 0.5*dDetPackageRadius, backStuffThick, 0., 2*M_PI);
  backStuff_log = new G4LogicalVolume(backStuffTube, SS304, Append(side, "backStuff_log_"));
  backStuff_phys = new G4PVPlacement(NULL, G4ThreeVector(0, 0, detFrameHalf_Z - 0.5*backStuffThick),
                                  backStuff_log, Append(side, "backStuff_phys_"), container_log, false, 0);
}
