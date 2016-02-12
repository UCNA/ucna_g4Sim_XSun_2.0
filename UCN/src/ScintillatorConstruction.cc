#include "ScintillatorConstruction.hh"

#include <G4Polycone.hh>

ScintillatorConstruction::ScintillatorConstruction()
: dScintRadius(7.5*cm), dBackingRadius(10*cm), dScintThick(3.5*mm), dDeadThick(3.0*um),
dBackingThick(1.*inch), dLightGuideThick(1.0*cm)
{
  // can do constructor initializations here but all class members have already been set (defaulted)
}

void ScintillatorConstruction::Build(int side)
{
  if((dBackingRadius < dScintRadius) || (dLightGuideThick < dScintThick))
        G4cout << "\n\nMajor geometry error! Scintillator measurements don't make sense! \n \n" << G4endl;

  dN2Volume_Z = dLightGuideThick + dBackingThick;
  dScintFace_PosZ = -dN2Volume_Z/2.;

  //---- Create the shapes used in the scintillator object
  // Overall container layer for the scintillator
  G4Tubs* N2VolTube = new G4Tubs("N2_vol_tube", 0., dBackingRadius, dN2Volume_Z/2., 0., 2*M_PI);

  // dead layer in scint
  G4Tubs* deadLayerTube = new G4Tubs("Dead_scint_tube", 0, dScintRadius, dDeadThick/2., 0., 2*M_PI);
  G4VisAttributes* visDScint= new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.5));

  // scintillator
  G4Tubs* scintTube = new G4Tubs("scint_tube", 0, dScintRadius, (dScintThick - dDeadThick)/2., 0., 2*M_PI);
  G4VisAttributes* visScint= new G4VisAttributes(G4Colour(0.0,1.0,1.0,0.2));

  // light guides around and behind detector
  G4double zPlane[] = {0., dScintThick, dScintThick, dLightGuideThick};
  G4double lightGuideRadius = dScintRadius - (dLightGuideThick - dScintThick);
  G4double innerRad[] = {dScintRadius, dScintRadius, lightGuideRadius, lightGuideRadius};
  G4double outerRad[] = {dBackingRadius, dBackingRadius, dBackingRadius, dBackingRadius};

  G4Polycone* lightGuidePoly = new G4Polycone("lightguide_polycone", 0., 2*M_PI, 4, zPlane, innerRad, outerRad);
  G4VisAttributes* visLG = new G4VisAttributes(G4Colour(0.0,1.0,0.5,0.2));

  // backing veto
  G4Tubs* backingTube = new G4Tubs("backing_tube", 0., dBackingRadius, dBackingThick/2., 0., 2*M_PI);
  G4VisAttributes* visBacking= new G4VisAttributes(G4Colour(0.0,0.0,1,0.2));

  // create logical volumes for each object and assign visualizations
  container_log = new G4LogicalVolume(N2VolTube, WCNitrogen, Append(side, "N2_Vol_log_"));
  container_log -> SetVisAttributes(G4VisAttributes::Invisible);
  deadLayer_log = new G4LogicalVolume(deadLayerTube, Sci, Append(side, "Dead_scint_log_"));
  deadLayer_log -> SetVisAttributes(visDScint);
  scintillator_log = new G4LogicalVolume(scintTube, Sci, Append(side, "scint_log_"));
  scintillator_log -> SetVisAttributes(visScint);
  lightGuide_log = new G4LogicalVolume(lightGuidePoly, Sci, Append(side, "light_guide_log_"));
  lightGuide_log -> SetVisAttributes(visLG);
  backing_log = new G4LogicalVolume(backingTube, Sci, Append(side, "backing_log_"));
  backing_log -> SetVisAttributes(visBacking);

  // placement of the physical geometries in the scintillator container_log
  deadLayer_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -(dN2Volume_Z - dDeadThick)/2.),
                                deadLayer_log, Append(side, "Dead_scint_phys_"), container_log, false, 0);
  scintillator_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, -dN2Volume_Z/2. + dDeadThick + (dScintThick - dDeadThick)/2.),
                                scintillator_log, Append(side, "scint_crystal_phys_"), container_log, false, 0);
  lightGuide_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, dN2Volume_Z/2.), lightGuide_log,
                                Append(side, "light_guide_phys_"), container_log, false, 0);
  backing_phys = new G4PVPlacement(NULL, G4ThreeVector(0,0, (dN2Volume_Z - dBackingThick)/2.),
                                backing_log, Append(side, "backing_phys_"), container_log, false, 0);
}
