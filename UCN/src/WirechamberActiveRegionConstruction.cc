#include "WirechamberActiveRegionConstruction.hh"

#include <iostream>

WirechamberActiveRegionConstruction::WirechamberActiveRegionConstruction()
: dAnodeRadius(5*um), dCathodeRadius(25*um), dPlatingThick(0.2*um), dWireSpacing(2*mm),
iNbOfWires(64), dPlaneSpacing(1*cm)
{
  // here is space to set initial values of class members as well.
  // But I think we've already set them in the pre-constructor call.
}

void WirechamberActiveRegionConstruction::Build(int side)
{
  G4Material* cathodeWireMaterial = Al;
  G4Material* anodeWireMaterial = Wu;
  G4Material* cathodePlateMaterial = Au;

  if(!mMWPCGas)
    std::cout << "ERROR: active region gas of mwpc is a NULL pointer right now." << std::endl;

  G4double wirePlaneWidth = iNbOfWires*dWireSpacing;

  // Create shapes and visualizations to be set later.
  // effective mwpc gas volume containing the cathodes and anodes
  G4Box* mwpcGasBox = new G4Box("mpwc_gas_box", wirePlaneWidth/2., wirePlaneWidth/2., dPlaneSpacing);

  // anode, cathode wire containers. These are "on their side" to allow wireplane parametrization. Rotated later.
  G4Box* cathContainerBox = new G4Box("cathContainer_Box", wirePlaneWidth/2., dCathodeRadius, wirePlaneWidth/2.);
  G4Box* anodeContainerBox = new G4Box("anodeContainer_Box", wirePlaneWidth/2., dAnodeRadius, wirePlaneWidth/2.);

  // anode, cathode wires and surrouding gas
  G4Tubs* cathPlateTube = new G4Tubs("cathplate_tube", dCathodeRadius - dPlatingThick, dCathodeRadius, wirePlaneWidth/2., 0., 2*M_PI);
  G4Tubs* cathodeTube = new G4Tubs("cathode_tube", 0, dCathodeRadius - dPlatingThick, wirePlaneWidth/2., 0., 2*M_PI);
  G4Tubs* anodeTube = new G4Tubs("anode_tube", 0, dAnodeRadius, wirePlaneWidth/2., 0, 2*M_PI);
  G4Box* cathSegBox = new G4Box("cathodeSegmentBox", dWireSpacing/2., dCathodeRadius, wirePlaneWidth/2.);
  G4Box* anodeSegBox = new G4Box("anodeSegmentBox", dWireSpacing/2., dAnodeRadius, wirePlaneWidth/2.);

  // create some rotation matrices and visualizations for use later
  G4RotationMatrix* xRot90 = new G4RotationMatrix;      // rotate 90 degrees around X axis
  xRot90->rotateX(M_PI/2.*rad);
  G4RotationMatrix* xzRot90 = new G4RotationMatrix;     // rotate 90 deg around X then Z (Y axis in local coordinates)
  xzRot90->rotateX(M_PI/2.*rad);
  xzRot90->rotateY(M_PI/2.*rad);

  G4VisAttributes* visCathWires = new G4VisAttributes(G4Colour(1,0.7,0,0.8));
  G4VisAttributes* visAnodeWires = new G4VisAttributes(G4Colour(1,0.3,0,0.8));

  // anode, cathode segments
  gas_log = new G4LogicalVolume(mwpcGasBox, mMWPCGas, Append(side, "mwpc_gas_log_"));
  gas_log -> SetVisAttributes(G4VisAttributes::Invisible);
  cathSeg_log = new G4LogicalVolume(cathSegBox, mMWPCGas, Append(side, "cathSeg_log_"));
  anodeSeg_log = new G4LogicalVolume(anodeSegBox, mMWPCGas, Append(side, "anodeSeg_log_"));
  cathWire_log = new G4LogicalVolume(cathodeTube, cathodeWireMaterial, Append(side, "cathodeWire_log_"));
  cathPlate_log = new G4LogicalVolume(cathPlateTube, cathodePlateMaterial, Append(side, "cathode_plate_log_"));
  anodeWire_log = new G4LogicalVolume(anodeTube, anodeWireMaterial, Append(side, "anodeWire_log_"));
  cathSeg_log -> SetVisAttributes(G4VisAttributes::Invisible);
  anodeSeg_log -> SetVisAttributes(G4VisAttributes::Invisible);
  cathWire_log -> SetVisAttributes(visCathWires);
  cathPlate_log -> SetVisAttributes(visCathWires);
  anodeWire_log -> SetVisAttributes(visAnodeWires);

  // anode, cathode plane container volumes
  G4LogicalVolume* cathContainer1_log = new G4LogicalVolume(cathContainerBox, mMWPCGas, Append(side, "cathContainer1_log_"));
  G4LogicalVolume* cathContainer2_log = new G4LogicalVolume(cathContainerBox, mMWPCGas, Append(side, "cathContainer1_log_"));
  G4LogicalVolume* anodeContainer_log = new G4LogicalVolume(anodeContainerBox, mMWPCGas, Append(side, "anodeContainer_log_"));

  // place all segments in gas_log, which will then be "placed" in DetectorConstruction
  new G4PVPlacement(NULL, G4ThreeVector(), cathWire_log, Append(side, "cathode_wire_phys_"), cathSeg_log, true, 0);
  new G4PVPlacement(NULL, G4ThreeVector(), cathPlate_log, Append(side, "cathode_plate_phys_"), cathSeg_log, true, 0);
  new G4PVPlacement(NULL, G4ThreeVector(), anodeWire_log, Append(side, "anode_wire_phys_"), anodeSeg_log, true, 0);

  new G4PVPlacement(xRot90, G4ThreeVector(0., 0., dCathodeRadius - dPlaneSpacing),
      cathContainer1_log, Append(side, "cathContainer1_phys_"), gas_log, false, 0);
  new G4PVPlacement(xzRot90, G4ThreeVector(0., 0., (-1)*(dCathodeRadius - dPlaneSpacing)),
      cathContainer2_log, Append(side, "cathContainer2_phys_"), gas_log, false, 0);
  new G4PVPlacement(xRot90, G4ThreeVector(0,0,0), anodeContainer_log, Append(side, "anodeContainer_phys_"), gas_log, false, 0);

  // replicate the segments defined above into cathode, anode arrays
  new G4PVReplica(Append(side, "CathodeArray1_"), cathSeg_log, cathContainer1_log, kXAxis, iNbOfWires, dWireSpacing);
  new G4PVReplica(Append(side, "CathodeArray2_"), cathSeg_log, cathContainer2_log, kXAxis, iNbOfWires, dWireSpacing);
  new G4PVReplica(Append(side, "AnodeArray_"), anodeSeg_log, anodeContainer_log, kXAxis, iNbOfWires, dWireSpacing);

}
