#include "SourceHolderConstruction.hh"

#include <G4SubtractionSolid.hh>


SourceHolderConstruction::SourceHolderConstruction()
: dSourceWindowThickness(4.7*um), dSourceCoatingThickness(0.1*um),
mWindowMat(Mylar), mCoatingMat(Al),
dSourceHolderThickness((3./16.)*inch)
{
  // sets a messenger class that can communicate with detector construction?
/*  pUIdir = new G4UIdirectory("/sourceholder/");
  pWindowThickCmd = new G4UIcmdWithADoubleAndUnit("/sourceholder/windowthickness",this);
  pWindowThickCmd->AvailableForStates(G4State_PreInit);
  pWindowThickCmd->SetGuidance("thickness of windows on either side of sealed source");
  pWindowThickCmd->SetDefaultValue(4.7*um); */
}

/*void SourceHolderConstruction::SetNewValue(G4UIcommand * command, G4String newValue)
{
  if(command == pWindowThickCmd)
  {
    dSourceWindowThickness = pWindowThickCmd -> GetNewDoubleValue(newValue);
    G4cout << "Setting source holder window thickness to " << dSourceWindowThickness/um << "um" << G4endl;
  }
  else
  {
    G4out << "Unknown command:" << command->GetCommandName() << " passed to SourceHolderConstruction::SetNewValue\n";
  }
}
*/

void SourceHolderConstruction::Build()
{
//  if(!dSourceWindowThickness) dSourceWindowThickness = 0.001*um;

  G4double ringRadius = 0.5*inch;        // not pre-initialization variables
  G4double windowRadius = ringRadius-3.0*mm;
  G4double ringThickness = 3.2*mm;
  G4double holderHeight = 1.5*inch;
  G4double holderWidth = 1.5*inch;

  // source holder container
  G4Box* holderBox = new G4Box("source_holder_box", 0.5*holderWidth, 0.5*holderHeight, 0.5*dSourceHolderThickness);
  sourceContainer_log = new G4LogicalVolume(holderBox, Vacuum, "source_container_log");

  // source holder paddle
  G4Tubs* holderHole = new G4Tubs("source_holder_hole", 0., ringRadius, dSourceHolderThickness, 0., 2*M_PI);
  G4SubtractionSolid* holder = new G4SubtractionSolid("source holder", holderBox, holderHole);
  G4LogicalVolume* source_holder_log = new G4LogicalVolume(holder, Brass, "source_holder_log");
  source_holder_log -> SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.7,0,0.5)));
  sourceHolder_phys = new G4PVPlacement(NULL, G4ThreeVector(), source_holder_log, "source_holder_phys", sourceContainer_log, false, 0);

  // sealed source foil
  G4Tubs* windowTube = new G4Tubs("window_tube", 0., windowRadius, dSourceWindowThickness, 0., 2*M_PI);
  sourceWindow_log = new G4LogicalVolume(windowTube, mWindowMat, "source_window_log");
  G4VisAttributes* visWindow = new G4VisAttributes(G4Colour(0,1.0,0,1));
  sourceWindow_log->SetVisAttributes(visWindow);
  sourceWindow_phys = new G4PVPlacement(NULL, G4ThreeVector(), sourceWindow_log, "source_window_phys", sourceContainer_log, false, 0);

  // source foil coating
  G4Tubs* coatingTube = new G4Tubs("source_coating_tube", 0., windowRadius, 0.5*dSourceCoatingThickness, 0., 2*M_PI);
  for(int i = 0; i <= 1; i++)   // 0 = EAST, 1 = WEST
  {
    sourceCoating_log[i] = new G4LogicalVolume(coatingTube, mCoatingMat, Append(i, "source_coating_log_"));
    sourceCoating_log[i] -> SetVisAttributes(new G4VisAttributes(G4Colour(0,1,0,0.5)));
  }

  sourceCoating_phys[0] = new G4PVPlacement(NULL, G4ThreeVector(0,0, (-1)*(dSourceWindowThickness + dSourceCoatingThickness*0.5)),
                                        sourceCoating_log[0], "source_coating_phys_EAST", sourceContainer_log, false, 0);
  sourceCoating_phys[1] = new G4PVPlacement(NULL, G4ThreeVector(0,0, dSourceWindowThickness + dSourceCoatingThickness*0.5),
                                        sourceCoating_log[1], "source_coating_phys_WEST", sourceContainer_log, false, 0);

  // source retaining ring
  G4Tubs* ringTube = new G4Tubs("source_ring_tube", windowRadius, ringRadius, ringThickness/2., 0., 2*M_PI);
  G4LogicalVolume* ring_log = new G4LogicalVolume(ringTube, Al, "source_ring_log");
  ring_log -> SetVisAttributes(new G4VisAttributes(G4Colour(0.7,0.7,0.7,0.5)));
  sourceRing_phys = new G4PVPlacement(NULL, G4ThreeVector(), ring_log, "source_ring_phys", sourceContainer_log, false, 0);



}
