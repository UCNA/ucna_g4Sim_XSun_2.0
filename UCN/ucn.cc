#include "DetectorConstruction.hh"
#include "PhysList495.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include <G4UnitsTable.hh>

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

int main(int argc,char** argv)
{
  // Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  G4int seed = time(NULL);
  G4Random::setTheEngine(new CLHEP::RanecuEngine);	// Choose the Random engine
  G4Random::setTheSeed(seed);

#ifdef G4MULTITHREADED	// Construct the default run manager
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysList495());

  // Chi Feng's optical photons physics list
//  G4String physicsList = "QGSP_BERT";
//  UCNbPhysicsList* physics = new UCNbPhysicsList(physicsList);
//  runManager->SetUserInitialization(physics);

  runManager->SetUserAction(new PrimaryGeneratorAction(detector));
  runManager->SetUserAction(new RunAction);
  EventAction* eventAction = new EventAction;
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(new SteppingAction(eventAction));

  new G4UnitDefinition("torr", "torr", "Pressure", atmosphere/760.);

  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( ! ui ) {		// batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {		// interactive mode
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
}
