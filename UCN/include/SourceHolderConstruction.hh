#ifndef SourceHolderConstruction_HH
#define SourceHolderConstruction_HH

#include "DetectorConstructionUtils.hh"

#include <G4SubtractionSolid.hh>	// comment out when you're done moving all source holder stuff

//#include <G4UImessenger.hh>		// used for Detector Messenger class
//#include <G4UIdirectory.hh>
//#include <G4UIcmdWithADoubleAndUnit.hh>


class SourceHolderConstruction: public MaterialUser
{
public:
  SourceHolderConstruction();	// constructor

  G4double GetSourceHolderThickness() { return dSourceHolderThickness; };

  G4double dSourceWindowThickness;
  G4double dSourceCoatingThickness;
  G4Material* mWindowMat;
  G4Material* mCoatingMat;

  G4LogicalVolume* sourceContainer_log;
  G4LogicalVolume* sourceWindow_log;
  G4LogicalVolume* sourceCoating_log[2];

  void Build();		// construct source holder logical volume

//  virtual void SetNewValue(G4UIcommand * command, G4String newValue);

protected:
  G4VPhysicalVolume* sourceWindow_phys;
  G4VPhysicalVolume* sourceCoating_phys[2];
  G4VPhysicalVolume* sourceHolder_phys;
  G4VPhysicalVolume* sourceRing_phys;

  G4double dSourceHolderThickness;

private:
//  G4UIdirectory* pUIdir;			// UI Directory for source holder construction
//  G4UIcmdWithADoubleAndUnit* pWindowThickCmd;	// source holder window thickness command


};

#endif
