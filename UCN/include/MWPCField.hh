#ifndef MWPCField_h
#define MWPCField_h 1

#include <vector>
#include <G4ElectroMagneticField.hh>
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;

using namespace std;

class MWPCField: public G4ElectroMagneticField
{
public:
  MWPCField();

  void LoadFieldMap();

  // get values in EM field. Note this method needs to be implemented or EM field doesn't work!
  void GetFieldValue( const G4double Point[4], G4double *Bfield ) const;
  // whether the field changes particle energy, also needs to be implemented
  G4bool DoesFieldChangeEnergy() const { return fE0 != 0; };


  void SetPotential(G4double Vanode);
  void SetActiveReg_d(G4double activeRegion_d) {d = activeRegion_d;};
  void SetActiveReg_L(G4double activeRegion_L) {L = activeRegion_L;};
  void SetActiveReg_r(G4double activeRegion_r) {r = activeRegion_r;};
  void SetSideRot(G4RotationMatrix* sideRot) {fChamberRot = sideRot;};
  void SetSideTrans(G4ThreeVector sideTrans) {fChamberTrans = sideTrans;};

protected:
  G4double fE0;		// apparently a field scaling constant

private:
  void AddPoint(G4double zPositions, G4double BValues);
  vector<G4double> Bpoints; ///< field profile B values
  vector<G4double> Zpoints; ///< field profile z positions

  G4double fSqOfMaxRadius;			// These two needed for making Mag field component of global field
  double fFieldScale;				// dimensionless scaling factor

  double d;
  double L;
  double r;
  G4RotationMatrix* fChamberRot;
  G4ThreeVector fChamberTrans;

};

#endif
