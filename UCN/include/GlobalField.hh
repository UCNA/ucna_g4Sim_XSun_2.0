#ifndef GlobalField_h
#define GlobalField_h 1

#include <vector>
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;

using namespace std;

class GlobalField: public G4MagneticField
{
public:
  GlobalField();

  void LoadFieldMap();
  void GetFieldValue( const G4double Point[3], G4double *Bfield ) const;
  void SetFieldScale(G4double val) { fFieldScale = val; }

private:
  void AddPoint(G4double zPositions, G4double BValues);
  vector<G4double> Bpoints; ///< field profile B values
  vector<G4double> Zpoints; ///< field profile z positions

  G4double fSqOfMaxRadius;
  double fFieldScale;				// dimensionless scaling factor

};

#endif
