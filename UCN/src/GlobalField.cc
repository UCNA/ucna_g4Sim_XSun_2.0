#include <cmath>

#include "GlobalField.hh"

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

GlobalField::GlobalField()
 : fSqOfMaxRadius((20*cm)*(20*cm)), fFieldScale(1.0)
{
  LoadFieldMap();
}

void GlobalField::LoadFieldMap()
{
  Bpoints.clear();
  Zpoints.clear();

  AddPoint(-3.0*m, 0.6*tesla);
  AddPoint(-2.2*m, 0.6*tesla);
  AddPoint(-1.5*m, 1.0*tesla);
  AddPoint(1.5*m, 1.0*tesla);
  AddPoint(2.2*m, 0.6*tesla);
  AddPoint(3.0*m, 0.6*tesla);
}

void GlobalField::AddPoint(G4double zPositions, G4double BValues)
{
  G4cout << "Setting field values: z = " << zPositions/m << " m, B = " << BValues/tesla << " tesla" << G4endl;
  Zpoints.push_back(zPositions);
  Bpoints.push_back(BValues);
}

void GlobalField::GetFieldValue(const G4double Point[3], G4double *Bfield) const
{
  G4double z = Point[2]; // point z
  unsigned int zindex = int(lower_bound(Zpoints.begin(), Zpoints.end(), z)-Zpoints.begin()); // location in points list

  if((zindex==0) || (zindex>=Zpoints.size()) || (Point[0]*Point[0]+Point[1]*Point[1]>fSqOfMaxRadius) || (!fFieldScale))
  {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;	// no field defined outside exp volume
  }
  else
  {
    // interpolate between defined regions
    G4double base = 0.5*(Bpoints[zindex-1]+Bpoints[zindex]);// midpoint value
    G4double amp = 0.5*(Bpoints[zindex-1]-Bpoints[zindex]); // variation amplitude between ends
    G4double dz = Zpoints[zindex]-Zpoints[zindex-1]; // z distance between ends
    G4double l = (z-Zpoints[zindex-1])/dz; // fractional distance between ends

    Bfield[2] = base*fFieldScale;

    if(amp)
    {
      Bfield[2] += amp*cos(l*M_PI)*fFieldScale; // interpolate B_z component with cosine
      // B_r component to obey Maxwell equation grad dot B = dB_z/dz + 1/r d(r B_r)/dr = 0
      G4double Brtemp = amp*M_PI*sin(l*M_PI)/(2*dz)*fFieldScale;
      Bfield[0] = Point[0]*Brtemp;
      Bfield[1] = Point[1]*Brtemp;
    }
    else
    {
      Bfield[0] = 0.0;
      Bfield[1] = 0.0;
    }
  }
}




