#include <cmath>

#include "MWPCField.hh"

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

MWPCField::MWPCField()
 : fE0(0), fSqOfMaxRadius((20*cm)*(20*cm)), fFieldScale(1.0)
{
  G4cout << "Creating MWPC electromagnetic field objects." << G4endl;
  fChamberRot = NULL;	// initialize some class members
  fChamberTrans = G4ThreeVector(0,0,0);

  LoadFieldMap();
}

void MWPCField::LoadFieldMap()
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

void MWPCField::AddPoint(G4double zPositions, G4double BValues)
{
//  G4cout << "Setting field values: z = " << zPositions/m << " m, B = " << BValues/tesla << " tesla" << G4endl;
  Zpoints.push_back(zPositions);
  Bpoints.push_back(BValues);
}

void MWPCField::SetPotential(G4double Vanode)
{
  fE0 = M_PI*Vanode/d/log(sinh(M_PI*L/d)/sinh(M_PI*r/d));
  G4cout << "Wirechamber voltage set to " << Vanode/volt <<" V => fE0 = " << fE0/(volt/cm) << " V/cm" << G4endl;
}

void MWPCField::GetFieldValue(const G4double Point[4], G4double *Bfield) const
{
  // magnetic field components computed below. Same code as in the GlobalField GetFieldValue(...)
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

  if(!fE0)	// if no electric potential, then set the E-field parts to 0 and leave
  {
    Bfield[3] = 0;
    Bfield[4] = 0;
    Bfield[5] = 0;
    return;
  }

  // create a local position vector to avoid offsets in computation
  G4ThreeVector localPos = G4ThreeVector(Point[0], Point[1], Point[2]) - fChamberTrans;
  if(fChamberRot != NULL)
  {
    localPos = (*fChamberRot)(localPos);
  }

  // compute the electric field components
  G4ThreeVector E(0,0,0);
  double l = localPos[2];
  if(fabs(l)<L)
  {
    double a = localPos[0]/d;
    a = (a-floor(a)-0.5)*d;
    if(a*a+l*l > r*r)
    {
      double denom = cosh(2*M_PI*l/d)-cos(2*M_PI*a/d);
      E[2] = fE0*sinh(2*M_PI*l/d)/denom;
      E[0] = fE0*sin(2*M_PI*a/d)/denom;
    }
  }

  if(fChamberRot != NULL)	// rotate back to global coordinate frame
  {
    E = fChamberRot->inverse()(E);
  }
  Bfield[3] = E[0];	// setting the 4,5,6th components of Bfield, which is a 6 entry array
  Bfield[4] = E[1];	// to the electric field values. This needs to be done for GEANT4 EM field.
  Bfield[5] = E[2];

}




