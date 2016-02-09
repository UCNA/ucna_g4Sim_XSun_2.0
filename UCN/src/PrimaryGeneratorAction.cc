#include "PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <string>
using   namespace       std;

#define	OUTPUT_FILE	"UCNASimOutput.txt"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fMyDetector(myDC),
//  fSourceRadius(3.*mm),	// this is default set. If you aren't using DiskRandom, don't care.
  fSourceRadius(0),
  fPosOffset(G4ThreeVector(0,0,0))	// if it's not a source, then set to 0
{
  G4int nPtcls = 1;
  fParticleGun = new G4ParticleGun(nPtcls);

  G4String base = getenv("UCNA_BASE");
  G4String path = base + "/UCN/EventGenTools/G4Sim_Ptcl_Input_Files/";
  G4String file = path + "big_initPtclInfo.txt";

  G4cout << "------> Path to primaries file: " << path << G4endl;
  G4cout << "Fetching initial particles info from file name: " << file << G4endl;

  LoadFile(file);
  // At the end of constructor, GEANT4 default calls GeneratePrimaries method
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // use GEANT4 event id to track which part of fEvtsArray we will use for generated event
  int nID = anEvent -> GetEventID();

  // use the particle species flag to set particle type
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("geantino");
  if(fEvtsArray[nID].event_speciesFlag == 11)
  {
    particle = particleTable->FindParticle("e-");
  }
  else if(fEvtsArray[nID].event_speciesFlag == 22)
  {
    particle = particleTable->FindParticle("gamma");
  }
  else
  {
    G4cout << "No matching particle species definition." << G4endl;
  }
  fParticleGun -> SetParticleDefinition(particle);

  fParticleGun -> SetParticleEnergy(fEvtsArray[nID].event_energy*keV);

  fParticleGun -> SetParticlePosition(G4ThreeVector(fEvtsArray[nID].event_xPos*m,
						fEvtsArray[nID].event_yPos*m,
						fEvtsArray[nID].event_zPos*m));

  fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(fEvtsArray[nID].event_xMo,
						fEvtsArray[nID].event_yMo,
						fEvtsArray[nID].event_zMo));


  fParticleGun -> SetParticleTime(fEvtsArray[nID].event_time*ns);

  // Call to method to save primary particle initial info.
  // Need to be super careful here. A Priori, there's no reason that the momentum vector
  // is normalized when I'm reading it from final.
  // And momentum direction is a normalized vector once GEANT4 gets a hold of it.
  // So when we print out what was read in vs. what we get from particle gun may not be the same.
  // Will need to check this. For now it is the same within float to double rounding.
  SavePrimPtclInfo(nID);

  fParticleGun -> GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction::LoadFile(G4String fileName)
{
  event eRead;
  int i = 0;

  string buf;
  ifstream infile;
  infile.open(fileName);

  //a check to make sure the file is open
  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  while(getline(infile, buf))
  {
    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> eRead.event_gen_id
		>> eRead.event_speciesFlag
		>> eRead.event_energy
		>> eRead.event_xPos >> eRead.event_yPos >> eRead.event_zPos
		>> eRead.event_xMo >> eRead.event_yMo >> eRead.event_zMo
		>> eRead.event_time
		>> eRead.event_weight;

      fEvtsArray[i] = eRead;
      i++;
    }
  }

}

void PrimaryGeneratorAction::DiskRandom(G4double radius, G4double& x, G4double& y)
{
  while(true)
  {
    x = (2.0*G4UniformRand()-1.)*radius;
    y = (2.0*G4UniformRand()-1.)*radius;
    if(x*x+y*y<=radius*radius) break;
  }
}

void PrimaryGeneratorAction::DisplayGunStatus()
{
  G4cout
  << fParticleGun->GetParticleDefinition()->GetParticleName()
  << " gun from " << fParticleGun->GetParticlePosition()/m
  << "m towards " << fParticleGun->GetParticleMomentumDirection()
  << " at " << fParticleGun->GetParticleTime()/ns
  << "ns : gun firing at energy " << fParticleGun->GetParticleEnergy()/keV
  << "keV" <<
  G4endl;
}

void PrimaryGeneratorAction::SavePrimPtclInfo(int index)
{
  ofstream outfile;     // output initial particle information into same file as final sim output
  outfile.open(OUTPUT_FILE, ios::app);
  outfile << fEvtsArray[index].event_gen_id << "\t"
	<< fEvtsArray[index].event_speciesFlag << "\t"
	<< fEvtsArray[index].event_energy << "\t"
	<< fEvtsArray[index].event_xPos << "\t"
	<< fEvtsArray[index].event_yPos << "\t"
	<< fEvtsArray[index].event_zPos << "\t"
	<< fEvtsArray[index].event_xMo << "\t"	// be careful here with momentum
	<< fEvtsArray[index].event_yMo << "\t"	// see comment back in GeneratePrimaries
	<< fEvtsArray[index].event_zMo << "\t"
	<< fEvtsArray[index].event_time << "\t"
	<< fEvtsArray[index].event_weight << "\t";	// has to be a \t since getting appended in EventAction

	// while debugging TrackerSD and TrackerHit, this is what we print.
//  outfile << fEvtsArray[index].event_gen_id << "\t"
//	<< fEvtsArray[index].event_energy << "\t";
  outfile.close();
}






























//***** OLD CODE I WROTE EARLIER TO "TEST" 113SN SOURCE *****//

void PrimaryGeneratorAction::Set_113SnSource()	// don't need additional arguments since we set the particle gun.
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  //----- Setting species and energy of particle decided by CE random sampling from nndc
  int r1;
  r1 = rand() % 1007246 + 1;    // This bound is # of digits I want to produce
                                // NOTE not set to 100 because on nndc we get 100.7% for 391 keV gammas.
  double percentage = r1/10000.;        // This gives us 0.001 precision.

  if((percentage >= 0) && (percentage <= 64.97))
  {
    fParticleGun -> SetParticleEnergy(391.698*keV);
    particle = particleTable->FindParticle(particleName="gamma");
    fParticleGun -> SetParticleMomentumDirection(G4ThreeVector(0,1,0));
  }
  else if((percentage > 64.97) && (percentage <= 93.77))
  {
    fParticleGun -> SetParticleEnergy(363.758*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 93.77) && (percentage <= 99.37))
  {
    fParticleGun -> SetParticleEnergy(387.461*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 99.37) && (percentage <= 100.507))
  {
    fParticleGun -> SetParticleEnergy(390.872*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 100.507) && (percentage <= 100.712))
  {
    fParticleGun -> SetParticleEnergy(391.576*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 100.712) && (percentage <= 100.7246))
  {
    fParticleGun -> SetParticleEnergy(391.697*keV);
    particle = particleTable -> FindParticle(particleName="e-");
  }
  else if((percentage > 100.7246) || (percentage < 0))
  {
    G4cout << "Random number sampled beyond the scope of the decay." << G4endl;
  }

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);        // Michael's has this line. Idk why.

  //----- Setting isotropic particle momentum direction
  G4double alphaMin = 0*deg;    // alpha is apex angle
  G4double alphaMax = 180*deg;  // 180* ensures cone -> full sphere
  G4double cosAlphaMin = cos(alphaMin);
  G4double cosAlphaMax = cos(alphaMax);
  G4double phiMin = 0*deg;      // phi in 0 to 2pi
  G4double phiMax = 360.*deg;   // Presumably the rotation angle. 2pi makes it a cone.

  G4double cosAlpha = cosAlphaMin-G4UniformRand()*(cosAlphaMin-cosAlphaMax);
  G4double sinAlpha = sqrt(1. - cosAlpha*cosAlpha);
  G4double phi = phiMin + G4UniformRand()*(phiMax - phiMin);

  G4double ux = sinAlpha*cos(phi);
  G4double uy = sinAlpha*sin(phi);
  G4double uz = cosAlpha;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));

/*  G4ThreeVector newUz;        // fires isotropic cone where cone axis can be rotated.
  G4double theta, phi, apex;

  G4double xCentre = 0*cm;
  G4double yCentre = 0*cm;
  G4double zCentre = -1*cm;

  phi = (180 + atan(yCentre/xCentre)*(180.0/3.1416))*deg; // basic geometry and then converted to degrees.
  theta = acos(zCentre/sqrt(xCentre*xCentre + yCentre*yCentre + zCentre*zCentre))*(180.0/3.1416)*deg;

  apex = 180*deg;

  newUz = G4ThreeVector(std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta));

  G4double cosAlpha = 1. - G4UniformRand()*(1.- std::cos(apex));
  G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
  G4double psi      = 2*3.1416*G4UniformRand();  //psi uniform in [0, 2*pi]
  G4ThreeVector dir(sinAlpha*std::cos(psi), sinAlpha*std::sin(psi), cosAlpha);
  dir.rotateUz(newUz);

  fParticleGun->SetParticleMomentumDirection(dir);
*/

  //----- Setting the particle generation position
  G4double x0 = 0;              // Says it is negligibly thin.
  G4double y0 = 0;              // Brad told me the source radius
  G4double z0 = 0;
  if(fSourceRadius != 0)
  {
    DiskRandom(fSourceRadius, x0, y0);
  }

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

}
