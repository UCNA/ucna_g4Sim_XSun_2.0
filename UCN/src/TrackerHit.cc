#include "TrackerHit.hh"

#include <G4SystemOfUnits.hh>
#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

G4Allocator<TrackerHit> TrackerHitAllocator;

TrackerHit::TrackerHit()
{
  fEdep = 0;
  fEdepQuenched = 0;
  fHitPos = G4ThreeVector();		// G4ThreeVector is initially set to zero vector
  fEdepWeightedPos = G4ThreeVector();	// this doesn't match Michael Mendenhall's code
  fEdepWeightedPos2 = G4ThreeVector();	// but hopefully accomplishes same thing
  fIncidentMo = G4ThreeVector();
  fExitMo = G4ThreeVector();
  fTrackVertex = G4ThreeVector();


  sEdep = 0;	// keep them here for now
}
/*eDepSoFar(0), eDepQuenchedSoFar(0), hitPosition(), edepWeightedPosition(),
edepWeightedPosition2(), incidentMomentum(), exitMomentum(), vertex() {} */

void TrackerHit::Print()
{
  G4cout << "Printing tracker hit info. There's nothing here." << G4endl;
/*  G4cout	<< "  trackID: " << trackID
			<< "  vertex: " << G4BestUnit(vertex,"Length")
			<< "  created in " << creatorVolumeName
			<< "  in " << volumeName
			<< "  incident energy " << G4BestUnit(incidentEnergy,"Energy")
			<< "  position: " << G4BestUnit(hitPosition,"Length")
			<< "  time: " << G4BestUnit(hitTime,"Time")
			<< "  edep: " << G4BestUnit(eDepSoFar,"Energy")
			<< "  edep quenched: " << G4BestUnit(eDepQuenchedSoFar,"Energy")
			<< "  incident momentm: " << G4BestUnit(incidentMomentum,"Energy")
			<< "  exit momentum " << G4BestUnit(exitMomentum,"Energy")
			<< G4endl;
*/
}


/*void TrackerHit::fillTrackInfo(TrackInfo h) const {
	h.trackID = GetTrackID();
	h.hitTime = GetHitTime()/ns;
	h.KE = GetIncidentEnergy()/keV;
	h.Edep = GetEdep()/keV;
	h.EdepQuenched = GetEdepQuenched()/keV;
	h.isEntering = (originEnergy==0);

	for(int d=0; d<=2; d++) {
		h.edepPos[d] = GetEdepPos()[d]/(cm*keV);
		h.edepPos2[d] = GetEdepPos2()[d]/(cm*cm*keV);
		h.vertexPos[d] = GetVertex()[d]/cm;
		h.inPos[d] = GetPos()[d]/cm;
		h.pIn[d] = GetIncidentMomentum()[d]/keV;
		h.pOut[d] = GetExitMomentum()[d]/keV;
	}

	h.pID = GetPID();
}*/

