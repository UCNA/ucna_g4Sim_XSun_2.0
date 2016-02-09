#ifndef TrackerHit_h
#define TrackerHit_h

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4Allocator.hh>
#include <G4ThreeVector.hh>
#include <G4String.hh>

/*struct TrackInfo
{
  int trackID;          // ID number for this track
  int hcID;             // SD ID number for this track
  double hitTime;       // time of track start [ns]
  double KE;            // KE at start of track [keV]
  double Edep;          // accumulated deposited energy along track (not countrying secondaries) [keV]
  double EdepQuenched;  // accumulated "quenched" energy along track [keV]
  double pIn[3];        // momentum at entry to volume
  double pOut[3];       // momentum at exit to volume

  double edepPos[3];    // hit position weighted by deposited energy [keV*cm]
  double edepPos2[3];   // hit position^2 weighted by deposited energy [keV*cm]

  double inPos[3];	// entry position to volume [cm]
  double vertexPos[3];	// track vertex position [cm]

  int pID;              // particle creating this track (PDG code)
  bool isEntering;      // whether this is initial track entering a volume

};
*/

// accumulates segment-by-segment information for a track in an SD
class TrackerHit : public G4VHit
{
  public:
    TrackerHit();
    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class. Won't bother really defining them. Not gonna use.
    virtual void Draw() { };
    virtual void Print();

    // Methods that add/accumulate private member variables
    void Add(G4double deltaE) {sEdep = sEdep + deltaE; };	// this was for the steppingAction accumulation (old)
    void AddEdep(G4double edep, G4ThreeVector xyz)
    {
      fEdep = fEdep + edep;
      fEdepWeightedPos = fEdepWeightedPos + xyz*edep;
      for(int i = 0; i < 3; i++)
      {
        fEdepWeightedPos2[i] = fEdepWeightedPos2[i] + xyz[i]*xyz[i]*edep;
      }
    };
    void AddEdepQuenched(G4double edep) {fEdepQuenched = fEdepQuenched + edep; };

    // Setter methods
    void SetTrackID(G4int trackid) { fTrackID = trackid; };
    void SetIncidentEnergy(G4double energy) { fIncidentE = energy; };
    void SetHitPos(G4ThreeVector xyz) { fHitPos = xyz; };
    void SetHitTime(G4double time) { fHitTime = time; };
    void SetIncidentMomentum(G4ThreeVector pIn) { fIncidentMo = pIn; };
    void SetExitMomentum(G4ThreeVector pOut) { fExitMo = pOut; };
    void SetPtclSpeciesID(G4int pID) { fPtclSpeciesID = pID; };
    void SetProcessName(G4String sProc) { fProcessName = sProc; };
    void SetVolumeName(G4String sVol) { fVolumeName = sVol; };
    void SetTrackVertex(G4ThreeVector xyz) { fTrackVertex = xyz; };
    void SetCreatorVolumeName(G4String sName) { fCreatorVolumeName = sName; };

    // Getter methods
    G4double GetStepEdep() { return sEdep; };	// won't need this when done

    G4int 		GetTrackID() { return fTrackID; };
    G4double 		GetIncidentEnergy() { return fIncidentE; };
    G4double 		GetEdep() { return fEdep; };
    G4double 		GetEdepQuenched() { return fEdepQuenched; };
    G4ThreeVector       GetHitPos() { return fHitPos; };
    G4double		GetHitTime() { return fHitTime; };
    G4ThreeVector	GetEdepPos() { return fEdepWeightedPos; };
    G4ThreeVector	GetEdepPos2() { return fEdepWeightedPos2; };
    G4ThreeVector	GetIncidentMomentum() { return fIncidentMo; };
    G4ThreeVector 	GetExitMomentum() { return fExitMo; };
    G4int		GetPtclSpeciesID() { return fPtclSpeciesID; };
    G4String		GetProcessName() { return fProcessName; };
    G4String		GetVolumeName() { return fVolumeName; };
    G4ThreeVector	GetTrackVertex() { return fTrackVertex; };
    G4String		GetCreatorVolumeName() { return fCreatorVolumeName; };

    G4double		fOriginEnergy;		// energy at split from "originating" track for EdepQ tracking
    unsigned int	fNbSecondaries;		// number of secondaries produced along track

  private:
    G4double sEdep;	// delete this once you complete the move over. Using 0th element of HC

    G4int		fTrackID;		// ID # for this track
    G4double		fIncidentE;		// incident kinetic energy at start of track
    G4double		fEdep;			// accumulator for energy deposition
    G4double		fEdepQuenched;		// accumulator for quenched energy
    G4double		fHitTime;		// time at entry into volume
    G4ThreeVector	fHitPos;		// position where this track entered volume
    G4ThreeVector	fEdepWeightedPos; 	// track position weighted by deposited energy
    G4ThreeVector	fEdepWeightedPos2;	// track pos^2 weighted by deposited energy
    G4ThreeVector	fIncidentMo;		// momentum when entering volume
    G4ThreeVector	fExitMo;		// momentum when exiting volume
    G4int		fPtclSpeciesID;		// particle creating track (PDG code)
    G4String		fProcessName;		// name of process when creating track
    G4String		fVolumeName;		// name of volume where track is
    G4ThreeVector	fTrackVertex;		// track vertex position
    G4String		fCreatorVolumeName;	// volume where track was created



};

typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;

extern G4Allocator<TrackerHit> TrackerHitAllocator;

inline void* TrackerHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) TrackerHitAllocator.MallocSingle();
  return aHit;
}

inline void TrackerHit::operator delete(void *aHit)
{
  TrackerHitAllocator.FreeSingle((TrackerHit*) aHit);
}

#endif
