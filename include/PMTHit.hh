#ifndef PMT_HIT_HH
#define PMT_HIT_HH

#include "G4VHit.hh"
#include "G4ThreeVector.hh"
#include "G4Allocator.hh"
#include "G4THitsCollection.hh"

class PMTHit : public G4VHit
{
public:
  PMTHit();                    // Default constructor
  virtual ~PMTHit();           // Destructor
  PMTHit(const PMTHit &right); // Copy constructor

  // Set and get hit position (local coordinates)
  void SetPosition(const G4ThreeVector &pos) { fPosition = pos; }
  G4ThreeVector GetPosition() const { return fPosition; }

  // Set and get hit position (world coordinates)
  void SetWorldPosition(const G4ThreeVector &pos) { fWorldPosition = pos; }
  G4ThreeVector GetWorldPosition() const { return fWorldPosition; }

  // Set and get time of flight (ToF)
  void SetTime(G4double t) { fTime = t; }
  G4double GetTime() const { return fTime; }

  // Set and get energy of the detected photon
  void SetEnergy(G4double e) { fEnergy = e; }
  G4double GetEnergy() const { return fEnergy; }

  // Set and get energy of the detected photon
  void SetWaveLength(G4double wl) { fWaveLength = wl; }
  G4double GetWaveLength() const { return fWaveLength; }

  // Set and get particle ID
  void SetParticleID(G4int pid) { fParticleID = pid; }
  G4int GetParticleID() const { return fParticleID; }

  // Set and get PMT copy number (which PMT was hit)
  void SetCopyNumber(G4int cn) { fCopyNumber = cn; }
  G4int GetCopyNumber() const { return fCopyNumber; }

  // Set and get event ID
  void SetEventID(G4int id) { fEventID = id; }
  G4int GetEventID() const { return fEventID; }

  // Set and get event ID
  void SetDetectFlag(G4int detectFlag) { fDetectFlag = detectFlag; }
  G4int GetDetectFlag() const { return fDetectFlag; }

  // Set and get delta-ray origin flag (1: from delta-ray, 0: otherwise)
  void SetDeltaFlag(G4int deltaFlag) { fDeltaFlag = deltaFlag; }
  G4int GetDeltaFlag() const { return fDeltaFlag; }

  // Set and get origin volume ID for photon generation
  void SetOriginID(G4int originID) { fOriginID = originID; }
  G4int GetOriginID() const { return fOriginID; }

  void Print() const; // Print hit details

private:
  G4ThreeVector fPosition;      // Local position of the hit
  G4ThreeVector fWorldPosition; // World position of the hit
  G4double fTime;               // Time of flight (ToF)
  G4double fEnergy;             // Photon energy
  G4double fWaveLength;         // wave length
  G4int fParticleID;            // Particle ID
  G4int fCopyNumber;            // PMT copy number
  G4int fEventID;               // Event ID
  G4int fDetectFlag;            // detect flag
  G4int fDeltaFlag;             // delta-ray origin flag
  G4int fOriginID;              // origin volume ID for photon generation
};

// Memory allocator for PMTHit objects
extern G4Allocator<PMTHit> PMTHitAllocator;

#endif
