#ifndef PRIMARYGENERATORACTION_HH
#define PRIMARYGENERATORACTION_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "TFile.h"
#include "TTree.h"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event *anEvent) override;

private:
  TFile *fBeamFile = nullptr;
  TTree *fBeamTree = nullptr;
  int fNEntries = 0;
  int fCurrentEntry = 0;
  double beam_x = 0.0;
  double beam_y = 0.0;
  G4ParticleGun *fParticleGun; // Particle gun
  void GenerateBeam(G4Event *anEvent);
};

#endif
