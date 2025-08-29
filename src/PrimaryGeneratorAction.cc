#include "PrimaryGeneratorAction.hh"
#include "AnaManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4UnitsTable.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"
#include "TFile.h"
#include "TTree.h"
#include "ConfManager.hh"

namespace
{
  using CLHEP::deg;
  using CLHEP::GeV;
  using CLHEP::mm;
  const auto particleTable = G4ParticleTable::GetParticleTable();
  auto &gAnaMan = AnaManager::GetInstance();
  auto &gConfMan = ConfManager::GetInstance();
}

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
{
  fParticleGun = new G4ParticleGun(1);

  static const G4String beamfile = "../conf/BeamProfile/" + gConfMan.Get("beamfile");
  fBeamFile = TFile::Open(beamfile);
  fBeamTree = (TTree *)fBeamFile->Get("beam");
  fBeamTree->SetBranchAddress("x", &beam_x);
  fBeamTree->SetBranchAddress("y", &beam_y);
  fNEntries = fBeamTree->GetEntries();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  GenerateBeam(anEvent);
  // GeneratePhoton(anEvent);
}

void PrimaryGeneratorAction::GenerateBeam(G4Event *anEvent)
{
  static const G4String particle_name = gConfMan.Get("particle");
  static const auto particle = particleTable->FindParticle(particle_name);
  fParticleGun->SetParticleDefinition(particle);

  // -----------------------
  // Momentum
  // -----------------------
  G4double p0 = gConfMan.GetDouble("momentum") * GeV;
  G4double sigma_p = p0 * 0.02 / 2.355;
  G4double momentum = G4RandGauss::shoot(p0, sigma_p);

  G4double mass = particle->GetPDGMass();
  G4double energy = std::sqrt(mass * mass + momentum * momentum);
  gAnaMan.SetBeamEnergy(energy);
  fParticleGun->SetParticleEnergy(energy);

  // -----------------------
  // Momentum direction
  // -----------------------
  G4ThreeVector direction(0., 0., 1.);

  // fParticleGun->SetParticleMomentum(momentum);
  fParticleGun->SetParticleMomentumDirection(direction);

  // -----------------------
  // Position
  // -----------------------
  if (fCurrentEntry < fNEntries)
  {
    fBeamTree->GetEntry(fCurrentEntry++);
  }
  else
  {
    G4cerr << "[PrimaryGeneratorAction] Error: fCurrentEntry (" << fCurrentEntry << ") >= fNEntries (" << fNEntries << ")" << G4endl;
    G4Exception("PrimaryGeneratorAction::GenerateBeam", "BeamEntryOverflow", FatalException, "Number of beam profile entries exceeded.");
  }

  G4double gel_z = gConfMan.GetDouble("gel_size_z") * mm;
  G4double teflon_thickness = gConfMan.GetDouble("teflon_thickness") * mm;
  G4double blacksheet_thickness = gConfMan.GetDouble("BlackSheet_thickness") * mm;

  G4double x = beam_x * mm;
  G4double y = beam_y * mm;
  G4double z = -gel_z / 2.0 - teflon_thickness - blacksheet_thickness - 10.0 * mm;

  G4ThreeVector position(x, y, z);
  fParticleGun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

  // -----------------------
  // Gus
  // -----------------------
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
