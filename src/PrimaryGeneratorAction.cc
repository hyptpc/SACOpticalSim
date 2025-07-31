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
  G4double gel_z = gConfMan.GetDouble("gel_size_z") * mm;
  G4double teflon_thickness = gConfMan.GetDouble("teflon_thickness") * mm;
  G4double blacksheet_thickness = gConfMan.GetDouble("BlackSheet_thickness") * mm;

  G4double x = 0.0 * mm;
  G4double y = 0.0 * mm;
  G4double z = -gel_z / 2.0 - teflon_thickness - blacksheet_thickness - 10.0 * mm;

  G4ThreeVector position(x, y, z);
  fParticleGun->SetParticlePosition(position);
  gAnaMan.SetBeamPosition(position);

  // -----------------------
  // Debug
  // -----------------------
  // G4cout << "Particle: " << particle->GetParticleName() << G4endl
  // 	 << " | Energy: " << energy / GeV << " GeV" << G4endl
  // 	 << " | Momentum: " << momentum / GeV << " GeV/c" << G4endl
  // 	 << " | Position: (" << x / mm << ", " << y / mm << ", " << z / mm << ") mm"  << G4endl
  // 	 << " | Direction: (" << direction.x() << ", " << direction.y() << ", " << direction.z() << ")"
  // 	 << G4endl;

  // -----------------------
  // Gus
  // -----------------------
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
