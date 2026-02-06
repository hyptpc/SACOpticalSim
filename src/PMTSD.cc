// Reference: https://kumaroot.readthedocs.io/ja/latest/geant4/geant4-qe.html

#include "PMTSD.hh"
#include "PMTHit.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "Randomize.hh"

#include "TGraph.h"
#include "TSpline.h"

PMTSD::PMTSD(const G4String &name)
    : G4VSensitiveDetector(name),
      m_qe_spline(nullptr),
      m_trans_spline(nullptr)
{
  collectionName.insert("PmtCollection");
  InitializeQESplines();
}

PMTSD::~PMTSD()
{
  delete m_qe_spline;
  delete m_trans_spline;
}

//_____________________________________________________________________________
void PMTSD::Initialize(G4HCofThisEvent *HCTE)
{
  m_hits_collection = new G4THitsCollection<PMTHit>(SensitiveDetectorName, collectionName[0]);
  HCTE->AddHitsCollection(GetCollectionID(0), m_hits_collection);
}

//_____________________________________________________________________________
G4bool PMTSD::ProcessHits(G4Step *aStep, G4TouchableHistory *)
{
  const auto preStepPoint = aStep->GetPreStepPoint();
  const auto aTrack = aStep->GetTrack();

  // Collect only Cherenkov photons
  if (aTrack->GetDefinition() != G4OpticalPhoton::Definition())
  {
    return false;
  }

  // Energy
  G4double energy = aTrack->GetKineticEnergy(); // exclude rest mass

  // Calculate Quantum Efficiency
  G4double eff_qe = 0.0;
  if (energy >= m_qe_spline->GetXmin() && energy <= m_qe_spline->GetXmax())
  {
    eff_qe = m_qe_spline->Eval(energy);
  }

  // Detection flag
  G4int detectFlag = 0;
  if (G4UniformRand() < eff_qe)
    detectFlag = 1;

  // Stop the track
  aTrack->SetTrackStatus(fStopAndKill);

  // Hit info
  G4ThreeVector worldPos = preStepPoint->GetPosition();
  G4ThreeVector pos = preStepPoint->GetTouchable()->GetHistory()->GetTopTransform().TransformPoint(worldPos);
  G4double hitTime = preStepPoint->GetGlobalTime();
  G4double waveLength = (CLHEP::h_Planck * CLHEP::c_light / energy) / CLHEP::nm;
  G4int copyNumber = preStepPoint->GetTouchableHandle()->GetCopyNumber();
  G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();
  G4int particleID = aTrack->GetDefinition()->GetPDGEncoding();
  G4int parentID = aTrack->GetParentID();
  G4int deltaFlag = (parentID > 1) ? 1 : 0;
  G4int originID = 0;
  const auto *origin_lv = aTrack->GetLogicalVolumeAtVertex();
  if (origin_lv)
  {
    const auto &name = origin_lv->GetName();
    if (name == "GelLV")
      originID = 1; // Aerogel
    else if (name == "TeflonSheetLV")
      originID = 2; // Teflon sheet
    else if (name == "TeflonFrameLV")
      originID = 3; // Teflon frame
    else if (name == "BlackSheetLV")
      originID = 4; // Black sheet
    else
      originID = 9; // Other
  }

  // Create hit
  auto *aHit = new PMTHit();
  aHit->SetPosition(pos);
  aHit->SetWorldPosition(worldPos);
  aHit->SetEnergy(energy);
  aHit->SetWaveLength(waveLength);
  aHit->SetTime(hitTime);
  aHit->SetCopyNumber(copyNumber);
  aHit->SetEventID(eventID);
  aHit->SetDetectFlag(detectFlag);
  aHit->SetParticleID(particleID);
  aHit->SetDeltaFlag(deltaFlag);
  aHit->SetOriginID(originID);

  m_hits_collection->insert(aHit);
  return true;
}

//_____________________________________________________________________________
void PMTSD::EndOfEvent(G4HCofThisEvent *) {}

//_____________________________________________________________________________
void PMTSD::InitializeQESplines()
{
  using namespace CLHEP;

  // ---- QE ----
  std::vector<double> photonEnergyQE = {1.74 * eV, 1.76 * eV, 1.77 * eV, 1.78 * eV, 1.79 * eV,
                                        1.81 * eV, 1.82 * eV, 1.83 * eV, 1.84 * eV, 1.85 * eV,
                                        1.86 * eV, 1.87 * eV, 1.89 * eV, 1.90 * eV, 1.92 * eV,
                                        1.93 * eV, 1.95 * eV, 1.97 * eV, 1.99 * eV, 2.02 * eV,
                                        2.05 * eV, 2.08 * eV, 2.11 * eV, 2.15 * eV, 2.20 * eV,
                                        2.24 * eV, 2.29 * eV, 2.35 * eV, 2.41 * eV, 2.49 * eV,
                                        2.58 * eV, 2.69 * eV, 2.82 * eV, 2.96 * eV, 3.12 * eV,
                                        3.30 * eV, 3.49 * eV, 3.68 * eV, 3.86 * eV, 3.99 * eV,
                                        4.09 * eV, 4.16 * eV, 4.22 * eV, 4.27 * eV, 4.31 * eV,
                                        4.37 * eV, 4.41 * eV, 4.45 * eV, 4.50 * eV, 4.56 * eV,
                                        4.60 * eV, 4.65 * eV, 4.69 * eV, 4.74 * eV};

  std::vector<double> QE = {0.0002, 0.0002, 0.0003, 0.0004, 0.0005,
                            0.0007, 0.0009, 0.0011, 0.0015, 0.0018,
                            0.0024, 0.0031, 0.0039, 0.0050, 0.0063,
                            0.0081, 0.0105, 0.0132, 0.0169, 0.0216,
                            0.0271, 0.0336, 0.0422, 0.0514, 0.0656,
                            0.0799, 0.0974, 0.1206, 0.1470, 0.1711,
                            0.1933, 0.2086, 0.2217, 0.2320, 0.2392,
                            0.2392, 0.2320, 0.2054, 0.1738, 0.1404,
                            0.1101, 0.0863, 0.0666, 0.0522, 0.0409,
                            0.0321, 0.0251, 0.0197, 0.0152, 0.0119,
                            0.0092, 0.0072, 0.0056, 0.0046};

  auto *qe_graph = new TGraph(photonEnergyQE.size(), &photonEnergyQE[0], &QE[0]);

  m_qe_spline = new TSpline3("qe_spline", qe_graph);

  delete qe_graph;
}
