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

  // Optical photon check
  if (aTrack->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
    return false;

  // Energy
  G4double energy = aTrack->GetKineticEnergy(); // exclude rest mass

  // Calculate the effective Quantum Efficiency
  G4double eff_qe = 0.0;
  if (energy >= m_qe_spline->GetXmin() && energy <= m_qe_spline->GetXmax() &&
      energy >= m_trans_spline->GetXmin() && energy <= m_trans_spline->GetXmax())
  {
    eff_qe = m_qe_spline->Eval(energy) * m_trans_spline->Eval(energy);
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
  std::vector<double> photonEnergyQE = {
      1.5000 * eV, 1.8000 * eV, 1.8300 * eV, 1.8600 * eV, 1.8900 * eV,
      1.9200 * eV, 1.9500 * eV, 1.9800 * eV, 2.0100 * eV, 2.0400 * eV,
      2.0700 * eV, 2.1000 * eV, 2.1300 * eV, 2.1600 * eV, 2.1900 * eV,
      2.2200 * eV, 2.2500 * eV, 2.2800 * eV, 2.3100 * eV, 2.3400 * eV,
      2.3700 * eV, 2.4000 * eV, 2.4300 * eV, 2.4600 * eV, 2.4900 * eV,
      2.5200 * eV, 2.5500 * eV, 2.5800 * eV, 2.6100 * eV, 2.6400 * eV,
      2.6700 * eV, 2.7000 * eV, 2.7300 * eV, 2.7600 * eV, 2.7900 * eV,
      2.8200 * eV, 2.8500 * eV, 2.8800 * eV, 2.9100 * eV, 2.9700 * eV,
      3.0000 * eV, 3.0300 * eV, 3.0600 * eV, 3.0900 * eV, 3.1200 * eV,
      3.1800 * eV, 3.2100 * eV, 3.2400 * eV, 3.3600 * eV, 3.3900 * eV,
      3.4200 * eV, 3.4800 * eV, 3.5100 * eV, 3.5400 * eV, 3.6000 * eV,
      3.6300 * eV, 3.6600 * eV, 3.7200 * eV, 3.7500 * eV, 3.8100 * eV,
      3.8400 * eV, 3.9000 * eV, 3.9300 * eV, 3.9900 * eV, 4.0200 * eV,
      4.0800 * eV, 4.1100 * eV, 4.1700 * eV, 4.2300 * eV, 4.2600 * eV,
      4.3200 * eV, 4.3500 * eV, 4.4100 * eV};
  std::vector<double> QE = {
      0.0007, 0.0017, 0.0024, 0.0057, 0.0081, 0.0106, 0.0156, 0.0205, 0.0253, 0.0358,
      0.0431, 0.0505, 0.0578, 0.0745, 0.0838, 0.0912, 0.0985, 0.1071, 0.1243, 0.1329,
      0.1427, 0.1526, 0.1619, 0.1712, 0.1784, 0.1855, 0.1933, 0.2011, 0.2108, 0.2205,
      0.2284, 0.2363, 0.2418, 0.2474, 0.2518, 0.2561, 0.2603, 0.2645, 0.2663, 0.2682,
      0.2697, 0.2713, 0.2719, 0.2725, 0.2744, 0.2763, 0.2770, 0.2776, 0.2744, 0.2713,
      0.2682, 0.2651, 0.2621, 0.2591, 0.2532, 0.2474, 0.2365, 0.2256, 0.2134, 0.2011,
      0.1923, 0.1834, 0.1698, 0.1561, 0.1415, 0.1269, 0.1116, 0.0963, 0.0792, 0.0621,
      0.0398, 0.0175, 0.0000};

  // ---- Transmittance ----
  std::vector<double> photonEnergyT = {
      1.5000 * eV, 3.1246 * eV, 3.1754 * eV, 3.2008 * eV, 3.2262 * eV,
      3.3531 * eV, 3.3785 * eV, 3.4292 * eV, 3.4800 * eV, 3.5054 * eV,
      3.5562 * eV, 3.5815 * eV, 3.6323 * eV, 3.6831 * eV, 3.7085 * eV,
      3.7592 * eV, 3.8100 * eV, 3.8354 * eV, 3.8862 * eV, 3.9369 * eV,
      3.9877 * eV, 4.0131 * eV, 4.0638 * eV, 4.1146 * eV, 4.1654 * eV,
      4.2162 * eV, 4.2669 * eV, 4.3177 * eV, 4.3685 * eV, 4.4192 * eV,
      4.4700 * eV, 4.5208 * eV, 4.5715 * eV, 4.6223 * eV, 4.6731 * eV,
      4.6985 * eV};
  std::vector<double> Trans = {
      0.9950, 0.9441, 0.9333, 0.9226, 0.9120, 0.8915, 0.8710, 0.8611, 0.8511, 0.8415,
      0.8318, 0.8130, 0.7943, 0.7678, 0.7413, 0.7246, 0.7079, 0.6843, 0.6607, 0.6316,
      0.6026, 0.5698, 0.5370, 0.5024, 0.4677, 0.4284, 0.3890, 0.3490, 0.3090, 0.2745,
      0.2399, 0.1992, 0.1585, 0.1094, 0.0603, 0.0000};

  auto *qe_graph = new TGraph(photonEnergyQE.size(), &photonEnergyQE[0], &QE[0]);
  auto *t_graph = new TGraph(photonEnergyT.size(), &photonEnergyT[0], &Trans[0]);

  m_qe_spline = new TSpline3("qe_spline", qe_graph);
  m_trans_spline = new TSpline3("trans_spline", t_graph);

  delete qe_graph;
  delete t_graph;
}