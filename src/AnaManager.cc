#include "AnaManager.hh"
#include "ConfManager.hh"
#include "G4Run.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"

#include "PMTHit.hh"

#include "Randomize.hh"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

#include <string>
#include <sstream>
#include <vector>

#include "G4ThreeVector.hh"

extern int gCerenkovCounter;
extern double decay_check;

AnaManager &AnaManager::GetInstance()
{
  static AnaManager instance;
  return instance;
}

AnaManager::AnaManager()
    : m_file(),
      m_output_rootfile_path("test.root"),
      m_tree(new TTree("tree", "GEANT4 optical simulation for SAC")),
      m_evnum(0),
      m_nhit_pmt(0),
      m_cerenkov_all(0),
      m_cerenkov_aerogel(0),
      m_beam_energy(0.),
      m_beam_mom_x(0.),
      m_beam_mom_y(0.),
      m_beam_mom_z(0.),
      m_beam_pos_x(0.),
      m_beam_pos_y(0.),
      m_beam_pos_z(0.)
{
}

AnaManager::~AnaManager()
{
}

//_____________________________________________________________________________
void AnaManager::BeginOfRunAction(const G4Run *)
{
  m_file = new TFile(m_output_rootfile_path, "RECREATE");
  m_tree->Reset();

  m_tree->Branch("evnum", &m_evnum, "evnum/I");
  m_tree->Branch("cerenkov_all", &m_cerenkov_all, "cerenkov_all/I");
  m_tree->Branch("cerenkov_aerogel", &m_cerenkov_aerogel, "cerenkov_aerogel/I");

  // -- beam -----
  m_tree->Branch("beam_energy", &m_beam_energy, "beam_energy/D");
  m_tree->Branch("beam_mom_x", &m_beam_mom_x, "beam_mom_x/D");
  m_tree->Branch("beam_mom_y", &m_beam_mom_y, "beam_mom_y/D");
  m_tree->Branch("beam_mom_z", &m_beam_mom_z, "beam_mom_z/D");
  m_tree->Branch("beam_pos_x", &m_beam_pos_x, "beam_pos_x/D");
  m_tree->Branch("beam_pos_y", &m_beam_pos_y, "beam_pos_y/D");
  m_tree->Branch("beam_pos_z", &m_beam_pos_z, "beam_pos_z/D");

  // -- PMT -----
  m_tree->Branch("nhit_pmt", &m_nhit_pmt, "nhit_pmt/I");
  m_tree->Branch("pos_x", &m_pos_x);
  m_tree->Branch("pos_y", &m_pos_y);
  m_tree->Branch("pos_z", &m_pos_z);
  m_tree->Branch("time", &m_time);
  m_tree->Branch("energy", &m_energy);
  m_tree->Branch("wave_length", &m_wave_length);
  m_tree->Branch("particle_id", &m_particle_id);
  m_tree->Branch("seg", &m_seg);
  m_tree->Branch("detect_flag", &m_detect_flag);
}

void AnaManager::BeginOfEventAction(const G4Event *anEvent)
{
}

void AnaManager::EndOfEventAction(const G4Event *anEvent)
{
  G4HCofThisEvent *HCTE = anEvent->GetHCofThisEvent();
  if (!HCTE)
    return;
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();

  m_nhit_pmt = 0;
  G4THitsCollection<PMTHit> *PMTHC;
  G4int ColIdPMT = SDMan->GetCollectionID("PmtCollection");
  if (ColIdPMT >= 0)
  {
    PMTHC = dynamic_cast<G4THitsCollection<PMTHit> *>(HCTE->GetHC(ColIdPMT));
    if (PMTHC)
    {
      m_nhit_pmt = PMTHC->entries();
    }
  }

  ResetContainer();
  for (int i = 0; i < m_nhit_pmt; i++)
  {
    PMTHit *aHit = (*PMTHC)[i];

    G4ThreeVector pos = aHit->GetPosition();
    // m_pos.push_back(TVector3(pos.x(), pos.y(), pos.z()));
    m_pos_x.push_back(pos.x());
    m_pos_y.push_back(pos.y());
    m_pos_z.push_back(pos.z());

    G4double time = aHit->GetTime();
    m_time.push_back(time);

    G4double energy = aHit->GetEnergy();
    m_energy.push_back(energy);

    G4double wave_length = aHit->GetWaveLength();
    m_wave_length.push_back(wave_length);

    G4int particle_id = aHit->GetParticleID();
    m_particle_id.push_back(particle_id);

    G4int seg = aHit->GetCopyNumber();
    m_seg.push_back(seg);

    G4int detect_flag = aHit->GetDetectFlag();
    m_detect_flag.push_back(detect_flag);
  }

  m_tree->Fill();
  m_evnum++;
  G4cout << m_evnum << ", " << m_nhit_pmt << G4endl;
}

void AnaManager::EndOfRunAction(const G4Run *aRun)
{
  if (m_file && m_file->IsOpen())
  {
    m_file->cd();
    m_tree->Write();
    m_file->Close();
  }
}

void AnaManager::ResetContainer()
{
  // m_pos.clear();
  m_pos_x.clear();
  m_pos_y.clear();
  m_pos_z.clear();
  m_time.clear();
  m_energy.clear();
  m_wave_length.clear();
  m_particle_id.clear();
  m_seg.clear();
  m_detect_flag.clear();
}

void AnaManager::SetNumOfCerenkovAll(G4int cerenkov_all)
{
  m_cerenkov_all = cerenkov_all;
}

void AnaManager::SetNumOfCerenkovAerogel(G4int cerenkov_aerogel)
{
  m_cerenkov_aerogel = cerenkov_aerogel;
}

void AnaManager::SetBeamEnergy(G4double beam_energy)
{
  m_beam_energy = beam_energy;
}

void AnaManager::SetBeamMomentum(G4ThreeVector beam_momentum)
{
  m_beam_mom_x = beam_momentum.x();
  m_beam_mom_y = beam_momentum.y();
  m_beam_mom_z = beam_momentum.z();
}

void AnaManager::SetBeamPosition(G4ThreeVector beam_position)
{
  m_beam_pos_x = beam_position.x();
  m_beam_pos_y = beam_position.y();
  m_beam_pos_z = beam_position.z();
}

void AnaManager::SetOutputRootfilePath(G4String output_rootfile_path)
{
  m_output_rootfile_path = output_rootfile_path;
}
G4String AnaManager::GetOutputRootfilePath()
{
  return m_output_rootfile_path;
}
