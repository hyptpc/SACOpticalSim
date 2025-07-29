#ifndef PMTSD_HH
#define PMTSD_HH

#include "G4VSensitiveDetector.hh"
#include "PMTHit.hh"

#include "TSpline.h"

class G4Step;
class G4TouchableHistory;
class G4HCofThisEvent;

class PMTSD : public G4VSensitiveDetector
{
public:
  PMTSD(const G4String &name);
  ~PMTSD() override;

  void Initialize(G4HCofThisEvent *HCE) override;
  G4bool ProcessHits(G4Step *step, G4TouchableHistory *history) override;
  void EndOfEvent(G4HCofThisEvent *HCE) override;

private:
  G4THitsCollection<PMTHit> *m_hits_collection;
  TSpline3 *m_qe_spline;
  TSpline3 *m_trans_spline;
  G4double m_range_min;
  G4double m_range_max;

  void InitializeQESplines();
};

#endif
