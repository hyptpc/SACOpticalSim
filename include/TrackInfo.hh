#ifndef TRACK_INFO_HH
#define TRACK_INFO_HH

#include <unordered_map>

#include "G4Track.hh"

class TrackInfo
{
public:
  static void ResetIfNewEvent(G4int event_id);
  static void RegisterTrack(const G4Track *track);
  static G4int GetParentPDG(G4int parent_id);

private:
  static std::unordered_map<G4int, G4int> s_track_pdg;
  static G4int s_last_event_id;
};

#endif
