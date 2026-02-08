#include "TrackInfo.hh"

#include "G4ParticleDefinition.hh"

std::unordered_map<G4int, G4int> TrackInfo::s_track_pdg;
G4int TrackInfo::s_last_event_id = -1;

void TrackInfo::ResetIfNewEvent(G4int event_id)
{
  if (event_id != s_last_event_id)
  {
    s_track_pdg.clear();
    s_last_event_id = event_id;
  }
}

void TrackInfo::RegisterTrack(const G4Track *track)
{
  if (!track)
    return;

  const G4int track_id = track->GetTrackID();
  if (s_track_pdg.find(track_id) != s_track_pdg.end())
    return;

  const G4ParticleDefinition *def = track->GetDefinition();
  if (!def)
    return;

  s_track_pdg.emplace(track_id, def->GetPDGEncoding());
}

G4int TrackInfo::GetParentPDG(G4int parent_id)
{
  auto it = s_track_pdg.find(parent_id);
  if (it == s_track_pdg.end())
    return 0;
  return it->second;
}
