#include "SteppingAction.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4TouchableHandle.hh"
#include "PMTSD.hh"
#include "TrackInfo.hh"

SteppingAction::SteppingAction()
{
}

SteppingAction::~SteppingAction()
{
}

// void SteppingAction::UserSteppingAction(const G4Step* step)
// {
// }

void SteppingAction::UserSteppingAction(const G4Step *step)
{
    const auto *track = step->GetTrack();
    const auto *event = G4RunManager::GetRunManager()->GetCurrentEvent();
    if (event)
    {
        TrackInfo::ResetIfNewEvent(event->GetEventID());
    }
    if (track && track->GetCurrentStepNumber() == 1)
    {
        TrackInfo::RegisterTrack(track);
    }
}
