#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"

class DetectorMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

private:
  std::map<G4String, G4Element *> m_element_map;
  std::map<G4String, G4Material *> m_material_map;
  G4LogicalVolume *m_world_lv;
  G4bool m_check_overlaps;
  G4OpticalSurface *gel_teflon_surf;

private:
  virtual G4VPhysicalVolume *Construct();
  void ConstructElements();
  void ConstructMaterials();
  void ConstructSAC();
  void AddOpticalProperties();
  void DumpMaterialProperties(G4Material *mat);

  void CheckOverlaps(G4bool flag) { m_check_overlaps = flag; }
};

#endif /*OpNoviceDetectorConstruction_h*/
