#include "DetectorConstruction.hh"
#include "PMTSD.hh"
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4SubtractionSolid.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "ConfManager.hh"
#include "G4Tubs.hh"

namespace
{
  auto &gConfMan = ConfManager::GetInstance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::DetectorConstruction()
    : G4VUserDetectorConstruction(), m_check_overlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4VPhysicalVolume *
DetectorConstruction::Construct()
{
  using CLHEP::m;

  ConstructElements();
  ConstructMaterials();
  AddOpticalProperties();

  auto world_solid = new G4Box("WorldSolid", 1. * m / 2, 1. * m / 2, 1. * m / 2);
  m_world_lv = new G4LogicalVolume(world_solid, m_material_map["Air"],
                                   "World");
  m_world_lv->SetVisAttributes(G4VisAttributes::GetInvisible());
  auto world_pv = new G4PVPlacement(nullptr, G4ThreeVector(), m_world_lv,
                                    "World", nullptr, false, 0, m_check_overlaps);

  ConstructSAC();

  return world_pv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructElements()
{
  using CLHEP::g;
  using CLHEP::mole;
  /* G4Element(name, symbol, Z, A) */
  G4String name, symbol;
  G4double Z, A;
  name = "Hydrogen";
  m_element_map[name] = new G4Element(name, symbol = "H", Z = 1.,
                                      A = 1.00794 * g / mole);
  name = "Boron";
  m_element_map[name] = new G4Element(name, symbol = "B", Z = 5.,
                                      A = 10.811 * g / mole);
  name = "Carbon";
  m_element_map[name] = new G4Element(name, symbol = "C", Z = 6.,
                                      A = 12.011 * g / mole);
  name = "Nitrogen";
  m_element_map[name] = new G4Element(name, symbol = "N", Z = 7.,
                                      A = 14.00674 * g / mole);
  name = "Oxygen";
  m_element_map[name] = new G4Element(name, symbol = "O", Z = 8.,
                                      A = 15.9994 * g / mole);
  name = "Fluorine";
  m_element_map[name] = new G4Element(name, symbol = "F", Z = 9.,
                                      A = 18.998 * g / mole);
  name = "Sodium";
  m_element_map[name] = new G4Element(name, symbol = "Na", Z = 11.,
                                      A = 22.989768 * g / mole);
  name = "Aluminum";
  m_element_map[name] = new G4Element(name, symbol = "Al", Z = 13.,
                                      A = 26.9815386 * g / mole);
  name = "Silicon";
  m_element_map[name] = new G4Element(name, symbol = "Si", Z = 14.,
                                      A = 28.0855 * g / mole);
  name = "Phoshorus";
  m_element_map[name] = new G4Element(name, symbol = "P", Z = 15.,
                                      A = 30.973762 * g / mole);
  name = "Sulfur";
  m_element_map[name] = new G4Element(name, symbol = "S", Z = 16.,
                                      A = 32.066 * g / mole);
  name = "Chlorine";
  m_element_map[name] = new G4Element(name, symbol = "Cl", Z = 17.,
                                      A = 35.453 * g / mole);
  name = "Argon";
  m_element_map[name] = new G4Element(name, symbol = "Ar", Z = 18.,
                                      A = 39.948 * g / mole);
  name = "Potassium";
  m_element_map[name] = new G4Element(name, symbol = "K", Z = 19.,
                                      A = 39.093 * g / mole);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::ConstructMaterials()
{
  using CLHEP::cm3;
  using CLHEP::g;
  using CLHEP::mg;
  using CLHEP::mole;
  using CLHEP::STP_Temperature;

  /*
    G4Material(name, density, nelement, state, temperature, pressure);
    G4Material(name, z, a, density, state, temperature, pressure);
  */
  G4String name;
  G4double Z, A, density, massfraction;
  G4int natoms, nel, ncomponents;
  const G4double room_temp = STP_Temperature + 20. * CLHEP::kelvin;

  // Vacuum
  name = "Vacuum";
  m_material_map[name] =
      new G4Material(name, density = CLHEP::universe_mean_density, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Nitrogen"], 0.7);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], 0.3);

  // Air
  name = "Air";
  m_material_map[name] = new G4Material(name, density = 1.2929e-03 * g / cm3,
                                        nel = 3, kStateGas, room_temp);
  G4double fracN = 75.47;
  G4double fracO = 23.20;
  G4double fracAr = 1.28;
  G4double denominator = fracN + fracO + fracAr;
  m_material_map[name]->AddElement(m_element_map["Nitrogen"],
                                   massfraction = fracN / denominator);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],
                                   massfraction = fracO / denominator);
  m_material_map[name]->AddElement(m_element_map["Argon"],
                                   massfraction = fracAr / denominator);

  // Kapton
  name = "Kapton";
  m_material_map[name] = new G4Material(name, density = 1.42 * g / cm3,
                                        ncomponents = 4);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"],
                                   massfraction = 0.0273);
  m_material_map[name]->AddElement(m_element_map["Carbon"],
                                   massfraction = 0.7213);
  m_material_map[name]->AddElement(m_element_map["Nitrogen"],
                                   massfraction = 0.0765);
  m_material_map[name]->AddElement(m_element_map["Oxygen"],
                                   massfraction = 0.1749);

  // Aerogel
  name = "Aerogel";
  m_material_map[name] = new G4Material(name, density = 0.2 * g / cm3, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms = 1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], natoms = 2);

  // blackSheet
  name = "BlackSheet";
  m_material_map[name] = new G4Material(name, density = 0.95 * g / cm3, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms = 1);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms = 2);

  // POM（Polyoxymethylene）
  name = "POM";
  m_material_map[name] = new G4Material(name, density = 1.41 * g / cm3, nel = 3);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms = 1);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms = 2);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], natoms = 1);

  // Epoxi
  name = "Epoxi";
  m_material_map[name] = new G4Material(name, density = 1.1 * g / cm3, nel = 4);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms = 21);
  m_material_map[name]->AddElement(m_element_map["Hydrogen"], natoms = 25);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], natoms = 5);
  m_material_map[name]->AddElement(m_element_map["Chlorine"], natoms = 1);

  // Teflon
  name = "Teflon";
  m_material_map[name] = new G4Material(name, density = 2.2 * g / cm3, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Carbon"], natoms = 2);
  m_material_map[name]->AddElement(m_element_map["Fluorine"], natoms = 4);

  // SiO2
  name = "SiO2";
  density = 2.2 * g / cm3;
  m_material_map[name] = new G4Material(name, density, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Silicon"], natoms = 1);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], natoms = 2);

  // B2O3
  name = "B2O3";
  density = 2.46 * g / cm3;
  m_material_map[name] = new G4Material(name, density, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Boron"], natoms = 2);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], natoms = 3);

  // Na2O
  name = "Na2O";
  density = 2.27 * g / cm3;
  m_material_map[name] = new G4Material(name, density, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Sodium"], natoms = 2);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], natoms = 1);

  // Al2O3
  name = "Al2O3";
  density = 3.97 * g / cm3;
  m_material_map[name] = new G4Material(name, density, nel = 2);
  m_material_map[name]->AddElement(m_element_map["Aluminum"], natoms = 2);
  m_material_map[name]->AddElement(m_element_map["Oxygen"], natoms = 3);

  // Borosilicate Glass (PMT Window)
  name = "Glass";
  density = 2.23 * g / cm3;
  m_material_map[name] = new G4Material(name, density, ncomponents = 4);
  m_material_map[name]->AddMaterial(m_material_map["SiO2"], massfraction = 0.806);
  m_material_map[name]->AddMaterial(m_material_map["B2O3"], massfraction = 0.130);
  m_material_map[name]->AddMaterial(m_material_map["Na2O"], massfraction = 0.040);
  m_material_map[name]->AddMaterial(m_material_map["Al2O3"], massfraction = 0.024);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::AddOpticalProperties()
{
  using CLHEP::eV;
  using CLHEP::m;
  using CLHEP::mm;
  using CLHEP::um;

  std::vector<G4double> photon_energy, refractive_index, absorption_length, reflectivity;
  G4int n_entries;

  // +-------------------------------------------------------------------------+
  // | Aerogel                                                                 |
  // | Ref: Hasegawa master thesis                                             |
  // | from web site http://refractiveindex.info/?group=CRYSTALS&material=SiO2 |
  // | source of the information on quartz seems to be Gorachand Ghosh,        |
  // | Dispersion-equation coefficients for the refractive index and           |
  // | birefringence of calcite and quartz crystals, Opt. Commun. 163,         |
  // | 95-102 (1999) doi:10.1016/S0030-4018(99)00091-7.                        |
  // |                                                                         |
  // | Translation (of density difference) from quartz (SiO2) to silicai       |
  // | aerogel was made by the equation:                                       |
  // | n(aerogel) = 1 + (n(SiO2)-1)*rho(aerogel)/rho(SiO2)                     |
  // +-------------------------------------------------------------------------+
  auto aerogel_prop = new G4MaterialPropertiesTable();

  photon_energy = {1.93 * eV, 1.99 * eV, 2.06 * eV, 2.13 * eV, 2.21 * eV,
                   2.29 * eV, 2.38 * eV, 2.47 * eV, 2.57 * eV, 2.69 * eV,
                   2.81 * eV, 2.94 * eV, 3.09 * eV, 3.25 * eV, 3.43 * eV,
                   3.64 * eV, 3.86 * eV, 4.12 * eV, 4.42 * eV, 4.76 * eV,
                   5.15 * eV, 5.62 * eV, 6.18 * eV, 6.98 * eV, 7.73 * eV};

  refractive_index = {1.049, 1.049, 1.049, 1.049, 1.050,
                      1.050, 1.050, 1.050, 1.050, 1.050,
                      1.050, 1.050, 1.051, 1.051, 1.051,
                      1.052, 1.052, 1.053, 1.053, 1.054,
                      1.055, 1.057, 1.059, 1.062, 1.068};

  absorption_length = {1.1750 * m, 1.1750 * m, 1.1750 * m, 1.0820 * m, 0.9704 * m,
                       0.8584 * m, 0.7592 * m, 0.6808 * m, 0.6022 * m, 0.5044 * m,
                       0.4224 * m, 0.3636 * m, 0.2884 * m, 0.2350 * m, 0.1916 * m,
                       0.1499 * m, 0.1173 * m, 0.0892 * m, 0.0698 * m, 0.0546 * m,
                       0.0439 * m, 0.0335 * m, 0.0269 * m, 0.0269 * m, 0.0269 * m};

  n_entries = photon_energy.size();
  aerogel_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);

  n_entries = photon_energy.size();
  aerogel_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Aerogel"]->SetMaterialPropertiesTable(aerogel_prop);

  // +--------------+
  // | Air Property |
  // +--------------+
  auto air_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.0, 1.0};
  n_entries = photon_energy.size();
  air_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  m_material_map["Air"]->SetMaterialPropertiesTable(air_prop);

  // +----------------------------------------------------------+
  // | Black sheet Property                                     |
  // | absorption length is set as a rough value just to ensure |
  // | that optical photon doesn't come out from the aerogel    |
  // +----------------------------------------------------------+
  auto blacksheet_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.6, 1.6};
  absorption_length = {1. * um, 1. * um};
  reflectivity = {0.0, 0.0};
  n_entries = photon_energy.size();
  blacksheet_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  blacksheet_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["BlackSheet"]->SetMaterialPropertiesTable(blacksheet_prop);

  // +----------------------------------------------------------------------------+
  // | Teflon Property                                                            |
  // | absorption length is set to a small value to ensure it doesn't             |
  // | emmit cherenkov light inside                                               |
  // | Reflectivity is based on measurment by Yamamoto-san (2015)                 |
  // | http://lambda.phys.tohoku.ac.jp/~db/human_resource/thesis/2009_B_3_M_1.pdf |
  // +----------------------------------------------------------------------------+
  auto teflon_prop = new G4MaterialPropertiesTable();

  photon_energy = {1.54 * eV, 1.76 * eV, 2.06 * eV, 2.47 * eV, 2.75 * eV,
                   2.94 * eV, 3.09 * eV, 3.17 * eV, 3.25 * eV, 3.34 * eV,
                   3.43 * eV, 3.53 * eV, 3.64 * eV, 3.75 * eV, 3.86 * eV};

  refractive_index = {1.35, 1.35, 1.35, 1.35, 1.35,
                      1.35, 1.35, 1.35, 1.35, 1.35,
                      1.35, 1.35, 1.35, 1.35, 1.35};

  absorption_length = {1. * um, 1. * um, 1. * um, 1. * um, 1. * um,
                       1. * um, 1. * um, 1. * um, 1. * um, 1. * um,
                       1. * um, 1. * um, 1. * um, 1. * um, 1. * um};

  n_entries = photon_energy.size();
  teflon_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  teflon_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Teflon"]->SetMaterialPropertiesTable(teflon_prop);

  // Optical surface settings: Aerogel ↔ Teflon sheet
  const G4int teflon_layer = gConfMan.GetInt("teflon_layer");
  const G4int SigmaAlpha = gConfMan.GetInt("SigmaAlpha");
  if (teflon_layer == 2)
  {
    reflectivity = {0.85, 0.90, 0.93, 0.97, 0.98,
                    0.99, 0.99, 1.00, 1.00, 1.00,
                    1.00, 1.00, 1.00, 1.00, 1.00};
  }
  else if (teflon_layer == 3)
  {
    reflectivity = {0.90, 0.93, 0.96, 0.98, 1.00,
                    0.99, 0.99, 1.00, 1.00, 1.00,
                    1.00, 1.00, 1.00, 1.00, 1.00};
  }

  auto gel_steflon_prop = new G4MaterialPropertiesTable();
  gel_steflon_prop->AddProperty("REFLECTIVITY", &photon_energy[0], &reflectivity[0], n_entries);
  gel_steflon_surf = new G4OpticalSurface("GelTeflonSheetSurface");
  gel_steflon_surf->SetType(dielectric_dielectric);
  gel_steflon_surf->SetModel(unified);
  gel_steflon_surf->SetFinish(groundbackpainted);
  gel_steflon_surf->SetSigmaAlpha(SigmaAlpha); // roughness
  gel_steflon_surf->SetMaterialPropertiesTable(gel_steflon_prop);

  // Optical surface settings: Aerogel ↔ Teflon frame
  reflectivity = {0.97, 0.98, 0.99, 1.00, 1.00,
                  0.99, 0.99, 1.00, 1.00, 1.00,
                  1.00, 1.00, 1.00, 1.00, 1.00};

  auto gel_fteflon_prop = new G4MaterialPropertiesTable();
  gel_fteflon_prop->AddProperty("REFLECTIVITY", &photon_energy[0], &reflectivity[0], n_entries);
  gel_fteflon_surf = new G4OpticalSurface("GelTeflonFrameSurface");
  gel_fteflon_surf->SetType(dielectric_dielectric);
  gel_fteflon_surf->SetModel(unified);
  gel_fteflon_surf->SetFinish(groundfrontpainted);
  gel_fteflon_surf->SetSigmaAlpha(SigmaAlpha); // roughness
  gel_fteflon_surf->SetMaterialPropertiesTable(gel_fteflon_prop);

  // +-----------------------------------------------------------------+
  // | PMT window (Glass) Property                                     |
  // | Ref: https://refractiveindex.info/?shelf=3d&book=glass&page=BK7 |
  // +-----------------------------------------------------------------+
  auto pmt_prop = new G4MaterialPropertiesTable();

  photon_energy = {1.76 * eV, 1.82 * eV, 1.87 * eV, 1.93 * eV, 1.99 * eV,
                   2.06 * eV, 2.13 * eV, 2.21 * eV, 2.29 * eV, 2.38 * eV,
                   2.47 * eV, 2.57 * eV, 2.69 * eV, 2.81 * eV, 2.94 * eV,
                   3.09 * eV, 3.25 * eV, 3.43 * eV, 3.64 * eV, 3.86 * eV,
                   4.12 * eV};

  refractive_index = {1.513, 1.514, 1.514, 1.515, 1.516,
                      1.516, 1.517, 1.518, 1.519, 1.520,
                      1.521, 1.523, 1.524, 1.526, 1.528,
                      1.531, 1.534, 1.537, 1.541, 1.546,
                      1.553};

  absorption_length = {6.24 * m, 5.02 * m, 4.15 * m, 4.15 * m, 4.15 * m,
                       4.52 * m, 4.99 * m, 5.64 * m, 5.88 * m, 4.90 * m,
                       4.15 * m, 3.85 * m, 3.56 * m, 3.18 * m, 3.56 * m,
                       3.11 * m, 1.46 * m, 5.33e-1 * m, 1.35e-1 * m, 3.82e-2 * m,
                       8.35e-3 * m};

  n_entries = photon_energy.size();
  pmt_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  pmt_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Glass"]->SetMaterialPropertiesTable(pmt_prop);

  // +------------------------------------------------------------------------------+
  // | PMT casing (POM) Property                                                    |
  // | Ref: https://www.sciencedirect.com/science/article/abs/pii/S1350449523003687 |
  // +------------------------------------------------------------------------------+
  auto pom_prop = new G4MaterialPropertiesTable();

  photon_energy = {1.79 * eV, 1.87 * eV, 1.96 * eV, 2.06 * eV, 2.17 * eV,
                   2.29 * eV, 2.43 * eV, 2.58 * eV, 2.75 * eV, 2.95 * eV,
                   3.18 * eV, 3.44 * eV, 3.77 * eV, 4.15 * eV, 4.61 * eV,
                   5.20 * eV, 5.97 * eV};

  refractive_index = {1.480, 1.480, 1.480, 1.480, 1.480,
                      1.481, 1.481, 1.482, 1.482, 1.482,
                      1.484, 1.485, 1.488, 1.492, 1.497,
                      1.509, 1.535};

  n_entries = photon_energy.size();
  pom_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);

  photon_energy = {1.79 * eV, 1.87 * eV, 1.97 * eV, 2.06 * eV, 2.17 * eV,
                   2.29 * eV, 2.43 * eV, 2.59 * eV, 2.76 * eV, 2.94 * eV,
                   3.18 * eV, 3.44 * eV, 3.78 * eV, 4.15 * eV, 4.61 * eV,
                   5.20 * eV, 5.89 * eV};

  absorption_length = {36.358 * nm, 34.729 * nm, 33.021 * nm, 31.504 * nm, 29.907 * nm,
                       28.293 * nm, 26.683 * nm, 24.989 * nm, 23.386 * nm, 21.959 * nm,
                       20.312 * nm, 18.733 * nm, 17.041 * nm, 15.519 * nm, 13.939 * nm,
                       12.362 * nm, 10.872 * nm};

  n_entries = photon_energy.size();
  pom_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["POM"]->SetMaterialPropertiesTable(pom_prop);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSAC()
{
  using CLHEP::deg;
  using CLHEP::mm;

  const G4ThreeVector gel_size(
      gConfMan.GetDouble("gel_size_x") * mm,
      gConfMan.GetDouble("gel_size_y") * mm,
      gConfMan.GetDouble("gel_size_z") * mm);

  const G4double teflon_thickness = gConfMan.GetDouble("teflon_thickness") * mm;
  const G4double BlackSheet_thickness = gConfMan.GetDouble("BlackSheet_thickness") * mm;
  const G4double frame_thickness = gConfMan.GetDouble("frame_thickness") * mm;
  const G4double pmt_x_spacing = gConfMan.GetDouble("pmt_x_spacing") * mm;
  const G4double pmt_y_spacing = gConfMan.GetDouble("pmt_y_spacing") * mm;
  const G4double pmt_casing_radius = gConfMan.GetDouble("pmt_casing_radius") * mm;
  const G4double pmt_window_radius = gConfMan.GetDouble("pmt_window_radius") * mm;
  const G4double pmt_thickness = gConfMan.GetDouble("pmt_thickness") * mm;
  const G4int pmt_channel = gConfMan.GetInt("pmt_channel");
  const G4ThreeVector origin(0, 0, 0);

  // ----------------------
  // Mother Volume（Air）
  // ----------------------
  auto mother_solid = new G4Box("SACMotherSolid",
                                gel_size.x() / 2 + 50 * mm,
                                gel_size.y() / 2 + 50 * mm,
                                gel_size.z() / 2 + 50 * mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "SACMotherLV");
  new G4PVPlacement(nullptr, origin, mother_lv, "SACMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  // ----------------------
  // Aerogel
  // ----------------------
  auto gel_solid = new G4Box("GelSolid", gel_size.x() / 2, gel_size.y() / 2, gel_size.z() / 2);
  auto gel_lv = new G4LogicalVolume(gel_solid, m_material_map["Aerogel"], "GelLV");
  auto gel_pv = new G4PVPlacement(nullptr, origin, gel_lv, "GelPV", mother_lv, false, 0, m_check_overlaps);
  // gel_lv->SetVisAttributes(G4Colour::Cyan());

  // ----------------------
  // Teflon Sheet
  // ----------------------
  auto sheet_solid = new G4Box("TeflonSheetSolid", gel_size.x() / 2, gel_size.y() / 2, teflon_thickness / 2);
  auto sheet_lv = new G4LogicalVolume(sheet_solid, m_material_map["Teflon"], "TeflonSheetLV");
  G4ThreeVector sheet_pos_top(0, 0, gel_size.z() / 2 + teflon_thickness / 2);
  G4ThreeVector sheet_pos_bot(0, 0, -gel_size.z() / 2 - teflon_thickness / 2);
  auto sheet_top_pv = new G4PVPlacement(nullptr, sheet_pos_top, sheet_lv, "TeflonSheetTopPV", mother_lv, false, 0, m_check_overlaps);
  auto sheet_bot_pv = new G4PVPlacement(nullptr, sheet_pos_bot, sheet_lv, "TeflonSheetBotPV", mother_lv, false, 1, m_check_overlaps);
  sheet_lv->SetVisAttributes(G4Colour::White());
  new G4LogicalBorderSurface("Gel_TeflonTop", gel_pv, sheet_top_pv, gel_steflon_surf);
  new G4LogicalBorderSurface("Gel_TeflonBot", gel_pv, sheet_bot_pv, gel_steflon_surf);

  // ----------------------
  // BlackSheet
  // ----------------------
  auto black_sheet_solid = new G4Box("BlackSheetSolid", gel_size.x() / 2, gel_size.y() / 2, BlackSheet_thickness / 2);
  auto black_sheet_lv = new G4LogicalVolume(black_sheet_solid, m_material_map["BlackSheet"], "BlackSheetLV");
  G4ThreeVector black_top_pos(0, 0, gel_size.z() / 2 + teflon_thickness + BlackSheet_thickness / 2);
  G4ThreeVector black_bot_pos(0, 0, -gel_size.z() / 2 - teflon_thickness - BlackSheet_thickness / 2);
  new G4PVPlacement(nullptr, black_top_pos, black_sheet_lv, "BlackSheetTopPV", mother_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(nullptr, black_bot_pos, black_sheet_lv, "BlackSheetBotPV", mother_lv, false, 1, m_check_overlaps);
  black_sheet_lv->SetVisAttributes(G4Colour::Gray());

  // ----------------------
  // Teflon Frame
  // ----------------------
  auto frame_outer = new G4Box("FrameOuter", gel_size.x() / 2 + frame_thickness, gel_size.y() / 2 + frame_thickness, gel_size.z() / 2);
  auto frame_inner = new G4Box("FrameInner", gel_size.x() / 2, gel_size.y() / 2, gel_size.z() / 2);
  G4SubtractionSolid *frame = new G4SubtractionSolid("Frame", frame_outer, frame_inner);

  // substract holes for PMTs
  auto hole = new G4Tubs("hole", 0.0, pmt_casing_radius / 2, frame_thickness / 2, 0.0, 360.0 * deg);

  auto rotX = new G4RotationMatrix();
  rotX->rotateX(90. * deg);
  auto rotY = new G4RotationMatrix();
  rotY->rotateY(90. * deg);

  // old SAC
  if (pmt_channel == 8)
  {
    for (int i = 0; i < 3; ++i)
    {
      double x_pos = (i - 1) * pmt_x_spacing;
      auto upper_trans = G4Transform3D(*rotX, G4ThreeVector(x_pos, gel_size.y() / 2 + frame_thickness / 2, 0));
      auto lower_trans = G4Transform3D(*rotX, G4ThreeVector(x_pos, -gel_size.y() / 2 - frame_thickness / 2, 0));
      frame = new G4SubtractionSolid("Frame", frame, hole, upper_trans);
      frame = new G4SubtractionSolid("Frame", frame, hole, lower_trans);
    }

    double y_pos = 0;
    auto left_trans = G4Transform3D(*rotY, G4ThreeVector(-gel_size.x() / 2 - frame_thickness / 2, y_pos, 0));
    auto right_trans = G4Transform3D(*rotY, G4ThreeVector(gel_size.x() / 2 + frame_thickness / 2, y_pos, 0));
    frame = new G4SubtractionSolid("Frame", frame, hole, left_trans);
    frame = new G4SubtractionSolid("Frame", frame, hole, right_trans);
  }

  // new SAC
  if (pmt_channel == 14)
  {
    for (int i = 0; i < 3; ++i)
    {
      double x_pos = (i - 1) * pmt_x_spacing;
      auto upper_trans = G4Transform3D(*rotX, G4ThreeVector(x_pos, gel_size.y() / 2 + frame_thickness / 2, 0));
      auto lower_trans = G4Transform3D(*rotX, G4ThreeVector(x_pos, -gel_size.y() / 2 - frame_thickness / 2, 0));
      frame = new G4SubtractionSolid("Frame", frame, hole, upper_trans);
      frame = new G4SubtractionSolid("Frame", frame, hole, lower_trans);
    }

    for (int i = 0; i < 4; ++i)
    {
      double y_pos = (1.5 - i) * pmt_y_spacing;
      auto left_trans = G4Transform3D(*rotY, G4ThreeVector(-gel_size.x() / 2 - frame_thickness / 2, y_pos, 0));
      auto right_trans = G4Transform3D(*rotY, G4ThreeVector(gel_size.x() / 2 + frame_thickness / 2, y_pos, 0));
      frame = new G4SubtractionSolid("Frame", frame, hole, left_trans);
      frame = new G4SubtractionSolid("Frame", frame, hole, right_trans);
    }
  }

  auto frame_lv = new G4LogicalVolume(frame, m_material_map["Teflon"], "TeflonFrameLV");
  auto frame_pv = new G4PVPlacement(nullptr, origin, frame_lv, "TeflonFramePV", mother_lv, false, 0, m_check_overlaps);
  frame_lv->SetVisAttributes(G4Colour::White());
  new G4LogicalBorderSurface("Gel_TeflonFrame", gel_pv, frame_pv, gel_fteflon_surf);

  // ----------------------
  // PMT (Window and Casing)
  // ----------------------
  auto pmt_casing_solid = new G4Tubs("PMTCasingSolid", pmt_window_radius / 2, pmt_casing_radius / 2, pmt_thickness / 2, 0.0, 360.0 * deg);
  auto pmt_window_solid = new G4Tubs("PMTWindowSolid", 0.0, pmt_window_radius / 2, pmt_thickness / 2, 0.0, 360.0 * deg);

  auto pmt_casing_lv = new G4LogicalVolume(pmt_casing_solid, m_material_map["POM"], "PMTCasingLV");
  auto pmt_window_lv = new G4LogicalVolume(pmt_window_solid, m_material_map["Glass"], "PMTWindowLV");

  pmt_window_lv->SetVisAttributes(G4VisAttributes(G4Colour::Yellow()));

  // old SAC
  if (pmt_channel == 8)
  {
    for (int i = 0; i < 3; ++i)
    {
      double x_pos = (i - 1) * pmt_x_spacing;
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2 + pmt_thickness / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i, m_check_overlaps);      // upper casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2 - pmt_thickness / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 3, m_check_overlaps); // lower casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2 + pmt_thickness / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i, m_check_overlaps);      // upper window
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2 - pmt_thickness / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 3, m_check_overlaps); // lower window
    }

    double y_pos = 0;
    new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2 - pmt_thickness / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, 6, m_check_overlaps); // left casing
    new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2 + pmt_thickness / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, 7, m_check_overlaps);  // right casing
    new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2 - pmt_thickness / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, 6, m_check_overlaps); // left window
    new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2 + pmt_thickness / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, 7, m_check_overlaps);  // right window
  }

  // new SAC
  if (pmt_channel == 14)
  {
    for (int i = 0; i < 3; ++i)
    {
      double x_pos = (i - 1) * pmt_x_spacing;
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2 + pmt_thickness / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i, m_check_overlaps);      // upper casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2 - pmt_thickness / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 3, m_check_overlaps); // lower casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2 + pmt_thickness / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i, m_check_overlaps);      // upper window
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2 - pmt_thickness / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 3, m_check_overlaps); // lower window
    }

    for (int i = 0; i < 4; ++i)
    {
      double y_pos = (1.5 - i) * pmt_y_spacing;
      new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2 - pmt_thickness / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 6, m_check_overlaps); // left casing
      new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2 + pmt_thickness / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 10, m_check_overlaps); // right casing
      new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2 - pmt_thickness / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 6, m_check_overlaps); // left window
      new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2 + pmt_thickness / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 10, m_check_overlaps); // right window
    }
  }

  // Set sensitive detector
  auto pmt_sd = new PMTSD("PMT_SD");
  G4SDManager::GetSDMpointer()->AddNewDetector(pmt_sd);
  pmt_window_lv->SetSensitiveDetector(pmt_sd);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void DetectorConstruction::DumpMaterialProperties(G4Material *mat)
{
  G4cout << "=== Material: " << mat->GetName() << " ===" << G4endl;

  auto matPropTable = mat->GetMaterialPropertiesTable();
  if (!matPropTable)
  {
    G4cout << "No material properties table found." << G4endl;
    return;
  }

  std::vector<G4String> propertyNames = {"RINDEX", "ABSLENGTH", "REFLECTIVITY"};

  for (const auto &prop : propertyNames)
  {
    if (matPropTable->ConstPropertyExists(prop))
    {
      G4cout << prop << ": " << matPropTable->GetConstProperty(prop) << G4endl;
    }
  }

  for (const auto &prop : propertyNames)
  {
    if (matPropTable->GetProperty(prop))
    {
      G4MaterialPropertyVector *mpv = matPropTable->GetProperty(prop);
      G4cout << prop << ":" << G4endl;
      for (size_t i = 0; i < mpv->GetVectorLength(); i++)
      {
        G4cout << "  Energy: " << mpv->Energy(i) / eV << " eV, "
               << " Value: " << (*mpv)[i] << G4endl;
      }
    }
  }
}
