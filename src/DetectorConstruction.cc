#include "DetectorConstruction.hh"
#include "PMTSD.hh"
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4GDMLParser.hh"
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

namespace
{
  auto &gConfMan = ConfManager::GetInstance();
}

//_____________________________________________________________________________
DetectorConstruction::DetectorConstruction()
    // : G4VUserDetectorConstruction(), m_check_overlaps(true) //debug
    : G4VUserDetectorConstruction(), m_check_overlaps(false)
{
}

//_____________________________________________________________________________
DetectorConstruction::~DetectorConstruction()
{
}

//_____________________________________________________________________________
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

//_____________________________________________________________________________
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

//_____________________________________________________________________________
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

//_____________________________________________________________________________
void DetectorConstruction::AddOpticalProperties()
{
  using CLHEP::eV;
  using CLHEP::m;
  using CLHEP::mm;

  std::vector<G4double> photon_energy, refractive_index, absorption_length;
  G4int n_entries;

  // +---------+
  // | Aerogel |
  // +---------+
  auto aerogel_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.05, 1.05};
  n_entries = photon_energy.size();
  aerogel_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  // for absorption_length
  photon_energy = {1.3 * eV, 1.56 * eV, 1.68 * eV, 1.84 * eV, 2.06 * eV, 2.26 * eV, 2.54 * eV, 2.90 * eV, 3.10 * eV, 3.28 * eV, 3.94 * eV, 4.94 * eV, 7.0 * eV};
  absorption_length = {500 * mm, 128 * mm, 120 * mm, 97 * mm, 77 * mm, 59 * mm, 41 * mm, 26 * mm, 20 * mm, 17 * mm, 8 * mm, 4 * mm, 1 * mm};
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

  // +----------------------+
  // | Black sheet Property |
  // +----------------------+
  auto blacksheet_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.6, 1.6};
  absorption_length = {1.0e-11 * m, 1.0e-11 * m};
  n_entries = photon_energy.size();
  blacksheet_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  blacksheet_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["BlackSheet"]->SetMaterialPropertiesTable(blacksheet_prop);

  // +-----------------+
  // | Teflon Property |
  // +-----------------+
  auto teflon_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.35, 1.35};
  absorption_length = {1.0e-11 * m, 1.0e-11 * m};
  n_entries = photon_energy.size();
  teflon_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  teflon_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Teflon"]->SetMaterialPropertiesTable(teflon_prop);

  // +-----------------------------+
  // | PMT window (Glass) Property |
  // +-----------------------------+
  auto pmt_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.52, 1.52};
  // absorption_length = {1.0 * cm, 1.0 * cm};
  n_entries = photon_energy.size();
  pmt_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  // pmt_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);
  m_material_map["Glass"]->SetMaterialPropertiesTable(pmt_prop);

  // +----------------------------+
  // | PMT casing (POM) Property |
  // +----------------------------+
  auto pom_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.48, 1.48};            // 一般的なPOMの屈折率
  absorption_length = {0.01 * mm, 0.01 * mm}; // 完全遮光に近い設定
  pom_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], 2);
  pom_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], 2);
  m_material_map["POM"]->SetMaterialPropertiesTable(pom_prop);
}

//_____________________________________________________________________________
void DetectorConstruction::ConstructSAC()
{
  using CLHEP::deg;
  using CLHEP::mm;

  // ----------------------
  // パラメータ定義
  // ----------------------
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
  new G4PVPlacement(nullptr, origin, gel_lv, "GelPV", mother_lv, false, 0, m_check_overlaps);
  gel_lv->SetVisAttributes(G4Colour::White());

  // ----------------------
  // Teflon Sheet
  // ----------------------
  auto sheet_solid = new G4Box("TeflonSheetSolid", gel_size.x() / 2, gel_size.y() / 2, teflon_thickness / 2);
  auto sheet_lv = new G4LogicalVolume(sheet_solid, m_material_map["Teflon"], "TeflonSheetLV");
  G4ThreeVector sheet_pos_top(0, 0, gel_size.z() / 2 + teflon_thickness / 2);
  G4ThreeVector sheet_pos_bot(0, 0, -gel_size.z() / 2 - teflon_thickness / 2);
  new G4PVPlacement(nullptr, sheet_pos_top, sheet_lv, "TeflonSheetTopPV", mother_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(nullptr, sheet_pos_bot, sheet_lv, "TeflonSheetBotPV", mother_lv, false, 1, m_check_overlaps);
  sheet_lv->SetVisAttributes(G4Colour::White());

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
  new G4PVPlacement(nullptr, origin, frame_lv, "TeflonFrameBotPV", mother_lv, false, 0, m_check_overlaps);
  frame_lv->SetVisAttributes(G4Colour::White());

  // ----------------------
  // PMT
  // ----------------------
  auto pmt_casing_solid = new G4Tubs("PMTCasingSolid", 0.0, pmt_casing_radius / 2, pmt_thickness / 2, 0.0, 360.0 * deg);
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
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i, m_check_overlaps);      // upper casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 3, m_check_overlaps); // lower casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i, m_check_overlaps);      // upper window
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 3, m_check_overlaps); // lower window
    }

    double y_pos = 0;
    new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, 6, m_check_overlaps); // left casing
    new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, 7, m_check_overlaps);  // right casing
    new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, 6, m_check_overlaps); // left window
    new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, 7, m_check_overlaps);  // right window
  }

  // new SAC
  if (pmt_channel == 14)
  {
    for (int i = 0; i < 3; ++i)
    {
      double x_pos = (i - 1) * pmt_x_spacing;
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i, m_check_overlaps);      // upper casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 3, m_check_overlaps); // lower casing
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, gel_size.y() / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i, m_check_overlaps);      // upper window
      new G4PVPlacement(rotX, G4ThreeVector(x_pos, -gel_size.y() / 2, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 3, m_check_overlaps); // lower window
    }

    for (int i = 0; i < 4; ++i)
    {
      double y_pos = (1.5 - i) * pmt_y_spacing;
      new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 6, m_check_overlaps); // left casing
      new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2, y_pos, 0), pmt_casing_lv, "PMTCasing", mother_lv, false, i + 10, m_check_overlaps); // right casing
      new G4PVPlacement(rotY, G4ThreeVector(-gel_size.x() / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 6, m_check_overlaps); // left window
      new G4PVPlacement(rotY, G4ThreeVector(gel_size.x() / 2, y_pos, 0), pmt_window_lv, "PMTWindow", mother_lv, false, i + 10, m_check_overlaps); // right window
    }
  }

  // Set sensitive detector
  auto pmt_sd = new PMTSD("PMT_SD");
  G4SDManager::GetSDMpointer()->AddNewDetector(pmt_sd);
  pmt_window_lv->SetSensitiveDetector(pmt_sd);
}

void DetectorConstruction::DumpMaterialProperties(G4Material *mat)
{
  G4cout << "=== Material: " << mat->GetName() << " ===" << G4endl;

  auto matPropTable = mat->GetMaterialPropertiesTable();
  if (!matPropTable)
  {
    G4cout << "No material properties table found." << G4endl;
    return;
  }

  std::vector<G4String> propertyNames = {"RINDEX", "ABSORPTION", "REFLECTIVITY"};

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
