#include "DetectorConstruction.hh"
#include "MPPCSD.hh"
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
    : G4VUserDetectorConstruction(), m_check_overlaps(true)
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

  // blacksheet
  name = "Blacksheet";
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

  // +-----------------+
  // | Quartz Property |
  // +-----------------+
  auto quartz_prop = new G4MaterialPropertiesTable();

  std::vector<G4double> photon_energy{
      1.0022 * eV, 1.0045 * eV, 1.0069 * eV, 1.0092 * eV, 1.0115 * eV, 1.0138 * eV,
      1.0162 * eV, 1.0185 * eV, 1.0209 * eV, 1.0232 * eV, 1.0256 * eV, 1.0279 * eV,
      1.0303 * eV, 1.0327 * eV, 1.0351 * eV, 1.0374 * eV, 1.0398 * eV, 1.0422 * eV,
      1.0446 * eV, 1.0470 * eV, 1.0495 * eV, 1.0519 * eV, 1.0543 * eV, 1.0567 * eV,
      1.0592 * eV, 1.0616 * eV, 1.0641 * eV, 1.0665 * eV, 1.0690 * eV, 1.0714 * eV,
      1.0739 * eV, 1.0764 * eV, 1.0789 * eV, 1.0813 * eV, 1.0838 * eV, 1.0863 * eV,
      1.0888 * eV, 1.0914 * eV, 1.0939 * eV, 1.0964 * eV, 1.0989 * eV, 1.1015 * eV,
      1.1040 * eV, 1.1065 * eV, 1.1091 * eV, 1.1116 * eV, 1.1142 * eV, 1.1168 * eV,
      1.1194 * eV, 1.1219 * eV, 1.1245 * eV, 1.1271 * eV, 1.1297 * eV, 1.1323 * eV,
      1.1349 * eV, 1.1375 * eV, 1.1402 * eV, 1.1428 * eV, 1.1454 * eV, 1.1481 * eV,
      1.1507 * eV, 1.1534 * eV, 1.1560 * eV, 1.1587 * eV, 1.1614 * eV, 1.1640 * eV,
      1.1667 * eV, 1.1694 * eV, 1.1721 * eV, 1.1748 * eV, 1.1775 * eV, 1.1802 * eV,
      1.1830 * eV, 1.1857 * eV, 1.1884 * eV, 1.1911 * eV, 1.1939 * eV, 1.1966 * eV,
      1.1994 * eV, 1.2022 * eV, 1.2049 * eV, 1.2077 * eV, 1.2105 * eV, 1.2133 * eV,
      1.2161 * eV, 1.2189 * eV, 1.2217 * eV, 1.2245 * eV, 1.2273 * eV, 1.2302 * eV,
      1.2330 * eV, 1.2359 * eV, 1.2387 * eV, 1.2416 * eV, 1.2444 * eV, 1.2473 * eV,
      1.2502 * eV, 1.2530 * eV, 1.2559 * eV, 1.2588 * eV, 1.2617 * eV, 1.2646 * eV,
      1.2676 * eV, 1.2705 * eV, 1.2734 * eV, 1.2763 * eV, 1.2793 * eV, 1.2822 * eV,
      1.2852 * eV, 1.2882 * eV, 1.2911 * eV, 1.2941 * eV, 1.2971 * eV, 1.3001 * eV,
      1.3031 * eV, 1.3061 * eV, 1.3091 * eV, 1.3121 * eV, 1.3151 * eV, 1.3182 * eV,
      1.3212 * eV, 1.3242 * eV, 1.3273 * eV, 1.3304 * eV, 1.3334 * eV, 1.3365 * eV,
      1.3396 * eV, 1.3427 * eV, 1.3458 * eV, 1.3489 * eV, 1.3520 * eV, 1.3551 * eV,
      1.3582 * eV, 1.3613 * eV, 1.3645 * eV, 1.3676 * eV, 1.3708 * eV, 1.3739 * eV,
      1.3771 * eV, 1.3803 * eV, 1.3835 * eV, 1.3866 * eV, 1.3898 * eV, 1.3930 * eV,
      1.3963 * eV, 1.3995 * eV, 1.4027 * eV, 1.4059 * eV, 1.4092 * eV, 1.4124 * eV,
      1.4157 * eV, 1.4189 * eV, 1.4222 * eV, 1.4255 * eV, 1.4288 * eV, 1.4321 * eV,
      1.4354 * eV, 1.4387 * eV, 1.4420 * eV, 1.4453 * eV, 1.4487 * eV, 1.4520 * eV,
      1.4553 * eV, 1.4587 * eV, 1.4621 * eV, 1.4654 * eV, 1.4688 * eV, 1.4722 * eV,
      1.4756 * eV, 1.4790 * eV, 1.4824 * eV, 1.4858 * eV, 1.4892 * eV, 1.4927 * eV,
      1.4961 * eV, 1.4996 * eV, 1.5030 * eV, 1.5065 * eV, 1.5100 * eV, 1.5134 * eV,
      1.5169 * eV, 1.5204 * eV, 1.5239 * eV, 1.5274 * eV, 1.5310 * eV, 1.5345 * eV,
      1.5380 * eV, 1.5416 * eV, 1.5451 * eV, 1.5487 * eV, 1.5523 * eV, 1.5558 * eV,
      1.5594 * eV, 1.5630 * eV, 1.5666 * eV, 1.5702 * eV, 1.5739 * eV, 1.5775 * eV,
      1.5811 * eV, 1.5848 * eV, 1.5884 * eV, 1.5921 * eV, 1.5958 * eV, 1.5994 * eV,
      1.6031 * eV, 1.6068 * eV, 1.6105 * eV, 1.6142 * eV, 1.6180 * eV, 1.6217 * eV,
      1.6254 * eV, 1.6292 * eV, 1.6329 * eV, 1.6367 * eV, 1.6405 * eV, 1.6442 * eV,
      1.6480 * eV, 1.6518 * eV, 1.6556 * eV, 1.6595 * eV, 1.6633 * eV, 1.6671 * eV,
      1.6710 * eV, 1.6748 * eV, 1.6787 * eV, 1.6825 * eV, 1.6864 * eV, 1.6903 * eV,
      1.6942 * eV, 1.6981 * eV, 1.7020 * eV, 1.7060 * eV, 1.7099 * eV, 1.7138 * eV,
      1.7178 * eV, 1.7217 * eV, 1.7257 * eV, 1.7297 * eV, 1.7337 * eV, 1.7377 * eV,
      1.7417 * eV, 1.7457 * eV, 1.7497 * eV, 1.7537 * eV, 1.7578 * eV, 1.7618 * eV,
      1.7659 * eV, 1.7700 * eV, 1.7741 * eV, 1.7781 * eV, 1.7822 * eV, 1.7863 * eV,
      1.7905 * eV, 1.7946 * eV, 1.7987 * eV, 1.8029 * eV, 1.8070 * eV, 1.8112 * eV,
      1.8154 * eV, 1.8196 * eV, 1.8238 * eV, 1.8280 * eV, 1.8322 * eV, 1.8364 * eV,
      1.8406 * eV, 1.8449 * eV, 1.8491 * eV, 1.8534 * eV, 1.8577 * eV, 1.8619 * eV,
      1.8662 * eV, 1.8705 * eV, 1.8748 * eV, 1.8792 * eV, 1.8835 * eV, 1.8878 * eV,
      1.8922 * eV, 1.8966 * eV, 1.9009 * eV, 1.9053 * eV, 1.9097 * eV, 1.9141 * eV,
      1.9185 * eV, 1.9229 * eV, 1.9274 * eV, 1.9318 * eV, 1.9363 * eV, 1.9407 * eV,
      1.9452 * eV, 1.9497 * eV, 1.9542 * eV, 1.9587 * eV, 1.9632 * eV, 1.9677 * eV,
      1.9723 * eV, 1.9768 * eV, 1.9814 * eV, 1.9859 * eV, 1.9905 * eV, 1.9951 * eV,
      1.9997 * eV, 2.0043 * eV, 2.0089 * eV, 2.0136 * eV, 2.0182 * eV, 2.0229 * eV,
      2.0275 * eV, 2.0322 * eV, 2.0369 * eV, 2.0416 * eV, 2.0463 * eV, 2.0510 * eV,
      2.0557 * eV, 2.0605 * eV, 2.0652 * eV, 2.0700 * eV, 2.0748 * eV, 2.0795 * eV,
      2.0843 * eV, 2.0891 * eV, 2.0939 * eV, 2.0988 * eV, 2.1036 * eV, 2.1085 * eV,
      2.1133 * eV, 2.1182 * eV, 2.1231 * eV, 2.1280 * eV, 2.1329 * eV, 2.1378 * eV,
      2.1427 * eV, 2.1477 * eV, 2.1526 * eV, 2.1576 * eV, 2.1626 * eV, 2.1675 * eV,
      2.1725 * eV, 2.1775 * eV, 2.1826 * eV, 2.1876 * eV, 2.1926 * eV, 2.1977 * eV,
      2.2028 * eV, 2.2078 * eV, 2.2129 * eV, 2.2180 * eV, 2.2231 * eV, 2.2283 * eV,
      2.2334 * eV, 2.2385 * eV, 2.2437 * eV, 2.2489 * eV, 2.2541 * eV, 2.2593 * eV,
      2.2645 * eV, 2.2697 * eV, 2.2749 * eV, 2.2802 * eV, 2.2854 * eV, 2.2907 * eV,
      2.2960 * eV, 2.3013 * eV, 2.3066 * eV, 2.3119 * eV, 2.3172 * eV, 2.3226 * eV,
      2.3279 * eV, 2.3333 * eV, 2.3387 * eV, 2.3440 * eV, 2.3494 * eV, 2.3549 * eV,
      2.3603 * eV, 2.3657 * eV, 2.3712 * eV, 2.3767 * eV, 2.3821 * eV, 2.3876 * eV,
      2.3931 * eV, 2.3986 * eV, 2.4042 * eV, 2.4097 * eV, 2.4153 * eV, 2.4208 * eV,
      2.4264 * eV, 2.4320 * eV, 2.4376 * eV, 2.4432 * eV, 2.4489 * eV, 2.4545 * eV,
      2.4602 * eV, 2.4659 * eV, 2.4715 * eV, 2.4772 * eV, 2.4829 * eV, 2.4887 * eV,
      2.4944 * eV, 2.5002 * eV, 2.5059 * eV, 2.5117 * eV, 2.5175 * eV, 2.5233 * eV,
      2.5291 * eV, 2.5349 * eV, 2.5408 * eV, 2.5466 * eV, 2.5525 * eV, 2.5584 * eV,
      2.5643 * eV, 2.5702 * eV, 2.5761 * eV, 2.5821 * eV, 2.5880 * eV, 2.5940 * eV,
      2.6000 * eV, 2.6060 * eV, 2.6120 * eV, 2.6180 * eV, 2.6240 * eV, 2.6301 * eV,
      2.6361 * eV, 2.6422 * eV, 2.6483 * eV, 2.6544 * eV, 2.6605 * eV, 2.6667 * eV,
      2.6728 * eV, 2.6790 * eV, 2.6851 * eV, 2.6913 * eV, 2.6975 * eV, 2.7037 * eV,
      2.7100 * eV, 2.7162 * eV, 2.7225 * eV, 2.7288 * eV, 2.7351 * eV, 2.7414 * eV,
      2.7477 * eV, 2.7540 * eV, 2.7604 * eV, 2.7667 * eV, 2.7731 * eV, 2.7795 * eV,
      2.7859 * eV, 2.7923 * eV, 2.7988 * eV, 2.8052 * eV, 2.8117 * eV, 2.8182 * eV,
      2.8247 * eV, 2.8312 * eV, 2.8377 * eV, 2.8442 * eV, 2.8508 * eV, 2.8574 * eV,
      2.8640 * eV, 2.8706 * eV, 2.8772 * eV, 2.8838 * eV, 2.8905 * eV, 2.8971 * eV,
      2.9038 * eV, 2.9105 * eV, 2.9172 * eV, 2.9239 * eV, 2.9307 * eV, 2.9374 * eV,
      2.9442 * eV, 2.9510 * eV, 2.9578 * eV, 2.9646 * eV, 2.9714 * eV, 2.9783 * eV,
      2.9852 * eV, 2.9920 * eV, 2.9989 * eV, 3.0058 * eV, 3.0128 * eV, 3.0197 * eV,
      3.0267 * eV, 3.0337 * eV, 3.0406 * eV, 3.0477 * eV, 3.0547 * eV, 3.0617 * eV,
      3.0688 * eV, 3.0759 * eV, 3.0829 * eV, 3.0901 * eV, 3.0972 * eV, 3.1043 * eV,
      3.1115 * eV, 3.1187 * eV, 3.1258 * eV, 3.1330 * eV, 3.1403 * eV, 3.1475 * eV,
      3.1548 * eV, 3.1620 * eV, 3.1693 * eV, 3.1766 * eV, 3.1839 * eV, 3.1913 * eV,
      3.1987 * eV, 3.2060 * eV, 3.2134 * eV, 3.2208 * eV, 3.2282 * eV, 3.2357 * eV,
      3.2431 * eV, 3.2506 * eV, 3.2581 * eV, 3.2656 * eV, 3.2732 * eV, 3.2807 * eV,
      3.2883 * eV, 3.2958 * eV, 3.3034 * eV, 3.3111 * eV, 3.3187 * eV, 3.3263 * eV,
      3.3340 * eV, 3.3417 * eV, 3.3494 * eV, 3.3571 * eV, 3.3649 * eV, 3.3726 * eV,
      3.3804 * eV, 3.3882 * eV, 3.3960 * eV, 3.4038 * eV, 3.4117 * eV, 3.4195 * eV,
      3.4274 * eV, 3.4353 * eV, 3.4432 * eV, 3.4512 * eV, 3.4591 * eV, 3.4671 * eV,
      3.4751 * eV, 3.4831 * eV, 3.4911 * eV, 3.4992 * eV, 3.5072 * eV, 3.5153 * eV,
      3.5234 * eV, 3.5316 * eV, 3.5397 * eV, 3.5479 * eV, 3.5560 * eV, 3.5642 * eV,
      3.5725 * eV, 3.5807 * eV, 3.5889 * eV, 3.5972 * eV, 3.6055 * eV, 3.6138 * eV,
      3.6222 * eV, 3.6305 * eV, 3.6389 * eV, 3.6473 * eV, 3.6557 * eV, 3.6641 * eV,
      3.6725 * eV, 3.6810 * eV, 3.6895 * eV, 3.6980 * eV, 3.7065 * eV, 3.7151 * eV,
      3.7236 * eV, 3.7322 * eV, 3.7408 * eV, 3.7494 * eV, 3.7581 * eV, 3.7667 * eV,
      3.7754 * eV, 3.7841 * eV, 3.7929 * eV, 3.8016 * eV, 3.8104 * eV, 3.8192 * eV,
      3.8279 * eV, 3.8368 * eV, 3.8456 * eV, 3.8545 * eV, 3.8634 * eV, 3.8723 * eV,
      3.8812 * eV, 3.8902 * eV, 3.8991 * eV, 3.9081 * eV, 3.9171 * eV, 3.9261 * eV,
      3.9352 * eV, 3.9443 * eV, 3.9534 * eV, 3.9625 * eV, 3.9716 * eV, 3.9808 * eV,
      3.9899 * eV, 3.9991 * eV, 4.0084 * eV, 4.0176 * eV, 4.0269 * eV, 4.0361 * eV,
      4.0455 * eV, 4.0548 * eV, 4.0641 * eV, 4.0735 * eV, 4.0829 * eV, 4.0923 * eV,
      4.1017 * eV, 4.1112 * eV, 4.1207 * eV, 4.1301 * eV, 4.1397 * eV, 4.1492 * eV,
      4.1588 * eV, 4.1684 * eV, 4.1780 * eV, 4.1876 * eV, 4.1973 * eV, 4.2069 * eV,
      4.2166 * eV, 4.2264 * eV, 4.2361 * eV, 4.2459 * eV, 4.2557 * eV, 4.2655 * eV,
      4.2753 * eV, 4.2852 * eV, 4.2950 * eV, 4.3049 * eV, 4.3149 * eV, 4.3248 * eV,
      4.3348 * eV, 4.3448 * eV, 4.3548 * eV, 4.3648 * eV, 4.3749 * eV, 4.3850 * eV,
      4.3951 * eV, 4.4052 * eV, 4.4154 * eV, 4.4255 * eV, 4.4357 * eV, 4.4460 * eV,
      4.4562 * eV, 4.4665 * eV, 4.4768 * eV, 4.4871 * eV, 4.4974 * eV, 4.5078 * eV,
      4.5182 * eV, 4.5286 * eV, 4.5391 * eV, 4.5495 * eV, 4.5600 * eV, 4.5705 * eV,
      4.5811 * eV, 4.5916 * eV, 4.6022 * eV, 4.6128 * eV, 4.6234 * eV, 4.6341 * eV,
      4.6448 * eV, 4.6555 * eV, 4.6662 * eV, 4.6770 * eV, 4.6878 * eV, 4.6986 * eV,
      4.7094 * eV, 4.7203 * eV, 4.7312 * eV, 4.7420 * eV, 4.7530 * eV, 4.7640 * eV,
      4.7749 * eV, 4.7859 * eV, 4.7970 * eV, 4.8080 * eV, 4.8191 * eV, 4.8302 * eV,
      4.8414 * eV, 4.8525 * eV, 4.8637 * eV, 4.8749 * eV, 4.8862 * eV, 4.8974 * eV,
      4.9087 * eV, 4.9200 * eV, 4.9314 * eV, 4.9427 * eV, 4.9541 * eV, 4.9655 * eV,
      4.9770 * eV, 4.9885 * eV, 5.0000 * eV, 5.0115 * eV, 5.0230 * eV, 5.0346 * eV,
      5.0462 * eV, 5.0579 * eV, 5.0695 * eV, 5.0812 * eV, 5.0929 * eV, 5.1046 * eV,
      5.1164 * eV, 5.1282 * eV, 5.1400 * eV, 5.1519 * eV, 5.1638 * eV, 5.1757 * eV,
      5.1876 * eV, 5.1996 * eV, 5.2115 * eV, 5.2236 * eV, 5.2356 * eV, 5.2477 * eV,
      5.2598 * eV, 5.2719 * eV, 5.2840 * eV, 5.2962 * eV, 5.3084 * eV, 5.3207 * eV,
      5.3329 * eV, 5.3452 * eV, 5.3576 * eV, 5.3699 * eV, 5.3823 * eV, 5.3947 * eV,
      5.4071 * eV, 5.4196 * eV, 5.4321 * eV, 5.4446 * eV, 5.4571 * eV, 5.4697 * eV,
      5.4823 * eV, 5.4950 * eV, 5.5076 * eV, 5.5203 * eV, 5.5331 * eV, 5.5458 * eV,
      5.5586 * eV, 5.5714 * eV, 5.5843 * eV, 5.5972 * eV, 5.6100 * eV, 5.6230 * eV,
      5.6360 * eV, 5.6489 * eV, 5.6619 * eV, 5.6750 * eV, 5.6881 * eV, 5.7012 * eV,
      5.7143 * eV, 5.7275 * eV, 5.7407 * eV, 5.7540 * eV, 5.7672 * eV, 5.7805 * eV,
      5.7938 * eV, 5.8072 * eV, 5.8206 * eV, 5.8340 * eV, 5.8475 * eV, 5.8609 * eV,
      5.8744 * eV, 5.8880 * eV, 5.9016 * eV, 5.9152 * eV, 5.9288 * eV, 5.9425 * eV,
      5.9562 * eV, 5.9699 * eV, 5.9836 * eV, 5.9975 * eV, 6.0113 * eV, 6.0251 * eV,
      6.0390 * eV, 6.0529 * eV, 6.0669 * eV, 6.0809 * eV, 6.0949 * eV, 6.1090 * eV,
      6.1230 * eV, 6.1371 * eV, 6.1513 * eV, 6.1655 * eV, 6.1797 * eV, 6.1939 * eV,
      6.2082 * eV, 6.2225 * eV, 6.2369 * eV, 6.2513 * eV, 6.2657 * eV, 6.2801 * eV,
      6.2946 * eV, 6.3091 * eV, 6.3236 * eV, 6.3382 * eV, 6.3528 * eV, 6.3675 * eV,
      6.3822 * eV, 6.3968 * eV, 6.4116 * eV, 6.4264 * eV, 6.4412 * eV, 6.4560 * eV,
      6.4709 * eV, 6.4859 * eV, 6.5008 * eV, 6.5158 * eV, 6.5308 * eV, 6.5458 * eV,
      6.5609 * eV, 6.5761 * eV, 6.5912 * eV, 6.6064 * eV, 6.6216 * eV, 6.6369 * eV,
      6.6522 * eV, 6.6675 * eV, 6.6829 * eV, 6.6983 * eV, 6.7138 * eV, 6.7292 * eV,
      6.7448 * eV, 6.7603 * eV, 6.7759 * eV, 6.7915 * eV, 6.8072 * eV, 6.8229 * eV,
      6.8386 * eV, 6.8543 * eV, 6.8701 * eV, 6.8860 * eV, 6.9019 * eV, 6.9178 * eV,
      6.9337 * eV, 6.9497 * eV, 6.9657 * eV, 6.9818 * eV, 6.9979 * eV};

  std::vector<G4double> refractive_index{
      1.4479, 1.4479, 1.4480, 1.4480, 1.4480, 1.4481, 1.4481, 1.4481,
      1.4482, 1.4482, 1.4482, 1.4482, 1.4483, 1.4483, 1.4483, 1.4484,
      1.4484, 1.4484, 1.4485, 1.4485, 1.4485, 1.4486, 1.4486, 1.4486,
      1.4486, 1.4487, 1.4487, 1.4487, 1.4488, 1.4488, 1.4488, 1.4489,
      1.4489, 1.4489, 1.4490, 1.4490, 1.4490, 1.4490, 1.4491, 1.4491,
      1.4491, 1.4492, 1.4492, 1.4492, 1.4493, 1.4493, 1.4493, 1.4493,
      1.4494, 1.4494, 1.4494, 1.4495, 1.4495, 1.4495, 1.4496, 1.4496,
      1.4496, 1.4496, 1.4497, 1.4497, 1.4497, 1.4498, 1.4498, 1.4498,
      1.4498, 1.4499, 1.4499, 1.4499, 1.4500, 1.4500, 1.4500, 1.4500,
      1.4501, 1.4501, 1.4501, 1.4502, 1.4502, 1.4502, 1.4503, 1.4503,
      1.4503, 1.4503, 1.4504, 1.4504, 1.4504, 1.4505, 1.4505, 1.4505,
      1.4505, 1.4506, 1.4506, 1.4506, 1.4507, 1.4507, 1.4507, 1.4507,
      1.4508, 1.4508, 1.4508, 1.4509, 1.4509, 1.4509, 1.4509, 1.4510,
      1.4510, 1.4510, 1.4511, 1.4511, 1.4511, 1.4511, 1.4512, 1.4512,
      1.4512, 1.4513, 1.4513, 1.4513, 1.4513, 1.4514, 1.4514, 1.4514,
      1.4515, 1.4515, 1.4515, 1.4516, 1.4516, 1.4516, 1.4516, 1.4517,
      1.4517, 1.4517, 1.4518, 1.4518, 1.4518, 1.4518, 1.4519, 1.4519,
      1.4519, 1.4520, 1.4520, 1.4520, 1.4520, 1.4521, 1.4521, 1.4521,
      1.4522, 1.4522, 1.4522, 1.4523, 1.4523, 1.4523, 1.4523, 1.4524,
      1.4524, 1.4524, 1.4525, 1.4525, 1.4525, 1.4526, 1.4526, 1.4526,
      1.4526, 1.4527, 1.4527, 1.4527, 1.4528, 1.4528, 1.4528, 1.4529,
      1.4529, 1.4529, 1.4529, 1.4530, 1.4530, 1.4530, 1.4531, 1.4531,
      1.4531, 1.4532, 1.4532, 1.4532, 1.4533, 1.4533, 1.4533, 1.4533,
      1.4534, 1.4534, 1.4534, 1.4535, 1.4535, 1.4535, 1.4536, 1.4536,
      1.4536, 1.4537, 1.4537, 1.4537, 1.4538, 1.4538, 1.4538, 1.4539,
      1.4539, 1.4539, 1.4540, 1.4540, 1.4540, 1.4541, 1.4541, 1.4541,
      1.4541, 1.4542, 1.4542, 1.4542, 1.4543, 1.4543, 1.4543, 1.4544,
      1.4544, 1.4544, 1.4545, 1.4545, 1.4545, 1.4546, 1.4546, 1.4547,
      1.4547, 1.4547, 1.4548, 1.4548, 1.4548, 1.4549, 1.4549, 1.4549,
      1.4550, 1.4550, 1.4550, 1.4551, 1.4551, 1.4551, 1.4552, 1.4552,
      1.4552, 1.4553, 1.4553, 1.4554, 1.4554, 1.4554, 1.4555, 1.4555,
      1.4555, 1.4556, 1.4556, 1.4556, 1.4557, 1.4557, 1.4558, 1.4558,
      1.4558, 1.4559, 1.4559, 1.4559, 1.4560, 1.4560, 1.4561, 1.4561,
      1.4561, 1.4562, 1.4562, 1.4562, 1.4563, 1.4563, 1.4564, 1.4564,
      1.4564, 1.4565, 1.4565, 1.4566, 1.4566, 1.4566, 1.4567, 1.4567,
      1.4568, 1.4568, 1.4568, 1.4569, 1.4569, 1.4570, 1.4570, 1.4570,
      1.4571, 1.4571, 1.4572, 1.4572, 1.4573, 1.4573, 1.4573, 1.4574,
      1.4574, 1.4575, 1.4575, 1.4576, 1.4576, 1.4576, 1.4577, 1.4577,
      1.4578, 1.4578, 1.4579, 1.4579, 1.4579, 1.4580, 1.4580, 1.4581,
      1.4581, 1.4582, 1.4582, 1.4583, 1.4583, 1.4584, 1.4584, 1.4584,
      1.4585, 1.4585, 1.4586, 1.4586, 1.4587, 1.4587, 1.4588, 1.4588,
      1.4589, 1.4589, 1.4590, 1.4590, 1.4591, 1.4591, 1.4592, 1.4592,
      1.4593, 1.4593, 1.4594, 1.4594, 1.4595, 1.4595, 1.4596, 1.4596,
      1.4597, 1.4597, 1.4598, 1.4598, 1.4599, 1.4599, 1.4600, 1.4600,
      1.4601, 1.4601, 1.4602, 1.4602, 1.4603, 1.4603, 1.4604, 1.4605,
      1.4605, 1.4606, 1.4606, 1.4607, 1.4607, 1.4608, 1.4608, 1.4609,
      1.4610, 1.4610, 1.4611, 1.4611, 1.4612, 1.4612, 1.4613, 1.4614,
      1.4614, 1.4615, 1.4615, 1.4616, 1.4617, 1.4617, 1.4618, 1.4618,
      1.4619, 1.4620, 1.4620, 1.4621, 1.4621, 1.4622, 1.4623, 1.4623,
      1.4624, 1.4624, 1.4625, 1.4626, 1.4626, 1.4627, 1.4628, 1.4628,
      1.4629, 1.4630, 1.4630, 1.4631, 1.4632, 1.4632, 1.4633, 1.4634,
      1.4634, 1.4635, 1.4636, 1.4636, 1.4637, 1.4638, 1.4638, 1.4639,
      1.4640, 1.4640, 1.4641, 1.4642, 1.4643, 1.4643, 1.4644, 1.4645,
      1.4645, 1.4646, 1.4647, 1.4648, 1.4648, 1.4649, 1.4650, 1.4651,
      1.4651, 1.4652, 1.4653, 1.4654, 1.4654, 1.4655, 1.4656, 1.4657,
      1.4657, 1.4658, 1.4659, 1.4660, 1.4661, 1.4661, 1.4662, 1.4663,
      1.4664, 1.4665, 1.4665, 1.4666, 1.4667, 1.4668, 1.4669, 1.4670,
      1.4670, 1.4671, 1.4672, 1.4673, 1.4674, 1.4675, 1.4676, 1.4676,
      1.4677, 1.4678, 1.4679, 1.4680, 1.4681, 1.4682, 1.4683, 1.4684,
      1.4684, 1.4685, 1.4686, 1.4687, 1.4688, 1.4689, 1.4690, 1.4691,
      1.4692, 1.4693, 1.4694, 1.4695, 1.4696, 1.4697, 1.4698, 1.4699,
      1.4700, 1.4701, 1.4702, 1.4703, 1.4704, 1.4705, 1.4706, 1.4707,
      1.4708, 1.4709, 1.4710, 1.4711, 1.4712, 1.4713, 1.4714, 1.4715,
      1.4716, 1.4717, 1.4718, 1.4719, 1.4720, 1.4721, 1.4722, 1.4724,
      1.4725, 1.4726, 1.4727, 1.4728, 1.4729, 1.4730, 1.4731, 1.4733,
      1.4734, 1.4735, 1.4736, 1.4737, 1.4738, 1.4740, 1.4741, 1.4742,
      1.4743, 1.4744, 1.4746, 1.4747, 1.4748, 1.4749, 1.4751, 1.4752,
      1.4753, 1.4754, 1.4756, 1.4757, 1.4758, 1.4759, 1.4761, 1.4762,
      1.4763, 1.4765, 1.4766, 1.4767, 1.4769, 1.4770, 1.4771, 1.4773,
      1.4774, 1.4775, 1.4777, 1.4778, 1.4779, 1.4781, 1.4782, 1.4784,
      1.4785, 1.4787, 1.4788, 1.4789, 1.4791, 1.4792, 1.4794, 1.4795,
      1.4797, 1.4798, 1.4800, 1.4801, 1.4803, 1.4804, 1.4806, 1.4807,
      1.4809, 1.4810, 1.4812, 1.4814, 1.4815, 1.4817, 1.4818, 1.4820,
      1.4822, 1.4823, 1.4825, 1.4827, 1.4828, 1.4830, 1.4832, 1.4833,
      1.4835, 1.4837, 1.4838, 1.4840, 1.4842, 1.4843, 1.4845, 1.4847,
      1.4849, 1.4851, 1.4852, 1.4854, 1.4856, 1.4858, 1.4860, 1.4861,
      1.4863, 1.4865, 1.4867, 1.4869, 1.4871, 1.4873, 1.4875, 1.4876,
      1.4878, 1.4880, 1.4882, 1.4884, 1.4886, 1.4888, 1.4890, 1.4892,
      1.4894, 1.4896, 1.4898, 1.4900, 1.4902, 1.4905, 1.4907, 1.4909,
      1.4911, 1.4913, 1.4915, 1.4917, 1.4919, 1.4922, 1.4924, 1.4926,
      1.4928, 1.4930, 1.4933, 1.4935, 1.4937, 1.4940, 1.4942, 1.4944,
      1.4946, 1.4949, 1.4951, 1.4954, 1.4956, 1.4958, 1.4961, 1.4963,
      1.4966, 1.4968, 1.4970, 1.4973, 1.4975, 1.4978, 1.4981, 1.4983,
      1.4986, 1.4988, 1.4991, 1.4993, 1.4996, 1.4999, 1.5001, 1.5004,
      1.5007, 1.5009, 1.5012, 1.5015, 1.5018, 1.5020, 1.5023, 1.5026,
      1.5029, 1.5032, 1.5034, 1.5037, 1.5040, 1.5043, 1.5046, 1.5049,
      1.5052, 1.5055, 1.5058, 1.5061, 1.5064, 1.5067, 1.5070, 1.5073,
      1.5076, 1.5080, 1.5083, 1.5086, 1.5089, 1.5092, 1.5096, 1.5099,
      1.5102, 1.5105, 1.5109, 1.5112, 1.5115, 1.5119, 1.5122, 1.5126,
      1.5129, 1.5133, 1.5136, 1.5140, 1.5143, 1.5147, 1.5150, 1.5154,
      1.5158, 1.5161, 1.5165, 1.5169, 1.5172, 1.5176, 1.5180, 1.5184,
      1.5188, 1.5192, 1.5195, 1.5199, 1.5203, 1.5207, 1.5211, 1.5215,
      1.5219, 1.5224, 1.5228, 1.5232, 1.5236, 1.5240, 1.5244, 1.5249,
      1.5253, 1.5257, 1.5262, 1.5266, 1.5271, 1.5275, 1.5279, 1.5284,
      1.5289, 1.5293, 1.5298, 1.5302, 1.5307, 1.5312, 1.5317, 1.5321,
      1.5326, 1.5331, 1.5336, 1.5341, 1.5346, 1.5351, 1.5356, 1.5361,
      1.5366, 1.5371, 1.5376, 1.5382, 1.5387, 1.5392, 1.5398, 1.5403,
      1.5408, 1.5414, 1.5419, 1.5425, 1.5430, 1.5436, 1.5442, 1.5448,
      1.5453, 1.5459, 1.5465, 1.5471, 1.5477, 1.5483, 1.5489, 1.5495,
      1.5501, 1.5507, 1.5514, 1.5520, 1.5526, 1.5533, 1.5539, 1.5546,
      1.5552, 1.5559, 1.5566, 1.5572, 1.5579, 1.5586, 1.5593, 1.5600,
      1.5607, 1.5614, 1.5621, 1.5628, 1.5636, 1.5643, 1.5650, 1.5658,
      1.5665, 1.5673, 1.5680, 1.5688, 1.5696, 1.5704, 1.5712, 1.5720,
      1.5728, 1.5736, 1.5744, 1.5752, 1.5761, 1.5769, 1.5778, 1.5786,
      1.5795, 1.5804, 1.5813, 1.5822, 1.5831, 1.5840, 1.5849, 1.5858,
      1.5868, 1.5877, 1.5887, 1.5896, 1.5906};

  G4int n_entries = photon_energy.size();
  quartz_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);

  // for absorption_length
  photon_energy = {
      0.9999 * eV,
      1.1922 * eV,
      1.3776 * eV,
      1.5694 * eV,
      1.7712 * eV,
      1.9680 * eV,
      2.1377 * eV,
      2.3393 * eV,
      2.5303 * eV,
      2.6953 * eV,
      2.9520 * eV,
      3.0996 * eV,
      3.2627 * eV,
      3.6466 * eV,
      4.1328 * eV,
      4.5920 * eV,
      5.1660 * eV,
      5.6356 * eV,
      5.9040 * eV,
      6.1992 * eV,
      6.8880 * eV,
  };

  std::vector<G4double> absorption_length{
      100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m,
      100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m,
      100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 100.0 * m, 9.4029 * m,
      2.3965 * m, 0.5321 * m, 0.0158 * m};

  n_entries = photon_energy.size();
  quartz_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);

  m_material_map["QuartzKVC"]->SetMaterialPropertiesTable(quartz_prop);

  // +---------+
  // | Aerogel |
  // +---------+
  auto aerogel_prop = new G4MaterialPropertiesTable();
  photon_energy = {1.3 * eV, 7.0 * eV};
  refractive_index = {1.1, 1.1};
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
  m_material_map["Blacksheet"]->SetMaterialPropertiesTable(blacksheet_prop);

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
  refractive_index = {1.52, 1.52};          // 一般的な光電子増倍管の窓材（硝子）
  absorption_length = {1.0 * cm, 1.0 * cm}; // 可視光域で透過性あり

  n_entries = photon_energy.size();
  pmt_prop->AddProperty("RINDEX", &photon_energy[0], &refractive_index[0], n_entries);
  // pmt_prop->AddProperty("ABSLENGTH", &photon_energy[0], &absorption_length[0], n_entries);

  m_material_map["Glass"]->SetMaterialPropertiesTable(pmt_prop);

  // +----------------------------+
  // | PMT surface (POM) Property |
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
  const G4ThreeVector sac_size = G4ThreeVector(113.8 * mm, 145.0 * mm, 33.0 * mm); // SAC aerogel size
  const G4double teflon_thickness = 3.0 * mm;
  const G4double frame_thickness = 10.0 * mm;
  const G4double pmt_spacing = 25.0 * mm;
  const G4double pmt_radius = 25.8 * mm;   // PMT window radius
  const G4double pmt_thickness = 1.0 * mm; // PMT window thickness

  const G4ThreeVector origin(0, 0, 0);

  // ----------------------
  // 母ボリューム（Air）
  // ----------------------
  auto mother_solid = new G4Box("SACMotherSolid",
                                sac_size.x() / 2 + 50 * mm,
                                sac_size.y() / 2 + 50 * mm,
                                sac_size.z() / 2 + 50 * mm);
  auto mother_lv = new G4LogicalVolume(mother_solid,
                                       m_material_map["Air"],
                                       "SACMotherLV");
  new G4PVPlacement(nullptr, origin, mother_lv, "SACMotherPV", m_world_lv, false, 0, m_check_overlaps);
  mother_lv->SetVisAttributes(G4VisAttributes::GetInvisible());

  // ----------------------
  // Aerogel
  // ----------------------
  auto sac_solid = new G4Box("SACSolid", sac_size.x() / 2, sac_size.y() / 2, sac_size.z() / 2);
  auto sac_lv = new G4LogicalVolume(sac_solid, m_material_map["Aerogel"], "SACLV");
  new G4PVPlacement(nullptr, origin, sac_lv, "SACPV", mother_lv, false, 0, m_check_overlaps);
  sac_lv->SetVisAttributes(G4Colour::Cyan());

  // ----------------------
  // 上下 Teflon シート
  // ----------------------
  auto sheet_solid = new G4Box("TeflonSheetSolid", sac_size.x() / 2, teflon_thickness / 2, sac_size.z() / 2);
  auto sheet_lv = new G4LogicalVolume(sheet_solid, m_material_map["Teflon"], "TeflonSheetLV");

  G4ThreeVector sheet_pos_top(0, sac_size.y() / 2 + teflon_thickness / 2, 0);
  G4ThreeVector sheet_pos_bot(0, -sac_size.y() / 2 - teflon_thickness / 2, 0);
  new G4PVPlacement(nullptr, sheet_pos_top, sheet_lv, "TeflonSheetTopPV", mother_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(nullptr, sheet_pos_bot, sheet_lv, "TeflonSheetBotPV", mother_lv, false, 1, m_check_overlaps);
  sheet_lv->SetVisAttributes(G4Colour::White());

  // ----------------------
  // 上下 BlackSheet 追加
  // ----------------------
  auto black_sheet_solid = new G4Box("BlackSheetSolid", sac_size.x() / 2, 0.1 * mm, sac_size.z() / 2);
  auto black_sheet_lv = new G4LogicalVolume(black_sheet_solid, m_material_map["BlackSheet"], "BlackSheetLV");

  G4ThreeVector black_top_pos(0, sac_size.y() / 2 + teflon_thickness + 0.1 * mm, 0);
  G4ThreeVector black_bot_pos(0, -sac_size.y() / 2 - teflon_thickness - 0.1 * mm, 0);
  new G4PVPlacement(nullptr, black_top_pos, black_sheet_lv, "BlackSheetTopPV", mother_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(nullptr, black_bot_pos, black_sheet_lv, "BlackSheetBotPV", mother_lv, false, 1, m_check_overlaps);
  black_sheet_lv->SetVisAttributes(G4Colour::Black());

  // ----------------------
  // 側面 Teflon フレーム (穴あき)
  // ----------------------
  auto frame_outer = new G4Box("FrameOuter", frame_thickness / 2, sac_size.y() / 2, sac_size.z() / 2);
  G4SubtractionSolid *frame_with_holes = frame_outer;

  for (int i = 0; i < 3; ++i)
  {
    double y_pos = (i - 1) * pmt_spacing;
    auto hole = new G4Tubs("Hole", 0.0, pmt_radius + 0.5 * mm, sac_size.z() / 2 + 1.0 * mm, 0.0, 360.0 * deg);
    auto trans = G4Transform3D(G4RotationMatrix(), G4ThreeVector(0, y_pos, 0));
    frame_with_holes = new G4SubtractionSolid("FrameWithHoles", frame_with_holes, hole, trans);
  }

  auto frame_lv = new G4LogicalVolume(frame_with_holes, m_material_map["Teflon"], "TeflonFrameLV");

  G4ThreeVector frame_pos_left(-sac_size.x() / 2 - frame_thickness / 2, 0, 0);
  G4ThreeVector frame_pos_right(+sac_size.x() / 2 + frame_thickness / 2, 0, 0);
  new G4PVPlacement(nullptr, frame_pos_left, frame_lv, "TeflonFrameLeftPV", mother_lv, false, 0, m_check_overlaps);
  new G4PVPlacement(nullptr, frame_pos_right, frame_lv, "TeflonFrameRightPV", mother_lv, false, 1, m_check_overlaps);

  // ----------------------
  // PMT配置（14チャンネル）
  // ----------------------
  auto pmt_lv = SetupPMT();
  G4RotationMatrix *rotX = new G4RotationMatrix();
  rotX->rotateX(90. * deg);
  G4RotationMatrix *rotY = new G4RotationMatrix();
  rotY->rotateY(90. * deg);

  // 短辺（左右）3chずつ
  for (int i = 0; i < 3; ++i)
  {
    G4double dy = (i - 1) * pmt_spacing;
    G4ThreeVector left_pos(-sac_size.x() / 2 - frame_thickness - pmt_thickness / 2, dy, 0);
    G4ThreeVector right_pos(+sac_size.x() / 2 + frame_thickness + pmt_thickness / 2, dy, 0);
    new G4PVPlacement(rotY, left_pos, pmt_lv, "PMTPV", mother_lv, false, i, m_check_overlaps);
    new G4PVPlacement(rotY, right_pos, pmt_lv, "PMTPV", mother_lv, false, i + 3, m_check_overlaps);
  }

  // 長辺（上下）4chずつ
  for (int i = 0; i < 4; ++i)
  {
    G4double dx = (i - 1.5) * pmt_spacing;
    G4ThreeVector top_pos(dx, sac_size.y() / 2 + teflon_thickness + pmt_thickness / 2, 0);
    G4ThreeVector bot_pos(dx, -sac_size.y() / 2 - teflon_thickness - pmt_thickness / 2, 0);
    new G4PVPlacement(nullptr, top_pos, pmt_lv, "PMTPV", mother_lv, false, i + 6, m_check_overlaps);
    new G4PVPlacement(nullptr, bot_pos, pmt_lv, "PMTPV", mother_lv, false, i + 10, m_check_overlaps);
  }
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
