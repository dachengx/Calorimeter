//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// 
/// \file B4cDetectorConstruction.cc
/// \brief Implementation of the B4cDetectorConstruction class

#include "B4cDetectorConstruction.hh"
#include "B4cCalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::B4cDetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cDetectorConstruction::~B4cDetectorConstruction()
{ 
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();

  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_Fe");
  nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");

  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4cDetectorConstruction::DefineVolumes()
{
  // Geometry parameters

  G4double layerThickness = absoThickness + gapThickness;
  G4double calorThickness = fNofLayers * layerThickness - absoThickness;
  // G4double totalThickness = calorThickness + 2 * (topThickness + surThickness);
  G4double worldSizeXY = 2 * calorSizeXY * fNofCells;
  G4double worldSizeZ  = 12.*m;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("G4_POLYSTYRENE");
  auto topMaterial = G4Material::GetMaterial("G4_POLYSTYRENE");
  auto surMaterial = G4Material::GetMaterial("G4_Fe");

  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial || ! topMaterial || ! surMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  

  //     
  // World
  //
  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //                               
  // Calorimeter
  //  
  auto xycalorimeterS
    = new G4Box("xyCalorimeter",     // its name
                 calorSizeXY/2 * fNofCells, calorSizeXY/2 * fNofCells, calorThickness/2); // its size

  auto xycalorLV
    = new G4LogicalVolume(
                 xycalorimeterS,     // its solid
                 absorberMaterial,  // its material
                 "xyCalorimeterLV");   // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 xycalorLV,          // its logical volume                         
                 "xyCalorimeter",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  //
  // Layer
  //
  auto xylayerS
    = new G4Box("xyLayer",           // its name
                 calorSizeXY/2, calorSizeXY/2 * fNofCells, calorThickness/2); //its size

  auto xylayerLV
    = new G4LogicalVolume(
                 xylayerS,           // its solid
                 absorberMaterial,  // its material
                 "xyLayerLV");         // its name

  new G4PVReplica(
                 "xyLayer",      // its name
                 xylayerLV,      // its logical volume
                 xycalorLV,      // its mother
                 kXAxis,           // axis of replication
                 fNofCells,        // number of replica
                 calorSizeXY);     // witdth of replica

  //
  // Layer
  //
  auto ylayerS
    = new G4Box("xLayer",           // its name
                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); //its size

  auto ylayerLV
    = new G4LogicalVolume(
                 ylayerS,           // its solid
                 absorberMaterial,  // its material
                 "xLayerLV");         // its name

  new G4PVReplica(
                 "xLayer",      // its name
                 ylayerLV,      // its logical volume
                 xylayerLV,      // its mother
                 kYAxis,           // axis of replication
                 fNofCells,        // number of replica
                 calorSizeXY);     // witdth of replica

  //
  // Layer
  //
  auto layerSgap
    = new G4Box("gapLayer",           // its name
                 calorSizeXY/2, calorSizeXY/2, gapThickness/2); //its size

  auto layerLVgap
    = new G4LogicalVolume(
                 layerSgap,      // its solid
                 absorberMaterial,    // its material
                 "gapLayerLV");    // its name

  new G4PVReplica(
                 "gapLayer",       // its name
                 layerLVgap,       // its logical volume
                 ylayerLV,       // its mother
                 kZAxis,           // axis of replication
                 fNofLayers,   // number of replica
                 layerThickness);  // witdth of replica

  //                               
  // Gap
  //
  auto gapS 
    = new G4Box("Gap",             // its name
                 calorSizeXY/2, calorSizeXY/2, gapThickness/2); // its size

  auto gapLV
    = new G4LogicalVolume(
                 gapS,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // its position
                 gapLV,            // its logical volume                         
                 "Gap",            // its name
                 layerLVgap,       // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //                               
  // Plastic
  //
  auto plasticSl 
    = new G4Box("Plasticl",             // its name
                 calorSizeXY/2 * fNofCells, calorSizeXY/2 * fNofCells, topThickness/2); // its size

  auto plasticLVl
    = new G4LogicalVolume(
                 plasticSl,             // its solid
                 topMaterial,      // its material
                 "PlasticLVl");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0, 0, calorThickness/2 + topThickness/2),  // its position
                 plasticLVl,            // its logical volume                         
                 "Plasticl",            // its name
                 worldLV,       // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //                               
  // Plastic
  //
  auto plasticSr
    = new G4Box("Plasticr",             // its name
                 calorSizeXY/2 * fNofCells, calorSizeXY/2 * fNofCells, topThickness/2); // its size

  auto plasticLVr
    = new G4LogicalVolume(
                 plasticSr,             // its solid
                 topMaterial,      // its material
                 "PlasticLVr");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0, 0, -(calorThickness/2 + topThickness/2)),  // its position
                 plasticLVr,            // its logical volume                         
                 "Plasticr",            // its name
                 worldLV,       // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //                               
  // Fe
  //
  auto feSl 
    = new G4Box("Fel",             // its name
                 calorSizeXY/2 * fNofCells, calorSizeXY/2 * fNofCells, surThickness/2); // its size

  auto feLVl
    = new G4LogicalVolume(
                 feSl,             // its solid
                 surMaterial,      // its material
                 "FeLVl");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0, 0, calorThickness/2 + topThickness + surThickness/2),  // its position
                 feLVl,            // its logical volume                         
                 "Fel",            // its name
                 worldLV,       // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //                               
  // Fe
  //
  auto feSr
    = new G4Box("Fer",             // its name
                 calorSizeXY/2 * fNofCells, calorSizeXY/2 * fNofCells, surThickness/2); // its size

  auto feLVr
    = new G4LogicalVolume(
                 feSr,             // its solid
                 surMaterial,      // its material
                 "FeLVr");         // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0, 0, -(calorThickness/2 + topThickness + surThickness/2)),  // its position
                 feLVr,            // its logical volume                         
                 "Fer",            // its name
                 worldLV,       // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName() 
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;

  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  xycalorLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // 
  // Sensitive detectors
  //

  auto gapSD
    = new B4cCalorimeterSD("GapSD", "GapHitsCollection", fNofLayers * fNofCells * fNofCells);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
