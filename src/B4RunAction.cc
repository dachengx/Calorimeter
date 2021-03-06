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
/// \file B4RunAction.cc
/// \brief Implementation of the B4RunAction class

#include "B4RunAction.hh"
#include "B4Analysis.hh"
#include "B4cDetectorConstruction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::B4RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  analysisManager->SetFileName("B5");
    // Note: merging ntuples is available only with Root output

  // Creating ntuple
  // G4int ntupleID = 
  // G4int eventID = analysisManager->CreateNtupleDColumn("Nevent");
  G4int nCells = B4cDetectorConstruction::fNofLayers * B4cDetectorConstruction::fNofCells * B4cDetectorConstruction::fNofCells;
  analysisManager->CreateNtuple("Energy3D", "Edep");
  for (G4int i=0; i<nCells; i++ ) {
    analysisManager->CreateNtupleFColumn(0, "e" + std::to_string(i));
  }
  analysisManager->FinishNtuple(0);

  analysisManager->CreateNtuple("XYZ", "Position");
  analysisManager->CreateNtupleFColumn(1, "x");
  analysisManager->CreateNtupleFColumn(1, "y");
  analysisManager->CreateNtupleFColumn(1, "z");
  analysisManager->FinishNtuple(1);

  analysisManager->CreateNtuple("Truth", "Truth");
  analysisManager->CreateNtupleFColumn(2, "energy");
  analysisManager->CreateNtupleFColumn(2, "px");
  analysisManager->CreateNtupleFColumn(2, "py");
  analysisManager->CreateNtupleFColumn(2, "pz");
  analysisManager->CreateNtupleFColumn(2, "x");
  analysisManager->CreateNtupleFColumn(2, "y");
  analysisManager->FinishNtuple(2);

  analysisManager->CreateNtuple("Energy2D", "Edep");
  for (G4int i=0; i<nCells; i++ ) {
    if ( i % B4cDetectorConstruction::fNofLayers == 0) {
      analysisManager->CreateNtupleFColumn(3, "e" + std::to_string(i));
    }
  }
  analysisManager->FinishNtuple(3);

  analysisManager->CreateNtuple("XY", "Position");
  analysisManager->CreateNtupleFColumn(4, "x");
  analysisManager->CreateNtupleFColumn(4, "y");
  analysisManager->FinishNtuple(4);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4RunAction::~B4RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  analysisManager->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  auto nCells = B4cDetectorConstruction::fNofLayers * B4cDetectorConstruction::fNofCells * B4cDetectorConstruction::fNofCells;
  auto calorThickness = B4cDetectorConstruction::fNofLayers * (B4cDetectorConstruction::absoThickness + B4cDetectorConstruction::gapThickness) - B4cDetectorConstruction::absoThickness;
  auto xorigin = -B4cDetectorConstruction::calorSizeXY * B4cDetectorConstruction::fNofCells / 2 + B4cDetectorConstruction::calorSizeXY / 2;
  auto yorigin = -B4cDetectorConstruction::calorSizeXY * B4cDetectorConstruction::fNofCells / 2 + B4cDetectorConstruction::calorSizeXY / 2;
  auto zorigin = -calorThickness / 2 + B4cDetectorConstruction::gapThickness / 2;

  for (G4int i=0; i<nCells; i++ ) {
    int ii = i / (B4cDetectorConstruction::fNofLayers * B4cDetectorConstruction::fNofCells);
    int ij = i % (B4cDetectorConstruction::fNofLayers * B4cDetectorConstruction::fNofCells) / B4cDetectorConstruction::fNofLayers;
    int ik = i % B4cDetectorConstruction::fNofLayers;
    analysisManager->FillNtupleFColumn(1, 0, (ii * B4cDetectorConstruction::calorSizeXY + xorigin) / mm);
    analysisManager->FillNtupleFColumn(1, 1, (ij * B4cDetectorConstruction::calorSizeXY + yorigin) / mm);
    analysisManager->FillNtupleFColumn(1, 2, (ik * (B4cDetectorConstruction::absoThickness + B4cDetectorConstruction::gapThickness) + zorigin) / mm);
    analysisManager->AddNtupleRow(1);
  }

  auto nCellsxy = B4cDetectorConstruction::fNofCells * B4cDetectorConstruction::fNofCells;
  for (G4int i=0; i<nCellsxy; i++ ) {
    int ii = i / B4cDetectorConstruction::fNofCells;
    int ij = i % B4cDetectorConstruction::fNofCells;
    analysisManager->FillNtupleFColumn(4, 0, (ii * B4cDetectorConstruction::calorSizeXY + xorigin) / mm);
    analysisManager->FillNtupleFColumn(4, 1, (ij * B4cDetectorConstruction::calorSizeXY + yorigin) / mm);
    analysisManager->AddNtupleRow(4);
  }

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
