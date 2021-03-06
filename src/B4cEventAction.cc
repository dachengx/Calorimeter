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
/// \file B4cEventAction.cc
/// \brief Implementation of the B4cEventAction class

#include "B4cEventAction.hh"
#include "B4cCalorimeterSD.hh"
#include "B4cCalorHit.hh"
#include "B4Analysis.hh"
#include "B4cDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::B4cEventAction()
 : G4UserEventAction(),
   fGapHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cEventAction::~B4cEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4cCalorHitsCollection* 
B4cEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<B4cCalorHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("B4cEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::BeginOfEventAction(const G4Event* event)
{
  auto analysisManager = G4AnalysisManager::Instance();
  G4PrimaryVertex* gpVertex = event->GetPrimaryVertex();
	G4PrimaryParticle* gpParticle = gpVertex->GetPrimary();
	G4ThreeVector pd = gpParticle->GetMomentumDirection();
  auto detectorThickness = B4cDetectorConstruction::fNofLayers * (B4cDetectorConstruction::absoThickness + B4cDetectorConstruction::gapThickness) - B4cDetectorConstruction::absoThickness;
  auto x = (-(gpVertex->GetZ0()) - detectorThickness / 2) * pd[0] / pd[2];
  auto y = (-(gpVertex->GetZ0()) - detectorThickness / 2) * pd[1] / pd[2];
  analysisManager->FillNtupleFColumn(2, 0, gpParticle->GetTotalEnergy()/MeV);
  analysisManager->FillNtupleFColumn(2, 1, pd[0]);
  analysisManager->FillNtupleFColumn(2, 2, pd[1]);
  analysisManager->FillNtupleFColumn(2, 3, pd[2]);
  analysisManager->FillNtupleFColumn(2, 4, x/mm);
  analysisManager->FillNtupleFColumn(2, 5, y/mm);
  analysisManager->AddNtupleRow(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4cEventAction::EndOfEventAction(const G4Event* event)
{  
  // Get hits collections IDs (only once)
  if ( fGapHCID == -1 ) {
    fGapHCID 
      = G4SDManager::GetSDMpointer()->GetCollectionID("GapHitsCollection");
  }

  // Get hits collections
  auto gapHC = GetHitsCollection(fGapHCID, event);

  // Print per event (modulo n)
  //
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;
  }  
  
  // Fill histograms, ntuple
  //

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
  // fill ntuple
  G4int nCells = gapHC->entries()-1;
  if (nCells != B4cDetectorConstruction::fNofLayers * B4cDetectorConstruction::fNofCells * B4cDetectorConstruction::fNofCells) {
    std::ostringstream message;
    message << "Sanity check: wrong solid extent." << G4endl
            << "        Replicated geometry, logical volume: ";
    G4Exception("G4SmartVoxelHeader::BuildReplicaVoxels", "GeomMgt0002", FatalException, message);
  }
  
  // G4cout << nCells << G4endl;
  for (G4int i=0; i<nCells; i++ ) {
    auto gapHiti = (*gapHC)[i];
    // if (gapHiti->GetEdep() > 0) { G4cout << i << "-" << gapHiti->GetEdep() << G4endl; }
    analysisManager->FillNtupleFColumn(0, i, gapHiti->GetEdep()/MeV);
    // analysisManager->FillNtupleDColumn(i, eventID);
  }

  analysisManager->AddNtupleRow(0);

  G4double xydeposit;
  for (G4int i=0; i<nCells; i++ ) {
    auto gapHiti = (*gapHC)[i];
    if ( i % B4cDetectorConstruction::fNofLayers == 0 ) {
      xydeposit = 0;
    }
    xydeposit += gapHiti->GetEdep();
    if ( (i + 1) % B4cDetectorConstruction::fNofLayers == 0 ) {
      analysisManager->FillNtupleFColumn(3, i / B4cDetectorConstruction::fNofLayers, xydeposit / MeV);
    }
  }
  analysisManager->AddNtupleRow(3);

}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
