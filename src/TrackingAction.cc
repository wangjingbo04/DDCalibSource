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
/// \file TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// $Id: TrackingAction.cc 69099 2013-04-18 12:25:19Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "TreeMaker.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TrackingAction::TrackingAction(DetectorConstruction* det, EventAction* EA)
:G4UserTrackingAction(),fEvent(EA),fDetector(det),
 fNbStep1(0),fNbStep2(0),fTrackLen1(0.),fTrackLen2(0.),fTime1(0.),fTime2(0.)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
	// original
  fNbStep1 = fNbStep2 = 0;
  fTrackLen1 = fTrackLen2 = 0.;
  fTime1 = fTime2 = 0.;
  
  // gamma cascade simulation
  Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  TreeMaker* treemaker = TreeMaker::Instance();
  //which volume ?
  G4LogicalVolume* lVolume = track->GetVolume()->GetLogicalVolume();
  G4int iVol = 0;
  if (lVolume == fDetector->GetLogicPool())   iVol = 1; //in LAr pool
  else iVol = 2; // elsewhere
  G4int trackID = track->GetTrackID();
  G4int parentID = track->GetParentID();
  G4ThreeVector vertex = track->GetPosition(); //get position
  G4double energy = track->GetKineticEnergy();
  G4double time   = track->GetGlobalTime();
  G4double weight = track->GetWeight();
  const G4ParticleDefinition* particle = track->GetParticleDefinition();  
  G4String name   = particle->GetParticleName();
  G4int pid       = particle->GetPDGEncoding();
  G4int Z         = particle->GetAtomicNumber();
  G4int A         = particle->GetAtomicMass();
  G4double charge = particle->GetPDGCharge(); 
  fCharge = particle->GetPDGCharge();
  fMass   = particle->GetPDGMass(); 
  G4bool condition = false;
  
  //Emission products for all processes
  analysisManager->FillH1(28, energy, weight);
  
  //primary particles
  if(parentID == 0) {
  	treemaker->AddPrimaryGamma(vertex.x()/cm, vertex.y()/cm, vertex.z()/cm, time/ns, energy/MeV);
//  	analysisManager->FillNtupleDColumn(0, 0, vertex.x());
//  	analysisManager->FillNtupleDColumn(0, 1, vertex.y());
//  	analysisManager->FillNtupleDColumn(0, 2, vertex.z());
//  	analysisManager->FillNtupleDColumn(0, 3, energy);
//  	analysisManager->AddNtupleRow(0);
  	
  	return; 
  }
  //secondary particles only
  //energy spectrum
  
  if (name == "e-" ) {
  	RecoElectron*  electron = new RecoElectron(vertex.x()/cm, vertex.y()/cm, vertex.z()/cm, time/ns, energy/eV, 0);
  	fEvent->AddElectron(electron);
//  	treemaker->AddElectron(vertex.x()/cm, vertex.y()/cm, vertex.z()/cm, time/ns, energy/keV, weight, 0, 0);
//  	analysisManager->FillH1(29, energy);
//  	analysisManager->FillNtupleDColumn(1, 0, double(pid));
//  	analysisManager->FillNtupleDColumn(1, 1, vertex.x()/cm);
//  	analysisManager->FillNtupleDColumn(1, 2, vertex.y()/cm);
//  	analysisManager->FillNtupleDColumn(1, 3, vertex.z()/cm);
//  	analysisManager->FillNtupleDColumn(1, 4, energy);
//  	analysisManager->AddNtupleRow(1);
  }
  if (name == "gamma") analysisManager->FillH1(30, energy);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::UpdateTrackInfo(G4double ekin,G4double trackl,
                                     G4double time)
{
  const G4double thermal = 1*eV;
  if (ekin > thermal) {
    fNbStep1++; fTrackLen1 = trackl; fTime1 = time;    
  } else {
    fNbStep2++; fTrackLen2 = trackl - fTrackLen1; fTime2 = time - fTime1;  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TrackingAction::PostUserTrackingAction(const G4Track* track)
{
 // keep only primary neutron
 //
 G4int trackID = track->GetTrackID();
 if (trackID > 1) return;
 
 Run* run 
    = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun()); 
 run->SumTrackLength(fNbStep1,fNbStep2,fTrackLen1,fTrackLen2,fTime1,fTime2);
 
 // histograms
 //
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
 analysisManager->FillH1(0,fNbStep1);
 analysisManager->FillH1(1,fTrackLen1);
 analysisManager->FillH1(2,fTime1); 
 analysisManager->FillH1(3,fNbStep2);
 analysisManager->FillH1(4,fTrackLen2);
 analysisManager->FillH1(5,fTime2);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

