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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// $Id: EventAction.cc 76293 2013-11-08 13:11:23Z gcosmo $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction()
:G4UserEventAction()
{  
                            
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdep1 = fEdep2 = fWeight1 = fWeight2 = 0.;
  
  fNNeutronExit_Generator = 0;
  
  fNNeutronEnter_Moderator = 0;  
  fNNeutronExit_Moderator = 0;   
                              
  fNNeutronEnter_Filter1 = 0;  
  fNNeutronExit_Filter1 = 0;   
                              
  fNNeutronEnter_Filter2 = 0;      
  fNNeutronExit_Filter2 = 0; 

  fNNeutronEnter_Filter3 = 0;      
  fNNeutronExit_Filter3 = 0;      
                              
  fNNeutronEnter_ThermalAbsorber = 0;  
  fNNeutronExit_ThermalAbsorber = 0; 
  
  fNNeutronEnter_Port = 0;  
  fNNeutronExit_Port = 0; 
  
  fNNeutronEnter_Cryostat = 0;
  fNNeutronExit_Cryostat = 0;
  
  fNNeutronEnter_ArBuffer = 0;  
  fNNeutronExit_ArBuffer = 0; 
  
  fNNeutronEnter_LArPool = 0;  
  fNNeutronExit_LArPool = 0;

  fGammaEnter_World = 0;
  fGammaEnter_Shield = 0;
  fGammaEnter_Reflector = 0;
  fGammaEnter_LArPool = 0;
  fGammaExit_LArPool = 0;
 
  fNNeutronExit_Shield = 0;
  fNNeutronEnter_World = 0;
  
//  if(fElectronList!=NULL){
//    for(G4int i=0; i<fElectronList->size(); i++) {
//      delete fElectronList->at(i);  fElectronList->at(i) = 0;
//      fElectronList->clear();
//    }
//  }
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
//----------------------------------------------------------------- 
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
   	
   //count electrons created in this event
   //G4cout<< "\n Number of electron tracks in this event = "<<fElectronList.size()<< G4endl;
   analysisManager->FillNtupleIColumn(3, 0, fElectronList.size());
   analysisManager->AddNtupleRow(3);   
}

void EventAction::AddEdep(G4int iVol, G4double edep,
                                      G4double time, G4double weight) //add energy doposition from all the tracks within 1us
{
  // initialize t0
  if (fTime0 < 0.) fTime0 = time;
  
  // out of time window ?
  const G4double TimeWindow (1*microsecond);
  
  if (std::fabs(time - fTime0) > TimeWindow) return;
  if (iVol == 1) { fEdep1 += edep; fWeight1 += edep*weight;}
  else {
  	fEdep2 += edep; fWeight2 += edep*weight;
  }  
  	
}

//std::vector<RecoCluster*>* EventAction::RecoClusters(std::vector<Electron*>* myElectronList) {
//  return 0;
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


