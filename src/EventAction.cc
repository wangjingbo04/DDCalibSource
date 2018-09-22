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
  
  fClusterRadius = 20.0*cm;    // clustering window (cm) //Ioana... initia 300.0
  fMinClusterElectrons = 50;    // minimum clustered digits
  fMinElectronsPerCluster = 0;
  fRecoElectronList = new std::vector<RecoElectron*>;
  fClusterList = new std::vector<RecoCluster*>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
//----------------------------------------------------------------- 
   G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

   // run clustering algorithm
   // ========================
   this->RecoClusters(fRecoElectronList);
   	
   //count electrons created in this event
   //G4cout<< "\n Number of electron tracks in this event = "<<fElectronList.size()<< G4endl;
   analysisManager->FillNtupleIColumn(3, 0, fRecoElectronList->size());
   analysisManager->FillNtupleIColumn(3, 1, fClusterList->size());
   int NCluster = 0;
   for(G4int i=0;i<fClusterList->size();i++) {
     if(fClusterList->at(i)->GetNElectrons()>420 ) NCluster;
   }
   analysisManager->FillNtupleIColumn(3, 2, fClusterList->size());
   analysisManager->AddNtupleRow(3);
   
   for(G4int i=0;i<fClusterList->size();i++) {
   	 analysisManager->FillNtupleIColumn(4, 0, fClusterList->at(i)->GetNElectrons());
     analysisManager->AddNtupleRow(4);
   }
   for(G4int i=0;i<fClusterList->size();i++) {
     if(fClusterList->at(i)->GetNElectrons()>420 ) {
     	 analysisManager->FillNtupleIColumn(4, 1, fClusterList->at(i)->GetNElectrons());
       analysisManager->AddNtupleRow(4);
     }
   }
  if(fRecoElectronList){
    fRecoElectronList->clear();
    delete fRecoElectronList; fRecoElectronList = 0;
  }
  if(fClusterList){
    fClusterList->clear();
    delete fClusterList; fClusterList = 0;
  }
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

std::vector<RecoElectron*>* EventAction::ResetElectrons(std::vector<RecoElectron*>* myElectronList)
{
  for( G4int ielectron=0; ielectron<myElectronList->size(); ielectron++ ){
    RecoElectron* recoElectron = (RecoElectron*)(myElectronList->at(ielectron));
  }

  return myElectronList;
}

std::vector<RecoCluster*>* EventAction::RecoClusters(std::vector<RecoElectron*>* myElectronList)
{  

  // delete cluster electrons
  // =====================
  for( G4int i=0; i<vClusterElectronList.size(); i++ ){
    delete (RecoClusterElectron*)(vClusterElectronList.at(i));
  }
  vClusterElectronList.clear();

  // delete clusters
  // ===============
  for( G4int i=0; i<vClusterList.size(); i++ ){
    delete (RecoCluster*)(vClusterList.at(i));
  }
  vClusterList.clear();

  // clear vector clusters
  // =====================
  fClusterList->clear();
  
  //histClusters->Fill(1.0);

  // make cluster electrons
  // ===================
  for( G4int ielectron=0; ielectron<myElectronList->size(); ielectron++ ){
    RecoElectron* recoElectron = (RecoElectron*)(myElectronList->at(ielectron));
    RecoClusterElectron* clusterElectron = new RecoClusterElectron(recoElectron);
    vClusterElectronList.push_back(clusterElectron);
  }
  
  //histClusters->Fill(1.0);

  // run clustering algorithm
  // ========================
  for( G4int ielectron1=0; ielectron1<vClusterElectronList.size(); ielectron1++ ){
    for( G4int ielectron2=ielectron1+1; ielectron2<vClusterElectronList.size(); ielectron2++ ){

      RecoClusterElectron* felectron1 = (RecoClusterElectron*)(vClusterElectronList.at(ielectron1));
      RecoClusterElectron* felectron2 = (RecoClusterElectron*)(vClusterElectronList.at(ielectron2));

      G4double dx = felectron1->GetX() - felectron2->GetX();
      G4double dy = felectron1->GetY() - felectron2->GetY();
      G4double dz = felectron1->GetZ() - felectron2->GetZ();
      G4double dt = felectron1->GetTime() - felectron2->GetTime();
      G4double drsq = dx*dx + dy*dy + dz*dz;

      if( drsq>0.0
       && drsq<fClusterRadius*fClusterRadius
       /*&& fabs(dt)<fTimeWindowC*/ ){
        felectron1->AddClusterElectron(felectron2);
        felectron2->AddClusterElectron(felectron1);
      }
    }
  }
  
  //histClusters->Fill(1.0);

  // collect up clusters
  // ===================
  G4bool carryon = 0;

  for( G4int ielectron=0; ielectron<vClusterElectronList.size(); ielectron++ ){
    RecoClusterElectron* felectron = (RecoClusterElectron*)(vClusterElectronList.at(ielectron));

    if( felectron->IsClustered()==0
     && felectron->GetNClusterElectrons()>0 ){
        
      vClusterElectronCollection.clear();
      vClusterElectronCollection.push_back(felectron);
      felectron->SetClustered();

      carryon = 1;
      while( carryon ){
        carryon = 0;
        for( G4int jelectron=0; jelectron<vClusterElectronCollection.size(); jelectron++ ){
          RecoClusterElectron* celectron = (RecoClusterElectron*)(vClusterElectronCollection.at(jelectron));
	  
	        //std::cout << "celectron->GetNClusterElectrons() = " << celectron->GetNClusterElectrons() << std::endl;
	        G4double nElectrons = celectron->GetNClusterElectrons();
	        //vNelectronsCluster.push_back(nElectrons);
	        
	        if (celectron->GetNClusterElectrons() > fMinElectronsPerCluster) {
                if( celectron->IsAllClustered()==0 ){
                  for( G4int kelectron=0; kelectron<celectron->GetNClusterElectrons(); kelectron++ ){
                    RecoClusterElectron* celectronnew = (RecoClusterElectron*)(celectron->GetClusterElectron(kelectron));
                     
                    if( celectronnew->IsClustered()==0 ){
                      vClusterElectronCollection.push_back(celectronnew);
                      celectronnew->SetClustered();
                      carryon = 1;
                    }
                  }
                }
	        
	        }//end if min of # of hits per cluster
        }
      } 
    
      if( (G4int)vClusterElectronCollection.size()>=fMinClusterElectrons ){
        RecoCluster* cluster = new RecoCluster();
        fClusterList->push_back(cluster);
        vClusterList.push_back(cluster);

        for( G4int jelectron=0; jelectron<vClusterElectronCollection.size(); jelectron++ ){
          RecoClusterElectron* celectron = (RecoClusterElectron*)(vClusterElectronCollection.at(jelectron));
          RecoElectron* recoelectron = (RecoElectron*)(celectron->GetRecoElectron());
          cluster->AddElectron(recoelectron);        
        }
      }
    }
  }

  // return vector of clusters
  // =========================
  return fClusterList;
}

//std::vector<RecoCluster*>* EventAction::RecoClusters(std::vector<Electron*>* myElectronList) {
//  return 0;
//}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


