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
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 71404 2013-06-14 16:56:38Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "TrackingAction.hh"
#include "EventAction.hh"
#include "HistoManager.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include<cmath>
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, TrackingAction* TrAct, EventAction* event)
: G4UserSteppingAction(),fDetector(det),fTrackingAction(TrAct), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  
  //track informations
  G4StepPoint* prePoint = aStep->GetPreStepPoint();   
  G4StepPoint* postPoint = aStep->GetPostStepPoint();
  G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition(); 
  
  G4VPhysicalVolume* preStepPhysical = prePoint->GetPhysicalVolume();
  G4VPhysicalVolume* postStepPhysical = postPoint->GetPhysicalVolume();

  // The track does not exist
  if(preStepPhysical == 0 || postStepPhysical == 0) return;
 
  // Both steps are in the World
  if(preStepPhysical->GetCopyNo() == -1 && postStepPhysical->GetCopyNo() == -1) return;
  
  G4LogicalVolume* preVolume = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4LogicalVolume* endVolume = postPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4Track* track         = aStep->GetTrack();
  G4double kinEnergy     = track->GetKineticEnergy();
  
  // incident neutron
  //
  if (aStep->GetTrack()->GetTrackID() == 1) { 
    G4double ekin  = postPoint->GetKineticEnergy();
    G4double trackl = aStep->GetTrack()->GetTrackLength();
    G4double time   = aStep->GetTrack()->GetLocalTime();           
    fTrackingAction->UpdateTrackInfo(ekin,trackl,time);
    G4AnalysisManager::Instance()->FillH1(6,ekin); // fill the energy of each step to the histogram
  }
  

  // neutrons from generator to 1st moderator
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
//  	&& aStep->GetTrack()->GetTrackID() == 1
  	&& postPoint->GetStepStatus() == fGeomBoundary 
  	&& preVolume == fDetector->GetLogicDDGenerator() 
  	&& endVolume == fDetector->GetLogicModerator1()) {
    //G4AnalysisManager::Instance()->FillH1(7,kinEnergy);
    fEventAction->AddNeutronEnter_Moderator1();
    if(fEventAction->GetNNeutronEnter_Moderator1() == 1) G4AnalysisManager::Instance()->FillH1(7,kinEnergy);
  }



  //Number of Collisions in 1st moderator
  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
   && fEventAction->GetNNeutronEnter_Moderator2()== 1
   && preVolume == fDetector->GetLogicModerator1() 
   || endVolume == fDetector->GetLogicModerator1()) {

         Run* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
         run->AddCollisionsMod1();
  }

  
  // neutrons from 1st moderator to 2nd moderator
  if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
//  	&& aStep->GetTrack()->GetTrackID() == 1
  	&& postPoint->GetStepStatus() == fGeomBoundary 
  	&& preVolume == fDetector->GetLogicModerator1() 
  	&& endVolume == fDetector->GetLogicModerator2()) {  
    //G4AnalysisManager::Instance()->FillH1(8,kinEnergy);
    fEventAction->AddNeutronEnter_Moderator2();
    if(fEventAction->GetNNeutronEnter_Moderator2() == 1) G4AnalysisManager::Instance()->FillH1(8,kinEnergy);
    Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    if(kinEnergy >= 57*keV && kinEnergy < 1*MeV && fEventAction->GetNNeutronEnter_Moderator2() == 1) run->AddARCount();
  }
  
  // neutrons from 2nd moderator to filter
  if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& postPoint->GetStepStatus() == fGeomBoundary 
  	&& preVolume == fDetector->GetLogicModerator2() 
  	&& endVolume == fDetector->GetLogicFilter()) {  
  	//G4AnalysisManager::Instance()->FillH1(9,kinEnergy);
  	fEventAction->AddNeutronEnter_Filter();
    if(fEventAction->GetNNeutronEnter_Filter() == 1) G4AnalysisManager::Instance()->FillH1(9,kinEnergy);
  }
  
  // neutrons from filter to thermal absorber
  if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& postPoint->GetStepStatus() == fGeomBoundary 
  	&& preVolume == fDetector->GetLogicFilter() 
  	&& endVolume == fDetector->GetLogicThermalAbsorber()) {
    //G4AnalysisManager::Instance()->FillH1(10,kinEnergy);
    fEventAction->AddNeutronEnter_ThermalAbsorber();
    if(fEventAction->GetNNeutronEnter_ThermalAbsorber() == 1) G4AnalysisManager::Instance()->FillH1(10,kinEnergy);
  }
  
  // neutrons from thermal neutron absorber to feedthrough port
  if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& postPoint->GetStepStatus() == fGeomBoundary 
  	&& preVolume == fDetector->GetLogicThermalAbsorber()
  	&& endVolume == fDetector->GetLogicPort()) {  
        //G4AnalysisManager::Instance()->FillH1(11,kinEnergy);

        G4ThreeVector p = aStep->GetTrack()->GetMomentumDirection();
        if (p.x() != 0 && p.z() != 0)
            G4AnalysisManager::Instance()->FillH1(14, GetTheta(p.x(), p.z()));
        if (p.x() != 0 && p.y() != 0)
            G4AnalysisManager::Instance()->FillH1(15, GetPhi(p.x(), p.y()));

        fEventAction->AddNeutronEnter_Port();
        if(fEventAction->GetNNeutronEnter_Port() == 1) G4AnalysisManager::Instance()->FillH1(11,kinEnergy);
  }
  
  // neutrons from feedthrough port to gas argon buffer
  if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& postPoint->GetStepStatus() == fGeomBoundary 
  	&& preVolume == fDetector->GetLogicPort()
  	&& endVolume == fDetector->GetLogicBuffer()) {  
    //G4AnalysisManager::Instance()->FillH1(12,kinEnergy);
    fEventAction->AddNeutronEnter_ArBuffer();
    if(fEventAction->GetNNeutronEnter_ArBuffer() == 1) G4AnalysisManager::Instance()->FillH1(12,kinEnergy);
  }
  
  // neutrons entering liquid Argon pool
  if ( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& postPoint->GetStepStatus() == fGeomBoundary 
  	&& preVolume == fDetector->GetLogicBuffer()
  	&& endVolume == fDetector->GetLogicPool()) {  
    //G4AnalysisManager::Instance()->FillH1(13,kinEnergy);
    fEventAction->AddNeutronEnter_LArPool();
    if(fEventAction->GetNNeutronEnter_LArPool() == 1) G4AnalysisManager::Instance()->FillH1(13,kinEnergy);	
  }
  
  // count processes
  // 
  const G4VProcess* process   = postPoint->GetProcessDefinedStep();
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->CountProcesses(process);
  
  G4String procName = process->GetProcessName();  
  if( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& procName == "nCapture" 
  	&& endVolume == fDetector->GetLogicModerator1()) {
  	G4ThreeVector postPoint  = aStep->GetPostStepPoint()->GetPosition();
  	G4double x = postPoint.x(), y = postPoint.y(), z = postPoint.z();
  	G4AnalysisManager::Instance()->FillH2(0,x, y);
  	G4AnalysisManager::Instance()->FillH2(1,x, z);	
  }
  
  procName = process->GetProcessName();  
  if( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& procName == "nCapture" 
  	&& endVolume == fDetector->GetLogicModerator2()) {
  	G4ThreeVector postPoint  = aStep->GetPostStepPoint()->GetPosition();
  	G4double x = postPoint.x(), y = postPoint.y(), z = postPoint.z();
  	G4AnalysisManager::Instance()->FillH2(2,x, y);
  	G4AnalysisManager::Instance()->FillH2(3,x, z);	
  }
  
  procName = process->GetProcessName();  
  if( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& procName == "nCapture" 
  	&& endVolume == fDetector->GetLogicFilter()) {
  	G4ThreeVector postPoint  = aStep->GetPostStepPoint()->GetPosition();
  	G4double x = postPoint.x(), y = postPoint.y(), z = postPoint.z();
  	G4AnalysisManager::Instance()->FillH2(4,x, y);
  	G4AnalysisManager::Instance()->FillH2(5,x, z);	
  }
  
  procName = process->GetProcessName();  
  if( aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  	&& procName == "nCapture" 
  	&& endVolume == fDetector->GetLogicPool()) {
  	G4ThreeVector postPoint  = aStep->GetPostStepPoint()->GetPosition();
  	G4double x = postPoint.x(), y = postPoint.y(), z = postPoint.z();
  	G4AnalysisManager::Instance()->FillH2(6,x, y);
  	G4AnalysisManager::Instance()->FillH2(7,x, z);	
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



G4double SteppingAction::GetTheta(G4double dx, G4double dz){

    const G4double pi = 3.14159265;
    G4double theta = std::atan(dz/dx) * (180/pi);


    if(dx > 0 && dz > 0)
        theta = theta + 90;
    else if(dx < 0 && dz < 0)
        theta = theta - 90;
    else if(dx > 0 && dz < 0)
        theta = theta + 90;
    else if(dx < 0 && dz > 0)
        theta = theta - 90;

    if(theta > 180 || theta < -180)
        G4cout << "Error - SteppingAction::GetTheta()  " << "dx = " << dx << ", dz = " << dz << G4endl;

    return theta;
}







G4double SteppingAction::GetPhi(G4double dx, G4double dy){
    
    const G4double pi = 3.14159265;
    G4double phi = std::atan(dy/dx) * (180/pi);

    G4cout << "Phi1: " << phi << G4endl;

    if(dx > 0 && dy > 0)
        phi = phi + 90;
    else if(dx < 0 && dy < 0)
        phi = phi - 90;
    else if(dx > 0 && dy < 0)
        phi = phi + 90;
    else if(dx < 0 && dy > 0)
        phi = phi - 90;
        
   if(phi > 180 || phi < -180)
        G4cout << "Error - SteppingAction::GetPhi() " << "dx = " << dx << ", dy = " << dy << G4endl;
    return phi;
}





