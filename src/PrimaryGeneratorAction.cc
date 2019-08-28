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
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorAction.cc 67268 2013-02-13 11:38:40Z ihrivnac $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include <math.h>
#include "g4root.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction(),fParticleGun(0), fDetector(det), fUseUserDefinedEnergy(false)
{
  // make angular distribution 
  probDist = new G4double[360]; 
  MakeDDAngularDistribution();
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle); 
  
  // default particle kinematic
  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));
  G4double x = fDetector->GetSourceCenterX();
  G4double y = fDetector->GetSourceCenterY();
  G4double z = fDetector->GenerateGunPosition();
  //G4double z = fDetector->GetCryostatHeight()/2 - 1*cm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
  fParticleGun->SetParticleEnergy(57*keV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,0,-1)); 

  fMessenger = new G4GenericMessenger(this,"/primary/", "...doc...");
  fMessenger->DeclareMethod("updateGunPosition", &PrimaryGeneratorAction::UpdateGunPosition, "...doc...");
  fMessenger->DeclareMethod("UseUserDefinedEnergy", &PrimaryGeneratorAction::UseUserDefinedEnergy, "...doc...");
  fMessenger->DeclareMethod("UseUserDefinedDirection", &PrimaryGeneratorAction::UseUserDefinedDirection, "...doc...");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun; fParticleGun = 0;
  delete [] probDist; probDist = 0;
  delete fMessenger; fMessenger = 0;
  if(fUseUserDefinedEnergy) {
    fNeutronEnergyFile->Close();
    delete fNeutronEnergyFile; fNeutronEnergyFile = 0;
  }
}

void PrimaryGeneratorAction::UseUserDefinedEnergy(TString histname)
{
	std::cout<<"PrimaryGenerator::UseUserDefinedEnergy(): Read external neutron beam energy from a ROOT file"<<std::endl;
	// Open the neutron energy file
	fUseUserDefinedEnergy = true;
  fNeutronEnergyFile = new TFile("NeutronEnergy.root");
  fNeutronEnergy = (TH1D*)fNeutronEnergyFile->Get(histname);
}

void PrimaryGeneratorAction::UseUserDefinedDirection()
{
	std::cout<<"PrimaryGenerator::UseUserDefinedDirection(): Use an analytical angular distribution "<<std::endl;
	fUseUserDefinedDirection = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{  
   // define particle energy 
   
   //fParticleGun->SetParticleEnergy(73*keV);  
   //fParticleGun->SetParticleEnergy(100.*G4UniformRand()*keV);
   // generate particle direction     
   if(fUseUserDefinedDirection) {
     G4double theta = DDrandom()*degree;
     G4double phi = twopi*G4UniformRand();
     G4double cosTheta = -std::cos(theta);
     G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
     G4double dx = sinTheta*std::cos(phi),
              dy = sinTheta*std::sin(phi),
              dz = cosTheta;
   	 fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx,dy,dz));  
   }
   if(fUseUserDefinedEnergy) {
   	 fParticleGun->SetParticleEnergy(fNeutronEnergy->GetRandom()*1000*keV);
   }
   //UpdateGunPosition();
   fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void PrimaryGeneratorAction::UpdateGunPosition() {
  G4double x = fDetector->GetSourceCenterX();
  G4double y = fDetector->GetSourceCenterY();
  G4double z = fDetector->GenerateGunPosition();
  //G4double z = fDetector->GetCryostatHeight()/2 - 1*cm;
  //G4double z = 0*cm;
  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));	
}                         
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::MakeDDAngularDistribution(){

    G4double thetaMin = -180;
    G4double thetaMax = 180;
    G4double A = 0.237;
    G4double B = -0.051;
    G4double C = -0.130;
    for(G4int i=0; i<360; i++){
      G4double theta = thetaMin + ((thetaMax - thetaMin)/360)*i;
      if(theta < -90 || theta > 90){                                               
        probDist[i] = A + B + C;
      }
      else{
        G4double cosine = std::cos(theta *(pi/180) + (pi/2));
        probDist[i] = A + B * std::pow(cosine,2) + C * std::pow(cosine,4);
      }
    }
}

G4double PrimaryGeneratorAction::DDrandom(){
    G4RandGeneral *randDistribution = new G4RandGeneral(probDist, 360);
    G4double random = randDistribution->shoot();
    G4double angle = -180 + (360 * random);
    return angle;
}
