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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* det)
: G4VUserPrimaryGeneratorAction(),fParticleGun(0), fDetector(det)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle); 
  
  // default particle kinematic

  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(2.5*MeV);
  G4double x = 0.*cm;
  G4double y = 0.*cm;
  G4double z = det->GetInsulatorHeight()/2 + det->GetThermalAbsorberHeight() 
               + det->GetFilterHeight() + det->GetModerator2Height() 
               + det->GetModerator1Height() - det->GetDDGeneratorHeight()/2;
  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
  fMessenger = new G4GenericMessenger(this,"/primary/", "...doc...");
  fMessenger->DeclareMethod("updateGunPosition", &PrimaryGeneratorAction::UpdateGunPosition, "...doc...");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of event
  //
  //distribution uniform in solid angle
  //
  G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double ux = sinTheta*std::cos(phi),
           uy = sinTheta*std::sin(phi),
           uz = std::abs(cosTheta);

//  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(ux,uy,uz));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));  
  
  //set energy
//  G4double eRand = G4UniformRand();
//  fParticleGun->SetParticleEnergy(eRand*MeV);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::UpdateGunPosition() {
  G4double x = 0.*cm;
  G4double y = 0.*cm;
  G4double z = fDetector->GetInsulatorHeight()/2 + fDetector->GetThermalAbsorberHeight() 
               + fDetector->GetFilterHeight() + fDetector->GetModerator2Height() 
               + fDetector->GetModerator1Height() - fDetector->GetDDGeneratorHeight()/2;
  fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));	
}

