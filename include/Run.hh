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
/// \file Run.hh
/// \brief Definition of the Run class
//
// $Id: Run.hh 71375 2013-06-14 07:39:33Z maire $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"
#include <map>
#include "g4root.hh"

class DetectorConstruction;
class G4ParticleDefinition;

static struct Electron {
     Electron(): fX(0), fY(0.), fZ(0.), fT(0.) {}
     Electron(G4double x, G4double y, G4double z, G4double t): fX(x), fY(y), fZ(z), fT(t) {}
     void SetValues(G4double x, G4double y, G4double z, G4double t) {
       fX = x;
       fY = y;
       fZ = z;
       fT = t;	
     }    
     G4double  fX;
     G4double  fY;
     G4double  fZ;
     G4double  fT;
};  

static struct RecoClusterElectron {
  	RecoClusterElectron(Electron* recoElectron) {
      fIsClustered = 0;
      fIsAllClustered = 0;   
      fRecoElectron = recoElectron;    
      fClusterElectronList = new std::vector<RecoClusterElectron*>;    
    }
    G4bool fIsClustered;
    G4bool fIsAllClustered;    
    Electron* fRecoElectron;   
    std::vector<RecoClusterElectron*>* fClusterElectronList;
};

static struct RecoCluster {
     RecoCluster() {}
     std::vector<Electron> fElectronList;
     G4int fNElectrons;
     G4double fX;
     G4double fY;
     G4double fZ;
     void AddElectron(Electron electron);
     G4int GetNElectron() {return fElectronList.size();}
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Run : public G4Run
{
  public:
    Run(DetectorConstruction*);
   ~Run();

  public:
    void CountProcesses(const G4VProcess* process);                  
    void ParticleCount(G4String, G4double);
    void SumTrackLength (G4int,G4int,G4double,G4double,G4double,G4double);
    G4int GetARCount();
    void AddARCount();
    void AddCollisionsMod1();
    G4int GetCollisionsMod1();
    void AddCollisionsF1();
    G4int GetCollisionsF1();
    void SetPrimary(G4ParticleDefinition* particle, G4double energy);    
    void EndOfRun();
            
    virtual void Merge(const G4Run*);
   
  private:
    struct ParticleData {
     ParticleData()
       : fCount(0), fEmean(0.), fEmin(0.), fEmax(0.) {}
     ParticleData(G4int count, G4double ekin, G4double emin, G4double emax)
       : fCount(count), fEmean(ekin), fEmin(emin), fEmax(emax) {}
     G4int     fCount;
     G4double  fEmean;
     G4double  fEmin;
     G4double  fEmax;
    };
     
  private:
    DetectorConstruction* fDetector;
    G4ParticleDefinition* fParticle;
    G4double              fEkin;
        
    std::map<G4String,G4int>        fProcCounter;            
    std::map<G4String,ParticleData> fParticleDataMap;
        
    G4int    fNbStep1, fNbStep2;
    G4double fTrackLen1, fTrackLen2;
    G4double fTime1, fTime2;
    G4int fNbColMod1;
    G4int fARCount; // potential anti-resonance neutron candidates
    G4int fNbColF1;
        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

