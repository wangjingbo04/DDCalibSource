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
/// \file EventAction.hh
/// \brief Definition of the EventAction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class EventAction : public G4UserEventAction
{
  public:
    EventAction();
   ~EventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
    void AddNeutronEnter_Moderator1() {fNNeutronEnter_Moderator1++;};
    void AddNeutronExit_Moderator1() {fNNeutronExit_Moderator1++;};
    
    void AddNeutronEnter_Moderator2() {fNNeutronEnter_Moderator2++;};
    void AddNeutronExit_Moderator2() {fNNeutronExit_Moderator2++;};
    
    void AddNeutronEnter_Filter() {fNNeutronEnter_Filter++;};
    void AddNeutronExit_Filter() {fNNeutronEnter_Filter++;};
    
    void AddNeutronEnter_ThermalAbsorber() {fNNeutronEnter_ThermalAbsorber++;};
    void AddNeutronExit_ThermalAbsorber() {fNNeutronEnter_ThermalAbsorber++;};
    
    void AddNeutronEnter_Port() {fNNeutronEnter_Port++;};
    void AddNeutronExit_Port() {fNNeutronEnter_Port++;};
    
    void AddNeutronEnter_ArBuffer() {fNNeutronEnter_ArBuffer++;};
    void AddNeutronExit_ArBuffer() {fNNeutronEnter_ArBuffer++;};
    
    void AddNeutronEnter_LArPool() {fNNeutronEnter_LArPool++;};
    void AddNeutronExit_LArPool() {fNNeutronEnter_LArPool++;};
    
    G4int GetNNeutronEnter_Moderator1() {return fNNeutronEnter_Moderator1;};
    G4int GetNNeutronExit_Moderator1() {return fNNeutronExit_Moderator1;};

    G4int GetNNeutronEnter_Moderator2() {return fNNeutronEnter_Moderator2;};
    G4int GetNNeutronExit_Moderator2() {return fNNeutronExit_Moderator2;};

    G4int GetNNeutronEnter_Filter() {return fNNeutronEnter_Filter;};
    G4int GetNNeutronExit_Filter() {return fNNeutronEnter_Filter;};

    G4int GetNNeutronEnter_ThermalAbsorber() {return fNNeutronEnter_ThermalAbsorber;};
    G4int GetNNeutronExit_ThermalAbsorber() {return fNNeutronEnter_ThermalAbsorber;};
    
    G4int GetNNeutronEnter_Port() {return fNNeutronEnter_Port;};
    G4int GetNNeutronExit_Port() {return fNNeutronEnter_Port;};
    
    G4int GetNNeutronEnter_ArBuffer() {return fNNeutronEnter_ArBuffer;};
    G4int GetNNeutronExit_ArBuffer() {return fNNeutronEnter_ArBuffer;};
    
    G4int GetNNeutronEnter_LArPool() {return fNNeutronEnter_LArPool;};
    G4int GetNNeutronExit_LArPool() {return fNNeutronEnter_LArPool;};
                
  private:                  
    G4double				fNNeutronEnter_Moderator1;
    G4double				fNNeutronExit_Moderator1;
    
    G4double				fNNeutronEnter_Moderator2;
    G4double				fNNeutronExit_Moderator2;
    
    G4double				fNNeutronEnter_Filter;
    G4double				fNNeutronExit_Filter;
    
    G4double				fNNeutronEnter_ThermalAbsorber;
    G4double				fNNeutronExit_ThermalAbsorber;  
    
    G4double				fNNeutronEnter_Port;
    G4double				fNNeutronExit_Port;  
    
    G4double				fNNeutronEnter_ArBuffer;
    G4double				fNNeutronExit_ArBuffer;  
    
    G4double				fNNeutronEnter_LArPool;
    G4double				fNNeutronExit_LArPool;  

    
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
