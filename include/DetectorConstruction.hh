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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();

    G4Material* 
    MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);
         
    void SetSize     (G4double);              
    void SetWorldMaterial (G4String);
    void SetFilterHeight(G4double);
    void SetModerator1Height(G4double);
    void SetModerator2Height(G4double);
    
    void GetMaterialTable();
    void SetModerator1Material (G4String);
    void GetModerator1Material ();
    void SetModerator2Material (G4String);
    void GetModerator2Material ();
    void SetFilterMaterial (G4String);
    void GetFilterMaterial ();
    void SetAbsorberMaterial (G4String);
    void GetAbsorberMaterial ();
    

  public:
  
     const
     G4VPhysicalVolume* GetWorld()      {return fPWorld;};   
     G4VPhysicalVolume* GetPhysiPool() {return fPhysiPool;};  
     G4VPhysicalVolume* GetPhysiBuffer() {return fPhysiBuffer;}; 
     G4VPhysicalVolume* GetPhysiPort() {return fPhysiPort;}; 
     G4VPhysicalVolume* GetPhysif() {return fPhysiNShield;};  
     G4VPhysicalVolume* GetPhysiDDGenerator() {return fPhysiDDGenerator;};
     G4VPhysicalVolume* GetPhysiReflector() {return fPhysiReflector;};  
     G4VPhysicalVolume* GetPhysiModerator1() {return fPhysiModerator1;};   
     G4VPhysicalVolume* GetPhysiModerator2() {return fPhysiModerator2;};  
     G4VPhysicalVolume* GetPhysiFilter() {return fPhysiFilter;};  
      
     
     G4LogicalVolume* GetLogicWorld() {return fLWorld;};   
     G4LogicalVolume* GetLogicPool() {return fLogicPool;}; 
     G4LogicalVolume* GetLogicBuffer() {return fLogicBuffer;}; 
     G4LogicalVolume* GetLogicPort() {return fLogicPort;}; 
     G4LogicalVolume* GetLogicShield() {return fLogicNShield;}; 
     G4LogicalVolume* GetLogicReflector() {return fLogicReflector;};  
     G4LogicalVolume* GetLogicDDGenerator() {return fLogicDDGenerator;}; 
     G4LogicalVolume* GetLogicModerator1() {return fLogicModerator1;};   
     G4LogicalVolume* GetLogicModerator2() {return fLogicModerator2;}; 
     G4LogicalVolume* GetLogicThermalAbsorber() {return fLogicThermalAbsorber;}; 
     G4LogicalVolume* GetLogicFilter() {return fLogicFilter;};              
                    
     G4double           GetSize()       {return fWorldSize;};      
     G4Material*        GetMaterial()   {return fMaterial;}; 
     G4double           GetPoolHeight()       {return fPoolHeight;};
     G4double           GetFilterHeight()       {return fFilterHeight;};
     G4double           GetModerator1Height()       {return fModerator1_Height;};
     G4double           GetModerator2Height()       {return fModerator2_Height;};
     G4double           GetDDGeneratorHeight()       {return fDDGeneratorHeight;};
     G4double           GetInsulatorHeight()       {return fInsulatorHeight;};
     G4double           GetThermalAbsorberHeight()       {return fThermalAbsorberHeight;};
     
     
     void               PrintParameters();
                       
  private:
  
     // world box
     G4VPhysicalVolume* fPWorld;
     G4LogicalVolume*   fLWorld;
     
     G4double           fWorldSize;
     G4Material*        fMaterial;   
     
     G4Material*        fWorldMater;  
     
     // polypropylene insulator
     G4double           fInsulatorThickness; 
     G4double 					fInsulatorLength;
     G4double           fInsulatorWidth;
     G4double 					fInsulatorHeight;
     G4LogicalVolume*   fLogicInsulator;
     G4VPhysicalVolume* fPhysiInsulator;
     G4Material*        fInsulatorMater;
     
     // Feedthrough port
     G4double						fPortHeight;
  	 G4double						fPortRadius;
  	 G4LogicalVolume*   fLogicPort;
     G4VPhysicalVolume* fPhysiPort;
     G4Material*        fPortMater;
     
     // stainless steel cryostat
     G4double           fCryostatThickness; 
     G4double fCryostatLength;
  	 G4double fCryostatWidth ;
  	 G4double fCryostatHeight;
     G4LogicalVolume*   fLogicCryostat;
     G4VPhysicalVolume* fPhysiCryostat;
     G4Material*        fCryostatMater;  
     
     // liquid argon pool
     G4double           fPoolHeight; 
     G4double           fPoolWidth;
     G4double           fPoolLength;
     G4LogicalVolume*   fLogicPool;
     G4VPhysicalVolume* fPhysiPool;
     G4Material*        fPoolMater; 
     
     // gas argon
     G4double           fBufferHeight; 
     G4double           fBufferWidth;
     G4double           fBufferLength;
     G4LogicalVolume*   fLogicBuffer;
     G4VPhysicalVolume* fPhysiBuffer;
     G4Material*        fBufferMater; 
     
     // neutron DD generator
     G4double           fDDGeneratorHeight; //UI
     G4double           fDDGeneratorRadius; //UI
     G4LogicalVolume*   fLogicDDGenerator;
     G4VPhysicalVolume* fPhysiDDGenerator;
     G4Material*        fDDGeneratorMater;  
     
     // Moderator1: 
     G4double           fModerator1_Height;
     G4double           fModerator1_Radius;
     G4LogicalVolume*   fLogicModerator1;
     G4VPhysicalVolume* fPhysiModerator1;
     G4Material*        fModerator1_Mater; 
     
     // Moderator2: 
     G4double           fModerator2_Height;
     G4double           fModerator2_comp1_Height; //UI
     G4double           fModerator2_comp1_Radius; //UI
     G4double           fModerator2_comp2_Height;
  	 G4double           fModerator2_comp2_Rbottom;
  	 G4double           fModerator2_comp2_Rtop;
     G4LogicalVolume*   fLogicModerator2;
     G4VPhysicalVolume* fPhysiModerator2;
     G4Material*        fModerator2_Mater;
     
     // Thermal neutron absorber
  	 G4double						fThermalAbsorberHeight;
  	 G4double						fThermalAbsorberRadius;
  	 G4LogicalVolume*   fLogicThermalAbsorber;
     G4VPhysicalVolume* fPhysiThermalAbsorber;
     G4Material*        fThermalAbsorberMater;
     
     // neutron reflector
     G4double           fReflectorThickness;
     G4double						fReflectorHeight;
     G4double						fReflectorRadius;
     G4LogicalVolume*   fLogicReflector;
     G4VPhysicalVolume* fPhysiReflector;
     G4Material*        fReflectorMater; 
     
     // neutron energy filter
     G4double           fFilterHeight;
  	 G4double           fFilterRadius_top;
  	 G4double           fFilterRadius_bottom; 	 
     G4LogicalVolume*   fLogicFilter;
     G4VPhysicalVolume* fPhysiFilter;
     G4Material*        fFilterMater;
     
     // neutron shield
     G4double						fNShieldThickness;
     G4double           fNShieldHeight;
  	 G4double           fNShieldRadius;
     G4LogicalVolume*   fLogicNShield;
     G4VPhysicalVolume* fPhysiNShield;
     G4Material*        fNShieldMater;    
     
     DetectorMessenger* fDetectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

