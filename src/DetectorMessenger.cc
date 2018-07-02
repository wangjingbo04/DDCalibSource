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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
// $Id: DetectorMessenger.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction * Det)
:G4UImessenger(), 
 fDetector(Det), fTestemDir(0), fDetDir(0), fMaterCmd(0), fSizeCmd(0),
 fIsotopeCmd(0)
{ 
  fTestemDir = new G4UIdirectory("/ddsource/");
  fTestemDir->SetGuidance("commands specific to this example");
  
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/ddsource/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
        
  fMaterCmd = new G4UIcmdWithAString("/ddsource/det/setMat",this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice",false);
  fMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);


  //Get and Set Material Functions
  
  fModerator1GetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetMaterialTable", this);
  fModerator1GetMaterCmd->SetGuidance("Parameterless, returns a table of materials used.");
  fModerator1GetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fModerator1SetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetModerator1Material",this);
  fModerator1SetMaterCmd->SetGuidance("Select material of the moderator1.");
  fModerator1SetMaterCmd->SetParameterName("choice",false);
  fModerator1SetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fModerator1GetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetModerator1Material", this);
  fModerator1GetMaterCmd->SetGuidance("Parameterless, returns Moderator1 material.");
  fModerator1GetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fModerator2SetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetModerator2Material",this);
  fModerator2SetMaterCmd->SetGuidance("Select material of the moderator2.");
  fModerator2SetMaterCmd->SetParameterName("choice",false);
  fModerator2SetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fModerator2GetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetModerator2Material", this);
  fModerator2GetMaterCmd->SetGuidance("Parameterless, returns Moderator2 material.");
  fModerator2GetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilterSetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetFilterMaterial",this);
  fFilterSetMaterCmd->SetGuidance("Select material of the Filter.");
  fFilterSetMaterCmd->SetParameterName("choice",false);
  fFilterSetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fFilterGetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetFilterMaterial", this);
  fFilterGetMaterCmd->SetGuidance("Parameterless, returns Filter material.");
  fFilterGetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAbsorberSetMaterCmd = new G4UIcmdWithAString("/ddsource/det/SetAbsorberMaterial",this);
  fAbsorberSetMaterCmd->SetGuidance("Select material of the Absorber.");
  fAbsorberSetMaterCmd->SetParameterName("choice",false);
  fAbsorberSetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fAbsorberGetMaterCmd = new G4UIcmdWithoutParameter("/ddsource/det/GetAbsorberMaterial", this);
  fAbsorberGetMaterCmd->SetGuidance("Parameterless, returns Absorber material.");
  fAbsorberGetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

 
  
  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/setSize",this);
  fSizeCmd->SetGuidance("Set size of the box");
  fSizeCmd->SetParameterName("Size",false);
  fSizeCmd->SetRange("Size>0.");
  fSizeCmd->SetUnitCategory("Length");
  fSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fFilterHeightCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/setFilterHeight",this);
  fFilterHeightCmd->SetGuidance("Set Height of the filter");
  fFilterHeightCmd->SetParameterName("FilterHeight",false);
  fFilterHeightCmd->SetRange("FilterHeight>0.");
  fFilterHeightCmd->SetUnitCategory("Length");
  fFilterHeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fModerator1HeightCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/setModerator1Height",this);
  fModerator1HeightCmd->SetGuidance("Set Height of the 1st moderator");
  fModerator1HeightCmd->SetParameterName("Moderator1Height",false);
  fModerator1HeightCmd->SetRange("Moderator1Height>0.");
  fModerator1HeightCmd->SetUnitCategory("Length");
  fModerator1HeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fModerator2HeightCmd = new G4UIcmdWithADoubleAndUnit("/ddsource/det/setModerator2Height",this);
  fModerator2HeightCmd->SetGuidance("Set Height of the 2nd moderator");
  fModerator2HeightCmd->SetParameterName("Moderator2Height",false);
  fModerator2HeightCmd->SetRange("Moderator2Height>0.");
  fModerator2HeightCmd->SetUnitCategory("Length");
  fModerator2HeightCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
       
  fIsotopeCmd = new G4UIcommand("/ddsource/det/setIsotopeMat",this);
  fIsotopeCmd->SetGuidance("Build and select a material with single isotope");
  fIsotopeCmd->SetGuidance("  symbol of isotope, Z, A, density of material");
  //
  G4UIparameter* symbPrm = new G4UIparameter("isotope",'s',false);
  symbPrm->SetGuidance("isotope symbol");
  fIsotopeCmd->SetParameter(symbPrm);
  //      
  G4UIparameter* ZPrm = new G4UIparameter("Z",'i',false);
  ZPrm->SetGuidance("Z");
  ZPrm->SetParameterRange("Z>0");
  fIsotopeCmd->SetParameter(ZPrm);
  //      
  G4UIparameter* APrm = new G4UIparameter("A",'i',false);
  APrm->SetGuidance("A");
  APrm->SetParameterRange("A>0");
  fIsotopeCmd->SetParameter(APrm);  
  //    
  G4UIparameter* densityPrm = new G4UIparameter("density",'d',false);
  densityPrm->SetGuidance("density of material");
  densityPrm->SetParameterRange("density>0.");
  fIsotopeCmd->SetParameter(densityPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of density");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm->SetParameterCandidates(unitList);
  fIsotopeCmd->SetParameter(unitPrm);
  //
  fIsotopeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fMaterCmd;

  //Material Commands
  delete fModerator1SetMaterCmd;
  delete fModerator1GetMaterCmd;
  delete fModerator2SetMaterCmd;
  delete fModerator2GetMaterCmd;
  delete fFilterSetMaterCmd;
  delete fFilterGetMaterCmd;
  delete fAbsorberSetMaterCmd;
  delete fAbsorberGetMaterCmd;

  delete fSizeCmd;
  delete fIsotopeCmd;
  delete fDetDir;
  delete fTestemDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fMaterCmd )
   { fDetector->SetWorldMaterial(newValue);}

  
  //Set and Get Material Commands
  if( command == fModerator1SetMaterCmd )
   { fDetector->SetModerator1Material(newValue);}

  if( command == fModerator1GetMaterCmd )
   { fDetector->GetModerator1Material();}

  if( command == fModerator2SetMaterCmd )
   { fDetector->SetModerator2Material(newValue);}

  if( command == fModerator2GetMaterCmd )
   { fDetector->GetModerator2Material();}

  if( command == fFilterSetMaterCmd )
   { fDetector->SetFilterMaterial(newValue);}

  if( command == fFilterGetMaterCmd )
   { fDetector->GetFilterMaterial();}

  if( command == fAbsorberSetMaterCmd )
   { fDetector->SetAbsorberMaterial(newValue);}

 if( command == fAbsorberGetMaterCmd )
   { fDetector->GetAbsorberMaterial();}



  if( command == fSizeCmd )
   { fDetector->SetSize(fSizeCmd->GetNewDoubleValue(newValue));}
   
  if( command == fFilterHeightCmd )
   { fDetector->SetFilterHeight(fSizeCmd->GetNewDoubleValue(newValue));}

  if( command == fModerator1HeightCmd )
   { fDetector->SetModerator1Height(fSizeCmd->GetNewDoubleValue(newValue));}
  
  if( command == fModerator2HeightCmd )
   { fDetector->SetModerator2Height(fSizeCmd->GetNewDoubleValue(newValue));}
        
  if (command == fIsotopeCmd)
   {
     G4int Z; G4int A; G4double dens;
     G4String name, unt;
     std::istringstream is(newValue);
     is >> name >> Z >> A >> dens >> unt;
     dens *= G4UIcommand::ValueOf(unt);
     fDetector->MaterialWithSingleIsotope (name,name,dens,Z,A);
     fDetector->SetWorldMaterial(name);    
   }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  fFilterGetMaterCmd->SetGuidance("Parameterless, returns Filter material.");
