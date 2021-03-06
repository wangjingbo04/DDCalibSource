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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Paraboloid.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fPWorld(0), fLWorld(0), fMaterial(0), 
 fPhysiPool(0), fLogicPool(0), fPoolMater(0), 
 fPhysiFilter(0), fLogicFilter(0), fFilterMater(0), 
 fPhysiReflector(0), fLogicReflector(0), fReflectorMater(0), 
 fPhysiDDGenerator(0), fLogicDDGenerator(0), fDDGeneratorMater(0),
 fPhysiNShield(0), fLogicNShield(0), fNShieldMater(0), 
 fDetectorMessenger(0)
{
  fWorldSize = 80*m;
  DefineMaterials();
  SetWorldMaterial("Air");
  
//  // ProtoDUNE Lar pool  
//  fPoolLength      = 8.9*m;
//  fPoolWidth      = 7.8*m;
//  fPoolHeight      = 7.3*m;
//  
//  // ProtoDUNE Gas argon buffer
//  fBufferLength      = 8.9*m; 
//  fBufferWidth      = 7.8*m; 
//  fBufferHeight      = 0.8*m; 
  
  // DUNE Lar pool  
  fPoolLength      = 58.0*m;
  fPoolWidth      = 14.5*m;
  fPoolHeight      = 12.0*m;
  
  // DUNE Gas argon buffer
  fBufferLength      = 58.0*m;
  fBufferWidth      = 14.5*m;
  fBufferHeight      = 0.8*m;
  
  // polypropylene insulator
  fInsulatorThickness = 90*cm;
  fInsulatorLength = fPoolLength + 2*fCryostatThickness + 2*fInsulatorThickness;
  fInsulatorWidth = fPoolWidth + 2*fCryostatThickness + 2*fInsulatorThickness;
  fInsulatorHeight= fPoolHeight + fBufferHeight + 2*fCryostatThickness + 2*fInsulatorThickness;
  
  // Feedthrough port
  fPortHeight = fInsulatorThickness;
  fPortRadius = 25.0*cm;
  
  // stainless steel cryostat
  fCryostatThickness = 1.0*cm;
  fCryostatLength = fPoolLength + 2*fCryostatThickness;
  fCryostatWidth = fPoolWidth + 2*fCryostatThickness;
  fCryostatHeight= fPoolHeight + fBufferHeight + 2*fCryostatThickness;
  
  // neutron DD generator
  fDDGeneratorHeight = 10.0*cm;
  fDDGeneratorRadius = 10.0*cm;
  
  // 1st moderator (Fe)
  fModerator1_Height = 25.0*cm; // must be larger than the DD generator height
  fModerator1_Radius = 20.0*cm;
  
  // 2nd moderator: Fluental 
  fModerator2_Height = 26.0*cm;
  fModerator2_comp1_Height = fModerator2_Height/2;
  fModerator2_comp2_Height = fModerator2_Height/2;
  fModerator2_comp1_Radius = 40.0*cm;
  fModerator2_comp2_Rtop = 40.0*cm;
  fModerator2_comp2_Rbottom = 30*cm;
  
  // Thermal neutron absorber
  fThermalAbsorberHeight = 0.5*cm;
  fThermalAbsorberRadius = 12.5*cm;
  
  // neutron energy filter
  fFilterHeight    = 10.0*cm;
  fFilterRadius_top    = 30.0*cm;
  fFilterRadius_bottom    = 12.5*cm;
  
  // neutron reflector
  fReflectorThickness = 20.0*cm;
  fReflectorHeight = fFilterHeight + fModerator2_comp1_Height + fModerator2_comp2_Height + fModerator1_Height + fReflectorThickness;
  fReflectorRadius = fModerator2_comp1_Radius + fReflectorThickness;
  
  // neutron shield
  fNShieldThickness = 20.0*cm;
  fNShieldHeight = fThermalAbsorberHeight + fReflectorHeight + fNShieldThickness;
  fNShieldRadius = fReflectorRadius + fNShieldThickness;
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // specific element name for thermal neutronHP
  // (see G4ParticleHPThermalScatteringNames.cc)

  G4int Z, A, a, ncomponents, natoms;
  G4double fractionmass, abundance;
 
  // world material
  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen",   "O", 8, 16.00*g/mole);     
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2, kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);


  
  // vacuum
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  G4Material* Vacuum = new G4Material("interGalactic", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);

 
  
  // Define elements for all materials not found in the NIST database  
  G4NistManager* man = G4NistManager::Instance();
  G4Element* Li = man->FindOrBuildElement("Li");
  G4Element* B = man->FindOrBuildElement("B");
  G4Element* C  = man->FindOrBuildElement("C");
  G4Element* F = man->FindOrBuildElement("F");
  G4Element* Na = man->FindOrBuildElement("Na");
  G4Element* Mg = man->FindOrBuildElement("Mg");
  G4Element* Al = man->FindOrBuildElement("Al");
  G4Element* Si = man->FindOrBuildElement("Si");
  G4Element* K = man->FindOrBuildElement("K");
  G4Element* Ca = man->FindOrBuildElement("Ca");
  G4Element* Ti = man->FindOrBuildElement("Ti");
  G4Element* Cr = man->FindOrBuildElement("Cr");
  G4Element* Mn = man->FindOrBuildElement("Mn");
  G4Element* Fe = man->FindOrBuildElement("Fe");
  G4Element* Ni = man->FindOrBuildElement("Ni");
  G4Element* Sb = man->FindOrBuildElement("Sb");
  G4Element* Xe = man->FindOrBuildElement("Xe");
  G4Element* Cs = man->FindOrBuildElement("Cs");
  G4Element* Bi = man->FindOrBuildElement("Bi");
  

  
  // stainless steel
  G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
      StainlessSteel->AddElement(C, fractionmass=0.015); 
      StainlessSteel->AddElement(Si, fractionmass=0.008);
      StainlessSteel->AddElement(Cr, fractionmass=0.18);
      StainlessSteel->AddElement(Mn, fractionmass=0.01);
      StainlessSteel->AddElement(Fe, fractionmass=0.697);
      StainlessSteel->AddElement(Ni, fractionmass=0.09);

	
  // MgF2
  G4Material* MgF2 = new G4Material("MgF2", 3.15*g/cm3, ncomponents=2, kStateSolid);
      MgF2->AddElement(Mg, natoms=1);
      MgF2->AddElement(F, natoms=2);

  
  // TiF3
  G4Material* TiF3 = new G4Material("TiF3", 3.4*g/cm3, ncomponents=2, kStateSolid);
      TiF3->AddElement(Ti, natoms=1);
      TiF3->AddElement(F, natoms=3);

	
  // Fe-56 isotope
  G4Isotope* iso_Fe = new G4Isotope("iso_Fe", Z=26, A=56, a=55.9349363*g/mole);
  G4Element* ele_Fe = new G4Element("ele_Fe", "Fe", ncomponents=1);
  ele_Fe->AddIsotope(iso_Fe,abundance=100.*perCent);
  G4Material* mat_Fe=new G4Material("mat_Fe",7.874*g/cm3, ncomponents = 1);
	mat_Fe->AddElement(ele_Fe, fractionmass = 1 );

	
  // Li-6 isotope
  G4Isotope* iso_Li = new G4Isotope("iso_Li", Z=3, A=6, a=6.015122795*g/mole);
  G4Element* ele_Li = new G4Element("ele_Li", "Li", ncomponents=1);
  ele_Li->AddIsotope(iso_Li,abundance=100.*perCent);
  G4Material* mat_Li=new G4Material("mat_Li",0.534*g/cm3, ncomponents = 1);
    mat_Li->AddElement(ele_Li, fractionmass = 1 );

	
  
  // Fluental
  G4Material* AlF3 = new G4Material("AlF3", 2.88*g/cm3, ncomponents=2, kStateSolid); //ALF3
      AlF3->AddElement(Al, natoms=1);
      AlF3->AddElement(F, natoms=3);
  
  G4Material* LiF = new G4Material("LiF", 2.64*g/cm3, ncomponents=2, kStateSolid);   //LiF
      LiF->AddElement(Li, natoms=1);
      LiF->AddElement(F, natoms=1);
  
  G4Material* Fluental = new G4Material ("Fluental", density=2.831, ncomponents=3);
      Fluental->AddElement (Al, fractionmass = 30*perCent);
      Fluental->AddMaterial (AlF3, fractionmass = 69*perCent);
      Fluental->AddMaterial (LiF, fractionmass = 1*perCent);



  //Teflon (CF2)
  G4Material* Teflon = new G4Material("Teflon", 2.20*g/cm3, ncomponents=2, kStateSolid);
      Teflon->AddElement(C, natoms=1);
      Teflon->AddElement(F, natoms=2);



  //Al2O3
  G4Material* Al2O3 = new G4Material("Al203", 3.987*g/cm3, ncomponents=2, kStateSolid);
      Al2O3->AddElement(Al, natoms=2);
      Al2O3->AddElement(O, natoms=3);



  //BiF3
  G4Material* BiF3 = new G4Material("BiF3", 5.32*g/cm3, ncomponents=2, kStateSolid);
      BiF3->AddElement(Bi, natoms=1);
      BiF3->AddElement(F, natoms=3);



  //BiF5
  G4Material* BiF5 = new G4Material("BiF5", 5.40*g/cm3, ncomponents=2, kStateSolid);
      BiF5->AddElement(Bi, natoms=1);
      BiF5->AddElement(F, natoms=5);



  //CaF2
  G4Material* CaF2 = new G4Material("CaF2", 3.18*g/cm3, ncomponents=2, kStateSolid);
      CaF2->AddElement(Ca, natoms=1);
      CaF2->AddElement(F, natoms=2);



  
  
  // world mater
  fWorldMater = Air20;
  
  // insulator
  fInsulatorMater = man->FindOrBuildMaterial("G4_POLYETHYLENE");

  // Feedthrough port
  fPortMater = Air20;
  
  // cryostat
  fCryostatMater = StainlessSteel;
  
  // liquid argon pool
  fPoolMater = man->FindOrBuildMaterial("G4_lAr");
  
  // gas argon buffer
  fBufferMater = man->FindOrBuildMaterial("G4_Ar");
  
  // neutron reflectorMater
  fReflectorMater = man->FindOrBuildMaterial("G4_Pb");
  
  // DD generator
  fDDGeneratorMater = Vacuum; 
  
  // Moderator1: 
  fModerator1_Mater = mat_Fe;
  
  // Moderator2: 
  fModerator2_Mater = man->FindOrBuildMaterial("G4_S");;
   
  // neutron energy filter
  fFilterMater = man->FindOrBuildMaterial("G4_S");;
  
  // neutron thermal absorber
  fThermalAbsorberMater = Vacuum;
  
  // neutron shield
  fNShieldMater = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  
 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

 
  
  // world box
  G4Box*
  sWorld = new G4Box("Container",                                 //its name
                     fWorldSize/2,fWorldSize/2,fWorldSize/2);     //its dimensions

  fLWorld = new G4LogicalVolume(sWorld,                           //its shape
                                fWorldMater,                      //its material
                                "World_l");                       //its name

  fPWorld = new G4PVPlacement(0,                                  //no rotation
                              G4ThreeVector(),                    //at (0,0,0)
                              fLWorld,                            //its logical volume
                              "World_p",                          //its name
                              0,                                  //its mother  volume
                              false,                              //no boolean operation
                              0);                                 //copy number




  // insulator
  
  G4Box*
  sInsulator = new G4Box("Insulator",                                                 //its name
                         fInsulatorLength/2,fInsulatorWidth/2, fInsulatorHeight/2);   //its dimensions
  
  fLogicInsulator = new G4LogicalVolume(sInsulator,                     //its shape
                                        fInsulatorMater,                //its material
                                        "Insulator_l");                 //its name

  fPhysiInsulator = new G4PVPlacement(0,                          //no rotation
                                      G4ThreeVector(0, 0, 0),     //its placement
                                      fLogicInsulator,            //its logical volume
                                      "Insulator_p",              //its name
                                      fLWorld,                    //its mother  volume
                                      false,                      //no boolean operation
                                      0);                         //copy number



 
  // Feedthrough port
  G4Tubs*
  sPort = new G4Tubs("Port",                                                                //its name
                     0, 2*fFilterRadius_bottom, fInsulatorThickness/2, 0.,CLHEP::twopi );   //its dimensions
  G4ThreeVector zTransInsulator(0, 0, fInsulatorHeight/2 - fInsulatorThickness/2); 
  fLogicPort = new G4LogicalVolume(sPort,                      //its shape
                                   fPortMater,                 //its material
                                   "Port_l");                  //its name

  fPhysiPort = new G4PVPlacement(0,                                                               //no rotation
                                 G4ThreeVector(0, 0, fInsulatorHeight/2 - fPortHeight/2),         //at (0,0,0)
                                 fLogicPort,                                                      //its logical volume
                                 "Port_p",                                                        //its name
                                 fLogicInsulator,                                                 //its mother  volume
                                 false,                                                           //no boolean operation
                                 0);                                                              //copy number



  
  // cryostat
  
  G4Box*
  sCryostat = new G4Box("Cryostat",                                                      //its name
                        fCryostatLength/2, fCryostatWidth/2, fCryostatHeight/2);         //its dimensions 
  fLogicCryostat = new G4LogicalVolume(sCryostat,                                        //its shape
                                       fCryostatMater,                                   //its material
                                       "Cryostat_l");                                    //its name

  fPhysiCryostat = new G4PVPlacement(0,                                                  //no rotation
                                     G4ThreeVector(0, 0, 0),                             //at (0,0,0)
                                     fLogicCryostat,                                     //its logical volume
                                     "Cryostat_p",                                       //its name
                                     fLogicInsulator,                                    //its mother  volume
                                     false,                                              //no boolean operation
                                     0);                                                 //copy number


 
                            
  // liquid argon pool
  G4Box*
  sLarPool = new G4Box("Larpool",                                   //its name
                       fPoolLength/2,fPoolWidth/2,fPoolHeight/2);   //its dimensions

  fLogicPool = new G4LogicalVolume(sLarPool,                        //its shape
                                   fPoolMater,                      //its material
                                   "LarPool_l");                    //its name

  fPhysiPool = new G4PVPlacement(0,                                                                             //no rotation
                                 G4ThreeVector(0, 0, -fCryostatHeight/2+fCryostatThickness+fPoolHeight/2),            
                                 fLogicPool,                                                                    //its logical volume
                                 "LarPool_p",                                                                   //its name
                                 fLogicCryostat,                                                                //its mother  volume
                                 false,                                                                         //no boolean operation
                                 0);                                                                            //copy number



  	                        
  
  // gas argon buffer
  G4Box*
  sGArBuffer = new G4Box("GArbuffer",                                       //its name
                         fBufferLength/2,fBufferWidth/2,fBufferHeight/2);   //its dimensions

  fLogicBuffer = new G4LogicalVolume(sGArBuffer,                            //its shape
                                     fBufferMater,                          //its material
                                     "GArbuffer_l");                        //its name

  fPhysiBuffer = new G4PVPlacement(0,                                                                                           //no rotation
                                   G4ThreeVector(0, 0, -fCryostatHeight/2+fCryostatThickness+fPoolHeight + fBufferHeight/2),    //at (0,0,0)
                                   fLogicBuffer,                                                                                //its logical volume
                                   "GArbuffer_p",                                                                               //its name
                                   fLogicCryostat,                                                                              //its mother  volume
                                   false,                                                                                       //no boolean operation
                                   0);                                                                                          //copy number   


    
                            
                            
  // neutron shield
  
  G4Tubs* 
  sNShield = new G4Tubs("NShield_s",                                                //its name
                        0, fNShieldRadius, 0.5*fNShieldHeight, 0.,CLHEP::twopi);    //its dimensions

  fLogicNShield = new G4LogicalVolume(sNShield,                                     //its shape
                                      fNShieldMater,                                //its material
                                      "NShield_l");                                 //its name

  fPhysiNShield = new G4PVPlacement(0,                                                              //no rotation
                                    G4ThreeVector(0, 0,  fInsulatorHeight/2 + fNShieldHeight/2),      
                                    fLogicNShield ,                                                 //its logical volume
                                    "NShield_p",                                                    //its name
                                    fLWorld,                                                        //its mother  volume
                                    false,                                                          //no boolean operation
                                    0);                                                             //copy number


          
  
  // neutron thermal absorber                        
  G4Tubs* 
  sThermalAbsorber = new G4Tubs("ThermalAbsorber_s",  
                                0, fThermalAbsorberRadius, 0.5*fThermalAbsorberHeight, 0.,CLHEP::twopi);


  fLogicThermalAbsorber = new G4LogicalVolume(sThermalAbsorber,                 //shape
                                              fThermalAbsorberMater,            //materiald
                                              "ThermalAbsorber_l");             //name
                               
  fPhysiThermalAbsorber = new G4PVPlacement(0,                                                                   //no rotation
                                            G4ThreeVector(0, 0, -fNShieldHeight/2 + fThermalAbsorberHeight/2),  
                                            fLogicThermalAbsorber,                                               //logical volume
                                            "ThermalAbsorber_p",                                                 //name
                                            fLogicNShield,                                                       //mother  volume
                                            false,                                                               //no boolean operation
                                            0);                                                                  //copy number



  
  // neutron reflector
  G4Tubs* 
  sReflector = new G4Tubs("Reflector_s",                                                 //its name
                          0, fReflectorRadius, 0.5*fReflectorHeight, 0.,CLHEP::twopi);   //its dimensions

  fLogicReflector = new G4LogicalVolume(sReflector,                                      //its shape
                                        fReflectorMater,                                 //its material
                                        "Reflector_l");                                  //its name

  fPhysiReflector = new G4PVPlacement(0,                                                 //no rotation
                                      G4ThreeVector(0, 0, -fNShieldHeight/2 + fThermalAbsorberHeight + fReflectorHeight/2),   
                                      fLogicReflector ,                                  //its logical volume
                                      "Reflector_p",                                     //its name
                                      fLogicNShield,                                     //its mother  volume
                                      false,                                             //no boolean operation
                                      0);                                                //copy number






  
  // neutron energy filter                          
  G4Cons* 
  sFilter = new G4Cons("Filter_s",  
                       0, fFilterRadius_bottom, 0, fFilterRadius_top, 0.5*fFilterHeight, 0.,CLHEP::twopi);


  fLogicFilter = new G4LogicalVolume(sFilter,                 //shape
                                     fFilterMater,            //materiald
                                     "Filter_l");             //name
                               
  fPhysiFilter = new G4PVPlacement(0,                         //no rotation
                                   G4ThreeVector(0, 0, -fReflectorHeight/2 + fFilterHeight/2),  //at (0,0,0)
                                   fLogicFilter,              //logical volume
                                   "Filter_p",                //name
                                   fLogicReflector,           //mother  volume
                                   false,                     //no boolean operation
                                   0);                        //copy number



                 
  // neutron Moderator2                         
  G4Tubs* 
  sModerator2_comp1 = new G4Tubs("Moderator2_comp1_s",  
                                 0, fModerator2_comp1_Radius, 0.5*fModerator2_comp1_Height, 0.,CLHEP::twopi);
  G4Cons* 
  sModerator2_comp2 = new G4Cons("Moderator2_comp2_s",  
                                 0, fModerator2_comp2_Rbottom, 0, fModerator2_comp2_Rtop, 0.5*fModerator2_comp2_Height, 0.,CLHEP::twopi);
                
  G4ThreeVector zTransModerator2(0, 0, - fModerator2_comp1_Height/2 - fModerator2_comp2_Height/2); 
  
  G4UnionSolid* sModerator2 = new G4UnionSolid("Moderator2_comp1_s + Moderator2_comp2_s", 
  											   sModerator2_comp1, sModerator2_comp2, 0, zTransModerator2); 

  fLogicModerator2 = new G4LogicalVolume(sModerator2,                   //shape
                                         fModerator2_Mater,             //materiald
                                         "Moderator2_l");               //name
                               
  fPhysiModerator2 = new G4PVPlacement(0,                               //no rotation
                                       G4ThreeVector(0, 0, -fReflectorHeight/2 + fFilterHeight + fModerator2_comp1_Height/2 + fModerator2_comp2_Height), 
                                       fLogicModerator2,                //logical volume
                                       "Moderator2_p",                  //name
                                       fLogicReflector,                 //mother  volume
                                       false,                           //no boolean operation
                                       0);                              //copy number




  
  // neutron Moderator1                         
  G4Tubs* 
  sModerator1 = new G4Tubs("Moderator1_s",  
                           0, fModerator1_Radius, 0.5*fModerator1_Height, 0.,CLHEP::twopi);


  fLogicModerator1 = new G4LogicalVolume(sModerator1,                   //shape
                                         fModerator1_Mater,             //materiald
                                         "Degrader_l");                 //name
                               
  fPhysiModerator1 = new G4PVPlacement(0,                               //no rotation
                                       G4ThreeVector(0, 0, -fReflectorHeight/2 + fFilterHeight + fModerator2_comp1_Height + fModerator2_comp2_Height + fModerator1_Height/2),             //at (0,0,0)
                                       fLogicModerator1,                //logical volume
                                       "Moderator1_p",                  //name
                                       fLogicReflector,                 //mother  volume
                                       false,                           //no boolean operation
                                       0);                              //copy number



                            
  // neutron DD generator                         
  G4Tubs* 
  sDDGenerator = new G4Tubs("DDGenerator_s",  
                            0, fDDGeneratorRadius, 0.5*fDDGeneratorHeight, 0.,CLHEP::twopi);

  fLogicDDGenerator = new G4LogicalVolume(sDDGenerator,                 //shape
                                          fDDGeneratorMater,            //materiald
                                          "DDGenerator_l");             //name
                               
  fPhysiDDGenerator = new G4PVPlacement(0,                              //no rotation
                                        G4ThreeVector(0, 0, fModerator1_Height/2 - fDDGeneratorHeight/2),
                                        fLogicDDGenerator,              //logical volume
                                        "DDGenerator_p",                //name
                                        fLogicModerator1,               //mother  volume
                                        false,                          //no boolean operation
                                        0);                             //copy number




  //--------- Visualization attributes -------------------------------
  fLWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
  	
  G4VisAttributes* VisAttPool= new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  fLogicPool->SetVisAttributes(VisAttPool);
  
  G4VisAttributes* VisAttBuffer= new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
  fLogicBuffer->SetVisAttributes(VisAttBuffer);
  
  G4VisAttributes* VisAttFilter= new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  fLogicFilter->SetVisAttributes(VisAttFilter);
  
  G4VisAttributes* VisAttReflector= new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  fLogicReflector->SetVisAttributes(VisAttReflector);

  G4VisAttributes* VisAttModerator1= new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  fLogicModerator1->SetVisAttributes(VisAttModerator1);
  
  G4VisAttributes* VisAttModerator2= new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
  fLogicModerator2->SetVisAttributes(VisAttModerator2);
  
  G4VisAttributes* VisAttDDGenerator = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0));
  fLogicDDGenerator->SetVisAttributes(VisAttDDGenerator);

  
//-------------------------------------------------------------------                         
  PrintParameters();
  
  //always return the root volume
  //
  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The World is " << G4BestUnit(fWorldSize,"Length")
         << " of " << fMaterial->GetName() 
         << "\n \n" << fMaterial << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fMaterial != pttoMaterial) {
      fMaterial = pttoMaterial;
      if(fLWorld) { fLWorld->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetWorldMaterial : "
           << materialChoice << " not found" << G4endl;
  }              
}





//SET and GET functions for moderator layers


void DetectorConstruction::GetMaterialTable(){

    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

void DetectorConstruction::SetModerator1Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fModerator1_Mater != pttoMaterial) {
      fModerator1_Mater = pttoMaterial;
      if(fLogicModerator1) { fLogicModerator1->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModerator1Material : "
           << materialChoice << " not found" << G4endl;
  }              
}




void DetectorConstruction::GetModerator1Material(){

    G4cout << "Moderator1 material is " << fModerator1_Mater->GetName() << G4endl;
}




void DetectorConstruction::SetModerator2Material(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fModerator2_Mater != pttoMaterial) {
      fModerator2_Mater = pttoMaterial;
      if(fLogicModerator2) { fLogicModerator2->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModerator1Material : "
           << materialChoice << " not found" << G4endl;
  }              
}



void DetectorConstruction::GetModerator2Material(){

    G4cout << "Moderator2 material is " << fModerator2_Mater->GetName() << G4endl;
}





void DetectorConstruction::SetFilterMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fFilterMater != pttoMaterial) {
      fFilterMater = pttoMaterial;
      if(fFilterMater) { fLogicFilter->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModerator1Material : "
           << materialChoice << " not found" << G4endl;
  }              
}




void DetectorConstruction::GetFilterMaterial(){

    G4cout << "Filter material is " << fFilterMater->GetName() << G4endl;
}




void DetectorConstruction::SetAbsorberMaterial(G4String materialChoice)
{
  // search the material by its name
  G4Material* pttoMaterial =
     G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
  
  if (pttoMaterial) { 
    if(fThermalAbsorberMater != pttoMaterial) {
      fThermalAbsorberMater = pttoMaterial;
      if(fThermalAbsorberMater) { fLogicThermalAbsorber->SetMaterial(pttoMaterial); }
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  } else {
    G4cout << "\n--> warning from DetectorConstruction::SetModerator1Material : "
           << materialChoice << " not found" << G4endl;
  }              
}




void DetectorConstruction::GetAbsorberMaterial(){

    G4cout << "Absorber material is " << fThermalAbsorberMater->GetName() << G4endl;
}



//END OF SET AND GET FUNCTIONS













//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fWorldSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetFilterHeight(G4double value)
{
  fFilterHeight = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetModerator1Height(G4double value)
{
  fModerator1_Height = value;
  fReflectorHeight = fFilterHeight + fModerator2_comp1_Height + fModerator2_comp2_Height + fModerator1_Height + fReflectorThickness;
  fNShieldHeight = fThermalAbsorberHeight + fReflectorHeight + fNShieldThickness;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

void DetectorConstruction::SetModerator2Height(G4double value)
{
  fModerator2_Height = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

