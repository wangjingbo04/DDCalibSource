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
/// \file HistoManager.cc
/// \brief Implementation of the HistoManager class
//
// $Id: HistoManager.cc 67909 2013-03-12 18:51:09Z vnivanch $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "HistoManager.hh"
#include "G4UnitsTable.hh"
#include "g4root.hh"
#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("DDsource.root")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms
  
/*  // Define 1D histograms start values
  const G4int kMaxHisto = 12;
  const G4String id[] = {"0","N collision","Track length","ToF","N collision","Track length","ToF","E at all collision","E degrader1","E degrader2","E filter","E LAr pool"};
  const G4String title[] = 
                { "dummy",                                           //0
                  "incident neutron: nb of collisions above 57 keV",   //1
                  "incident neutron: total track length above 57 keV", //2
                  "incident neutron: time of flight above 57 keV",     //3
                  "incident neutron: nb of collisions below 57 keV",   //4
                  "incident neutron: total track length below 57 keV", //5
                  "incident neutron: time of flight below 57 keV",     //6
                  "incident neutron: energy distribution for all steps", //7
                  "neutron entering degrader1",                				 //8
                  "neutron entering degrader2",                				 //9
                  "neutron entering LAr filter",                				 //10
                  "neutron entering LAr pool"                				 //11
                 };  

  // Default values (to be reset via /analysis/h1/set command)               
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  // Create all histograms as inactivated 
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto; k++) {
    G4int ih = analysisManager->CreateH1(id[k], title[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, true);
  }*/
  
  
  // histos 1D
  //
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 100.;

  G4double xmin = -180.;
  G4double xmax = 180.; 

  //Histogram 0 - nb of collisions above 57 keV
  G4int ih = analysisManager->CreateH1("h1.0", "incident neutron: nb of collisions above 57 keV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 1 - total track length above 57 keV
  ih = analysisManager->CreateH1("h1.1", "incident neutron: total track length above 57 keV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 2 - time of light above 57 keV
  ih = analysisManager->CreateH1("h1.2", "incident neutron: time of flight above 57 keV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 3 - nb of collisions bellow 57 keV
  ih = analysisManager->CreateH1("h1.3", "incident neutron: nb of collisions below 57 keV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 4 - total track length below 57 keV
  ih = analysisManager->CreateH1("h1.4", "incident neutron: total track length below 57 keV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 5 - time of flight below 57 keV
  ih = analysisManager->CreateH1("h1.5", "incident neutron: time of flight below 57 keV", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 6 - energy distribution for all steps
  ih = analysisManager->CreateH1("h1.6", "Neutron energy distribution for all steps", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 7 - energy spectrum of neutrons from generator
  ih = analysisManager->CreateH1("h1.7", "Energy Spectrum from Generator", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 8 - energy spectrum of neutrons from moderator
  ih = analysisManager->CreateH1("h1.8", "Energy Spectrum from Moderator", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 9 - energy spectrum of neutrons from filter 1
  ih = analysisManager->CreateH1("h1.9",  "Energy Spectrum from Filter 1", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 10 - energy spectrum of neutrons from filter 2
  ih = analysisManager->CreateH1("h1.10", "Energy Spectrum from Filter 2", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 11 - energy spectrum of neutrons from filter 3
  ih = analysisManager->CreateH1("h1.11", "Energy Spetrum from Fileter 3", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 12 - energy spectrum of neutrons from thermal absorber
  ih = analysisManager->CreateH1("h1.12", "Energy Spectrum from Thermal Absorber", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 13 - energy spectrum of neutrons entering liquid argon pool (from gas buffer)
  ih = analysisManager->CreateH1("h1.13", "Neutrons entering liquid Argon pool", nbins, vmin, vmax);
  analysisManager->SetH1Activation(ih, true);

  //Histogram 14 - Theta angle for neutrons entering the port
  ih = analysisManager->CreateH1("h1.14", "Theta angle for neutrons entering the gas buffer", 360, xmin, xmax, "degree");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 15 - Phi angle for neutrons entering the port
  ih = analysisManager->CreateH1("h1.15", "Phi angle for neutrons entering the gas buffer", 360, xmin, xmax, "degree");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 16 - DD gun neutron initial momentum angle
  ih = analysisManager->CreateH1("h1.16", "DD gun neutrons initial momentum angle", 360, xmin, xmax, "degree");
  analysisManager->SetH1Activation(ih, true);




  
  // histos 2D
  //
  ih = analysisManager->CreateH2("h2.0","neutron capture position in the moderator (top view, y:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  ih = analysisManager->CreateH2("h2.1","neutron capture position in the moderator (side view, z:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  
  ih = analysisManager->CreateH2("h2.2","neutron capture position in filter 1 (top view, y:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  ih = analysisManager->CreateH2("h2.3","neutron capture position in filter 1 (side view, z:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  
  ih = analysisManager->CreateH2("h2.4","neutron capture position in filter 2 (top view, y:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  ih = analysisManager->CreateH2("h2.5","neutron capture position in filter 2 (side view, z:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);

  ih = analysisManager->CreateH2("h2.6","neutron capture position in filter 3 (top view, y:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  ih = analysisManager->CreateH2("h2.7","neutron capture position in filter 3 (side view, z:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  
  ih = analysisManager->CreateH2("h2.8","neutron capture position in  (top view, y:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);
  ih = analysisManager->CreateH2("h2.9","neutron capture position in LAr TPC (side view, z:x)",nbins,vmin,vmax, nbins,vmin,vmax);
  analysisManager->SetH2Activation(ih, true);






  //Gamma Histograms
    
  G4double maxE = 3;
  G4double minE = 0;
  
  //Histogram 17 - Total gamma energy spectrum
  ih = analysisManager->CreateH1("h1.17", "Total gamma energy spectrum", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 18 - Energy Spectrum of Escaping Gammas
  ih = analysisManager->CreateH1("h1.18", "Energy Spectrum of Escaped Gammas", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 19 - Energy Spectrum of Gammas Entering Shield
  ih = analysisManager->CreateH1("h1.19", "Energy Spectrum of Gammas Entering Shield", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 20 - Energy Spectrum of Gammas Escaping Shield
  ih = analysisManager->CreateH1("h1.20", "Energy Spectrum of Gammas Escaping Shield", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 21 - Energy Spectrum of Gammas Entering LAr Pool
  ih = analysisManager->CreateH1("h1.21", "Energy Spectrum of Gammas Entering LAr Pool", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 22 - Energy Spectrum of Gammas Exiting LAr Pool
  ih = analysisManager->CreateH1("h1.22", "Energy Spectrum of Gammas Exiting LAr Pool", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);

    
  //Neutron Radiation Histograms

  //Histogram 23 - Energy Spectrum of Neutrons Entering World
  ih = analysisManager->CreateH1("h1.23", "Energy Spectrum of Escaping Shield", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 24 - Energy Spectrum of Neutrons Exiting Shield
  ih = analysisManager->CreateH1("h1.24", "Energy Spectrum of Neutrons Entering World", nbins, minE, maxE, "MeV");
  analysisManager->SetH1Activation(ih, true);


  //Neutron Angles into LAr Pool
  //Histogram 25 - Theta angle for neutrons entering LAr Pool
  ih = analysisManager->CreateH1("h1.25", "Theta angle for neutrons entering LAr Pool", 360, -180, 180, "degree");
  analysisManager->SetH1Activation(ih, true);

  //Histogram 26 - Phi angle for neutrons entering LAr Pool
  ih = analysisManager->CreateH1("h1.26", "Phi angle for neutrons entering LAr Pool", 360, -180, 180, "degree");
  analysisManager->SetH1Activation(ih, true);

  

   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
