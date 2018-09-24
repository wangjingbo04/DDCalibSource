
#ifndef TreeMaker_HH
#define TreeMaker_HH

#include "G4Track.hh"
#include "globals.hh"
#include "Run.hh"
#include "TFile.h"
#include "TTree.h"
#include "RecoElectron.hh"
#include <string>
#include <iostream>
#include <vector>


class TreeMaker {

 public:
 	static TreeMaker* Instance();
  TreeMaker();
  ~TreeMaker();
  
  /// \brief Reset all variables. 
 	void ResetVariables();
 	
 	/// \brief ROOT TFile that will be used to store the output from this tool
  TFile* fOutput_tfile = nullptr;
  
  /// \brief TTree that will be used to store output
  TTree* fRecoTree = nullptr;
  
  /// \brief Branch variables
  G4int fEventNumber;
  /// primary gammas
  std::vector<G4double> fGamma_x;
  std::vector<G4double> fGamma_y;
  std::vector<G4double> fGamma_z;
  std::vector<G4double> fGamma_t;
  std::vector<G4double> fGamma_energy;
  	
  /// created electrons
  std::vector<G4double> fElectron_x;
  std::vector<G4double> fElectron_y;
  std::vector<G4double> fElectron_z;
  std::vector<G4double> fElectron_t;
  std::vector<G4double> fElectron_dx;
  std::vector<G4double> fElectron_dy;
  std::vector<G4double> fElectron_dz;
  std::vector<G4double> fElectron_energy;
  std::vector<G4double> fElectron_clusterID;
  std::vector<G4int> fElectron_nsecondaries_last_step;
  std::vector<G4double> fElectron_dE_last_step;
  std::vector<G4double> fElectron_weight;
  std::vector<G4double> fElectron_weight_modified;
  	
  // clusters
  G4int fNclusters;
  std::vector<G4int> fNelectron_in_cluster;
  std::vector<G4double> fClusteredE_x;
  std::vector<G4double> fClusteredE_y;
  std::vector<G4double> fClusteredE_z;
  	
  // event
  G4int fNElectrons;
  
  // functions
  void Initialize(std::string output_filename);
  void Fill() {fRecoTree->Fill();}
  void CloseFile() {fRecoTree->Write(); fOutput_tfile->Close();}
  
  void AddPrimaryGamma(G4double x, G4double y, G4double z, G4double t, G4double e);
  void AddElectron(G4double x, G4double y, G4double z, G4double t, G4double e, G4int cluster, G4double weight1, G4double weight2, G4int Nsecondary);
  void SetNumberOfElectrons(G4int n) {fNElectrons = n;}
  void AddNelectronsInCluster(G4int n) {fNelectron_in_cluster.push_back(n);}
  void SetNumberOfClusters(G4int n) {fNclusters = n;}
  void AddEventNumber() {fEventNumber++;}
  
  // 
  
  
  // Reset
  void Reset();

 private:


};

#endif







