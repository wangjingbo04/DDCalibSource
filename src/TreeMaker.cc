#include "TreeMaker.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

static TreeMaker* fgTreeMaker = 0;

TreeMaker* TreeMaker::Instance()
{
  if( !fgTreeMaker ){
    fgTreeMaker = new TreeMaker();
  }

  return fgTreeMaker;
}

TreeMaker::TreeMaker() {
	
}

TreeMaker::~TreeMaker() {
}

void TreeMaker::AddElectron(G4double x, G4double y, G4double z, G4double t, G4double e, G4int cluster, G4double weight1, G4double weight2, G4int Nsecondary) {
  fElectron_x.push_back(x);
  fElectron_y.push_back(y);   
  fElectron_z.push_back(z);   
  fElectron_t.push_back(t);   
  fElectron_energy.push_back(e);	
  fElectron_clusterID.push_back(cluster);
  fElectron_weight.push_back(weight1);
  fElectron_weight_modified.push_back(weight2);
  fElectron_nsecondaries_last_step.push_back(Nsecondary);
}

void TreeMaker::AddPrimaryGamma(G4double x, G4double y, G4double z, G4double t, G4double e) {
  fGamma_x.push_back(x);
  fGamma_y.push_back(y);
  fGamma_z.push_back(z);
  fGamma_t.push_back(t);
  fGamma_energy.push_back(e);	
}

void TreeMaker::Reset() {
  fGamma_x.clear();         
  fGamma_y.clear();         
  fGamma_z.clear();         
  fGamma_t.clear();         
  fGamma_energy.clear();                                          
  fElectron_x.clear();      
  fElectron_y.clear();      
  fElectron_z.clear();      
  fElectron_t.clear();      
  fElectron_dx.clear();     
  fElectron_dy.clear();     
  fElectron_dz.clear();     
  fElectron_energy.clear(); 
  fElectron_clusterID.clear();
  fElectron_nsecondaries_last_step.clear();
  fElectron_dE_last_step.clear();
  fElectron_weight.clear();
  fElectron_weight_modified.clear();
  fNelectron_in_cluster.clear();
	
}

void TreeMaker::Initialize(std::string output_filename) {
	fOutput_tfile = new TFile(output_filename.c_str(), "recreate");
	fRecoTree = new TTree("gamma", "Neutron capture gamma cascade");
	fRecoTree->Branch("eventNo", &fEventNumber, "fEventNumber/I");
	fRecoTree->Branch("gamma_x",&fGamma_x);
	fRecoTree->Branch("gamma_y",&fGamma_y);
	fRecoTree->Branch("gamma_z",&fGamma_z);
	fRecoTree->Branch("gamma_t",&fGamma_t);
	fRecoTree->Branch("gamma_energy",&fGamma_energy);
	fRecoTree->Branch("electron_x",&fElectron_x);
	fRecoTree->Branch("electron_y",&fElectron_y);
	fRecoTree->Branch("electron_z",&fElectron_z);
	fRecoTree->Branch("electron_t",&fElectron_t);
	fRecoTree->Branch("electron_dx",&fElectron_dx);
	fRecoTree->Branch("electron_dy",&fElectron_dy);
	fRecoTree->Branch("electron_dz",&fElectron_dz);
	fRecoTree->Branch("electron_energy",&fElectron_energy);
	fRecoTree->Branch("electron_clusterID",&fElectron_clusterID);
	fRecoTree->Branch("electron_nsecondary",&fElectron_nsecondaries_last_step);
	fRecoTree->Branch("electron_dE_last",&fElectron_dE_last_step);
	fRecoTree->Branch("electron_weight",&fElectron_weight);
	
	fRecoTree->Branch("electron_weight_modified",&fElectron_weight_modified);
	fRecoTree->Branch("Nelectron_in_cluster", &fNelectron_in_cluster);
  fRecoTree->Branch("Nclusters",&fNclusters, "fNclusters/I");
  fRecoTree->Branch("Nelectrons",&fNElectrons, "fElectrons/I");

  
}

