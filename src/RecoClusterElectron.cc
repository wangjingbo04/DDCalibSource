
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RecoClusterElectron.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RecoClusterElectron::RecoClusterElectron(RecoElectron* recoElectron)
{
  fIsClustered = 0;
  fIsAllClustered = 0;

  fRecoElectron = recoElectron;

  fClusterElectronList = new std::vector<RecoClusterElectron*>;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RecoClusterElectron::~RecoClusterElectron(){
  delete fClusterElectronList;	
}
	
void RecoClusterElectron::AddClusterElectron(RecoClusterElectron* clusterElectron)
{ 
  fClusterElectronList->push_back(clusterElectron); 
}

G4int RecoClusterElectron::GetNClusterElectrons() 
{ 
  return fClusterElectronList->size(); 
}
  
RecoClusterElectron* RecoClusterElectron::GetClusterElectron(G4int idigit) 
{ 
  return (RecoClusterElectron*)(fClusterElectronList->at(idigit)); 
}
  
std::vector<RecoClusterElectron*>* RecoClusterElectron::GetClusterElectronList()
{ 
  return fClusterElectronList; 
}

G4bool RecoClusterElectron::IsAllClustered()
{
  if( fIsAllClustered==0 ){
    fIsAllClustered = 1;
    for( G4int n=0; n<fClusterElectronList->size(); n++ ){
      RecoClusterElectron* clusterElectron = (RecoClusterElectron*)(fClusterElectronList->at(n));
      if( clusterElectron->IsClustered()==0 ) fIsAllClustered = 0;
    }
  }

  return fIsAllClustered;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


