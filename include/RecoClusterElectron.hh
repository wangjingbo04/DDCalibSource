
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RecoClusterElectron_h
#define RecoClusterElectron_h 1

#include "G4Track.hh"
#include "globals.hh"
#include "Run.hh"
#include "RecoElectron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RecoClusterElectron
{
  public:
    RecoClusterElectron(RecoElectron* recoElectron);
   ~RecoClusterElectron();
   void AddClusterElectron(RecoClusterElectron* clusterElectron);
   G4int GetNClusterElectrons(); 
   RecoClusterElectron* GetClusterElectron(G4int idigit); 
   std::vector<RecoClusterElectron*>* GetClusterElectronList();
   G4double GetX(){ return fRecoElectron->GetX(); }
   G4double GetY(){ return fRecoElectron->GetY(); }
   G4double GetZ(){ return fRecoElectron->GetZ(); }
   G4double GetTime(){ return fRecoElectron->GetTime(); }
   G4int SetClusterID(G4int id) {fClusterID = id;}

  public:
  	
    void SetValues(G4double x, G4double y, G4double z, G4double t);
    
    void SetClustered( G4bool yesno = 1 ){ fIsClustered = yesno; }
    G4bool IsClustered(){ return fIsClustered; }
    G4bool IsAllClustered();
    RecoElectron* GetRecoElectron(){ return fRecoElectron; }
               
  private:                  
  	G4bool fIsClustered;
    G4bool fIsAllClustered;
    G4int			fClusterID = -1;

    RecoElectron* fRecoElectron;

    std::vector<RecoClusterElectron*>* fClusterElectronList;
          
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
