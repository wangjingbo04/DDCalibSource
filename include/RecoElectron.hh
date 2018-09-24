
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef RecoElectron_h
#define RecoElectron_h 1

#include "G4Track.hh"
#include "globals.hh"
#include "Run.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class RecoElectron
{
  public:
    RecoElectron();
    RecoElectron(G4double x, G4double y, G4double z, G4double t, G4double e, G4int nsec): fX(x), fY(y), fZ(z), fT(t),fEnergy(e), fNSecondary(nsec){}
   ~RecoElectron();
   G4double GetX() { return fX; }
   G4double GetY() { return fY; }
   G4double GetZ() { return fZ; }
   G4double GetTime() {return fT;}
   G4double GetEnergy() {return fEnergy;}
   G4int GetNumberOfSecondaries() {return fNSecondary;}
   G4int SetClusterID(G4int id) {fClusterID = id;}
   G4int GetClusterID() {return fClusterID;}

  public:
  	
    void SetValues(G4double x, G4double y, G4double z, G4double t, G4double e);
               
  private:                  
  	G4double  fX;
    G4double  fY;
    G4double  fZ;
    G4double  fT; 
    G4double	fEnergy;
    G4int			fNSecondary;
    G4int			fClusterID = -1;
          
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
