
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
    RecoElectron(G4double x, G4double y, G4double z, G4double t, G4int nsec): fX(x), fY(y), fZ(z), fT(t), fNSecondary(nsec){}
   ~RecoElectron();
   G4double GetX() { return fX; }
   G4double GetY() { return fY; }
   G4double GetZ() { return fZ; }
   G4double GetTime() {return fT;}
   G4int GetNumberOfSecondaries() {return fNSecondary;}

  public:
  	
    void SetValues(G4double x, G4double y, G4double z, G4double t);
               
  private:                  
  	G4double  fX;
    G4double  fY;
    G4double  fZ;
    G4double  fT; 
    G4int			fNSecondary;
          
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
