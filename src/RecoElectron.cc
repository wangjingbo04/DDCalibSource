
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RecoElectron.hh"

#include "Run.hh"
#include "HistoManager.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RecoElectron::RecoElectron(): fX(0), fY(0.), fZ(0.), fT(0.) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RecoElectron::~RecoElectron(){}
	
void RecoElectron::SetValues(G4double x, G4double y, G4double z, G4double t) {
       fX = x;
       fY = y;
       fZ = z;
       fT = t;	
} 



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


