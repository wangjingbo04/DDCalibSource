
#ifndef RECOCLUSTER_HH
#define RECOCLUSTER_HH

#include "G4Track.hh"
#include "globals.hh"
#include "Run.hh"
#include "RecoElectron.hh"


class RecoCluster {

 public:
  RecoCluster();
  ~RecoCluster();

  void Reset();
  void SortCluster();

  void AddElectron(RecoElectron* digit);

  RecoElectron* GetElectron(G4int n);
  G4int GetNElectrons();

 private:

  std::vector<RecoElectron*> fElectronList;

};

#endif







