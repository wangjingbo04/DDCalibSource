#include "RecoCluster.hh"
#include "RecoElectron.hh"

RecoCluster::RecoCluster()
{
}

RecoCluster::~RecoCluster()
{
  Reset();
}

void RecoCluster::Reset()
{
  int j=fElectronList.size();
  for(int i=0;i<j;i++){
    delete fElectronList[i];
  }
  fElectronList.clear();
}

static bool CompareTimes(RecoElectron *rd1, RecoElectron *rd2)
{
  return ( rd1->GetTime() > rd2->GetTime() );
}

void RecoCluster::SortCluster()
{
  sort(fElectronList.begin(), fElectronList.end(), CompareTimes);
}

void RecoCluster::AddElectron(RecoElectron* electron)
{
  fElectronList.push_back(electron);
}

RecoElectron* RecoCluster::GetElectron(G4int n)
{
  return (RecoElectron*)(fElectronList.at(n));
}

G4int RecoCluster::GetNElectrons()
{
  return fElectronList.size();
}
