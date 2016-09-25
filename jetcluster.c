#include "TFile.h"
#include "vector"


void jetcluster(){
  
  TFile *f= new TFile("/home/cmsumd/Desktop/exercise/Jetcluster/for_jet_clustering_exercise.root");
  TTree *tree = (TTree*)f->Get("tree");
  
  vector<float> *pt, *eta, *phi, *mass;
  tree->SetBranchAddress("pt",   &pt);
  tree->SetBranchAddress("eta",  &eta);
  tree->SetBranchAddress("phi",  &phi);
  tree->SetBranchAddress("mass", &mass);

  Int_t nentries = (Int_t)tree->GetEntries();
  for(Int_t i=0; i<nentries; i++){
     tree->GetEntry(i);   // i.th event

      
     
  

  }
  
   

}
