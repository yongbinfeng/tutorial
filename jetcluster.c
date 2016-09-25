#include "TFile.h"
#include "vector"

int antikT(vector<float> *pt, vector<float> *eta, vector<float>*phi, vector<float> *mass){
    return 1;
}

void jetcluster(){
  
  TFile *f= new TFile("for_jet_clustering_exercise.root");
  TTree *t1 = (TTree*)f->Get("tree");
  
  vector<float> *pt = 0;// pointer muast be initialized before passed as branch address
  vector<float> *eta = 0;
  vector<float> *phi = 0;
  vector<float> *mass = 0;

  t1->SetBranchAddress("pt",   &pt);
  t1->SetBranchAddress("eta",  &eta);
  t1->SetBranchAddress("phi",  &phi);
  t1->SetBranchAddress("mass", &mass);

  Int_t nentries = (Int_t)t1->GetEntries();
  cout<<nentries<<endl;
  for(Int_t i=0; i<nentries; i++){
     t1->GetEntry(i);   // i.th event
      
     
  

  }
  
   

}
