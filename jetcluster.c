#include "TFile.h"
#include "vector"
#include <math.h>
#include <algorithm>
#include <iostream>
using namespace std;

struct p4info{
  double pT;
  double eta;
  double phi;
  double mass;
};


double dij(p4info particle1, p4info particle2){//calculate d_ij defined in the anti-kT algorithm
  return std::min(1/pow(particle1.pT,2),1/pow(particle2.pT,2)) * sqrt(pow(particle1.phi-particle2.phi,2) + pow(particle1.eta-particle2.eta, 2));
}

p4info merge(p4info particle1, p4info particle2){//combine particle1 and particle2 into one particle
  p4info particle12;

  TLorentzVector v, v1, v2;
  v1.SetPtEtaPhiM(particle1.pT, particle1.eta, particle1.phi, particle1.mass);
  v2.SetPtEtaPhiM(particle2.pT, particle2.eta, particle2.phi, particle2.mass);

  v = v1 + v2;
  particle12.pT = v.Pt();
  particle12.eta = v.Eta();
  particle12.phi = v.Phi();
  particle12.mass = v.M();  
  
  return particle12;
}


void jetcluster(){
  
  TFile *f= new TFile("for_jet_clustering_exercise.root");
  TTree *t1 = (TTree*)f->Get("tree");
  
  vector<float> *pt = 0;// pointer must be initialized before passed as branch address
  vector<float> *eta = 0;
  vector<float> *phi = 0;
  vector<float> *mass = 0;

  t1->SetBranchAddress("pt",   &pt);
  t1->SetBranchAddress("eta",  &eta);
  t1->SetBranchAddress("phi",  &phi);
  t1->SetBranchAddress("mass", &mass);

  //new root file, store info 
  TFile *fnew = new TFile("jetresult.root","RECREATE");
  TTree t2("t2", "new tree");

  double jetpt, jeteta, jetphi, jetmass;
  t2.Branch("jetpt",   &jetpt,   "jetpt/D");
  t2.Branch("jeteta",  &jeteta,  "jeteta/D");
  t2.Branch("jetphi",  &jetphi,  "jetphi/D");
  t2.Branch("jetmass", &jetmass, "jetmass/D");

  TH1D *histo1 = new TH1D("histo1", "number of jets with pT > 20GeV", 10, 0, 10);

  Int_t nevents = (Int_t)t1->GetEntries();
 
  for(Int_t i=0; i<nevents; i++){
     t1->GetEntry(i);   // i.th event
     int nparticles = (*pt).size(); // number of particles in i.th event 
     int num = 0; //number of jets with pT > 20GeV

     double d2[nparticles][nparticles];
     double d1[nparticles];
     p4info particle[nparticles];
     bool single[nparticles];

     //--------------**-----------------
     //initialize
     //-------------**------------------

     for(Int_t j=0; j<nparticles; j++){
        particle[j].pT = (*pt)[j];
        particle[j].eta = (*eta)[j];
        particle[j].phi = (*phi)[j];
        particle[j].mass = (*mass)[j];

        single[nparticles] = true;
     }

     for(Int_t j=0; j<nparticles; j++)
        d1[j] = 1/pow(particle[j].pT,2);

     for(Int_t j=0; j<nparticles-1; j++)
        for(Int_t k=j+1; k<nparticles; k++)
           d2[j][k] = dij(particle[j], particle[k]);
 
     //-------------**------------------
     // anti-kT algorithm
     //-------------**------------------

     int nleft = nparticles;
     while(nleft){
        double d2min = std::numeric_limits<double>::max();
	double d1min = std::numeric_limits<double>::max();
	int i1, i2, i3;

        for(Int_t j=0; j<nparticles-1; j++){ // find the minimum of dij
           if(!single[j]) continue;
           for(Int_t k=j+1; k<nparticles; k++){
              if(!single[k]) continue;
              if(d2min > d2[j][k]){
                 i1 = j;
                 i2 = k;
                 d2min = d2[j][k];
              }
           } 
        }
 
        for(Int_t j=0; j<nparticles; j++){ // find the minimum of di
           if(!single[j]) continue;
           if(d1min > d1[j]){
              i3 = j;
              d1min = d1[j];
           }
        }

        if(d2min < d1min){//combine i1 and i2, eliminate i2
           particle[i1] = merge(particle[i1], particle[i2]);
	   single[i2] = false;
	   nleft--;

	   d1[i1] = 1/pow(particle[i1].pT,2);

           for(Int_t k=0; k<i1; k++){
              if(!single[k]) continue; 
              d2[k][i1] = dij(particle[k], particle[i1]);
           
           }
           for(Int_t k=i1+1; k<nparticles; k++){
              if(!single[k]) continue;
	      d2[i1][k] = dij(particle[i1], particle[k]);
	   }   
        }

	else { //one jet found
	   single[i3] = false;
	   nleft--;

	   jetpt = particle[i3].pT;
	   jeteta = particle[i3].eta;
	   jetphi = particle[i3].phi;
	   jetmass = particle[i3].mass;

	   t2.Fill();

	   if(jetpt > 20) num++;
	}

     }//finish one event

     histo1->Fill(num);

  }//finish

  fnew->Write();   
}
