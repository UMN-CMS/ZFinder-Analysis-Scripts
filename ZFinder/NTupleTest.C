#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TChain.h"

#include "TMatrixDSparse.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TArray.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <sstream>
#include <string>

#include <iostream>
#include <fstream>
void NTupleTest(const std::string infile, const std::string cteqfile, const bool reco){
  TFile f(infile.c_str(), "update");

  Float_t cteq6;
  std::string reco_name="Combined Single Reco";
  if (!reco) reco_name="Combined Gen Cuts Reco";
  
  TTree *t = (TTree*)f.Get(reco_name.c_str());
  TBranch *newBranch = t->Branch("weight_cteq6", &cteq6, "weight_cteq6/F");
  
  Long64_t nentries = t->GetEntries();
  TBranch *b_event=t->GetBranch("event_info");
  TLeaf *l_en=b_event->GetLeaf("event_number");
  int nwcteq;
  t->SetBranchAddress("weight_cteq_size",&nwcteq);
  t->GetEntry(0);
  double weights_cteq[nwcteq];
  t->SetBranchAddress("weights_cteq",&weights_cteq);
  
  TChain* t2 = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles2;
  nfiles2=t2->Add(cteqfile.c_str());
  t2->SetBranchStatus("reco",0); //to disable all branches
  t2->SetBranchStatus("truth",0); //to disable all branches
  TBranch *b_event2=t2->GetBranch("event_info");
  TLeaf *l_en2=b_event2->GetLeaf("event_number");
  double weights_cteq6[1];
  t2->SetBranchAddress("weights_cteq",&weights_cteq6);
  
  vector<int > eventnvec;
  vector<double > wcteqvec;
  int nentries2= t2->GetEntries();
  cout<<nentries2<<endl;
  for (Long64_t i=0; i < nentries2; i++){
    if (i%10000==0)cout<<i<<endl;
    t2->GetEntry(i); 
    eventnvec.push_back(l_en2->GetValue());
    wcteqvec.push_back(weights_cteq6[0]);
    if (i%10000==0)cout<<eventnvec[i]<<endl;
  }
  int idx=0;
  int event_not_found=0;
  for (Long64_t i = 0; i < nentries; i++){
    bool found =0;
    if (i%10000==0)cout<<i<<" "<<idx<<endl;
    // cout<<i<<endl;
    t->GetEntry(i);
    cteq6= weights_cteq[0];
    double eventnumber=l_en->GetValue();
    if (idx<nentries2 && eventnvec[idx]==eventnumber){
      cteq6=wcteqvec[idx];
      found=1;
      idx++;
    }
    else {
      idx++;
      if (idx<nentries2 && eventnvec[idx]==eventnumber){
	cteq6=wcteqvec[idx];
	found=1;
	idx++;
      }
    }
    if (!found){
      for (Long64_t j = 0; j < nentries2; j++){
	if (eventnvec[j]==eventnumber){
	  cteq6=wcteqvec[6];
	  found=1;
	  cout<<i<<" "<<j<<" "<<idx<<endl;
	  idx=j+1;
	  continue;
	}
      }
    }
    if (!found) {event_not_found++; idx=idx-1;}
    newBranch->Fill();
  }
  cout<<event_not_found<<endl;
  // save only the new version of the tree
  t->Write("", TObject::kOverwrite);
}

