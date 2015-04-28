#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
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
#include "RooUnfold.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldInvert.h"
#include "RooUnfoldSvd.h"
#include <sstream>
#include <string>

#include <iostream>
#include <fstream>

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_reco.root";
std::string File_Powheg_reco="/afs/cern.ch/work/r/ruckstuh/public/Powheg_reco.root";
std::string File_Signal_gen="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_gen.root";
std::string File_Powheg_gen="/afs/cern.ch/work/r/ruckstuh/public/Powheg_gen.root";

std::string reco_name="Combined Single Reco";
std::string gen_name ="Combined Gen Cuts Reco";

TH1D* GetPhistar(bool madgraph, bool gen, int elec=0){
  std::string name=reco_name;
  if (gen) name=gen_name;
  TChain* t = new TChain(name.c_str(),name.c_str());
  int nfiles;
  if (!gen){
    if (madgraph)  nfiles=t->Add(File_Signal_reco.c_str());
    else nfiles=t->Add(File_Powheg_reco.c_str());
  }
  else{
    if (madgraph)  nfiles=t->Add(File_Signal_gen.c_str());
    else nfiles=t->Add(File_Powheg_gen.c_str());
  }
  TBranch *b_truth=t->GetBranch("truth");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  if (elec==1) l_phistar_true=b_truth->GetLeaf("z_phistar_born");
  if (elec==2) l_phistar_true=b_truth->GetLeaf("z_phistar_naked");
  TLeaf *l_y=b_truth->GetLeaf("z_y");
  int nweights;
  t->SetBranchAddress("weight_size",&nweights);
  t->GetEntry(0);
  cout<<"The sample has nweights: "<<nweights<<endl;
  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);

  TH1D* ph=new TH1D("phistar","phistar",nphistar,phistarBins);
  ph->Sumw2();
  for (int i=0; i<t->GetEntries();i++){
    t->GetEntry(i);
    //if (fabs(l_y->GetValue())<1.5)continue;
    double phistar_true=l_phistar_true->GetValue();
    double weight =1;
    for (int w=0; w<nweights;w++){
      if (!gen){
	if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
      }
      else{
	if (weightid[w]==1 || weightid[w]==2) {weight=weight*weights[w];}
      }
    }
    ph->Fill(phistar_true,weight);
  }
  return ph;
}

void NTupleEventEff(){
  for (int elec=0; elec<1; elec++){
    TH1D *MG_reco= GetPhistar(1, 0, elec);
    TH1D *MG_gen = GetPhistar(1, 1, elec);
    TH1D *PH_reco= GetPhistar(0, 0, elec);
    TH1D *PH_gen = GetPhistar(0, 1, elec);

    cout<<"done reading"<<endl; 
    TH1D *MG_eff=(TH1D*)MG_reco->Clone();
    TH1D *PH_eff=(TH1D*)PH_reco->Clone();
    cout<<"done cloning"<<endl; 
    MG_eff->Divide(MG_reco,MG_gen,1.,1.,"B");
    PH_eff->Divide(PH_reco,PH_gen,1.,1.,"B");
    cout<<"done dividing"<<endl; 

    TCanvas* Efficiency = new TCanvas("Efficiency","Efficiency",800,900);
    Efficiency->cd();
    Efficiency->SetLogx();
    MG_eff->GetXaxis()->SetRangeUser(0.001,3.2);
    MG_eff->GetXaxis()->SetTitle("#phi*_{generated}");
    MG_eff->GetXaxis()->SetTitleOffset(0.8);
    MG_eff->GetXaxis()->SetTitleSize(0.04);
    MG_eff->GetXaxis()->SetLabelOffset(-0.01);
    MG_eff->GetXaxis()->SetLabelSize(0.04);
    MG_eff->GetYaxis()->SetTitle("N_{reconstructed}/N_{generated}");
    MG_eff->GetYaxis()->SetTitleOffset(1.2);
    MG_eff->GetYaxis()->SetTitleSize(0.04);
    MG_eff->GetYaxis()->SetLabelSize(0.04);
    //    MG_eff->GetYaxis()->SetRangeUser(0.5,0.57);
    MG_eff->SetStats(0);
    MG_eff->SetBit( TH1::kNoTitle, true );
    MG_eff->SetLineColor(1);
    MG_eff->SetMarkerColor(1);
    MG_eff->SetMarkerStyle(20);
    MG_eff->Draw();
    PH_eff->SetMarkerStyle(21);
    PH_eff->SetLineColor(2);
    PH_eff->SetMarkerColor(2);
    PH_eff->Draw("same");
    MG_eff->Draw("same");
    TLegend* leg_eff = new TLegend(0.45,0.77,0.85,0.91);
    leg_eff->SetFillStyle(0);
    leg_eff->SetBorderSize(0);
    leg_eff->SetLineWidth(1);
    leg_eff->SetNColumns(1);
    leg_eff->SetTextFont(42);
    leg_eff->SetTextSize(0.04);
    leg_eff->AddEntry(MG_eff,"MadGraph","PL");
    leg_eff->AddEntry(PH_eff,"Powheg","PL");
    leg_eff->Draw();
    return; 
  }
}
