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

std::string reco_name="Combined Single Reco";

void GetBinM(RooUnfoldResponse* &BinM, TH1D* &h_gen, TH1D* &h_reco, bool madgraph=1, int elec=0){
  TH2D* BinMigration=new TH2D("BinMigration","BinMigration",nphistar,phistarBins,nphistar,phistarBins);
  BinMigration->Sumw2();
  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  if (madgraph)  nfiles=t->Add(File_Signal_reco.c_str());
  else nfiles=t->Add(File_Powheg_reco.c_str());

  TBranch *b_reco=t->GetBranch("reco");
  TBranch *b_truth=t->GetBranch("truth");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  if (elec==1) l_phistar_true=b_truth->GetLeaf("z_phistar_born");
  if (elec==2) l_phistar_true=b_truth->GetLeaf("z_phistar_naked");
  int nweights;
  t->SetBranchAddress("weight_size",&nweights);
  t->GetEntry(0);
  cout<<"The sample has nweights: "<<nweights<<endl;
  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);

  cout<<"Entries: "<<t->GetEntries()<<endl;
  h_gen   = new TH1D("phistar","phistar",nphistar,phistarBins);
  h_gen->Sumw2();
  h_reco  = new TH1D("phistar","phistar",nphistar,phistarBins);
  h_reco->Sumw2();
  BinM  =new RooUnfoldResponse (h_reco,h_gen);

  TH1D *h_w=new TH1D("h_w","h_w",10000,0,10);

  if (madgraph){
    for (int i=0; i<t->GetEntries();i++){
      t->GetEntry(i);
      double phistar_true=l_phistar_true->GetValue();
      double weight =1;
      for (int w=0; w<nweights;w++){
	if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
      }
      h_w ->Fill(phistar_true,weight);
    }
  }

  cout<<"reading signal reco "<<endl;
  for (int i=0; i<t->GetEntries();i++){
    // if (i>20000) continue;
    //for (int i=0; i<50000;i++){
    t->GetEntry(i);
    double phistar=l_phistar->GetValue();
    double phistar_true=l_phistar_true->GetValue();
    double weight =1;
    for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
    }
    if (madgraph){
      // cout<<i<<" "<<phistar_true<<endl;
      int idx=h_w->GetXaxis()->FindBin(phistar_true);
      // cout<<idx<<" "<<h_w->GetBinContent(idx)<<endl;
      if (h_w->GetBinContent(idx)==0) {
	cout<<"error bin content is 0:"<<weight<<endl;
	weight=0;
      }
      else weight=10000*weight/h_w->GetBinContent(idx);
      // cout<<weight<<endl;
    }
    h_gen ->Fill(phistar_true,weight);
    h_reco->Fill(phistar,weight);
    BinM  ->Fill(phistar,phistar_true,weight); 
    BinMigration->Fill(phistar,phistar_true,weight);
  }
 
  TLatex mark;
  mark.SetTextSize(0.04);
  mark.SetTextFont(42);
  mark.SetNDC(true);

  BinMigration->GetXaxis()->SetRangeUser(0.001,3.2);
  BinMigration->GetXaxis()->SetTitle("#phi*(reconstructed)");
  BinMigration->GetXaxis()->SetTitleOffset(0.8);
  BinMigration->GetXaxis()->SetTitleSize(0.04);
  BinMigration->GetXaxis()->SetLabelOffset(-0.01);
  BinMigration->GetXaxis()->SetLabelSize(0.04);
  BinMigration->GetYaxis()->SetTitleOffset(1.05);
  BinMigration->GetYaxis()->SetTitleSize(0.04);
  BinMigration->GetYaxis()->SetLabelSize(0.04);
  BinMigration->GetYaxis()->SetRangeUser(0.001,3.2);
  BinMigration->GetYaxis()->SetTitle("#phi*(generated)");
  BinMigration->GetZaxis()->SetTitleOffset(-0.004);
  BinMigration->SetStats(0);
  BinMigration->SetBit( TH2::kNoTitle, true );

  if (madgraph){
    TCanvas* BinM_M = new TCanvas("BinM_M","BinM_M",800,900);
    BinM_M->cd();
    BinM_M->SetLogx();
    BinM_M->SetLogy();
    BinM_M->SetLogz();
    BinMigration->Draw("COLZ");
    mark.DrawLatex(0.19,0.80,"reweighed MadGraph");
  }
  return;
}

void NTupleBinM(){
  for (int elec=0; elec<1; elec++){
    RooUnfoldResponse *m_BinM;
    TH1D *m_gen, *m_reco;
    GetBinM(m_BinM, m_gen, m_reco, 1, elec);

    RooUnfoldResponse *p_BinM;
    TH1D *p_gen, *p_reco;
    GetBinM(p_BinM, p_gen, p_reco, 0, elec);

    for (uint i=0; i<nphistar;i++){
      m_gen->SetBinError(i+1,0);
      p_gen->SetBinError(i+1,0);
    }

    RooUnfoldBayes  unfold_MP(m_BinM,p_reco, 4);
    TH1D* h_BinM_MP  = (TH1D*) unfold_MP.Hreco();  
    TH1D* h_BinM_MP_t= (TH1D*) unfold_MP.Hreco();
    h_BinM_MP->Divide(h_BinM_MP_t,p_gen);
 
    TLatex mark;
    mark.SetTextSize(0.04);
    mark.SetTextFont(42);
    mark.SetNDC(true);

    TCanvas* BinMigration_MP = new TCanvas("BinMigration_MP","BinMigration_MP",800,900);
    BinMigration_MP->cd();
    BinMigration_MP->SetLogx();
    h_BinM_MP->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_MP->GetXaxis()->SetTitle("#phi*");
    h_BinM_MP->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_MP->GetXaxis()->SetTitleSize(0.04);
    h_BinM_MP->GetXaxis()->SetLabelOffset(-0.01);
    h_BinM_MP->GetXaxis()->SetLabelSize(0.04);
    h_BinM_MP->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_MP->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_MP->GetYaxis()->SetTitleSize(0.04);
    h_BinM_MP->GetYaxis()->SetLabelSize(0.03);
    h_BinM_MP->GetYaxis()->SetRangeUser(0.9,1.1);
    h_BinM_MP->SetMarkerStyle(20);
    h_BinM_MP->SetStats(0);
    h_BinM_MP->SetBit( TH1::kNoTitle, true );
    h_BinM_MP->SetLineColor(1);
    h_BinM_MP->Draw();
    TLegend* leg_MP = new TLegend(0.45,0.80,0.85,0.86);
    leg_MP->SetFillStyle(0);
    leg_MP->SetBorderSize(0);
    leg_MP->SetLineWidth(1);
    leg_MP->SetNColumns(1);
    leg_MP->SetTextFont(42);
    leg_MP->SetTextSize(0.04);
    leg_MP->AddEntry(h_BinM_MP,"Powheg","PL");
    leg_MP->Draw();
    mark.DrawLatex(0.19,0.17,"Unfolded using reweighed MadGraph");
   }
}
