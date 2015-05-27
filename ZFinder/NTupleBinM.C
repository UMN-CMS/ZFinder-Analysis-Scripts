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

int N_MC=-1;

void GetBinM(RooUnfoldResponse* &BinM, TH1D* &h_gen, TH1D* &h_reco, RooUnfoldResponse* &BinM_1, TH1D* &h_gen_1, TH1D* &h_reco_1, RooUnfoldResponse* &BinM_2, TH1D* &h_gen_2, TH1D* &h_reco_2, bool madgraph=1, int elec=0){
  TH2D* BinMigration=new TH2D("BinMigration","BinMigration",nphistar,phistarBins,nphistar,phistarBins);
  BinMigration->Sumw2();
  TH1D* Dphistar=new TH1D("Dphistar","Dphistar",2000,-1,1);
  Dphistar->Sumw2();
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
  h_gen_1 = new TH1D("phistar","phistar",nphistar,phistarBins);
  h_gen_1->Sumw2();
  h_reco_1= new TH1D("phistar","phistar",nphistar,phistarBins);
  h_reco_1->Sumw2();
  h_gen_2 = new TH1D("phistar","phistar",nphistar,phistarBins);
  h_gen_2->Sumw2();
  h_reco_2= new TH1D("phistar","phistar",nphistar,phistarBins);
  h_reco_2->Sumw2();
  BinM  =new RooUnfoldResponse (h_reco,h_gen);
  BinM_1=new RooUnfoldResponse (h_reco,h_gen);
  BinM_2=new RooUnfoldResponse (h_reco,h_gen);

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
    h_gen ->Fill(phistar_true,weight);
    h_reco->Fill(phistar,weight);
    if (N_MC==-1 || i<N_MC) BinM  ->Fill(phistar,phistar_true,weight); 
    BinMigration->Fill(phistar,phistar_true,weight);
    Dphistar->Fill((phistar_true-phistar)/phistar_true,weight);
    if (i%2==0){
      h_gen_1 ->Fill(phistar_true,weight);
      h_reco_1->Fill(phistar,weight);
      if (N_MC==-1 || i<N_MC) BinM_1  ->Fill(phistar,phistar_true,weight); 
    }
    else{
      h_gen_2 ->Fill(phistar_true,weight);
      h_reco_2->Fill(phistar,weight);
      if (N_MC==-1 || i<N_MC) BinM_2  ->Fill(phistar,phistar_true,weight); 
    }
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
  BinMigration->GetYaxis()->SetTitle("#phi* (generated)");
  BinMigration->GetZaxis()->SetTitleOffset(-0.004);
  BinMigration->SetStats(0);
  BinMigration->SetBit( TH2::kNoTitle, true );

  Dphistar->GetXaxis()->SetTitle("(#phi*(generated) - #phi*(reconstructed)/#phi*(generated)");
  Dphistar->GetXaxis()->SetTitleOffset(1.0);
  Dphistar->GetXaxis()->SetTitleSize(0.04);
  Dphistar->GetXaxis()->SetLabelOffset(0.0);
  Dphistar->GetXaxis()->SetLabelSize(0.04);
  Dphistar->GetYaxis()->SetTitleOffset(1.05);
  Dphistar->GetYaxis()->SetTitleSize(0.04);
  Dphistar->GetYaxis()->SetLabelSize(0.04);
  Dphistar->GetYaxis()->SetTitle("a.u.");
  Dphistar->SetLineColor(1);
  Dphistar->SetStats(0);
  Dphistar->SetBit( TH2::kNoTitle, true );

  if (madgraph){
    TCanvas* BinM_M = new TCanvas("BinM_M","BinM_M",800,900);
    BinM_M->cd();
    BinM_M->SetLogx();
    BinM_M->SetLogy();
    BinM_M->SetLogz();
    BinMigration->Draw("COLZ");
    mark.DrawLatex(0.19,0.80,"MadGraph");
    TCanvas* Dp_M = new TCanvas("Dp_M","Dp_M",800,900);
    Dp_M->cd();
    Dp_M->SetLogy();
    Dphistar->Draw();
    mark.DrawLatex(0.19,0.80,"MadGraph");
  }
  else{
    TCanvas* BinM_P = new TCanvas("BinM_P","BinM_P",800,900);
    BinM_P->cd();
    BinM_P->SetLogx();
    BinM_P->SetLogy();
    BinM_P->SetLogz();
    BinMigration->Draw("COLZ");
    mark.DrawLatex(0.19,0.80,"Powheg");
    TCanvas* Dp_P = new TCanvas("Dp_P","Dp_P",800,900);
    Dp_P->cd();
    Dp_P->SetLogy();
    Dphistar->Draw();
    mark.DrawLatex(0.19,0.80,"Powheg");
  }

  double weighted_events= BinMigration->GetSumOfWeights();  
  double diagonal=0;
  for (int i=0; i<nphistar; i++){
    diagonal+=BinMigration->GetBinContent(i+1,i+1);
  }
  cout<<weighted_events<<" of which diagonal "<<diagonal<<" non diagonal faction:"<<(weighted_events-diagonal)/weighted_events<<endl;
  cout<<"done reading data for "<<File_Signal_reco<<"  "<<reco_name<<endl;

  for (int i=0; i<nphistar;i++){
    h_gen->SetBinError(i+1,0);
    h_gen_1->SetBinError(i+1,0);
    h_gen_2->SetBinError(i+1,0);
  }
  return;
}

void NTupleBinM(){
  for (int elec=0; elec<1; elec++){
    RooUnfoldResponse *m_BinM, *m_BinM_1, *m_BinM_2;
    TH1D *m_gen, *m_reco, *m_gen_1, *m_reco_1, *m_gen_2, *m_reco_2;
    GetBinM(m_BinM, m_gen, m_reco, m_BinM_1, m_gen_1, m_reco_1, m_BinM_2, m_gen_2, m_reco_2, 1, elec);

    RooUnfoldResponse *p_BinM, *p_BinM_1, *p_BinM_2;
    TH1D *p_gen, *p_reco, *p_gen_1, *p_reco_1, *p_gen_2, *p_reco_2;
    GetBinM(p_BinM, p_gen, p_reco, p_BinM_1, p_gen_1, p_reco_1, p_BinM_2, p_gen_2, p_reco_2, 0, elec);

    RooUnfold::ErrorTreatment witherror=RooUnfold::kCovToy;
    //RooUnfold::ErrorTreatment witherror=3;

    RooUnfoldBayes  unfold_MM(m_BinM,m_reco, 4);
    TH1D* h_BinM_MM  = (TH1D*) unfold_MM.Hreco(witherror);  
    TH1D* h_BinM_MM_t= (TH1D*) unfold_MM.Hreco(witherror);
    h_BinM_MM->Divide(h_BinM_MM_t,m_gen);
 
    RooUnfoldBayes  unfold_M1M2(m_BinM_1,m_reco_2, 4);
    TH1D* h_BinM_M1M2  = (TH1D*) unfold_M1M2.Hreco(witherror);  
    TH1D* h_BinM_M1M2_t= (TH1D*) unfold_M1M2.Hreco(witherror);
    h_BinM_M1M2->Divide(h_BinM_M1M2_t,m_gen_2);
 
    RooUnfoldBayes  unfold_M2M1(m_BinM_2,m_reco_1, 4);
    TH1D* h_BinM_M2M1  = (TH1D*) unfold_M2M1.Hreco(witherror);  
    TH1D* h_BinM_M2M1_t= (TH1D*) unfold_M2M1.Hreco(witherror);
    h_BinM_M2M1->Divide(h_BinM_M2M1_t,m_gen_1);
 
    RooUnfoldBayes  unfold_PP(p_BinM,p_reco, 4);
    TH1D* h_BinM_PP  = (TH1D*) unfold_PP.Hreco(witherror);  
    TH1D* h_BinM_PP_t= (TH1D*) unfold_PP.Hreco(witherror);
    h_BinM_PP->Divide(h_BinM_PP_t,p_gen);
 
    RooUnfoldBayes  unfold_P1P2(p_BinM_1,p_reco_2, 4);
    TH1D* h_BinM_P1P2  = (TH1D*) unfold_P1P2.Hreco(witherror);  
    TH1D* h_BinM_P1P2_t= (TH1D*) unfold_P1P2.Hreco(witherror);
    h_BinM_P1P2->Divide(h_BinM_P1P2_t,p_gen_2);
 
    RooUnfoldBayes  unfold_P2P1(p_BinM_2,p_reco_1, 4);
    TH1D* h_BinM_P2P1  = (TH1D*) unfold_P2P1.Hreco(witherror);  
    TH1D* h_BinM_P2P1_t= (TH1D*) unfold_P2P1.Hreco(witherror);
    h_BinM_P2P1->Divide(h_BinM_P2P1_t,p_gen_1);
 
    RooUnfoldBayes  unfold_MP(m_BinM,p_reco, 7);
    TH1D* h_BinM_MP  = (TH1D*) unfold_MP.Hreco(witherror);  
    TH1D* h_BinM_MP_t= (TH1D*) unfold_MP.Hreco(witherror);
    h_BinM_MP->Divide(h_BinM_MP_t,p_gen);
 
    RooUnfoldBayes  unfold_PM(p_BinM,m_reco, 7);
    TH1D* h_BinM_PM  = (TH1D*) unfold_PM.Hreco(witherror);  
    TH1D* h_BinM_PM_t= (TH1D*) unfold_PM.Hreco(witherror);
    h_BinM_PM->Divide(h_BinM_PM_t,m_gen);
 
    TLatex mark;
    mark.SetTextSize(0.04);
    mark.SetTextFont(42);
    mark.SetNDC(true);

    TCanvas* BinMigration_MM = new TCanvas("BinMigration_MM","BinMigration_MM",800,900);
    BinMigration_MM->cd();
    BinMigration_MM->SetLogx();
    h_BinM_MM->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_MM->GetXaxis()->SetTitle("#phi*");
    h_BinM_MM->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_MM->GetXaxis()->SetTitleSize(0.04);
    h_BinM_MM->GetXaxis()->SetLabelOffset(-0.01);
    h_BinM_MM->GetXaxis()->SetLabelSize(0.04);
    h_BinM_MM->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_MM->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_MM->GetYaxis()->SetTitleSize(0.04);
    h_BinM_MM->GetYaxis()->SetLabelSize(0.03);
    h_BinM_MM->GetYaxis()->SetRangeUser(0.9,1.1);
    h_BinM_MM->SetMarkerStyle(20);
    h_BinM_MM->SetStats(0);
    h_BinM_MM->SetBit( TH1::kNoTitle, true );
    h_BinM_MM->SetLineColor(1);
    h_BinM_MM->Draw();
    TLegend* leg_MM = new TLegend(0.45,0.80,0.85,0.86);
    leg_MM->SetFillStyle(0);
    leg_MM->SetBorderSize(0);
    leg_MM->SetLineWidth(1);
    leg_MM->SetNColumns(1);
    leg_MM->SetTextFont(42);
    leg_MM->SetTextSize(0.04);
    leg_MM->AddEntry(h_BinM_MM,"MadGraph","PL");
    leg_MM->Draw();
    mark.DrawLatex(0.19,0.17,"Unfolded using MadGraph");
    BinMigration_MM->SaveAs("BinM_MM.pdf");
    BinMigration_MM->SaveAs("BinM_MM.C");

    TCanvas* BinMigration_M1M2 = new TCanvas("BinMigration_M1M2","BinMigration_M1M2",800,900);
    BinMigration_M1M2->cd();
    BinMigration_M1M2->SetLogx();
    h_BinM_M1M2->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_M1M2->GetXaxis()->SetTitle("#phi*");
    h_BinM_M1M2->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_M1M2->GetXaxis()->SetTitleSize(0.04);
    h_BinM_M1M2->GetXaxis()->SetLabelOffset(-0.01);
    h_BinM_M1M2->GetXaxis()->SetLabelSize(0.04);
    h_BinM_M1M2->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_M1M2->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_M1M2->GetYaxis()->SetTitleSize(0.04);
    h_BinM_M1M2->GetYaxis()->SetLabelSize(0.03);
    h_BinM_M1M2->GetYaxis()->SetRangeUser(0.95,1.05);
    h_BinM_M1M2->SetMarkerStyle(21);
    h_BinM_M1M2->SetStats(0);
    h_BinM_M1M2->SetBit( TH1::kNoTitle, true );
    h_BinM_M1M2->SetLineColor(2);
    h_BinM_M1M2->SetMarkerColor(2);
    h_BinM_M1M2->Draw();
    h_BinM_M2M1->SetMarkerStyle(20);
    h_BinM_M2M1->SetLineColor(1);
    h_BinM_M2M1->Draw("same");
    TLegend* leg_M1M2 = new TLegend(0.45,0.77,0.85,0.91);
    leg_M1M2->SetFillStyle(0);
    leg_M1M2->SetBorderSize(0);
    leg_M1M2->SetLineWidth(1);
    leg_M1M2->SetNColumns(1);
    leg_M1M2->SetTextFont(42);
    leg_M1M2->SetTextSize(0.04);
    leg_M1M2->AddEntry(h_BinM_M2M1,"MadGraph part 1","PL");
    leg_M1M2->AddEntry(h_BinM_M1M2,"MadGraph part 2","PL");
    leg_M1M2->Draw();
    mark.DrawLatex(0.19,0.17,"Unfolded using MadGraph");
    BinMigration_M1M2->SaveAs("BinM_M1M2.pdf");
    BinMigration_M1M2->SaveAs("BinM_M1M2.C");

    TCanvas* BinMigration_PP = new TCanvas("BinMigration_PP","BinMigration_PP",800,900);
    BinMigration_PP->cd();
    BinMigration_PP->SetLogx();
    h_BinM_PP->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_PP->GetXaxis()->SetTitle("#phi*");
    h_BinM_PP->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_PP->GetXaxis()->SetTitleSize(0.04);
    h_BinM_PP->GetXaxis()->SetLabelOffset(-0.01);
    h_BinM_PP->GetXaxis()->SetLabelSize(0.04);
    h_BinM_PP->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_PP->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_PP->GetYaxis()->SetTitleSize(0.04);
    h_BinM_PP->GetYaxis()->SetLabelSize(0.03);
    h_BinM_PP->GetYaxis()->SetRangeUser(0.95,1.05);
    h_BinM_PP->SetMarkerStyle(20);
    h_BinM_PP->SetStats(0);
    h_BinM_PP->SetBit( TH1::kNoTitle, true );
    h_BinM_PP->SetLineColor(1);
    h_BinM_PP->Draw();
    TLegend* leg_PP = new TLegend(0.45,0.80,0.85,0.86);
    leg_PP->SetFillStyle(0);
    leg_PP->SetBorderSize(0);
    leg_PP->SetLineWidth(1);
    leg_PP->SetNColumns(1);
    leg_PP->SetTextFont(42);
    leg_PP->SetTextSize(0.04);
    leg_PP->AddEntry(h_BinM_PP,"Powheg","PL");
    leg_PP->Draw();
    mark.DrawLatex(0.19,0.17,"Unfolded using Powheg");
    BinMigration_PP->SaveAs("BinM_PP.pdf");
    BinMigration_PP->SaveAs("BinM_PP.C");

    TCanvas* BinMigration_P1P2 = new TCanvas("BinMigration_P1P2","BinMigration_P1P2",800,900);
    BinMigration_P1P2->cd();
    BinMigration_P1P2->SetLogx();
    h_BinM_P1P2->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_P1P2->GetXaxis()->SetTitle("#phi*");
    h_BinM_P1P2->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_P1P2->GetXaxis()->SetTitleSize(0.04);
    h_BinM_P1P2->GetXaxis()->SetLabelOffset(-0.01);
    h_BinM_P1P2->GetXaxis()->SetLabelSize(0.04);
    h_BinM_P1P2->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_P1P2->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_P1P2->GetYaxis()->SetTitleSize(0.04);
    h_BinM_P1P2->GetYaxis()->SetLabelSize(0.03);
    h_BinM_P1P2->GetYaxis()->SetRangeUser(0.95,1.05);
    h_BinM_P1P2->SetMarkerStyle(21);
    h_BinM_P1P2->SetStats(0);
    h_BinM_P1P2->SetBit( TH1::kNoTitle, true );
    h_BinM_P1P2->SetLineColor(2);
    h_BinM_P1P2->SetMarkerColor(2);
    h_BinM_P1P2->Draw();
    h_BinM_P2P1->SetMarkerStyle(20);
    h_BinM_P2P1->SetLineColor(1);
    h_BinM_P2P1->Draw("same");
    TLegend* leg_P1P2 = new TLegend(0.45,0.77,0.85,0.91);
    leg_P1P2->SetFillStyle(0);
    leg_P1P2->SetBorderSize(0);
    leg_P1P2->SetLineWidth(1);
    leg_P1P2->SetNColumns(1);
    leg_P1P2->SetTextSize(0.04);
    leg_P1P2->SetTextFont(42);
    leg_P1P2->AddEntry(h_BinM_P2P1,"Powheg part 1","PL");
    leg_P1P2->AddEntry(h_BinM_P1P2,"Powheg part 2","PL");
    leg_P1P2->Draw();
    mark.DrawLatex(0.19,0.17,"Unfolded using Powheg");
    BinMigration_P1P2->SaveAs("BinM_P1P2.pdf");
    BinMigration_P1P2->SaveAs("BinM_P1P2.C");

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
    h_BinM_MP->GetYaxis()->SetRangeUser(0.95,1.05);
    h_BinM_MP->SetStats(0);
    h_BinM_MP->SetBit( TH1::kNoTitle, true );
    h_BinM_MP->SetLineColor(1);
    h_BinM_MP->SetMarkerStyle(20);
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
    mark.DrawLatex(0.19,0.17,"Unfolded using MadGraph");
    BinMigration_MP->SaveAs("BinM_MP.pdf");
    BinMigration_MP->SaveAs("BinM_MP.C");

    TCanvas* BinMigration_PM = new TCanvas("BinMigration_PM","BinMigration_PM",800,900);
    BinMigration_PM->cd();
    BinMigration_PM->SetLogx();
    h_BinM_PM->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_PM->GetXaxis()->SetTitle("#phi*");
    h_BinM_PM->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_PM->GetXaxis()->SetLabelOffset(-0.01);
    h_BinM_PM->GetXaxis()->SetTitleSize(0.04);
    h_BinM_PM->GetXaxis()->SetLabelSize(0.04);
    h_BinM_PM->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_PM->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_PM->GetYaxis()->SetTitleSize(0.04);
    h_BinM_PM->GetYaxis()->SetLabelSize(0.03);
    h_BinM_PM->GetYaxis()->SetRangeUser(0.95,1.05);
    h_BinM_PM->SetStats(0);
    h_BinM_PM->SetBit( TH1::kNoTitle, true );
    h_BinM_PM->SetLineColor(1);
    h_BinM_PM->SetMarkerStyle(20);
    h_BinM_PM->Draw();
    TLegend* leg_PM = new TLegend(0.45,0.80,0.85,0.86);
    leg_PM->SetFillStyle(0);
    leg_PM->SetBorderSize(0);
    leg_PM->SetLineWidth(1);
    leg_PM->SetNColumns(1);
    leg_PM->SetTextFont(42);
    leg_PM->SetTextSize(0.04);
    leg_PM->AddEntry(h_BinM_PM,"MadGraph","PL");
    leg_PM->Draw();
    mark.DrawLatex(0.19,0.17,"Unfolded using Powheg");
    BinMigration_PM->SaveAs("BinM_PM.pdf");
    BinMigration_PM->SaveAs("BinM_PM.C");
  }
}
