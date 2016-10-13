#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TMatrixDSparse.h"
#include "TMatrix.h"
#include "TVector.h"
#include "TArray.h"
#include "TAxis.h"
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
using namespace std;

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,5.0,10.0};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
const double yBins[] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};
size_t ny=(sizeof(yBins)/sizeof(yBins[0]))-1;
size_t nbins=nphistar*ny;
std::string FileType = ".pdf";
std::string File_Signal_reco_b="/data/whybee0a/user/lesko_2/PhiStar/NicolesrootFiles/MG_born_reco.root";
std::string File_Powheg_reco_b="/data/whybee0a/user/lesko_2/PhiStar/NicolesrootFiles/PH_born_reco.root";
std::string File_Signal_reco_n="/data/whybee0a/user/lesko_2/PhiStar/NicolesrootFiles/MG_bare_reco.root";
std::string File_Powheg_reco_n="/data/whybee0a/user/lesko_2/PhiStar/NicolesrootFiles/PH_bare_reco.root";
std::string File_Signal_reco_d="/data/whybee0a/user/lesko_2/PhiStar/NicolesrootFiles/MG_dressed_reco.root";
std::string File_Powheg_reco_d="/data/whybee0a/user/lesko_2/PhiStar/NicolesrootFiles/PH_dressed_reco.root";

std::string reco_name="Combined Single Reco";

int N_MC=-1;

int GetBin(double phistar, double y){
  TAxis* A_phistar=new TAxis(nphistar,phistarBins);
  TAxis* A_y=new TAxis(ny,yBins);
  int bin_phistar=A_phistar->FindBin(phistar)-1;
  int bin_y=A_y->FindBin(fabs(y))-1;
  int bin=bin_y*nphistar+bin_phistar;
  if (bin_phistar<0 || bin_y<0) bin=-1;
  if (bin_phistar>=nphistar || bin_y>=ny) bin=nbins;
  // cout<<phistar<<" "<<y<<" "<<bin<<endl;
  return bin;
}

void GetBinM(RooUnfoldResponse* &BinM, TH1D* &h_gen, TH1D* &h_reco, RooUnfoldResponse* &BinM_1, TH1D* &h_gen_1, TH1D* &h_reco_1, RooUnfoldResponse* &BinM_2, TH1D* &h_gen_2, TH1D* &h_reco_2, bool madgraph=1, int elec=0){

  TH2D* BinMigration=new TH2D("BinMigration","BinMigration",nbins,0,nbins,nbins,0,nbins);
  BinMigration->Sumw2();
  TH1D* Dphistar=new TH1D("Dphistar","Dphistar",2000,-1,1);
  Dphistar->Sumw2();
  TH1D* Dy=new TH1D("Dy","Dy",2000,-1,1);
  Dy->Sumw2();
  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  if (madgraph){
    if (elec==0) nfiles=t->Add(File_Signal_reco_d.c_str());
    if (elec==1) nfiles=t->Add(File_Signal_reco_b.c_str());
    if (elec==2) nfiles=t->Add(File_Signal_reco_n.c_str());
  }
  else {
    if (elec==0) nfiles=t->Add(File_Powheg_reco_d.c_str());
    if (elec==1) nfiles=t->Add(File_Powheg_reco_b.c_str());
    if (elec==2) nfiles=t->Add(File_Powheg_reco_n.c_str());
  }

  TBranch *b_reco=t->GetBranch("reco");
  TBranch *b_truth=t->GetBranch("truth");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  if (elec==1) l_phistar_true=b_truth->GetLeaf("z_phistar_born");
  if (elec==2) l_phistar_true=b_truth->GetLeaf("z_phistar_naked");
  TLeaf *l_y=b_reco->GetLeaf("z_y");
  TLeaf *l_y_true=b_truth->GetLeaf("z_y");
  if (elec==1) l_y_true=b_reco->GetLeaf("z_yBorn");
  if (elec==2) l_y_true=b_reco->GetLeaf("z_yNaked");
  int nweights;
  t->SetBranchAddress("weight_size",&nweights);
  t->GetEntry(0);
  cout<<"The sample has nweights: "<<nweights<<endl;
  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);

  cout<<"Entries: "<<t->GetEntries()<<endl;
  h_gen   = new TH1D("phistar","phistar",nbins,0,nbins);
  h_gen->Sumw2();
  h_reco  = new TH1D("phistar","phistar",nbins,0,nbins);
  h_reco->Sumw2();
  h_gen_1 = new TH1D("phistar","phistar",nbins,0,nbins);
  h_gen_1->Sumw2();
  h_reco_1= new TH1D("phistar","phistar",nbins,0,nbins);
  h_reco_1->Sumw2();
  h_gen_2 = new TH1D("phistar","phistar",nbins,0,nbins);
  h_gen_2->Sumw2();
  h_reco_2= new TH1D("phistar","phistar",nbins,0,nbins);
  h_reco_2->Sumw2();
  BinM  =new RooUnfoldResponse (h_reco,h_gen);
  BinM_1=new RooUnfoldResponse (h_reco,h_gen);
  BinM_2=new RooUnfoldResponse (h_reco,h_gen);

  cout<<"reading signal reco "<<endl;
  for (int i=0; i<t->GetEntries();i++){
    // if (i>20000) continue;
    //for (int i=0; i<50000;i++){
    t->GetEntry(i);
    double weight =1;
    for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
    }
    double phistar=l_phistar->GetValue();
    double phistar_true=l_phistar_true->GetValue();
    Dphistar->Fill((phistar_true-phistar)/phistar_true,weight);
    double y=l_y->GetValue();
    double y_true=l_y_true->GetValue();
    Dy->Fill((y_true-y)/y_true,weight);

    int bin=GetBin(phistar,y);
    int bin_true=GetBin(phistar_true,y_true);

    h_gen ->Fill(bin_true,weight);
    h_reco->Fill(bin,weight);
    if (N_MC==-1 || i<N_MC) BinM  ->Fill(bin,bin_true,weight); 
    BinMigration->Fill(bin,bin_true,weight);
    if (i%2==0){
      h_gen_1 ->Fill(bin_true,weight);
      h_reco_1->Fill(bin,weight);
      if (N_MC==-1 || i<N_MC) BinM_1  ->Fill(bin,bin_true,weight); 
    }
    else{
      h_gen_2 ->Fill(bin_true,weight);
      h_reco_2->Fill(bin,weight);
      if (N_MC==-1 || i<N_MC) BinM_2  ->Fill(bin,bin_true,weight); 
    }
  }
////////////////////////////here is where I left.
  TLatex mark;
  mark.SetTextSize(0.04);
  mark.SetTextFont(42);
  mark.SetNDC(true);

//  BinMigration->GetXaxis()->SetRangeUser(0.001,3.2);
  BinMigration->GetXaxis()->SetTitle("(#phi*,y) bin (reconstructed)");
  BinMigration->GetXaxis()->SetTitleOffset(0.8);
  BinMigration->GetXaxis()->SetTitleSize(0.04);
  BinMigration->GetXaxis()->SetLabelOffset(0);
  BinMigration->GetXaxis()->SetLabelSize(0.03);
  BinMigration->GetYaxis()->SetTitleOffset(1.05);
  BinMigration->GetYaxis()->SetTitleSize(0.04);
  BinMigration->GetYaxis()->SetLabelSize(0.03);
//  BinMigration->GetYaxis()->SetRangeUser(0.001,3.2);
  BinMigration->GetYaxis()->SetTitle("(#phi*,y) bin (generated)");
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
  Dphistar->SetMarkerStyle(20);
  Dphistar->SetMarkerColor(1);
  Dphistar->SetLineColor(1);
  Dphistar->SetStats(0);
  Dphistar->SetBit( TH2::kNoTitle, true );

  Dy->GetXaxis()->SetTitle("(Z(y)(generated) - Z(y)(reconstructed)/Z(y) (generated)");
  Dy->GetXaxis()->SetTitleOffset(1.0);
  Dy->GetXaxis()->SetTitleSize(0.04);
  Dy->GetXaxis()->SetLabelOffset(0.0);
  Dy->GetXaxis()->SetLabelSize(0.04);
  Dy->GetYaxis()->SetTitleOffset(1.05);
  Dy->GetYaxis()->SetTitleSize(0.04);
  Dy->GetYaxis()->SetLabelSize(0.04);
  Dy->GetYaxis()->SetTitle("a.u.");
  Dy->SetMarkerStyle(20);
  Dy->SetMarkerColor(1);
  Dy->SetLineColor(1);
  Dy->SetStats(0);
  Dy->SetBit( TH2::kNoTitle, true );
   
  std::string plotname = "Plots/BinM_";
  if (elec==0) plotname+="Dressed_";
  if (elec==1) plotname+="Born_";
  if (elec==2) plotname+="Naked_";
  if (madgraph)plotname+="MG_";
  else         plotname+="PH_";

  if (madgraph){
    TCanvas* BinM_M = new TCanvas("BinM_M","BinM_M",800,900);
    BinM_M->cd();
    // BinM_M->SetLogx();
    // BinM_M->SetLogy();
    BinM_M->SetLogz();
    BinMigration->Draw("COLZ");
    mark.DrawLatex(0.19,0.80,"MadGraph");
    BinM_M->SaveAs((plotname+"res"+FileType).c_str());
    TCanvas* Dp_M = new TCanvas("Dp_M","Dp_M",800,900);
    Dp_M->cd();
    Dp_M->SetLogy();
    Dphistar->Draw();
    mark.DrawLatex(0.19,0.80,"MadGraph");
    Dp_M->SaveAs((plotname+"Dphistar"+FileType).c_str());
    TCanvas* Dy_M = new TCanvas("Dy_M","Dy_M",800,900);
    Dy_M->cd();
    Dy_M->SetLogy();
    Dy->Draw();
    mark.DrawLatex(0.19,0.80,"MadGraph");
    Dy_M->SaveAs((plotname+"Dy"+FileType).c_str());
  }
  else{
    TCanvas* BinM_P = new TCanvas("BinM_P","BinM_P",800,900);
    BinM_P->cd();
    // BinM_P->SetLogx();
    // BinM_P->SetLogy();
    BinM_P->SetLogz();
    BinMigration->Draw("COLZ");
    mark.DrawLatex(0.19,0.80,"Powheg");
    BinM_P->SaveAs((plotname+"res"+FileType).c_str());
    TCanvas* Dp_P = new TCanvas("Dp_P","Dp_P",800,900);
    Dp_P->cd();
    Dp_P->SetLogy();
    Dphistar->Draw();
    mark.DrawLatex(0.19,0.80,"Powheg");
    Dp_P->SaveAs((plotname+"Dphistar"+FileType).c_str());
    TCanvas* Dy_P = new TCanvas("Dy_P","Dy_P",800,900);
    Dy_P->cd();
    Dy_P->SetLogy();
    Dy->Draw();
    mark.DrawLatex(0.19,0.80,"Powheg");
    Dy_P->SaveAs((plotname+"Dy"+FileType).c_str());
  }

  double weighted_events= BinMigration->GetSumOfWeights();  
  double diagonal=0;
  for (uint i=0; i<nbins; i++){
    diagonal+=BinMigration->GetBinContent(i+1,i+1);
  }
  cout<<weighted_events<<" of which diagonal "<<diagonal<<" non diagonal faction:"<<(weighted_events-diagonal)/weighted_events<<endl;
  //Find amount of bin migration of Y
  double InTheY=0;//Gets values in the Same Y region 
  for (uint i=1; i<=nphistar; i++){
      for (uint j=1; j<=nphistar; j++){
          InTheY+=BinMigration->GetBinContent(i,j);
          InTheY+=BinMigration->GetBinContent(i+nphistar,j+nphistar);
          InTheY+=BinMigration->GetBinContent(i+nphistar*2,j+nphistar*2);
          InTheY+=BinMigration->GetBinContent(i+nphistar*3,j+nphistar*3);
          InTheY+=BinMigration->GetBinContent(i+nphistar*4,j+nphistar*4);
          InTheY+=BinMigration->GetBinContent(i+nphistar*5,j+nphistar*5);
      }
      
  }
  cout<<"Percent mig Y axis: "<<100*(weighted_events-InTheY)/weighted_events<<std::endl<<std::endl;
  
  cout<<"done reading data "<<endl;

  for (uint i=0; i<nbins;i++){
    h_gen->SetBinError(i+1,0);
    h_gen_1->SetBinError(i+1,0);
    h_gen_2->SetBinError(i+1,0);
  }
  return;
}

void NTupleBinM(){
  for (int elec=1; elec<2; elec++){
    RooUnfoldResponse *m_BinM, *m_BinM_1, *m_BinM_2;
    TH1D *m_gen, *m_reco, *m_gen_1, *m_reco_1, *m_gen_2, *m_reco_2;
    GetBinM(m_BinM, m_gen, m_reco, m_BinM_1, m_gen_1, m_reco_1, m_BinM_2, m_gen_2, m_reco_2, 1, elec);

    RooUnfoldResponse *p_BinM, *p_BinM_1, *p_BinM_2;
    TH1D *p_gen, *p_reco, *p_gen_1, *p_reco_1, *p_gen_2, *p_reco_2;
    GetBinM(p_BinM, p_gen, p_reco, p_BinM_1, p_gen_1, p_reco_1, p_BinM_2, p_gen_2, p_reco_2, 0, elec);

    // RooUnfold::ErrorTreatment witherror=RooUnfold::kCovToy;
    // //RooUnfold::ErrorTreatment witherror=3;

    RooUnfoldBayes  unfold_MM(m_BinM,m_reco, 4);
    TH1D* h_BinM_MM  = (TH1D*) unfold_MM.Hreco();  
    TH1D* h_BinM_MM_t= (TH1D*) unfold_MM.Hreco();
    h_BinM_MM->Divide(h_BinM_MM_t,m_gen);
    // double chi2_MM=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_MM+=pow((h_BinM_MM->GetBinContent(i+1)-1)/h_BinM_MM->GetBinError(i+1),2);
    // }
    // chi2_MM=chi2_MM/double(nbins);

    RooUnfoldBayes  unfold_M1M2(m_BinM_1,m_reco_2, 4);
    TH1D* h_BinM_M1M2  = (TH1D*) unfold_M1M2.Hreco();  
    TH1D* h_BinM_M1M2_t= (TH1D*) unfold_M1M2.Hreco();
    h_BinM_M1M2->Divide(h_BinM_M1M2_t,m_gen_2);
    // double chi2_M1M2=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_M1M2+=pow((h_BinM_M1M2->GetBinContent(i+1)-1)/h_BinM_M1M2->GetBinError(i+1),2);
    // }
    // chi2_M1M2=chi2_M1M2/double(nbins);
 
    RooUnfoldBayes  unfold_M2M1(m_BinM_2,m_reco_1, 4);
    TH1D* h_BinM_M2M1  = (TH1D*) unfold_M2M1.Hreco();  
    TH1D* h_BinM_M2M1_t= (TH1D*) unfold_M2M1.Hreco();
    h_BinM_M2M1->Divide(h_BinM_M2M1_t,m_gen_1);
    // double chi2_M2M1=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_M2M1+=pow((h_BinM_M2M1->GetBinContent(i+1)-1)/h_BinM_M2M1->GetBinError(i+1),2);
    // }
    // chi2_M2M1=chi2_M2M1/double(nbins);
 
    RooUnfoldBayes  unfold_PP(p_BinM,p_reco, 4);
    TH1D* h_BinM_PP  = (TH1D*) unfold_PP.Hreco();  
    TH1D* h_BinM_PP_t= (TH1D*) unfold_PP.Hreco();
    h_BinM_PP->Divide(h_BinM_PP_t,p_gen);
    // double chi2_PP=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_PP+=pow((h_BinM_PP->GetBinContent(i+1)-1)/h_BinM_PP->GetBinError(i+1),2);
    // }
    // chi2_PP=chi2_PP/double(nbins);
 
    RooUnfoldBayes  unfold_P1P2(p_BinM_1,p_reco_2, 4);
    TH1D* h_BinM_P1P2  = (TH1D*) unfold_P1P2.Hreco();  
    TH1D* h_BinM_P1P2_t= (TH1D*) unfold_P1P2.Hreco();
    h_BinM_P1P2->Divide(h_BinM_P1P2_t,p_gen_2);
    // double chi2_P1P2=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_P1P2+=pow((h_BinM_P1P2->GetBinContent(i+1)-1)/h_BinM_P1P2->GetBinError(i+1),2);
    // }
    // chi2_P1P2=chi2_P1P2/double(nbins);
 
    RooUnfoldBayes  unfold_P2P1(p_BinM_2,p_reco_1, 4);
    TH1D* h_BinM_P2P1  = (TH1D*) unfold_P2P1.Hreco();  
    TH1D* h_BinM_P2P1_t= (TH1D*) unfold_P2P1.Hreco();
    h_BinM_P2P1->Divide(h_BinM_P2P1_t,p_gen_1);
    // double chi2_P2P1=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_P2P1+=pow((h_BinM_P2P1->GetBinContent(i+1)-1)/h_BinM_P2P1->GetBinError(i+1),2);
    // }
    // chi2_P2P1=chi2_P2P1/double(nbins);
 
    RooUnfoldBayes  unfold_MP(m_BinM,p_reco, 7);
    TH1D* h_BinM_MP  = (TH1D*) unfold_MP.Hreco();  
    TH1D* h_BinM_MP_t= (TH1D*) unfold_MP.Hreco();
    h_BinM_MP->Divide(h_BinM_MP_t,p_gen);
    // double chi2_MP=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_MP+=pow((h_BinM_MP->GetBinContent(i+1)-1)/h_BinM_MP->GetBinError(i+1),2);
    // }
    // chi2_MP=chi2_MP/double(nbins);
 
    RooUnfoldBayes  unfold_PM(p_BinM,m_reco, 7);
    TH1D* h_BinM_PM  = (TH1D*) unfold_PM.Hreco();  
    TH1D* h_BinM_PM_t= (TH1D*) unfold_PM.Hreco();
    h_BinM_PM->Divide(h_BinM_PM_t,m_gen);
    // double chi2_PM=0;
    // for (int i=0; i<nbins;i++){
    //   chi2_PM+=pow((h_BinM_PM->GetBinContent(i+1)-1)/h_BinM_PM->GetBinError(i+1),2);
    //   std::cout<<i<<" "<<h_BinM_PM->GetBinContent(i+1)-1<<"  "<<h_BinM_PM->GetBinError(i+1)<<std::endl;
    // }
    // chi2_PM=chi2_PM/double(nbins);
 
    TLatex mark;
    mark.SetTextSize(0.04);
    mark.SetTextFont(42);
    mark.SetNDC(true);

    std::string plotname = "Plots/BinM_";
    if (elec==0) plotname+="Dressed_";
    if (elec==1) plotname+="Born_";
    if (elec==2) plotname+="Naked_";

    TCanvas* BinMigration_MM = new TCanvas("BinMigration_MM","BinMigration_MM",1600,900);
    BinMigration_MM->cd();
    // BinMigration_MM->SetLogx();
    // h_BinM_MM->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_MM->GetXaxis()->SetTitle("(#phi*,y) bin");
    h_BinM_MM->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_MM->GetXaxis()->SetTitleSize(0.04);
    h_BinM_MM->GetXaxis()->SetLabelOffset(0);
    h_BinM_MM->GetXaxis()->SetLabelSize(0.03);
    h_BinM_MM->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_MM->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_MM->GetYaxis()->SetTitleSize(0.04);
    h_BinM_MM->GetYaxis()->SetLabelSize(0.03);
    h_BinM_MM->GetYaxis()->SetRangeUser(0.85,1.15);
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
    // std::string chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_MM;
    // strs_MM << chi2_MM;
    // chi=chi+strs_MM.str();
    // mark.DrawLatex(0.19,0.80,chi.c_str());
    mark.DrawLatex(0.19,0.17,"Unfolded using MadGraph");
    BinMigration_MM->SaveAs((plotname+"MM"+FileType).c_str());
    // BinMigration_MM->SaveAs((plotname+"MM.C").c_str());

    TCanvas* BinMigration_M1M2 = new TCanvas("BinMigration_M1M2","BinMigration_M1M2",1600,900);
    BinMigration_M1M2->cd();
    // BinMigration_M1M2->SetLogx();
    // h_BinM_M1M2->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_M1M2->GetXaxis()->SetTitle("(#phi*,y) bin");
    h_BinM_M1M2->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_M1M2->GetXaxis()->SetTitleSize(0.04);
    h_BinM_M1M2->GetXaxis()->SetLabelOffset(0);
    h_BinM_M1M2->GetXaxis()->SetLabelSize(0.03);
    h_BinM_M1M2->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_M1M2->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_M1M2->GetYaxis()->SetTitleSize(0.04);
    h_BinM_M1M2->GetYaxis()->SetLabelSize(0.03);
    h_BinM_M1M2->GetYaxis()->SetRangeUser(0.85,1.15);
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
    // chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_M1M2;
    // strs_M1M2 << chi2_M1M2;
    // chi=chi+strs_M1M2.str();
    // // mark.DrawLatex(0.19,0.80,chi.c_str());
    // chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_M2M1;
    // strs_M2M1 << chi2_M2M1;
    // chi=chi+strs_M2M1.str();
    // // mark.SetTextColor(2);
    // // mark.DrawLatex(0.19,0.74,chi.c_str());
    // // mark.SetTextColor(1);
    mark.DrawLatex(0.19,0.17,"Unfolded using MadGraph");
    BinMigration_M1M2->SaveAs((plotname+"M1M2"+FileType).c_str());
    // BinMigration_M1M2->SaveAs((plotname+"M1M2.C").c_str());

    TCanvas* BinMigration_PP = new TCanvas("BinMigration_PP","BinMigration_PP",1600,900);
    BinMigration_PP->cd();
    // BinMigration_PP->SetLogx();
    // h_BinM_PP->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_PP->GetXaxis()->SetTitle("(#phi*,y) bin");
    h_BinM_PP->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_PP->GetXaxis()->SetTitleSize(0.04);
    h_BinM_PP->GetXaxis()->SetLabelOffset(0);
    h_BinM_PP->GetXaxis()->SetLabelSize(0.03);
    h_BinM_PP->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_PP->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_PP->GetYaxis()->SetTitleSize(0.04);
    h_BinM_PP->GetYaxis()->SetLabelSize(0.03);
    h_BinM_PP->GetYaxis()->SetRangeUser(0.85,1.15);
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
    // chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_PP;
    // strs_PP << chi2_PP;
    // chi=chi+strs_PP.str();
    // // mark.DrawLatex(0.19,0.80,chi.c_str());
    mark.DrawLatex(0.19,0.17,"Unfolded using Powheg");
    BinMigration_PP->SaveAs((plotname+"PP"+FileType).c_str());
    // BinMigration_PP->SaveAs((plotname+"PP.C").c_str());

    TCanvas* BinMigration_P1P2 = new TCanvas("BinMigration_P1P2","BinMigration_P1P2",1600,900);
    BinMigration_P1P2->cd();
    // BinMigration_P1P2->SetLogx();
    // h_BinM_P1P2->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_P1P2->GetXaxis()->SetTitle("(#phi*,y) bin");
    h_BinM_P1P2->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_P1P2->GetXaxis()->SetTitleSize(0.04);
    h_BinM_P1P2->GetXaxis()->SetLabelOffset(0);
    h_BinM_P1P2->GetXaxis()->SetLabelSize(0.03);
    h_BinM_P1P2->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_P1P2->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_P1P2->GetYaxis()->SetTitleSize(0.04);
    h_BinM_P1P2->GetYaxis()->SetLabelSize(0.03);
    h_BinM_P1P2->GetYaxis()->SetRangeUser(0.85,1.15);
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
    // chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_P1P2;
    // strs_P1P2 << chi2_P1P2;
    // chi=chi+strs_P1P2.str();
    // // mark.DrawLatex(0.19,0.80,chi.c_str());
    // chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_P2P1;
    // strs_P2P1 << chi2_P2P1;
    // chi=chi+strs_P2P1.str();
    // // mark.SetTextColor(2);
    // // mark.DrawLatex(0.19,0.74,chi.c_str());
    // // mark.SetTextColor(1);
    mark.DrawLatex(0.19,0.17,"Unfolded using Powheg");
    BinMigration_P1P2->SaveAs((plotname+"P1P2"+FileType).c_str());
    // BinMigration_P1P2->SaveAs((plotname+"P1P2.C").c_str());

    TCanvas* BinMigration_MP = new TCanvas("BinMigration_MP","BinMigration_MP",1600,900);
    BinMigration_MP->cd();
    // BinMigration_MP->SetLogx();
    // h_BinM_MP->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_MP->GetXaxis()->SetTitle("(#phi*,y) bin");
    h_BinM_MP->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_MP->GetXaxis()->SetTitleSize(0.04);
    h_BinM_MP->GetXaxis()->SetLabelOffset(0);
    h_BinM_MP->GetXaxis()->SetLabelSize(0.03);
    h_BinM_MP->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_MP->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_MP->GetYaxis()->SetTitleSize(0.04);
    h_BinM_MP->GetYaxis()->SetLabelSize(0.03);
    h_BinM_MP->GetYaxis()->SetRangeUser(0.85,1.15);
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
    // chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_MP;
    // strs_MP << chi2_MP;
    // chi=chi+strs_MP.str();
    // mark.DrawLatex(0.19,0.80,chi.c_str());
    mark.DrawLatex(0.19,0.17,"Unfolded using MadGraph");
    BinMigration_MP->SaveAs((plotname+"MP"+FileType).c_str());
    // BinMigration_MP->SaveAs((plotname+"MP.C").c_str());

    TCanvas* BinMigration_PM = new TCanvas("BinMigration_PM","BinMigration_PM",1600,900);
    BinMigration_PM->cd();
    // BinMigration_PM->SetLogx();
    // h_BinM_PM->GetXaxis()->SetRangeUser(0.001,3.2);
    h_BinM_PM->GetXaxis()->SetTitle("(#phi*,y) bin");
    h_BinM_PM->GetXaxis()->SetTitleOffset(0.8);
    h_BinM_PM->GetXaxis()->SetLabelOffset(0);
    h_BinM_PM->GetXaxis()->SetTitleSize(0.04);
    h_BinM_PM->GetXaxis()->SetLabelSize(0.03);
    h_BinM_PM->GetYaxis()->SetTitle("Unfolded/Generated");
    h_BinM_PM->GetYaxis()->SetTitleOffset(1.2);
    h_BinM_PM->GetYaxis()->SetTitleSize(0.04);
    h_BinM_PM->GetYaxis()->SetLabelSize(0.03);
    h_BinM_PM->GetYaxis()->SetRangeUser(0.85,1.15);
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
    // chi="#Chi^{2}/D.o.F=";
    // std::ostringstream strs_PM;
    // strs_PM << chi2_PM;
    // chi=chi+strs_PM.str();
    // mark.DrawLatex(0.19,0.80,chi.c_str());
    mark.DrawLatex(0.19,0.17,"Unfolded using Powheg");
    BinMigration_PM->SaveAs((plotname+"PM"+FileType).c_str());
    // BinMigration_PM->SaveAs((plotname+"PM.C").c_str());
  }
}
