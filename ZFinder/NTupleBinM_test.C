#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGenericPdf.h"
#include "RooDataSet.h"
#include "RooAbsData.h"
#include "RooArgSet.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooKeysPdf.h"
#include "RooHistPdf.h"
#include "RooCategory.h"
#include "RooBinning.h"
#include "TMatrixDSparse.h"
#include "TRandom.h"
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
using namespace RooFit;
//using namespace RooStats;

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
const double phistarBinsPlot[] = {0.0012,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
RooBinning phistarbin(nphistar,phistarBins);

std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/Signal_reco.root";
std::string File_Powheg_reco="/afs/cern.ch/work/r/ruckstuh/Powheg_cteq_reco.root";

std::string reco_name="Combined Single Reco";

int N_MC=5000;

int N_TOY=500;

TGraphAsymmErrors* ConvertToTGraph(TH1D* h){
  TGraphAsymmErrors* g=new TGraphAsymmErrors(nphistar);
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    g->SetPoint(iphistar, (phistarBins[iphistar]+phistarBins[iphistar+1])/2., h->GetBinContent(iphistar+1));
    g->SetPointError(iphistar, 0, 0, h->GetBinError(iphistar+1), h->GetBinError(iphistar+1));
  }
  return g;
}

void GetBinM(vector<RooUnfoldResponse*> &BinM, vector<TH1D*> &h_gen, vector<TH1D*> &h_reco, bool madgraph=1, int elec=0){

  TH2D* BinMigration=new TH2D("BinMigration","BinMigration",nphistar,phistarBins,nphistar,phistarBins);
  BinMigration->Sumw2();
  TH1D* Gen=new TH1D("Gen","gen",nphistar,phistarBins);
  Gen->Sumw2();
  TH1D* Reco=new TH1D("Reco","Reco",nphistar,phistarBins);
  Reco->Sumw2();

  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  if (madgraph)  nfiles=t->Add(File_Signal_reco.c_str());
  else nfiles=t->Add(File_Powheg_reco.c_str());

  TBranch *b_reco=t->GetBranch("reco");
  TBranch *b_truth=t->GetBranch("truth");
  TBranch *b_event=t->GetBranch("event_info");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  TLeaf *l_en=b_event->GetLeaf("event_number");
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
  cout<<"reading signal reco "<<endl;
  for (int i=0; i<t->GetEntries();i++){
    // if (!madgraph && (i>N_MC+50000 || i<50000) &&N_MC!=-1) continue;
    if (!madgraph && i>N_MC &&N_MC!=-1) continue;
    //for (int i=0; i<50000;i++){
    t->GetEntry(i);
    double phistar=l_phistar->GetValue();
    double phistar_true=l_phistar_true->GetValue();
    double weight =1;
    double eventn=l_en->GetValue();
    // else en=eventn;
    for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
    }
    Gen ->Fill(phistar_true,weight);
    Reco->Fill(phistar,weight);
    BinMigration->Fill(phistar,phistar_true,weight); 
  }
  if (!madgraph){
    TCanvas* BM = new TCanvas("BM","BM",800,900);
    BinMigration->GetXaxis()->SetRangeUser(0.001,3.2);
    BinMigration->GetXaxis()->SetTitle("#phi^{*}(reconstructed)");
    BinMigration->GetXaxis()->SetTitleOffset(0.8);
    BinMigration->GetXaxis()->SetTitleSize(0.04);
    BinMigration->GetXaxis()->SetLabelOffset(-0.01);
    BinMigration->GetXaxis()->SetLabelSize(0.04);
    BinMigration->GetYaxis()->SetTitleOffset(1.05);
    BinMigration->GetYaxis()->SetTitleSize(0.04);
    BinMigration->GetYaxis()->SetLabelSize(0.04);
    BinMigration->GetYaxis()->SetRangeUser(0.001,3.2);
    BinMigration->GetYaxis()->SetTitle("#phi^{*}(generated)");
    BinMigration->GetZaxis()->SetTitleOffset(-0.004);
    BinMigration->SetStats(0);
    BinMigration->SetBit( TH2::kNoTitle, true );
    BM->SetLogx();
    BM->SetLogy();
    BM->SetLogz();
    BinMigration->Draw("COLZ");
    TLatex mark;
    mark.SetTextSize(0.04);
    mark.SetTextFont(42);
    mark.SetNDC(true);
    mark.DrawLatex(0.19,0.80,"Powheg");
  }
  RooUnfoldResponse* BinM1  =new RooUnfoldResponse (Reco,Gen,BinMigration);
  BinM.push_back(BinM1);
  h_gen.push_back(Gen);
  h_reco.push_back(Reco);
  for (int i=0; i<N_TOY; i++){
    TH1D* Recotemp=new TH1D("Reco","Reco",nphistar,phistarBins);
    Recotemp->Sumw2();
    TH2D* BinMigrationtemp=new TH2D("BinMigration","BinMigration",nphistar,phistarBins,nphistar,phistarBins);
    BinMigrationtemp->Sumw2();
    for (int j=0; j<nphistar; j++){
      double x=gRandom->Gaus(Reco->GetBinContent(j+1),Reco->GetBinError(j+1));
      Recotemp->SetBinContent(j+1,x);
      Recotemp->SetBinError(j+1,Reco->GetBinError(j+1));
      for (int k=0; k<nphistar; k++){
	double mean=BinMigration->GetBinContent(j+1,k+1);
	if (mean==0) continue;
	double error=BinMigration->GetBinError(j+1,k+1);
	if (mean/error<5){
	  x=gRandom->Poisson(mean);
	}
	else{
	  x=gRandom->Gaus(mean,error);
	}
	BinMigrationtemp->SetBinContent(j+1,k+1,x);
	BinMigrationtemp->SetBinError(j+1,k+1,error);
      }
    }
    TH1D* recotemp=BinMigrationtemp->ProjectionX();
    TH1D* gentemp=BinMigrationtemp->ProjectionY();
    RooUnfoldResponse* BinM1temp  =new RooUnfoldResponse (recotemp,gentemp,BinMigrationtemp);
    BinM.push_back(BinM1temp);
    h_reco.push_back(Recotemp);
  }

  cout<<"done reading data for "<<File_Signal_reco<<"  "<<reco_name<<endl;
  return;
}

void NTupleBinM_test(){
  for (int elec=0; elec<1; elec++){
    vector<RooUnfoldResponse *> m_BinM;
    vector<TH1D *> m_gen, m_reco;
    GetBinM(m_BinM, m_gen, m_reco, 1, elec);

    vector<RooUnfoldResponse *>p_BinM;
    vector<TH1D *> p_gen, p_reco;
    GetBinM(p_BinM, p_gen, p_reco, 0, elec);
 
    RooUnfoldBayes  unfold_PM(p_BinM[0],m_reco[0], 4);
    TH1D* h_BinM_PM  = (TH1D*) unfold_PM.Hreco();  
    TH1D* h_BinM_PM_t= (TH1D*) unfold_PM.Hreco();
    h_BinM_PM->Divide(h_BinM_PM_t,m_gen[0]);

    RooUnfoldBayes  unfold_MP(m_BinM[0],p_reco[0], 4);
    TH1D* h_BinM_MP  = (TH1D*) unfold_MP.Hreco();  
    TH1D* h_BinM_MP_t= (TH1D*) unfold_MP.Hreco();
    h_BinM_MP->Divide(h_BinM_MP_t,p_gen[0]);

    TGraphAsymmErrors* g_BinM_PM=ConvertToTGraph(h_BinM_PM); 
    TGraphAsymmErrors* g_BinM_MP=ConvertToTGraph(h_BinM_MP); 

    vector<TH1D *> h_Data, h_MC;
    for (int i=0; i<N_TOY+1; i++){
      cout<<i<<endl;
      RooUnfoldBayes  unfold_MP_Data(m_BinM[0],p_reco[i], 4);
      TH1D* h_BinM_MP_Data  = (TH1D*) unfold_MP_Data.Hreco();  
      TH1D* h_BinM_MP_Data_t= (TH1D*) unfold_MP_Data.Hreco();
      h_BinM_MP_Data->Divide(h_BinM_MP_Data_t,p_gen[0]);

      RooUnfoldBayes  unfold_PM_MC(p_BinM[i],m_reco[0], 4);
      TH1D* h_BinM_PM_MC  = (TH1D*) unfold_PM_MC.Hreco();  
      TH1D* h_BinM_PM_MC_t= (TH1D*) unfold_PM_MC.Hreco();
      h_BinM_PM_MC->Divide(h_BinM_PM_MC_t,m_gen[0]);

      h_Data.push_back(h_BinM_MP_Data);
      h_MC.push_back(h_BinM_PM_MC);

      cout<<i<<"done"<<endl;
    }

    TGraphAsymmErrors* MC_68=new TGraphAsymmErrors(nphistar);
    TGraphAsymmErrors* Data_68=new TGraphAsymmErrors(nphistar);

    TH2D* Toy_MC=new TH2D("Toy_MC","Toy_MC",nphistar,phistarBinsPlot,60,0.7,1.3);

    for (int j=0; j<nphistar; j++){
      vector<double> recoun;
      vector<double> binmun;
      for (int i=0; i<N_TOY+1; i++){
	// cout<<h_Data[i]->GetBinContent(j)<<endl;
	recoun.push_back(h_Data[i]->GetBinContent(j+1));
	binmun.push_back(h_MC[i]->GetBinContent(j+1));
	Toy_MC->Fill((phistarBins[j]+phistarBins[j+1])/2.,h_MC[i]->GetBinContent(j+1));
      }
      std::sort(binmun.begin(),binmun.end());
      std::sort(recoun.begin(),recoun.end());
      int idx68=0.1587*(N_TOY+1)-0.5;
      // int idx96=0.02275*(N_TOY+1)-0.5;
      double b_med=binmun[N_TOY/2];
      double b_min68=binmun[idx68];
      double b_max68=binmun[(N_TOY+1)-idx68];
      MC_68->SetPoint(j,(phistarBins[j]+phistarBins[j+1])/2.,b_med);
      MC_68->SetPointError(j,0,0,b_med-b_min68,b_max68-b_med);
      double r_med=recoun[N_TOY/2];
      double r_min68=recoun[idx68];
      double r_max68=recoun[(N_TOY+1)-idx68];
      Data_68->SetPoint(j,(phistarBins[j]+phistarBins[j+1])/2.,r_med);
      Data_68->SetPointError(j,0,0,r_med-r_min68,r_max68-r_med);
    }
  
    TLatex mark;
    mark.SetTextSize(0.04);
    mark.SetTextFont(42);
    mark.SetNDC(true);

    TCanvas* BinMigration_PM = new TCanvas("BinMigration_PM","BinMigration_PM",800,900);
    BinMigration_PM->cd();
    BinMigration_PM->SetLogx();
    g_BinM_PM->GetXaxis()->SetRangeUser(0.001,3.2);
    g_BinM_PM->GetXaxis()->SetTitle("#phi^{*}");
    g_BinM_PM->GetXaxis()->SetTitleOffset(0.8);
    g_BinM_PM->GetXaxis()->SetLabelOffset(-0.01);
    g_BinM_PM->GetXaxis()->SetTitleSize(0.04);
    g_BinM_PM->GetXaxis()->SetLabelSize(0.04);
    g_BinM_PM->GetYaxis()->SetTitle("Unfolded/Generated");
    g_BinM_PM->SetTitle("");
    g_BinM_PM->GetYaxis()->SetTitleOffset(1.2);
    g_BinM_PM->GetYaxis()->SetTitleSize(0.04);
    g_BinM_PM->GetYaxis()->SetLabelSize(0.03);
    g_BinM_PM->GetYaxis()->SetRangeUser(0.9,1.1);
    //    g_BinM_PM->SetStats(0);
    g_BinM_PM->SetBit( TH1::kNoTitle, true );
    g_BinM_PM->SetLineColor(1);
    g_BinM_PM->SetMarkerColor(1);
    g_BinM_PM->SetMarkerStyle(20);
    g_BinM_PM->Draw("PA");
    MC_68->SetLineColor(2);
    MC_68->SetMarkerColor(2);
    MC_68->SetMarkerStyle(21);
    MC_68->Draw("Psame");
    if (N_MC!=-1) g_BinM_PM->Draw("Psame");
    TLegend* leg_PM = new TLegend(0.45,0.77,0.85,0.91);
    leg_PM->SetFillStyle(0);
    leg_PM->SetBorderSize(0);
    leg_PM->SetLineWidth(1);
    leg_PM->SetNColumns(1);
    leg_PM->SetTextFont(42);
    leg_PM->SetTextSize(0.04);
    leg_PM->AddEntry(g_BinM_PM,"RooUnfold","PL");
    leg_PM->AddEntry(MC_68,"Toy MC","PL");
    leg_PM->Draw();
    mark.DrawLatex(0.19,0.21,"MadGraph");
    if (N_MC==-1) mark.DrawLatex(0.19,0.17,"unfolded using full Powheg sample");
    if (N_MC==5000) mark.DrawLatex(0.19,0.17,"unfolded using 5000 Powheg events");
    if (N_MC==50000) mark.DrawLatex(0.19,0.17,"unfolded using 50000 Powheg events");

    TCanvas* BinMigration_MP = new TCanvas("BinMigration_MP","BinMigration_MP",800,900);
    BinMigration_MP->cd();
    BinMigration_MP->SetLogx();
    g_BinM_MP->GetXaxis()->SetRangeUser(0.001,3.2);
    g_BinM_MP->GetXaxis()->SetTitle("#phi^{*}");
    g_BinM_MP->GetXaxis()->SetTitleOffset(0.8);
    g_BinM_MP->GetXaxis()->SetLabelOffset(-0.01);
    g_BinM_MP->GetXaxis()->SetTitleSize(0.04);
    g_BinM_MP->GetXaxis()->SetLabelSize(0.04);
    g_BinM_MP->GetYaxis()->SetTitle("Unfolded/Generated");
    g_BinM_MP->GetYaxis()->SetTitleOffset(1.2);
    g_BinM_MP->GetYaxis()->SetTitleSize(0.04);
    g_BinM_MP->GetYaxis()->SetLabelSize(0.03);
    g_BinM_MP->GetYaxis()->SetRangeUser(0.9,1.1);
    g_BinM_MP->SetTitle("");
    //    g_BinM_MP->SetStats(0);
    g_BinM_MP->SetBit( TH1::kNoTitle, true );
    g_BinM_MP->SetLineColor(1);
    g_BinM_MP->SetMarkerColor(1);
    g_BinM_MP->SetMarkerStyle(20);
    g_BinM_MP->Draw("PA");
    Data_68->SetLineColor(2);
    Data_68->SetMarkerColor(2);
    Data_68->SetMarkerStyle(21);
    Data_68->Draw("Psame");
    TLegend* leg_MP = new TLegend(0.45,0.77,0.85,0.91);
    leg_MP->SetFillStyle(0);
    leg_MP->SetBorderSize(0);
    leg_MP->SetLineWidth(1);
    leg_MP->SetNColumns(1);
    leg_MP->SetTextFont(42);
    leg_MP->SetTextSize(0.04);
    leg_MP->AddEntry(g_BinM_MP,"RooUnfold","PL");
    leg_MP->AddEntry(Data_68,"Toy MC","PL");
    leg_MP->Draw();
    if (N_MC==-1) mark.DrawLatex(0.19,0.21,"Full Powheg sample");
    if (N_MC==5000) mark.DrawLatex(0.19,0.21,"5000 Powheg events");
    if (N_MC==50000) mark.DrawLatex(0.19,0.21,"50000 Powheg events");
    mark.DrawLatex(0.19,0.17,"unfolded using MadGraph");

    TCanvas* ToyMC_PM = new TCanvas("ToyMC_PM","ToyMC_PM",800,900);
    ToyMC_PM->cd();
    ToyMC_PM->SetLogx();
    ToyMC_PM->SetLogz();
    Toy_MC->GetXaxis()->SetRangeUser(0.001,3.2);
    Toy_MC->GetXaxis()->SetTitle("#phi^{*}");
    Toy_MC->GetXaxis()->SetTitleOffset(0.8);
    Toy_MC->GetXaxis()->SetLabelOffset(-0.01);
    Toy_MC->GetXaxis()->SetTitleSize(0.04);
    Toy_MC->GetXaxis()->SetLabelSize(0.04);
    Toy_MC->GetYaxis()->SetTitle("Unfolded/Generated");
    Toy_MC->GetYaxis()->SetTitleOffset(1.2);
    Toy_MC->GetYaxis()->SetTitleSize(0.04);
    Toy_MC->GetYaxis()->SetLabelSize(0.03);
    //Toy_MC->GetYaxis()->SetRangeUser(0.8,1.2);
    Toy_MC->SetStats(0);
    Toy_MC->SetBit( TH1::kNoTitle, true );
    Toy_MC->Draw("COLZ");
    mark.DrawLatex(0.19,0.21,"MadGraph");
    if (N_MC==-1) mark.DrawLatex(0.19,0.17,"unfolded using full Powheg sample");
    if (N_MC==5000) mark.DrawLatex(0.19,0.17,"unfolded using 5000 Powheg events");
    if (N_MC==50000) mark.DrawLatex(0.19,0.17,"unfolded using 50000 Powheg events");
  }
}
