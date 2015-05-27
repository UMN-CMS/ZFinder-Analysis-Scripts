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

//bool doGen=1;
bool doMass=1;

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_reco.root";
std::string File_Signal_gen="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_gen.root";

std::string reco_name="Combined Single Reco";
std::string gen_name ="Combined Gen Cuts Reco";

void readmcsample(TH1D &nom, TH1D &rw, TH1D &nomratio, bool doGen){
 cout<<"reading data for "<<File_Signal_reco<<endl;
  TChain* t;
  if (!doGen) t= new TChain(reco_name.c_str(),reco_name.c_str());
  else t= new TChain(gen_name.c_str(),gen_name.c_str());
  int nfiles;
  if (!doGen) nfiles=t->Add(File_Signal_reco.c_str());
  else nfiles=t->Add(File_Signal_gen.c_str());

  TBranch *b_reco=t->GetBranch("reco");
  TBranch *b_truth=t->GetBranch("truth");
  TLeaf *l_mass=b_reco->GetLeaf("z_y");
  if (doMass) l_mass=b_reco->GetLeaf("z_m");
  if (doGen) {
    l_mass=b_truth->GetLeaf("z_y");
    if (doMass) l_mass=b_truth->GetLeaf("z_m");
  }
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  int nweights;
  t->SetBranchAddress("weight_size",&nweights);
  t->GetEntry(0);
  cout<<"The sample has nweights: "<<nweights<<endl;
  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);
  
  // TFile f_zy("ratio_y.root");
  // ratio = (TH1F*)f_zy.Get("ratio");
  // f_zy.Close();
  TFile f_zm("ratio_zmass.root");
  TH1F *ratio = (TH1F*)f_zm.Get("ratio_histo");
  cout<<"Entries: "<<t->GetEntries()<<endl;
  for (int i=0; i<t->GetEntries();i++){
    t->GetEntry(i);
    double phistar=l_phistar->GetValue();
    double mass=l_mass->GetValue();
    double phistar_true=l_phistar_true->GetValue();
    double weight =1;
    for (int w=0; w<nweights;w++){
      if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
    }
    int idx=ratio->GetXaxis()->FindBin(mass);
    if (ratio->GetBinContent(idx)!=0){
     if (!doGen){
	nom.Fill(phistar,weight);
	nomratio.Fill(phistar,weight);
	weight=weight*ratio->GetBinContent(idx);
	rw.Fill(phistar,weight);
     }
      else{
	nom.Fill(phistar_true,weight);
	nomratio.Fill(phistar_true,weight);
	weight=weight*ratio->GetBinContent(idx);
	rw.Fill(phistar_true,weight);
      }
    }
  }
}


void NTupleZMass(){
  TH1D *phistar_reco= new TH1D("phistar","phistar",nphistar,phistarBins);
  TH1D *phistar_reco_ratio= new TH1D("phistar","phistar",nphistar,phistarBins);
  TH1D *phistar_reco_rw= new TH1D("phistar","phistar",nphistar,phistarBins);
  TH1D *phistar_gen= new TH1D("phistar","phistar",nphistar,phistarBins);
  TH1D *phistar_gen_ratio= new TH1D("phistar","phistar",nphistar,phistarBins);
  TH1D *phistar_gen_rw= new TH1D("phistar","phistar",nphistar,phistarBins);
  phistar_reco->Sumw2();
  phistar_reco_ratio->Sumw2();
  phistar_reco_rw->Sumw2();
  phistar_gen->Sumw2();
  phistar_gen_ratio->Sumw2();
  phistar_gen_rw->Sumw2();

  gErrorIgnoreLevel = kError;
  readmcsample(*phistar_reco, *phistar_reco_rw, *phistar_reco_ratio, 0);
  readmcsample(*phistar_gen, *phistar_gen_rw, *phistar_gen_ratio, 1); 
  double norm1=phistar_reco_rw->Integral();
  double norm2=phistar_reco->Integral();
  double norm3=phistar_gen_rw->Integral();
  double norm4=phistar_gen->Integral();
  phistar_reco_rw->Scale(1./norm1);
  phistar_reco->Scale(1./norm2);
  phistar_gen_rw->Scale(1./norm3);
  phistar_gen->Scale(1./norm4);
  phistar_reco_ratio->Divide(phistar_reco,phistar_reco_rw,1.,1.,"B");
  phistar_gen_ratio->Divide(phistar_gen,phistar_gen_rw,1.,1.,"B");

  TLatex mark;
  mark.SetTextSize(0.04);
  mark.SetTextFont(42);
  mark.SetNDC(true);
  TCanvas* zmassrw = new TCanvas("zmassrw","zmassrw",800,900);
  zmassrw->cd();
  zmassrw->SetLogx();
  phistar_reco_ratio->GetXaxis()->SetRangeUser(0.001,3.2);
  phistar_reco_ratio->GetXaxis()->SetTitle("#phi*");
  phistar_reco_ratio->GetXaxis()->SetTitleOffset(0.8);
  phistar_reco_ratio->GetXaxis()->SetTitleSize(0.04);
  phistar_reco_ratio->GetXaxis()->SetLabelOffset(-0.01);
  phistar_reco_ratio->GetXaxis()->SetLabelSize(0.04);
  phistar_reco_ratio->GetYaxis()->SetTitle("Nominal/Reweighed");
  phistar_reco_ratio->GetYaxis()->SetTitleOffset(1.2);
  phistar_reco_ratio->GetYaxis()->SetTitleSize(0.04);
  phistar_reco_ratio->GetYaxis()->SetLabelSize(0.03);
  // phistar_ratio->GetYaxis()->SetRangeUser(0.9,1.1);
  phistar_reco_ratio->SetMarkerStyle(20);
  phistar_reco_ratio->SetStats(0);
  phistar_reco_ratio->SetBit( TH1::kNoTitle, true );
  phistar_reco_ratio->Draw();
  phistar_gen_ratio->SetLineColor(2);
  phistar_gen_ratio->SetMarkerColor(2);
  phistar_gen_ratio->SetMarkerStyle(21);
  phistar_gen_ratio->Draw("same");
  TLegend* leg_mass = new TLegend(0.45,0.77,0.85,0.91);
  leg_mass->SetFillStyle(0);
  leg_mass->SetBorderSize(0);
  leg_mass->SetLineWidth(1);
  leg_mass->SetNColumns(1);
  leg_mass->SetTextFont(42);
  leg_mass->SetTextSize(0.04);
  leg_mass->AddEntry(phistar_reco_ratio,"#phi*_{reconstructed}","PL");
  leg_mass->AddEntry(phistar_gen_ratio,"#phi*_{generated}","PL");
  leg_mass->Draw();
  mark.DrawLatex(0.19,0.17,"MadGraph");
 }
