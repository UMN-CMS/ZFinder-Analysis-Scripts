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

using namespace std;

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,10};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

std::string File_MadGraphNoIso_reco="/data/whybee0a/user/lesko_2/fermi/NoIsoMC/MadReco.root";
std::string File_PowhegNoIso_reco="/data/whybee0a/user/lesko_2/fermi/NoIsoMC/PowReco.root";
std::string File_MadGraphNoIso_gen="/data/whybee0a/user/lesko_2/fermi/NoIsoMC/MadGen.root";
std::string File_PowhegNoIso_gen="/data/whybee0a/user/lesko_2/fermi/NoIsoMC/PowGen.root";



//data/whybee0a/user/lesko_2/fermi/NoIso/ZfinderPOWHEGoIso2015-12-21/NoIsoPow.root
//data/whybee0a/user/lesko_2/fermi/NoIso/ZfinderMadgraphNoIso2015-12-21/NoIsoMad.root
        
std::string File_MadGraph_reco="/data/whybee0a/user/lesko_2/fermi/WithIsoForEventEff/MadReco.root";
std::string File_Powheg_reco="/data/whybee0a/user/lesko_2/fermi/WithIsoForEventEff/PowReco.root";
std::string File_MadGraph_gen="/data/whybee0a/user/lesko_2/fermi/WithIsoForEventEff/MadGen.root";
std::string File_Powheg_gen="/data/whybee0a/user/lesko_2/fermi/WithIsoForEventEff/PowGen.root";




std::string reco_name="Combined Single Reco";
std::string gen_name ="Combined Gen Cuts Reco";

TH1D* GetPhistar(bool madgraph=false, bool gen=false, int elec=0, bool NoIso=false){
  std::string name=reco_name;
  if (gen) name=gen_name;
  TChain* t = new TChain(name.c_str(),name.c_str());
  int nfiles;//WHY????
  if (NoIso) {
    if (!gen) {
      if (madgraph) nfiles = t->Add(File_MadGraphNoIso_reco.c_str());
      else nfiles = t->Add(File_PowhegNoIso_reco.c_str());
    } else {
      if (madgraph) nfiles = t->Add(File_MadGraphNoIso_gen.c_str());
      else nfiles = t->Add(File_PowhegNoIso_gen.c_str());
    }
  } else {
    if (!gen) {
      if (madgraph) nfiles = t->Add(File_MadGraph_reco.c_str());
      else nfiles = t->Add(File_Powheg_reco.c_str());
    } else {
      if (madgraph) nfiles = t->Add(File_MadGraph_gen.c_str());
      else nfiles = t->Add(File_Powheg_gen.c_str());
    }
    
  }
  cout<<"test 1"<<nfiles<<endl;
  TBranch *b_truth=t->GetBranch("truth");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  if (elec==1) l_phistar_true=b_truth->GetLeaf("z_phistar_born");
  if (elec==2) l_phistar_true=b_truth->GetLeaf("z_phistar_naked");
  //TLeaf *l_y=b_truth->GetLeaf("z_y");//maybe include later
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

void NTupleEventEff(int elec=0){
    TH1D *MGNoIso_reco= GetPhistar(1, 0, elec, true);
    TH1D *MGNoIso_gen = GetPhistar(1, 1, elec, true);
    TH1D *PHNoIso_reco= GetPhistar(0, 0, elec, true);
    TH1D *PHNoIso_gen = GetPhistar(0, 1, elec, true);

    cout<<"done reading"<<endl; 
    TH1D *MGNoIso_eff=(TH1D*)MGNoIso_reco->Clone();
    TH1D *PHNoIso_eff=(TH1D*)PHNoIso_reco->Clone();
    cout<<"done cloning"<<endl; 
    MGNoIso_eff->Divide(MGNoIso_reco,MGNoIso_gen,1.,1.,"B");
    PHNoIso_eff->Divide(PHNoIso_reco,PHNoIso_gen,1.,1.,"B");
    cout<<"done dividing"<<endl; 

    TCanvas* EfficiencyNoIso = new TCanvas("EfficiencyNoIso","EfficiencyNoIso",800,900);
    EfficiencyNoIso->cd();
    EfficiencyNoIso->SetLogx();
    MGNoIso_eff->GetXaxis()->SetRangeUser(0.001,3.2);
    MGNoIso_eff->GetXaxis()->SetTitle("#phi*_{generated}");
    MGNoIso_eff->GetXaxis()->SetTitleOffset(0.8);
    MGNoIso_eff->GetXaxis()->SetTitleSize(0.04);
    MGNoIso_eff->GetXaxis()->SetLabelOffset(-0.01);
    MGNoIso_eff->GetXaxis()->SetLabelSize(0.04);
    MGNoIso_eff->GetYaxis()->SetTitle("N_{reconstructed}/N_{generated}");
    MGNoIso_eff->GetYaxis()->SetTitleOffset(1.2);
    MGNoIso_eff->GetYaxis()->SetTitleSize(0.04);
    MGNoIso_eff->GetYaxis()->SetLabelSize(0.04);
    //    MGNoIso_eff->GetYaxis()->SetRangeUser(0.5,0.57);
    MGNoIso_eff->SetStats(0);
    MGNoIso_eff->SetBit( TH1::kNoTitle, true );
    MGNoIso_eff->SetLineColor(1);
    MGNoIso_eff->SetMarkerColor(1);
    MGNoIso_eff->SetMarkerStyle(20);
    MGNoIso_eff->Draw();
    PHNoIso_eff->SetMarkerStyle(21);
    PHNoIso_eff->SetLineColor(2);
    PHNoIso_eff->SetMarkerColor(2);
    PHNoIso_eff->Draw("same");
    MGNoIso_eff->Draw("same");
    TLegend* legNoIso_eff = new TLegend(0.45,0.77,0.85,0.91);
    legNoIso_eff->SetFillStyle(0);
    legNoIso_eff->SetBorderSize(0);
    legNoIso_eff->SetLineWidth(1);
    legNoIso_eff->SetNColumns(1);
    legNoIso_eff->SetTextFont(42);
    legNoIso_eff->SetTextSize(0.04);
    legNoIso_eff->AddEntry(MGNoIso_eff,"MadGraph","PL");
    legNoIso_eff->AddEntry(PHNoIso_eff,"Powheg","PL");
    legNoIso_eff->Draw();
    EfficiencyNoIso->Print("/home/user1/lesko/work/Phistar/RandomwMacros/RootFiles/EffNoIso.pdf");
    EfficiencyNoIso->Print("/home/user1/lesko/work/Phistar/RandomwMacros/RootFiles/EffNoIso.png");
    
  //SECOND PART
  
    TH1D *MG_reco= GetPhistar(1, 0, elec, false);
    TH1D *MG_gen = GetPhistar(1, 1, elec, false);
    TH1D *PH_reco= GetPhistar(0, 0, elec, false);
    TH1D *PH_gen = GetPhistar(0, 1, elec, false);

    cout<<"done reading"<<endl; 
    TH1D *MG_eff=(TH1D*)MG_reco->Clone();
    TH1D *PH_eff=(TH1D*)PH_reco->Clone();
    cout<<"done cloning"<<endl; 
    MG_eff->Divide(MG_reco,MG_gen,1.,1.,"B");
    PH_eff->Divide(PH_reco,PH_gen,1.,1.,"B");
    cout<<"done dividing"<<endl; 

    TCanvas* EfficiencyWithIso = new TCanvas("EfficiencyWithIso","EfficiencyWithIso",800,900);
    EfficiencyWithIso->cd();
    EfficiencyWithIso->SetLogx();
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
    TLegend* legWithIso_eff = new TLegend(0.45,0.77,0.85,0.91);
    legWithIso_eff->SetFillStyle(0);
    legWithIso_eff->SetBorderSize(0);
    legWithIso_eff->SetLineWidth(1);
    legWithIso_eff->SetNColumns(1);
    legWithIso_eff->SetTextFont(42);
    legWithIso_eff->SetTextSize(0.04);
    legWithIso_eff->AddEntry(MG_eff,"MadGraph","PL");
    legWithIso_eff->AddEntry(PH_eff,"Powheg","PL");
    legWithIso_eff->Draw();
    EfficiencyWithIso->Print("/home/user1/lesko/work/Phistar/RandomwMacros/RootFiles/EffWithIso.pdf");
    EfficiencyWithIso->Print("/home/user1/lesko/work/Phistar/RandomwMacros/RootFiles/EffWithIso.png");
    TH1D * NoIsoRatioTH1;//=new TH1D("NoIsophistar","NoIsophistar",nphistar,phistarBins);
    TH1D * WithIsoRatioTH1;//=new TH1D("NoIsophistar","NoIsophistar",nphistar,phistarBins);

    NoIsoRatioTH1=(TH1D*)MGNoIso_eff->Clone("MyCloneNoIso");
    WithIsoRatioTH1=(TH1D*)MG_eff->Clone("MyCloneWITHIso");
    
  NoIsoRatioTH1->Divide( PHNoIso_eff);
  WithIsoRatioTH1->Divide(PH_eff);

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 900);
  c1->cd();
  c1->SetLogx();
  WithIsoRatioTH1->GetXaxis()->SetRangeUser(0.001, 3.2);
  WithIsoRatioTH1->GetXaxis()->SetTitle("#phi*_{generated}");
  WithIsoRatioTH1->GetYaxis()->SetTitle("MadGraphEff/PowHeggEff");
  WithIsoRatioTH1->SetLineColor(1);
  WithIsoRatioTH1->SetMarkerColor(1);
  WithIsoRatioTH1->SetMarkerStyle(20);
  WithIsoRatioTH1->Draw();
  NoIsoRatioTH1->SetMarkerStyle(21);
  NoIsoRatioTH1->SetLineColor(2);
  NoIsoRatioTH1->SetMarkerColor(2);
  NoIsoRatioTH1->Draw("same");


  TLegend* leg = new TLegend(0.45, 0.77, 0.85, 0.91);
  leg->AddEntry(WithIsoRatioTH1, "With Iso Cuts", "PL");
  leg->AddEntry(NoIsoRatioTH1, "No Iso Cuts", "PL");
  leg->Draw();
  c1->Print("/home/user1/lesko/work/Phistar/RandomwMacros/RootFiles/RatioPlotsINeed.pdf");
  c1->Print("/home/user1/lesko/work/Phistar/RandomwMacros/RootFiles/RatioPlotsINeed.png");
  return; 
}
