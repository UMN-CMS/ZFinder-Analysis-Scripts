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
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
RooBinning phistarbin(nphistar,phistarBins);

vector<std::string> Get_File_Data(){
  vector<std::string> vec;
  vec.push_back("/afs/cern.ch/work/r/ruckstuh/public/Data.root");
  return vec;
}
std::string File_Signal_reco="/afs/cern.ch/work/r/ruckstuh/public/MadGraph_reco.root";
std::string File_Signal_gen= "/afs/cern.ch/work/r/ruckstuh/public/MadGraph_gen.root";
std::string File_Powheg_reco="/afs/cern.ch/work/r/ruckstuh/public/Powheg_reco.root";
std::string File_Powheg_gen= "/afs/cern.ch/work/r/ruckstuh/public/Powheg_gen.root";
std::string File_tt=         "/afs/cern.ch/work/r/ruckstuh/public/TTbar.root";
std::string File_tautau=     "/afs/cern.ch/work/r/ruckstuh/public/TauTau.root";
std::string File_tbarw=      "/afs/cern.ch/work/r/ruckstuh/public/TbarW.root";
std::string File_tw=         "/afs/cern.ch/work/r/ruckstuh/public/TW.root";
std::string File_ww=         "/afs/cern.ch/work/r/ruckstuh/public/WW.root";
std::string File_wz=         "/afs/cern.ch/work/r/ruckstuh/public/WZ.root";
std::string File_zz=         "/afs/cern.ch/work/r/ruckstuh/public/ZZ.root";

std::string reco_name="Combined Single Reco";
std::string gen_name ="Combined Gen Cuts Reco";

TH1D* GetBGPhiStar(std::string FileName, double sampleweight){
  //TH1D *phistartemp= new TH1D("phistar","phistar",nphistar,phistarBins);
  TH1D *phistartemp= new TH1D("phistar","phistar",160,50,130);
  phistartemp->Sumw2();

  gErrorIgnoreLevel = kError;
  cout<<"reading data for "<<FileName<<endl;
  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  nfiles=t->Add(FileName.c_str());

  TBranch *b_reco=t->GetBranch("reco");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_m=b_reco->GetLeaf("z_m");
  TLeaf *l_e0_q =b_reco->GetLeaf("e_charge0");
  TLeaf *l_e1_q =b_reco->GetLeaf("e_charge1");
  int nweights;
  t->SetBranchAddress("weight_size",&nweights);
  t->GetEntry(0);
  cout<<"The sample has nweights: "<<nweights<<endl;
  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);

  cout<<"reading data for "<<FileName<<endl;
  cout<<"Entries: "<<t->GetEntries()<<endl;
  for (int i=0; i<t->GetEntries();i++){
    t->GetEntry(i);
    double E0_q =l_e0_q->GetValue();
    double E1_q =l_e1_q->GetValue();
    if (E0_q!=E1_q) continue;
    double phistar=l_phistar->GetValue();
    double weight =sampleweight;
    for (int w=0; w<nweights;w++){
      // if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
    }
    //phistartemp->Fill(phistar,weight);
    phistartemp->Fill(l_m->GetValue(),weight);
  }
  cout<<"done reading data for "<<FileName<<endl;
  return phistartemp;
}

TH1D* GetDataPhiStar(){
  cout<<"reading data"<<endl;
  vector<std::string> File_Data=Get_File_Data();
  // TH1D *h_phistar= new TH1D("phistar","phistar",nphistar,phistarBins);
  TH1D *h_phistar= new TH1D("phistar","phistar",160,50,130);
  h_phistar->Sumw2();
  TChain* t = new TChain(reco_name.c_str(),reco_name.c_str());
  int nfiles;
  for (uint i=0; i<File_Data.size(); i++){
    nfiles=t->Add(File_Data[i].c_str());
  }
  TBranch *b_reco=t->GetBranch("reco");
  TLeaf *l_m=b_reco->GetLeaf("z_m");
  TLeaf *l_phistar=b_reco->GetLeaf("z_phistar_dressed");
  TLeaf *l_e0_q =b_reco->GetLeaf("e_charge0");
  TLeaf *l_e1_q =b_reco->GetLeaf("e_charge1");
  cout<<"Entries: "<<t->GetEntries()<<endl;
  for (int i=0; i<t->GetEntries();i++){
    t->GetEntry(i);
    double E0_q =l_e0_q->GetValue();
    double E1_q =l_e1_q->GetValue();
    if (E0_q!=E1_q) continue;
    //h_phistar->Fill(l_phistar->GetValue(),1);
    h_phistar->Fill(l_m->GetValue(),1);
  }
  cout<<"filled data phistar histogram"<<endl;
  return h_phistar;
}

void NTupleSS(){
  double Lumi = 19712.;
  double ttbar_weight  = Lumi*(23.64/4246440.);
  double tautau_weight = Lumi*(1966.7/47271600.);
  double tbarw_weight  = Lumi*(11.1/493460.);
  double tw_weight     = Lumi*(11.1/497658.);
  double ww_weight     = Lumi*(54.84/10000430.);
  double wz_weight     = Lumi*(33.21/10000280.);
  double zz_weight     = Lumi*(17.0/9799908.);
  double signal_weight = Lumi*(3531.89/30459500.);//3504
 
  TH1D* Data=GetDataPhiStar();
  TH1D* h_ee    =GetBGPhiStar(File_Signal_reco,    signal_weight);
  TH1D* h_zz    =GetBGPhiStar(File_zz,    zz_weight);
  TH1D* h_tt    =GetBGPhiStar(File_tt,    ttbar_weight);
  TH1D* h_tautau=GetBGPhiStar(File_tautau,tautau_weight);
  TH1D* h_tbarw =GetBGPhiStar(File_tbarw, tbarw_weight);
  TH1D* h_tw    =GetBGPhiStar(File_tw,    tw_weight);
  TH1D* h_ww    =GetBGPhiStar(File_ww,    ww_weight);
  TH1D* h_wz    =GetBGPhiStar(File_wz,    wz_weight);
  //  TH1D* h_zz    =GetBGPhiStar(File_zz,    zz_weight);

  double data_sel=Data->GetSumOfWeights();
  double ee_sel=h_ee->GetSumOfWeights();
  double tt_sel=h_tt->GetSumOfWeights();
  double tautau_sel=h_tautau->GetSumOfWeights();
  double tbarw_sel=h_tbarw->GetSumOfWeights();
  double tw_sel=h_tw->GetSumOfWeights();
  double ww_sel=h_ww->GetSumOfWeights();
  double wz_sel=h_wz->GetSumOfWeights();
  double zz_sel=h_zz->GetSumOfWeights();

  // double t_bg=ee_sel+tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel+wz_sel+zz_sel;
  double t_bg=ee_sel+tt_sel+tautau_sel+tbarw_sel+tw_sel+ww_sel+wz_sel+zz_sel;
  cout<<"total left:"<<data_sel-t_bg;
  cout<<"Data: "<<data_sel<<"  ee: "<<ee_sel<<"  tt: "<<tt_sel<<"  tautau: "<<tautau_sel<<"  tbarw: "<<tbarw_sel<<"  tw: "<<tw_sel<<" singletop: " <<tbarw_sel+tw_sel<<"  ww: "<<ww_sel<<"  wz: "<<wz_sel<<"  zz: "<<zz_sel<< endl;
  cout<<"ratio:"          <<"  ee: "<<ee_sel/data_sel<<"  tt: "<<tt_sel/data_sel<<"  tautau: "<<tautau_sel/data_sel<<"  tbarw: "<<tbarw_sel/data_sel<<"  tw: "<<tw_sel/data_sel<<" singletop: " <<(tbarw_sel+tw_sel)/data_sel<< "  ww: "<<ww_sel/data_sel<<"  wz: "<<wz_sel/data_sel<<"  zz: "<<zz_sel/data_sel<< "data: "<< (data_sel-t_bg)/data_sel<<endl;

  TH1D* h_qcd  =(TH1D*)Data->Clone();
  h_qcd->Add(h_ee,-1.0);
  h_qcd->Add(h_tt,-1.0);
  h_qcd->Add(h_tautau,-1.0);
  h_qcd->Add(h_tbarw,-1.0);
  h_qcd->Add(h_tw,-1.0);
  h_qcd->Add(h_ww,-1.0);
  h_qcd->Add(h_wz,-1.0);
  h_qcd->Add(h_zz,-1.0);

  h_qcd->Draw();
  cout<<"norm="<<h_qcd->GetSumOfWeights()<<endl;;
}
