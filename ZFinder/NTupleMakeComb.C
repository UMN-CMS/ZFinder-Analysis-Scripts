#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
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
#include <sstream>
#include <string>

#include <iomanip>
#include <iostream>
#include <fstream>

const bool doNorm=1;
const bool doMG=0;
const std::string Tag="";

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277};
//const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;

TH1D* GetHistMuon(std::string name, bool gen=0){
  TH1D* h_Muon=new TH1D("h_Comb","h_Comb",nphistar,phistarBins);
  h_Muon->Sumw2();
  TFile fm(name.c_str());
  TH1F *h_temp;
  if (!gen) h_temp= (TH1F*)fm.Get("hRecoClone1");
  else h_temp= (TH1F*)fm.Get("Gen_Dressed_phistar_NoPU");
  for (uint i=0; i<nphistar; i++){
    h_Muon->SetBinContent(i+1,h_temp->GetBinContent(i+1));
    h_Muon->SetBinError(i+1,h_temp->GetBinError(i+1));
  }
  return h_Muon;
} 

TH1D* GetHistElec(std::string name, std::string histname){
  TH1D* h_Elec=new TH1D("h_Elec","h_Elec",nphistar,phistarBins);
  h_Elec->Sumw2();
  TFile fm(name.c_str());
  TH1D *h_temp = (TH1D*)fm.Get(histname.c_str());
  for (uint i=0; i<nphistar; i++){
    h_Elec->SetBinContent(i+1,h_temp->GetBinContent(i+1));
    h_Elec->SetBinError(i+1,h_temp->GetBinError(i+1));
  }
  return h_Elec;
} 

TH1D* ConvertToHist(TGraphAsymmErrors* g){
  TH1D* h_temp;
  h_temp=new TH1D("h_tot","h_tot",nphistar,phistarBins);
  h_temp->Sumw2();
  for (size_t iphistar=0; iphistar<nphistar; iphistar++){
    double x,y;
    g->GetPoint(iphistar,x,y);
    double error =g->GetErrorYhigh(iphistar);
    h_temp->SetBinContent(iphistar+1,y);
    h_temp->SetBinError(iphistar+1,error);
  }
  return h_temp;
}

void NTupleMakeComb(){
  std::string muon_name="Muon_CrossSection_For_Combination/Madgraph_Dressed_";
  if (doNorm) muon_name+="Normal_";
  else        muon_name+="Absolute_";
  std::string muon_lep=muon_name+"Central_UnCorrelated_Errors_No_BGND.root";
  std::string muon_bkg=muon_name+"Central_UnCorrelated_Errors_Only_BGND.root";
  // if (doNorm){
  //  muon_lep=muon_name+"UnCorrelated_Errors_No_BGND.root";
  //  muon_bkg=muon_name+"UnCorrelated_Errors_Only_BGND.root";
  // }
  std::string muon_pup=muon_name+"PU_UP.root";
  std::string muon_pum=muon_name+"PU_DOWN.root";

  TH1D* h_Muon_lep=GetHistMuon(muon_lep);
  TH1D* h_Muon_bkg=GetHistMuon(muon_bkg);
  TH1D* h_Muon_pup=GetHistMuon(muon_pup);
  TH1D* h_Muon_pum=GetHistMuon(muon_pum);

  if (!doNorm){
    h_Muon_bkg->Scale(1./1000.);
    h_Muon_pup->Scale(1./1000.);
    h_Muon_pum->Scale(1./1000.);
    h_Muon_lep->Scale(1./1000.);
  }
  std::string elec_name="Final_Hist_";
  elec_name+=Tag;
  if (doNorm) elec_name+="Norm_";
  else        elec_name+="Abs_";
  elec_name+="MG_";
  elec_name+="Dressed.root";

  TH1D* h_Elec_lep=GetHistElec(elec_name,"h_data_elec");
  TH1D* h_Elec_bkg=GetHistElec(elec_name,"h_data_bgnd");
  TH1D* h_Elec_pup=GetHistElec(elec_name,"h_data_pup");
  TH1D* h_Elec_pum=GetHistElec(elec_name,"h_data_pum");

  TH1D* h_Data_lep=(TH1D*)h_Muon_lep->Clone();
  TH1D* h_Data_pup=(TH1D*)h_Muon_pup->Clone();
  TH1D* h_Data_pum=(TH1D*)h_Muon_pum->Clone();

  h_Data_lep->Add(h_Elec_lep,1);
  h_Data_pup->Add(h_Elec_pup,1);
  h_Data_pum->Add(h_Elec_pum,1);
  for (uint i=0; i<nphistar; i++){
    double error_elec=h_Elec_lep->GetBinError(i+1)/h_Data_lep->GetBinContent(i+1);
    double error_muon=h_Muon_lep->GetBinError(i+1)/h_Data_lep->GetBinContent(i+1);
    double error_lep=sqrt(error_elec*error_elec+error_muon*error_muon);
    double error_bkg=(h_Elec_bkg->GetBinError(i+1)+h_Muon_bkg->GetBinError(i+1))/h_Data_lep->GetBinContent(i+1);
    double error_lumi=0;
    if (!doNorm) error_lumi=0.026;
    double error_tot=sqrt(error_lep*error_lep+error_bkg*error_bkg+error_lumi*error_lumi)*h_Data_lep->GetBinContent(i+1);
    h_Data_lep->SetBinError(i+1,error_tot);
  }
  h_Data_pup->Scale(0.5);
  h_Data_pum->Scale(0.5);
  h_Data_lep->Scale(0.5);

  for (uint i=0; i<nphistar; i++){
    //    double error_fsr=(h_Data_fsr->GetBinContent(i+1)-h_Data_lep->GetBinContent(i+1))/h_Data_lep->GetBinContent(i+1);
    double error_pum=fabs((h_Data_pum->GetBinContent(i+1)-h_Data_lep->GetBinContent(i+1))/h_Data_lep->GetBinContent(i+1));
    double error_tot=h_Data_lep->GetBinError(i+1)/h_Data_lep->GetBinContent(i+1);
    error_tot=sqrt(error_tot*error_tot+error_pum*error_pum)*h_Data_lep->GetBinContent(i+1);
    h_Data_lep->SetBinError(i+1,error_tot);
  }

  std::string textn="Comb_Hist_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  textn+="MG_";
  textn+="Dressed.root";
  TFile tr(textn.c_str(),"RECREATE");
  h_Data_lep->Write();
}

void NTupleMakeCombMC(){
  std::string muon_name="Muon_CrossSection_For_Combination/Gen_";
  if (doMG)   muon_name+="Madgraph_";
  else        muon_name+="Powheg_";
  muon_name+="Dressed_";
  if (doNorm) muon_name+="Normal_";
  else        muon_name+="Absolute_";
  muon_name+="Central_Stat_Errors.root";

  cout<<muon_name<<endl;
  TH1D* h_Muon=GetHistMuon(muon_name);

  if (!doNorm){
    h_Muon->Scale(1./1000.);
  }

  std::string textn="Final_Hist_MC_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  textn+="Dressed.root";
  TH1D* h_Elec=GetHistElec(textn,"h_mc");

  textn="Data_Graph_MC_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  textn+="Dressed.root";
  TFile fm(textn.c_str());  
  TGraphAsymmErrors* g_temp = (TGraphAsymmErrors*)fm.Get("Graph");
  TH1D* h_Elec_Tot=ConvertToHist(g_temp);
 
  TH1D* h_MC=(TH1D*)h_Muon->Clone();
  h_MC->Add(h_Elec,1);
  for (uint i=0; i<nphistar; i++){
    double error_elec=h_Elec->GetBinError(i+1)/h_MC->GetBinContent(i+1);
    double error_muon=h_Muon->GetBinError(i+1)/h_MC->GetBinContent(i+1);
    double error_elec2=sqrt(h_Elec_Tot->GetBinError(i+1)*h_Elec_Tot->GetBinError(i+1) - h_Elec->GetBinError(i+1)*h_Elec->GetBinError(i+1))/h_Elec->GetBinContent(i+1);
    // double error_tot=sqrt(error_elec*error_elec+error_muon*error_muon+error_elec2*error_elec2)*h_MC->GetBinContent(i+1);
    double error_tot=sqrt(error_elec*error_elec+error_elec*error_elec+error_elec2*error_elec2)*h_MC->GetBinContent(i+1);
    // cout<<error_elec<<" "<<error_muon<<" "<<error_elec2<<" "<<error_tot<<endl;
    cout<<h_Elec->GetBinError(i+1)<<" "<<h_Muon->GetBinError(i+1)<<endl;
    h_MC->SetBinError(i+1,error_tot);
  }
  h_MC->Scale(0.5);
  textn="Comb_Hist_MC_";
  textn+=Tag;
  if (doNorm) textn+="Norm_";
  else        textn+="Abs_";
  if (doMG)   textn+="MG_";
  else        textn+="PH_";
  textn+="Dressed.root";
  TFile tr(textn.c_str(),"RECREATE");
  h_MC->Write();
}
void MakeComp(){
  std::string muon_name="Muon_CrossSection_For_Combination/Madgraph_Dressed_";
  if (doNorm) muon_name+="Normal_";
  else        muon_name+="Absolute_";
  std::string muon_lep=muon_name+"Central_UnCorrelated_Errors_No_BGND.root";
  cout<<"Hi1"<<endl;
  cout<<muon_lep<<endl;
  TH1D* h_Muon_lep=GetHistMuon(muon_lep);
  cout<<"Hi2"<<endl;

  if (!doNorm){
    h_Muon_lep->Scale(1./1000.);
  }
  std::string elec_name="Final_Hist_";
  elec_name+=Tag;
  if (doNorm) elec_name+="Norm_";
  else        elec_name+="Abs_";
  elec_name+="MG_";
  elec_name+="Dressed.root";
  cout<<"Hi3"<<endl;

  TH1D* h_Elec_lep=GetHistElec(elec_name,"h_data_elec");

  cout<<"Hi"<<endl;
  TH1D* h_Data_lep=(TH1D*)h_Muon_lep->Clone();
  h_Data_lep->Sumw2();
  h_Data_lep->Divide(h_Elec_lep);
  cout<<"Hi4"<<endl;
  h_Data_lep->Draw();
  double chi2=0;
  for (int i=0; i<nphistar;i++){
    chi2+=pow((h_Data_lep->GetBinContent(i+1)-1)/h_Data_lep->GetBinError(i+1),2);
    cout<<h_Data_lep->GetBinContent(i+1)<<" "<<h_Data_lep->GetBinError(i+1)<<" "<<pow((h_Data_lep->GetBinContent(i+1)-1)/h_Data_lep->GetBinError(i+1),2)<<endl;
  }
  cout<<chi2/nphistar<<endl;
}
