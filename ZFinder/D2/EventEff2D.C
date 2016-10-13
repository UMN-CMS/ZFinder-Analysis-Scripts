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
#include <iomanip>

#include <iostream>
#include <fstream>

const double phistarBins[] = {0.000,0.004,0.008,0.012,0.016,0.020,0.024,0.029,0.034,0.039,0.045,0.052,0.057,0.064,0.072,0.081,0.091,0.102,0.114,0.128,0.145,0.165,0.189,0.219,0.258,0.312,0.391,0.524,0.695,0.918,1.153,1.496,1.947,2.522,3.277,5.0,10.0};
size_t nphistar=(sizeof(phistarBins)/sizeof(phistarBins[0]))-1;
const double yBins[] = {0.0,0.4,0.8,1.2,1.6,2.0,2.4};
size_t ny=(sizeof(yBins)/sizeof(yBins[0]))-1;
size_t nbins=nphistar*ny;

std::string File_Signal_reco_b="/afs/cern.ch/work/r/ruckstuh/public/MG_born_reco.root";
std::string File_Powheg_reco_b="/afs/cern.ch/work/r/ruckstuh/public/PH_born_reco.root";
std::string File_Signal_gen_b="/afs/cern.ch/work/r/ruckstuh/public/MG_born_gen.root";
std::string File_Powheg_gen_b="/afs/cern.ch/work/r/ruckstuh/public/PH_born_gen.root";
std::string File_Signal_reco_n="/afs/cern.ch/work/r/ruckstuh/public/MG_bare_reco.root";
std::string File_Powheg_reco_n="/afs/cern.ch/work/r/ruckstuh/public/PH_bare_reco.root";
std::string File_Signal_gen_n="/afs/cern.ch/work/r/ruckstuh/public/MG_bare_gen.root";
std::string File_Powheg_gen_n="/afs/cern.ch/work/r/ruckstuh/public/PH_bare_gen.root";
std::string File_Signal_reco_d="/afs/cern.ch/work/r/ruckstuh/public/MG_dressed_reco.root";
std::string File_Powheg_reco_d="/afs/cern.ch/work/r/ruckstuh/public/PH_dressed_reco.root";
std::string File_Signal_gen_d="/afs/cern.ch/work/r/ruckstuh/public/MG_dressed_gen.root";
std::string File_Powheg_gen_d="/afs/cern.ch/work/r/ruckstuh/public/PH_dressed_gen.root";

std::string reco_name="Combined Single Reco";
std::string gen_name ="Combined Gen Cuts Reco";

vector<TH1D*> SplitHist(TH1D* hist){
  vector<TH1D*> v;
  for (uint i=0; i<ny; i++){
    TH1D* h=new TH1D("phistar","phistar",nphistar,phistarBins);
    h->Sumw2();
    for (uint j=0; j<nphistar; j++){
      int bin=i*nphistar+j;
      h->SetBinContent(j+1,hist->GetBinContent(bin+1));
      h->SetBinError(j+1,hist->GetBinError(bin+1));
    }
    v.push_back(h);
  }
  return v;
}

int GetBin(double phistar, double y){
  TAxis* A_phistar=new TAxis(nphistar,phistarBins);
  TAxis* A_y=new TAxis(ny,yBins);
  int bin_phistar=A_phistar->FindBin(phistar)-1;
  int bin_y=A_y->FindBin(fabs(y))-1;
  int bin=bin_y*nphistar+bin_phistar;
  if (bin_phistar<0 || bin_y<0) bin=-1;
  if (bin_phistar>=nphistar || bin_y>=ny) bin=nbins;
  return bin;
}

TH1D* GetPhistar(bool madgraph, bool gen, int elec=0, bool ref=0, TH1D *rew_m=0, TH1D *rew_p=0){
  std::string name=reco_name;
  if (gen) name=gen_name;
  TChain* t = new TChain(name.c_str(),name.c_str());
  int nfiles;
  if (!gen){
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
  }
  else{
    if (madgraph) {
      if (elec==0) nfiles=t->Add(File_Signal_gen_d.c_str());
      if (elec==1) nfiles=t->Add(File_Signal_gen_b.c_str());
      if (elec==2) nfiles=t->Add(File_Signal_gen_n.c_str());
    }
    else {
      if (elec==0) nfiles=t->Add(File_Powheg_gen_d.c_str());
      if (elec==1) nfiles=t->Add(File_Powheg_gen_b.c_str());
      if (elec==2) nfiles=t->Add(File_Powheg_gen_n.c_str());
    }
  }
  TBranch *b_truth=t->GetBranch("truth");
  TBranch *b_reco=t->GetBranch("reco");
  TLeaf *l_phistar_true=b_truth->GetLeaf("z_phistar_dressed");
  if (elec==1) l_phistar_true=b_truth->GetLeaf("z_phistar_born");
  if (elec==2) l_phistar_true=b_truth->GetLeaf("z_phistar_naked");
  TLeaf *l_y_true=b_truth->GetLeaf("z_y");
  if (elec==1) l_y_true=b_reco->GetLeaf("z_yBorn");
  if (elec==2) l_y_true=b_reco->GetLeaf("z_yNaked");
  TLeaf *l_pt=b_truth->GetLeaf("z_pt");
  int nweights;
  t->SetBranchAddress("weight_size",&nweights);
  t->GetEntry(0);
  cout<<"The sample has nweights: "<<nweights<<endl;
  double weights[nweights];
  int weightid[nweights];
  t->SetBranchAddress("weights",&weights);
  t->SetBranchAddress("weight_ids",&weightid);

  TH1D* ph;
  if (!ref) ph=new TH1D("bin","bin",nbins,0,nbins);
  else ph=new TH1D("ref","ref",100000,0,1000);
  ph->Sumw2();
  cout<<gen<<" "<<madgraph<<"  "<<t->GetEntries()<<endl;
  for (int i=0; i<t->GetEntries();i++){
    t->GetEntry(i);
    //if (fabs(l_y->GetValue())<1.5)continue;
    double weight =1;
    for (int w=0; w<nweights;w++){
      if (!gen){
	if (weightid[w]==1 || weightid[w]==2 || weightid[w]==12 || weightid[w]==13 || weightid[w]==20 || weightid[w]==30) {weight=weight*weights[w];}
      }
      else{
	if (weightid[w]==1 || weightid[w]==2) {weight=weight*weights[w];}
      }
    }
    double phistar_true=l_phistar_true->GetValue();
    double y_true=l_y_true->GetValue();
    int bin_true=GetBin(phistar_true,y_true);
    double zpt=l_pt->GetValue();
    if (ref){ph->Fill(zpt,weight);}
    // if (!ref && madgraph) {
    //   int idx=rew_m->GetXaxis()->FindBin(zpt);
    //   if (rew_m->GetBinContent(idx)!=0) {
    // 	   	weight=weight*rew_p->GetBinContent(idx)/rew_m->GetBinContent(idx);
    // 	//	weight=weight;
    //   }
    //   else {
    // 	weight=0;
    // 	cout<<rew_p->GetBinContent(idx)<<" "<<rew_m->GetBinContent(idx)<<"  "<<zpt<<"  "<<idx<<endl;
    //   }
    // }
    if (!ref) {
      ph->Fill(bin_true,weight);
      // if (i%100000==0) cout<<phistar_true<<" "<<weight<<endl;
    }
  }
  return ph;
}

void NTupleEventEff(){
  for (int elec=1; elec<2; elec++){
    // TH1D *MG_ref = GetPhistar(1, 1, elec, 1);
    // TH1D *PH_ref = GetPhistar(0, 1, elec, 1);
    // TH1D *MG_reco= GetPhistar(1, 0, elec, 0, MG_ref, PH_ref);
    // TH1D *MG_gen = GetPhistar(1, 1, elec, 0, MG_ref, PH_ref);
    TH1D *MG_reco= GetPhistar(1, 0, elec, 0);
    TH1D *MG_gen = GetPhistar(1, 1, elec, 0);
    TH1D *PH_reco= GetPhistar(0, 0, elec, 0);
    TH1D *PH_gen = GetPhistar(0, 1, elec, 0);

    cout<<"done reading"<<endl; 
    TH1D *MG_eff=(TH1D*)MG_reco->Clone();
    TH1D *PH_eff=(TH1D*)PH_reco->Clone();
    cout<<"done cloning"<<endl; 
    MG_eff->Divide(MG_reco,MG_gen,1.,1.,"B");
    PH_eff->Divide(PH_reco,PH_gen,1.,1.,"B");
    cout<<"done dividing"<<endl; 

    std::string plotname = "Plots/EventEff_";
    if (elec==0) plotname+="Dressed_";
    if (elec==1) plotname+="Born_";
    if (elec==2) plotname+="Naked_";

    TH1D *ratio=(TH1D*)MG_reco->Clone();
    ratio->Divide(MG_eff,PH_eff,1.,1.);
    TCanvas* Efficiency = new TCanvas("Efficiency","Efficiency",1600,900);
    Efficiency->cd();
    //   Efficiency->SetLogx();
    //    MG_eff->GetXaxis()->SetRangeUser(0.001,3.2);
    MG_eff->GetXaxis()->SetTitle("(#phi*,y)-bin_{generated}");
    MG_eff->GetXaxis()->SetTitleOffset(0.8);
    MG_eff->GetXaxis()->SetTitleSize(0.04);
    MG_eff->GetXaxis()->SetLabelOffset(0);
    MG_eff->GetXaxis()->SetLabelSize(0.03);
    MG_eff->GetYaxis()->SetTitle("N_{reconstructed}/N_{generated}");
    MG_eff->GetYaxis()->SetTitleOffset(1.2);
    MG_eff->GetYaxis()->SetTitleSize(0.04);
    MG_eff->GetYaxis()->SetLabelSize(0.03);
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
    Efficiency->SaveAs((plotname+"AllBins.pdf").c_str());


    TCanvas* Correction = new TCanvas("Correction","Correction",1600,900);
    Correction->cd();
    // Correction->SetLogx();
    // ratio->GetXaxis()->SetRangeUser(0.001,3.2);
    ratio->GetXaxis()->SetTitle("(#phi*,y)-bin_{generated}");
    ratio->GetXaxis()->SetTitleOffset(0.8);
    ratio->GetXaxis()->SetTitleSize(0.04);
    ratio->GetXaxis()->SetLabelOffset(0);
    ratio->GetXaxis()->SetLabelSize(0.03);
    ratio->GetYaxis()->SetTitle("MadGraph/Powheg");
    ratio->GetYaxis()->SetTitleOffset(1.2);
    ratio->GetYaxis()->SetTitleSize(0.04);
    ratio->GetYaxis()->SetLabelSize(0.03);
    //    ratio->GetYaxis()->SetRangeUser(0.5,0.57);
    ratio->SetStats(0);
    ratio->SetBit( TH1::kNoTitle, true );
    ratio->SetLineColor(1);
    ratio->SetMarkerColor(1);
    ratio->SetMarkerStyle(20);
    ratio->Draw();
    ratio->SaveAs((plotname+"Ratio_AllBins.pdf").c_str());

    vector<TH1D*> v_MG_eff=SplitHist(MG_eff);
    vector<TH1D*> v_PH_eff=SplitHist(PH_eff);
    vector<TH1D*> v_ratio=SplitHist(ratio);

    for (int i=0; i<ny; i++){ 
      std::ostringstream strs;
      strs << i;
      std::string Canvasname="EventEff_Bin"+strs.str();
      std::string plotname2=plotname+"_Bin"+strs.str();
      TCanvas* FinalPhiTot = new TCanvas(Canvasname.c_str(),Canvasname.c_str(),800,900);
      FinalPhiTot->Divide(1,2);
      FinalPhiTot->cd(1);
      gPad->SetPad("p1","p1",0,2.5/9.0, 1,1,kWhite,0,0);
      gPad->SetBottomMargin(0.01);
      gPad->SetTopMargin(0.06);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.06);
      gPad->SetLogx(1);
      v_MG_eff[i]->GetXaxis()->SetRangeUser(0.001,10);
      v_MG_eff[i]->GetXaxis()->SetTitleOffset(1.05); 
      v_MG_eff[i]->GetXaxis()->SetTitle(0);
      v_MG_eff[i]->GetXaxis()->SetTitleSize(0.05);
      v_MG_eff[i]->GetXaxis()->SetLabelSize(0.05);
      v_MG_eff[i]->GetYaxis()->SetTitleOffset(1.05); 
      v_MG_eff[i]->GetYaxis()->SetTitleSize(0.05);
      v_MG_eff[i]->GetYaxis()->SetLabelSize(0.04);
      //      v_MG_eff[i]->GetXaxis()->SetRangeUser(0.001,3.0);
      //      v_MG_eff[i]->GetXaxis()->SetTitleOffset(0.8);
      //      v_MG_eff[i]->GetXaxis()->SetTitleSize(0.04);
      //      v_MG_eff[i]->GetXaxis()->SetLabelOffset(0);
      //      v_MG_eff[i]->GetXaxis()->SetLabelSize(0.03);
      v_MG_eff[i]->GetYaxis()->SetTitle("N_{reconstructed}/N_{generated}");
      //      v_MG_eff[i]->GetYaxis()->SetTitleOffset(1.2);
      //      v_MG_eff[i]->GetYaxis()->SetTitleSize(0.04);
      //      v_MG_eff[i]->GetYaxis()->SetLabelSize(0.03);
      if (i==0)  v_MG_eff[i]->GetYaxis()->SetRangeUser(0.53,0.65);
      if (i==1)  v_MG_eff[i]->GetYaxis()->SetRangeUser(0.49,0.60);
      if (i==2)  v_MG_eff[i]->GetYaxis()->SetRangeUser(0.43,0.55);
      if (i==3)  v_MG_eff[i]->GetYaxis()->SetRangeUser(0.38,0.55);
      if (i==4)  v_MG_eff[i]->GetYaxis()->SetRangeUser(0.30,0.55);
      if (i==5)  v_MG_eff[i]->GetYaxis()->SetRangeUser(0.26,0.50);
      v_MG_eff[i]->SetStats(0);
      v_MG_eff[i]->SetBit( TH1::kNoTitle, true );
      v_MG_eff[i]->SetLineColor(1);
      v_MG_eff[i]->SetMarkerColor(1);
      v_MG_eff[i]->SetMarkerStyle(20);
      v_MG_eff[i]->Draw();
      v_PH_eff[i]->SetMarkerStyle(21);
      v_PH_eff[i]->SetLineColor(2);
      v_PH_eff[i]->SetMarkerColor(2);
      v_PH_eff[i]->Draw("same");
      v_MG_eff[i]->Draw("same");
      TLegend* leg= new TLegend(0.65,0.83,0.94,0.93);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->SetLineWidth(1);
      leg->SetNColumns(1);
      leg->SetTextFont(42);
      leg->AddEntry(MG_eff,"MadGraph","PL");
      leg->AddEntry(PH_eff,"Powheg","PL");
      leg->Draw();

      TLatex mark;
      mark.SetTextSize(0.043);
      mark.SetTextFont(42);
      //mark.SetTextSize(0.035);
      mark.SetNDC(true);
      std::ostringstream ystrs;
      ystrs << std::fixed << std::setprecision(1) << yBins[i];
      ystrs<<" #leq |y_{ee}| < ";
      ystrs<<std::setprecision(1) <<yBins[i+1];
      mark.DrawLatex(0.19,0.88,ystrs.str().c_str());
      FinalPhiTot->cd(2);
      gPad->SetPad("p2","p2",0,0,1,2.5/9.0,kWhite,0,0);
      gPad->SetBottomMargin(0.37);
      gPad->SetTopMargin(0.01);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.06);
      gPad->SetLogx(1);
      v_ratio[i]->GetXaxis()->SetRangeUser(0.001,10);
      v_ratio[i]->GetXaxis()->SetTitleOffset(1.05);
      //      v_ratio[i]->GetXaxis()->SetTitleSize(0.12);
      v_ratio[i]->GetXaxis()->SetLabelSize(0.12);
      v_ratio[i]->GetXaxis()->SetTitleSize(0.15);
      v_ratio[i]->GetXaxis()->CenterTitle();
      v_ratio[i]->GetXaxis()->SetTitle("#phi*(generated)");
      //      v_ratio[i]->GetXaxis()->SetTitleOffset(0.8);
      //      v_ratio[i]->GetXaxis()->SetTitleSize(0.04);
      v_ratio[i]->GetXaxis()->SetLabelOffset(0);
      //      v_ratio[i]->GetXaxis()->SetLabelSize(0.03);
     
      v_ratio[i]->GetYaxis()->SetRangeUser(0.8,1.2);
      //      v_ratio[i]->GetYaxis()->SetTitleOffset(0.32);
      v_ratio[i]->GetYaxis()->SetTitleSize(0.1);
      v_ratio[i]->GetYaxis()->SetLabelSize(0.1);
      v_ratio[i]->GetYaxis()->SetNdivisions(3,5,0);
      v_ratio[i]->GetYaxis()->SetTitleOffset(0.45);
      v_ratio[i]->GetYaxis()->CenterTitle();
      v_ratio[i]->GetYaxis()->SetTitle("MG/PH");
      if (i==0)  v_ratio[i]->GetYaxis()->SetRangeUser(0.94,1.06);
      if (i==1)  v_ratio[i]->GetYaxis()->SetRangeUser(0.94,1.06);
      if (i==2)  v_ratio[i]->GetYaxis()->SetRangeUser(0.94,1.06);
      if (i==3)  v_ratio[i]->GetYaxis()->SetRangeUser(0.88,1.12);
      if (i==4)  v_ratio[i]->GetYaxis()->SetRangeUser(0.76,1.24);
      if (i==5)  v_ratio[i]->GetYaxis()->SetRangeUser(0.41,1.59);
       //      v_ratio[i]->GetYaxis()->SetTitleOffset(1.2);
      //      v_ratio[i]->GetYaxis()->SetTitleSize(0.04);
      //      v_ratio[i]->GetYaxis()->SetLabelSize(0.03);
      v_ratio[i]->SetStats(0);
      v_ratio[i]->SetBit( TH1::kNoTitle, true );
      v_ratio[i]->SetLineColor(1);
      v_ratio[i]->SetMarkerColor(1);
      v_ratio[i]->SetMarkerStyle(20);
      v_ratio[i]->Draw();

      FinalPhiTot->SaveAs((plotname2+".pdf").c_str());
    }      
  }
  return; 
}
