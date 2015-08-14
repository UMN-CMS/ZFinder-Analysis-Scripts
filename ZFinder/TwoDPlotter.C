// Standard Library
#include <algorithm>  //std::equal_range, std::sort, std::unique
#include <cmath>  // std::abs
#include <iostream>
#include <string>

// ROOT
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TFile.h>
#include <TStyle.h>
#include "TGraphAsymmErrors.h"

const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277, 10};
const Int_t nphistar = ((sizeof (phistarBins) / sizeof (phistarBins[0])) - 1);

using namespace std;

void TwoDPlotter() {


    TFile* YDistributionFile = new TFile("Data_Graph_Norm_PH_YSeperated-0.10-0.20-0.30-0.40-0.50-0.60-0.80-1.00-1.50Dressed.root");
    TGraphAsymmErrors* AllGraph;
    AllGraph = (TGraphAsymmErrors*) YDistributionFile->Get("Graph");
    if (!AllGraph) {
        cout << "Didn't open properly " << endl;
        return;
    }
    vector <TH1D*> AllHistos;
    string Names[10];
    Names[0] = "Y_0-.5";
    AllHistos.push_back(new TH1D(Names[0].c_str(), Names[0].c_str(), nphistar, phistarBins));
    Names[1] = "Y_.5-1";
    AllHistos.push_back(new TH1D(Names[1].c_str(), Names[1].c_str(), nphistar, phistarBins));
    Names[2] = "Y_1-1.5";
    AllHistos.push_back(new TH1D(Names[2].c_str(), Names[2].c_str(), nphistar, phistarBins));
    Names[3] = "Y_1.5-inf";
    AllHistos.push_back(new TH1D(Names[3].c_str(), Names[3].c_str(), nphistar, phistarBins));
   // Names[4] = "Y_.4-.5";
   // AllHistos.push_back(new TH1D(Names[4].c_str(), Names[4].c_str(), nphistar, phistarBins));
   // Names[5] = "Y_.5-.6";
   // AllHistos.push_back(new TH1D(Names[5].c_str(), Names[5].c_str(), nphistar, phistarBins));
   // Names[6] = "Y_.6-.8";
   // AllHistos.push_back(new TH1D(Names[6].c_str(), Names[6].c_str(), nphistar, phistarBins));
   // Names[7] = "Y_.8-1";
   // AllHistos.push_back(new TH1D(Names[7].c_str(), Names[7].c_str(), nphistar, phistarBins));
   // Names[8] = "Y_1-1.5";
   // AllHistos.push_back(new TH1D(Names[8].c_str(), Names[8].c_str(), nphistar, phistarBins));
   // Names[9] = "Y_1.5-inf";
   // AllHistos.push_back(new TH1D(Names[9].c_str(), Names[9].c_str(), nphistar, phistarBins));

    for (Int_t iphistar = 0; iphistar < (AllGraph->GetN()); iphistar++) {
        double x = 0;
        double y = 0;
        if (iphistar < nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[0]->Fill(x, y);
        }
        if (iphistar > nphistar && iphistar < 2 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[0]->Fill(x - (10 + 1), y);
        }
        if (iphistar > 2 * nphistar && iphistar < 3 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[0]->Fill(x - 2 * (10 + 1), y);
        }
        if (iphistar > 3 * nphistar && iphistar < 4 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[0]->Fill(x - 3 * (10 + 1), y);
        }
        if (iphistar > 4 * nphistar && iphistar < 5 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[0]->Fill(x - 4 * (10 + 1), y);
        }
        if (iphistar > 5 * nphistar && iphistar < 6 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[1]->Fill(x - 5 * (10 + 1), y);
        }
        if (iphistar > 6 * nphistar && iphistar < 7 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[1]->Fill(x - 6 * (10 + 1), y);
        }
        if (iphistar > 7 * nphistar && iphistar < 8 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[1]->Fill(x - 7 * (10 + 1), y);
        }
        if (iphistar > 8 * nphistar && iphistar < 9 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[2]->Fill(x - 8 * (10 + 1), y);
        }
        if (iphistar > 9 * nphistar && iphistar < 10 * nphistar) {
            AllGraph->GetPoint(iphistar, x, y);
            AllHistos[3]->Fill(x - 9 * (10 + 1), y);
        }
    }
    double Normalizer = AllHistos[0]->Integral();
    if (Normalizer)AllHistos[0]->Scale(1 / Normalizer);

    Normalizer = AllHistos[1]->Integral();
    if (Normalizer) AllHistos[1]->Scale(1 / Normalizer);

    Normalizer = AllHistos[2]->Integral();
    if (Normalizer) AllHistos[2]->Scale(1 / Normalizer);

    Normalizer = AllHistos[3]->Integral();
    if (Normalizer) AllHistos[3]->Scale(1 / Normalizer);

   // Normalizer = AllHistos[4]->Integral();
   // if (Normalizer)AllHistos[4]->Scale(1 / Normalizer);
//
   // Normalizer = AllHistos[5]->Integral();
   // if (Normalizer) AllHistos[5]->Scale(1 / Normalizer);
//
   // Normalizer = AllHistos[6]->Integral();
   // if (Normalizer) AllHistos[6]->Scale(1 / Normalizer);
//
   // Normalizer = AllHistos[7]->Integral();
   // if (Normalizer)AllHistos[7]->Scale(1 / Normalizer);
//
   // Normalizer = AllHistos[8]->Integral();
   // if (Normalizer) AllHistos[8]->Scale(1 / Normalizer);
//
   // Normalizer = AllHistos[9]->Integral();
   // if (Normalizer)AllHistos[9]->Scale(1 / Normalizer);


    TCanvas* c1 = new TCanvas("c1", "", 800, 700);
    c1->SetFillColor(10);
    c1->SetFillColor(10);
    c1->cd();
    c1->SetLogy();
    c1->SetLogx();
    gStyle->SetOptStat("");
    gStyle->SetCanvasColor(0);
    gStyle->SetStatBorderSize(1);
    AllHistos[0]->SetLineColor(kBlack);
    AllHistos[1]->SetLineColor(kBlue);
    AllHistos[2]->SetLineColor(kRed);
    AllHistos[3]->SetLineColor(kGreen);
   // AllHistos[4]->SetLineColor(kOrange + 3);
   // AllHistos[6]->SetLineColor(kCyan);
   // AllHistos[7]->SetLineColor(kViolet + 6);
   // AllHistos[8]->SetLineColor(kMagenta + 1);
   // AllHistos[9]->SetLineColor(kYellow+2);

    AllHistos[0]->SetTitle("Seperated by rapidity");
    AllHistos[0]->Draw();
    AllHistos[1]->Draw("same");
    AllHistos[2]->Draw("same");
    AllHistos[3]->Draw("same");
    //AllHistos[4]->Draw("same");
    //AllHistos[5]->Draw("same");
    //AllHistos[6]->Draw("same");
    //AllHistos[7]->Draw("same");
    //AllHistos[8]->Draw("same");
    //AllHistos[9]->Draw("same");

    TLegend* leg = new TLegend(0.65, 0.65, 0.9, 0.9);
    leg-> SetFillColor(kWhite);
    leg-> SetShadowColor(kWhite);
    leg->AddEntry(AllHistos[0], Names[0].c_str());
    leg->AddEntry(AllHistos[1], Names[1].c_str());
    leg->AddEntry(AllHistos[2], Names[2].c_str());
    leg->AddEntry(AllHistos[3], Names[3].c_str());
   // leg->AddEntry(AllHistos[4], Names[4].c_str());
   // leg->AddEntry(AllHistos[5], Names[5].c_str());
   // leg->AddEntry(AllHistos[6], Names[6].c_str());
   // leg->AddEntry(AllHistos[7], Names[7].c_str());
   // leg->AddEntry(AllHistos[8], Names[8].c_str());
   // leg->AddEntry(AllHistos[9], Names[9].c_str());
    leg->Draw();

    c1->Print("YSeperated-0.50-1.00-1.50Dressed.png");
    delete c1;

}