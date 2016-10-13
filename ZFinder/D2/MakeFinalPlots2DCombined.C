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

const bool doNorm = false;
const int elec = 1;
const int doMG = 0;
const std::string Tag = "";
const std::string Type = "elec"; //elec, muon or combined
const std::string OutType = "png";
const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277, 5, 10};
size_t nphistar = (sizeof (phistarBins) / sizeof (phistarBins[0])) - 1;
const double yBins[] = {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};
size_t ny = (sizeof (yBins) / sizeof (yBins[0])) - 1;
size_t nbins = nphistar*ny;

void PrintValue(ofstream &outputfile, int n, int p, double value, double error) {
    if (p > 1 || p<-1) {
        if (n + p > 1)outputfile << std::setprecision(n + p) << "(" << value / double(pow(10, p)) << "&" << error / double(pow(10, p)) << ")$\\times 10^{" << p << "}$";
        else outputfile << std::setprecision(2) << "(" << value / double(pow(10, p)) << "&" << error / double(pow(10, p)) << ")$\\times 10^{" << p << "}$";
    } else {
        if (n > 1) outputfile << std::setprecision(n) << value << "&" << error;
        else outputfile << std::setprecision(2) << value << "&" << error;
    }
}

TGraphAsymmErrors* OffSetter(vector<TGraphAsymmErrors*> Orginal, bool Data = false) {
    size_t Newnphistar = nphistar - 2;
    TGraphAsymmErrors* NewPlot = new TGraphAsymmErrors(Newnphistar * ny);
    for (size_t ybin = 0; ybin < ny; ybin++) {
        for (size_t j = 0; j < Newnphistar - 2; j++) {
            double x, y;
            Orginal[ybin]->GetPoint(j, x, y);
            double Error = Orginal[ybin]->GetErrorYhigh(j);
            NewPlot->SetPoint(j + ybin * Newnphistar, x, y + .5 * (ny - ybin - 1));
            if (Data)NewPlot->SetPointError(j + ybin * Newnphistar, (phistarBins[j + 1] - phistarBins[j]) / 2., (phistarBins[j + 1] - phistarBins[j]) / 2., Error, Error);
            else NewPlot->SetPointError(j + ybin*Newnphistar, 0, 0, Error, Error);
        }
    }
    return NewPlot;
}

vector<TGraphAsymmErrors*> SplitGraph(TGraphAsymmErrors* graph, bool doXerrors = 0) {
    vector<TGraphAsymmErrors*> v;
    for (uint i = 0; i < ny; i++) {
        TGraphAsymmErrors* g = new TGraphAsymmErrors(nphistar);
        for (uint j = 0; j < nphistar; j++) {
            int bin = i * nphistar + j;
            double x, y;
            graph->GetPoint(bin, x, y);
            g->SetPoint(j, (phistarBins[j] + phistarBins[j + 1]) / 2., y);
            g->SetPointError(j, 0, 0, graph->GetErrorYlow(bin), graph->GetErrorYhigh(bin));
            if (doXerrors) g->SetPointError(j, (phistarBins[j + 1] - phistarBins[j]) / 2., (phistarBins[j + 1] - phistarBins[j]) / 2., graph->GetErrorYlow(bin), graph->GetErrorYhigh(bin));
            //cout << i << " " << j << " " << bin << " " << (phistarBins[j] + phistarBins[j + 1]) / 2. << " " << (phistarBins[j + 1] - phistarBins[j]) / 2. << " " << y << " " << graph->GetErrorYlow(bin) << " " << graph->GetErrorYhigh(bin) << endl;
        }
        v.push_back(g);
    }
    return v;
}

TGraphAsymmErrors* ResbosFromRaj(int FType = 0) {
    TGraphAsymmErrors* ResHolder = new TGraphAsymmErrors(nphistar);
    string FName = "";
    if (doNorm) {
        if (FType == 0)FName = "Resbos_0_1000_8TeV_Normalized.root";
        if (FType == 1)FName = "PowhegPythia8_NLO_AODSIM_8TeV_Normalized.root";
        if (FType == 2)FName = "AMCAT_NLO_AODSIM_8TeV_Normalized.root";
    } else {
        if (FType == 0)FName = "Resbos_0_1000_8TeV_Absolute.root";
        else if (FType == 1)FName = "PhiStar_PowhegPythia8_phistar_Status3_Absolute.root";
        else if (FType == 2) FName = "AMCAT_NLO_AODSIM_8TeV_Normalized.root";
    }
    TFile ResFile(FName.c_str());
    TH1D* Bin0 = (TH1D*) ResFile.Get("PhiStar_YBin_0");
    TH1D* Bin1 = (TH1D*) ResFile.Get("PhiStar_YBin_1");
    TH1D* Bin2 = (TH1D*) ResFile.Get("PhiStar_YBin_2");
    TH1D* Bin3 = (TH1D*) ResFile.Get("PhiStar_YBin_3");
    TH1D* Bin4 = (TH1D*) ResFile.Get("PhiStar_YBin_4");
    TH1D* Bin5 = (TH1D*) ResFile.Get("PhiStar_YBin_5");
    //TH1D* Allstuff = (TH1D*) ResFile.Get("PhiStar")
    if (Bin0 == 0)cout << "missing Bin0" << endl;
    if (Bin1 == 0)cout << "missing Bin1" << endl;
    if (Bin2 == 0)cout << "missing Bin2" << endl;
    if (Bin3 == 0)cout << "missing Bin3" << endl;
    if (Bin4 == 0)cout << "missing Bin4" << endl;
    if (Bin5 == 0)cout << "missing Bin5" << endl;

    if (!doNorm) {
        Bin0->Scale(1 / .4);
        Bin1->Scale(1 / .4);
        Bin2->Scale(1 / .4);
        Bin3->Scale(1 / .4);
        Bin4->Scale(1 / .4);
        Bin5->Scale(1 / .4);
    }

    for (uint i = 1; i <= nphistar; i++) {
        uint iphistar = i - 1;
        ResHolder->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2, Bin0->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 1 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 1 * phistarBins[nphistar], Bin1->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 2 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 2 * phistarBins[nphistar], Bin2->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 3 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 3 * phistarBins[nphistar], Bin3->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 4 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 4 * phistarBins[nphistar], Bin4->GetBinContent(i));
        ResHolder->SetPoint(iphistar + 5 * nphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2 + 5 * phistarBins[nphistar], Bin5->GetBinContent(i));

        ResHolder->SetPointEYhigh(iphistar, Bin0->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 1 * nphistar, Bin1->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 2 * nphistar, Bin2->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 3 * nphistar, Bin3->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 4 * nphistar, Bin4->GetBinError(i));
        ResHolder->SetPointEYhigh(iphistar + 5 * nphistar, Bin5->GetBinError(i));

        ResHolder->SetPointEYlow(iphistar, Bin0->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 1 * nphistar, Bin1->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 2 * nphistar, Bin2->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 3 * nphistar, Bin3->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 4 * nphistar, Bin4->GetBinError(i));
        ResHolder->SetPointEYlow(iphistar + 5 * nphistar, Bin5->GetBinError(i));

        ResHolder->SetPointEYhigh(iphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 1 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 2 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 3 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 4 * nphistar, .000001);
        ResHolder->SetPointEYhigh(iphistar + 5 * nphistar, .000001);

        ResHolder->SetPointEYlow(iphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 1 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 2 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 3 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 4 * nphistar, .000001);
        ResHolder->SetPointEYlow(iphistar + 5 * nphistar, .000001);
    }
    return ResHolder;
}

TGraphAsymmErrors* CreateDummy(TGraphAsymmErrors* graph) {
    vector<TGraphAsymmErrors*> g_dummyvec;
    // double x,y,errorl,errorh;
    double x, y;
    TGraphAsymmErrors* g_dummy = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        g_dummy->SetPoint(iphistar, x, y);
        g_dummy->SetPointError(iphistar, x - phistarBins[iphistar], x - phistarBins[iphistar], 0, 0);
    }
    g_dummy->GetXaxis()->SetRangeUser(0.001, 10);
    g_dummy->GetXaxis()->SetTitleOffset(1.05);
    g_dummy->GetXaxis()->SetTitle(0);
    g_dummy->GetXaxis()->SetTitleSize(0.05);
    g_dummy->GetXaxis()->SetLabelSize(0.05);
    g_dummy->GetYaxis()->SetTitleOffset(1.05);
    //  if (!doNorm)g_dummy->GetYaxis()->SetTitleOffset(1.05); 
    g_dummy->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi*dy ");
    if (!doNorm) g_dummy->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi*dy (pb)");
    g_dummy->GetYaxis()->SetTitleSize(0.05);
    g_dummy->GetYaxis()->SetLabelSize(0.05);
    g_dummy->SetFillColor(kWhite);
    g_dummy->SetTitle(0);
    return g_dummy;
}

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData) {
    double x, y, errorl, errorh, xmc, ymc, errorlmc, errorhmc;
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nbins);
    cout << "test 0" << endl;
    for (size_t ibin = 0; ibin < nbins; ibin++) {
        graph->GetPoint(ibin, x, y);
        if (!isData) {
            graphmc->GetPoint(ibin, xmc, ymc);
            errorlmc = graphmc->GetErrorYlow(ibin);
            errorhmc = graphmc->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, ymc / y);
            g_ratio->SetPointError(ibin, 0, 0, errorlmc / y, errorhmc / y);
        } else {
            errorl = graph->GetErrorYlow(ibin);
            errorh = graph->GetErrorYhigh(ibin);
            g_ratio->SetPoint(ibin, x, 1);
            g_ratio->SetPointError(ibin, 0, 0, errorl / y, errorh / y);
        }
    }
    return g_ratio;
}

vector<TGraphAsymmErrors*> CreateDummy(vector<TGraphAsymmErrors*> graphs) {
    vector<TGraphAsymmErrors*> v;
    for (uint i = 0; i < graphs.size(); i++) {
        TGraphAsymmErrors* d = CreateDummy(graphs[i]);
        v.push_back(d);
    }
    return v;
}

double IntValue(TGraphAsymmErrors* Graph) {
    double Output = 0;
    for (size_t i = 0; i < nphistar; i++) {
        double x, y;
        Graph->GetPoint(i, x, y);
        Output += y * (phistarBins[i + 1] - phistarBins[i]);
    }
    return Output;
}

double YErrorValue(TGraphAsymmErrors* Graph) {
    double YError = 0;
    for (size_t i = 0; i < nphistar; i++) {
        double BinSizeSquared = (phistarBins[i + 1] - phistarBins[i])*(phistarBins[i + 1] - phistarBins[i]);
        YError = sqrt(YError * YError + Graph->GetErrorYhigh(i) * Graph->GetErrorYhigh(i) * BinSizeSquared);
    }
    return YError;
}

void RatioValue(TGraphAsymmErrors* NewPlot, TGraphAsymmErrors* GraphDem, TGraphAsymmErrors* GraphNum) {

    for (int i = 0; i < (int) ny; i++) {
        double x1, y1, x2, y2, Ratio;
        double Error;
        GraphDem->GetPoint(i, x1, y1);
        GraphNum->GetPoint(i, x2, y2);
        Ratio = y2 / y1;
        NewPlot->SetPoint(i, x1, Ratio);
        Error = GraphNum->GetErrorY(i) / y1;

        NewPlot->SetPointError(i, 0, 0, Error, Error);
        if (GraphDem == GraphNum) {
            cout << "our error is " << Error << endl;
            NewPlot->SetPointError(i, .2, .2, Error, Error);
        }
    }
}

void OneDYPlot(vector<TGraphAsymmErrors*> g_data, vector<TGraphAsymmErrors*> g_mg, vector<TGraphAsymmErrors*> g_ph, vector<TGraphAsymmErrors*> g_re, vector<TGraphAsymmErrors*> g_ANlo, vector<TGraphAsymmErrors*> g_Pyth8) {
    TGraphAsymmErrors* g_data_Y = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_mg_Y = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_ph_Y = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_re_Y = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_ANlo_Y = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_Pyh8_Y = new TGraphAsymmErrors(ny);

    TGraphAsymmErrors* g_Data_Y_ratio = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_mg_Y_ratio = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_ph_Y_ratio = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_re_Y_ratio = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_ANlo_Y_ratio = new TGraphAsymmErrors(ny);
    TGraphAsymmErrors* g_Pyh8_Y_ratio = new TGraphAsymmErrors(ny);

    for (size_t i = 0; i < ny; i++) {
        cout << "our error for the first bin is " << g_data[i]->GetErrorYhigh(0) << endl;
    }

    double YError = 0;
    size_t YBin = 5;
    for (size_t i = 0; i < nphistar; i++) {

        double BinSizeSquared = (phistarBins[i + 1] - phistarBins[i])*(phistarBins[i + 1] - phistarBins[i]);
        double x, y;

        g_data[YBin]->GetPoint(i, x, y);
        cout << "Error For Phistar bin " << i << " is: " << sqrt(g_data[YBin]->GetErrorYhigh(i) * g_data[YBin]->GetErrorYhigh(i)) << endl;
        YError = sqrt(YError * YError + g_data[YBin]->GetErrorYhigh(i) * g_data[YBin]->GetErrorYhigh(i) * BinSizeSquared);
    }

    cout << "OOOKAy our int value is " << IntValue(g_data[YBin]) << endl;
    cout << "And our error is " << YError << endl;

    for (int i = 0; i < (int) ny; i++) {


        g_data_Y->SetPoint(i, .2 + .4 * i, IntValue(g_data[i]));
        g_mg_Y->SetPoint(i, .2 + .4 * i, IntValue(g_mg[i]));
        g_ph_Y->SetPoint(i, .2 + .4 * i, IntValue(g_ph[i]));
        g_re_Y->SetPoint(i, .2 + .4 * i, IntValue(g_re[i]));
        g_ANlo_Y->SetPoint(i, .2 + .4 * i, IntValue(g_ANlo[i]));
        g_Pyh8_Y->SetPoint(i, .2 + .4 * i, IntValue(g_Pyth8[i]));

        g_data_Y->SetPointError(i, 0, 0, YErrorValue(g_data[i]), YErrorValue(g_data[i]));
        g_mg_Y->SetPointError(i, 0, 0, YErrorValue(g_mg[i]), YErrorValue(g_mg[i]));
        g_ph_Y->SetPointError(i, 0, 0, YErrorValue(g_ph[i]), YErrorValue(g_ph[i]));
        g_re_Y->SetPointError(i, 0, 0, YErrorValue(g_re[i]), YErrorValue(g_re[i]));
        g_ANlo_Y->SetPointError(i, 0, 0, YErrorValue(g_ANlo[i]), YErrorValue(g_ANlo[i]));
        g_Pyh8_Y->SetPointError(i, 0, 0, YErrorValue(g_Pyth8[i]), YErrorValue(g_Pyth8[i]));

    }



    RatioValue(g_Data_Y_ratio, g_data_Y, g_data_Y);
    RatioValue(g_mg_Y_ratio, g_data_Y, g_mg_Y);
    RatioValue(g_ph_Y_ratio, g_data_Y, g_ph_Y);
    RatioValue(g_re_Y_ratio, g_data_Y, g_re_Y);
    RatioValue(g_ANlo_Y_ratio, g_data_Y, g_ANlo_Y);
    RatioValue(g_Pyh8_Y_ratio, g_data_Y, g_Pyh8_Y);


    TCanvas* FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot", 800, 900);
    FinalPhiTot->Divide(1, 2);
    FinalPhiTot->cd(1);
    gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
    gPad->SetBottomMargin(0.01);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    //gPad->SetLogy(1);
    g_data_Y->SetTitle("");
    g_data_Y->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/dy");
    g_data_Y->GetYaxis()->SetTitleSize(.065);
    g_data_Y->GetYaxis()->SetTitleOffset(1.);
    g_data_Y->GetYaxis()->CenterTitle();
    g_data_Y->GetYaxis()->SetLabelSize(0.05);
    g_data_Y->GetXaxis()->SetRangeUser(0.0, 2.4);
    g_data_Y->GetXaxis()->SetLabelSize(.0);
    g_data_Y->SetFillColor(kGray);
    if (doNorm) g_data_Y->GetYaxis()->SetRangeUser(.01, .8);
    else g_data_Y->GetYaxis()->SetRangeUser(1000, 400000.0);
    g_data_Y->Draw("A2");
    int Color1 = kBlue;
    int Color2 = kRed;
    g_mg_Y->SetMarkerColor(Color1);
    g_mg_Y->SetLineColor(Color1);
    g_mg_Y->SetMarkerSize(1);
    g_mg_Y->SetLineWidth(2);
    g_mg_Y->SetMarkerStyle(4);
    g_mg_Y->Draw("PEsame");
    g_ph_Y->SetMarkerColor(Color2);
    g_ph_Y->SetLineColor(Color2);
    g_ph_Y->SetMarkerSize(1);
    g_ph_Y->SetLineWidth(2);
    g_ph_Y->SetMarkerStyle(20);
    g_ph_Y->Draw("PEsame");
    int ColorError = kYellow;
    g_data_Y->SetFillColor(ColorError);
    g_data_Y->SetMarkerSize(1);
    g_data_Y->SetLineWidth(2);
    g_data_Y->SetMarkerStyle(20);
    g_data_Y->Draw("PEsame");
    g_data_Y->SetFillColor(ColorError);

    TLegend* leg = new TLegend(0.23, 0.76, 0.95, 0.94);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);

    leg->AddEntry(g_data_Y, "2012 data", "PEF");

    leg->AddEntry(g_mg_Y, "Z #rightarrow ee MadGraph+Pythia6 (Z2star)", "P");
    leg->AddEntry(g_ph_Y, "Z #rightarrow ee POWHEG+Pythia6 (Z2star)", "P");
    leg->Draw();

    TLatex mark;
    // mark.SetTextSize(0.05);
    mark.SetTextSize(0.035);
    mark.SetNDC(true);
    mark.DrawLatex(0.745, 0.95, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.19, 0.955, "CMS Preliminary");
    //if (Type == "elec") {
    mark.DrawLatex(0.19, 0.20, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
    mark.DrawLatex(0.19, 0.13, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
    mark.DrawLatex(0.19, 0.06, "60 GeV < M_{ll} < 120 GeV");
    //}

    FinalPhiTot->cd(2);
    gPad->SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, kWhite, 0, 0);
    gPad->SetBottomMargin(0.37);
    gPad->SetTopMargin(0.01);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);


    g_Data_Y_ratio->SetTitle("");
    g_Data_Y_ratio->GetYaxis()->SetTitle("MC/Data ");
    g_Data_Y_ratio->GetYaxis()->SetTitleOffset(0.4);
    g_Data_Y_ratio->GetYaxis()->SetTitleSize(.15);
    g_Data_Y_ratio->SetFillColor(kYellow);
    g_Data_Y_ratio->GetYaxis()->SetRangeUser(0.9, 1.1);



    g_Data_Y_ratio->GetYaxis()->SetLabelSize(0.04);
    g_Data_Y_ratio->GetXaxis()->SetLabelSize(0.13);
    g_Data_Y_ratio->GetXaxis()->SetLabelOffset(0.02);
    g_Data_Y_ratio->GetXaxis()->SetRangeUser(0, 2.4);
    g_Data_Y_ratio->GetXaxis()->SetNdivisions(510);
    g_Data_Y_ratio->GetXaxis()->SetTitle("y_{ll}");
    g_Data_Y_ratio->GetXaxis()->SetTitleSize(0.15);
    g_Data_Y_ratio->GetXaxis()->SetTitleOffset(1);
    g_Data_Y_ratio->GetXaxis()->CenterTitle();
    g_Data_Y_ratio->GetYaxis()->SetNdivisions(505);
    g_Data_Y_ratio->GetYaxis()->SetLabelSize(0.13);
    g_Data_Y_ratio->Draw("AE2");
    g_mg_Y_ratio->SetMarkerSize(1);
    g_mg_Y_ratio->SetLineWidth(2);
    g_mg_Y_ratio->SetMarkerStyle(4);
    g_mg_Y_ratio->SetMarkerColor(Color1);
    g_mg_Y_ratio->SetLineColor(Color1);
    g_mg_Y_ratio->Draw("PEsame");
    g_ph_Y_ratio->SetMarkerSize(1);
    g_ph_Y_ratio->SetLineWidth(2);
    g_ph_Y_ratio->SetMarkerStyle(20);
    g_ph_Y_ratio->SetMarkerColor(Color2);
    g_ph_Y_ratio->SetLineColor(Color2);
    g_ph_Y_ratio->Draw("PEsame");
    std::string plotname = "Plots/YFull";
    if (!doMG) plotname += "PH_";
    if (doNorm) plotname += "Norm_";
    else plotname += "Abs_";
    if (elec == 0)plotname += "Dressed.";
    if (elec == 1)plotname += "Born.";
    if (elec == 2)plotname += "Naked.";
    FinalPhiTot->SaveAs((plotname + "pdf").c_str());
    FinalPhiTot->SaveAs((plotname + "png").c_str());
    FinalPhiTot->SaveAs((plotname + "C").c_str());
    //RatioPlotStart
    TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio", 800, 900);
    FinalPhiRatio->cd();
    g_Data_Y_ratio->SetTitle("");
    g_Data_Y_ratio->GetYaxis()->SetTitle("MC/Data");
    g_Data_Y_ratio->GetYaxis()->CenterTitle();
    g_Data_Y_ratio->SetFillColor(kYellow);
    g_Data_Y_ratio->GetYaxis()->SetRangeUser(0.5, 1.25);
    g_Data_Y_ratio->GetYaxis()->SetTitleOffset(1.2);
    g_Data_Y_ratio->GetYaxis()->SetTitleSize(0.04);
    g_Data_Y_ratio->GetYaxis()->SetLabelSize(0.04);
    g_Data_Y_ratio->GetXaxis()->SetRangeUser(0.0, 2.4);
    g_Data_Y_ratio->GetXaxis()->SetTitleSize(0.043);
    g_Data_Y_ratio->GetXaxis()->SetLabelSize(0.04);
    g_Data_Y_ratio->GetXaxis()->SetLabelOffset(0.012);
    g_Data_Y_ratio->Draw("AE2");

    g_mg_Y_ratio->SetMarkerStyle(kFullTriangleUp);
    g_mg_Y_ratio->SetMarkerColor(kBlue - 7);
    g_mg_Y_ratio->SetLineColor(kBlue - 7);

    g_ph_Y_ratio->SetMarkerSize(1);
    g_ph_Y_ratio->SetLineWidth(2);
    g_ph_Y_ratio->SetMarkerStyle(kFullSquare);
    g_ph_Y_ratio->SetMarkerColor(kRed);
    g_ph_Y_ratio->SetLineColor(kRed);
    g_ph_Y_ratio->Draw("PEsame");

    g_mg_Y_ratio->Draw("PEsame");

    g_re_Y_ratio->SetMarkerColor(kGreen + 1);
    g_re_Y_ratio->SetLineColor(kGreen + 1);
    g_re_Y_ratio->SetMarkerSize(1);
    g_re_Y_ratio->SetLineWidth(2);
    g_re_Y_ratio->SetMarkerStyle(kStar);
    g_re_Y_ratio->Draw("PEsame");

    g_ANlo_Y_ratio->SetMarkerSize(1);
    g_ANlo_Y_ratio->SetLineWidth(2);
    g_ANlo_Y_ratio->SetMarkerStyle(kOpenCircle);
    g_ANlo_Y_ratio->SetMarkerColor(kCyan + 2);
    g_ANlo_Y_ratio->SetLineColor(kCyan + 2);
    g_ANlo_Y_ratio->Draw("PEsame");

    g_Pyh8_Y_ratio->SetMarkerSize(1);
    g_Pyh8_Y_ratio->SetLineWidth(2);
    g_Pyh8_Y_ratio->SetMarkerStyle(kOpenSquare);
    g_Pyh8_Y_ratio->SetMarkerColor(kRed);
    g_Pyh8_Y_ratio->SetLineColor(kRed);
    g_Pyh8_Y_ratio->Draw("PEsame");

    //r_data_WO->Draw("PEsame");

    mark.SetTextSize(0.03);
    mark.DrawLatex(0.7, 0.907, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.19, 0.907, "CMS Preliminary");

    if (Type == "elec") {
        mark.DrawLatex(0.15, 0.25, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
        mark.DrawLatex(0.15, 0.20, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ll} < 120 GeV");
    }


    TLegend* leg2 = new TLegend(0.13, 0.72, 0.7, 0.9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetLineWidth(1);
    leg2->SetNColumns(1);
    leg2->SetTextFont(22);

    leg2->AddEntry(g_Data_Y_ratio, "2012 data", "F");
    if (Type == "elec") {
        leg2->AddEntry(g_mg_Y_ratio, " MadGraph+Pythia6 (Z2*)", "P");
        leg2->AddEntry(g_ph_Y_ratio, "POWHEG+Pythia6 (Z2*)", "P");
        leg2->AddEntry(g_re_Y_ratio, "Resbos", "P");
        leg2->AddEntry(g_ANlo_Y_ratio, "AMC@nlo+Pythia8(CUETP8M1)", "P");
        leg2->AddEntry(g_Pyh8_Y_ratio, "POWHEG+Pythia8 (CT10)", "P");
    }

    leg2->Draw();

    plotname = "Plots/Ratio_YRatio";
    plotname += Tag;
    //if (isPlot2 && Type == "combined")plotname += "MuEl";
    //else if (isPlot2 && Type == "elec")plotname += "PHMG";
    plotname += Type;
    plotname += "_";
    if ((Type == "elec") && !doMG) plotname += "PH_";
    if ((Type == "elec") && doMG) plotname += "MG_";
    if (doNorm) plotname += "Norm_";
    else plotname += "Abs_";
    if (elec == 0)plotname += "Dressed.";
    if (elec == 1)plotname += "Born.";
    if (elec == 2)plotname += "Naked.";
    FinalPhiRatio->RedrawAxis();
    FinalPhiRatio->SaveAs((plotname + "pdf").c_str());
    FinalPhiRatio->SaveAs((plotname + "png").c_str());
    FinalPhiRatio->SaveAs((plotname + "C").c_str());



}

void PlotFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, TGraphAsymmErrors* g_ratio_phistar, TGraphAsymmErrors* g_ratio_mg_phistar, TGraphAsymmErrors* g_ratio_ph_phistar, bool isPlot2 = 0, TGraphAsymmErrors* g_re_final = 0, TGraphAsymmErrors* g_ratio_re_phistar = 0) {

    vector<TGraphAsymmErrors*> g_data = SplitGraph(g_data_final);
    vector<TGraphAsymmErrors*> g_mg = SplitGraph(g_mg_final);
    vector<TGraphAsymmErrors*> g_ph = SplitGraph(g_ph_final);
    vector<TGraphAsymmErrors*> r_data = SplitGraph(g_ratio_phistar, 1);
    vector<TGraphAsymmErrors*> r_mg = SplitGraph(g_ratio_mg_phistar);
    vector<TGraphAsymmErrors*> r_ph = SplitGraph(g_ratio_ph_phistar);
    vector<TGraphAsymmErrors*> g_re = SplitGraph(g_re_final);
    vector<TGraphAsymmErrors*> r_re = SplitGraph(g_ratio_re_phistar);
    //New Graphs
    TGraphAsymmErrors* g_ANlo_final = ResbosFromRaj(2);
    TGraphAsymmErrors* r_ANlo_Ratio_final = CreateRatio(g_data_final, g_ANlo_final, 0);
    TGraphAsymmErrors* g_PowPyth8_final = ResbosFromRaj(1);
    TGraphAsymmErrors* r_PowPyth8_Ratio_final = CreateRatio(g_data_final, g_PowPyth8_final, 0);

    vector<TGraphAsymmErrors*> g_ANlo = SplitGraph(g_ANlo_final);
    vector<TGraphAsymmErrors*> r_ANlo = SplitGraph(r_ANlo_Ratio_final);
    vector<TGraphAsymmErrors*> g_Pyth8 = SplitGraph(g_PowPyth8_final);
    vector<TGraphAsymmErrors*> r_Pyth8 = SplitGraph(r_PowPyth8_Ratio_final);

    vector<TGraphAsymmErrors*> g_dummy = CreateDummy(g_data);
    OneDYPlot(g_data, g_mg, g_ph, g_re, g_ANlo, g_Pyth8);
    for (uint i = 0; i < ny; i++) {
        //  for (uint i=0; i<5; i++){ 
        std::ostringstream strs;
        strs << i;
        //TO HERE
        std::string Canvasname = "EventEff_Bin" + strs.str();

        TCanvas* FinalPhiTot = new TCanvas(Canvasname.c_str(), Canvasname.c_str(), 800, 900);
        FinalPhiTot->Divide(1, 2);
        FinalPhiTot->cd(1);
        gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
        gPad->SetBottomMargin(0.01);
        gPad->SetTopMargin(0.06);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.06);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        g_dummy[i]->GetXaxis()->SetRangeUser(0.001, 10.0);
        if (doNorm) g_dummy[i]->GetYaxis()->SetRangeUser(0.00000005, 1000.0);
        else g_dummy[i]->GetYaxis()->SetRangeUser(0.00005, 1000000.0);
        g_dummy[i]->GetXaxis()->CenterTitle();
        g_dummy[i]->GetYaxis()->CenterTitle();
        g_dummy[i]->Draw("A2");
        g_mg[i]->SetMarkerColor(kBlue - 7);
        g_mg[i]->SetLineColor(kBlue - 7);
        g_mg[i]->SetMarkerSize(1);
        g_mg[i]->SetLineWidth(2);
        g_mg[i]->SetMarkerStyle(21);
        g_mg[i]->Draw("PEsame");
        g_ph[i]->SetMarkerColor(kRed);
        g_ph[i]->SetLineColor(kRed);
        g_ph[i]->SetMarkerSize(1);
        g_ph[i]->SetLineWidth(2);
        g_ph[i]->SetMarkerStyle(22);
        g_ph[i]->Draw("PEsame");
        if (!isPlot2 && elec == 1 && g_re_final) {
            g_re[i]->SetMarkerColor(kGreen + 1);
            g_re[i]->SetLineColor(kGreen + 1);
            g_re[i]->SetMarkerSize(1);
            g_re[i]->SetLineWidth(2);
            g_re[i]->SetMarkerStyle(23);
            g_re[i]->Draw("PEsame");
        }
        g_data[i]->SetFillColor(kYellow);
        g_data[i]->SetMarkerSize(1);
        g_data[i]->SetLineWidth(2);
        g_data[i]->SetMarkerStyle(20);
        g_data[i]->Draw("PEsame");
        g_data[i]->SetFillColor(kYellow);

        TLegend* leg;
        if (isPlot2) leg = new TLegend(0.15, 0.06, 0.80, 0.27);
        else leg = new TLegend(0.15, 0.06, 0.80, 0.31); //TLegend(0.45,0.73,0.94,0.93);//.19 0.06
        leg->SetFillStyle(0);
        leg->SetBorderSize(0);
        leg->SetLineWidth(1);
        leg->SetNColumns(1);
        leg->SetTextFont(42);

        if (!isPlot2) {
            leg->AddEntry(g_data[i], "2012 data", "PEF");
            if (Type == "elec") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow ee (MadGraph+Pythia6 Z2*)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow ee (Powheg+Pythia6 Z2*)", "P");
                //ToDo AMCAt decisions
                //if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ee (Resbos)", "P");
                if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ee (POWHEG+Pythia8)", "P");
            }
            if (Type == "muon") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow #mu#mu (MadGraph+Pythia6 Z2*)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow #mu#mu (Powheg+Pythia6 Z2*)", "P");
                if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow #mu#mu (Resbos)", "P");
                //if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow #mu#mu (POWHEG+Pythia8 CUETP8M1)", "P");
            }
            if (Type == "combined") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "P");
                if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ll (Resbos)", "P");
                //if (elec == 1 && g_re_final) leg->AddEntry(g_re[i], "#gamma*/Z #rightarrow ll (POWHEG+Pythia8 CUETP8M1)", "P");
            }
        } else {
            if (Type == "combined") {
                leg->AddEntry(g_mg[i], "#gamma*/Z #rightarrow ee (2012 data)", "P");
                leg->AddEntry(g_ph[i], "#gamma*/Z #rightarrow #mu#mu (2012 data)", "P");
                if (doMG)leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "PEF");
                else leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "PEF");
            }
            if (Type == "elec") {
                leg->AddEntry(g_mg[i], "data (unfolded with Powheg)", "P");
                leg->AddEntry(g_ph[i], "data (unfolded with MadGraph)", "P");
                if (doMG)leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (MadGraph+Pythia6 Z2*)", "PEF");
                else leg->AddEntry(g_data[i], "#gamma*/Z #rightarrow ll (Powheg+Pythia6 Z2*)", "PEF");
            }
        }
        leg->Draw();

        TLatex mark3;
        mark3.SetTextSize(0.043);
        mark3.SetTextFont(42);
        mark3.SetNDC(true);
        mark3.DrawLatex(0.71, 0.955, "19.7 fb^{-1} (8 TeV)");
        TLatex mark;
        mark.SetTextSize(0.043);
        mark.SetTextFont(42);
        mark.SetNDC(true);
        if (Type == "elec" && !isPlot2) {
            mark.DrawLatex(0.53, 0.88, "|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
            mark.DrawLatex(0.53, 0.81, "p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
            mark.DrawLatex(0.53, 0.74, "60 GeV < M_{ee} < 120 GeV");
        }
        if (Type == "combined" || isPlot2) {
            mark.DrawLatex(0.53, 0.88, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
            mark.DrawLatex(0.53, 0.81, "p_{T}^{l_{0}} > 30 GeV,  p_{T}^{l_{1}} > 20 GeV");
            mark.DrawLatex(0.53, 0.74, "60 GeV < M_{ll} < 120 GeV");
        }
        if (i == 0) mark.DrawLatex(0.2, 0.88, "0.0 < |y_{ee}| < 0.4");
        if (i == 1) mark.DrawLatex(0.2, 0.88, "0.4 < |y_{ee}| < 0.8");
        if (i == 2) mark.DrawLatex(0.2, 0.88, "0.8 < |y_{ee}| < 1.2");
        if (i == 3) mark.DrawLatex(0.2, 0.88, "1.2 < |y_{ee}| < 1.6");
        if (i == 4) mark.DrawLatex(0.2, 0.88, "1.6 < |y_{ee}| < 2.0");
        if (i == 5) mark.DrawLatex(0.2, 0.88, "2.0 < |y_{ee}| < 2.4");
        FinalPhiTot->cd(2);
        gPad->SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, kWhite, 0, 0);
        gPad->SetBottomMargin(0.37);
        gPad->SetTopMargin(0.01);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.06);
        gPad->SetLogx(1);

        r_data[i]->SetLineWidth(2);
        r_data[i]->GetXaxis()->SetRangeUser(0.001, 10);
        r_data[i]->GetXaxis()->SetTitle("#phi*");
        r_data[i]->GetXaxis()->SetTitleOffset(1.05);
        r_data[i]->GetXaxis()->SetTitleSize(0.12);
        r_data[i]->GetXaxis()->SetLabelSize(0.12);
        r_data[i]->GetYaxis()->SetTitle("MC/Data  ");
        //   r_data[i]->GetYaxis()->SetRangeUser(0.8,1.2);
        r_data[i]->GetYaxis()->SetTitleOffset(0.32);
        r_data[i]->GetYaxis()->SetTitleSize(0.12);
        r_data[i]->GetYaxis()->SetLabelSize(0.12);
        r_data[i]->SetTitle(0);
        if (isPlot2) r_data[i]->GetYaxis()->SetTitle("Data/MC  ");
        else r_data[i]->GetYaxis()->SetTitle("MC/Data   ");
        r_data[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
        r_data[i]->GetYaxis()->SetNdivisions(2, 5, 0);
        if (isPlot2) r_data[i]->GetYaxis()->SetRangeUser(0.3, 1.35);
        if (isPlot2 && (!doNorm)) r_data[i]->GetYaxis()->SetRangeUser(0.7, 1.3);
        r_data[i]->GetYaxis()->SetNdivisions(2, 5, 0);
        r_data[i]->GetYaxis()->SetTitleOffset(0.45);
        r_data[i]->SetFillColor(kYellow);
        r_data[i]->GetXaxis()->SetTitleSize(0.15);
        r_data[i]->GetXaxis()->CenterTitle();
        r_data[i]->GetYaxis()->CenterTitle();
        r_data[i]->Draw("APE2");
        r_mg[i]->SetMarkerSize(1);
        r_mg[i]->SetLineWidth(2);
        r_mg[i]->SetMarkerStyle(21);
        r_mg[i]->SetMarkerColor(kBlue - 7);
        r_mg[i]->SetLineColor(kBlue - 7);
        r_mg[i]->Draw("PEsame");
        r_ph[i]->SetMarkerSize(1);
        r_ph[i]->SetLineWidth(2);
        r_ph[i]->SetMarkerStyle(22);
        r_ph[i]->SetMarkerColor(kRed);
        r_ph[i]->SetLineColor(kRed);
        r_ph[i]->Draw("PEsame");
        r_re[i]->SetMarkerSize(1);
        r_re[i]->SetLineWidth(2);
        r_re[i]->SetMarkerStyle(23);
        r_re[i]->SetMarkerColor(kGreen + 1);
        r_re[i]->SetLineColor(kGreen + 1);
        r_re[i]->Draw("PEsame");

        std::string plotname = "Plots/ZShape_2D_";
        plotname += "Bin" + strs.str() + "_";
        plotname += Tag;
        if (isPlot2 && Type == "combined")plotname += "MuEl";
        else if (isPlot2 && Type == "elec")plotname += "PHMG";
        else plotname += Type;
        plotname += "_";
        if ((Type == "elec" || isPlot2) && !doMG) plotname += "PH_";
        if ((Type == "elec" || isPlot2) && doMG) plotname += "MG_";
        if (doNorm) plotname += "Norm_";
        else plotname += "Abs_";
        if (elec == 0)plotname += "Dressed.";
        if (elec == 1)plotname += "Born.";
        if (elec == 2)plotname += "Naked.";
        //FinalPhiTot->SaveAs((plotname + OutType).c_str());
        delete FinalPhiTot;
        //TO HERE ZACH


        TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio", 800, 900);
        FinalPhiRatio->cd();
        FinalPhiRatio->SetLogx();
        r_data[i]->GetYaxis()->SetTitle("MC/Data");
        r_data[i]->SetFillColor(kYellow);
        r_data[i]->GetYaxis()->SetRangeUser(0.5, 1.5);
        r_data[i]->GetYaxis()->SetTitleOffset(1.2);
        r_data[i]->GetYaxis()->SetTitleSize(0.04);
        r_data[i]->GetYaxis()->SetLabelSize(0.04);
        r_data[i]->GetXaxis()->SetTitleSize(0.04);
        r_data[i]->GetXaxis()->SetLabelSize(0.04);
        r_data[i]->GetXaxis()->SetLabelOffset(-0.01);
        r_data[i]->Draw("AE2");

        r_mg[i]->SetMarkerStyle(kFullTriangleUp);
        r_mg[i]->SetMarkerColor(kBlue - 7);
        r_mg[i]->SetLineColor(kBlue - 7);

        r_ph[i]->SetMarkerSize(1);
        r_ph[i]->SetLineWidth(2);
        r_ph[i]->SetMarkerStyle(kFullSquare);
        r_ph[i]->SetMarkerColor(kRed);
        r_ph[i]->SetLineColor(kRed);
        r_ph[i]->Draw("PEsame");

        r_mg[i]->Draw("PEsame");

        r_re[i]->SetMarkerColor(kGreen + 1);
        r_re[i]->SetLineColor(kGreen + 1);
        r_re[i]->SetMarkerSize(1);
        r_re[i]->SetLineWidth(2);
        r_re[i]->SetMarkerStyle(kStar);
        r_re[i]->Draw("PEsame");

        r_ANlo[i]->SetMarkerSize(1);
        r_ANlo[i]->SetLineWidth(2);
        r_ANlo[i]->SetMarkerStyle(kOpenCircle);
        r_ANlo[i]->SetMarkerColor(kCyan + 2);
        r_ANlo[i]->SetLineColor(kCyan + 2);
        r_ANlo[i]->Draw("PEsame");

        r_Pyth8[i]->SetMarkerSize(1);
        r_Pyth8[i]->SetLineWidth(2);
        r_Pyth8[i]->SetMarkerStyle(kOpenSquare);
        r_Pyth8[i]->SetMarkerColor(kRed);
        r_Pyth8[i]->SetLineColor(kRed);
        r_Pyth8[i]->Draw("PEsame");

        //r_data[i]->Draw("PEsame");

        mark.SetTextSize(0.03);
        mark.DrawLatex(0.7, 0.907, "19.7 fb^{-1} (8 TeV)");
        mark.DrawLatex(0.19, 0.907, "CMS Preliminary");
        if (i == 0) mark.DrawLatex(0.15, 0.3, "0.0 < |y_{ll}| < 0.4");
        if (i == 1) mark.DrawLatex(0.15, 0.3, "0.4 < |y_{ll}| < 0.8");
        if (i == 2) mark.DrawLatex(0.15, 0.3, "0.8 < |y_{ll}| < 1.2");
        if (i == 3) mark.DrawLatex(0.15, 0.3, "1.2 < |y_{ll}| < 1.6");
        if (i == 4) mark.DrawLatex(0.15, 0.3, "1.6 < |y_{ll}| < 2.0");
        if (i == 5) mark.DrawLatex(0.15, 0.35, "2.0 < |y_{ll}| < 2.4");
        if (Type == "elec" && !isPlot2) {
            mark.DrawLatex(0.15, 0.25, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
            mark.DrawLatex(0.15, 0.20, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
            if (i != 5)mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ll} < 120 GeV");
            else mark.DrawLatex(0.15, 0.30, "60 GeV < M_{ll} < 120 GeV");
        }
        if (Type == "muon" && !isPlot2) {
            mark.DrawLatex(0.15, 0.25, "|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
            mark.DrawLatex(0.15, 0.20, "p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
            mark.DrawLatex(0.15, 0.15, "60 GeV < M_{#mu#mu} < 120 GeV");
        }
        if (Type == "combined" || isPlot2) {
            mark.DrawLatex(0.15, 0.25, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
            mark.DrawLatex(0.15, 0.20, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
            mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ll} < 120 GeV");
        }

        TLegend* leg2 = new TLegend(0.13, 0.72, 0.7, 0.9);
        leg2->SetFillStyle(0);
        leg2->SetBorderSize(0);
        leg2->SetLineWidth(1);
        leg2->SetNColumns(1);
        leg2->SetTextFont(22);

        leg2->AddEntry(r_data[i], "2012 data", "F");
        if (Type == "elec") {
            leg2->AddEntry(r_mg[i], " MadGraph+Pythia6 (Z2*)", "P");
            leg2->AddEntry(r_ph[i], "POWHEG+Pythia6 (Z2*)", "P");
            leg2->AddEntry(r_re[i], "Resbos", "P");
            leg2->AddEntry(r_ANlo[i], "AMC@nlo+Pythia8(CUETP8M1)", "P");
            leg2->AddEntry(r_Pyth8[i], "POWHEG+Pythia8 (CT10)", "P");
        }
        if (Type == "muon") {
            leg2->AddEntry(r_mg[i], "Z #rightarrow #mu#mu MadGraph", "P");
            leg2->AddEntry(r_ph[i], "Z #rightarrow #mu#mu Powheg", "P");
        }
        if (Type == "combined") {
            leg2->AddEntry(r_mg[i], "Z #rightarrow ll MadGraph", "P");
            leg2->AddEntry(r_ph[i], "Z #rightarrow ll Powheg", "P");
        }
        leg2->Draw();
        plotname = "ZShape_Ratio";

        plotname += Tag;
        if (isPlot2)plotname += "MuEl";
        else plotname += Type;
        plotname = "Plots/Ratio_ZShape_2D_";
        plotname += "Bin" + strs.str() + "_";
        plotname += Tag;
        //if (isPlot2 && Type == "combined")plotname += "MuEl";
        //else if (isPlot2 && Type == "elec")plotname += "PHMG";
        plotname += Type;
        plotname += "_";
        if ((Type == "elec" || isPlot2) && !doMG) plotname += "PH_";
        if ((Type == "elec" || isPlot2) && doMG) plotname += "MG_";
        if (doNorm) plotname += "Norm_";
        else plotname += "Abs_";
        if (elec == 0)plotname += "Dressed.";
        if (elec == 1)plotname += "Born.";
        if (elec == 2)plotname += "Naked.";
        FinalPhiRatio->RedrawAxis();
        FinalPhiRatio->SaveAs((plotname + "pdf").c_str());
        FinalPhiRatio->SaveAs((plotname + "png").c_str());
        FinalPhiRatio->SaveAs((plotname + "C").c_str());
    }


    //TOTAL PLOT FROM HERE
    TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio", 800, 900);
    FinalPhiRatio->SetBottomMargin(0.1);

    TPad* Info = new TPad("p1", "p1", 0, .03, 1, .9);
    Info->Draw();
    Info->cd();
    Info->Divide(1, 6, 0, 0);


    for (uint i = 0; i < ny; i++) {
        Info->cd(i + 1);
        std::ostringstream strs;
        strs << i;
        std::string gPadName = "p" + strs.str();
        //gPad->SetPad(gPadName.c_str(), gPadName.c_str(), 0, i / 6, 1, (1 + i) / 6, kWhite, 0, 0);
        //gPad->SetLeftMargin(0.15);
        //gPad->SetRightMargin(0.06);
        if (i == 5)gPad->SetBottomMargin(0.2);
        gPad->SetLogx(1);

        if (i == 0)r_data[i]->SetTitle("");
        else r_data[i]->SetTitle("");

        if (i == 0)r_data[i]->GetYaxis()->SetTitle("");
        else if (i == 1)r_data[i]->GetYaxis()->SetTitle("");
        else if (i == 2)r_data[i]->GetYaxis()->SetTitle("MC/Data");
        else if (i == 3)r_data[i]->GetYaxis()->SetTitle("");
        else if (i == 4)r_data[i]->GetYaxis()->SetTitle("");
        else if (i == 5)r_data[i]->GetYaxis()->SetTitle("");

        r_data[i]->GetYaxis()->CenterTitle();
        r_data[i]->SetFillColor(kYellow);
        r_data[i]->GetYaxis()->SetNdivisions(503);
        r_data[i]->GetYaxis()->SetRangeUser(0.865, 1.135);
        if (i == 4)r_data[i]->GetYaxis()->SetRangeUser(0.77, 1.23);
        if (i == 5)r_data[i]->GetYaxis()->SetRangeUser(0.4, 1.6);
        r_data[i]->GetYaxis()->SetTitleOffset(0.2); //OFFSET
        r_data[i]->GetYaxis()->SetTitleSize(0.2); //Y TITLE SIZE
        r_data[i]->GetYaxis()->SetLabelSize(0.15);
        r_data[i]->GetXaxis()->SetTitleSize(0);
        r_data[i]->GetXaxis()->SetLabelSize(0);
        r_data[i]->GetXaxis()->SetTitle("");
        r_data[i]->GetXaxis()->CenterTitle();
        if (i == ny - 1) {
            r_data[i]->GetXaxis()->SetTitleOffset(0.42); //OFFSET
            r_data[i]->GetXaxis()->SetTitleSize(0.2); //X TITLE SIZE
            r_data[i]->GetXaxis()->SetLabelSize(0.15);
            r_data[i]->GetXaxis()->SetLabelOffset(-0.01);
            r_data[i]->GetXaxis()->SetTitle("#phi*");
            r_data[i]->GetXaxis()->CenterTitle();
        }
        r_data[i]->GetXaxis()->SetRangeUser(.015, 1.95);
        r_data[i]->Draw("AE2");

        r_mg[i]->SetMarkerStyle(kFullTriangleUp);
        r_mg[i]->SetMarkerColor(kBlue - 7);
        r_mg[i]->SetLineColor(kBlue - 7);

        r_ph[i]->SetMarkerSize(1);
        r_ph[i]->SetLineWidth(2);
        r_ph[i]->SetMarkerStyle(kFullSquare);
        r_ph[i]->SetMarkerColor(kRed);
        r_ph[i]->SetLineColor(kRed);
        r_ph[i]->Draw("PEsame");

        r_mg[i]->Draw("PEsame");

        r_re[i]->SetMarkerColor(kGreen + 1);
        r_re[i]->SetLineColor(kGreen + 1);
        r_re[i]->SetMarkerSize(1);
        r_re[i]->SetLineWidth(2);
        r_re[i]->SetMarkerStyle(kStar);
        r_re[i]->Draw("PEsame");

        r_ANlo[i]->SetMarkerSize(1);
        r_ANlo[i]->SetLineWidth(2);
        r_ANlo[i]->SetMarkerStyle(kOpenCircle);
        r_ANlo[i]->SetMarkerColor(kCyan + 2);
        r_ANlo[i]->SetLineColor(kCyan + 2);
        r_ANlo[i]->Draw("PEsame");

        r_Pyth8[i]->SetMarkerSize(1);
        r_Pyth8[i]->SetLineWidth(2);
        r_Pyth8[i]->SetMarkerStyle(kOpenSquare);
        r_Pyth8[i]->SetMarkerColor(kRed);
        r_Pyth8[i]->SetLineColor(kRed);
        r_Pyth8[i]->Draw("PEsame");
    }
    //r_data_WO->Draw("PEsame");
    FinalPhiRatio->cd(0);
    TLatex mark;
    mark.SetTextSize(0.02);
    mark.SetNDC(kTRUE);
    mark.DrawLatex(.75, .95, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.1, .95, "CMS Preliminary");
    TLatex mark2;
    mark2.SetTextSize(0.02);
    mark2.SetNDC(kTRUE);
    mark2.DrawLatex(.12, .84, "0.0 < |y_{ll}| < 0.4");
    mark2.DrawLatex(.12, .7, "0.4 < |y_{ll}| < 0.8");
    mark2.DrawLatex(.12, .55, "0.8 < |y_{ll}| < 1.2");
    mark2.DrawLatex(.12, .42, "1.2 < |y_{ll}| < 1.6");
    mark2.DrawLatex(.12, .28, "1.6 < |y_{ll}| < 2.0");
    mark2.DrawLatex(.12, .15, "2.0 < |y_{ll}| < 2.4");
    //mark.DrawLatex(-0.15, 0.25, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
    //mark.DrawLatex(0.15, 0.20, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
    //if (i != 5)mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ll} < 120 GeV");
    //else mark.DrawLatex(0.15, 0.30, "60 GeV < M_{ll} < 120 GeV");


    TLegend* leg2 = new TLegend(0.09, 0.888, 0.9, 0.94);
    leg2->SetNColumns(3);
    leg2->SetFillStyle(0);
    //leg2->SetBorderSize(1);
    leg2->SetLineWidth(1);
    leg2->SetTextFont(22);

    leg2->AddEntry(r_data[1], "2012 data", "F");
    leg2->AddEntry(r_mg[1], " MadGraph+Pythia6 (Z2*)", "P");
    leg2->AddEntry(r_ph[1], "POWHEG+Pythia6 (Z2*)", "P");
    leg2->AddEntry(r_re[1], "Resbos", "P");
    leg2->AddEntry(r_ANlo[1], "AMC@nlo+Pythia8(CUETP8M1)", "P");
    leg2->AddEntry(r_Pyth8[1], "POWHEG+Pythia8 (CT10)", "P");

    leg2->Draw();
    std::string plotname = "Plots/Ratio_ZShape_2D_ALL";
    plotname += Tag;
    //if (isPlot2 && Type == "combined")plotname += "MuEl";
    //else if (isPlot2 && Type == "elec")plotname += "PHMG";
    plotname += Type;
    plotname += "_";
    if ((Type == "elec" || isPlot2) && !doMG) plotname += "PH_";
    if (doNorm) plotname += "Norm_";
    else plotname += "Abs_";
    if (elec == 0)plotname += "Dressed.";
    if (elec == 1)plotname += "Born.";

    if (elec == 2)plotname += "Naked.";
    FinalPhiRatio->RedrawAxis();
    FinalPhiRatio->SaveAs((plotname + "pdf").c_str());
    FinalPhiRatio->SaveAs((plotname + "png").c_str());
    FinalPhiRatio->SaveAs((plotname + "C").c_str());


    //TO HEREs


}

TGraphAsymmErrors* CreateCopy(TGraphAsymmErrors* graph, double scale = 1.0) {
    double x, y, errorl, errorh;
    TGraphAsymmErrors* g_copy = new TGraphAsymmErrors(nbins);
    for (size_t ibin = 0; ibin < nbins; ibin++) {

        graph->GetPoint(ibin, x, y);
        errorl = graph->GetErrorYlow(ibin);
        errorh = graph->GetErrorYhigh(ibin);
        g_copy->SetPoint(ibin, x, y * scale);
        g_copy->SetPointError(ibin, 0, 0, errorl*scale, errorh * scale);
    }
    return g_copy;
}

void MakeFinalPlots() {

    TGraphAsymmErrors* g_data_final;
    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_ph_final;
    TGraphAsymmErrors* g_re_final = 0;
    TGraphAsymmErrors* g_ratio_re_phistar = 0;

    if (Type == "elec") {
        std::string textn = "Output/Data_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        if (doMG) textn += "MG_";
        else textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "BornTT.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile dg(textn.c_str());
        TGraphAsymmErrors * g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
        g_data_final = CreateCopy(g_data_temp);

        textn = "Output/MC_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "MG_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile mg(textn.c_str());
        TGraphAsymmErrors * g_mg_temp = (TGraphAsymmErrors*) mg.Get("Graph");
        g_mg_final = CreateCopy(g_mg_temp);

        textn = "Output/MC_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile fph(textn.c_str());
        TGraphAsymmErrors * g_ph_temp = (TGraphAsymmErrors*) fph.Get("Graph");
        g_ph_final = CreateCopy(g_ph_temp);
    }

    // if (Type=="combined"){
    // }
    // if (elec==1){
    // }

    g_re_final = ResbosFromRaj();
    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_ratio_phistar = CreateRatio(g_data_final, g_data_final, 1);
    TGraphAsymmErrors* g_ratio_mg_phistar = CreateRatio(g_data_final, g_mg_final, 0);
    TGraphAsymmErrors* g_ratio_ph_phistar = CreateRatio(g_data_final, g_ph_final, 0);
    if (g_re_final != 0)g_ratio_re_phistar = CreateRatio(g_data_final, g_re_final, 0);


    cout << "test 1" << endl;
    //PrintFinal(g_data_final, g_mg_final, g_ph_final, 0, g_re_final);
    cout << "test 2" << endl;
    for (size_t i = 0; i < g_data_final->GetN(); i++) {
        cout << "our Error is  for bin " << i << " is :" << g_data_final->GetErrorYhigh(i) << endl;
    }

    PlotFinal(g_data_final, g_mg_final, g_ph_final, g_ratio_phistar, g_ratio_mg_phistar, g_ratio_ph_phistar, 0, g_re_final, g_ratio_re_phistar);
}
