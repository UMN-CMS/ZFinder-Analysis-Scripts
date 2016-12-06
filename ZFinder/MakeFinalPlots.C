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

#include <TMatrixD.h>

#include <iomanip>

#include <iostream>
#include <fstream>
#include <TAxis.h>

const bool AlexPlots = false;
const bool doNorm = false;
const bool debug = true;
const int elec = 1;
const int doMG = 0;
const std::string Tag = "";
const std::string Type = "combined"; //elec, muon or combined
bool AverageData = false;

const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277};
size_t nphistar = (sizeof (phistarBins) / sizeof (phistarBins[0])) - 1;
const double PHDataUnfoldUncertainty[] = {.0082, .0093, .0099, .0103, .0089, .0096, .0091, .0084, .0088, .0077, .0077, .0099, .0081, .0078, .0077, .0075, .0076,
    .0078, .0076, .0076, .0081, .0076, .0076, .0076, .0074, .0075, .0081, .0096, .0122, .0166, .0188, .0234, .0299, .0364};

TGraphAsymmErrors* SplittGraph(bool Muons, bool RemoveLumi) {//Seperates the muon+Electron into 
    string FileName;
    if (doNorm)FileName = "/home/user1/lesko/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFile/Comb_ForBlue_Norm_1D_Born.root";
    else FileName = "/home/user1/lesko/work/HomeWork/Phistar/CombineElectWithMu/EPlusMuFile/Comb_ForBlue_Abs_1D_Born.root";
    TFile AllStuff(FileName.c_str());
    cout << "test 2" << endl;
    TGraphAsymmErrors* graph = (TGraphAsymmErrors*) AllStuff.Get("Nominal");
    TGraphAsymmErrors* Output = new TGraphAsymmErrors(nphistar);
    TMatrixD* ErrorMatrix = (TMatrixD*) AllStuff.Get("CovM_tot");
    TMatrixD* ErrorLumiMatrix = (TMatrixD*) AllStuff.Get("CovM_lumi");
    if (RemoveLumi)(*ErrorMatrix) = (*ErrorMatrix)-(*ErrorLumiMatrix);
    if (!graph)cout << "Couldn't find the graph  " << endl;
    if (!ErrorMatrix)cout << "Couldn't find the marrix  " << endl;
    cout << "Begin loop" << endl;

    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, y, yerror;
        size_t BigBin = iphistar + nphistar * (Muons);
        graph->GetPoint(iphistar + nphistar * (Muons), x, y);
        x = (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2;
        Output->SetPoint(iphistar, x, y);
        yerror = sqrt((*ErrorMatrix)(BigBin, BigBin));
        Output->SetPointError(iphistar, 0, 0, yerror, yerror);
    }
    cout << "End loop" << endl;
    return Output;
}

TGraphAsymmErrors * ResbosFromRaj(int FType = 0) {//Only plotting Ratio soooo why not?
    TGraphAsymmErrors* ResHolder = new TGraphAsymmErrors(nphistar);
    string FName = "";
    if (doNorm) {
        if (FType == 0)FName = "D2/PhiStar_Resbos_phistar_Status3_Normalized.root";
        if (FType == 1)FName = "D2/PhiStar_PowhegPythia8_NLO_phistar_Status3_Normalized.root";
        if (FType == 2)FName = "D2/PhiStar_AMCAT_NLO_phistar_Status3_Normalized.root";
    } else {
        if (FType == 0)FName = "D2/PhiStar_Resbos_phistar_Status3_Absolute.root";
        else if (FType == 1)FName = "D2/PhiStar_PowhegPythia8_NLO_phistar_Status3_Absolute.root";
        else if (FType == 2) FName = "D2/PhiStar_AMCAT_NLO_phistar_Status3_Absolute.root";
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

    if (!doNorm && FType != 2) {
        Bin0->Scale(1 / .4);
        Bin1->Scale(1 / .4);
        Bin2->Scale(1 / .4);
        Bin3->Scale(1 / .4);
        Bin4->Scale(1 / .4);
        Bin5->Scale(1 / .4);
    }


    TH1D* phiStarHist = (TH1D*) ResFile.Get("Phi_Star");
    if(FType == 1)cout<<"Our errors are "<<phiStarHist->GetBinError(1)<<"  and bin value is "<<phiStarHist->GetBinContent(1)<<endl;;
    if (doNorm && FType != 0&&false) {
        for (uint i = 1; i <= nphistar; i++) {
            phiStarHist->SetBinContent(i, phiStarHist->GetBinContent(i) / (phistarBins[i] - phistarBins[i - 1]));
        }
    }
if(FType == 1)cout<<"Second Our errors are "<<phiStarHist->GetBinError(1)<<"  and bin value is "<<phiStarHist->GetBinContent(1)<<endl;;
    //phiStarHist->Add(Bin1);
    //phiStarHist->Add(Bin2);
    //phiStarHist->Add(Bin3);
    //phiStarHist->Add(Bin4);
    //phiStarHist->Add(Bin5);

    double Normilizer = 0;

    if (false) {
        for (uint i = 1; i <= nphistar; i++) {
            Normilizer += phiStarHist->GetBinContent(i)*(phistarBins[i] - phistarBins[i - 1]);
        }
        phiStarHist->Scale(1 / Normilizer);
    }
    cout << "Test Bin 1 is " << phiStarHist->GetBinContent(1) << endl;
    //TH1D* Allstuff = (TH1D*) ResFile.Get("PhiStar")
    if (phiStarHist == 0)cout << "missing phiStarHist" << endl;



    for (uint i = 1; i <= nphistar; i++) {
        uint iphistar = i - 1;
        ResHolder->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2, phiStarHist->GetBinContent(i));

        ResHolder->SetPointEYhigh(iphistar, phiStarHist->GetBinError(i));

        ResHolder->SetPointEYlow(iphistar, phiStarHist->GetBinError(i));

        if (FType == 0) {
            ResHolder->SetPointEYhigh(iphistar, .000001);

            ResHolder->SetPointEYlow(iphistar, .000001);
        }
    }
    return ResHolder;
}

TGraphAsymmErrors* ConvertToTGraph(TH1D* h) {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        g->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2., h->GetBinContent(iphistar + 1));
        g->SetPointError(iphistar, 0, 0, h->GetBinError(iphistar + 1), h->GetBinError(iphistar + 1));
    }
    return g;
}

TGraphAsymmErrors* CreateRatio(TGraphAsymmErrors* graph, TGraphAsymmErrors* graphmc, bool isData) {
    double x, y, errorl, errorh, xmc, ymc, errorlmc, errorhmc;
    TGraphAsymmErrors* g_ratio = new TGraphAsymmErrors(nphistar);

    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        errorl = graph->GetErrorYlow(iphistar);
        errorh = graph->GetErrorYhigh(iphistar);
        if (!isData) {
            graphmc->GetPoint(iphistar, xmc, ymc);
            errorlmc = graphmc->GetErrorYlow(iphistar);
            errorhmc = graphmc->GetErrorYhigh(iphistar);
            g_ratio->SetPoint(iphistar, x, ymc / y);
            g_ratio->SetPointError(iphistar, 0, 0, errorlmc / y, errorhmc / y);
        } else {
            g_ratio->SetPoint(iphistar, x, 1);
            g_ratio->SetPointError(iphistar, x - phistarBins[iphistar], x - phistarBins[iphistar], errorl / y, errorh / y);
        }
    }
    if (isData) {
        g_ratio->SetLineWidth(2);
        g_ratio->GetXaxis()->SetRangeUser(0.001, 3.2);
        g_ratio->GetXaxis()->SetTitle("#phi*");
        g_ratio->GetXaxis()->SetTitleOffset(1.05);
        g_ratio->GetXaxis()->SetTitleSize(0.12);
        g_ratio->GetXaxis()->SetLabelSize(0.12);
        g_ratio->GetYaxis()->SetTitle("MC/Data  ");
        g_ratio->GetYaxis()->SetRangeUser(0.8, 1.12);
        g_ratio->GetYaxis()->SetTitleOffset(0.32);
        g_ratio->GetYaxis()->SetTitleSize(0.12);
        g_ratio->GetYaxis()->SetLabelSize(0.12);
        g_ratio->GetYaxis()->SetNdivisions(3, 5, 0);
        g_ratio->SetTitle(0);
    }
    return g_ratio;
}

void AlexsPlots(TGraphAsymmErrors* g_ph2_finalAlex, TGraphAsymmErrors* g_ph1_finalAlex) {



    const double PowErrorsPer[] = {0.63, 0.63, 0.65, 0.65, 0.65, 0.65, 0.60, 0.60, 0.62, 0.58, 0.55,
        0.68, 0.60, 0.59, 0.57, 0.57, 0.61, 0.60, 0.61, 0.60, 0.62, 0.64, 0.65, 0.63, 0.63, 0.65, 0.67, 0.81, 1.05, 1.47, 1.74, 2.15, 2.76, 3.26}; //powheg percentage error]
    const std::string textn = "/home/user1/lesko/work/Downloads/output.root";
    TFile tr2(textn.c_str());
    TH1D* h_PH1_temp = (TH1D*) tr2.Get("phi_1"); //Powheg stuff
    TH1D* h_PH2_temp = (TH1D*) tr2.Get("phi_2"); //Powheg stuff
    double Normph1 = h_PH1_temp->Integral();
    h_PH1_temp->Scale(1 / Normph1, "width");

    double Normph2 = h_PH2_temp->Integral();
    h_PH2_temp->Scale(1 / Normph2, "width");
    if (!h_PH1_temp) {
        cout << "couldn't get phi 1" << endl;
    }
    if (!h_PH2_temp) {
        cout << "couldn't get phi 2" << endl;
    }

    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        g_ph1_finalAlex->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2., h_PH1_temp->GetBinContent(iphistar + 1));
        g_ph1_finalAlex->SetPointError(iphistar, 0, 0, PowErrorsPer[iphistar] * h_PH1_temp->GetBinContent(iphistar + 1) / 100, PowErrorsPer[iphistar] * h_PH1_temp->GetBinContent(iphistar + 1) / 100);
    }
    //g_ratio_ph1_phistarAlex = CreateRatio(g_data_final, g_ph1_finalAlex, 0);

    //g_ph2_finalAlex = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        g_ph2_finalAlex->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2., h_PH2_temp->GetBinContent(iphistar + 1));
        g_ph2_finalAlex->SetPointError(iphistar, 0, 0, PowErrorsPer[iphistar] * h_PH2_temp->GetBinContent(iphistar + 1) / 100, PowErrorsPer[iphistar] * h_PH2_temp->GetBinContent(iphistar + 1) / 100);
    }
    // cout<<"making phi 2 ratio"<<endl;

    if (!g_ph2_finalAlex)cout << "g_ph2_finalAlex is inside null" << endl;
    //  TGraphAsymmErrors* g_un = new TGraphAsymmErrors(nphistar);


    //FinalPhiRatio->SaveAs((plotname + "C").c_str());
}

void PrintFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, bool dataonly = 0) {
    ofstream outputfile;
    std::string textname = "Table_Values_";
    textname += Tag;
    if (!dataonly) {
        textname += Type;
        textname += "_";
    }
    if (dataonly) textname += "MuEl_";
    if (Type == "elec" && !doMG) textname += "PH_";
    if (doNorm) textname += "Norm_";
    else textname += "Abs_";
    if (elec == 0)textname += "Dressed.txt";
    if (elec == 1)textname += "Born.txt";
    if (elec == 2)textname += "Naked.txt";
    outputfile.open(textname.c_str());
    std::string tableheader = "";
    if (!dataonly) {
        if (!doNorm) tableheader = "$\\phi^*$ range & \\multicolumn{2}{c}{Data \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{MadGraph \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{Powheg \\frac{d\\sigma^{fid}}{d\\phi*}(pb)}\\\\ \\hline";
        else tableheader = "$\\phi^*$ range & \\multicolumn{2}{c}{Data \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{MadGraph \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{Powheg \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}}\\\\ \\hline";
    } else {
        if (!doNorm) tableheader = "$\\phi^*$ range & \\multicolumn{2}{c}{Data: Z \\rightarrow ll \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{Data: Z \\rightarrow ee \\frac{d\\sigma^{fid}}{d\\phi*}(pb)} & \\multicolumn{2}{c}{Data: Z \\rightarrow \\mu\\mu \\frac{d\\sigma^{fid}}{d\\phi*}(pb)}\\\\ \\hline";
        else tableheader = "$\\phi^*$ range & \\multicolumn{2}{c}{Data: Z \\rightarrow ll \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{Data: Z \\rightarrow ee \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}} & \\multicolumn{2}{c}{Data: Z \\rightarrow \\mu\\mu \\frac{1}{\\sigma^{fid}}\\frac{d\\sigma^{fid}}{d\\phi*}}\\\\ \\hline";
    }
    outputfile << tableheader << "\n";
    std::cout << tableheader << endl;
    for (size_t i = 0; i < nphistar; i++) {
        double x, y, xmg, ymg, xph, yph;
        g_data_final->GetPoint(i, x, y);
        g_mg_final->GetPoint(i, xmg, ymg);
        g_ph_final->GetPoint(i, xph, yph);
        double temp_d = g_data_final->GetErrorYhigh(i);
        double temp_m = g_mg_final->GetErrorYhigh(i);
        double temp_p = g_ph_final->GetErrorYhigh(i);
        int n_d = 0;
        int n_m = 0;
        int n_p = 0;
        while (temp_d < 1) {
            temp_d = temp_d * 10.;
            n_d++;
        }
        while (temp_m < 1) {
            temp_m = temp_m * 10.;
            n_m++;
        }
        while (temp_p < 1) {
            temp_p = temp_p * 10.;
            n_p++;
        }
        outputfile << std::fixed << std::setprecision(3) << phistarBins[i] << "-" << phistarBins[i + 1] << " & " << std::setprecision(n_d) << y << " & " << g_data_final->GetErrorYhigh(i) << " & " << std::setprecision(n_m) << ymg << " & " << g_mg_final->GetErrorYhigh(i) << " & " << std::setprecision(n_p) << yph << " & " << g_ph_final->GetErrorYhigh(i) << " \\\\ \\hline" << "\n";
    }
    outputfile.close();
}

void PlotFinal(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final, TGraphAsymmErrors* g_dummy_phistar, TGraphAsymmErrors* g_ratio_phistar, TGraphAsymmErrors* g_ratio_mg_phistar, TGraphAsymmErrors* g_ratio_ph_phistar, bool isPlot2 = 0) {

    double WordSize = .075;
    TCanvas* FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot", 800, 900);
    FinalPhiTot->Divide(1, 2);
    FinalPhiTot->cd(1);
    gPad->SetPad("p1", "p1", 0, .5, 1, 1, kWhite, 0, 0);
    gPad->SetBottomMargin(0.01);
    gPad->SetTopMargin(0.09);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    TGraphAsymmErrors* g_ph2_finalAlex;
    TGraphAsymmErrors* g_ratio_ph2_phistarAlex;
    TGraphAsymmErrors* g_ph1_finalAlex;
    TGraphAsymmErrors* g_ratio_ph1_phistarAlex;


    std::string textn = "Comb_Hist_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "PH_";
    textn += "Born.root";
    TFile tr(textn.c_str());
    TH1D* h_data_temp = (TH1D*) tr.Get("h_Comb");
    TGraphAsymmErrors* AveragedCombined = ConvertToTGraph(h_data_temp);
    TGraphAsymmErrors* AveragedRatio = new TGraphAsymmErrors(nphistar);
    AverageData = true;
    cout << " Error ? " << g_data_final->GetErrorYhigh(1);
    AveragedRatio = CreateRatio(g_data_final, AveragedCombined, 0);



    //make  Set Plot stuff
    TGraphAsymmErrors* g_Resbos_phistar = ResbosFromRaj(0);
    TGraphAsymmErrors*g_PowPyth8_phistar = ResbosFromRaj(1);
    TGraphAsymmErrors*g_AMCatnlo_phistar = ResbosFromRaj(2);

    TGraphAsymmErrors* g_ratio_Resbos_phistar = CreateRatio(g_data_final, ResbosFromRaj(0), 0);
    TGraphAsymmErrors* g_ratio_PowPyth8_phistar = CreateRatio(g_data_final, ResbosFromRaj(1), 0);
    TGraphAsymmErrors* g_ratio_AMCatnlo_phistar = CreateRatio(g_data_final, ResbosFromRaj(2), 0);


    g_ratio_mg_phistar->SetMarkerSize(1);
    g_ratio_mg_phistar->SetLineWidth(2);
    g_ratio_mg_phistar->SetMarkerStyle(kFullTriangleUp);
    g_ratio_mg_phistar->SetMarkerColor(kBlue - 7);
    g_ratio_mg_phistar->SetLineColor(kBlue - 7);


    g_ratio_ph_phistar->SetMarkerSize(1);
    g_ratio_ph_phistar->SetLineWidth(2);
    g_ratio_ph_phistar->SetMarkerStyle(kFullSquare);
    g_ratio_ph_phistar->SetMarkerColor(kRed);
    g_ratio_ph_phistar->SetLineColor(kRed);

    g_ratio_PowPyth8_phistar->SetMarkerSize(1);
    g_ratio_PowPyth8_phistar->SetLineWidth(2);
    g_ratio_PowPyth8_phistar->SetMarkerStyle(kOpenSquare);
    g_ratio_PowPyth8_phistar->SetMarkerColor(kRed);
    g_ratio_PowPyth8_phistar->SetLineColor(kRed);

    g_ratio_Resbos_phistar->SetMarkerColor(kGreen + 1);
    g_ratio_Resbos_phistar->SetLineColor(kGreen + 1);
    g_ratio_Resbos_phistar->SetMarkerSize(1);
    g_ratio_Resbos_phistar->SetLineWidth(2);
    g_ratio_Resbos_phistar->SetMarkerStyle(kStar);


    g_ratio_AMCatnlo_phistar->SetMarkerSize(1);
    g_ratio_AMCatnlo_phistar->SetLineWidth(2);
    g_ratio_AMCatnlo_phistar->SetMarkerStyle(kOpenCircle);
    g_ratio_AMCatnlo_phistar->SetMarkerColor(kCyan + 2);
    g_ratio_AMCatnlo_phistar->SetLineColor(kCyan + 2);



    g_mg_final->SetMarkerSize(1);
    g_mg_final->SetLineWidth(2);
    g_mg_final->SetMarkerStyle(kFullTriangleUp);
    g_mg_final->SetMarkerColor(kBlue - 7);
    g_mg_final->SetLineColor(kBlue - 7);


    g_ph_final->SetMarkerSize(1);
    g_ph_final->SetLineWidth(2);
    g_ph_final->SetMarkerStyle(kFullSquare);
    g_ph_final->SetMarkerColor(kRed);
    g_ph_final->SetLineColor(kRed);

    g_data_final->SetMarkerSize(1);
    g_data_final->SetLineWidth(2);
    g_data_final->SetMarkerStyle(kFullCircle);
    g_data_final->SetMarkerColor(kBlack);
    g_data_final->SetLineColor(kBlack);

    g_PowPyth8_phistar->SetMarkerSize(1);
    g_PowPyth8_phistar->SetLineWidth(2);
    g_PowPyth8_phistar->SetMarkerStyle(kOpenSquare);
    g_PowPyth8_phistar->SetMarkerColor(kRed);
    g_PowPyth8_phistar->SetLineColor(kRed);

    g_Resbos_phistar->SetMarkerColor(kGreen + 1);
    g_Resbos_phistar->SetLineColor(kGreen + 1);
    g_Resbos_phistar->SetMarkerSize(1);
    g_Resbos_phistar->SetLineWidth(2);
    g_Resbos_phistar->SetMarkerStyle(kStar);


    g_AMCatnlo_phistar->SetMarkerSize(1);
    g_AMCatnlo_phistar->SetLineWidth(2);
    g_AMCatnlo_phistar->SetMarkerStyle(kOpenCircle);
    g_AMCatnlo_phistar->SetMarkerColor(kCyan + 2);
    g_AMCatnlo_phistar->SetLineColor(kCyan + 2);

    //End


    if (AlexPlots && doNorm) {
        
        g_ph2_finalAlex = new TGraphAsymmErrors(nphistar);
        g_ratio_ph2_phistarAlex = new TGraphAsymmErrors(nphistar);
        g_ph1_finalAlex = new TGraphAsymmErrors(nphistar);
        g_ratio_ph1_phistarAlex = new TGraphAsymmErrors(nphistar);
        AlexsPlots(g_ph2_finalAlex, g_ph1_finalAlex);
        g_ratio_ph2_phistarAlex = CreateRatio(g_data_final, g_ph2_finalAlex, 0);

        g_ratio_ph1_phistarAlex = CreateRatio(g_data_final, g_ph1_finalAlex, 0);
        if (!g_ph2_finalAlex)cout << "g_ph2_finalAlex after null" << endl;
    }
    g_dummy_phistar->GetXaxis()->SetRangeUser(0.001, 3.0);
    if (doNorm) g_dummy_phistar->GetYaxis()->SetRangeUser(0.0012, 100.0);
    else g_dummy_phistar->GetYaxis()->SetRangeUser(0.8, 40000.0);
    g_dummy_phistar->GetXaxis()->CenterTitle();
    g_dummy_phistar->GetXaxis()->SetTitleOffset(.8);
    g_dummy_phistar->GetYaxis()->CenterTitle();
    g_dummy_phistar->GetYaxis()->SetTitleSize(1.2*WordSize);
    g_dummy_phistar->GetYaxis()->SetLabelSize(WordSize);
    g_dummy_phistar->GetYaxis()->SetTitleOffset(0.7);
    g_dummy_phistar->Draw("A2");
    int Color1 = kBlue;
    int Color2 = kRed;
    int Color3 = kGreen + 2;
    int Color4 = kOrange + 7;
    g_mg_final->Draw("PEsame");
    g_ph_final->Draw("PEsame");
    g_AMCatnlo_phistar->Draw("PESame");
    g_Resbos_phistar->Draw("PESame");
    g_PowPyth8_phistar->Draw("PESame");
    if (AlexPlots && doNorm) {
        if (debug) {
            cout << "test 2" << endl;
        }
        g_ph2_finalAlex->SetMarkerColor(Color3);
        if (debug) {
            cout << "test 3.1" << endl;
        }
        g_ph2_finalAlex->SetLineColor(Color3);
        g_ph2_finalAlex->SetMarkerSize(1);
        g_ph2_finalAlex->SetLineWidth(2);
        g_ph2_finalAlex->SetMarkerStyle(21);
        g_ph2_finalAlex->Draw("PEsame");


        g_ph1_finalAlex->SetMarkerColor(Color4);
        g_ph1_finalAlex->SetLineColor(Color4);
        g_ph1_finalAlex->SetMarkerSize(1);
        g_ph1_finalAlex->SetLineWidth(2);
        g_ph1_finalAlex->SetMarkerStyle(22);
        g_ph1_finalAlex->Draw("PEsame");
        //this will make the table of values and uncertainties. 
        ofstream VAndU;
        VAndU.open("ValuesAndUnvertaities.txt");
        VAndU << "phi* range, Data, MadGraph 2ZStar, POWHEG 2ZStar, POWHEG Tunepp 14, POWHEG Tunepp5 \n";
        for (size_t i = 0; i < nphistar; i++) {
            VAndU << phistarBins[i] << "-" << phistarBins[i + 1] << ", ";
            double x, y;
            double yerrhigh, yerrlow, yerrAve;
            g_data_final->GetPoint(i, x, y);
            VAndU << y << "+/-";
            if (doMG == 0)VAndU << y * PHDataUnfoldUncertainty[i] << ", ";
            g_mg_final->GetPoint(i, x, y);
            yerrhigh = g_mg_final->GetErrorYhigh(i);
            yerrlow = g_mg_final->GetErrorYlow(i);
            yerrAve = (yerrhigh + yerrlow) / 2;
            VAndU << y << "+/-" << yerrAve << ", ";
            g_ph_final->GetPoint(i, x, y);
            yerrhigh = g_ph_final->GetErrorYhigh(i);
            yerrlow = g_ph_final->GetErrorYlow(i);
            yerrAve = (yerrhigh + yerrlow) / 2;
            VAndU << y << "+/-" << yerrAve << ", ";

            g_ph1_finalAlex->GetPoint(i, x, y);
            yerrhigh = g_ph1_finalAlex->GetErrorYhigh(i);
            yerrlow = g_ph1_finalAlex->GetErrorYlow(i);
            yerrAve = (yerrhigh + yerrlow) / 2;
            VAndU << y << "+/-" << yerrAve << ", ";

            g_ph2_finalAlex->GetPoint(i, x, y);
            yerrhigh = g_ph2_finalAlex->GetErrorYhigh(i);
            yerrlow = g_ph2_finalAlex->GetErrorYlow(i);
            yerrAve = (yerrhigh + yerrlow) / 2;
            VAndU << y << "+/-" << yerrAve << "\n ";

        }
        VAndU.close();
    }
    int ColorError = kYellow;
    g_data_final->Draw("PEsame");
    g_data_final->SetFillColor(ColorError);

    double legendOffset = -.4;
    TLegend* leg = new TLegend(0.16, 0.45 + legendOffset, 0.4, 0.88 + legendOffset);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);
    leg->SetTextSize(WordSize);

    if (!isPlot2) {
        leg->AddEntry(g_data_final, "Combined Data", "P");
        if (Type == "elec") {
            //leg->AddEntry(g_mg_final, "Z #rightarrow ee MadGraph+PYTHIA6 (Z2*)", "P");
            leg->AddEntry(g_ph_final, "Z #rightarrow ee POWHEG+PYTHIA6 (Z2*)", "P");
            if (AlexPlots && doNorm) {
                leg->AddEntry(g_ph2_finalAlex, "Z #rightarrow ee POWHEG+PYTHIA8 (Tunepp 5)", "P");
                leg->AddEntry(g_ph1_finalAlex, "Z #rightarrow ee POWHEG+PYTHIA8 (Tunepp 14)", "P");
            }
        }
        if (Type == "muon") {
            //leg->AddEntry(g_mg_final, "Z #rightarrow #mu#mu MadGraph", "P");
            leg->AddEntry(g_ph_final, "Z #rightarrow #mu#mu Powheg", "P");
        }
        if (Type == "combined") {
            leg->AddEntry(g_ratio_mg_phistar, "MadGraph+PYTHIA6 (Z2*)", "P");
            leg->AddEntry(g_ratio_ph_phistar, "POWHEG+PYTHIA6 (Z2*)", "P");
            leg->AddEntry(g_ratio_PowPyth8_phistar, "POWHEG+PYTHIA8 (CT10)", "P");
            leg->AddEntry(g_ratio_Resbos_phistar, "ResBos", "P");
            leg->AddEntry(g_ratio_AMCatnlo_phistar, "aMC@NLO+PYTHIA8(CUETP8M1)", "P");
        }
    } else {
        //leg->AddEntry(g_mg_final, "data: Z #rightarrow ee", "P");
        leg->AddEntry(g_ph_final, "Data: Z #rightarrow #mu#mu", "P");
        leg->AddEntry(g_data_final, "Blue Combined: Z #rightarrow ll", "PEF");
    }
    leg->Draw();

    TLatex mark;
    // mark.SetTextSize(0.05);
    mark.SetTextSize(WordSize);
    mark.SetNDC(true);
    mark.DrawLatex(0.66, 0.93, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.15, 0.93, "CMS");
    if (Type == "elec" && !isPlot2) {
        mark.DrawLatex(0.19, 0.24, "|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
        mark.DrawLatex(0.19, 0.19, "p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
        mark.DrawLatex(0.19, 0.14, "60 GeV < M_{ee} < 120 GeV");
    }
    if (Type == "muon" && !isPlot2) {
        mark.DrawLatex(0.19, 0.28, "|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
        mark.DrawLatex(0.19, 0.22, "p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
        mark.DrawLatex(0.19, 0.14, "60 GeV < M_{#mu#mu} < 120 GeV");
    }
    //if (Type == "combined" || isPlot2) {
    //    double Yoffset = .6;
    //    mark.DrawLatex(0.58, 0.28 + Yoffset, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
    //    mark.DrawLatex(0.58, 0.2 + Yoffset, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
    //    mark.DrawLatex(0.58, 0.12 + Yoffset, "60 GeV < M_{ll} < 120 GeV");
    //}
    if (true) {
        FinalPhiTot->cd(2);
        gPad->SetPad("p2", "p2", 0, 0, 1, .5, kWhite, 0, 0);
        gPad->SetBottomMargin(0.2);
        gPad->SetTopMargin(0.01);
        gPad->SetLeftMargin(0.15);
        gPad->SetRightMargin(0.06);
        gPad->SetLogx(1);

        if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("Data/MC");
        else g_ratio_phistar->GetYaxis()->SetTitle("Theory/Data");
        g_ratio_phistar->GetYaxis()->CenterTitle();
        g_ratio_phistar->GetXaxis()->SetRangeUser(0.001, 3.0);
        g_ratio_phistar->GetYaxis()->SetRangeUser(0.76, 1.24);
        if (isPlot2) g_ratio_phistar->GetYaxis()->SetRangeUser(0.88, 1.12);
        if (isPlot2 && (!doNorm)) g_ratio_phistar->GetYaxis()->SetRangeUser(0.93, 1.07);
        g_ratio_phistar->GetYaxis()->SetTitleOffset(.7);
        g_ratio_phistar->SetFillColor(ColorError);
        g_ratio_phistar->GetYaxis()->SetRangeUser(.78, 1.22);
        g_ratio_phistar->GetXaxis()->CenterTitle();
        g_ratio_phistar->GetXaxis()->SetTitleOffset(.65);
        g_ratio_phistar->GetXaxis()->SetTitleSize(1.2*WordSize); 
        g_ratio_phistar->GetYaxis()->SetTitleSize(1.2*WordSize);
        g_ratio_phistar->GetXaxis()->SetLabelSize(WordSize);
        g_ratio_phistar->GetYaxis()->SetLabelSize(WordSize);
        g_ratio_phistar->Draw("AE2");
        g_ratio_mg_phistar->Draw("PEsame");
        g_ratio_ph_phistar->Draw("PEsame");
        g_ratio_PowPyth8_phistar->Draw("PEsame");
        g_ratio_Resbos_phistar->Draw("PEsame");
        g_ratio_AMCatnlo_phistar->Draw("PEsame");
        TLine Test(0, 1, 3.12, 1);
        Test.Draw();
    }
    std::string plotname = "ZShape_";
    plotname += Tag;
    if (isPlot2)plotname += "MuEl";
    else plotname += Type;
    plotname += "_";
    if (Type == "elec" && !doMG) plotname += "PH_";
    if (doNorm) plotname += "Norm_";
    else plotname += "Abs_";
    if (elec == 0)plotname += "Dressed.";
    if (elec == 1)plotname += "Born.";
    if (elec == 2)plotname += "Naked.";
    FinalPhiTot->SaveAs((plotname + "pdf").c_str());
    FinalPhiTot->SaveAs((plotname + "png").c_str());
    FinalPhiTot->SaveAs((plotname + "C").c_str());
    
    if(doNorm)FinalPhiTot->SaveAs("/home/user1/lesko/work/Phistar/papers/SMP-15-002/trunk/NormAll.pdf");
    else FinalPhiTot->SaveAs("/home/user1/lesko/work/Phistar/papers/SMP-15-002/trunk/AbsAll.pdf");

    return;
    TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio", 800, 900);
    FinalPhiRatio->cd();
    FinalPhiRatio->SetRightMargin(.01);
    FinalPhiRatio->SetLogx();
    if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("MC/Data");
    else g_ratio_phistar->GetYaxis()->SetTitle("MC/Data");
    g_ratio_phistar->SetFillColor(kYellow);
    g_ratio_phistar->GetYaxis()->CenterTitle();
    g_ratio_phistar->GetXaxis()->CenterTitle();
    g_ratio_phistar->SetFillColor(ColorError);
    g_ratio_phistar->GetYaxis()->SetRangeUser(0.7, 1.3);

    g_ratio_phistar->GetYaxis()->SetTitleOffset(1.2);
    g_ratio_phistar->GetYaxis()->SetTitleSize(0.04);
    g_ratio_phistar->GetYaxis()->SetLabelSize(0.04);
    g_ratio_phistar->GetXaxis()->SetTitleOffset(.65);
    g_ratio_phistar->GetXaxis()->SetTitleSize(0.05);
    g_ratio_phistar->GetXaxis()->SetLabelSize(0.04);
    g_ratio_phistar->GetYaxis()->SetLabelSize(0.04);
    g_ratio_phistar->GetXaxis()->SetLabelOffset(-0.01);
    g_ratio_phistar->GetYaxis()->SetNdivisions(510);
    g_ratio_phistar->Draw("AE2");


    g_ratio_mg_phistar->SetMarkerSize(1);
    g_ratio_mg_phistar->SetLineWidth(2);
    g_ratio_mg_phistar->SetMarkerStyle(kOpenTriangleUp);
    g_ratio_mg_phistar->SetMarkerColor(kBlue - 7);
    g_ratio_mg_phistar->SetLineColor(kBlue - 7);
    g_ratio_mg_phistar->Draw("PEsame");



    g_ratio_ph_phistar->Draw("PEsame");


    g_ratio_PowPyth8_phistar->Draw("PEsame");
    g_ratio_Resbos_phistar->Draw("PEsame");
    g_ratio_AMCatnlo_phistar->Draw("PEsame");
    if (AlexPlots && doNorm) {
        g_ratio_ph2_phistarAlex->Draw("PEsame");
        g_ratio_ph1_phistarAlex->Draw("PEsame");
    }
    mark.SetTextSize(0.04);
    mark.DrawLatex(0.71, 0.907, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.101, 0.907, "CMS");
    if (Type == "elec" && !isPlot2) {
        mark.DrawLatex(0.15, 0.25, "|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
        mark.DrawLatex(0.15, 0.20, "p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ee} < 120 GeV");
    }
    if (Type == "muon" && !isPlot2) {
        mark.DrawLatex(0.15, 0.25, "|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
        mark.DrawLatex(0.15, 0.20, "p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{#mu#mu} < 120 GeV");
    }
    if (Type == "combined" || isPlot2) {
        mark.DrawLatex(0.15, 0.25, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
        mark.DrawLatex(0.15, 0.2, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ll} < 120 GeV");

        //AveragedRatio->SetMarkerSize(1);
        //AveragedRatio->SetLineWidth(2);
        //AveragedRatio->SetMarkerStyle(kOpenCircle);
        //AveragedRatio->SetMarkerColor(kCyan + 2);
        //AveragedRatio->SetLineColor(kCyan + 2);
        //AveragedRatio->Draw("PEsame");

    }
    TLegend* leg2 = new TLegend(0.1, 0.66, 0.65, 0.9);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetLineWidth(1);
    leg2->SetNColumns(1);
    leg2->SetTextFont(42);
    leg2->SetTextSize(.04);
    if (!isPlot2) {
        //g_ratio_phistar->SetMarkerColor(kYellow);
        //leg2->AddEntry(g_ratio_phistar, "Data", "EF");
        if (Type == "elec") {
            leg2->AddEntry(g_ratio_mg_phistar, "MadGraph+PYTHIA6 (Z2*)", "P");
            leg2->AddEntry(g_ratio_ph_phistar, "POWHEG+PYTHIA6 (Z2*)", "P");
            leg2->AddEntry(g_ratio_PowPyth8_phistar, "POWHEG+PYTHIA8 (CT10)", "P");
            leg2->AddEntry(g_ratio_Resbos_phistar, "Resbos", "P");
            leg2->AddEntry(g_ratio_AMCatnlo_phistar, "aMC@NLO+PYTHIA8(CUETP8M1)", "P");
            if (AlexPlots && doNorm) {
                leg2->AddEntry(g_ratio_ph2_phistarAlex, "POWHEG+PYTHIA8  (Tunepp 5)", "P");
                leg2->AddEntry(g_ratio_ph1_phistarAlex, "POWHEG+PYTHIA8 (Tunepp 14)", "P");
            }
        }
        if (Type == "muon") {
            leg2->AddEntry(g_mg_final, "Z #rightarrow #mu#mu MadGraph", "P");
            leg2->AddEntry(g_ph_final, "Z #rightarrow #mu#mu Powheg", "P");
        }
        if (Type == "combined") {
            leg2->AddEntry(g_ratio_mg_phistar, "MadGraph+PYTHIA6 (Z2*)", "P");
            leg2->AddEntry(g_ratio_ph_phistar, "POWHEG+PYTHIA6 (Z2*)", "P");
            leg2->AddEntry(g_ratio_PowPyth8_phistar, "POWHEG+PYTHIA8 (CT10)", "P");
            leg2->AddEntry(g_ratio_Resbos_phistar, "Resbos", "P");
            leg2->AddEntry(g_ratio_AMCatnlo_phistar, "aMC@NLO+PYTHIA8(CUETP8M1)", "P");

            //leg2->AddEntry(AveragedRatio, "Z #rightarrow ll Averaged Combined", "P");
        }
    } else {
        leg2->AddEntry(g_mg_final, "Data: Z #rightarrow ee", "P");
        leg2->AddEntry(g_ph_final, "Data: Z #rightarrow #mu#mu", "P");
        leg2->AddEntry(g_data_final, "Blue Combined: Z #rightarrow ll", "PEF");
    }
    leg2->Draw();
    plotname = "ZShape_Ratio_Blue";
    if (AlexPlots && doNorm) {
        plotname += "Alex";
    }
    plotname += Tag;
    if (isPlot2)plotname += "MuEl";
    else plotname += Type;
    plotname += "_";
    if (Type == "elec" && !doMG) plotname += "PH_";
    if (doNorm) plotname += "Norm_";
    else plotname += "Abs_";

    if (elec == 0)plotname += "Dressed.";
    if (elec == 1)plotname += "Born.";
    if (elec == 2)plotname += "Naked.";
    if (AlexPlots && doNorm) {

    }
    FinalPhiRatio->SaveAs((plotname + "pdf").c_str());
    FinalPhiRatio->SaveAs((plotname + "png").c_str());
    FinalPhiRatio->SaveAs((plotname + "C").c_str());
}

TGraphAsymmErrors* CreateDummy(TGraphAsymmErrors* graph) {
    vector<TGraphAsymmErrors*> g_dummyvec;
    // double x,y,errorl,errorh;
    double x, y;
    TGraphAsymmErrors* g_dummy = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        // errorl=graph->GetErrorYlow(iphistar);
        // errorh=graph->GetErrorYhigh(iphistar);
        g_dummy->SetPoint(iphistar, x, y);
        g_dummy->SetPointError(iphistar, x - phistarBins[iphistar], x - phistarBins[iphistar], 0, 0);
    }
    g_dummy->GetXaxis()->SetRangeUser(0.001, 3.2);
    g_dummy->GetXaxis()->SetTitleOffset(1.05);
    g_dummy->GetXaxis()->SetTitle("#phi*");
    g_dummy->GetXaxis()->SetTitleSize(0.05);
    g_dummy->GetXaxis()->SetLabelSize(0.05);
    g_dummy->GetYaxis()->SetTitleOffset(1.05);
    //  if (!doNorm)g_dummy->GetYaxis()->SetTitleOffset(1.05); 
    g_dummy->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/d#phi* ");
    if (!doNorm) g_dummy->GetYaxis()->SetTitle("d#sigma^{fid}/d#phi* (pb)");
    g_dummy->GetYaxis()->SetTitleSize(0.05);
    g_dummy->GetYaxis()->SetLabelSize(0.05);
    g_dummy->SetFillColor(kWhite);
    g_dummy->SetTitle(0);
    return g_dummy;
}

TGraphAsymmErrors* CreateCopy(TGraphAsymmErrors* graph, double scale = 1.0) {
    double x, y, errorl, errorh;
    TGraphAsymmErrors* g_copy = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        graph->GetPoint(iphistar, x, y);
        errorl = graph->GetErrorYlow(iphistar);
        errorh = graph->GetErrorYhigh(iphistar);
        g_copy->SetPoint(iphistar, x, y * scale);
        g_copy->SetPointError(iphistar, 0, 0, errorl*scale, errorh * scale);
    }
    return g_copy;
}

void MakeFinalPlots() {

    AverageData = false;
    TGraphAsymmErrors* g_data_final;
    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_ph_final;

    if (Type == "elec") {
        std::string textn = "Data_Graph_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        if (doMG) textn += "MG_";
        else textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile dg(textn.c_str());
        TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
        g_data_final = CreateCopy(g_data_temp);
        textn = "Data_Graph_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "MG_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile mg(textn.c_str());
        TGraphAsymmErrors* g_mg_temp = (TGraphAsymmErrors*) mg.Get("Graph");
        g_mg_final = CreateCopy(g_mg_temp);
        textn = "Data_Graph_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "PH_";
        if (elec == 0)textn += "Dressed.root";
        if (elec == 1)textn += "Born.root";
        if (elec == 2)textn += "Naked.root";
        cout << textn.c_str() << endl;
        TFile fph(textn.c_str());
        TGraphAsymmErrors* g_ph_temp = (TGraphAsymmErrors*) fph.Get("Graph");
        g_ph_final = CreateCopy(g_ph_temp);
    }

    if (Type == "combined") {

        std::string textn;
        if (doNorm)textn = "/home/user1/lesko/work/HomeWork/Phistar/CombineElectWithMu/Results/Comb_Norm_UsingBlue1D.root";
        else textn = "/home/user1/lesko/work/HomeWork/Phistar/CombineElectWithMu/Results/Comb_Abs_UsingBlue1D.root";

        TFile tr(textn.c_str());

        TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) tr.Get("h_Comb");
        if (!g_data_temp) {
            cout << "couldn't find BLUE results" << endl;
            return;
        }

        g_data_final = CreateCopy(g_data_temp);
        textn = "Comb_Hist_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "MG_";
        textn += "Dressed.root";
        TFile tr2(textn.c_str());
        TH1D* h_mg_temp = (TH1D*) tr2.Get("h_Comb");
        TGraphAsymmErrors* g_mg_temp = ConvertToTGraph(h_mg_temp);
        g_mg_final = CreateCopy(g_mg_temp);
        textn = "Comb_Hist_MC_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "PH_";
        textn += "Dressed.root";
        TFile tr3(textn.c_str());
        TH1D* h_ph_temp = (TH1D*) tr3.Get("h_Comb");
        TGraphAsymmErrors* g_ph_temp = ConvertToTGraph(h_ph_temp);
        g_ph_final = CreateCopy(g_ph_temp);
    }

    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_dummy_phistar = CreateDummy(g_data_final);
    TGraphAsymmErrors* g_ratio_phistar = CreateRatio(g_data_final, g_data_final, 1);
    TGraphAsymmErrors* g_ratio_mg_phistar = CreateRatio(g_data_final, g_mg_final, 0);
    TGraphAsymmErrors* g_ratio_ph_phistar = CreateRatio(g_data_final, g_ph_final, 0);

    PrintFinal(g_data_final, g_mg_final, g_ph_final);
    PlotFinal(g_data_final, g_mg_final, g_ph_final, g_dummy_phistar, g_ratio_phistar, g_ratio_mg_phistar, g_ratio_ph_phistar);
}

TH1D* GetHistMuon(std::string name) {
    TH1D* h_Muon = new TH1D("h_Comb", "h_Comb", nphistar, phistarBins);
    h_Muon->Sumw2();
    TFile fm(name.c_str());
    TH1F *h_temp = (TH1F*) fm.Get("hRecoClone1");
    for (uint i = 0; i < nphistar; i++) {
        h_Muon->SetBinContent(i + 1, h_temp->GetBinContent(i + 1));
        h_Muon->SetBinError(i + 1, h_temp->GetBinError(i + 1));
    }
    return h_Muon;
}

void MakeFinalPlots2(bool RemoveLumiUncertainty = 0) {

    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_data_muon;
    TGraphAsymmErrors* g_data_elec;

    std::string textn = "Data_Graph_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "PH_";
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    cout << textn.c_str() << endl;
    TFile dg(textn.c_str());
    TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
    g_data_elec = CreateCopy(g_data_temp);
    g_data_elec = SplittGraph(false, RemoveLumiUncertainty); //Overriting for testing purposes
    textn = "Data_Graph_MC_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn = "Comb_Hist_MC_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "MG_";
    textn += "Dressed.root";
    TFile tr2(textn.c_str());
    TH1D* h_mg_temp = (TH1D*) tr2.Get("h_Comb");
    TGraphAsymmErrors* g_mg_temp = ConvertToTGraph(h_mg_temp);
    g_mg_final = CreateCopy(g_mg_temp);
    //   textn+="MG_";
    // if (elec==0)textn+="Dressed.root";
    // if (elec==1)textn+="Born.root";
    // if (elec==2)textn+="Naked.root";
    // cout<<textn.c_str()<<endl;
    // TFile mg(textn.c_str());  
    // TGraphAsymmErrors* g_mg_temp= (TGraphAsymmErrors*)mg.Get("Graph");
    // g_mg_final=CreateCopy(g_mg_temp);

    //std::string muon_name = "Muon_CrossSection_For_Combination/Madgraph_Dressed_";
    //if (doNorm) muon_name += "Normal_";
    //else muon_name += "Absolute_";
    //muon_name += "Central_Full_Errors.root";
    //if (doNorm);
    //else muon_name = "Muon_CrossSection_For_Combination/Powheg_Born_Absolute_Central_Full_Errors.root";
    TFile MuonFile("Muon_CrossSection_For_Combination/1D_unf_dist.root");
    TH1D* h_Muon_temp;
    if (doNorm)h_Muon_temp = (TH1D*) MuonFile.Get("hnorm");
    else h_Muon_temp = (TH1D*) MuonFile.Get("habs");

    if (!doNorm) {
        h_Muon_temp->Scale(1. / 1000.);
    }
    g_data_muon = ConvertToTGraph(h_Muon_temp);
    g_data_muon = SplittGraph(true, RemoveLumiUncertainty);

    if (RemoveLumiUncertainty) {
        for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
            double x, y;
            g_data_muon->GetPoint(iphistar, x, y);
            g_data_elec->GetPoint(iphistar, x, y);

        }
    }
    cout << "going to make ratio plots" << endl;
    if (doNorm)textn = "~/work/HomeWork/Phistar/CombineElectWithMu/Comb_Norm_UsingBlue1D.root";
    else textn = "~/work/HomeWork/Phistar/CombineElectWithMu/Comb_Abs_UsingBlue1D.root";
    TFile tr_comb(textn.c_str());
    TGraphAsymmErrors* g_comb_temp = (TGraphAsymmErrors*) tr_comb.Get("h_Comb");
    TGraphAsymmErrors* g_comb_final = CreateCopy(g_comb_temp);


    TGraphAsymmErrors* g_dummy_phistar = CreateDummy(g_comb_final);
    TGraphAsymmErrors* g_ratio_elec = CreateRatio(g_comb_final, g_data_elec, 0);
    TGraphAsymmErrors* g_ratio_mg = CreateRatio(g_comb_final, g_comb_final, 1);
    TGraphAsymmErrors* g_ratio_muon = CreateRatio(g_comb_final, g_data_muon, 0);
    //
    //  PrintFinal(g_data_final,g_mg_final,g_ph_final);

    PlotFinal(g_comb_final, g_data_elec, g_data_muon, g_dummy_phistar, g_ratio_mg, g_ratio_elec, g_ratio_muon, 1);



    for (size_t phibinN = 0; phibinN < nphistar; phibinN++) {
        double xe, ye;
        double xc, yc;
        g_data_elec->GetPoint(phibinN, xe, ye);
        g_comb_final->GetPoint(phibinN, xc, yc);
        cout << "AND FOR OUR RATIO IN BIN " << phibinN << " we have " << ye / yc << endl;
    }
    PrintFinal(g_comb_final, g_data_elec, g_data_muon, 1);
}
