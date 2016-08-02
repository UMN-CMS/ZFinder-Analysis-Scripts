#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TMatrixDSparse.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "TRandom.h"
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

const bool AlexPlots = false;
const bool doNorm = 1;
const bool debug = true;
const int elec = 0;
const int doMG = 0;
const std::string Tag = "";
const std::string Type = "elec"; //elec, muon or combined

const double phistarBins[] = {0.000, 0.004, 0.008, 0.012, 0.016, 0.020, 0.024, 0.029, 0.034, 0.039, 0.045, 0.052, 0.057, 0.064, 0.072, 0.081, 0.091, 0.102, 0.114, 0.128, 0.145, 0.165, 0.189, 0.219, 0.258, 0.312, 0.391, 0.524, 0.695, 0.918, 1.153, 1.496, 1.947, 2.522, 3.277};
size_t nphistar = (sizeof (phistarBins) / sizeof (phistarBins[0])) - 1;
const double PHDataUnfoldUncertainty[] = {.0082, .0093, .0099, .0103, .0089, .0096, .0091, .0084, .0088, .0077, .0077, .0099, .0081, .0078, .0077, .0075, .0076,
    .0078, .0076, .0076, .0081, .0076, .0076, .0076, .0074, .0075, .0081, .0096, .0122, .0166, .0188, .0234, .0299, .0364};

void ChiSquared(TGraphAsymmErrors* g_data_final, TGraphAsymmErrors* g_mg_final, TGraphAsymmErrors* g_ph_final) {

    std::string CovarianceMatrixEff = "CovarianceMatrixAbs_MG_Dressed.root";
    TFile* theFile = new TFile(CovarianceMatrixEff.c_str());

    TMatrixD* UnFoldMatrix = (TMatrixD*) theFile->Get("Unfolded");
    TMatrixD* EffMatrix = (TMatrixD*) theFile->Get("Efficency");
    TMatrixD* MCstatMatrix = (TMatrixD*) theFile->Get("MCstat");
    TMatrixD* BGMatrix = (TMatrixD*) theFile->Get("BG");
    TMatrixD* PileupMatrix = (TMatrixD*) theFile->Get("PileUp");
    TMatrixD* LumiMatrix = (TMatrixD*) theFile->Get("Lumi");
    TMatrixD* PtMatrix = (TMatrixD*) theFile->Get("Pt");

    TMatrixD PilPlusLumi(*LumiMatrix + *PileupMatrix);
    cout << "before :" << PilPlusLumi(1, 1) << endl;
    PilPlusLumi.Invert();
    cout << "after :" << PilPlusLumi(1, 1) << endl;
    TMatrixD PilPlusLumiInverted(PilPlusLumi);
    double chiSquared = 0;
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        double x, data1, MC1;
        double hold = 0;
        for (size_t iphistar2 = 0; iphistar2 < nphistar; iphistar2++) {
            double data2, MC2;
            g_data_final->GetPoint(iphistar2, x, data2);
            g_mg_final->GetPoint(iphistar2, x, MC2);
            hold += PilPlusLumiInverted(iphistar, iphistar2)*(data2 - MC2);
        }
        g_data_final->GetPoint(iphistar, x, data1);
        g_mg_final->GetPoint(iphistar, x, MC1);
        chiSquared += hold * (data1 - MC1);
    }
    cout << " and we have a chiSquared of " << chiSquared << endl;
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
            errorlmc = graph->GetErrorYlow(iphistar);
            errorhmc = graph->GetErrorYhigh(iphistar);

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
        g_ratio->GetYaxis()->SetRangeUser(0.8, 1.2);
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

    TCanvas* FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot", 800, 900);
    FinalPhiTot->Divide(1, 2);
    FinalPhiTot->cd(1);
    gPad->SetPad("p1", "p1", 0, 2.5 / 9.0, 1, 1, kWhite, 0, 0);
    gPad->SetBottomMargin(0.01);
    gPad->SetTopMargin(0.06);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetLogx(1);
    gPad->SetLogy(1);
    TGraphAsymmErrors* g_ph2_finalAlex;
    TGraphAsymmErrors* g_ratio_ph2_phistarAlex;
    TGraphAsymmErrors* g_ph1_finalAlex;
    TGraphAsymmErrors* g_ratio_ph1_phistarAlex;
    if (AlexPlots && doNorm) {
        if (debug) {
            cout << "test 1" << endl;
        }
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
    g_dummy_phistar->Draw("A2");
    int Color1 = kBlue;
    int Color2 = kRed;
    int Color3 = kGreen + 2;
    int Color4 = kOrange + 7;
    g_mg_final->SetMarkerColor(Color1);
    g_mg_final->SetLineColor(Color1);
    g_mg_final->SetMarkerSize(1);
    g_mg_final->SetLineWidth(2);
    g_mg_final->SetMarkerStyle(4);
    g_mg_final->Draw("PEsame");
    g_ph_final->SetMarkerColor(Color2);
    g_ph_final->SetLineColor(Color2);
    g_ph_final->SetMarkerSize(1);
    g_ph_final->SetLineWidth(2);
    g_ph_final->SetMarkerStyle(20);
    g_ph_final->Draw("PEsame");
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
    int ColorError = kGray + 1;
    g_data_final->SetFillColor(ColorError);
    g_data_final->SetMarkerSize(1);
    g_data_final->SetLineWidth(2);
    g_data_final->SetMarkerStyle(20);
    g_data_final->Draw("PEsame");
    g_data_final->SetFillColor(ColorError);

    TLegend* leg = new TLegend(0.23, 0.76, 0.95, 0.94);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetLineWidth(1);
    leg->SetNColumns(1);
    leg->SetTextFont(42);

    if (!isPlot2) {
        leg->AddEntry(g_data_final, "2012 data", "PEF");
        if (Type == "elec") {
            leg->AddEntry(g_mg_final, "Z #rightarrow ee MadGraph+Pythia6 (Z2star)", "P");
            leg->AddEntry(g_ph_final, "Z #rightarrow ee POWHEG+Pythia6 (Z2star)", "P");
            if (AlexPlots && doNorm) {
                leg->AddEntry(g_ph2_finalAlex, "Z #rightarrow ee POWHEG+Pythia8 (Tunepp 5)", "P");
                leg->AddEntry(g_ph1_finalAlex, "Z #rightarrow ee POWHEG+Pythia8 (Tunepp 14)", "P");
            }
        }
        if (Type == "muon") {
            leg->AddEntry(g_mg_final, "Z #rightarrow #mu#mu MadGraph", "P");
            leg->AddEntry(g_ph_final, "Z #rightarrow #mu#mu Powheg", "P");
        }
        if (Type == "combined") {
            leg->AddEntry(g_mg_final, "Z #rightarrow ll MadGraph", "P");
            leg->AddEntry(g_ph_final, "Z #rightarrow ll Powheg", "P");
        }
    } else {
        leg->AddEntry(g_mg_final, "2012 data: Z #rightarrow ee", "P");
        leg->AddEntry(g_ph_final, "2012 data: Z #rightarrow #mu#mu", "P");
        leg->AddEntry(g_data_final, "MadGraph: Z #rightarrow ll", "PEF");
    }
    leg->Draw();

    TLatex mark;
    // mark.SetTextSize(0.05);
    mark.SetTextSize(0.035);
    mark.SetNDC(true);
    mark.DrawLatex(0.745, 0.95, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.19, 0.955, "CMS Preliminary");
    if (Type == "elec" && !isPlot2) {
        mark.DrawLatex(0.19, 0.20, "|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
        mark.DrawLatex(0.19, 0.13, "p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
        mark.DrawLatex(0.19, 0.06, "60 GeV < M_{ee} < 120 GeV");
    }
    if (Type == "muon" && !isPlot2) {
        mark.DrawLatex(0.19, 0.20, "|#eta^{#mu_{0}}| < 2.1,        |#eta^{#mu_{1}}| < 2.4");
        mark.DrawLatex(0.19, 0.13, "p_{T}^{#mu_{0}} > 30 GeV,   p_{T}^{#mu_{1}} > 20 GeV");
        mark.DrawLatex(0.19, 0.06, "60 GeV < M_{#mu#mu} < 120 GeV");
    }
    if (Type == "combined" || isPlot2) {
        mark.DrawLatex(0.19, 0.20, "|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
        mark.DrawLatex(0.19, 0.13, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
        mark.DrawLatex(0.19, 0.06, "60 GeV < M_{ll} < 120 GeV");
    }
    FinalPhiTot->cd(2);
    gPad->SetPad("p2", "p2", 0, 0, 1, 2.5 / 9.0, kWhite, 0, 0);
    gPad->SetBottomMargin(0.37);
    gPad->SetTopMargin(0.01);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.06);
    gPad->SetLogx(1);

    if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("Data/MC  ");
    else g_ratio_phistar->GetYaxis()->SetTitle("MC/Data   ");
    g_ratio_phistar->GetXaxis()->SetRangeUser(0.001, 3.0);
    g_ratio_phistar->GetYaxis()->SetRangeUser(0.76, 1.24);
    if (isPlot2) g_ratio_phistar->GetYaxis()->SetRangeUser(0.88, 1.12);
    if (isPlot2 && (!doNorm)) g_ratio_phistar->GetYaxis()->SetRangeUser(0.88, 1.22);
    g_ratio_phistar->GetYaxis()->SetTitleOffset(0.45);
    g_ratio_phistar->SetFillColor(ColorError);
    g_ratio_phistar->Draw("AE2leg->AddEn");
    g_ratio_mg_phistar->SetMarkerSize(1);
    g_ratio_mg_phistar->SetLineWidth(2);
    g_ratio_mg_phistar->SetMarkerStyle(4);
    g_ratio_mg_phistar->SetMarkerColor(Color1);
    g_ratio_mg_phistar->SetLineColor(Color1);
    g_ratio_mg_phistar->Draw("PEsame");
    g_ratio_ph_phistar->SetMarkerSize(1);
    g_ratio_ph_phistar->SetLineWidth(2);
    g_ratio_ph_phistar->SetMarkerStyle(20);
    g_ratio_ph_phistar->SetMarkerColor(Color2);
    g_ratio_ph_phistar->SetLineColor(Color2);
    g_ratio_ph_phistar->Draw("PEsame");
    if (AlexPlots && doNorm) {
        g_ratio_ph2_phistarAlex->SetMarkerColor(Color3);
        g_ratio_ph2_phistarAlex->SetLineColor(Color3);
        g_ratio_ph2_phistarAlex->SetMarkerSize(1);
        g_ratio_ph2_phistarAlex->SetLineWidth(2);
        g_ratio_ph2_phistarAlex->SetMarkerStyle(21);
        g_ratio_ph2_phistarAlex->Draw("PEsame");


        g_ratio_ph1_phistarAlex->SetMarkerColor(Color4);
        g_ratio_ph1_phistarAlex->SetLineColor(Color4);
        g_ratio_ph1_phistarAlex->SetMarkerSize(1);
        g_ratio_ph1_phistarAlex->SetLineWidth(2);
        g_ratio_ph1_phistarAlex->SetMarkerStyle(22);
        g_ratio_ph1_phistarAlex->Draw("PEsame");
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

    TCanvas* FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio", 800, 900);
    FinalPhiRatio->cd();
    FinalPhiRatio->SetLogx();
    if (isPlot2) g_ratio_phistar->GetYaxis()->SetTitle("Data/MC");
    else g_ratio_phistar->GetYaxis()->SetTitle("MC/Data");
    g_ratio_phistar->SetFillColor(ColorError);
    g_ratio_phistar->GetYaxis()->SetRangeUser(0.7, 1.3);
    g_ratio_phistar->GetYaxis()->SetTitleOffset(1.2);
    g_ratio_phistar->GetYaxis()->SetTitleSize(0.04);
    g_ratio_phistar->GetYaxis()->SetLabelSize(0.04);
    g_ratio_phistar->GetXaxis()->SetTitleSize(0.04);
    g_ratio_phistar->GetXaxis()->SetLabelSize(0.04);
    g_ratio_phistar->GetXaxis()->SetLabelOffset(-0.01);
    g_ratio_phistar->Draw("AE2");

    g_ratio_mg_phistar->Draw("PEsame");
    g_ratio_ph_phistar->Draw("PEsame");
    if (AlexPlots && doNorm) {
        g_ratio_ph2_phistarAlex->Draw("PEsame");
        g_ratio_ph1_phistarAlex->Draw("PEsame");
    }
    mark.SetTextSize(0.03);
    mark.DrawLatex(0.7, 0.907, "19.7 fb^{-1} (8 TeV)");
    mark.DrawLatex(0.19, 0.87, "CMS Preliminary");
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
        mark.DrawLatex(0.15, 0.20, "p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
        mark.DrawLatex(0.15, 0.15, "60 GeV < M_{ll} < 120 GeV");
    }
    TLegend* leg2 = new TLegend(0.13, 0.62, 0.95, 0.84);
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->SetLineWidth(1);
    leg2->SetNColumns(1);
    leg2->SetTextFont(42);

    if (!isPlot2) {
        leg2->AddEntry(g_ratio_phistar, "2012 data", "F");
        if (Type == "elec") {
            leg2->AddEntry(g_mg_final, "Z #rightarrow ee MadGraph+Pythia6 (Z2star)", "P");
            leg2->AddEntry(g_ph_final, "Z #rightarrow ee POWHEG+Pythia6 (Z2star)", "P");
            if (AlexPlots && doNorm) {
                leg2->AddEntry(g_ratio_ph2_phistarAlex, "Z #rightarrow ee POWHEG+Pythia8  (Tunepp 5)", "P");
                leg2->AddEntry(g_ratio_ph1_phistarAlex, "Z #rightarrow ee POWHEG+Pythia8 (Tunepp 14)", "P");
            }
        }
        if (Type == "muon") {
            leg2->AddEntry(g_mg_final, "Z #rightarrow #mu#mu MadGraph", "P");
            leg2->AddEntry(g_ph_final, "Z #rightarrow #mu#mu Powheg", "P");
        }
        if (Type == "combined") {
            leg2->AddEntry(g_mg_final, "Z #rightarrow ll MadGraph", "P");
            leg2->AddEntry(g_ph_final, "Z #rightarrow ll Powheg", "P");
        }
    } else {
        leg2->AddEntry(g_mg_final, "2012 data: Z #rightarrow ee", "P");
        leg2->AddEntry(g_ph_final, "2012 data: Z #rightarrow #mu#mu", "P");
        leg2->AddEntry(g_data_final, "MadGraph: Z #rightarrow ll", "PEF");
    }
    leg2->Draw();
    plotname = "ZShape_Ratio";
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
    g_dummy->GetXaxis()->SetTitle(0);
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

TGraphAsymmErrors* ConvertToTGraph(TH1D* h) {
    TGraphAsymmErrors* g = new TGraphAsymmErrors(nphistar);
    for (size_t iphistar = 0; iphistar < nphistar; iphistar++) {
        g->SetPoint(iphistar, (phistarBins[iphistar] + phistarBins[iphistar + 1]) / 2., h->GetBinContent(iphistar + 1));
        g->SetPointError(iphistar, 0, 0, h->GetBinError(iphistar + 1), h->GetBinError(iphistar + 1));
    }
    return g;
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
        std::string textn = "Comb_Hist_";
        textn += Tag;
        if (doNorm) textn += "Norm_";
        else textn += "Abs_";
        textn += "MG_";
        textn += "Dressed.root";
        TFile tr(textn.c_str());
        TH1D* h_data_temp = (TH1D*) tr.Get("h_Comb");
        TGraphAsymmErrors* g_data_temp = ConvertToTGraph(h_data_temp);
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
    ChiSquared(g_data_final, g_mg_final, g_ph_final);
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

void MakeFinalPlots2() {

    TGraphAsymmErrors* g_mg_final;
    TGraphAsymmErrors* g_data_muon;
    TGraphAsymmErrors* g_data_elec;

    std::string textn = "Data_Graph_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "MG_";
    if (elec == 0)textn += "Dressed.root";
    if (elec == 1)textn += "Born.root";
    if (elec == 2)textn += "Naked.root";
    cout << textn.c_str() << endl;
    TFile dg(textn.c_str());
    TGraphAsymmErrors* g_data_temp = (TGraphAsymmErrors*) dg.Get("Graph");
    g_data_elec = CreateCopy(g_data_temp);
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
    std::string muon_name = "Muon_CrossSection_For_Combination/Madgraph_Dressed_";
    if (doNorm) muon_name += "Normal_";
    else muon_name += "Absolute_";
    muon_name += "Central_Full_Errors.root";
    TH1D* h_Muon_temp = GetHistMuon(muon_name);
    if (!doNorm) {
        h_Muon_temp->Scale(1. / 1000.);
    }
    g_data_muon = ConvertToTGraph(h_Muon_temp);

    cout << "going to make ratio plots" << endl;
    TGraphAsymmErrors* g_dummy_phistar = CreateDummy(g_mg_final);
    TGraphAsymmErrors* g_ratio_elec = CreateRatio(g_mg_final, g_data_elec, 0);
    TGraphAsymmErrors* g_ratio_mg = CreateRatio(g_mg_final, g_mg_final, 1);
    TGraphAsymmErrors* g_ratio_muon = CreateRatio(g_mg_final, g_data_muon, 0);
    //
    //  PrintFinal(g_data_final,g_mg_final,g_ph_final);
    PlotFinal(g_mg_final, g_data_elec, g_data_muon, g_dummy_phistar, g_ratio_mg, g_ratio_elec, g_ratio_muon, 1);
    textn = "Comb_Hist_";
    textn += Tag;
    if (doNorm) textn += "Norm_";
    else textn += "Abs_";
    textn += "MG_";
    textn += "Dressed.root";
    TFile tr_comb(textn.c_str());
    TH1D* h_comb_temp = (TH1D*) tr_comb.Get("h_Comb");
    TGraphAsymmErrors* g_comb_temp = ConvertToTGraph(h_comb_temp);
    TGraphAsymmErrors* g_comb_final = CreateCopy(g_comb_temp);

    PrintFinal(g_comb_final, g_data_elec, g_data_muon, 1);
}

