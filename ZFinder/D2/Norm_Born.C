{
//=========Macro generated from canvas: FinalPhiTot/FinalPhiTot
//=========  (Mon Jul 18 12:04:03 2016) by ROOT version5.34/18
   TCanvas *FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot",1684,48,800,900);
   FinalPhiTot->Range(0,0,1,1);
   FinalPhiTot->SetFillColor(0);
   FinalPhiTot->SetBorderMode(0);
   FinalPhiTot->SetBorderSize(2);
   FinalPhiTot->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: p1
   TPad *p1 = new TPad("p1", "p1",0,0.2777778,1,1);
   p1->Draw();
   p1->cd();
   p1->Range(-0.4556962,0.6817436,2.582278,2.404389);
   p1->SetFillColor(0);
   p1->SetBorderMode(0);
   p1->SetBorderSize(0);
   p1->SetLogy();
   p1->SetLeftMargin(0.15);
   p1->SetRightMargin(0.06);
   p1->SetTopMargin(0.06);
   p1->SetBottomMargin(0.01);
   p1->SetFrameBorderMode(0);
   p1->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.2,70.30767);
   grae->SetPointError(0,0,0,0.1324974,0.1324974);
   grae->SetPoint(1,0.6,69.33778);
   grae->SetPointError(1,0,0,0.1306982,0.1306982);
   grae->SetPoint(2,1,67.85886);
   grae->SetPointError(2,0,0,0.1412105,0.1412105);
   grae->SetPoint(3,1.4,59.04005);
   grae->SetPointError(3,0,0,0.1418763,0.1418763);
   grae->SetPoint(4,1.8,38.86572);
   grae->SetPointError(4,0,0,0.135296,0.135296);
   grae->SetPoint(5,2.2,9.376096);
   grae->SetPointError(5,0,0,0.3236447,0.3236447);
   
   TH1F *Graph_Graph959 = new TH1F("Graph_Graph959","Graph",100,0,2.4);
   Graph_Graph959->SetMinimum(5);
   Graph_Graph959->SetMaximum(200);
   Graph_Graph959->SetDirectory(0);
   Graph_Graph959->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph959->SetLineColor(ci);
   Graph_Graph959->GetXaxis()->SetRange(1,100);
   Graph_Graph959->GetXaxis()->SetLabelFont(42);
   Graph_Graph959->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph959->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph959->GetXaxis()->SetTitleFont(42);
   Graph_Graph959->GetYaxis()->SetLabelFont(42);
   Graph_Graph959->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph959->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph959->GetYaxis()->SetTitleFont(42);
   Graph_Graph959->GetZaxis()->SetLabelFont(42);
   Graph_Graph959->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph959->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph959->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph959);
   
   grae->Draw("a2");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(4);
   grae->SetPoint(0,0.2,68.44178);
   grae->SetPointError(0,0,0,0.09457353,0.09457353);
   grae->SetPoint(1,0.6,68.48141);
   grae->SetPointError(1,0,0,0.09464392,0.09464392);
   grae->SetPoint(2,1,67.67099);
   grae->SetPointError(2,0,0,0.09412069,0.09412069);
   grae->SetPoint(3,1.4,59.9608);
   grae->SetPointError(3,0,0,0.08887153,0.08887153);
   grae->SetPoint(4,1.8,41.24982);
   grae->SetPointError(4,0,0,0.07383172,0.07383172);
   grae->SetPoint(5,2.2,10.25012);
   grae->SetPointError(5,0,0,0.03692826,0.03692826);
   
   TH1F *Graph_Graph960 = new TH1F("Graph_Graph960","Graph",100,0,2.4);
   Graph_Graph960->SetMinimum(4.376904);
   Graph_Graph960->SetMaximum(74.41234);
   Graph_Graph960->SetDirectory(0);
   Graph_Graph960->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph960->SetLineColor(ci);
   Graph_Graph960->GetXaxis()->SetLabelFont(42);
   Graph_Graph960->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph960->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph960->GetXaxis()->SetTitleFont(42);
   Graph_Graph960->GetYaxis()->SetLabelFont(42);
   Graph_Graph960->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph960->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph960->GetYaxis()->SetTitleFont(42);
   Graph_Graph960->GetZaxis()->SetLabelFont(42);
   Graph_Graph960->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph960->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph960->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph960);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.2,68.74086);
   grae->SetPointError(0,0,0,0.2070766,0.2070766);
   grae->SetPoint(1,0.6,68.25162);
   grae->SetPointError(1,0,0,0.1441424,0.1441424);
   grae->SetPoint(2,1,66.92132);
   grae->SetPointError(2,0,0,0.06730436,0.06730436);
   grae->SetPoint(3,1.4,58.06196);
   grae->SetPointError(3,0,0,0.0934479,0.0934479);
   grae->SetPoint(4,1.8,39.24395);
   grae->SetPointError(4,0,0,0.1163179,0.1163179);
   grae->SetPoint(5,2.2,9.449127);
   grae->SetPointError(5,0,0,0.04228468,0.04228468);
   
   TH1F *Graph_Graph961 = new TH1F("Graph_Graph961","Graph",100,0,2.4);
   Graph_Graph961->SetMinimum(3.452733);
   Graph_Graph961->SetMaximum(74.90204);
   Graph_Graph961->SetDirectory(0);
   Graph_Graph961->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph961->SetLineColor(ci);
   Graph_Graph961->GetXaxis()->SetLabelFont(42);
   Graph_Graph961->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph961->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph961->GetXaxis()->SetTitleFont(42);
   Graph_Graph961->GetYaxis()->SetLabelFont(42);
   Graph_Graph961->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph961->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph961->GetYaxis()->SetTitleFont(42);
   Graph_Graph961->GetZaxis()->SetLabelFont(42);
   Graph_Graph961->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph961->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph961->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph961);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.2,70.30767);
   grae->SetPointError(0,0,0,0.1324974,0.1324974);
   grae->SetPoint(1,0.6,69.33778);
   grae->SetPointError(1,0,0,0.1306982,0.1306982);
   grae->SetPoint(2,1,67.85886);
   grae->SetPointError(2,0,0,0.1412105,0.1412105);
   grae->SetPoint(3,1.4,59.04005);
   grae->SetPointError(3,0,0,0.1418763,0.1418763);
   grae->SetPoint(4,1.8,38.86572);
   grae->SetPointError(4,0,0,0.135296,0.135296);
   grae->SetPoint(5,2.2,9.376096);
   grae->SetPointError(5,0,0,0.3236447,0.3236447);
   
   TH1F *Graph_Graph_Graph959962 = new TH1F("Graph_Graph_Graph959962","Graph",100,0,2.4);
   Graph_Graph_Graph959962->SetMinimum(5);
   Graph_Graph_Graph959962->SetMaximum(200);
   Graph_Graph_Graph959962->SetDirectory(0);
   Graph_Graph_Graph959962->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph959962->SetLineColor(ci);
   Graph_Graph_Graph959962->GetXaxis()->SetRange(1,100);
   Graph_Graph_Graph959962->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph959962->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph959962->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph959962->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph959962->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph959962->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph959962->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph959962->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph959962->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph959962->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph959962->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph959962->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph959962);
   
   grae->Draw("pe");
   
   TLegend *leg = new TLegend(0.23,0.76,0.95,0.94,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph0","2012 data","PEF");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph1","Z #rightarrow ee MadGraph+Pythia6 (Z2star)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(4);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","Z #rightarrow ee POWHEG+Pythia6 (Z2star)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   TLatex *   tex = new TLatex(0.745,0.95,"19.7 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.955,"CMS Preliminary");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.2,"|#eta^{e_{0}}| < 2.1,        |#eta^{e_{1}}| < 2.4");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.13,"p_{T}^{e_{0}} > 30 GeV,   p_{T}^{e_{1}} > 20 GeV");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.06,"60 GeV < M_{ee} < 120 GeV");
tex->SetNDC();
   tex->SetTextSize(0.035);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TPaveText *pt = new TPaveText(0.4315829,0.9359316,0.5684171,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("Graph");
   pt->Draw();
   p1->Modified();
   FinalPhiTot->cd();
  
// ------------>Primitives in pad: p2
   p2 = new TPad("p2", "p2",0,0,1,0.2777778);
   p2->Draw();
   p2->cd();
   p2->Range(-0.4561519,0.5612903,2.584861,1.206452);
   p2->SetFillColor(0);
   p2->SetBorderMode(0);
   p2->SetBorderSize(0);
   p2->SetLeftMargin(0.15);
   p2->SetRightMargin(0.06);
   p2->SetTopMargin(0.01);
   p2->SetBottomMargin(0.37);
   p2->SetFrameBorderMode(0);
   p2->SetFrameBorderMode(0);
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetPoint(0,0.2,1);
   grae->SetPointError(0,0.2,0.2,0.001884537,0.001884537);
   grae->SetPoint(1,0.6,1);
   grae->SetPointError(1,0.2,0.2,0.00188495,0.00188495);
   grae->SetPoint(2,1,1);
   grae->SetPointError(2,0.2,0.2,0.002080945,0.002080945);
   grae->SetPoint(3,1.4,1);
   grae->SetPointError(3,0.2,0.2,0.002403052,0.002403052);
   grae->SetPoint(4,1.8,1);
   grae->SetPointError(4,0.2,0.2,0.003481113,0.003481113);
   grae->SetPoint(5,2.2,1);
   grae->SetPointError(5,0.2,0.2,0.03451806,0.03451806);
   
   TH1F *Graph_Graph963 = new TH1F("Graph_Graph963","Graph",100,0,2.64);
   Graph_Graph963->SetMinimum(0.8);
   Graph_Graph963->SetMaximum(1.2);
   Graph_Graph963->SetDirectory(0);
   Graph_Graph963->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph963->SetLineColor(ci);
   Graph_Graph963->GetXaxis()->SetRange(1,91);
   Graph_Graph963->GetXaxis()->SetLabelFont(42);
   Graph_Graph963->GetXaxis()->SetLabelOffset(-0.01);
   Graph_Graph963->GetXaxis()->SetTitleFont(42);
   Graph_Graph963->GetYaxis()->SetTitle("MC/Data");
   Graph_Graph963->GetYaxis()->SetLabelFont(42);
   Graph_Graph963->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph963->GetYaxis()->SetTitleFont(42);
   Graph_Graph963->GetZaxis()->SetLabelFont(42);
   Graph_Graph963->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph963->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph963->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph963);
   
   grae->Draw("ae2");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(4);
   grae->SetPoint(0,0.2,0.973461);
   grae->SetPointError(0,0,0,0.001345138,0.001345138);
   grae->SetPoint(1,0.6,0.9876494);
   grae->SetPointError(1,0,0,0.001364969,0.001364969);
   grae->SetPoint(2,1,0.9972316);
   grae->SetPointError(2,0,0,0.001387007,0.001387007);
   grae->SetPoint(3,1.4,1.015595);
   grae->SetPointError(3,0,0,0.001505275,0.001505275);
   grae->SetPoint(4,1.8,1.061342);
   grae->SetPointError(4,0,0,0.001899662,0.001899662);
   grae->SetPoint(5,2.2,1.093218);
   grae->SetPointError(5,0,0,0.003938553,0.003938553);
   
   TH1F *Graph_Graph964 = new TH1F("Graph_Graph964","Graph",100,0,2.4);
   Graph_Graph964->SetMinimum(0.9596118);
   Graph_Graph964->SetMaximum(1.109661);
   Graph_Graph964->SetDirectory(0);
   Graph_Graph964->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph964->SetLineColor(ci);
   Graph_Graph964->GetXaxis()->SetLabelFont(42);
   Graph_Graph964->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph964->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph964->GetXaxis()->SetTitleFont(42);
   Graph_Graph964->GetYaxis()->SetLabelFont(42);
   Graph_Graph964->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph964->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph964->GetYaxis()->SetTitleFont(42);
   Graph_Graph964->GetZaxis()->SetLabelFont(42);
   Graph_Graph964->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph964->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph964->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph964);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.2,0.9777149);
   grae->SetPointError(0,0,0,0.002945292,0.002945292);
   grae->SetPoint(1,0.6,0.9843353);
   grae->SetPointError(1,0,0,0.002078844,0.002078844);
   grae->SetPoint(2,1,0.986184);
   grae->SetPointError(2,0,0,0.0009918287,0.0009918287);
   grae->SetPoint(3,1.4,0.9834334);
   grae->SetPointError(3,0,0,0.001582788,0.001582788);
   grae->SetPoint(4,1.8,1.009732);
   grae->SetPointError(4,0,0,0.002992815,0.002992815);
   grae->SetPoint(5,2.2,1.007789);
   grae->SetPointError(5,0,0,0.004509838,0.004509838);
   
   TH1F *Graph_Graph965 = new TH1F("Graph_Graph965","Graph",100,0,2.4);
   Graph_Graph965->SetMinimum(0.9709741);
   Graph_Graph965->SetMaximum(1.01652);
   Graph_Graph965->SetDirectory(0);
   Graph_Graph965->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph965->SetLineColor(ci);
   Graph_Graph965->GetXaxis()->SetLabelFont(42);
   Graph_Graph965->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph965->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph965->GetXaxis()->SetTitleFont(42);
   Graph_Graph965->GetYaxis()->SetLabelFont(42);
   Graph_Graph965->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph965->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph965->GetYaxis()->SetTitleFont(42);
   Graph_Graph965->GetZaxis()->SetLabelFont(42);
   Graph_Graph965->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph965->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph965->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph965);
   
   grae->Draw("pe");
   
   pt = new TPaveText(0.4680151,0.9355505,0.5319849,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("Graph");
   pt->Draw();
   p2->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->SetSelected(FinalPhiTot);
}
