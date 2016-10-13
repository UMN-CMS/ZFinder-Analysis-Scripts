{
//=========Macro generated from canvas: FinalPhiTot/FinalPhiTot
//=========  (Tue Jul 19 10:29:54 2016) by ROOT version5.34/18
   TCanvas *FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot",1684,48,625,898);
   FinalPhiTot->Range(0,0,1,1);
   FinalPhiTot->SetFillColor(0);
   FinalPhiTot->SetBorderMode(0);
   FinalPhiTot->SetBorderSize(2);
   FinalPhiTot->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: p1
   TPad *p1 = new TPad("p1", "p1",0,0.2777778,1,1);
   p1->Draw();
   p1->cd();
   p1->Range(-0.4556962,-0.00860215,2.582278,0.8516129);
   p1->SetFillColor(0);
   p1->SetBorderMode(0);
   p1->SetBorderSize(0);
   p1->SetLeftMargin(0.15);
   p1->SetRightMargin(0.06);
   p1->SetTopMargin(0.06);
   p1->SetBottomMargin(0.01);
   p1->SetFrameBorderMode(0);
   p1->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.2,0.5642977);
   grae->SetPointError(0,0,0,0.0008302697,0.0008302697);
   grae->SetPoint(1,0.6,0.5559633);
   grae->SetPointError(1,0,0,0.0008096385,0.0008096385);
   grae->SetPoint(2,1,0.5417405);
   grae->SetPointError(2,0,0,0.0008599714,0.0008599714);
   grae->SetPoint(3,1.4,0.4624889);
   grae->SetPointError(3,0,0,0.0008575522,0.0008575522);
   grae->SetPoint(4,1.8,0.3028549);
   grae->SetPointError(4,0,0,0.0008078161,0.0008078161);
   grae->SetPoint(5,2.2,0.07265468);
   grae->SetPointError(5,0,0,0.0004232915,0.0004232915);
   
   TH1F *Graph_Graph428 = new TH1F("Graph_Graph428","",100,0,2.4);
   Graph_Graph428->SetMinimum(0);
   Graph_Graph428->SetMaximum(0.8);
   Graph_Graph428->SetDirectory(0);
   Graph_Graph428->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph428->SetLineColor(ci);
   Graph_Graph428->GetXaxis()->SetRange(1,100);
   Graph_Graph428->GetXaxis()->SetLabelFont(42);
   Graph_Graph428->GetXaxis()->SetLabelSize(0);
   Graph_Graph428->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph428->GetXaxis()->SetTitleFont(42);
   Graph_Graph428->GetYaxis()->SetLabelFont(42);
   Graph_Graph428->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph428->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph428->GetYaxis()->SetTitleFont(42);
   Graph_Graph428->GetZaxis()->SetLabelFont(42);
   Graph_Graph428->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph428->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph428->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph428);
   
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
   grae->SetPoint(0,0.2,0.5489086);
   grae->SetPointError(0,0,0,0.0006210304,0.0006210304);
   grae->SetPoint(1,0.6,0.5471653);
   grae->SetPointError(1,0,0,0.0006200069,0.0006200069);
   grae->SetPoint(2,1,0.537466);
   grae->SetPointError(2,0,0,0.0006142905,0.0006142905);
   grae->SetPoint(3,1.4,0.4688073);
   grae->SetPointError(3,0,0,0.0005740266,0.0005740266);
   grae->SetPoint(4,1.8,0.3192421);
   grae->SetPointError(4,0,0,0.0004742752,0.0004742752);
   grae->SetPoint(5,2.2,0.07841084);
   grae->SetPointError(5,0,0,0.0002353245,0.0002353245);
   
   TH1F *Graph_Graph429 = new TH1F("Graph_Graph429","Graph",100,0,2.4);
   Graph_Graph429->SetMinimum(0.03104011);
   Graph_Graph429->SetMaximum(0.596665);
   Graph_Graph429->SetDirectory(0);
   Graph_Graph429->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph429->SetLineColor(ci);
   Graph_Graph429->GetXaxis()->SetLabelFont(42);
   Graph_Graph429->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph429->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph429->GetXaxis()->SetTitleFont(42);
   Graph_Graph429->GetYaxis()->SetLabelFont(42);
   Graph_Graph429->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph429->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph429->GetYaxis()->SetTitleFont(42);
   Graph_Graph429->GetZaxis()->SetLabelFont(42);
   Graph_Graph429->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph429->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph429->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph429);
   
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
   grae->SetPoint(0,0.2,0.5627517);
   grae->SetPointError(0,0,0,0.001427191,0.001427191);
   grae->SetPoint(1,0.6,0.5569279);
   grae->SetPointError(1,0,0,0.0009846963,0.0009846963);
   grae->SetPoint(2,1,0.5403595);
   grae->SetPointError(2,0,0,0.0004463375,0.0004463375);
   grae->SetPoint(3,1.4,0.4603819);
   grae->SetPointError(3,0,0,0.0006118243,0.0006118243);
   grae->SetPoint(4,1.8,0.3064496);
   grae->SetPointError(4,0,0,0.0007501454,0.0007501454);
   grae->SetPoint(5,2.2,0.07312933);
   grae->SetPointError(5,0,0,0.0002714195,0.0002714195);
   
   TH1F *Graph_Graph430 = new TH1F("Graph_Graph430","Graph",100,0,2.4);
   Graph_Graph430->SetMinimum(0.02372581);
   Graph_Graph430->SetMaximum(0.613311);
   Graph_Graph430->SetDirectory(0);
   Graph_Graph430->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph430->SetLineColor(ci);
   Graph_Graph430->GetXaxis()->SetLabelFont(42);
   Graph_Graph430->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph430->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph430->GetXaxis()->SetTitleFont(42);
   Graph_Graph430->GetYaxis()->SetLabelFont(42);
   Graph_Graph430->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph430->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph430->GetYaxis()->SetTitleFont(42);
   Graph_Graph430->GetZaxis()->SetLabelFont(42);
   Graph_Graph430->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph430->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph430->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph430);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("");

   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.2,0.5642977);
   grae->SetPointError(0,0,0,0.0008302697,0.0008302697);
   grae->SetPoint(1,0.6,0.5559633);
   grae->SetPointError(1,0,0,0.0008096385,0.0008096385);
   grae->SetPoint(2,1,0.5417405);
   grae->SetPointError(2,0,0,0.0008599714,0.0008599714);
   grae->SetPoint(3,1.4,0.4624889);
   grae->SetPointError(3,0,0,0.0008575522,0.0008575522);
   grae->SetPoint(4,1.8,0.3028549);
   grae->SetPointError(4,0,0,0.0008078161,0.0008078161);
   grae->SetPoint(5,2.2,0.07265468);
   grae->SetPointError(5,0,0,0.0004232915,0.0004232915);
   
   TH1F *Graph_Graph_Graph428431 = new TH1F("Graph_Graph_Graph428431","",100,0,2.4);
   Graph_Graph_Graph428431->SetMinimum(0);
   Graph_Graph_Graph428431->SetMaximum(0.8);
   Graph_Graph_Graph428431->SetDirectory(0);
   Graph_Graph_Graph428431->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph428431->SetLineColor(ci);
   Graph_Graph_Graph428431->GetXaxis()->SetRange(1,100);
   Graph_Graph_Graph428431->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph428431->GetXaxis()->SetLabelSize(0);
   Graph_Graph_Graph428431->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph428431->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph428431->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph428431->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph428431->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph428431->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph428431->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph428431->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph428431->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph428431->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph428431);
   
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
   p1->Modified();
   FinalPhiTot->cd();
  
// ------------>Primitives in pad: p2
   p2 = new TPad("p2", "p2",0,0,1,0.2777778);
   p2->Draw();
   p2->cd();
   p2->Range(-0.4561519,0.7806452,2.584861,1.103226);
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
   grae->SetTitle("");

   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetPoint(0,0.2,1);
   grae->SetPointError(0,0.2,0.2,0.001471333,0.001471333);
   grae->SetPoint(1,0.6,1);
   grae->SetPointError(1,0.2,0.2,0.00145628,0.00145628);
   grae->SetPoint(2,1,1);
   grae->SetPointError(2,0.2,0.2,0.001587423,0.001587423);
   grae->SetPoint(3,1.4,1);
   grae->SetPointError(3,0.2,0.2,0.001854211,0.001854211);
   grae->SetPoint(4,1.8,1);
   grae->SetPointError(4,0.2,0.2,0.002667337,0.002667337);
   grae->SetPoint(5,2.2,1);
   grae->SetPointError(5,0.2,0.2,0.005826074,0.005826074);
   
   TH1F *Graph_Graph432 = new TH1F("Graph_Graph432","",100,0,2.64);
   Graph_Graph432->SetMinimum(0.9);
   Graph_Graph432->SetMaximum(1.1);
   Graph_Graph432->SetDirectory(0);
   Graph_Graph432->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph432->SetLineColor(ci);
   Graph_Graph432->GetXaxis()->SetRange(1,91);
   Graph_Graph432->GetXaxis()->SetLabelFont(42);
   Graph_Graph432->GetXaxis()->SetLabelOffset(0.2);
   Graph_Graph432->GetXaxis()->SetLabelSize(0.1);
   Graph_Graph432->GetXaxis()->SetTitleFont(42);
   Graph_Graph432->GetYaxis()->SetTitle("MC/Data     ");
   Graph_Graph432->GetYaxis()->SetNdivisions(505);
   Graph_Graph432->GetYaxis()->SetLabelFont(42);
   Graph_Graph432->GetYaxis()->SetTitleSize(0.07);
   Graph_Graph432->GetYaxis()->SetTitleOffset(0.5);
   Graph_Graph432->GetYaxis()->SetTitleFont(42);
   Graph_Graph432->GetZaxis()->SetLabelFont(42);
   Graph_Graph432->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph432->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph432->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph432);
   
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
   grae->SetPoint(0,0.2,0.9727288);
   grae->SetPointError(0,0,0,0.001100537,0.001100537);
   grae->SetPoint(1,0.6,0.9841751);
   grae->SetPointError(1,0,0,0.001115194,0.001115194);
   grae->SetPoint(2,1,0.9921097);
   grae->SetPointError(2,0,0,0.00113392,0.00113392);
   grae->SetPoint(3,1.4,1.013662);
   grae->SetPointError(3,0,0,0.001241168,0.001241168);
   grae->SetPoint(4,1.8,1.054109);
   grae->SetPointError(4,0,0,0.001566015,0.001566015);
   grae->SetPoint(5,2.2,1.079226);
   grae->SetPointError(5,0,0,0.003238944,0.003238944);
   
   TH1F *Graph_Graph433 = new TH1F("Graph_Graph433","Graph",100,0,2.4);
   Graph_Graph433->SetMinimum(0.9605445);
   Graph_Graph433->SetMaximum(1.093549);
   Graph_Graph433->SetDirectory(0);
   Graph_Graph433->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph433->SetLineColor(ci);
   Graph_Graph433->GetXaxis()->SetLabelFont(42);
   Graph_Graph433->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph433->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph433->GetXaxis()->SetTitleFont(42);
   Graph_Graph433->GetYaxis()->SetLabelFont(42);
   Graph_Graph433->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph433->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph433->GetYaxis()->SetTitleFont(42);
   Graph_Graph433->GetZaxis()->SetLabelFont(42);
   Graph_Graph433->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph433->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph433->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph433);
   
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
   grae->SetPoint(0,0.2,0.9972603);
   grae->SetPointError(0,0,0,0.002529146,0.002529146);
   grae->SetPoint(1,0.6,1.001735);
   grae->SetPointError(1,0,0,0.001771153,0.001771153);
   grae->SetPoint(2,1,0.9974508);
   grae->SetPointError(2,0,0,0.0008238954,0.0008238954);
   grae->SetPoint(3,1.4,0.9954442);
   grae->SetPointError(3,0,0,0.001322895,0.001322895);
   grae->SetPoint(4,1.8,1.01187);
   grae->SetPointError(4,0,0,0.002476914,0.002476914);
   grae->SetPoint(5,2.2,1.006533);
   grae->SetPointError(5,0,0,0.003735747,0.003735747);
   
   TH1F *Graph_Graph434 = new TH1F("Graph_Graph434","Graph",100,0,2.4);
   Graph_Graph434->SetMinimum(0.9920988);
   Graph_Graph434->SetMaximum(1.016369);
   Graph_Graph434->SetDirectory(0);
   Graph_Graph434->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph434->SetLineColor(ci);
   Graph_Graph434->GetXaxis()->SetLabelFont(42);
   Graph_Graph434->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph434->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph434->GetXaxis()->SetTitleFont(42);
   Graph_Graph434->GetYaxis()->SetLabelFont(42);
   Graph_Graph434->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph434->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph434->GetYaxis()->SetTitleFont(42);
   Graph_Graph434->GetZaxis()->SetLabelFont(42);
   Graph_Graph434->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph434->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph434->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph434);
   
   grae->Draw("pe");
   p2->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->SetSelected(FinalPhiTot);
}
