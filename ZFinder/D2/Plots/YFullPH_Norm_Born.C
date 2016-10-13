{
//=========Macro generated from canvas: FinalPhiTot/FinalPhiTot
//=========  (Wed Oct  5 11:18:57 2016) by ROOT version5.34/18
   TCanvas *FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot",6,128,450,896);
   FinalPhiTot->Range(0,0,1,1);
   FinalPhiTot->SetFillColor(0);
   FinalPhiTot->SetBorderMode(0);
   FinalPhiTot->SetBorderSize(2);
   FinalPhiTot->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: p1
   TPad *p1 = new TPad("p1", "p1",0,0.2777778,1,1);
   p1->Draw();
   p1->cd();
   p1->Range(-0.4556962,0.001505377,2.582278,0.8509677);
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
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",100,0,2.4);
   Graph_Graph1->SetMinimum(0.01);
   Graph_Graph1->SetMaximum(0.8);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph1->SetLineColor(ci);
   Graph_Graph1->GetXaxis()->SetRange(1,100);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelSize(0);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/dy");
   Graph_Graph1->GetYaxis()->CenterTitle(true);
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1);
   
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
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,2.4);
   Graph_Graph2->SetMinimum(0.03104011);
   Graph_Graph2->SetMaximum(0.596665);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph2->SetLineColor(ci);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph2);
   
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
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,0,2.4);
   Graph_Graph3->SetMinimum(0.02372581);
   Graph_Graph3->SetMaximum(0.613311);
   Graph_Graph3->SetDirectory(0);
   Graph_Graph3->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph3->SetLineColor(ci);
   Graph_Graph3->GetXaxis()->SetLabelFont(42);
   Graph_Graph3->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetXaxis()->SetTitleFont(42);
   Graph_Graph3->GetYaxis()->SetLabelFont(42);
   Graph_Graph3->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetYaxis()->SetTitleFont(42);
   Graph_Graph3->GetZaxis()->SetLabelFont(42);
   Graph_Graph3->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph3->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph3);
   
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
   
   TH1F *Graph_Graph_Graph14 = new TH1F("Graph_Graph_Graph14","",100,0,2.4);
   Graph_Graph_Graph14->SetMinimum(0.01);
   Graph_Graph_Graph14->SetMaximum(0.8);
   Graph_Graph_Graph14->SetDirectory(0);
   Graph_Graph_Graph14->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph14->SetLineColor(ci);
   Graph_Graph_Graph14->GetXaxis()->SetRange(1,100);
   Graph_Graph_Graph14->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph14->GetXaxis()->SetLabelSize(0);
   Graph_Graph_Graph14->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph14->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph14->GetYaxis()->SetTitle("1/#sigma^{fid} #bullet d#sigma^{fid}/dy");
   Graph_Graph_Graph14->GetYaxis()->CenterTitle(true);
   Graph_Graph_Graph14->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph14->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph_Graph14->GetYaxis()->SetTitleSize(0.065);
   Graph_Graph_Graph14->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph14->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph14->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph14->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph14->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph14);
   
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
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","",100,0,2.64);
   Graph_Graph5->SetMinimum(0.9);
   Graph_Graph5->SetMaximum(1.1);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph5->SetLineColor(ci);
   Graph_Graph5->GetXaxis()->SetTitle("y_{ll}");
   Graph_Graph5->GetXaxis()->SetRange(1,91);
   Graph_Graph5->GetXaxis()->CenterTitle(true);
   Graph_Graph5->GetXaxis()->SetLabelFont(42);
   Graph_Graph5->GetXaxis()->SetLabelOffset(0.02);
   Graph_Graph5->GetXaxis()->SetLabelSize(0.13);
   Graph_Graph5->GetXaxis()->SetTitleSize(0.15);
   Graph_Graph5->GetXaxis()->SetTitleFont(42);
   Graph_Graph5->GetYaxis()->SetTitle("MC/Data ");
   Graph_Graph5->GetYaxis()->SetNdivisions(505);
   Graph_Graph5->GetYaxis()->SetLabelFont(42);
   Graph_Graph5->GetYaxis()->SetLabelSize(0.13);
   Graph_Graph5->GetYaxis()->SetTitleSize(0.15);
   Graph_Graph5->GetYaxis()->SetTitleOffset(0.4);
   Graph_Graph5->GetYaxis()->SetTitleFont(42);
   Graph_Graph5->GetZaxis()->SetLabelFont(42);
   Graph_Graph5->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph5);
   
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
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,0,2.4);
   Graph_Graph6->SetMinimum(0.9605445);
   Graph_Graph6->SetMaximum(1.093549);
   Graph_Graph6->SetDirectory(0);
   Graph_Graph6->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph6->SetLineColor(ci);
   Graph_Graph6->GetXaxis()->SetLabelFont(42);
   Graph_Graph6->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetXaxis()->SetTitleFont(42);
   Graph_Graph6->GetYaxis()->SetLabelFont(42);
   Graph_Graph6->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetYaxis()->SetTitleFont(42);
   Graph_Graph6->GetZaxis()->SetLabelFont(42);
   Graph_Graph6->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph6->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph6->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph6);
   
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
   
   TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","Graph",100,0,2.4);
   Graph_Graph7->SetMinimum(0.9920988);
   Graph_Graph7->SetMaximum(1.016369);
   Graph_Graph7->SetDirectory(0);
   Graph_Graph7->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph7->SetLineColor(ci);
   Graph_Graph7->GetXaxis()->SetLabelFont(42);
   Graph_Graph7->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetXaxis()->SetTitleFont(42);
   Graph_Graph7->GetYaxis()->SetLabelFont(42);
   Graph_Graph7->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetYaxis()->SetTitleFont(42);
   Graph_Graph7->GetZaxis()->SetLabelFont(42);
   Graph_Graph7->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph7->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph7);
   
   grae->Draw("pe");
   p2->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->Modified();
   FinalPhiTot->cd();
   FinalPhiTot->SetSelected(FinalPhiTot);
}
