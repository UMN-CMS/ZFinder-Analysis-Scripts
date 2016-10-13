{
//=========Macro generated from canvas: FinalPhiTot/FinalPhiTot
//=========  (Wed Oct  5 11:26:58 2016) by ROOT version5.34/18
   TCanvas *FinalPhiTot = new TCanvas("FinalPhiTot", "FinalPhiTot",52,80,450,896);
   FinalPhiTot->Range(0,0,1,1);
   FinalPhiTot->SetFillColor(0);
   FinalPhiTot->SetBorderMode(0);
   FinalPhiTot->SetBorderSize(2);
   FinalPhiTot->SetFrameBorderMode(0);
  
// ------------>Primitives in pad: p1
   TPad *p1 = new TPad("p1", "p1",0,0.2777778,1,1);
   p1->Draw();
   p1->cd();
   p1->Range(-0.4556962,-3290.322,2.582278,425741.9);
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
   grae->SetPoint(0,0.2,270.7929);
   grae->SetPointError(0,0,0,1.391488,1.391488);
   grae->SetPoint(1,0.6,266.7934);
   grae->SetPointError(1,0,0,1.362168,1.362168);
   grae->SetPoint(2,1,259.9682);
   grae->SetPointError(2,0,0,1.332492,1.332492);
   grae->SetPoint(3,1.4,221.9373);
   grae->SetPointError(3,0,0,1.16218,1.16218);
   grae->SetPoint(4,1.8,145.3328);
   grae->SetPointError(4,0,0,0.8118371,0.8118371);
   grae->SetPoint(5,2.2,34.86524);
   grae->SetPointError(5,0,0,0.2687413,0.2687413);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","",100,0,2.4);
   Graph_Graph1->SetMinimum(1000);
   Graph_Graph1->SetMaximum(400000);
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
   grae->SetPoint(0,0.2,253.9099);
   grae->SetPointError(0,0,0,1.559823,1.559823);
   grae->SetPoint(1,0.6,253.1035);
   grae->SetPointError(1,0,0,1.555536,1.555536);
   grae->SetPoint(2,1,248.6169);
   grae->SetPointError(2,0,0,1.531514,1.531514);
   grae->SetPoint(3,1.4,216.8573);
   grae->SetPointError(3,0,0,1.34455,1.34455);
   grae->SetPoint(4,1.8,147.6725);
   grae->SetPointError(4,0,0,0.9275748,0.9275748);
   grae->SetPoint(5,2.2,36.27068);
   grae->SetPointError(5,0,0,0.2475475,0.2475475);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,2.4);
   Graph_Graph2->SetMinimum(14.07848);
   Graph_Graph2->SetMaximum(277.4144);
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
   grae->SetPoint(0,0.2,265.7569);
   grae->SetPointError(0,0,0,1.744915,1.744915);
   grae->SetPoint(1,0.6,263.0067);
   grae->SetPointError(1,0,0,1.573809,1.573809);
   grae->SetPoint(2,1,255.1823);
   grae->SetPointError(2,0,0,1.287629,1.287629);
   grae->SetPoint(3,1.4,217.4133);
   grae->SetPointError(3,0,0,0.9129161,0.9129161);
   grae->SetPoint(4,1.8,144.7194);
   grae->SetPointError(4,0,0,0.5488226,0.5488226);
   grae->SetPoint(5,2.2,34.53499);
   grae->SetPointError(5,0,0,0.1452997,0.1452997);
   
   TH1F *Graph_Graph3 = new TH1F("Graph_Graph3","Graph",100,0,2.4);
   Graph_Graph3->SetMinimum(11.07848);
   Graph_Graph3->SetMaximum(290.8131);
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
   grae->SetPoint(0,0.2,270.7929);
   grae->SetPointError(0,0,0,1.391488,1.391488);
   grae->SetPoint(1,0.6,266.7934);
   grae->SetPointError(1,0,0,1.362168,1.362168);
   grae->SetPoint(2,1,259.9682);
   grae->SetPointError(2,0,0,1.332492,1.332492);
   grae->SetPoint(3,1.4,221.9373);
   grae->SetPointError(3,0,0,1.16218,1.16218);
   grae->SetPoint(4,1.8,145.3328);
   grae->SetPointError(4,0,0,0.8118371,0.8118371);
   grae->SetPoint(5,2.2,34.86524);
   grae->SetPointError(5,0,0,0.2687413,0.2687413);
   
   TH1F *Graph_Graph_Graph14 = new TH1F("Graph_Graph_Graph14","",100,0,2.4);
   Graph_Graph_Graph14->SetMinimum(1000);
   Graph_Graph_Graph14->SetMaximum(400000);
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
   grae->SetPointError(0,0.2,0.2,0.005138568,0.005138568);
   grae->SetPoint(1,0.6,1);
   grae->SetPointError(1,0.2,0.2,0.005105705,0.005105705);
   grae->SetPoint(2,1,1);
   grae->SetPointError(2,0.2,0.2,0.005125594,0.005125594);
   grae->SetPoint(3,1.4,1);
   grae->SetPointError(3,0.2,0.2,0.005236523,0.005236523);
   grae->SetPoint(4,1.8,1);
   grae->SetPointError(4,0.2,0.2,0.005586056,0.005586056);
   grae->SetPoint(5,2.2,1);
   grae->SetPointError(5,0.2,0.2,0.007708001,0.007708001);
   
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
   grae->SetPoint(0,0.2,0.9376534);
   grae->SetPointError(0,0,0,0.005760207,0.005760207);
   grae->SetPoint(1,0.6,0.948687);
   grae->SetPointError(1,0,0,0.00583049,0.00583049);
   grae->SetPoint(2,1,0.9563355);
   grae->SetPointError(2,0,0,0.005891158,0.005891158);
   grae->SetPoint(3,1.4,0.9771103);
   grae->SetPointError(3,0,0,0.006058244,0.006058244);
   grae->SetPoint(4,1.8,1.016099);
   grae->SetPointError(4,0,0,0.00638242,0.00638242);
   grae->SetPoint(5,2.2,1.040311);
   grae->SetPointError(5,0,0,0.007100123,0.007100123);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,0,2.4);
   Graph_Graph6->SetMinimum(0.9203415);
   Graph_Graph6->SetMaximum(1.058963);
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
   grae->SetPoint(0,0.2,0.9814029);
   grae->SetPointError(0,0,0,0.006443725,0.006443725);
   grae->SetPoint(1,0.6,0.9858064);
   grae->SetPointError(1,0,0,0.005898979,0.005898979);
   grae->SetPoint(2,1,0.9815903);
   grae->SetPointError(2,0,0,0.004953025,0.004953025);
   grae->SetPoint(3,1.4,0.9796156);
   grae->SetPointError(3,0,0,0.004113396,0.004113396);
   grae->SetPoint(4,1.8,0.9957798);
   grae->SetPointError(4,0,0,0.003776317,0.003776317);
   grae->SetPoint(5,2.2,0.990528);
   grae->SetPointError(5,0,0,0.004167468,0.004167468);
   
   TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","Graph",100,0,2.4);
   Graph_Graph7->SetMinimum(0.9724994);
   Graph_Graph7->SetMaximum(1.002016);
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
