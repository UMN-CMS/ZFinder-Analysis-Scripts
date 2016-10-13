{
//=========Macro generated from canvas: FinalPhiRatio/FinalPhiRatio
//=========  (Wed Oct  5 11:26:58 2016) by ROOT version5.34/18
   TCanvas *FinalPhiRatio = new TCanvas("FinalPhiRatio", "FinalPhiRatio",54,104,450,896);
   FinalPhiRatio->Range(-0.3003,0.40625,2.7027,1.34375);
   FinalPhiRatio->SetFillColor(0);
   FinalPhiRatio->SetBorderMode(0);
   FinalPhiRatio->SetBorderSize(2);
   FinalPhiRatio->SetFrameBorderMode(0);
   FinalPhiRatio->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("");

   Int_t ci;   // for color index setting
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
   
   TH1F *Graph_Graph_Graph58 = new TH1F("Graph_Graph_Graph58","",100,0,2.64);
   Graph_Graph_Graph58->SetMinimum(0.5);
   Graph_Graph_Graph58->SetMaximum(1.25);
   Graph_Graph_Graph58->SetDirectory(0);
   Graph_Graph_Graph58->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph58->SetLineColor(ci);
   Graph_Graph_Graph58->GetXaxis()->SetTitle("y_{ll}");
   Graph_Graph_Graph58->GetXaxis()->SetRange(1,91);
   Graph_Graph_Graph58->GetXaxis()->CenterTitle(true);
   Graph_Graph_Graph58->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph58->GetXaxis()->SetLabelOffset(0.012);
   Graph_Graph_Graph58->GetXaxis()->SetTitleSize(0.043);
   Graph_Graph_Graph58->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph58->GetYaxis()->SetTitle("MC/Data");
   Graph_Graph_Graph58->GetYaxis()->CenterTitle(true);
   Graph_Graph_Graph58->GetYaxis()->SetNdivisions(505);
   Graph_Graph_Graph58->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph58->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph_Graph58->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph58->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph58->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph58->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph58->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph58);
   
   grae->Draw("ae2");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#ff0000");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(21);
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
   
   TH1F *Graph_Graph_Graph79 = new TH1F("Graph_Graph_Graph79","Graph",100,0,2.4);
   Graph_Graph_Graph79->SetMinimum(0.9724994);
   Graph_Graph_Graph79->SetMaximum(1.002016);
   Graph_Graph_Graph79->SetDirectory(0);
   Graph_Graph_Graph79->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph79->SetLineColor(ci);
   Graph_Graph_Graph79->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph79->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph79->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph79->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph79->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph79->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph79->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph79->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph79->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph79->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph79->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph79->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph79);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#6666ff");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#6666ff");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(22);
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
   
   TH1F *Graph_Graph_Graph610 = new TH1F("Graph_Graph_Graph610","Graph",100,0,2.4);
   Graph_Graph_Graph610->SetMinimum(0.9203415);
   Graph_Graph_Graph610->SetMaximum(1.058963);
   Graph_Graph_Graph610->SetDirectory(0);
   Graph_Graph_Graph610->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph_Graph610->SetLineColor(ci);
   Graph_Graph_Graph610->GetXaxis()->SetLabelFont(42);
   Graph_Graph_Graph610->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph610->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph610->GetXaxis()->SetTitleFont(42);
   Graph_Graph_Graph610->GetYaxis()->SetLabelFont(42);
   Graph_Graph_Graph610->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph610->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph610->GetYaxis()->SetTitleFont(42);
   Graph_Graph_Graph610->GetZaxis()->SetLabelFont(42);
   Graph_Graph_Graph610->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph_Graph610->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph_Graph610->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph_Graph610);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#00cc00");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#00cc00");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(3);
   grae->SetPoint(0,0.2,1.002653);
   grae->SetPointError(0,0,0,2.000841e-08,2.000841e-08);
   grae->SetPoint(1,0.6,1.0057);
   grae->SetPointError(1,0,0,2.030835e-08,2.030835e-08);
   grae->SetPoint(2,1,1.000945);
   grae->SetPointError(2,0,0,2.084153e-08,2.084153e-08);
   grae->SetPoint(3,1.4,0.9457246);
   grae->SetPointError(3,0,0,2.441291e-08,2.441291e-08);
   grae->SetPoint(4,1.8,0.7882002);
   grae->SetPointError(4,0,0,3.728089e-08,3.728089e-08);
   grae->SetPoint(5,2.2,0.5088992);
   grae->SetPointError(5,0,0,1.554022e-07,1.554022e-07);
   
   TH1F *Graph_Graph11 = new TH1F("Graph_Graph11","Graph",100,0,2.4);
   Graph_Graph11->SetMinimum(0.4592189);
   Graph_Graph11->SetMaximum(1.055381);
   Graph_Graph11->SetDirectory(0);
   Graph_Graph11->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph11->SetLineColor(ci);
   Graph_Graph11->GetXaxis()->SetLabelFont(42);
   Graph_Graph11->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetXaxis()->SetTitleFont(42);
   Graph_Graph11->GetYaxis()->SetLabelFont(42);
   Graph_Graph11->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetYaxis()->SetTitleFont(42);
   Graph_Graph11->GetZaxis()->SetLabelFont(42);
   Graph_Graph11->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph11->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph11->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph11);
   
   grae->Draw("pe");
   
   grae = new TGraphAsymmErrors(6);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");
   grae->SetFillColor(1);

   ci = TColor::GetColor("#009999");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);

   ci = TColor::GetColor("#009999");
   grae->SetMarkerColor(ci);
   grae->SetMarkerStyle(24);
   grae->SetPoint(0,0.2,0.005150031);
   grae->SetPointError(0,0,0,2.000841e-08,2.000841e-08);
   grae->SetPoint(1,0.6,0.005184999);
   grae->SetPointError(1,0,0,2.030835e-08,2.030835e-08);
   grae->SetPoint(2,1,0.005190294);
   grae->SetPointError(2,0,0,2.084153e-08,2.084153e-08);
   grae->SetPoint(3,1.4,0.005199018);
   grae->SetPointError(3,0,0,2.441291e-08,2.441291e-08);
   grae->SetPoint(4,1.8,0.005359058);
   grae->SetPointError(4,0,0,3.728089e-08,3.728089e-08);
   grae->SetPoint(5,2.2,0.005400582);
   grae->SetPointError(5,0,0,1.554022e-07,1.554022e-07);
   
   TH1F *Graph_Graph12 = new TH1F("Graph_Graph12","Graph",100,0,2.4);
   Graph_Graph12->SetMinimum(0.005124938);
   Graph_Graph12->SetMaximum(0.00542581);
   Graph_Graph12->SetDirectory(0);
   Graph_Graph12->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph12->SetLineColor(ci);
   Graph_Graph12->GetXaxis()->SetLabelFont(42);
   Graph_Graph12->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetXaxis()->SetTitleFont(42);
   Graph_Graph12->GetYaxis()->SetLabelFont(42);
   Graph_Graph12->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetYaxis()->SetTitleFont(42);
   Graph_Graph12->GetZaxis()->SetLabelFont(42);
   Graph_Graph12->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph12->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph12);
   
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
   grae->SetMarkerStyle(25);
   grae->SetPoint(0,0.2,0.5419867);
   grae->SetPointError(0,0,0,2.000841e-08,2.000841e-08);
   grae->SetPoint(1,0.6,0.5447232);
   grae->SetPointError(1,0,0,2.030835e-08,2.030835e-08);
   grae->SetPoint(2,1,0.5424681);
   grae->SetPointError(2,0,0,2.084153e-08,2.084153e-08);
   grae->SetPoint(3,1.4,0.5427314);
   grae->SetPointError(3,0,0,2.441291e-08,2.441291e-08);
   grae->SetPoint(4,1.8,0.5514604);
   grae->SetPointError(4,0,0,3.728089e-08,3.728089e-08);
   grae->SetPoint(5,2.2,0.5479487);
   grae->SetPointError(5,0,0,1.554022e-07,1.554022e-07);
   
   TH1F *Graph_Graph13 = new TH1F("Graph_Graph13","Graph",100,0,2.4);
   Graph_Graph13->SetMinimum(0.5410393);
   Graph_Graph13->SetMaximum(0.5524078);
   Graph_Graph13->SetDirectory(0);
   Graph_Graph13->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph13->SetLineColor(ci);
   Graph_Graph13->GetXaxis()->SetLabelFont(42);
   Graph_Graph13->GetXaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetXaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetXaxis()->SetTitleFont(42);
   Graph_Graph13->GetYaxis()->SetLabelFont(42);
   Graph_Graph13->GetYaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetYaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetYaxis()->SetTitleFont(42);
   Graph_Graph13->GetZaxis()->SetLabelFont(42);
   Graph_Graph13->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph13->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph13->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph13);
   
   grae->Draw("pe");
   TLatex *   tex = new TLatex(0.7,0.907,"19.7 fb^{-1} (8 TeV)");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.907,"CMS Preliminary");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.25,"|#eta^{l_{0}}| < 2.1,        |#eta^{l_{1}}| < 2.4");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.2,"p_{T}^{l_{0}} > 30 GeV,   p_{T}^{l_{1}} > 20 GeV");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.15,0.15,"60 GeV < M_{ll} < 120 GeV");
tex->SetNDC();
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->Draw();
   
   TLegend *leg = new TLegend(0.13,0.72,0.7,0.9,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(22);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(19);
   leg->SetFillStyle(0);
   TLegendEntry *entry=leg->AddEntry("Graph0","2012 data","F");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph1"," MadGraph+Pythia6 (Z2*)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#6666ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(22);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph2","POWHEG+Pythia6 (Z2*)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph0","Resbos","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#00cc00");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(3);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph1","AMC@nlo+Pythia8(CUETP8M1)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#009999");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(24);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   entry=leg->AddEntry("Graph2","POWHEG+Pythia8 (CT10)","P");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);

   ci = TColor::GetColor("#ff0000");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(25);
   entry->SetMarkerSize(1);
   entry->SetTextFont(22);
   leg->Draw();
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","",100,0,2.64);
   Graph_Graph5->SetMinimum(0.5);
   Graph_Graph5->SetMaximum(1.25);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);

   ci = TColor::GetColor("#000099");
   Graph_Graph5->SetLineColor(ci);
   Graph_Graph5->GetXaxis()->SetTitle("y_{ll}");
   Graph_Graph5->GetXaxis()->SetRange(1,91);
   Graph_Graph5->GetXaxis()->CenterTitle(true);
   Graph_Graph5->GetXaxis()->SetLabelFont(42);
   Graph_Graph5->GetXaxis()->SetLabelOffset(0.012);
   Graph_Graph5->GetXaxis()->SetTitleSize(0.043);
   Graph_Graph5->GetXaxis()->SetTitleFont(42);
   Graph_Graph5->GetYaxis()->SetTitle("MC/Data");
   Graph_Graph5->GetYaxis()->CenterTitle(true);
   Graph_Graph5->GetYaxis()->SetNdivisions(505);
   Graph_Graph5->GetYaxis()->SetLabelFont(42);
   Graph_Graph5->GetYaxis()->SetTitleOffset(1.2);
   Graph_Graph5->GetYaxis()->SetTitleFont(42);
   Graph_Graph5->GetZaxis()->SetLabelFont(42);
   Graph_Graph5->GetZaxis()->SetLabelSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleSize(0.035);
   Graph_Graph5->GetZaxis()->SetTitleFont(42);
   Graph_Graph5->Draw("sameaxis");
   FinalPhiRatio->Modified();
   FinalPhiRatio->cd();
   FinalPhiRatio->SetSelected(FinalPhiRatio);
}
