{
//=========Macro generated from canvas: CPull/CPull
//=========  (Wed Apr  1 12:38:34 2015) by ROOT version5.34/18
   TCanvas *CPull = new TCanvas("CPull", "CPull",52,27,800,900);
   CPull->Range(-6.25,-1.875,6.25,16.875);
   CPull->SetFillColor(0);
   CPull->SetBorderMode(0);
   CPull->SetBorderSize(2);
   CPull->SetFrameBorderMode(0);
   CPull->SetFrameBorderMode(0);
   
   TH1D *Pull = new TH1D("Pull","Pull",20,-5,5);
   Pull->SetBinContent(7,1);
   Pull->SetBinContent(8,2);
   Pull->SetBinContent(9,7);
   Pull->SetBinContent(10,10);
   Pull->SetBinContent(11,4);
   Pull->SetBinContent(12,6);
   Pull->SetBinContent(13,5);
   Pull->SetBinError(7,1);
   Pull->SetBinError(8,1.414214);
   Pull->SetBinError(9,2.645751);
   Pull->SetBinError(10,3.162278);
   Pull->SetBinError(11,2);
   Pull->SetBinError(12,2.44949);
   Pull->SetBinError(13,2.236068);
   Pull->SetMinimum(0);
   Pull->SetMaximum(15);
   Pull->SetEntries(35);
   Pull->SetStats(0);
   Pull->GetXaxis()->SetTitle("(Unfolded/Generated-1)/#sigma(Unfolded/Generated)");
   Pull->GetXaxis()->SetLabelFont(42);
   Pull->GetXaxis()->SetLabelOffset(0);
   Pull->GetXaxis()->SetTitleFont(42);
   Pull->GetYaxis()->SetTitle("#phi* bins");
   Pull->GetYaxis()->SetLabelFont(42);
   Pull->GetYaxis()->SetTitleFont(42);
   Pull->GetZaxis()->SetLabelFont(42);
   Pull->GetZaxis()->SetLabelSize(0.035);
   Pull->GetZaxis()->SetTitleSize(0.035);
   Pull->GetZaxis()->SetTitleFont(42);
   Pull->Draw("");
   TLatex *   tex = new TLatex(0.19,0.84,"MadGraph");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
      tex = new TLatex(0.19,0.8,"unfolded using 50000 Powheg events");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   CPull->Modified();
   CPull->cd();
   CPull->SetSelected(CPull);
}
