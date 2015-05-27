{
//=========Macro generated from canvas: CPull/CPull
//=========  (Wed Apr  1 12:46:12 2015) by ROOT version5.34/18
   TCanvas *CPull = new TCanvas("CPull", "CPull",8,131,800,900);
   CPull->Range(-6.25,-2.5,6.25,22.5);
   CPull->SetFillColor(0);
   CPull->SetBorderMode(0);
   CPull->SetBorderSize(2);
   CPull->SetFrameBorderMode(0);
   CPull->SetFrameBorderMode(0);
   
   TH1D *Pull = new TH1D("Pull","Pull",20,-5,5);
   Pull->SetBinContent(8,2);
   Pull->SetBinContent(9,3);
   Pull->SetBinContent(10,12);
   Pull->SetBinContent(11,11);
   Pull->SetBinContent(12,5);
   Pull->SetBinContent(13,1);
   Pull->SetBinContent(14,1);
   Pull->SetBinError(8,1.414214);
   Pull->SetBinError(9,1.732051);
   Pull->SetBinError(10,3.464102);
   Pull->SetBinError(11,3.316625);
   Pull->SetBinError(12,2.236068);
   Pull->SetBinError(13,1);
   Pull->SetBinError(14,1);
   Pull->SetMinimum(0);
   Pull->SetMaximum(20);
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
      tex = new TLatex(0.19,0.8,"unfolded using full Powheg sample");
tex->SetNDC();
   tex->SetTextFont(42);
   tex->SetTextSize(0.04);
   tex->SetLineWidth(2);
   tex->Draw();
   CPull->Modified();
   CPull->cd();
   CPull->SetSelected(CPull);
}
