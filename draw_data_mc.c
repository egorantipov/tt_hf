///////////////////////////////
// Draw Data/Bkgd Comparison //
///////////////////////////////
int draw_data_mc_plot(TH1 *h_data, TH1 *h_mc, TString title, TString savename)
{
  TCanvas *c = new TCanvas(title, title, 1600, 1200);
  
  gStyle->SetOptStat(0);

  h_data->SetMarkerStyle(20);
  h_data->SetMarkerColor(1);
  h_mc->SetLineColor(4);
  
  TPad *tPad = new TPad("tPad", "tPad", 0.005, 0.4, 0.995, 0.995);
  TPad *bPad = new TPad("bPad", "bPad", 0.005, 0.1, 0.995, 0.4);
  tPad->Draw();
  bPad->Draw("same");

  tPad->cd();
  h_data->Draw("P");
  h_data->GetXaxis()->SetLabelSize(0);
  h_data->GetYaxis()->SetTitle("events");
  h_mc->Draw("Lsame");
  
  bPad->cd();
  TH1 *h_ratio = (TH1*)h_data->Clone();
  h_ratio->Divide(h_mc);
  h_ratio->Draw("E1");
  
  c->Print("Plots/data_mc_comparison/"+savename+".png");
  return 0;
}


//////////////
//   MAIN   //
//////////////
void draw_data_mc()
{
  TFile *mc_hists_file = new TFile("hists_mc.root");
  TH1 *h_met_mc = (TH1*)mc_hists_file->Get("2b_emu_OS_met");
  
  TFile *data_hists_file = new TFile("hists_data.root");
  TH1 *h_met_data = (TH1*)data_hists_file->Get("2b_emu_OS_met");

  int draw_met = draw_data_mc_plot(h_met_data, h_met_mc, "met", "2b_emu_OS_met");
}
