///////////////////////////////
// Draw Data/Bkgd Comparison //
///////////////////////////////
int draw_data_mc_plot(TH1 *h_data, TH1 *h_mc, TString title, TString savename)
{
  cout << "Drawing " << title << endl;

  TCanvas *c = new TCanvas("c", "c", 1600, 1200);
  
  gStyle->SetOptStat(0);

  h_data->SetMarkerStyle(20);
  h_data->SetMarkerColor(1);
  h_mc->SetLineColor(4);
  
  TPad *tPad = new TPad("tPad", "tPad", 0.005, 0.4, 0.995, 0.995);
  TPad *bPad = new TPad("bPad", "bPad", 0.005, 0.1, 0.995, 0.4);
  tPad->Draw();
  bPad->Draw("same");


  // Top pad: hists
  tPad->cd();
  tPad->SetGrid();
  tPad->SetRightMargin(0.05);
  tPad->SetLeftMargin(0.07);
  tPad->SetBottomMargin(0);
  tPad->SetTopMargin(0.03);
  gPad->SetLogy();
  
  h_data->Draw("E1P");
  h_data->SetMarkerSize(1.5);
  h_data->SetTitle("");
  h_data->GetXaxis()->SetLabelSize(0);
  h_data->GetYaxis()->SetTitle("#bf{N Events}");
  h_mc->Draw("e2 hist same");
  h_mc->SetLineWidth(4);
  h_mc->SetLineColor(2);
  

  // Bottom pad: ratio
  bPad->cd();
  bPad->SetRightMargin(0.05);
  bPad->SetLeftMargin(0.07);
  bPad->SetTopMargin(0);
  bPad->SetBottomMargin(0.4);
  bPad->SetGrid();
  TH1 *h_ratio = (TH1*)h_data->Clone();
  
  h_ratio->Divide(h_mc);
  h_ratio->SetTitle("");
  
  h_ratio->GetXaxis()->SetLabelSize(0.10);
  h_ratio->GetXaxis()->SetTitle(title);
  h_ratio->GetXaxis()->SetTitleOffset(1.3);
  h_ratio->GetXaxis()->SetTitleSize(0.12);
  
  h_ratio->GetYaxis()->SetLabelSize(0.05);  
  h_ratio->GetYaxis()->SetTitle("#bf{Data/mc}");
  h_ratio->GetYaxis()->SetTitleOffset(0.5);
  h_ratio->GetYaxis()->SetTitleSize(0.07);
  h_ratio->GetYaxis()->CenterTitle();
  h_ratio->GetYaxis()->SetRangeUser(0,2);
  h_ratio->GetYaxis()->SetNdivisions(8);
  h_ratio->Draw("E1");
  
  
  // Save the plot and print yileds
  c->Print("Plots/data_mc_comparison/"+savename+".png");

  double data_int = h_data->Integral(0, h_data->GetNbinsX());
  double mc_int = h_mc->Integral(0, h_mc->GetNbinsX());
  cout << "Data integral = " << data_int << "\nmc integral =   " << mc_int << "\n\n" << endl;
  
  return 0;
}


//////////////
//   MAIN   //
//////////////
void draw_data_mc()
{
  // Get mc hists
  TFile *mc_hists_file = new TFile("hists_mc.root");
  TH1 *h_mc_jets_n = (TH1*)mc_hists_file->Get("2b_emu_OS_jets_n");
  TH1 *h_mc_bjets_n = (TH1*)mc_hists_file->Get("2b_emu_OS_bjets_n");
  TH1 *h_met_mc = (TH1*)mc_hists_file->Get("2b_emu_OS_met");
  TH1 *h_jet_pt_mc[3];
  for (int i=0; i<3; i++) {
    TString title = "2b_emu_OS_jet_pt_" + to_string(i);
    h_jet_pt_mc[i] = (TH1*)mc_hists_file->Get(title); }
  TH1 *h_mc_NN__2b_emu_OS_m_bjet_lep_min_dR = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_m_bjet_lep_min_dR");
  TH1 *h_mc_NN__2b_emu_OS_m_bjet_lep_min = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_m_bjet_lep_min");
  TH1 *h_mc_NN__2b_emu_OS_m_bjet_lep_max = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_m_bjet_lep_max");
  TH1 *h_mc_NN__2b_emu_OS_m_bjet_el = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_m_bjet_el");
  TH1 *h_mc_NN__2b_emu_OS_m_bjet_mu = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_m_bjet_mu");
  TH1 *h_mc_NN__2b_emu_OS_dR_bjet_lep0 = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_dR_bjet_lep0");
  TH1 *h_mc_NN__2b_emu_OS_dR_bjet_lep1 = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_dR_bjet_lep1");
  TH1 *h_mc_NN__2b_emu_OS_min_dR_bjet_lep = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_min_dR_bjet_lep");
  TH1 *h_mc_NN__2b_emu_OS_min_dR_jet_bjet = (TH1*)mc_hists_file->Get("NN__2b_emu_OS_min_dR_jet_bjet");
  

  // Get data hists
  TFile *data_hists_file = new TFile("hists_data.root");
  TH1 *h_data_jets_n = (TH1*)data_hists_file->Get("2b_emu_OS_jets_n");
  TH1 *h_data_bjets_n = (TH1*)data_hists_file->Get("2b_emu_OS_bjets_n");
  TH1 *h_met_data = (TH1*)data_hists_file->Get("2b_emu_OS_met");
  TH1 *h_jet_pt_data[3];
  for (int i=0; i<3; i++) {
    TString title = "2b_emu_OS_jet_pt_" + to_string(i);
    h_jet_pt_data[i] = (TH1*)data_hists_file->Get(title); }
  TH1 *h_data_NN__2b_emu_OS_m_bjet_lep_min_dR = (TH1*)data_hists_file->Get("NN__2b_emu_OS_m_bjet_lep_min_dR");
  TH1 *h_data_NN__2b_emu_OS_m_bjet_lep_min = (TH1*)data_hists_file->Get("NN__2b_emu_OS_m_bjet_lep_min");
  TH1 *h_data_NN__2b_emu_OS_m_bjet_lep_max = (TH1*)data_hists_file->Get("NN__2b_emu_OS_m_bjet_lep_max");
  TH1 *h_data_NN__2b_emu_OS_m_bjet_el = (TH1*)data_hists_file->Get("NN__2b_emu_OS_m_bjet_el");
  TH1 *h_data_NN__2b_emu_OS_m_bjet_mu = (TH1*)data_hists_file->Get("NN__2b_emu_OS_m_bjet_mu");
  TH1 *h_data_NN__2b_emu_OS_dR_bjet_lep0 = (TH1*)data_hists_file->Get("NN__2b_emu_OS_dR_bjet_lep0");
  TH1 *h_data_NN__2b_emu_OS_dR_bjet_lep1 = (TH1*)data_hists_file->Get("NN__2b_emu_OS_dR_bjet_lep1");
  TH1 *h_data_NN__2b_emu_OS_min_dR_bjet_lep = (TH1*)data_hists_file->Get("NN__2b_emu_OS_min_dR_bjet_lep");
  TH1 *h_data_NN__2b_emu_OS_min_dR_jet_bjet = (TH1*)data_hists_file->Get("NN__2b_emu_OS_min_dR_jet_bjet");

  
  // Draw different hists
  
  int draw_met = draw_data_mc_plot(h_met_data, h_met_mc, "#bf{E_{T}^{miss}, GeV}", "2b_emu_OS_met");
  
  for (int i=0; i<3; i++) {
    TString title = "#bf{p_{T}^{jet-" + to_string(i) + "}, GeV}";
    TString savename = "2b_emu_OS_jet_pt_" + to_string(i);
    int draw_jet_pt_0 = draw_data_mc_plot(h_jet_pt_data[i], h_jet_pt_mc[i], title, savename); }
  
  int draw_jets_n = draw_data_mc_plot(h_data_jets_n, h_mc_jets_n, "#bf{N_{jets}}", "2b_emu_OS_jets_n");
  int draw_bjets_n = draw_data_mc_plot(h_data_bjets_n, h_mc_bjets_n, "#bf{N_{bjets}}", "2b_emu_OS_bjets_n");
  

  int draw_NN__2b_emu_OS_m_bjet_lep_min_dR = draw_data_mc_plot(h_data_NN__2b_emu_OS_m_bjet_lep_min_dR, h_mc_NN__2b_emu_OS_m_bjet_lep_min_dR, "#bf{m_{inv}^{min dR}(bjet, lep)}", "2b_emu_OS_m_bjet_lep_min_dR");
  int draw_NN__2b_emu_OS_m_bjet_lep_min = draw_data_mc_plot(h_data_NN__2b_emu_OS_m_bjet_lep_min, h_mc_NN__2b_emu_OS_m_bjet_lep_min, "#bf{m_{inv}^{min}(bjet, lep)}", "2b_emu_OS_m_bjet_lep_min");
  int draw_NN__2b_emu_OS_m_bjet_lep_max = draw_data_mc_plot(h_data_NN__2b_emu_OS_m_bjet_lep_max, h_mc_NN__2b_emu_OS_m_bjet_lep_max, "#bf{m_{inv}^{max}(bjet, lep)}", "2b_emu_OS_m_bjet_lep_max");
  int draw_NN__2b_emu_OS_m_bjet_el = draw_data_mc_plot(h_data_NN__2b_emu_OS_m_bjet_el, h_mc_NN__2b_emu_OS_m_bjet_el, "#bf{m_{inv}(bjet, el)}", "2b_emu_OS_m_bjet_el");
  int draw_NN__2b_emu_OS_m_bjet_mu = draw_data_mc_plot(h_data_NN__2b_emu_OS_m_bjet_mu, h_mc_NN__2b_emu_OS_m_bjet_mu, "#bf{m_{inv}(bjet, mu)}", "2b_emu_OS_m_bjet_mu");
  int draw_NN__2b_emu_OS_dR_bjet_lep0 = draw_data_mc_plot(h_data_NN__2b_emu_OS_dR_bjet_lep0, h_mc_NN__2b_emu_OS_dR_bjet_lep0, "#bf{dR(bjet, lep0)}", "2b_emu_OS_dR_bjet_lep0");
  int draw_NN__2b_emu_OS_dR_bjet_lep1 = draw_data_mc_plot(h_data_NN__2b_emu_OS_dR_bjet_lep1, h_mc_NN__2b_emu_OS_dR_bjet_lep1, "#bf{dR(bjet, lep1)}", "2b_emu_OS_dR_bjet_lep1");
  int draw_NN__2b_emu_OS_min_dR_bjet_lep = draw_data_mc_plot(h_data_NN__2b_emu_OS_min_dR_bjet_lep, h_mc_NN__2b_emu_OS_min_dR_bjet_lep, "#bf{dR_{min}(bjet-fix, lep)}", "2b_emu_OS_min_dR_bjet_lep");
  int draw_NN__2b_emu_OS_min_dR_jet_bjet = draw_data_mc_plot(h_data_NN__2b_emu_OS_min_dR_jet_bjet, h_mc_NN__2b_emu_OS_min_dR_jet_bjet, "#bf{dR_{min}(jet-fix, bjet)}", "2b_emu_OS_min_dR_jet_bjet");
}
