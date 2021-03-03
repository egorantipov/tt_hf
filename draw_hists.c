#include <TH2.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TMath.h>
#include <TColor.h>

#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

//Colors: red , green, blue, orange, cyan, magenta, yellow
vector<Int_t> colors = {632, 416+1, 600, 800-3, 432+2, 616+1, 400+1};



// ###########################
// ## Draw a few histograms ## 
// ###########################
int draw_n_histos(vector<TH1*> h_vec, vector<TString> h_title, TString x_axis_title, TString title, bool normalize=false, Double_t y_min=0, Double_t y_max=10000)
{
  // Draws N histogram on one canvas (not stacked)
  cout << "Start drawing " << title << " !" << endl;
  if (h_vec.size()==0) { cout << "h_vec is emmpty, aborting!!!" << endl; return 0; }
  
  
  TCanvas *c = new TCanvas(h_title[0], h_title[0], 1600, 1200);
  gStyle->SetOptStat(0);
  gPad->SetGrid();
  if (normalize==false) gPad->SetLogy();
  double legend_height = 0.09*h_vec.size();
  double legend_y1 = 0.90 - legend_height;
  TLegend *legend = new TLegend(0.70, legend_y1, 0.90, 0.90);

  for (int i=0; i<h_vec.size(); i++){
    // TEST 
    if (!h_vec[i]) cout << "Requested object TH_[" << i << "] wasn't found!" << endl;
    
    double h_int = h_vec[i]->Integral(1, h_vec[i]->GetNbinsX());
    double sf = 1/h_int;
    cout << "Hist: " << h_title[i] << " , Integral: " << h_int << endl;
    
    h_vec[i]->SetMarkerStyle(20);
    h_vec[i]->SetMarkerSize(2);
    h_vec[i]->SetMarkerColor(colors[i]);
    //h_vec[i]->SetLineColor(i+2);
    h_vec[i]->SetLineColor(colors[i]);
    h_vec[i]->SetLineWidth(4);
    if (normalize==true) h_vec[i]->Scale(sf);
    
    if (i==0) {
      h_vec[i]->Draw("hist");
      h_vec[i]->SetTitle(title);
      
      if (normalize==false) {
	h_vec[i]->GetYaxis()->SetRangeUser(y_min, y_max);
	h_vec[i]->GetYaxis()->SetTitle("#bf{Events}"); }
      
      else {
	if (y_min==0) y_min = 1;
	h_vec[i]->GetYaxis()->SetRangeUser(y_min, y_max);
	h_vec[i]->GetYaxis()->SetTitle("#bf{Events norm to 1}"); 
	gPad->SetLogy(); }
      
      h_vec[i]->GetXaxis()->SetTitle(x_axis_title); }
    
    else { h_vec[i]->Draw("hist same"); }
    legend->AddEntry(h_vec[i], h_title[i]); }
  legend->Draw("same");
  
  c->Print("Plots/" + title + ".png");
  cout << "Drawn " + title + " !\n\n" << endl;
  
  return 0;
}



// ################
// ##    MAIN    ##
// ################
void draw_hists()
{
  // Open the file with histograms
  TFile *hists_file_mc = TFile::Open("hists_mc.root");
  
  
  
  // Get the histograms

  // DL1r tag weights
  vector<TString> processes = {"2b1l", "4b", "3b", "2b1c"};
  TH1 *mc16_tag0_DL1r[4];
  TH1 *mc16_tag1_DL1r[4];
  TH1 *mc16_tag2_DL1r[4];
  TH1 *mc16_bjets_n_2b[4];
  for (int topHFFF_i=0; topHFFF_i<4; topHFFF_i++) {
    mc16_tag0_DL1r[topHFFF_i] = (TH1*)hists_file_mc->Get("DL1r_templates_"+processes[topHFFF_i]+"_1st_tag");
    mc16_tag1_DL1r[topHFFF_i] = (TH1*)hists_file_mc->Get("DL1r_templates_"+processes[topHFFF_i]+"_2nd_tag");
    mc16_tag2_DL1r[topHFFF_i] = (TH1*)hists_file_mc->Get("DL1r_templates_"+processes[topHFFF_i]+"_3rd_tag");
    
    mc16_bjets_n_2b[topHFFF_i] = (TH1*)hists_file_mc->Get("2b_emu_OS_"+processes[topHFFF_i]+"bjets_n"); }

  // jet_pT
  TH1 *mc16_jet_pT0_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_jet_pT0_from_top");
  TH1 *mc16_jet_pT1_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_jet_pT1_from_top");
  TH1 *mc16_jet_pT2_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_jet_pT2_from_top");
  TH1 *mc16_jet_pT0_not_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_jet_pT0_not_from_top");
  TH1 *mc16_jet_pT1_not_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_jet_pT1_not_from_top");
  TH1 *mc16_jet_pT2_not_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_jet_pT2_not_from_top");
  TH1 *mc16_jet_pT0 = (TH1*)hists_file_mc->Get("2b_emu_OS_jet_pt0");
  TH1 *mc16_jet_pT1 = (TH1*)hists_file_mc->Get("2b_emu_OS_jet_pt1");
  TH1 *mc16_jet_pT2 = (TH1*)hists_file_mc->Get("2b_emu_OS_jet_pt2");
  TH1 *mc16_jet_pT3 = (TH1*)hists_file_mc->Get("2b_emu_OS_jet_pt3");
  TH1 *mc16_jet_pT4 = (TH1*)hists_file_mc->Get("2b_emu_OS_jet_pt4");
  TH1 *mc16_jet_pT5 = (TH1*)hists_file_mc->Get("2b_emu_OS_jet_pt5");

  // MET
  TH1 *mc16_met = (TH1*)hists_file_mc->Get("2b_emu_OS_met");
  TH1 *mc16_met_phi = (TH1*)hists_file_mc->Get("2b_emu_OS_met_phi");

  // leptons variables
  TH1 *mc16_lep0_pt = (TH1*)hists_file_mc->Get("2b_emu_OS_lep0_pt");
  TH1 *mc16_lep1_pt = (TH1*)hists_file_mc->Get("2b_emu_OS_lep1_pt");
  TH1 *mc16_lep_pt = (TH1*)hists_file_mc->Get("2b_emu_OS_lep_pt");
  TH1 *mc16_lep0_eta = (TH1*)hists_file_mc->Get("2b_emu_OS_lep0_eta");
  TH1 *mc16_lep1_eta = (TH1*)hists_file_mc->Get("2b_emu_OS_lep1_eta");
  TH1 *mc16_lep_eta = (TH1*)hists_file_mc->Get("2b_emu_OS_lep_eta");
  TH1 *mc16_lep0_phi = (TH1*)hists_file_mc->Get("2b_emu_OS_lep0_phi");
  TH1 *mc16_lep1_phi = (TH1*)hists_file_mc->Get("2b_emu_OS_lep1_phi");
  TH1 *mc16_lep_phi = (TH1*)hists_file_mc->Get("2b_emu_OS_lep_phi");
  TH1 *mc16_dR_lep0_lep1 = (TH1*)hists_file_mc->Get("2b_emu_OS_dR_lep0_lep1");
  
  // dR bN and bM
  TH1 *mc16_dR_b0_b1 = (TH1*)hists_file_mc->Get("2b_emu_OS_dR_b0_b1");
  TH1 *mc16_dR_b0_b2 = (TH1*)hists_file_mc->Get("2b_emu_OS_dR_b0_b2");
  TH1 *mc16_dR_b1_b2 = (TH1*)hists_file_mc->Get("2b_emu_OS_dR_b1_b2");
  TH1 *mc16_minDeltaR_b01_b2_from_top = (TH1*)hists_file_mc->Get("2b_emu_OS_min_dR_b01_b2_from_top");
  TH1 *mc16_minDeltaR_b01_b2_not_from_top = (TH1*)hists_file_mc->Get("2b_emu_OS_min_dR_b01_b2_not_from_top");
  
  // min_dR b and leptons
  TH1 *mc16_dR_min_lep0_b_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_min_dR_lep0_b_from_top");
  TH1 *mc16_dR_min_lep0_b_not_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_min_dR_lep0_b_not_from_top");
  TH1 *mc16_dR_min_lep1_b_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_min_dR_lep1_b_from_top");
  TH1 *mc16_dR_min_lep1_b_not_from_top = (TH1*)hists_file_mc->Get("3b_emu_OS_min_dR_lep1_b_not_from_top");
  TH1 *mc16_dR_min_lep0_b_from_top_tags[3];
  TH1 *mc16_dR_min_lep0_b_not_from_top_tags[3];
  TH1 *mc16_dR_min_lep1_b_from_top_tags[3];
  TH1 *mc16_dR_min_lep1_b_not_from_top_tags[3];
  for (int i=0; i<3; i++) {
    TString title1 = "2b_emu_OS_min_dR_lep0_b_from_top_" + to_string(i+1) + "tag";
    TString title2 = "2b_emu_OS_min_dR_lep0_b_not_from_top_" + to_string(i+1) + "tag";
    TString title3 = "2b_emu_OS_min_dR_lep1_b_from_top_" + to_string(i+1) + "tag";
    TString title4 = "2b_emu_OS_min_dR_lep1_b_not_from_top_" + to_string(i+1) + "tag"; 
    mc16_dR_min_lep0_b_from_top_tags[i] = (TH1*)hists_file_mc->Get(title1);
    mc16_dR_min_lep0_b_not_from_top_tags[i] = (TH1*)hists_file_mc->Get(title2);
    mc16_dR_min_lep1_b_from_top_tags[i] = (TH1*)hists_file_mc->Get(title3);
    mc16_dR_min_lep1_b_not_from_top_tags[i] = (TH1*)hists_file_mc->Get(title4);
}

  // mid_dR(obj, obj) proposed by Sasha
  TH1 *mc16_minDeltaR_b_from_top_to_b = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_b_from_top_to_b");
  TH1 *mc16_minDeltaR_b_not_from_top_to_b = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_b_not_from_top_to_b");
  TH1 *mc16_minDeltaR_not_b_to_b = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_not_b_to_b");
  TH1 *mc16_minDeltaR_b_from_top_to_jet = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_b_from_top_to_jet");
  TH1 *mc16_minDeltaR_b_not_from_top_to_jet = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_b_not_from_top_to_jet");
  TH1 *mc16_minDeltaR_not_b_to_jet = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_not_b_to_jet");
  TH1 *mc16_minDeltaR_b_from_top_to_lep = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_b_from_top_to_lep");
  TH1 *mc16_minDeltaR_b_not_from_top_to_lep = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_b_not_from_top_to_lep");
  TH1 *mc16_minDeltaR_not_b_to_lep = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_minDeltaR_not_b_to_lep");

  // Invariant mass
  TH1 *mc16_inv_mass_lep_bjet_from_top_min_dR = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_inv_mass_lep_bjet_from_top_min_dR");
  TH1 *mc16_inv_mass_lep_bjet_not_from_top_min_dR = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_inv_mass_lep_bjet_not_from_top_min_dR");
  TH1 *mc16_inv_mass_lep_btag_from_top_min_dR = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_inv_mass_lep_btag_from_top_min_dR");
  TH1 *mc16_inv_mass_lep_btag_not_from_top_min_dR = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_inv_mass_lep_btag_not_from_top_min_dR");
  TH1 *mc16_min_inv_mass_lep_bjet_from_top = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_min_inv_mass_lep_bjet_from_top");
  TH1 *mc16_max_inv_mass_lep_bjet_from_top = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_max_inv_mass_lep_bjet_from_top");
  TH1 *mc16_min_inv_mass_lep_bjet_not_from_top = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_min_inv_mass_lep_bjet_not_from_top");
  TH1 *mc16_max_inv_mass_lep_bjet_not_from_top = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_max_inv_mass_lep_bjet_not_from_top");
  TH1 *mc16_min_inv_mass_lep_other_jet = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_min_inv_mass_lep_other_jet");
  TH1 *mc16_max_inv_mass_lep_other_jet = (TH1*)hists_file_mc->Get("2b_emu_OS_mc16_max_inv_mass_lep_other_jet");

  // topHFFF overlap removal studies
  TH1 *mc16_bjets_n_411076 = (TH1*)hists_file_mc->Get("topHFFF_study_2b_emu_OS_bjets_n_411076");
  TH1 *mc16_bjets_n_411077 = (TH1*)hists_file_mc->Get("topHFFF_study_2b_emu_OS_bjets_n_411077");
  TH1 *mc16_bjets_n_411078 = (TH1*)hists_file_mc->Get("topHFFF_study_2b_emu_OS_bjets_n_411078");
  TH1 *mc16_bjets_n_410472 = (TH1*)hists_file_mc->Get("topHFFF_study_2b_emu_OS_bjets_n_410472");

  
  // Draw the plots in two steps:
  // (1) prepare different combinations of hists (vectors) to draw them on one canvas
  // (2) call a function to draw the collection of hists on one canvas
  /*
  // First three tags for same process
  vector<TString> mc16_DL1r_topHFFF_3tags_title = {"1st tag", "2nd tag", "3rd tag"};
  vector<TH1*> mc16_DL1r_topHFFF0_3tags = {mc16_tag0_DL1r[0], mc16_tag1_DL1r[0], mc16_tag2_DL1r[0]};
  vector<TH1*> mc16_DL1r_topHFFF1_3tags = {mc16_tag0_DL1r[1], mc16_tag1_DL1r[1], mc16_tag2_DL1r[1]};
  vector<TH1*> mc16_DL1r_topHFFF2_3tags = {mc16_tag0_DL1r[2], mc16_tag1_DL1r[2], mc16_tag2_DL1r[2]};
  vector<TH1*> mc16_DL1r_topHFFF3_3tags = {mc16_tag0_DL1r[3], mc16_tag1_DL1r[3], mc16_tag2_DL1r[3]};
  int mc16_DL1r_topHFFF0_3tags_draw = draw_n_histos(mc16_DL1r_topHFFF0_3tags, mc16_DL1r_topHFFF_3tags_title, "#bf{DL1r tag weight}", "2b1l_first_three_tags", true, 0, 0.6);
  int mc16_DL1r_topHFFF1_3tags_draw = draw_n_histos(mc16_DL1r_topHFFF1_3tags, mc16_DL1r_topHFFF_3tags_title, "#bf{DL1r tag weight}", "4b_first_three_tags", true, 0, 0.6);
  int mc16_DL1r_topHFFF2_3tags_draw = draw_n_histos(mc16_DL1r_topHFFF2_3tags, mc16_DL1r_topHFFF_3tags_title, "#bf{DL1r tag weight}", "3b_first_three_tags", true, 0, 0.6);
  int mc16_DL1r_topHFFF3_3tags_draw = draw_n_histos(mc16_DL1r_topHFFF3_3tags, mc16_DL1r_topHFFF_3tags_title, "#bf{DL1r tag weight}", "2b1c_first_three_tag", true, 0, 0.6);


  // 1st, 2nd and 3rd tag for different processes
  vector<TH1*> mc16_DL1r_1st_tag_all = {mc16_tag0_DL1r[0], mc16_tag0_DL1r[1], mc16_tag0_DL1r[2], mc16_tag0_DL1r[3]};
  vector<TH1*> mc16_DL1r_2nd_tag_all = {mc16_tag1_DL1r[0], mc16_tag1_DL1r[1], mc16_tag1_DL1r[2], mc16_tag1_DL1r[3]};
  vector<TH1*> mc16_DL1r_3rd_tag_all = {mc16_tag2_DL1r[0], mc16_tag2_DL1r[1], mc16_tag2_DL1r[2], mc16_tag2_DL1r[3]};
  int mc16_DL1r_1st_tag_all_draw = draw_n_histos(mc16_DL1r_1st_tag_all, processes, "#bf{DL1r, 1st tag weight}", "1st_tag", true, 0, 0.6);
  int mc16_DL1r_2nd_tag_all_draw = draw_n_histos(mc16_DL1r_2nd_tag_all, processes, "#bf{DL1r, 2nd tag weight}", "2nd_tag", true, 0, 0.6);
  int mc16_DL1r_3rd_tag_all_draw = draw_n_histos(mc16_DL1r_3rd_tag_all, processes, "#bf{DL1r, 3rd tag weight}", "3rd_tag", true, 0, 0.6);


  // bjets_n for the 2+b channel
  vector<TH1*> mc16_bjets_n_2b_vec = {mc16_bjets_n_2b[0], mc16_bjets_n_2b[1], mc16_bjets_n_2b[2], mc16_bjets_n_2b[3]};
  int mc16_bjets_n_2b_draw = draw_n_histos(mc16_bjets_n_2b_vec, processes, "#bf{N_{bjets}}", "bjets_n", true, 0, 1.1);


  // Regular dR_min between leptons and bjets
  vector<TString> dR_min_title = {"From top", "Not from top"};
  vector<TH1*> mc16_min_dR_lep0_b = {mc16_dR_min_lep0_b_from_top, mc16_dR_min_lep0_b_not_from_top};
  vector<TH1*> mc16_min_dR_lep1_b = {mc16_dR_min_lep1_b_from_top, mc16_dR_min_lep1_b_not_from_top};
  int mc16_min_dR_lep0_b_draw = draw_n_histos(mc16_min_dR_lep0_b, dR_min_title, "#bf{min_dR_lep0_bjets}", "min_dR_lep0_bjets", true, 0, 0.6);
  int mc16_min_dR_lep1_b_draw = draw_n_histos(mc16_min_dR_lep1_b, dR_min_title, "#bf{min_dR_lep1_bjets}", "min_dR_lep1_bjets", true, 0, 0.6);


  // dR_min for the first three tag weights
  vector<TH1*> mc16_min_dR_lep0_b_1st_tag = {mc16_dR_min_lep0_b_from_top_tags[0], mc16_dR_min_lep0_b_not_from_top_tags[0]};
  vector<TH1*> mc16_min_dR_lep0_b_2nd_tag = {mc16_dR_min_lep0_b_from_top_tags[1], mc16_dR_min_lep0_b_not_from_top_tags[1]};
  vector<TH1*> mc16_min_dR_lep0_b_3rd_tag = {mc16_dR_min_lep0_b_from_top_tags[2], mc16_dR_min_lep0_b_not_from_top_tags[2]};
  vector<TH1*> mc16_min_dR_lep1_b_1st_tag = {mc16_dR_min_lep1_b_from_top_tags[0], mc16_dR_min_lep1_b_not_from_top_tags[0]};
  vector<TH1*> mc16_min_dR_lep1_b_2nd_tag = {mc16_dR_min_lep1_b_from_top_tags[1], mc16_dR_min_lep1_b_not_from_top_tags[1]};
  vector<TH1*> mc16_min_dR_lep1_b_3rd_tag = {mc16_dR_min_lep1_b_from_top_tags[2], mc16_dR_min_lep1_b_not_from_top_tags[2]};
  int mc16_min_dR_lep0_b_1st_tag_draw = draw_n_histos(mc16_min_dR_lep0_b_1st_tag, dR_min_title, "#bf{min_dR_lep0_bjets}", "min_dR_lep0_bjets_1st_tag", true, 0, 0.6);
  int mc16_min_dR_lep0_b_2rd_tag_draw = draw_n_histos(mc16_min_dR_lep0_b_2nd_tag, dR_min_title, "#bf{min_dR_lep0_bjets}", "min_dR_lep0_bjets_2nd_tag", true, 0, 0.6);
  int mc16_min_dR_lep0_b_3rd_tag_draw = draw_n_histos(mc16_min_dR_lep0_b_3rd_tag, dR_min_title, "#bf{min_dR_lep0_bjets}", "min_dR_lep0_bjets_3rd_tag", true, 0, 0.6);
  int mc16_min_dR_lep1_b_1st_tag_draw = draw_n_histos(mc16_min_dR_lep0_b_1st_tag, dR_min_title, "#bf{min_dR_lep1_bjets}", "min_dR_lep1_bjets_1st_tag", true, 0, 0.6);
  int mc16_min_dR_lep1_b_2nd_tag_draw = draw_n_histos(mc16_min_dR_lep0_b_2nd_tag, dR_min_title, "#bf{min_dR_lep1_bjets}", "min_dR_lep1_bjets_2nd_tag", true, 0, 0.6);
  int mc16_min_dR_lep1_b_3rd_tag_draw = draw_n_histos(mc16_min_dR_lep0_b_3rd_tag, dR_min_title, "#bf{min_dR_lep1_bjets}", "min_dR_lep1_bjets_3rd_tag", true, 0, 0.6);
  
  vector<TString> dR_min_tags_title = {"1st tag", "2nd tag", "3rd_tag"};
  vector<TH1*> mc16_min_dR_lep0_bjets_from_top_three_tags = {mc16_dR_min_lep0_b_from_top_tags[0], mc16_dR_min_lep0_b_from_top_tags[1], mc16_dR_min_lep0_b_from_top_tags[2]};
  vector<TH1*> mc16_min_dR_lep0_bjets_not_from_top_three_tags = {mc16_dR_min_lep0_b_not_from_top_tags[0], mc16_dR_min_lep0_b_not_from_top_tags[1], mc16_dR_min_lep0_b_not_from_top_tags[2]};
  vector<TH1*> mc16_min_dR_lep1_bjets_from_top_three_tags = {mc16_dR_min_lep1_b_from_top_tags[0], mc16_dR_min_lep1_b_from_top_tags[1], mc16_dR_min_lep1_b_from_top_tags[2]};
  vector<TH1*> mc16_min_dR_lep1_bjets_not_from_top_three_tags = {mc16_dR_min_lep1_b_not_from_top_tags[0], mc16_dR_min_lep1_b_not_from_top_tags[1], mc16_dR_min_lep1_b_not_from_top_tags[2]};
  int mc16_min_dR_lep0_bjets_from_top_3t_draw = draw_n_histos(mc16_min_dR_lep0_bjets_from_top_three_tags, dR_min_tags_title, "#bf{min_dR_lep0_bjets}", "min_dR_lep0_bjets_from_top_three_tags", true, 0, 0.6);
  int mc16_min_dR_lep0_bjets_not_from_top_3t_draw = draw_n_histos(mc16_min_dR_lep0_bjets_not_from_top_three_tags, dR_min_tags_title, "#bf{min_dR_lep0_bjets}", "min_dR_lep0_bjets_not_from_top_three_tags", true, 0, 0.6);
  int mc16_min_dR_lep1_bjets_from_top_3t_draw = draw_n_histos(mc16_min_dR_lep1_bjets_from_top_three_tags, dR_min_tags_title, "#bf{min_dR_lep1_bjets}", "min_dR_lep1_bjets_from_top_three_tags", true, 0, 0.6);
  int mc16_min_dR_lep1_bjets_not_from_top_3t_draw = draw_n_histos(mc16_min_dR_lep1_bjets_not_from_top_three_tags, dR_min_tags_title, "#bf{min_dR_lep1_bjets}", "min_dR_lep1_bjets_not_from_top_tree_tags", true, 0, 0.6);


  // dR between bN and bM
  vector<TString> dR_bN_bM_title = {"1st-2nd tags", "1st-3rd tags", "2nd-3rd tags"};
  vector<TH1*> mc16_dR_bN_bM = {mc16_dR_b0_b1, mc16_dR_b0_b2, mc16_dR_b1_b2};
  int mc16_dR_bN_bM_draw = draw_n_histos(mc16_dR_bN_bM, dR_bN_bM_title, "#bf{dR for two tags (DL1r)}", "dR_bN_bM", true, 0, 0.6);


  // min_dR between the 3rd tag and the (1st/2nd) tag, (not)from top:
  vector<TString> min_dR_b01_b2_title = {"from top", "not from top"};
  vector<TH1*> mc16_minDeltaR_b01_b2_collection = {mc16_minDeltaR_b01_b2_from_top, mc16_minDeltaR_b01_b2_not_from_top};
  int mc16_minDeltaR_b01_b2_collection_draw = draw_n_histos(mc16_minDeltaR_b01_b2_collection, min_dR_b01_b2_title, "#bf{min dR 3rd tag to 1st/2nd tag}", "min_dR_b01_b2", true, 0, 0.6);


  // MET
  vector<TString> met_title = {"MET"};
  vector<TH1*> mc16_met_collection = {mc16_met};
  int mc16_met_draw = draw_n_histos(mc16_met_collection, met_title, "#bf{met, #it{GeV}}", "met", true, 0, 0.6);
  vector<TString> met_phi_title = {"MET Phi"};
  vector<TH1*> mc16_met_phi_collection = {mc16_met_phi};
  int mc16_met_phi_draw = draw_n_histos(mc16_met_phi_collection, met_phi_title, "#bf{met phi}", "met_phi", true, 0, 0.6);

  
  // jet pT in the reverce (pT) order
  vector<TString> jet_pt_title = {"1st", "2nd", "3rd", "4th", "5th", "6th"};
  vector<TH1*> mc16_jet_pt_collection = {mc16_jet_pT0, mc16_jet_pT1, mc16_jet_pT2, mc16_jet_pT3, mc16_jet_pT4, mc16_jet_pT5};
  int mc16_jet_pTN_draw = draw_n_histos(mc16_jet_pt_collection, jet_pt_title, "#bf{jet pT, #it{GeV}}", "jet_pt_collection", true, 0, 0.6);


  // Leptons
  vector<TString> lep_all_title = {"1st", "2nd", "combined"};
  vector<TH1*> mc16_lepN_pt_collection = {mc16_lep0_pt, mc16_lep1_pt, mc16_lep_pt};
  int mc16_lepN_pt_draw = draw_n_histos(mc16_lepN_pt_collection, lep_all_title, "#bf{lep. pT, #it{GeV}}", "lep_pT", true, 0, 0.6);
  
  vector<TH1*> mc16_lepN_eta_collection = {mc16_lep0_eta, mc16_lep1_eta, mc16_lep_eta};
  int mc16_lepN_eta_draw = draw_n_histos(mc16_lepN_eta_collection, lep_all_title, "#bf{lep. eta}", "lep_eta", true, 0, 0.6);
  
  vector<TH1*> mc16_lepN_phi_collection = {mc16_lep0_phi, mc16_lep1_phi, mc16_lep_phi};
  int mc16_lepN_phi_draw = draw_n_histos(mc16_lepN_phi_collection, lep_all_title, "#bf{lep. phi}", "lep_phi", true, 0, 0.6);

  vector<TString> lep_dR_title = {"dR lep0 - lep1"};
  vector<TH1*> mc16_dR_lep0_lep1_collection = {mc16_dR_lep0_lep1};
  int mc16_dR_lep0_lep1_draw = draw_n_histos(mc16_dR_lep0_lep1_collection, lep_dR_title, "#bf{dR lep0-lep1}","dR_lep0_lep1",  true, 0, 0.6);


  // mid_dR(obj, obj) proposed by Sasha
  vector<TString> dR_jet_obj_title = {"other non-b-jet", "other b-jet", "b-jet from top"};
  
  vector<TH1*> mc16_min_dR_jet_bjet_collection = {mc16_minDeltaR_not_b_to_b, mc16_minDeltaR_b_not_from_top_to_b, mc16_minDeltaR_b_from_top_to_b};
  int  mc16_min_dR_jet_bjet_draw = draw_n_histos(mc16_min_dR_jet_bjet_collection, dR_jet_obj_title, "#bf{#DeltaR_{min}(jet, b-jet)}", "min_dR_jet_bjet", true, 0, 0.2);
  
  vector<TH1*> mc16_min_dR_jet_jet_collection = {mc16_minDeltaR_not_b_to_jet, mc16_minDeltaR_b_not_from_top_to_jet, mc16_minDeltaR_b_from_top_to_jet};
  int mc16_min_dR_jet_jet_draw = draw_n_histos(mc16_min_dR_jet_jet_collection, dR_jet_obj_title, "#bf{#DeltaR_{min}(jet, jet)}", "min_dR_jet_jet", true, 0, 0.2);

  vector<TH1*> mc16_min_dR_je_lep_collection = {mc16_minDeltaR_not_b_to_lep, mc16_minDeltaR_b_not_from_top_to_lep, mc16_minDeltaR_b_from_top_to_lep};
  int mc16_min_dR_jet_lep_draw = draw_n_histos(mc16_min_dR_je_lep_collection, dR_jet_obj_title, "#bf{#DeltaR_{min}(jet, lepton)}", "min_dR_jet_lep", true, 0, 0.2);


  // Invariant mass
  vector<TString> mc16_inv_mass_lep_bjet_min_dR_title = {"from top", "not from top"};
  vector<TH1*> mc16_inv_mass_lep_bjet_min_dR_collection = {mc16_inv_mass_lep_bjet_from_top_min_dR, mc16_inv_mass_lep_bjet_not_from_top_min_dR};
  int mc16_inv_mass_lep_bjet_min_dR_draw = draw_n_histos(mc16_inv_mass_lep_bjet_min_dR_collection, mc16_inv_mass_lep_bjet_min_dR_title, "#bf{M(bjet - closest lep.)^{inv}}", "inv_mass_lep_bjet_min_dR", true, 0, 0.02);

  vector<TH1*> mc16_inv_mass_lep_btag_min_dR_collection = {mc16_inv_mass_lep_btag_from_top_min_dR, mc16_inv_mass_lep_btag_not_from_top_min_dR};
  int mc16_inv_mass_lep_btag_min_dR_draw = draw_n_histos(mc16_inv_mass_lep_btag_min_dR_collection, mc16_inv_mass_lep_bjet_min_dR_title, "#bf{M(btag - closest lep.)^{inv}}", "inv_mass_lep_btag_min_dR", true, 0, 0.02);

  vector<TString> mc16_inv_mass_lep_obj_title = {"other non-b-jet", "other b-jet", "b-jet from top"};
  vector<TH1*> mc16_min_inv_mass_lep_obj_collection = {mc16_min_inv_mass_lep_other_jet, mc16_min_inv_mass_lep_bjet_not_from_top, mc16_min_inv_mass_lep_bjet_from_top};
  int mc16_min_inv_mass_lep_obj_draw = draw_n_histos(mc16_min_inv_mass_lep_obj_collection, mc16_inv_mass_lep_obj_title, "#bf{M(jet - lep.)^{inv}_{min}}", "inv_mass_min_lep_obj", true, 0, 0.02);

  vector<TH1*> mc16_max_inv_mass_lep_obj_collection = {mc16_max_inv_mass_lep_other_jet, mc16_max_inv_mass_lep_bjet_not_from_top, mc16_max_inv_mass_lep_bjet_from_top};
  int mc16_max_inv_mass_lep_obj_draw = draw_n_histos(mc16_max_inv_mass_lep_obj_collection, mc16_inv_mass_lep_obj_title, "#bf{M(jet - lep.)^{inv}_{max}}", "inv_mass_max_lep_obj", true, 0, 0.02);
  */

  // topHFFF overlar removal study
  vector<TString> mc16_bjets_n_4_title = {"411076 - 4b", "411077 - 3b", "411078 - 2b1c", "410472 - tt incl."};
  vector<TH1*> mc16_bjets_n_4_collection = {mc16_bjets_n_411076, mc16_bjets_n_411077, mc16_bjets_n_411078, mc16_bjets_n_410472};
  int mc16_bjets_n_4_draw = draw_n_histos(mc16_bjets_n_4_collection, mc16_bjets_n_4_title, "#bf{N_{bjets}}", "bjets_n_topHFFF_study", false, 1, pow(10,6));


  // Close the hists file
  hists_file_mc->Close();
}


