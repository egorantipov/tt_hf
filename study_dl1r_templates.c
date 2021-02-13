#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TPad.h>
#include <TMath.h>
#include <TFractionFitter.h>

#include <iostream>
#include <sstream>
#include <vector>
using namespace std;



// ###########################
// ## Draw a few histograms ##
// ###########################
int draw_n_histos(vector<TH1*> h_vec, vector<TString> h_title, TString x_axis_title, TString title, bool normalize=false, double y_min=0, double y_max=10000)
{
  // Draws N histogram on one canvas (not stacked)
  cout << "Start drawing " << title << " !" << endl;
  if (h_vec.size()==0) { cout << "h_vec is emmpty, aborting!!!" << endl; return 0; }


  TCanvas *c = new TCanvas(h_title[0], h_title[0], 1600, 900);
  gStyle->SetOptStat(0);
  gPad->SetGrid();
  if (normalize==false) gPad->SetLogy();
  double legend_height = 0.09*h_vec.size();
  double legend_y1 = 0.90 - legend_height;
  TLegend *legend = new TLegend(0.80, legend_y1, 0.90, 0.90);

  for (int i=0; i<h_vec.size(); i++) {
    // TEST                                                                                        
    if (!h_vec[i]) cout << "Requested object TH_[" << i << "] wasn't found!" << endl;
    
    double h_int = h_vec[i]->Integral(0, h_vec[i]->GetNbinsX());
    double sf = 1/h_int;
    
    h_vec[i]->SetMarkerStyle(20);
    h_vec[i]->SetMarkerSize(2);
    h_vec[i]->SetMarkerColor(i+1);
    h_vec[i]->SetLineColor(i+1);
    h_vec[i]->SetLineWidth(2);
    if (normalize==true) h_vec[i]->Scale(sf);
    if (i==0) {
      h_vec[i]->Draw("C");
      h_vec[i]->SetTitle(title);
      
      if (normalize==false){
	h_vec[i]->GetYaxis()->SetRangeUser(y_min, y_max);
	h_vec[i]->GetYaxis()->SetTitle("#bf{Events}"); }
      else {
	h_vec[i]->GetYaxis()->SetRangeUser(0, 0.6);
	h_vec[i]->GetYaxis()->SetTitle("#bf{Events norm to 1}"); }
      
      h_vec[i]->GetXaxis()->SetTitle(x_axis_title); }
    else { h_vec[i]->Draw("same C"); }
    
    legend->AddEntry(h_vec[i], h_title[i]); }
  
  legend->Draw("same");
  
  c->Print("Plots/fits/" + title + ".png");
  cout << "Drawn " + title + " !" << endl;
  
  return 0;
}



// ##############
// ##   MAIN   ##
// ##############
void study_dl1r_templates()
{
  // OPen the file with histograms
  TFile *hists_file_mc = TFile::Open("hists_mc.root");


  
  // Get required histograms
  // Scale them to unity right away
  vector<TString> processes = {"2b1l", "4b", "3b", "2b1c"};
  TH1 *mc16_tag0_DL1r[4];
  TH1 *mc16_tag1_DL1r[4];
  TH1 *mc16_tag2_DL1r[4];
  TH1 *mc16_tag2_DL1r_tmp[4];
  for (int topHFFF_i=0; topHFFF_i<4; topHFFF_i++) {
    mc16_tag0_DL1r[topHFFF_i] = (TH1*)hists_file_mc->Get("DL1r_templates_"+processes[topHFFF_i]+"_1st_tag");
    double tag0_sf = mc16_tag0_DL1r[topHFFF_i]->Integral(0,  mc16_tag0_DL1r[topHFFF_i]->GetNbinsX() + 1);
    mc16_tag0_DL1r[topHFFF_i]->Scale(1/tag0_sf);

    mc16_tag1_DL1r[topHFFF_i] = (TH1*)hists_file_mc->Get("DL1r_templates_"+processes[topHFFF_i]+"_2nd_tag");
    double tag1_sf = mc16_tag1_DL1r[topHFFF_i]->Integral(0, mc16_tag1_DL1r[topHFFF_i]->GetNbinsX() + 1);
    mc16_tag1_DL1r[topHFFF_i]->Scale(1/tag1_sf);
    
    mc16_tag2_DL1r[topHFFF_i] = (TH1*)hists_file_mc->Get("DL1r_templates_"+processes[topHFFF_i]+"_3rd_tag");
    double tag2_sf = mc16_tag2_DL1r[topHFFF_i]->Integral(0, mc16_tag2_DL1r[topHFFF_i]->GetNbinsX() + 1);
    mc16_tag2_DL1r[topHFFF_i]->Scale(1/tag2_sf); 
  }


  
  // Make mixtures of the taggers in known ratios
  // Reference order: 2b1l, 4b, 3b, 2b1c
  TH1 *mc16_tag0_DL1r_mix[9];
  TH1 *mc16_tag1_DL1r_mix[9];
  TH1 *mc16_tag2_DL1r_mix[9];
  vector<double> fraction_2b1l = {0.30, 0.30, 0.30, 0.35, 0.35, 0.35, 0.40, 0.40, 0.40};
  vector<double> fraction_4b   = {0.15, 0.15, 0.10, 0.15, 0.10, 0.10, 0.10, 0.10, 0.10};
  vector<double> fraction_3b   = {0.30, 0.25, 0.25, 0.25, 0.25, 0.20, 0.25, 0.20, 0.15};
  vector<double> fraction_2b1c = {0.25, 0.30, 0.35, 0.25, 0.30, 0.35, 0.25, 0.30, 0.35};
  for (int i=0; i<fraction_2b1l.size(); i++) {
    mc16_tag0_DL1r_mix[i] = (TH1F*)mc16_tag0_DL1r[0]->Clone();
    mc16_tag0_DL1r_mix[i]->Scale(fraction_2b1l[i]);
    mc16_tag1_DL1r_mix[i] = (TH1F*)mc16_tag1_DL1r[0]->Clone();
    mc16_tag1_DL1r_mix[i]->Scale(fraction_2b1l[i]);
    mc16_tag2_DL1r_mix[i] = (TH1F*)mc16_tag2_DL1r[0]->Clone();
    mc16_tag2_DL1r_mix[i]->Scale(fraction_2b1l[i]);

    mc16_tag0_DL1r_mix[i]->Add(mc16_tag0_DL1r[1], fraction_4b[i]);
    mc16_tag0_DL1r_mix[i]->Add(mc16_tag0_DL1r[2], fraction_3b[i]);
    mc16_tag0_DL1r_mix[i]->Add(mc16_tag0_DL1r[3], fraction_2b1c[i]);

    mc16_tag1_DL1r_mix[i]->Add(mc16_tag1_DL1r[1], fraction_4b[i]);
    mc16_tag1_DL1r_mix[i]->Add(mc16_tag1_DL1r[2], fraction_3b[i]);
    mc16_tag1_DL1r_mix[i]->Add(mc16_tag1_DL1r[3], fraction_2b1c[i]);

    mc16_tag2_DL1r_mix[i]->Add(mc16_tag2_DL1r[1], fraction_4b[i]);
    mc16_tag2_DL1r_mix[i]->Add(mc16_tag2_DL1r[2], fraction_3b[i]);
    mc16_tag2_DL1r_mix[i]->Add(mc16_tag2_DL1r[3], fraction_2b1c[i]); }


  
  // Perforn fit of the mixture with histograms
  vector<TH1*> tag2_fit_results;
  TH1 *empty_hists = new TH1F("empty_hist", "empty_hist", 30, -15, 15);
  TObjArray *tag2_4_templates = new TObjArray(4);
  for (int topHFFF_i=0; topHFFF_i<4; topHFFF_i++) {
    tag2_4_templates->Add(mc16_tag2_DL1r[topHFFF_i]); }

  // Loop over predefined fractions
  for (int i=0; i<fraction_2b1l.size(); i++) {
  
    // Create combined templates
    // (1) 3b+2b1l+2b1c, (2) 4b
    TH1 *combined_2b1l_3b_2b1c_template = new TH1F("combined_2b+3b", "combined_2b+3b", 30, -15, 15);
    combined_2b1l_3b_2b1c_template->Add(mc16_tag2_DL1r[0], fraction_2b1l[i]); // 2b1l
    combined_2b1l_3b_2b1c_template->Add(mc16_tag2_DL1r[2], fraction_3b[i]); // 3b
    combined_2b1l_3b_2b1c_template->Add(mc16_tag2_DL1r[3], fraction_2b1c[i]); // 2b1c
    combined_2b1l_3b_2b1c_template->Scale(1/combined_2b1l_3b_2b1c_template->Integral(0, combined_2b1l_3b_2b1c_template->GetNbinsX()+1));
    TObjArray *tag2_2_templates = new TObjArray(2);
    tag2_2_templates->Add(combined_2b1l_3b_2b1c_template);
    tag2_2_templates->Add(mc16_tag2_DL1r[1]);
    
    // (1) 4b+3b, (2) 2b1l, (3) 2b1c
    TH1 *combined_4b_3b_template = new TH1F("combined_extra_b", "combined_extra_b", 30, -15, 15);
    combined_4b_3b_template->Add(mc16_tag2_DL1r[1], fraction_4b[i]); // 4b
    combined_4b_3b_template->Add(mc16_tag2_DL1r[2], fraction_3b[i]); // 3b
    combined_4b_3b_template->Scale(1/combined_4b_3b_template->Integral(0, combined_4b_3b_template->GetNbinsX()+1));
    TObjArray *tag2_3_templates = new TObjArray(3);
    tag2_3_templates->Add(mc16_tag2_DL1r[0]);
    tag2_3_templates->Add(combined_4b_3b_template);
    tag2_3_templates->Add(mc16_tag2_DL1r[3]);
    
    // Collect templates and run the fit
    cout << "\n\n\n\n" << endl;
    //cout << "Performing fit [" << i << "] for\n\t2b1l frac. = " << fraction_2b1l[i] << ";\n\t4b frac. = " << fraction_4b[i] << ";\n\t3b frac. = " << fraction_3b[i] << ";\n\t2b1c frac. = " << fraction_2b1c[i] << ";" << endl;
    cout << "Performing fit [" << i << "] for\n\t2b1l frac. = " << fraction_2b1l[i] <<";\n\t3b4b frac. = " << fraction_4b[i]+fraction_3b[i] << ";\n\t2b1c frac. = " << fraction_2b1c[i] << ";" << endl;
    TObjArray *tag2_self_template = new TObjArray(1);
    tag2_self_template->Add(mc16_tag2_DL1r_mix[i]);
    tag2_self_template->Add(mc16_tag2_DL1r_mix[i]);
    //TFractionFitter *fit = new TFractionFitter(mc16_tag2_DL1r_mix[i], tag2_4_templates);
    //TFractionFitter *fit = new TFractionFitter(mc16_tag2_DL1r_mix[i], tag2_self_template);
    //TFractionFitter *fit = new TFractionFitter(mc16_tag2_DL1r_mix[i], tag2_2_templates);
    TFractionFitter *fit = new TFractionFitter(mc16_tag2_DL1r_mix[i], tag2_3_templates);
    fit->Constrain(0, 0.0, 100.0); // constrain fraction 1 to be between 0 and 100
    fit->Constrain(1, 0.0, 100.0);
    fit->Constrain(2, 0.0, 100.0);
    //fit->Constrain(3, 0.0, 100.0);
    //fit->SetRangeX(10, 25); // use bins from V1 through V2 in the fit
    Int_t status = fit->Fit(); // perform the fit
    
    if (status==0) { 
      cout << "Fit [" << i << "] success! STATUS = " << status << endl;
      tag2_fit_results.push_back( (TH1F*) fit->GetPlot() ); } 
    else { 
      cout << "Fit [" << i << "] error! STATUS = " << status << endl; 
      tag2_fit_results.push_back(empty_hists); } }
  

  
  // Draw mixtures and fits on one canvas
  for (int i=0; i<fraction_2b1l.size(); i++) {
    TString save_name = "mixrute_fit_result_" + to_string(i);
    vector<TString> fit_titles = {"mixture", "fit"};
    vector<TH1*> pair_of_hists  = {mc16_tag2_DL1r_mix[i], tag2_fit_results[i]};
    int draw_fit_results = draw_n_histos(pair_of_hists, fit_titles, "#bf{3^{rd} DL1r tag weight}", save_name, true); }
  


  // Work is done, close the hists file
  hists_file_mc->Close();

}
