// ##############################################
// ###   Draw 4 derived fractions from fits   ###
// ##############################################
int draw_funcs(vector<double> vec_2b1l, vector<double> vec_4b, vector<double> vec_3b, vector<double> vec_2b1c, vector<double> norm_factors, vector<double> initial_fractions, TString savename_ext)
{
  Int_t len = norm_factors.size();

  TLine *aim_2b1l = new TLine(norm_factors[0], initial_fractions[0], norm_factors[len-1], initial_fractions[0]);
  aim_2b1l->SetLineColor(1);
  aim_2b1l->SetLineStyle(9);
  aim_2b1l->SetLineWidth(3);
  TLine *aim_4b = new TLine(norm_factors[0], initial_fractions[1], norm_factors[len-1], initial_fractions[1]);
  aim_4b->SetLineColor(2);
  aim_4b->SetLineStyle(9);
  aim_4b->SetLineWidth(3);
  TLine *aim_3b = new TLine(norm_factors[0], initial_fractions[2], norm_factors[len-1], initial_fractions[2]);
  aim_3b->SetLineColor(3);
  aim_3b->SetLineStyle(9);
  aim_3b->SetLineWidth(3);
  TLine *aim_2b1c = new TLine(norm_factors[0], initial_fractions[3], norm_factors[len-1], initial_fractions[3]);
  aim_2b1c->SetLineColor(4);
  aim_2b1c->SetLineStyle(9);
  aim_2b1c->SetLineWidth(3);
  
  double arr_norm_factors[norm_factors.size()];
  copy(norm_factors.begin(), norm_factors.end(), arr_norm_factors);
  
  double arr_2b1l[vec_2b1l.size()];
  copy(vec_2b1l.begin(), vec_2b1l.end(), arr_2b1l);
  TGraph *gr_2b1l = new TGraph(len, arr_norm_factors, arr_2b1l);
  gr_2b1l->SetMarkerStyle(20);
  gr_2b1l->SetMarkerSize(2);
  gr_2b1l->SetLineColor(1);
  gr_2b1l->SetLineWidth(4);
  gr_2b1l->SetMarkerColor(1);
  
  double arr_4b[vec_4b.size()];
  copy(vec_4b.begin(), vec_4b.end(), arr_4b);
  TGraph *gr_4b = new TGraph(len, arr_norm_factors, arr_4b);
  gr_4b->SetMarkerStyle(20);
  gr_4b->SetMarkerSize(2);
  gr_4b->SetLineColor(2);
  gr_4b->SetLineWidth(4);
  gr_4b->SetMarkerColor(2);
  
  double arr_3b[vec_3b.size()];
  copy(vec_3b.begin(), vec_3b.end(), arr_3b);
  TGraph *gr_3b = new TGraph(len, arr_norm_factors, arr_3b);
  gr_3b->SetMarkerStyle(20);
  gr_3b->SetMarkerSize(2);
  gr_3b->SetLineColor(3);
  gr_3b->SetLineWidth(4);
  gr_3b->SetMarkerColor(3);

  double arr_2b1c[vec_2b1c.size()];
  copy(vec_2b1c.begin(), vec_2b1c.end(), arr_2b1c);
  TGraph *gr_2b1c = new TGraph(len, arr_norm_factors, arr_2b1c);
  gr_2b1c->SetMarkerStyle(20);
  gr_2b1c->SetMarkerSize(2);
  gr_2b1c->SetLineColor(4);
  gr_2b1c->SetLineWidth(4);
  gr_2b1c->SetMarkerColor(4);
  
  TCanvas *c1 = new TCanvas("fit results", "fit results", 2000, 1500);
  gStyle->SetOptStat(0);
  gPad->SetGrid();
  
  TLegend *legend = new TLegend(0.70, 0.70, 0.90, 0.90);
  legend->AddEntry(gr_2b1l, "2b1l");
  legend->AddEntry(gr_2b1l, "4b");
  legend->AddEntry(gr_2b1l, "3b");
  legend->AddEntry(gr_2b1l, "2b1c");

  gr_2b1l->Draw("ACP");
  gr_4b->Draw("CPsame");
  gr_3b->Draw("CPsame");
  gr_2b1c->Draw("CPsame");
  legend->Draw("same");
  aim_2b1l->Draw("same");
  aim_4b->Draw("same");
  aim_3b->Draw("same");
  aim_2b1c->Draw("same");

  gr_2b1l->GetYaxis()->SetRangeUser(0, 1);
  gr_2b1l->SetTitle("Fit results");
  gr_2b1l->GetXaxis()->SetTitle("#bf{DL1r tag weight, 3rd tag}");
  gr_2b1l->GetYaxis()->SetTitle("#bf{Fit fractions}");
  
  c1->Print("Plots/fit_results_plot_" + savename_ext + ".png");

  return 0;
}




// ####################
// #####   MAIN   #####
// ####################
void plot_fit_results()
{
  vector<double> initial_fractions_1= {0.60, 0.10, 0.10, 0.20};
  vector<double> norm_factors1 = {   1,    3,    5,    7,   10,   15};
  vector<double> f_2b1l  =       {0.40, 0.47, 0.54, 0.54, 0.58, 0.58};
  vector<double> f_4b    =       {0.09, 0.09, 0.09, 0.10, 0.10, 0.09};
  vector<double> f_3b    =       {0.26, 0.22, 0.17, 0.17, 0.14, 0.14};
  vector<double> f_2b1c  =       {0.24, 0.22, 0.20, 0.20, 0.18, 0.18};

  vector<double> initial_fractions_2 = {0.60, 0.15, 0.15, 0.10};
  vector<double> norm_factors2 = {   1,    3,    5,    7,   10,   15,   20,   30,   40,   50};
  vector<double> s_2b1l  =       {0.47, 0.46, 0.54, 0.54, 0.56, 0.58, 0.58, 0.60, 0.60, 0.59};
  vector<double> s_4b    =       {0.09, 0.12, 0.14, 0.14, 0.14, 0.14, 0.14, 0.15, 0.15, 0.15};
  vector<double> s_3b    =       {0.22, 0.23, 0.18, 0.18, 0.17, 0.17, 0.18, 0.16, 0.16, 0.16};
  vector<double> s_2b1c  =       {0.22, 0.18, 0.15, 0.15, 0.13, 0.11, 0.11, 0.10, 0.10, 0.10};
  
  
  int plot_s = draw_funcs(s_2b1l, s_4b, s_3b, s_2b1c, norm_factors2, initial_fractions_2, "60_15_15_10");
  int plot_f = draw_funcs(f_2b1l, f_4b, f_3b, f_2b1c, norm_factors1, initial_fractions_1, "60_10_10_20");
}
