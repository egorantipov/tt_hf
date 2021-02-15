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
#include <TLorentzVector.h>

#include <iostream>
#include <sstream>
#include <vector>

#include "KLFitter/DetectorAtlas_8TeV.h"
#include "KLFitter/Fitter.h"
#include "KLFitter/LikelihoodTopDilepton.h"
#include "KLFitter/Permutations.h"

using namespace std;



// ##################################
// ## Split string into components ##
// ##################################
vector<TString> split(TString split_string, char delimiter)
{
  stringstream ss;
  ss << split_string;
  string component;

  vector<TString> container;
  while(getline(ss, component, delimiter))
    {
      container.push_back(component);
    }

  return container;
}



// #################################################
// ## Make a list of files in the given directory ##
// #################################################
vector<TString> get_list_of_files(TString dirname, vector<TString> container = {})
{
  TSystemDirectory dir(dirname, dirname);
  TList *files = dir.GetListOfFiles();
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (fname != "." && fname != "..") {
	if (fname.EndsWith(".root")) { container.push_back(dirname + fname); }
	else { container.push_back(dirname + fname + "/"); }
      }
    }
  }
  return container;
}



// #####################
// ## Compute delta R ##
// #####################
double dR(float phi_1st, float eta_1st, float phi_2nd, float eta_2nd)
{
  double dR_val = 0;
  
  double dPhi = phi_1st - phi_2nd;
  double dEta = eta_1st - eta_2nd;
  
  if (dPhi > TMath::Pi()) dPhi = 2*TMath::Pi() - dPhi;

  dR_val = sqrt( pow(dPhi,2) + pow(dEta,2) );
  
  return dR_val;
}



// ##############
// ##   MAIN   ##
// ##############
void prepare_histograms()
{
  // Create a list of directories with ntuples
  TString path_to_ntuples = "/eos/user/e/eantipov/Files/tt_hf/";
  vector<TString> dir_paths = get_list_of_files(path_to_ntuples);

  
  // Declare histograms
  
  // dR_min between all bjets/btags and leptons
  TH1 *mc16_minDeltaR_lep0_bjets_from_top = new TH1F("minDeltaR_lep1_bjets_fromTop", "minDeltaR_lep1_bjets_fromTop", 20, 0, 5);
  TH1 *mc16_minDeltaR_lep0_bjets_not_from_top = new TH1F("minDeltaR_lep1_bjets_notFromTop", "minDeltaR_lep1_bjets_notFromTop", 20, 0, 5);
  TH1 *mc16_minDeltaR_lep1_bjets_from_top = new TH1F("minDeltaR_lep2_bjets_fromTop", "minDeltaR_lep2_bjets_fromTop", 20, 0, 5);
  TH1 *mc16_minDeltaR_lep1_bjets_not_from_top = new TH1F("minDeltaR_lep2_bjets_notFromTop", "minDeltaR_lep2_bjets_notFromTop", 20, 0, 5);
  
  // dR_min between the 1st/2nd/3rd bjets/btags and leptons
  TH1 *mc16_minDeltaR_lep0_bjets_from_top_tags[3];
  TH1 *mc16_minDeltaR_lep0_bjets_not_from_top_tags[3];
  TH1 *mc16_minDeltaR_lep1_bjets_from_top_tags[3];
  TH1 *mc16_minDeltaR_lep1_bjets_not_from_top_tags[3];
  for (int i=0; i<3; i++) {
      TString title1 = "minDeltaR_lep1_bjets_fromTop_" + to_string(i) + "_tag";
      TString title2 = "minDeltaR_lep1_bjets_notFromTop_" + to_string(i) + "_tag";
      TString title3 = "minDeltaR_lep2_bjets_fromTop" + to_string(i) + "_tag";
      TString title4 = "minDeltaR_lep2_bjets_notFromTop" + to_string(i) + "_tag";
      mc16_minDeltaR_lep0_bjets_from_top_tags[i] = new TH1F(title1, title1, 20, 0, 5);
      mc16_minDeltaR_lep0_bjets_not_from_top_tags[i] = new TH1F(title2, title2, 20, 0, 5);
      mc16_minDeltaR_lep1_bjets_from_top_tags[i] = new TH1F(title3, title3, 20, 0, 5);
      mc16_minDeltaR_lep1_bjets_not_from_top_tags[i] = new TH1F(title4, title4, 20, 0, 5); }
  
  // dR_min between bjets/btags
  TH1 *mc16_dR_b0_b1 = new TH1F("mc16_dR_b0_b1", "mc16_dR_b0_b1", 20, 0, 5);
  TH1 *mc16_dR_b0_b2 = new TH1F("mc16_dR_b0_b2", "mc16_dR_b0_b2", 20, 0, 5);
  TH1 *mc16_dR_b1_b2 = new TH1F("mc16_dR_b1_b2", "mc16_dR_b1_b2", 20, 0, 5);
  TH1 *mc16_minDeltaR_b01_b2_from_top = new TH1F("mc16_minDeltaR_b0_b2_from_top", "mc16_minDeltaR_b0_b2_from_top", 20, 0, 5);
  TH1 *mc16_minDeltaR_b01_b2_not_from_top = new TH1F("mc16_minDeltaR_b01_b2_not_from_top", "mc16_minDeltaR_b01_b2_not_from_top", 20, 0, 5);
  TH1 *mc16_minDeltaR_b_from_top_to_b = new TH1F("mc16_minDeltaR_b_from_top_to_b", "mc16_minDeltaR_b_from_top_to_b", 20, 0, 5);
  TH1 *mc16_minDeltaR_b_not_from_top_to_b = new TH1F("mc16_minDeltaR_b_not_from_top_to_b", "mc16_minDeltaR_b_not_from_top_to_b", 20, 0, 5);
  TH1 *mc16_minDeltaR_not_b_to_b = new TH1F("mc16_minDeltaR_not_b_to_b", "mc16_minDeltaR_not_b_to_b", 20, 0, 5);
  TH1 *mc16_minDeltaR_b_from_top_to_jet = new TH1F("mc16_minDeltaR_b_from_top_to_jet", "mc16_minDeltaR_b_from_top_to_jet", 20, 0, 5);
  TH1 *mc16_minDeltaR_b_not_from_top_to_jet = new TH1F("mc16_minDeltaR_b_not_from_top_to_jet", "mc16_minDeltaR_b_not_from_top_to_jet", 20, 0, 5);
  TH1 *mc16_minDeltaR_not_b_to_jet = new TH1F("mc16_minDeltaR_not_b_to_jet", "mc16_minDeltaR_not_b_to_jet", 20, 0, 5);
  TH1 *mc16_minDeltaR_b_from_top_to_lep = new TH1F("mc16_minDeltaR_b_from_top_to_lep", "mc16_minDeltaR_b_from_top_to_lep", 20, 0, 5);
  TH1 *mc16_minDeltaR_b_not_from_top_to_lep = new TH1F("mc16_minDeltaR_b_not_from_top_to_lep", "mc16_minDeltaR_b_not_from_top_to_lep", 20, 0, 5);
  TH1 *mc16_minDeltaR_not_b_to_lep = new TH1F("mc16_minDeltaR_not_b_to_lep", "mc16_minDeltaR_not_b_to_lep", 20, 0, 5);
  
  
  // pT of the first three bjets/btags (not)from top
  TH1 *mc16_jet_pT0_from_top = new TH1F("mc16_jet_pT0_from_top", "mc16_jet_pT0_from_top", 10, 0, 500);
  TH1 *mc16_jet_pT1_from_top = new TH1F("mc16_jet_pT1_from_top", "mc16_jet_pT1_from_top", 10, 0, 500);
  TH1 *mc16_jet_pT2_from_top = new TH1F("mc16_jet_pT2_from_top", "mc16_jet_pT2_from_top", 10, 0, 500);
  TH1 *mc16_jet_pT0_not_from_top = new TH1F("mc16_jet_pT0_not_from_top", "mc16_jet_pT0_not_from_top", 10, 0, 500);
  TH1 *mc16_jet_pT1_not_from_top = new TH1F("mc16_jet_pT1_not_from_top", "mc16_jet_pT1_not_from_top", 10, 0, 500);
  TH1 *mc16_jet_pT2_not_from_top = new TH1F("mc16_jet_pT2_not_from_top", "mc16_jet_pT2_not_from_top", 10, 0, 500);
 
  // the first three DL1r tags distributions for 2b1l / 4b / 3b / 2b1c
  TH1 *mc16_bjets_n_2b[4];
  TH1 *mc16_tag0_DL1r[4];
  TH1 *mc16_tag1_DL1r[4];
  TH1 *mc16_tag2_DL1r[4];
  for (int i=0; i<4; i++) {
      TString h_title0 = "mc16_tag0_DL1r_TopHFFF" + to_string(i);
      TString h_title1 = "mc16_tag1_DL1r_TopHFFF" + to_string(i);
      TString h_title2 = "mc16_tag2_DL1r_TopHFFF" + to_string(i);
      mc16_tag0_DL1r[i] = new TH1F(h_title0, h_title0, 30, -15, 15);
      mc16_tag1_DL1r[i] = new TH1F(h_title1, h_title1, 30, -15, 15);
      mc16_tag2_DL1r[i] = new TH1F(h_title2, h_title2, 30, -15, 15);

      TString h_bjn_title = "mc16_bjets_n_TopHFFF"+to_string(i);
      mc16_bjets_n_2b[i] = new TH1F(h_bjn_title, h_bjn_title, 5, 0, 5); }

  // MET
  TH1 *mc16_met = new TH1F("mc16_met", "mc16_met", 20, 0, 1000);
  TH1 *mc16_met_phi = new TH1F("mc16_met_phi", "mc16_met_phi", 40, -4, 4);

  // jet pT in high-low order
  TH1 *mc16_jet_pT0 = new TH1F("mc16_jet_pT0", "mc16_jet_pT0", 20, 0, 1000);
  TH1 *mc16_jet_pT1 = new TH1F("mc16_jet_pT1", "mc16_jet_pT1", 20, 0, 1000);
  TH1 *mc16_jet_pT2 = new TH1F("mc16_jet_pT2", "mc16_jet_pT2", 20, 0, 1000);
  TH1 *mc16_jet_pT3 = new TH1F("mc16_jet_pT3", "mc16_jet_pT3", 20, 0, 1000);
  TH1 *mc16_jet_pT4 = new TH1F("mc16_jet_pT4", "mc16_jet_pT4", 20, 0, 1000);
  TH1 *mc16_jet_pT5 = new TH1F("mc16_jet_pT5", "mc16_jet_pT5", 20, 0, 1000);
  
  // leptons
  TH1 *mc16_lep0_pt = new TH1F("mc16_lep0_pt", "mc16_lep0_pt", 20, 0, 1000);
  TH1 *mc16_lep1_pt = new TH1F("mc16_lep1_pt", "mc16_lep1_pt", 20, 0, 1000);
  TH1 *mc16_lep_pt = new TH1F("mc16_lep_pt", "mc16_lep_pt", 20, 0, 1000);
  TH1 *mc16_lep0_eta = new TH1F("mc16_lep0_eta", "mc16_let0_eta", 20, -5, 5);
  TH1 *mc16_lep1_eta = new TH1F("mc16_lep1_eta", "mc16_lep1_eta", 20, -5, 5);
  TH1 *mc16_lep_eta = new TH1F("mc16_lep_eta", "mc16_lep_eta", 20, -5, 5);
  TH1 *mc16_lep0_phi = new TH1F("mc16_lep0_phi", "mc16_lep0_phi", 40, -4, 4);
  TH1 *mc16_lep1_phi = new TH1F("mc16_lep1_phi", "mc16_lep1_phi", 40, -4, 4);
  TH1 *mc16_lep_phi = new TH1F("mc16_lep_phi", "mc16_lep_phi", 40, -4, 4);
  TH1 *mc16_dR_lep0_lep1 = new TH1F("mc16_dR_lep0_lep1", "mc16_dR_lep0_lep1", 20, 0, 5);

  // Different jets multiplicity 
  TH1 *mc16_njets_from_top = new TH1F("m16_njets_from_top", "m16_njets_from_top", 4, 0, 4);
  TH1 *mc16_nbjets_from_top = new TH1F("mc16_nbjets_from_top", "mc16_nbjets_from_top", 4, 0, 4);
  
  // Invariant mass
  TH1 *mc16_inv_mass_lep_bjet_from_top_min_dR = new TH1F("mc16_inv_mass_lep_bjet_from_top_min_dR", "mc16_inv_mass_lep_bjet_from_top_min_dR", 1000, 0, 1000);
  TH1 *mc16_inv_mass_lep_bjet_not_from_top_min_dR = new TH1F("mc16_inv_mass_lep_bjet_not_from_top_min_dR", "mc16_inv_mass_lep_bjet_not_from_top_min_dR", 1000, 0, 1000);
  TH1 *mc16_inv_mass_lep_btag_from_top_min_dR = new TH1F("mc16_inv_mass_lep_btag_from_top_min_dR", "mc16_inv_mass_lep_btag_from_top_min_dR", 1000, 0, 1000);
  TH1 *mc16_inv_mass_lep_btag_not_from_top_min_dR = new TH1F("mc16_inv_mass_lep_btag_not_from_top_min_dR", "mc16_inv_mass_lep_btag_not_from_top_min_dR", 1000, 0, 1000);
  TH1 *mc16_min_inv_mass_lep_bjet_from_top = new TH1F("mc16_min_inv_mass_lep_bjet_from_top", "mc16_min_inv_mass_lep_bjet_from_top", 1000, 0, 1000);
  TH1 *mc16_max_inv_mass_lep_bjet_from_top = new TH1F("mc16_max_inv_mass_lep_bjet_from_top", "mc16_max_inv_mass_lep_bjet_from_top", 1000, 0, 1000);
  TH1 *mc16_min_inv_mass_lep_bjet_not_from_top = new TH1F("mc16_min_inv_mass_lep_bjet_not_from_top", "mc16_min_inv_mass_lep_bjet_not_from_top", 1000, 0, 1000);
  TH1 *mc16_max_inv_mass_lep_bjet_not_from_top = new TH1F("mc16_max_inv_mass_lep_bjet_not_from_top", "mc16_max_inv_mass_lep_bjet_not_from_top", 1000, 0, 1000);
  TH1 *mc16_min_inv_mass_lep_other_jet = new TH1F("mc16_min_inv_mass_lep_other_jet", "mc16_min_inv_mass_lep_other_jet", 1000, 0, 1000);
  TH1 *mc16_max_inv_mass_lep_other_jet = new TH1F("mc16_max_inv_mass_lep_other_jet", "mc16_max_inv_mass_lep_other_jet", 1000, 0, 1000);

  // Initialize KLFitter
  KLFitter::Fitter fitter{};
  
  KLFitter::DetectorAtlas_8TeV detector{"/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/AnalysisTop/KLFitterTFs/mc15c/akt4_EMtopo_PP6"};
  fitter.SetDetector(&detector);

  KLFitter::LikelihoodTopDilepton likelihood{};
  likelihood.PhysicsConstants()->SetMassTop(172.5);
  likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kNotag);
  likelihood.SetFlagTopMassFixed(true);
  fitter.SetLikelihood(&likelihood);


  // CReate a ROOT file for NN
  vector<vector<int>> NN_tHOF_v, NN_jet_truthflav_v;
  vector<float> NN_min_dR_b_from_top_to_lep_v, NN_min_dR_b_not_from_top_to_lep_v, NN_min_dR_not_b_to_lep_v;
  vector<float> NN_M_lep_bjet_min_dR_from_top_v, NN_M_lep_bjet_min_dR_not_from_top_v;
  vector<float> NN_min_M_lep_bjet_from_top_v, NN_max_M_lep_bjet_from_top_v;
  vector<float> NN_min_M_lep_bjet_not_from_top_v, NN_max_M_lep_bjet_not_from_top_v;
  vector<float> NN_min_M_lep_other_jet_v, NN_max_M_lep_other_jet_v;


  // Loop over directories with ntuples collections
  for (int dir_counter=0; dir_counter<dir_paths.size(); dir_counter++)
    {
      // Announce current directory
      cout << "\n\n\n" << dir_paths[dir_counter] << endl;
      

      // Check for the content: data/mc? which campaign?
      vector<TString> dir_path_components = split(dir_paths[dir_counter], '/');
      int last_element_index = dir_path_components.size();
      vector<TString> dir_name_components = split(dir_path_components[last_element_index-1], '_');
      bool is_data = false;
      bool is_2015 = false;
      bool is_2016 = false;
      bool is_2017 = false;
      bool is_2018 = false;
      bool is_mc16a = false;
      bool is_mc16d = false;
      bool is_mc16e = false;
      for (int i=0; i<dir_name_components.size(); i++)
	{
	  if (dir_name_components[i] == "data") is_data = true;
          if (dir_name_components[i] == "2015") is_2015 = true;
          if (dir_name_components[i] == "2016") is_2016 = true;
          if (dir_name_components[i] == "2017") is_2017 = true;
          if (dir_name_components[i] == "2018") is_2018 = true;
          if (dir_name_components[i] == "mc16a") is_mc16a = true;
          if (dir_name_components[i] == "mc16d") is_mc16d = true;
          if (dir_name_components[i] == "mc16e") is_mc16e = true;
	}
      
      
      // We work with MC only
      if (is_data == true) continue;
      

      // Testing option: run over mc16a campaign only to save time
      //if (is_mc16a != true) continue;
      

      // Make a list of paths to jobs/DIDs outputs (pieces of a full ntuple)
      vector<TString> paths_to_jobs = get_list_of_files(dir_paths[dir_counter]);


      // Loop over jobs/DIDs
      for (int job_number=0; job_number<paths_to_jobs.size(); job_number++)
	{
	  // Get info about the job/DID from its name
	  vector<TString> path_to_jobs_components = split(paths_to_jobs[job_number], '/');
	  TString job_name = path_to_jobs_components[path_to_jobs_components.size() - 1];
	  vector<TString> job_name_components = split(job_name, '.');
	  TString job_DID = job_name_components[2];
	  vector<TString> campaign_info = split(job_name_components[5], '_');

	  
	  // Select only jobs/physics_processes of our interest: 
	  // (1) regular (not alternamtive) samples
	  // (2) tt+any, ttbb, ttb, ttc
	  if (campaign_info[1]!="s3126") continue;
	  if (job_DID!="410472" && job_DID!="411076" && job_DID!="411077" && job_DID!="411078") { continue; }
	  else { cout << "\n\nDID: " << job_DID << endl; }


	  // Testing option: keep only tt+all if true
	  bool only_410472 = false;
	  //if (job_DID=="410472") { only_410472=true; } else { continue; }


	  // Make a list of paths to ntuples of the given job/DID
	  vector<TString> paths_to_ntuples = get_list_of_files(paths_to_jobs[job_number]);


	  // Loop over ntuples of one job/DID
	  for (int ntuple_number=0; ntuple_number<paths_to_ntuples.size(); ntuple_number++)
	    {
	      // Open ntuple
	      cout << paths_to_ntuples[ntuple_number] << endl;
	      TFile *ntuple = new TFile (paths_to_ntuples[ntuple_number]);
	      TTree *tree_nominal = (TTree*)ntuple->Get("nominal");


	      // Set all the needed branches
	      vector<Float_t> *jet_pt, *jet_DL1r, *jet_eta, *jet_phi, *jet_e;
	      vector<Float_t> *el_pt, *el_eta, *el_cl_eta, *el_phi, *el_charge, *el_e;
	      vector<Float_t> *mu_pt, *mu_eta, *mu_phi, *mu_charge, *mu_e;
	      vector<int> *topHadronOriginFlag, *jet_truthflav;
	      vector<char> *jet_DL1r_77;
	      jet_pt = jet_DL1r = jet_eta = jet_phi = jet_e = 0;
	      el_pt = el_eta = el_cl_eta = el_phi = el_charge = el_e = 0;
	      mu_pt = mu_eta = mu_phi = mu_charge = mu_e = 0;
	      topHadronOriginFlag = jet_truthflav = 0;
	      jet_DL1r_77 = 0;
	      Float_t met, met_phi;
	      tree_nominal->SetBranchAddress("jet_pt", &jet_pt);
              tree_nominal->SetBranchAddress("jet_eta", &jet_eta);
              tree_nominal->SetBranchAddress("jet_phi", &jet_phi);
	      tree_nominal->SetBranchAddress("jet_e", &jet_e);
              tree_nominal->SetBranchAddress("jet_DL1r", &jet_DL1r);
	      tree_nominal->SetBranchAddress("jet_isbtagged_DL1r_77", &jet_DL1r_77);
              tree_nominal->SetBranchAddress("jet_truthflav", &jet_truthflav);
	      tree_nominal->SetBranchAddress("el_pt", &el_pt);
              tree_nominal->SetBranchAddress("el_eta", &el_eta);
	      tree_nominal->SetBranchAddress("el_cl_eta", &el_cl_eta);
              tree_nominal->SetBranchAddress("el_phi", &el_phi);
              tree_nominal->SetBranchAddress("el_charge", &el_charge);
	      tree_nominal->SetBranchAddress("el_e", &el_e);
              tree_nominal->SetBranchAddress("mu_pt", &mu_pt);
              tree_nominal->SetBranchAddress("mu_eta", &mu_eta);
              tree_nominal->SetBranchAddress("mu_phi", &mu_phi);
              tree_nominal->SetBranchAddress("mu_charge", &mu_charge);
	      tree_nominal->SetBranchAddress("mu_e", &mu_e);
              tree_nominal->SetBranchAddress("jet_GBHInit_topHadronOriginFlag", &topHadronOriginFlag); // https://gitlab.cern.ch/TTJ/Ntuple/-/blob/master/TTJNtuple/TTJNtuple/EventSaver.h#L55 
	      tree_nominal->SetBranchAddress("met_met", &met);
	      tree_nominal->SetBranchAddress("met_phi", &met_phi);
	      

	      // Weights
	      float w_mc, w_pu, w_leptonSF, w_DL1r_77, w_jvt;
              UInt_t runNumber;
              tree_nominal->SetBranchAddress("weight_mc", &w_mc);
              tree_nominal->SetBranchAddress("weight_pileup", &w_pu);
              tree_nominal->SetBranchAddress("weight_leptonSF", &w_leptonSF);
              tree_nominal->SetBranchAddress("weight_bTagSF_DL1r_77", &w_DL1r_77);
              tree_nominal->SetBranchAddress("weight_jvt", &w_jvt);
              tree_nominal->SetBranchAddress("runNumber", &runNumber);


	      // Top flavor filter flag
	      int topHFFF;
              tree_nominal->SetBranchAddress("topHeavyFlavorFilterFlag", &topHFFF);


	      // Ignore the "ReadStreamerInfo, class:string, illegal uid=-2" erro


	      // Loop over entries
	      Int_t nEntries = tree_nominal->GetEntries();
	      cout << "\tEntries = " << nEntries << endl;
	      for (int entry=0; entry<nEntries; entry++)
		{
		  // Show events counter
		  if (entry%1000==0) { cout << "\t" << entry << "\r"; cout.flush(); }
		  tree_nominal->GetEntry(entry);
		  

		  // Compute weights
		  double weight_lumi = 1;
                  double sumWeights = 1;
                  double campaign_lumi = 1;
                  double campaign_xsection = 1;
                  double campaign_genFiltEff = 1;
                  double kFactor = 1;
                  double total_lumi = 3.21956 + 32.9881 + 44.3074 + 58.4501;
		  if (runNumber==284500)
                    {
                      campaign_lumi = 3.21956 + 32.9881;
                      if (job_DID=="411076")
                        {
                          sumWeights = 3.33006*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.008814;
                          kFactor = 1.1397;
                        }
                      if (job_DID=="411077")
                        {
                          sumWeights = 3.61088*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.046655;
                          kFactor = 1.1398;
                        }
                      if (job_DID=="411078")
                        {
                          sumWeights = 3.61598*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.039503;
                          kFactor = 1.1397;
                        }
                      if (job_DID=="410472")
                        {
                          sumWeights = 5.82869*pow(10, 10);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.10547;
                          kFactor = 1.13975636159;
                        }
                    }
		  if (runNumber==300000)
                    {
                      campaign_lumi = 44.3074;
                      if (job_DID=="411076")
                        {
                          sumWeights = 4.21891*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.008814;
                          kFactor = 1.1397;
                        }
                      if (job_DID=="411077")
                        {
                          sumWeights = 4.49595*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.046655;
                          kFactor = 1.1398;
                        }
                      if (job_DID=="411078")
                        {
                          sumWeights = 4.49400*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.039503;
                          kFactor = 1.1397;
                        }
                      if (job_DID=="410472")
                        {
                          sumWeights = 7.26510*pow(10, 10);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.10547;
                          kFactor = 1.13975636159;
                        }
                    }
		  if (runNumber==310000)
                    {
                      campaign_lumi = 58.4501;
                      if (job_DID=="411076")
                        {
                          sumWeights = 5.47811*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.008814;
                          kFactor = 1.1397;
                        }
                      if (job_DID=="411077")
                        {
                          sumWeights = 5.94763*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.046655;
                          kFactor = 1.1398;
                        }
                      if (job_DID=="411078")
                        {
                          sumWeights = 5.94190*pow(10, 9);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.039503;
                          kFactor = 1.1397;
                        }
                      if (job_DID=="410472")
                        {
                          sumWeights = 1.01641*pow(10, 11);
                          campaign_xsection = 0.72977;
                          campaign_genFiltEff = 0.10547;
                          kFactor = 1.13975636159;
                        }
                    }
		  // Actual computation:
		  weight_lumi = campaign_lumi * campaign_xsection * pow(10,6) * campaign_genFiltEff * kFactor / sumWeights;
		  double weights = w_mc * w_pu * w_leptonSF * w_DL1r_77 * w_jvt * weight_lumi;
		  //double weights = w_mc * weight_lumi; // Particle level
		  //double weights = 1; // Data

		  
		  // Initiate cuts names
		  bool emu_cut = false;
                  bool OS_cut = false;
                  bool jets_n_cut = false;
                  bool bjets_n2_cut = false;
                  bool bjets_n3_cut = false;
                  bool lep_pt_cut = false;
                  bool jet_pt_cut = true;
                  bool all_cuts = false;
                  bool topHFFF_cut = false;
		  

		  // Define cuts themselves
		  if ((*el_pt).size()==1 && (*mu_pt).size()==1) emu_cut = true;
		  if ((*el_charge)[0]!=(*mu_charge)[0]) OS_cut = true;
		  int bjets_n = 0;
		  for (int i=0; i<(*jet_pt).size(); i++) { if ( int((*jet_DL1r_77)[i]==1) ) bjets_n++; }
		  if (bjets_n==3) bjets_n3_cut = true;
                  if (bjets_n>=2) bjets_n2_cut = true;
                  int jets_n = (*jet_pt).size();
                  if (jets_n >=3) jets_n_cut = true;
                  float pt_value_cut = 25000;
                  if ((*el_pt)[0]>=pt_value_cut && (*mu_pt)[0]>=pt_value_cut) lep_pt_cut = true;
                  for (int i=0; i<(*jet_pt).size(); i++) { if ((*jet_pt)[i]<pt_value_cut) jet_pt_cut = false; }
		  if ( only_410472==true || ( (topHFFF==1 && job_DID=="411076") || (topHFFF==2 && job_DID=="411077") || (topHFFF==3 && job_DID=="411078") || (topHFFF==0 && job_DID=="410472") ) ) topHFFF_cut = true;


		  // 3b, emu, OS channel: Draw min_dR
		  if (emu_cut*OS_cut*bjets_n3_cut*lep_pt_cut*jet_pt_cut*topHFFF_cut*jets_n_cut == true)
		    {
		      // Define min(dR) as some lagre value to begin with
                      double min_dR1_top = 999999.; // leading lepton
                      double min_dR2_top = 999999.; // subleading lepton
                      double min_dR1_not_top = 999999.; // leading lepton
                      double min_dR2_not_top = 999999.; // subleading lepton

		      
		      // Set iterator over b-jets
		      int bjet_iter = 0;
		      
		      
		      // Loop over all jets, and select only b-tagged in the loop
		      for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++)
			{
			  if ( int((*jet_DL1r_77)[jet_i]==1) )
			  //if ( (*topHadronOriginFlag)[i]!=-99 ) // Particle level
			    {
			      // Define initial dR's
			      double dR1 = 0;
			      double dR2 = 0;

			      // Assign dR1 to the leading lep and dR2 to the subleading 
			      if ((*mu_pt)[0]>(*el_pt)[0]) {
				dR1 = dR((*mu_phi)[0], (*mu_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]);
				dR2 = dR((*el_phi)[0], (*el_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]); }
			      else {
				dR1 = dR((*el_phi)[0], (*el_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]);
				dR2 = dR((*mu_phi)[0], (*mu_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]); }

			      // Sort into two cases: from top and not from top
			      // Fill pT histograms
			      if ((*topHadronOriginFlag)[jet_i]==4) {
				min_dR1_top = min(min_dR1_top, dR1); min_dR2_top = min(min_dR2_top, dR2);
				if (bjet_iter==0) mc16_jet_pT0_from_top->Fill((*jet_pt)[jet_i]/1000, weights);
				if (bjet_iter==1) mc16_jet_pT1_from_top->Fill((*jet_pt)[jet_i]/1000, weights);
				if (bjet_iter==2) mc16_jet_pT2_from_top->Fill((*jet_pt)[jet_i]/1000, weights); }
			      else {
				min_dR1_not_top = min(min_dR1_not_top, dR1); min_dR2_not_top = min(min_dR2_not_top, dR2);
				if (bjet_iter==0) mc16_jet_pT0_not_from_top->Fill((*jet_pt)[jet_i]/1000, weights);
				if (bjet_iter==1) mc16_jet_pT1_not_from_top->Fill((*jet_pt)[jet_i]/1000, weights);
				if (bjet_iter==2) mc16_jet_pT2_not_from_top->Fill((*jet_pt)[jet_i]/1000, weights); }


			      // Increment number of bjets by one
			      bjet_iter++;

			    } // selection of b-jets from all jets
			  
			} // [jet_i] - loop over jets, sorted wrt pT


		      // Fill min_dR histograms
		      mc16_minDeltaR_lep0_bjets_from_top->Fill(min_dR1_top, weights);
                      mc16_minDeltaR_lep0_bjets_not_from_top->Fill(min_dR1_not_top, weights);
                      mc16_minDeltaR_lep1_bjets_from_top->Fill(min_dR2_top, weights);
                      mc16_minDeltaR_lep1_bjets_not_from_top->Fill(min_dR2_not_top, weights);
		      
		    } // 3b, emu, OS cuts



		  // 2+b, emu, OS channel
		  if (emu_cut*OS_cut*bjets_n2_cut*lep_pt_cut*jet_pt_cut*topHFFF_cut*jets_n_cut == true)
		    {
		      // MET plots:
		      mc16_met->Fill(met, 0.001);
		      mc16_met_phi->Fill(met_phi);


		      // jet pT plots in pT order:
		      mc16_jet_pT0->Fill((*jet_pt)[0], 0.001);
                      mc16_jet_pT1->Fill((*jet_pt)[1], 0.001);
                      mc16_jet_pT2->Fill((*jet_pt)[2], 0.001);
                      if ((*jet_pt).size()>=4) mc16_jet_pT3->Fill((*jet_pt)[3], 0.001);
                      if ((*jet_pt).size()>=5) mc16_jet_pT4->Fill((*jet_pt)[4], 0.001);
                      if ((*jet_pt).size()>=6) mc16_jet_pT5->Fill((*jet_pt)[5], 0.001);
		      

		      // leptons kinematic hists
		      if ( (*el_pt)[0] > (*mu_pt)[0] ) {		
			mc16_lep0_pt->Fill((*el_pt)[0], 0.001);
			mc16_lep1_pt->Fill((*mu_pt)[0], 0.001);
			mc16_lep0_eta->Fill((*el_eta)[0]);
			mc16_lep1_eta->Fill((*mu_eta)[0]);
			mc16_lep0_phi->Fill((*el_phi)[0]);
			mc16_lep1_phi->Fill((*mu_phi)[0]); }
		      else {
			mc16_lep0_pt->Fill((*mu_pt)[0], 0.001);
                        mc16_lep1_pt->Fill((*el_pt)[0], 0.001);
                        mc16_lep0_eta->Fill((*mu_eta)[0]);
                        mc16_lep1_eta->Fill((*el_eta)[0]);
                        mc16_lep0_phi->Fill((*mu_phi)[0]);
                        mc16_lep1_phi->Fill((*el_phi)[0]); }
		      
		      mc16_lep_pt->Fill((*el_pt)[0], 0.001);
		      mc16_lep_pt->Fill((*mu_pt)[0], 0.001);
		      mc16_lep_eta->Fill((*el_eta)[0]);
		      mc16_lep_eta->Fill((*mu_eta)[0]);
		      mc16_lep_phi->Fill((*el_phi)[0]);
		      mc16_lep_phi->Fill((*mu_phi)[0]);
		      
		      
		      // dR_lep0_lep1 hist:
		      double dR_lep0_lep1 = dR((*el_phi)[0], (*el_eta)[0], (*mu_phi)[0], (*mu_eta)[0]);
		      mc16_dR_lep0_lep1->Fill(dR_lep0_lep1);
		      

		      // Declare a value for n jets from top
		      int njets_from_top = 0;
		      int nbjets_from_top = 0;
		      

		      // Compute min dR for different jet-obj combinations
		      double min_dR_b_from_top_to_b = 999999.;
		      double min_dR_b_not_from_top_to_b = 999999.;
		      double min_dR_not_b_to_b = 999999.;
		      double min_dR_b_from_top_to_jet = 999999.;
		      double min_dR_b_not_from_top_to_jet = 999999.;
		      double min_dR_not_b_to_jet = 999999.;
		      double min_dR_b_from_top_to_lep = 999999.;
		      double min_dR_b_not_from_top_to_lep = 999999.;
		      double min_dR_not_b_to_lep = 999999.;
		      
		      for (int i=0; i<(*jet_pt).size(); i++) {
			
			for (int j=0; j<(*jet_pt).size(); j++) {
			  if (i==j) continue;
			  
			  // dR 
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]==4 && (*jet_truthflav)[j]==5) {
			    double dR_b_from_top_to_b = dR((*jet_phi)[i], (*jet_eta)[i], (*jet_phi)[j], (*jet_eta)[j]);
			    if (dR_b_from_top_to_b < min_dR_b_from_top_to_b) min_dR_b_from_top_to_b = dR_b_from_top_to_b; }
			  
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]!=4 && (*jet_truthflav)[j]==5) {
			    double dR_b_not_from_top_to_b = dR((*jet_phi)[i], (*jet_eta)[i], (*jet_phi)[j], (*jet_eta)[j]);
			    if (dR_b_not_from_top_to_b < min_dR_b_not_from_top_to_b) min_dR_b_not_from_top_to_b = dR_b_not_from_top_to_b; }
			  
			  if ((*jet_truthflav)[i]!=5 && (*topHadronOriginFlag)[i]!=4 && (*jet_truthflav)[j]==5) {
			    double dR_not_b_to_b = dR((*jet_phi)[i], (*jet_eta)[i], (*jet_phi)[j], (*jet_eta)[j]);
			    if (dR_not_b_to_b < min_dR_not_b_to_b) min_dR_not_b_to_b = dR_not_b_to_b; }
			  
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]==4) {
			    double dR_b_from_top_to_jet = dR((*jet_phi)[i], (*jet_eta)[i], (*jet_phi)[j], (*jet_eta)[j]);
			    if (dR_b_from_top_to_jet < min_dR_b_from_top_to_jet) min_dR_b_from_top_to_jet = dR_b_from_top_to_jet; }
			  
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]!=4) {
			    double dR_b_not_from_top_to_jet = dR((*jet_phi)[i], (*jet_eta)[i], (*jet_phi)[j], (*jet_eta)[j]);
			    if (dR_b_not_from_top_to_jet < min_dR_b_not_from_top_to_jet) min_dR_b_not_from_top_to_jet = dR_b_not_from_top_to_jet; }
			  
			  if ((*jet_truthflav)[i]!=5 && (*topHadronOriginFlag)[i]!=4) {
			    double dR_not_b_to_jet = dR((*jet_phi)[i], (*jet_eta)[i], (*jet_phi)[j], (*jet_eta)[j]);
			    if (dR_not_b_to_jet < min_dR_not_b_to_jet) min_dR_not_b_to_jet = dR_not_b_to_jet; }
			} // loop over jet[j]
			
			if ((*jet_truthflav)[i]==1 && (*topHadronOriginFlag)[i]==4) {
			  double dR_b_from_top_to_el = dR((*jet_phi)[i], (*jet_eta)[i], (*el_phi)[0], (*el_eta)[0]);
			  double dR_b_from_top_to_mu = dR((*jet_phi)[i], (*jet_eta)[i], (*mu_phi)[0], (*mu_eta)[0]);
			  double dR_b_from_top_to_lep = min(dR_b_from_top_to_el, dR_b_from_top_to_mu);
			  if (dR_b_from_top_to_lep < min_dR_b_from_top_to_lep) min_dR_b_from_top_to_lep = dR_b_from_top_to_lep; }
			
			if ((*jet_truthflav)[i]==1 && (*topHadronOriginFlag)[i]!=4) {
			  double dR_b_not_from_top_to_el = dR((*jet_phi)[i], (*jet_eta)[i], (*el_phi)[0], (*el_eta)[0]);
			  double dR_b_not_from_top_to_mu = dR((*jet_phi)[i], (*jet_eta)[i], (*mu_phi)[0], (*mu_eta)[0]);
			  double dR_b_not_from_top_to_lep = min(dR_b_not_from_top_to_el, dR_b_not_from_top_to_mu);
			  if (dR_b_not_from_top_to_lep < min_dR_b_not_from_top_to_lep) min_dR_b_not_from_top_to_lep = dR_b_not_from_top_to_lep; }
			
			if ((*jet_truthflav)[i]!=1 && (*topHadronOriginFlag)[i]!=4) {
			  double dR_not_b_to_el = dR((*jet_phi)[i], (*jet_eta)[i], (*el_phi)[0], (*el_eta)[0]);
			  double dR_not_b_to_mu = dR((*jet_phi)[i], (*jet_eta)[i], (*mu_phi)[0], (*mu_eta)[0]);
			  double dR_not_b_to_lep = min(dR_not_b_to_el, dR_not_b_to_mu);
			  if (dR_not_b_to_lep < min_dR_not_b_to_lep) min_dR_not_b_to_lep = dR_not_b_to_lep; }

			// Count number of jets from top
			if ((*topHadronOriginFlag)[i]==4) njets_from_top++;
			if ((*topHadronOriginFlag)[i]==4 && (*jet_truthflav)[i]==5) nbjets_from_top++;
			
		      } // loop over jet[i]

		      mc16_njets_from_top->Fill(njets_from_top, weights);
		      mc16_nbjets_from_top->Fill(nbjets_from_top, weights);
		      
		      mc16_minDeltaR_b_from_top_to_b->Fill(min_dR_b_from_top_to_b, weights);
		      mc16_minDeltaR_b_not_from_top_to_b->Fill(min_dR_b_not_from_top_to_b, weights);
		      mc16_minDeltaR_not_b_to_b->Fill(min_dR_not_b_to_b, weights);
		      mc16_minDeltaR_b_from_top_to_jet->Fill(min_dR_b_from_top_to_jet, weights);
		      mc16_minDeltaR_b_not_from_top_to_jet->Fill(min_dR_b_not_from_top_to_jet, weights);
		      mc16_minDeltaR_not_b_to_jet->Fill(min_dR_not_b_to_jet, weights);
		      mc16_minDeltaR_b_from_top_to_lep->Fill(min_dR_b_from_top_to_lep, weights);
		      mc16_minDeltaR_b_not_from_top_to_lep->Fill(min_dR_b_not_from_top_to_lep, weights);
		      mc16_minDeltaR_not_b_to_lep->Fill(min_dR_not_b_to_lep, weights);


		      // KLFitter invariant mass calculations:
		      //KLFitter::Particles particles{};
                      //likelihood.SetLeptonType(KLFitter::LikelihoodTopDilepton::kElectron, KLFitter::LikelihoodTopDilepton::kMuon);
		      // Add leptons
                      TLorentzVector el_lvec;
                      TLorentzVector mu_lvec;
                      el_lvec.SetPtEtaPhiE((*el_pt)[0]*0.001, (*el_eta)[0], (*el_phi)[0], (*el_e)[0]*0.001);
                      mu_lvec.SetPtEtaPhiE((*mu_pt)[0]*0.001, (*mu_eta)[0], (*mu_phi)[0], (*mu_e)[0]*0.001);
                      //particles.AddParticle(el_lvec, (*el_cl_eta)[0], (*el_charge)[0], KLFitter::Particles::kElectron);
		      //particles.AddParticle(mu_lvec, mu_lvec.Eta(), (*mu_charge)[0], KLFitter::Particles::kMuon);
                      // Add two leading pT jets
		      //for (int jet_i=0; jet_i<3; jet_i++) {
		      //TLorentzVector jet_lvec;
		      //jet_lvec.SetPtEtaPhiE((*jet_pt)[jet_i]*0.001, (*jet_eta)[jet_i], (*jet_phi)[jet_i], (*jet_e)[jet_i]*0.001); 
		      //particles.AddParticle(jet_lvec, jet_lvec.Eta(), KLFitter::Particles::kParton, "", jet_i); 
		      //}
		      //fitter.SetParticles(&particles);
		      // Add MET
		      //fitter.SetET_miss_XY_SumET(met*0.001*cos(met_phi), met*0.001*sin(met_phi), met*0.001);


		      //TEST
		      //cout << "\n" << entry << endl;
                      //cout << el_lvec.Pt() << " , " << el_lvec.Eta() << " , " << el_lvec.Phi() << " , " << el_lvec.E() << endl;
                      //cout << mu_lvec.Pt() << " , " << mu_lvec.Eta() << " , " << mu_lvec.Phi() << " , " << mu_lvec.E() << endl;
		      
		      // Loop over permutations
		      //int n_perm = fitter.Permutations()->NPermutations();
		      //for (int perm_i=0; perm_i < n_perm; perm_i++) {
		      //fitter.Fit(perm_i);
		      //auto permutedParticles = fitter.Likelihood()->PParticlesPermuted();
		      //double llh = fitter.Likelihood()->LogLikelihood(fitter.Likelihood()->GetBestFitParameters()); }
		      
		      

		      // Invariant mass of bjets/btags and leptons
		      vector<TLorentzVector> bjet_from_top_lvec;
		      vector<TLorentzVector> bjet_not_from_top_lvec;
		      vector<TLorentzVector> btag_from_top_lvec;
		      vector<TLorentzVector> btag_not_from_top_lvec;
		      
		      double inv_mass_lep_bjet_from_top_min_dR = 0;
		      double inv_mass_lep_bjet_not_from_top_min_dR = 0;
		      double min_inv_mass_lep_bjet_from_top = 999999;
		      double max_inv_mass_lep_bjet_from_top = 0;
		      double min_inv_mass_lep_bjet_not_from_top = 999999;
		      double max_inv_mass_lep_bjet_not_from_top = 0;
		      double min_inv_mass_lep_other_jet = 999999;
		      double max_inv_mass_lep_other_jet = 0;
		      
		      for (int jet_i=0; jet_i<jet_pt->size(); jet_i++) {
			TLorentzVector lvec;
			lvec.SetPtEtaPhiE((*jet_pt)[jet_i]*0.001, (*jet_eta)[jet_i], (*jet_phi)[jet_i], (*jet_e)[jet_i]*0.001);
			
			// bjets and the clesest leptons
			if ((*jet_truthflav)[jet_i]==5) {
			  if ((*topHadronOriginFlag)[jet_i]==4) {
			    bjet_from_top_lvec.push_back(lvec);
			    double dr_j_el = lvec.DeltaR(el_lvec);
			    double dr_j_mu = lvec.DeltaR(mu_lvec);
			    double inv_mass_j_lep = 0;
			    if (dr_j_el<= dr_j_mu) { inv_mass_j_lep = (lvec + el_lvec).M(); }
			    else { inv_mass_j_lep = (lvec + mu_lvec).M(); }
			    inv_mass_lep_bjet_from_top_min_dR = inv_mass_j_lep;
			    mc16_inv_mass_lep_bjet_from_top_min_dR->Fill(inv_mass_j_lep, weights); }
			  else { 
			    bjet_not_from_top_lvec.push_back(lvec);
			    double dr_j_el = lvec.DeltaR(el_lvec);
			    double dr_j_mu = lvec.DeltaR(mu_lvec);
			    double inv_mass_j_lep = 0;
			    if (dr_j_el<= dr_j_mu) { inv_mass_j_lep = (lvec + el_lvec).M(); }
			    else { inv_mass_j_lep = (lvec + mu_lvec).M(); }
			    inv_mass_lep_bjet_not_from_top_min_dR = inv_mass_j_lep;
			    mc16_inv_mass_lep_bjet_not_from_top_min_dR->Fill(inv_mass_j_lep, weights); } }
			
			// btags and the closest leptons
			if ((*jet_DL1r_77)[jet_i]==1) {
			  if ((*topHadronOriginFlag)[jet_i]==4) { 
			    btag_from_top_lvec.push_back(lvec);
			    double dr_j_el = lvec.DeltaR(el_lvec);
			    double dr_j_mu = lvec.DeltaR(mu_lvec);
			    double inv_mass_j_lep = 0;
			    if (dr_j_el<= dr_j_mu) { inv_mass_j_lep = (lvec + el_lvec).M(); }
			    else { inv_mass_j_lep = (lvec + mu_lvec).M(); }
			    mc16_inv_mass_lep_btag_from_top_min_dR->Fill(inv_mass_j_lep, weights); }
			  else { 
			    btag_not_from_top_lvec.push_back(lvec);
			    double dr_j_el = lvec.DeltaR(el_lvec);
			    double dr_j_mu = lvec.DeltaR(mu_lvec);
			    double inv_mass_j_lep = 0;
			    if (dr_j_el<= dr_j_mu) { inv_mass_j_lep = (lvec + el_lvec).M(); }
			    else { inv_mass_j_lep = (lvec + mu_lvec).M(); }
			    mc16_inv_mass_lep_btag_not_from_top_min_dR->Fill(inv_mass_j_lep, weights); } }
			
			// bjets and leptons - min and max invarinat masses
			if ((*jet_truthflav)[jet_i]==5) {
			  if ((*topHadronOriginFlag)[jet_i]==4) {
			    double min_inv_mass_lep_bjet_from_top_tmp = min( (lvec + el_lvec).M(), (lvec + mu_lvec).M() );
			    double max_inv_mass_lep_bjet_from_top_tmp = max( (lvec + mu_lvec).M(), (lvec + mu_lvec).M() );
			    min_inv_mass_lep_bjet_from_top = min(min_inv_mass_lep_bjet_from_top_tmp, min_inv_mass_lep_bjet_from_top);
			    max_inv_mass_lep_bjet_from_top = max(max_inv_mass_lep_bjet_from_top_tmp, max_inv_mass_lep_bjet_from_top_tmp); }
			  else {
			    double min_inv_mass_lep_bjet_not_from_top_tmp = min( (lvec + el_lvec).M(), (lvec + mu_lvec).M() );
			    double max_inv_mass_lep_bjet_not_from_top_tmp = max( (lvec + el_lvec).M(), (lvec + mu_lvec).M() );
			    min_inv_mass_lep_bjet_not_from_top = min(min_inv_mass_lep_bjet_not_from_top_tmp, min_inv_mass_lep_bjet_not_from_top);
			    max_inv_mass_lep_bjet_not_from_top = max(max_inv_mass_lep_bjet_not_from_top_tmp, max_inv_mass_lep_bjet_not_from_top); } }
			// other than bjets and leptons - min and max invariant masses
			else {
			  double min_inv_mass_lep_other_jet_tmp = min( (lvec + el_lvec).M(), (lvec + mu_lvec).M() );
			  double max_inv_mass_lep_other_jet_tmp = max( (lvec + el_lvec).M(), (lvec + mu_lvec).M() );
			  min_inv_mass_lep_other_jet = min( (lvec + el_lvec).M(), (lvec + mu_lvec).M() );
			  max_inv_mass_lep_other_jet = max( (lvec + el_lvec).M(), (lvec + mu_lvec).M() ); }
		      }
		      
		      // Fill the min/max invariant mass hists
		      if (min_inv_mass_lep_bjet_from_top!=0) mc16_min_inv_mass_lep_bjet_from_top->Fill(min_inv_mass_lep_bjet_from_top, weights);
		      if (max_inv_mass_lep_bjet_from_top!=999999) mc16_max_inv_mass_lep_bjet_from_top->Fill(max_inv_mass_lep_bjet_from_top, weights);
		      if (min_inv_mass_lep_bjet_not_from_top!=0) mc16_min_inv_mass_lep_bjet_not_from_top->Fill(min_inv_mass_lep_bjet_not_from_top, weights);
		      if (max_inv_mass_lep_bjet_not_from_top!=999999) mc16_max_inv_mass_lep_bjet_not_from_top->Fill(max_inv_mass_lep_bjet_not_from_top, weights);
		      if (min_inv_mass_lep_other_jet!=0) mc16_min_inv_mass_lep_other_jet->Fill(min_inv_mass_lep_other_jet, weights);
		      if (max_inv_mass_lep_other_jet!=999999) mc16_max_inv_mass_lep_other_jet->Fill(max_inv_mass_lep_other_jet, weights);


		      // Fill the NN vector variables
		      NN_min_dR_b_from_top_to_lep_v.push_back(min_dR_b_from_top_to_lep);
		      NN_min_dR_b_not_from_top_to_lep_v.push_back(min_dR_b_not_from_top_to_lep);
		      NN_min_dR_not_b_to_lep_v.push_back(min_dR_not_b_to_lep);
		      NN_M_lep_bjet_min_dR_from_top_v.push_back(inv_mass_lep_bjet_from_top_min_dR);
		      NN_M_lep_bjet_min_dR_not_from_top_v.push_back(inv_mass_lep_bjet_not_from_top_min_dR);
		      NN_min_M_lep_bjet_from_top_v.push_back(min_inv_mass_lep_bjet_from_top);
		      NN_max_M_lep_bjet_from_top_v.push_back(max_inv_mass_lep_bjet_from_top);
		      NN_min_M_lep_bjet_not_from_top_v.push_back(min_inv_mass_lep_bjet_not_from_top);
		      NN_max_M_lep_bjet_not_from_top_v.push_back(max_inv_mass_lep_bjet_not_from_top);
		      NN_min_M_lep_other_jet_v.push_back(min_inv_mass_lep_other_jet);
		      NN_max_M_lep_other_jet_v.push_back(max_inv_mass_lep_other_jet);
		      NN_tHOF_v.push_back((*topHadronOriginFlag));
		      NN_jet_truthflav_v.push_back((*jet_truthflav));	
		      
		      

		      // Sort jets wrt DL1r tag weights
		      sort (jet_DL1r->begin(), jet_DL1r->end(), greater<int>());


		      // Fill DL1r tag weight histos for the first three tags, also sort wrt topHFFF
		      mc16_tag0_DL1r[topHFFF]->Fill((*jet_DL1r)[0], weights);
                      mc16_tag1_DL1r[topHFFF]->Fill((*jet_DL1r)[1], weights);
                      mc16_tag2_DL1r[topHFFF]->Fill((*jet_DL1r)[2], weights);


		      // And fill the bjets_n histogram (need it for testing purposes)
		      mc16_bjets_n_2b[topHFFF]->Fill(bjets_n, weights);

		      
		      // min_dR for the first three tags
		      // Define min(dR) as some lagre value to begin with
		      double min_dR0_top = 999999.; // leading lepton 
		      double min_dR1_top = 999999.; // subleading lepton
		      double min_dR0_not_top = 999999.; // leading lepton
		      double min_dR1_not_top = 999999.; // subleading lepton
		      

		      // Indeces for the b_jets beign selected for min_dR value
		      int min_dR0_top_idx = 0;
		      int min_dR1_top_idx = 0;
		      int min_dR0_not_top_idx = 0;
		      int min_dR1_not_top_idx = 0;
		      
		      // Loop over jets to compute min_dR
		      for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++)
			{
			  // Compute dR and assign to the corresponding leptons
			  double dR0 = 0;
			  double dR1 = 0;
			  if ((*mu_pt)[0]>(*el_pt)[0]) {
			    dR0 = dR((*mu_phi)[0], (*mu_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]);
			    dR1 = dR((*el_phi)[0], (*el_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]); }
			  else {
			    dR0 = dR((*el_phi)[0], (*el_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]);
			    dR1 = dR((*mu_phi)[0], (*mu_eta)[0], (*jet_phi)[jet_i], (*jet_eta)[jet_i]); }


			  // Sort into two cases: from top and not from top
			  if ((*topHadronOriginFlag)[jet_i]==4) {
			    double min_dR0_top_old = min_dR0_top; double min_dR1_top_old = min_dR1_top;
			    min_dR0_top = min(min_dR0_top, dR0); min_dR1_top = min(min_dR1_top, dR1); 
			    if (min_dR0_top_old!=min_dR0_top) min_dR0_top_idx = jet_i;
			    if (min_dR1_top_old!=min_dR1_top) min_dR1_top_idx = jet_i; }
			  else {
			    double min_dR0_not_top_old = min_dR0_not_top; double min_dR1_not_top_old = min_dR1_not_top;
			    min_dR0_not_top = min(min_dR0_not_top, dR0); min_dR1_not_top = min(min_dR1_not_top, dR1); 
			    if (min_dR0_not_top_old!=min_dR0_not_top) min_dR0_not_top_idx = jet_i; 
			    if (min_dR1_not_top_old!=min_dR1_not_top) min_dR1_not_top_idx = jet_i; }
       
			} //[jet_i] - loop over jets t ocompute min_dR


		      // Fill the min_dR histos for the first three tag weights only
		      if (min_dR0_top_idx<3) mc16_minDeltaR_lep0_bjets_from_top_tags[min_dR0_top_idx]->Fill(min_dR0_top, weights);
		      if (min_dR0_not_top_idx<3) mc16_minDeltaR_lep0_bjets_not_from_top_tags[min_dR0_not_top_idx]->Fill(min_dR0_not_top, weights);
		      if (min_dR1_top_idx<3) mc16_minDeltaR_lep1_bjets_from_top_tags[min_dR1_top_idx]->Fill(min_dR1_top, weights);
		      if (min_dR1_not_top_idx<3) mc16_minDeltaR_lep1_bjets_not_from_top_tags[min_dR1_not_top_idx]->Fill(min_dR1_not_top, weights); 

		      
		      // Initiate variables to compute dR between bjets/btags:
		      double dR_b0_b1 = 0;
		      double dR_b0_b2 = 0;
		      double dR_b1_b2 = 0;
		      double min_dR_b01_b2_from_top = 0;
		      double min_dR_b01_b2_not_from_top = 0;
		      
		      // Compute dR between bjets/btags and fill hists:
		      dR_b0_b1 = dR((*jet_phi)[0], (*jet_eta)[0], (*jet_phi)[1], (*jet_eta)[1]);
		      dR_b0_b2 = dR((*jet_phi)[0], (*jet_eta)[0], (*jet_phi)[2], (*jet_eta)[2]);
		      dR_b1_b2 = dR((*jet_phi)[1], (*jet_eta)[1], (*jet_phi)[2], (*jet_eta)[2]);
		      mc16_dR_b0_b1->Fill(dR_b0_b1);
		      mc16_dR_b0_b2->Fill(dR_b0_b2);
		      mc16_dR_b1_b2->Fill(dR_b1_b2);
		      
		      // Assuming the 3rd btag is additional:
		      if ((*topHadronOriginFlag)[2]==4) { mc16_minDeltaR_b01_b2_from_top->Fill(min(dR_b0_b2, dR_b1_b2)); }
		      else if ((*topHadronOriginFlag)[2]!=-99) { mc16_minDeltaR_b01_b2_not_from_top->Fill(min(dR_b0_b2, dR_b1_b2)); }


		    } // 2+b, emu, OS cuts 

		} // [entry] - loop over entries
	      
	      // Close ntuple after we're done with it
	      ntuple->Close();
	      
	    } // [ntuple_number] - loop over ntuples of a particular job
	  
	} // [job_number] - loop over jobs (pieces) of a collection.

    } // [dir_counter] - loop over directories names with jobs folders: mc16a, mc16d, mc16e, data



  // Save histograms
  TFile *hists_file = new TFile("hists_mc.root", "RECREATE");
  for (int topHFFF_i=0; topHFFF_i<4; topHFFF_i++) {
    TString process = "";
    if (topHFFF_i==0) process = "2b1l";
    if (topHFFF_i==1) process = "4b";
    if (topHFFF_i==2) process = "3b";
    if (topHFFF_i==3) process = "2b1c";
    mc16_tag0_DL1r[topHFFF_i]->Write("DL1r_templates_"+process+"_1st_tag");
    mc16_tag1_DL1r[topHFFF_i]->Write("DL1r_templates_"+process+"_2nd_tag");
    mc16_tag2_DL1r[topHFFF_i]->Write("DL1r_templates_"+process+"_3rd_tag");
    mc16_bjets_n_2b[topHFFF_i]->Write("2b_emu_OS_"+process+"bjets_n"); }
  
  mc16_jet_pT0_from_top->Write("3b_emu_OS_jet_pT0_from_top");
  mc16_jet_pT1_from_top->Write("3b_emu_OS_jet_pT1_from_top");
  mc16_jet_pT2_from_top->Write("3b_emu_OS_jet_pT2_from_top");
  mc16_jet_pT0_not_from_top->Write("3b_emu_OS_jet_pT0_not_from_top");
  mc16_jet_pT1_not_from_top->Write("3b_emu_OS_jet_pT0_not_from_top");
  mc16_jet_pT2_not_from_top->Write("3b_emu_OS_jet_pT0_not_from_top");

  mc16_minDeltaR_lep0_bjets_from_top->Write("3b_emu_OS_min_dR_lep0_b_from_top");
  mc16_minDeltaR_lep1_bjets_from_top->Write("3b_emu_OS_min_dR_lep1_b_from_top");
  mc16_minDeltaR_lep0_bjets_not_from_top->Write("3b_emu_OS_min_dR_lep0_b_not_from_top");
  mc16_minDeltaR_lep1_bjets_not_from_top->Write("3b_emu_OS_min_dR_lep1_b_not_from_top");
      
  for (int i=0; i<3; i++) {
    TString savename1 = "2b_emu_OS_min_dR_lep0_b_from_top_" + to_string(i+1) + "tag";
    TString savename2 = "2b_emu_OS_min_dR_lep0_b_not_from_top_" + to_string(i+1) + "tag";
    TString savename3 = "2b_emu_OS_min_dR_lep1_b_from_top_" + to_string(i+1) + "tag";
    TString savename4 = "2b_emu_OS_min_dR_lep1_b_not_from_top_" + to_string(i+1) + "tag";
    mc16_minDeltaR_lep0_bjets_from_top_tags[i]->Write(savename1);
    mc16_minDeltaR_lep0_bjets_not_from_top_tags[i]->Write(savename2);
    mc16_minDeltaR_lep1_bjets_from_top_tags[i]->Write(savename3);
    mc16_minDeltaR_lep1_bjets_not_from_top_tags[i]->Write(savename4); }

  mc16_jet_pT0->Write("2b_emu_OS_jet_pt0");
  mc16_jet_pT1->Write("2b_emu_OS_jet_pt1");
  mc16_jet_pT2->Write("2b_emu_OS_jet_pt2");
  mc16_jet_pT3->Write("2b_emu_OS_jet_pt3");
  mc16_jet_pT4->Write("2b_emu_OS_jet_pt4");
  mc16_jet_pT5->Write("2b_emu_OS_jet_pt5");
  
  mc16_dR_b0_b1->Write("2b_emu_OS_dR_b0_b1");
  mc16_dR_b0_b2->Write("2b_emu_OS_dR_b0_b2");
  mc16_dR_b1_b2->Write("2b_emu_OS_dR_b1_b2");
  mc16_minDeltaR_b01_b2_from_top->Write("2b_emu_OS_min_dR_b01_b2_from_top");
  mc16_minDeltaR_b01_b2_not_from_top->Write("2b_emu_OS_min_dR_b01_b2_not_from_top");
  
  mc16_met->Write("2b_emu_OS_met");
  mc16_met_phi->Write("2b_emu_OS_met_phi");
  
  mc16_lep0_pt->Write("2b_emu_OS_lep0_pt");
  mc16_lep1_pt->Write("2b_emu_OS_lep1_pt");
  mc16_lep_pt->Write("2b_emu_OS_lep_pt");
  mc16_lep0_eta->Write("2b_emu_OS_lep0_eta");
  mc16_lep1_eta->Write("2b_emu_OS_lep1_eta");
  mc16_lep_eta->Write("2b_emu_OS_lep_eta");
  mc16_lep0_phi->Write("2b_emu_OS_lep0_phi");
  mc16_lep1_phi->Write("2b_emu_OS_lep1_phi");
  mc16_lep_phi->Write("2b_emu_OS_lep_phi");
  mc16_dR_lep0_lep1->Write("2b_emu_OS_dR_lep0_lep1");
  
  mc16_minDeltaR_b_from_top_to_b->Write("2b_emu_OS_mc16_minDeltaR_b_from_top_to_b");
  mc16_minDeltaR_b_not_from_top_to_b->Write("2b_emu_OS_mc16_minDeltaR_b_not_from_top_to_b");
  mc16_minDeltaR_not_b_to_b->Write("2b_emu_OS_mc16_minDeltaR_not_b_to_b");
  mc16_minDeltaR_b_from_top_to_jet->Write("2b_emu_OS_mc16_minDeltaR_b_from_top_to_jet");
  mc16_minDeltaR_b_not_from_top_to_jet->Write("2b_emu_OS_mc16_minDeltaR_b_not_from_top_to_jet");
  mc16_minDeltaR_not_b_to_jet->Write("2b_emu_OS_mc16_minDeltaR_not_b_to_jet");
  mc16_minDeltaR_b_from_top_to_lep->Write("2b_emu_OS_mc16_minDeltaR_b_from_top_to_lep");
  mc16_minDeltaR_b_not_from_top_to_lep->Write("2b_emu_OS_mc16_minDeltaR_b_not_from_top_to_lep");
  mc16_minDeltaR_not_b_to_lep->Write("2b_emu_OS_mc16_minDeltaR_not_b_to_lep");
  
  mc16_njets_from_top->Write("2b_emu_OS_njets_from_top");
  mc16_njets_from_top->Write("2b_emu_OS_nbjets_from_top");
  
  // Invariant mass
  mc16_inv_mass_lep_bjet_from_top_min_dR->Write("2b_emu_OS_mc16_inv_mass_lep_bjet_from_top_min_dR");
  mc16_inv_mass_lep_bjet_not_from_top_min_dR->Write("2b_emu_OS_mc16_inv_mass_lep_bjet_not_from_top_min_dR");
  mc16_inv_mass_lep_btag_from_top_min_dR->Write("2b_emu_OS_mc16_inv_mass_lep_btag_from_top_min_dR");
  mc16_inv_mass_lep_btag_not_from_top_min_dR->Write("2b_emu_OS_mc16_inv_mass_lep_btag_not_from_top_min_dR");
  mc16_min_inv_mass_lep_bjet_from_top->Write("2b_emu_OS_mc16_min_inv_mass_lep_bjet_from_top");
  mc16_max_inv_mass_lep_bjet_from_top->Write("2b_emu_OS_mc16_max_inv_mass_lep_bjet_from_top");
  mc16_min_inv_mass_lep_bjet_not_from_top->Write("2b_emu_OS_mc16_min_inv_mass_lep_bjet_not_from_top");
  mc16_max_inv_mass_lep_bjet_not_from_top->Write("2b_emu_OS_mc16_max_inv_mass_lep_bjet_not_from_top");
  mc16_min_inv_mass_lep_other_jet->Write("2b_emu_OS_mc16_min_inv_mass_lep_other_jet");
  mc16_max_inv_mass_lep_other_jet->Write("2b_emu_OS_mc16_max_inv_mass_lep_other_jet");

  // Close the hists file
  hists_file->Close();


  // Fill NN ROOT file
  TFile *NN_tfile= new TFile("tt_jets_NN_input.root", "RECREATE");
  TTree *NN_ttree = new TTree("nominal", "NN_input");
  vector<int> *NN_tHOF, *NN_jet_truthflav;
  float *NN_min_dR_b_from_top_to_lep, *NN_min_dR_b_not_from_top_to_lep, *NN_min_dR_not_b_to_lep;
  float *NN_M_lep_bjet_min_dR_from_top, *NN_M_lep_bjet_min_dR_not_from_top;
  float *NN_min_M_lep_bjet_from_top, *NN_max_M_lep_bjet_from_top;
  float *NN_min_M_lep_bjet_not_from_top, *NN_max_M_lep_bjet_not_from_top;
  float *NN_min_M_lep_other_jet, *NN_max_M_lep_other_jet;
  TBranch *NN_topHadronOriginFlag_br = NN_ttree->Branch("topHadronOriginFlag", &NN_tHOF, "topHadronOriginFlag/I");
  TBranch *NN_jet_truthflav_br = NN_ttree->Branch("jet_truthflav", &NN_jet_truthflav, "jet_truthflav/I");
  TBranch *NN_min_dR_b_from_top_to_lep_br = NN_ttree->Branch("min_dR_b_from_top_to_lep", &NN_min_dR_b_from_top_to_lep, "min_dR_b_from_top_to_lep/F");
  TBranch *NN_min_dR_b_not_from_top_to_lep_br = NN_ttree->Branch("min_dR_b_not_from_top_to_lep", &NN_min_dR_b_not_from_top_to_lep, "min_dR_b_not_from_top_to_lep/F");
  TBranch *NN_min_dR_not_b_to_lep_br = NN_ttree->Branch("min_dR_not_b_to_lep", &NN_min_dR_not_b_to_lep, "min_dR_not_b_to_lep/F");
  TBranch *NN_M_lep_bjet_min_dR_from_top_br = NN_ttree->Branch("M_lep_bjet_min_dR_br_from_top", &NN_M_lep_bjet_min_dR_from_top, "M_lep_bjet_min_dR_br_from_top/F");
  TBranch *NN_M_lep_bjet_min_dR_not_from_top_br = NN_ttree->Branch("M_lep_bjet_min_dR_not_from_top", &NN_M_lep_bjet_min_dR_not_from_top, "M_lep_bjet_min_dR_not_from_top/F");
  TBranch *NN_min_M_lep_bjet_from_top_br = NN_ttree->Branch("min_M_lep_bjet_from_top", &NN_min_M_lep_bjet_from_top, "min_M_lep_bjet_from_top/F");
  TBranch *NN_max_M_lep_bjet_from_top_br = NN_ttree->Branch("max_M_lep_bjet_from_top", &NN_max_M_lep_bjet_from_top, "max_M_lep_bjet_from_top/F");
  TBranch *NN_min_M_lep_bjet_not_from_top_br = NN_ttree->Branch("min_M_lep_bjet_not_from_top", &NN_min_M_lep_bjet_not_from_top, "min_M_lep_bjet_not_from_top/F");
  TBranch *NN_max_M_lep_bjet_not_from_top_br = NN_ttree->Branch("max_M_lep_bjet_not_from_top", &NN_max_M_lep_bjet_not_from_top, "max_M_lep_bjet_not_from_top/F");
  TBranch *NN_min_M_lep_other_jet_br = NN_ttree->Branch("min_M_lep_other_jet", &NN_min_M_lep_other_jet, "min_M_lep_other_jet/F");
  TBranch *NN_max_M_lep_other_jet_br = NN_ttree->Branch("max_M_lep_other_jet", &NN_max_M_lep_other_jet, "max_M_lep_other_jet/F");
  for (int entry=0; entry<NN_tHOF_v.size(); entry++) {
    NN_min_dR_b_from_top_to_lep = &NN_min_dR_b_from_top_to_lep_v[entry];
    NN_min_dR_b_not_from_top_to_lep = &NN_min_dR_b_not_from_top_to_lep_v[entry];
    NN_min_dR_not_b_to_lep = &NN_min_dR_not_b_to_lep[entry];
    NN_M_lep_bjet_min_dR_from_top = &NN_M_lep_bjet_min_dR_from_top_v[entry];
    NN_M_lep_bjet_min_dR_not_from_top = &NN_M_lep_bjet_min_dR_not_from_top_v[entry];
    NN_min_M_lep_bjet_from_top = &NN_min_M_lep_bjet_from_top[entry];
    NN_max_M_lep_bjet_from_top = &NN_max_M_lep_bjet_from_top_v[entry];
    NN_min_M_lep_bjet_not_from_top = &NN_min_M_lep_bjet_not_from_top_v[entry];
    NN_max_M_lep_bjet_not_from_top = &NN_max_M_lep_bjet_not_from_top_v[entry];
    NN_min_M_lep_other_jet = &NN_min_M_lep_other_jet_v[entry];
    NN_tHOF = &NN_tHOF_v[entry];
    NN_jet_truthflav = &NN_jet_truthflav_v[entry];
    
    NN_topHadronOriginFlag_br->Fill();
    NN_jet_truthflav_br->Fill();
    NN_min_dR_b_from_top_to_lep_br->Fill();
    NN_min_dR_b_not_from_top_to_lep_br->Fill();
    NN_min_dR_not_b_to_lep_br->Fill();
    NN_M_lep_bjet_min_dR_from_top_br->Fill();
    NN_M_lep_bjet_min_dR_not_from_top_br->Fill();
    NN_min_M_lep_bjet_from_top_br->Fill();
    NN_max_M_lep_bjet_from_top_br->Fill();
    NN_min_M_lep_bjet_not_from_top_br->Fill();
    NN_max_M_lep_bjet_not_from_top_br->Fill();
    NN_min_M_lep_other_jet_br->Fill();
    NN_max_M_lep_other_jet_br->Fill();
    NN_ttree->Fill();
  }
  NN_ttree->Write("nominal", TTree::kOverwrite);
  NN_tfile->Close();
}
