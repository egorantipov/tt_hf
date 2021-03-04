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
	if (!(file->IsDirectory())) { container.push_back(dirname + fname); }
	//if (fname.EndsWith(".root")) { container.push_back(dirname + fname); }
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
void prepare_hists_mc()
{
  // Create a list of directories with ntuples
  TString path_to_ntuples = "/eos/atlas/atlascerngroupdisk/phys-top/ttjets/v4/";
  vector<TString> dir_paths = get_list_of_files(path_to_ntuples);

  
  // Declare histograms
  
  // dR_min between bjets and leptons, 3b channel
  TH1 *h_minDeltaR_lep0_bjets_from_top = new TH1F("minDeltaR_lep1_bjets_fromTop", "minDeltaR_lep1_bjets_fromTop", 20, 0, 5);
  TH1 *h_minDeltaR_lep0_bjets_not_from_top = new TH1F("minDeltaR_lep1_bjets_notFromTop", "minDeltaR_lep1_bjets_notFromTop", 20, 0, 5);
  TH1 *h_minDeltaR_lep1_bjets_from_top = new TH1F("minDeltaR_lep2_bjets_fromTop", "minDeltaR_lep2_bjets_fromTop", 20, 0, 5);
  TH1 *h_minDeltaR_lep1_bjets_not_from_top = new TH1F("minDeltaR_lep2_bjets_notFromTop", "minDeltaR_lep2_bjets_notFromTop", 20, 0, 5);
  
  // dR_min between btags and leptons, 2b channel
  TH1 *h_minDeltaR_lep0_btags_from_top = new TH1F("h_minDeltaR_lep0_btags_from_top", "h_minDeltaR_lep0_btags_from_top", 20, 0, 5);
  TH1 *h_minDeltaR_lep0_btags_not_from_top = new TH1F("h_minDeltaR_lep0_btags_not_from_top", "h_minDeltaR_lep0_btags_not_from_top", 20, 0, 5);;
  TH1 *h_minDeltaR_lep1_btags_from_top = new TH1F("h_minDeltaR_lep1_btags_from_top", "h_minDeltaR_lep1_btags_from_top", 20, 0, 5);
  TH1 *h_minDeltaR_lep1_btags_not_from_top = new TH1F("h_minDeltaR_lep1_btags_not_from_top", "h_minDeltaR_lep1_btags_not_from_top", 20, 0, 5);
  
  // dR_min, 2b channel 
  TH1 *h_minDeltaR_b_from_top_to_b = new TH1F("h_minDeltaR_b_from_top_to_b", "h_minDeltaR_b_from_top_to_b", 20, 0, 5);
  TH1 *h_minDeltaR_b_not_from_top_to_b = new TH1F("h_minDeltaR_b_not_from_top_to_b", "h_minDeltaR_b_not_from_top_to_b", 20, 0, 5);
  TH1 *h_minDeltaR_not_b_to_b = new TH1F("h_minDeltaR_not_b_to_b", "h_minDeltaR_not_b_to_b", 20, 0, 5);
  TH1 *h_minDeltaR_b_from_top_to_jet = new TH1F("h_minDeltaR_b_from_top_to_jet", "h_minDeltaR_b_from_top_to_jet", 20, 0, 5);
  TH1 *h_minDeltaR_b_not_from_top_to_jet = new TH1F("h_minDeltaR_b_not_from_top_to_jet", "h_minDeltaR_b_not_from_top_to_jet", 20, 0, 5);
  TH1 *h_minDeltaR_not_b_to_jet = new TH1F("h_minDeltaR_not_b_to_jet", "h_minDeltaR_not_b_to_jet", 20, 0, 5);
  TH1 *h_minDeltaR_b_from_top_to_lep = new TH1F("h_minDeltaR_b_from_top_to_lep", "h_minDeltaR_b_from_top_to_lep", 20, 0, 5);
  TH1 *h_minDeltaR_b_not_from_top_to_lep = new TH1F("h_minDeltaR_b_not_from_top_to_lep", "h_minDeltaR_b_not_from_top_to_lep", 20, 0, 5);
  TH1 *h_minDeltaR_not_b_to_lep = new TH1F("h_minDeltaR_not_b_to_lep", "h_minDeltaR_not_b_to_lep", 20, 0, 5);
  
  // pT of the three leading jets, 2b channel
  TH1 *h_jet_pt[3];
  for (int i=0; i<3; i++) {
      TString title = "jet_pt_" + to_string(i);
      h_jet_pt[i] = new TH1F(title, title, 10, 0, 400); }
 
  // the first three DL1r tag distributions for 2b1l / 4b / 3b / 2b1c, 2b channel
  TH1 *h_tag0_DL1r[4];
  TH1 *h_tag1_DL1r[4];
  TH1 *h_tag2_DL1r[4];
  for (int i=0; i<4; i++) {
      TString h_title0 = "h_tag0_DL1r_TopHFFF" + to_string(i);
      TString h_title1 = "h_tag1_DL1r_TopHFFF" + to_string(i);
      TString h_title2 = "h_tag2_DL1r_TopHFFF" + to_string(i);
      h_tag0_DL1r[i] = new TH1F(h_title0, h_title0, 30, -15, 15);
      h_tag1_DL1r[i] = new TH1F(h_title1, h_title1, 30, -15, 15);
      h_tag2_DL1r[i] = new TH1F(h_title2, h_title2, 30, -15, 15); }

  // MET, 2b channel
  TH1 *h_met = new TH1F("h_met", "h_met", 20, 0, 1000);
  TH1 *h_met_phi = new TH1F("h_met_phi", "h_met_phi", 40, -4, 4);

  // jets_n, 2b channel
  TH1 *h_jets_n = new TH1F("h_jets_n", "h_jets_n", 6, 2, 8);
  TH1 *h_bjets_n = new TH1F("h_bjets_n", "h_bjets_n", 6, 2, 8);
  
  // leptons, 2b channel
  TH1 *h_lep0_pt = new TH1F("h_lep0_pt", "h_lep0_pt", 20, 0, 1000);
  TH1 *h_lep1_pt = new TH1F("h_lep1_pt", "h_lep1_pt", 20, 0, 1000);
  TH1 *h_lep_pt = new TH1F("h_lep_pt", "h_lep_pt", 20, 0, 1000);
  TH1 *h_lep0_eta = new TH1F("h_lep0_eta", "h_let0_eta", 20, -5, 5);
  TH1 *h_lep1_eta = new TH1F("h_lep1_eta", "h_lep1_eta", 20, -5, 5);
  TH1 *h_lep_eta = new TH1F("h_lep_eta", "h_lep_eta", 20, -5, 5);
  TH1 *h_lep0_phi = new TH1F("h_lep0_phi", "h_lep0_phi", 40, -4, 4);
  TH1 *h_lep1_phi = new TH1F("h_lep1_phi", "h_lep1_phi", 40, -4, 4);
  TH1 *h_lep_phi = new TH1F("h_lep_phi", "h_lep_phi", 40, -4, 4);
  TH1 *h_dR_lep0_lep1 = new TH1F("h_dR_lep0_lep1", "h_dR_lep0_lep1", 20, 0, 5);
  
  // Invariant mass, 2b channel
  TH1 *h_inv_mass_lep_bjet_from_top_min_dR = new TH1F("h_inv_mass_lep_bjet_from_top_min_dR", "h_inv_mass_lep_bjet_from_top_min_dR", 1000, 0, 1000);
  TH1 *h_inv_mass_lep_bjet_not_from_top_min_dR = new TH1F("h_inv_mass_lep_bjet_not_from_top_min_dR", "h_inv_mass_lep_bjet_not_from_top_min_dR", 1000, 0, 1000);
  TH1 *h_inv_mass_lep_btag_from_top_min_dR = new TH1F("h_inv_mass_lep_btag_from_top_min_dR", "h_inv_mass_lep_btag_from_top_min_dR", 1000, 0, 1000);
  TH1 *h_inv_mass_lep_btag_not_from_top_min_dR = new TH1F("h_inv_mass_lep_btag_not_from_top_min_dR", "h_inv_mass_lep_btag_not_from_top_min_dR", 1000, 0, 1000);
  TH1 *h_min_inv_mass_lep_bjet_from_top = new TH1F("h_min_inv_mass_lep_bjet_from_top", "h_min_inv_mass_lep_bjet_from_top", 1000, 0, 1000);
  TH1 *h_max_inv_mass_lep_bjet_from_top = new TH1F("h_max_inv_mass_lep_bjet_from_top", "h_max_inv_mass_lep_bjet_from_top", 1000, 0, 1000);
  TH1 *h_min_inv_mass_lep_bjet_not_from_top = new TH1F("h_min_inv_mass_lep_bjet_not_from_top", "h_min_inv_mass_lep_bjet_not_from_top", 1000, 0, 1000);
  TH1 *h_max_inv_mass_lep_bjet_not_from_top = new TH1F("h_max_inv_mass_lep_bjet_not_from_top", "h_max_inv_mass_lep_bjet_not_from_top", 1000, 0, 1000);
  TH1 *h_min_inv_mass_lep_other_jet = new TH1F("h_min_inv_mass_lep_other_jet", "h_min_inv_mass_lep_other_jet", 1000, 0, 1000);
  TH1 *h_max_inv_mass_lep_other_jet = new TH1F("h_max_inv_mass_lep_other_jet", "h_max_inv_mass_lep_other_jet", 1000, 0, 1000);
  
  // NN variables hists
  TH1 *h_m_bjet_lep_min_dR = new TH1F("h_m_bjet_lep_min_dR", "h_m_bjet_lep_min_dR", 120, 0, 700);
  TH1 *h_m_bjet_lep_min = new TH1F("h_m_bjet_lep_min", "h_m_bjet_lep_min", 120, 0, 700);
  TH1 *h_m_bjet_lep_max = new TH1F("h_m_bjet_lep_max", "h_m_bjet_lep_max", 120, 0, 700);
  TH1 *h_m_bjet_el = new TH1F("h_m_bjet_el", "h_m_bjet_el", 120, 0, 700);
  TH1 *h_m_bjet_mu = new TH1F("h_m_bjet_mu", "h_m_bjet_mu", 120, 0, 700);
  TH1 *h_dR_bjet_lep0 = new TH1F("h_min_dR_bjet_lep0", "h_min_dR_bjet_lep0", 20, 0, 5);
  TH1 *h_dR_bjet_lep1 = new TH1F("h_min_dR_bjet_lep1", "h_min_dR_bjet_lep1", 20, 0, 5);
  TH1 *h_min_dR_bjet_lep = new TH1F("h_min_dR_bjet_lep", "h_min_dR_bjet_lep", 20, 0, 5);
  TH1 *h_min_dR_jet_bjet = new TH1F("h_min_dR_jet_bjet", "h_min_dR_jet_bjet", 20, 0, 5);

  // Study topHFF overlap removal
  TH1 *h_bjets_n_410472 = new TH1F("bjets_n_410472", "bjets_n_410472", 6, 2, 8);
  TH1 *h_bjets_n_411076 = new TH1F("bjets_n_411076", "bjets_n_411076", 6, 2, 8);
  TH1 *h_bjets_n_411077 = new TH1F("bjets_n_411077", "bjets_n_411077", 6, 2, 8);
  TH1 *h_bjets_n_411078 = new TH1F("bjets_n_411078", "bjets_n_411078", 6, 2, 8);
  TH1 *h_topHFFF_410472 = new TH1F("topHFFF_410472", "topHFFF_410472", 5, 0, 5);
  TH1 *h_topHFFF_411076 = new TH1F("topHFFF_411076", "topHFFF_411076", 5, 0 ,5);
  TH1 *h_topHFFF_411077 = new TH1F("topHFFF_411077", "topHFFF_411077", 5, 0 ,5);
  TH1 *h_topHFFF_411078 = new TH1F("topHFFF_411078", "topHFFF_411078", 5, 0 ,5);
  

  // Initialize KLFitter
  KLFitter::Fitter fitter{};
  
  KLFitter::DetectorAtlas_8TeV detector{"/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/AnalysisTop/KLFitterTFs/mc15c/akt4_EMtopo_PP6"};
  fitter.SetDetector(&detector);

  KLFitter::LikelihoodTopDilepton likelihood{};
  likelihood.PhysicsConstants()->SetMassTop(172.5);
  likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kNotag);
  likelihood.SetFlagTopMassFixed(true);
  fitter.SetLikelihood(&likelihood);


  // Create vars for NN
  vector<vector<int>> NN_tHOF_v, NN_jet_DL1r_77_v;


  // Loop over directories with ntuples collections
  for (int dir_counter=0; dir_counter<dir_paths.size(); dir_counter++)
    {
      // Announce current directory
      TSystemFile dir(dir_paths[dir_counter], dir_paths[dir_counter]);
      if (!(dir.IsDirectory())) continue; 
      cout << "\n\n\n" << dir_paths[dir_counter] << endl;
      
      
      // Check for the content: data/mc? which campaign?
      vector<TString> dir_path_components = split(dir_paths[dir_counter], '/');
      int last_element_index = dir_path_components.size();
      vector<TString> dir_name_components = split(dir_path_components[last_element_index-1], '_');
      bool is_data = false;
      bool is_mc16a = false;
      bool is_mc16d = false;
      bool is_mc16e = false;
      for (int i=0; i<dir_name_components.size(); i++) {
	  if (dir_name_components[i] == "data") is_data = true;
          if (dir_name_components[i] == "mc16a") is_mc16a = true;
          if (dir_name_components[i] == "mc16d") is_mc16d = true;
          if (dir_name_components[i] == "mc16e") is_mc16e = true; }
      
      
      // We work with MC only
      if (is_data == true) continue;

      
      // Only nominal trees
      dir_paths[dir_counter] += "nominal/";
      

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
		  
		  if (runNumber==284500) {
		    campaign_lumi = 3.21956 + 32.9881;
		    if (job_DID=="411076") {
		      sumWeights = 3.33006*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.008814;
		      kFactor = 1.1397; }
		    if (job_DID=="411077") {
		      sumWeights = 3.61088*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.046655;
		      kFactor = 1.1398; }
		    if (job_DID=="411078") {
		      sumWeights = 3.61598*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.039503;
		      kFactor = 1.1397; }
		    if (job_DID=="410472") {
		      sumWeights = 5.82869*pow(10, 10);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.10547;
		      kFactor = 1.13975636159; } }
		  if (runNumber==300000) {
		    campaign_lumi = 44.3074;
		    if (job_DID=="411076") {
		      sumWeights = 4.21891*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.008814;
		      kFactor = 1.1397; }
		    if (job_DID=="411077") {
		      sumWeights = 4.49595*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.046655;
		      kFactor = 1.1398; }
		    if (job_DID=="411078") {
		      sumWeights = 4.49400*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.039503;
		      kFactor = 1.1397; }
		    if (job_DID=="410472") {
		      sumWeights = 7.26510*pow(10, 10);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.10547;
		      kFactor = 1.13975636159; } }
		  if (runNumber==310000) {
		    campaign_lumi = 58.4501;
		    if (job_DID=="411076") {
		      sumWeights = 5.47811*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.008814;
		      kFactor = 1.1397; }
		    if (job_DID=="411077") {
		      sumWeights = 5.94763*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.046655;
		      kFactor = 1.1398; }
		    if (job_DID=="411078") {
		      sumWeights = 5.94190*pow(10, 9);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.039503;
		      kFactor = 1.1397; }
		    if (job_DID=="410472") {
		      sumWeights = 1.01641*pow(10, 11);
		      campaign_xsection = 0.72977;
		      campaign_genFiltEff = 0.10547;
		      kFactor = 1.13975636159; } }
		  
		  // Actual computation:
		  weight_lumi = campaign_lumi * campaign_xsection * pow(10,6) * campaign_genFiltEff * kFactor / sumWeights;
		  double weights = w_mc * w_pu * w_leptonSF * w_DL1r_77 * w_jvt * weight_lumi;

		  
		  // Initiate cuts names
		  bool emu_cut = false;
                  bool OS_cut = false;
                  bool jets_n_cut = false;
		  bool btags_n2_cut = false;
                  bool bjets_n2_cut = false;
                  bool bjets_n3_cut = false;
                  bool topHFFF_cut = false;
		  

		  // Define cuts themselves
		  if ((*el_pt).size()==1 && (*mu_pt).size()==1) emu_cut = true;
		  if ((*el_charge)[0]!=(*mu_charge)[0]) OS_cut = true;
		  
		  int bjets_n = 0;
		  for (int i=0; i<(*jet_pt).size(); i++) { if ( int((*jet_truthflav)[i]==1) ) bjets_n++; }
		  if (bjets_n==3) bjets_n3_cut = true;
                  if (bjets_n>=2) bjets_n2_cut = true;
                  
		  int jets_n = (*jet_pt).size();
                  if (jets_n >=3) jets_n_cut = true;
		  
		  int btags_n = 0;
		  for (int i=0; i<(*jet_pt).size(); i++) { if ((*jet_DL1r_77)[i]==1) btags_n++; }
		  if (btags_n >=2) btags_n2_cut = true;
                  
		  if ( only_410472==true || ( (topHFFF==1 && job_DID=="411076") || (topHFFF==2 && job_DID=="411077") || (topHFFF==3 && job_DID=="411078") || (topHFFF==0 && job_DID=="410472") ) ) topHFFF_cut = true;


		  
		  // TLorentzVectors for leptons and jets
		  TLorentzVector el_lvec;
		  TLorentzVector mu_lvec;
		  vector<TLorentzVector> jets_lvec;
		  el_lvec.SetPtEtaPhiE((*el_pt)[0]*0.001, (*el_eta)[0], (*el_phi)[0], (*el_e)[0]*0.001);
		  mu_lvec.SetPtEtaPhiE((*mu_pt)[0]*0.001, (*mu_eta)[0], (*mu_phi)[0], (*mu_e)[0]*0.001);
		  for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++) {
		    TLorentzVector lvec;
		    lvec.SetPtEtaPhiE((*jet_pt)[jet_i]*0.001, (*jet_eta)[jet_i], (*jet_phi)[jet_i], (*jet_e)[jet_i]*0.001);
		    jets_lvec.push_back(lvec); }


		  
		  // 3b, emu, OS channel: Draw min_dR
		  if (emu_cut*OS_cut*bjets_n3_cut*topHFFF_cut*jets_n_cut == true)
		    {
		      
		      // Define min_dR as some lagre value to begin with
                      double min_dR1_top = 999999.; // leading lepton
                      double min_dR2_top = 999999.; // subleading lepton
                      double min_dR1_not_top = 999999.; // leading lepton
                      double min_dR2_not_top = 999999.; // subleading lepton
		      
		      
		      // Loop over all jets, and select only b-tagged in the loop
		      for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++) {
			if ( (*jet_truthflav)[jet_i]==5) {
			      
			    // Define initial dR's
			    double dR1 = 0;
			    double dR2 = 0;
			    
			    // Assign dR1 to the leading lep and dR2 to the subleading 
			    if ((*mu_pt)[0]>(*el_pt)[0]) {
			      dR1 = mu_lvec.DeltaR(jets_lvec[jet_i]);
			      dR2 = el_lvec.DeltaR(jets_lvec[jet_i]); }
			    else {
			      dR1 = el_lvec.DeltaR(jets_lvec[jet_i]);
			      dR2 = mu_lvec.DeltaR(jets_lvec[jet_i]); }
			    
			    // Sort wrt origin
			    if ((*topHadronOriginFlag)[jet_i]==4) { 
			      min_dR1_top = min(min_dR1_top, dR1); min_dR2_top = min(min_dR2_top, dR2); }
			    else {
			      min_dR1_not_top = min(min_dR1_not_top, dR1); min_dR2_not_top = min(min_dR2_not_top, dR2); }
			  
			  }  // selection of b-jets from all jets
			
		      } // [jet_i] - loop over jets
		      
		      
		      // Fill min_dR histograms
		      h_minDeltaR_lep0_bjets_from_top->Fill(min_dR1_top, weights);
		      h_minDeltaR_lep0_bjets_not_from_top->Fill(min_dR1_not_top, weights);
		      h_minDeltaR_lep1_bjets_from_top->Fill(min_dR2_top, weights);
		      h_minDeltaR_lep1_bjets_not_from_top->Fill(min_dR2_not_top, weights);
		      
		    } // 3b, emu, OS cuts
		  
		  

		  // 2+b (tags), emu, OS
		  if (emu_cut*OS_cut*btags_n2_cut*topHFFF_cut*jets_n_cut == true) {

		    // MET hists:
		    h_met->Fill(met*0.001, weights);
		    h_met_phi->Fill(met_phi, weights);
		    
		    
		    // jet pt hists:
		    for (int i=0; i<3; i++) { h_jet_pt[i]->Fill((*jet_pt)[i]*0.001, weights); }
		    

		    // jets_n hist:
		    h_jets_n->Fill((*jet_pt).size(), weights);
		    h_bjets_n->Fill(btags_n, weights);
		    

		    // leptons hists:
		    if ( (*el_pt)[0] > (*mu_pt)[0] ) {
		      h_lep0_pt->Fill((*el_pt)[0]*0.001, weights);
		      h_lep1_pt->Fill((*mu_pt)[0]*0.001, weights);
		      h_lep0_eta->Fill((*el_eta)[0], weights);
		      h_lep1_eta->Fill((*mu_eta)[0], weights);
		      h_lep0_phi->Fill((*el_phi)[0], weights);
		      h_lep1_phi->Fill((*mu_phi)[0], weights); }
		    else {
		      h_lep0_pt->Fill((*mu_pt)[0]*0.001, weights);
		      h_lep1_pt->Fill((*el_pt)[0]*0.001, weights);
		      h_lep0_eta->Fill((*mu_eta)[0], weights);
		      h_lep1_eta->Fill((*el_eta)[0], weights);
		      h_lep0_phi->Fill((*mu_phi)[0], weights);
		      h_lep1_phi->Fill((*el_phi)[0], weights); }
		    h_lep_pt->Fill((*el_pt)[0]*0.001, weights);
		    h_lep_pt->Fill((*mu_pt)[0]*0.001, weights);
		    h_lep_eta->Fill((*el_eta)[0], weights);
		    h_lep_eta->Fill((*mu_eta)[0], weights);
		    h_lep_phi->Fill((*el_phi)[0], weights);
		    h_lep_phi->Fill((*mu_phi)[0], weights);
		  
		    
		    // dR(lep0, lep1) hist:
		    h_dR_lep0_lep1->Fill(el_lvec.DeltaR(mu_lvec));
		    
		    
		    // dR_min lep0/1 btags
		    double min_dR1_top = 999999.; // leading lepton
		    double min_dR2_top = 999999.; // subleading lepton
		    double min_dR1_not_top = 999999.; // leading lepton
		    double min_dR2_not_top = 999999.; // subleading lepton
		    
		    for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++) {
		      if ((*jet_DL1r_77)[jet_i]==1) {

			double dR1 = 0;
			double dR2 = 0;

			if ((*mu_pt)[0]>(*el_pt)[0]) {
			  dR1 = mu_lvec.DeltaR(jets_lvec[jet_i]);
			  dR2 = el_lvec.DeltaR(jets_lvec[jet_i]); }
			else {
			  dR1 = mu_lvec.DeltaR(jets_lvec[jet_i]);
			  dR2 = mu_lvec.DeltaR(jets_lvec[jet_i]); }
		      
			if ((*topHadronOriginFlag)[jet_i]==4) {
			  min_dR1_top = min(min_dR1_top, dR1); 
			  min_dR2_top = min(min_dR2_top, dR2);}
			else {
			  min_dR1_not_top = min(min_dR1_not_top, dR1); 
			  min_dR2_not_top = min(min_dR2_not_top, dR2); }
			
		      } // [if] DL1r tag 
		    } // [jet_i] - loop over jets
		    
		    h_minDeltaR_lep0_btags_from_top->Fill(min_dR1_top, weights);
		    h_minDeltaR_lep0_btags_not_from_top->Fill(min_dR1_not_top, weights);
		    h_minDeltaR_lep1_btags_from_top->Fill(min_dR2_top, weights);
		    h_minDeltaR_lep1_btags_not_from_top->Fill(min_dR2_not_top, weights);
		
		    

		    // NN variables
		    vector<int> NN_tHOF_v_one_event;
		    vector<int> NN_jet_DL1r_77_v_one_event;
		  
		    double min_dR0 = 999999.;
		    double min_dR1 = 999999.;

		    for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++) {

		      NN_tHOF_v_one_event.push_back((*topHadronOriginFlag)[jet_i]);
		      NN_jet_DL1r_77_v_one_event.push_back((*jet_DL1r_77)[jet_i]);
		      
		      if ((*jet_DL1r_77)[jet_i]==1) {

			// Compute min_dR bjet-lep0/1
			double dR0 = 0;
			double dR1 = 0;

			if ((*mu_pt)[0]>(*el_pt)[0]) {
			  dR0 = mu_lvec.DeltaR(jets_lvec[jet_i]);
			  dR1 = el_lvec.DeltaR(jets_lvec[jet_i]); }
			else {
			  dR0 = mu_lvec.DeltaR(jets_lvec[jet_i]);
			  dR1 = mu_lvec.DeltaR(jets_lvec[jet_i]); }

			h_dR_bjet_lep0->Fill(dR0, weights);
			h_dR_bjet_lep1->Fill(dR1, weights);
			h_min_dR_bjet_lep->Fill( min(dR0, dR1), weights );

			// Compute inv masses for bjet-lep pairs
			// bjet and the closest lepton
			double m_bjet_lep = 0;
			double dr_bjet_el = jets_lvec[jet_i].DeltaR(el_lvec);
			double dr_bjet_mu = jets_lvec[jet_i].DeltaR(mu_lvec);
			if (dr_bjet_el <= dr_bjet_mu) { m_bjet_lep = (jets_lvec[jet_i] + el_lvec).M(); }
			else { m_bjet_lep = (jets_lvec[jet_i] + mu_lvec).M(); }
			if (m_bjet_lep!=0) h_m_bjet_lep_min_dR->Fill(m_bjet_lep, weights);
			// bjet and el/mu
			h_m_bjet_el->Fill( (jets_lvec[jet_i] + el_lvec).M() , weights);
			h_m_bjet_mu->Fill( (jets_lvec[jet_i] + mu_lvec).M() , weights);
			// bjet and lepton to min/max inv mass
			double m_max_bjet_lep = max((jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M());
			double m_min_bjet_lep = min((jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M());
			h_m_bjet_lep_max->Fill(m_max_bjet_lep, weights);
			h_m_bjet_lep_min->Fill(m_min_bjet_lep, weights);

		      } // [if] - DL1r tagget jet_i

		      // Compute min inv mass for jet-bjet pairs
		      for (int jet_j=0; jet_j<(*jet_pt).size(); jet_j++) {
			if (jet_i==jet_j) continue;
			double min_dR_jet_bjet = 999999.;
			if ( (*jet_DL1r_77)[jet_j]==1) {
			  double dR_jet_bjet = jets_lvec[jet_i].DeltaR(jets_lvec[jet_j]);
			  min_dR_jet_bjet = min(min_dR_jet_bjet, dR_jet_bjet);
			} // [if] - DL1r tagger jet_j

			h_min_dR_jet_bjet->Fill(min_dR_jet_bjet, weights);

		      } // [jet_j] - loop over jets

		    } // [jet_i] - loop over jets
		    
		    NN_tHOF_v.push_back(NN_tHOF_v_one_event);
                    NN_jet_DL1r_77_v.push_back(NN_jet_DL1r_77_v_one_event);

		    
		    
		    // topHFFF overlap removal studies
		    if (job_DID=="411076") { h_bjets_n_411076->Fill(btags_n, weights); h_topHFFF_411076->Fill(topHFFF, weights); }
		    if (job_DID=="411077") { h_bjets_n_411077->Fill(btags_n, weights); h_topHFFF_411077->Fill(topHFFF, weights); }
		    if (job_DID=="411078") { h_bjets_n_411078->Fill(btags_n, weights); h_topHFFF_411078->Fill(topHFFF, weights); }
		    if (job_DID=="410472") { h_bjets_n_410472->Fill(btags_n, weights); h_topHFFF_410472->Fill(topHFFF, weights); }

		  } // 2+b (tags) selection

		  
		  
		  // 2+b (jets), emu, OS channel
		  if (emu_cut*OS_cut*bjets_n2_cut*topHFFF_cut*jets_n_cut == true) {
		      
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
			  
			  // dR_min 
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]==4 && (*jet_truthflav)[j]==5) {
			    double dR_b_from_top_to_b = jets_lvec[i].DeltaR(jets_lvec[j]);
			    if (dR_b_from_top_to_b < min_dR_b_from_top_to_b) min_dR_b_from_top_to_b = dR_b_from_top_to_b; }
			  
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]!=4 && (*jet_truthflav)[j]==5) {
			    double dR_b_not_from_top_to_b = jets_lvec[i].DeltaR(jets_lvec[j]);
			    if (dR_b_not_from_top_to_b < min_dR_b_not_from_top_to_b) min_dR_b_not_from_top_to_b = dR_b_not_from_top_to_b; }
			  
			  if ((*jet_truthflav)[i]!=5 && (*topHadronOriginFlag)[i]!=4 && (*jet_truthflav)[j]==5) {
			    double dR_not_b_to_b = jets_lvec[i].DeltaR(jets_lvec[j]);
			    if (dR_not_b_to_b < min_dR_not_b_to_b) min_dR_not_b_to_b = dR_not_b_to_b; }
			  
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]==4 && (*jet_truthflav)[j]!=5) {
			    double dR_b_from_top_to_jet = jets_lvec[i].DeltaR(jets_lvec[j]);
			    if (dR_b_from_top_to_jet < min_dR_b_from_top_to_jet) min_dR_b_from_top_to_jet = dR_b_from_top_to_jet; }
			  
			  if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]!=4 && (*jet_truthflav)[j]!=5) {
			    double dR_b_not_from_top_to_jet = jets_lvec[i].DeltaR(jets_lvec[j]);
			    if (dR_b_not_from_top_to_jet < min_dR_b_not_from_top_to_jet) min_dR_b_not_from_top_to_jet = dR_b_not_from_top_to_jet; }
			  
			  if ((*jet_truthflav)[i]!=5 && (*topHadronOriginFlag)[i]!=4 && (*jet_truthflav)[j]!=5) {
			    double dR_not_b_to_jet = jets_lvec[i].DeltaR(jets_lvec[j]);
			    if (dR_not_b_to_jet < min_dR_not_b_to_jet) min_dR_not_b_to_jet = dR_not_b_to_jet; }
			} // loop over jet[j]
			
			if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]==4) {
			  double dR_b_from_top_to_el = jets_lvec[i].DeltaR(el_lvec);
			  double dR_b_from_top_to_mu = jets_lvec[i].DeltaR(mu_lvec);
			  double dR_b_from_top_to_lep = min(dR_b_from_top_to_el, dR_b_from_top_to_mu);
			  if (dR_b_from_top_to_lep < min_dR_b_from_top_to_lep) min_dR_b_from_top_to_lep = dR_b_from_top_to_lep; }
			
			if ((*jet_truthflav)[i]==5 && (*topHadronOriginFlag)[i]!=4) {
			  double dR_b_not_from_top_to_el = jets_lvec[i].DeltaR(el_lvec);
			  double dR_b_not_from_top_to_mu = jets_lvec[i].DeltaR(mu_lvec);
			  double dR_b_not_from_top_to_lep = min(dR_b_not_from_top_to_el, dR_b_not_from_top_to_mu);
			  if (dR_b_not_from_top_to_lep < min_dR_b_not_from_top_to_lep) min_dR_b_not_from_top_to_lep = dR_b_not_from_top_to_lep; }
			
			if ((*jet_truthflav)[i]!=5 && (*topHadronOriginFlag)[i]!=4) {
			  double dR_not_b_to_el = jets_lvec[i].DeltaR(el_lvec);
			  double dR_not_b_to_mu = jets_lvec[i].DeltaR(mu_lvec);;
			  double dR_not_b_to_lep = min(dR_not_b_to_el, dR_not_b_to_mu);
			  if (dR_not_b_to_lep < min_dR_not_b_to_lep) min_dR_not_b_to_lep = dR_not_b_to_lep; }
			
		      } // loop over jet[i]
		      
		      h_minDeltaR_b_from_top_to_b->Fill(min_dR_b_from_top_to_b, weights);
		      h_minDeltaR_b_not_from_top_to_b->Fill(min_dR_b_not_from_top_to_b, weights);
		      h_minDeltaR_not_b_to_b->Fill(min_dR_not_b_to_b, weights);
		      h_minDeltaR_b_from_top_to_jet->Fill(min_dR_b_from_top_to_jet, weights);
		      h_minDeltaR_b_not_from_top_to_jet->Fill(min_dR_b_not_from_top_to_jet, weights);
		      h_minDeltaR_not_b_to_jet->Fill(min_dR_not_b_to_jet, weights);
		      h_minDeltaR_b_from_top_to_lep->Fill(min_dR_b_from_top_to_lep, weights);
		      h_minDeltaR_b_not_from_top_to_lep->Fill(min_dR_b_not_from_top_to_lep, weights);
		      h_minDeltaR_not_b_to_lep->Fill(min_dR_not_b_to_lep, weights);


		      // KLFitter invariant mass calculations:
		      //KLFitter::Particles particles{};
                      //likelihood.SetLeptonType(KLFitter::LikelihoodTopDilepton::kElectron, KLFitter::LikelihoodTopDilepton::kMuon);
		      // Add leptons
                      //particles.AddParticle(el_lvec, (*el_cl_eta)[0], (*el_charge)[0], KLFitter::Particles::kElectron);
		      //particles.AddParticle(mu_lvec, mu_lvec.Eta(), (*mu_charge)[0], KLFitter::Particles::kMuon);
                      // Add two leading pT jets
		      //for (int jet_i=0; jet_i<3; jet_i++) { particles.AddParticle(jet_lvec, jet_lvec.Eta(), KLFitter::Particles::kParton, "", jet_i); }
		      //fitter.SetParticles(&particles);
		      // Add MET
		      //fitter.SetET_miss_XY_SumET(met*0.001*cos(met_phi), met*0.001*sin(met_phi), met*0.001);

		      
		      // Loop over permutations
		      //int n_perm = fitter.Permutations()->NPermutations();
		      //for (int perm_i=0; perm_i < n_perm; perm_i++) {
		      //fitter.Fit(perm_i);
		      //auto permutedParticles = fitter.Likelihood()->PParticlesPermuted();
		      //double llh = fitter.Likelihood()->LogLikelihood(fitter.Likelihood()->GetBestFitParameters()); }
		      
		      

		      // Invariant mass of bjet-lepton pairs		      
		      double min_inv_mass_lep_bjet_from_top = 999999;
		      double max_inv_mass_lep_bjet_from_top = 0;
		      double min_inv_mass_lep_bjet_not_from_top = 999999;
		      double max_inv_mass_lep_bjet_not_from_top = 0;
		      double min_inv_mass_lep_other_jet = 999999;
		      double max_inv_mass_lep_other_jet = 0;
		      
		      for (int jet_i=0; jet_i<jet_pt->size(); jet_i++) {
			
			// bjets and the closest leptons
			if ((*jet_truthflav)[jet_i]==5) {
			  if ((*topHadronOriginFlag)[jet_i]==4) {
			    double dr_j_el = jets_lvec[jet_i].DeltaR(el_lvec);
			    double dr_j_mu = jets_lvec[jet_i].DeltaR(mu_lvec);
			    double inv_mass_j_lep = 0;
			    if (dr_j_el <= dr_j_mu) { inv_mass_j_lep = (jets_lvec[jet_i] + el_lvec).M(); }
			    else { inv_mass_j_lep = (jets_lvec[jet_i] + mu_lvec).M(); }
			    if (inv_mass_j_lep != 0) h_inv_mass_lep_bjet_from_top_min_dR->Fill(inv_mass_j_lep, weights); }
			  else { 
			    double dr_j_el = jets_lvec[jet_i].DeltaR(el_lvec);
			    double dr_j_mu = jets_lvec[jet_i].DeltaR(mu_lvec);
			    double inv_mass_j_lep = 0;
			    if (dr_j_el <= dr_j_mu) { inv_mass_j_lep = (jets_lvec[jet_i] + el_lvec).M(); }
			    else { inv_mass_j_lep = (jets_lvec[jet_i] + mu_lvec).M(); }
			    if (inv_mass_j_lep != 0) h_inv_mass_lep_bjet_not_from_top_min_dR->Fill(inv_mass_j_lep, weights); } }
			
			// bjets and leptons - min and max invarinat masses
			if ((*jet_truthflav)[jet_i]==5) {
			  if ((*topHadronOriginFlag)[jet_i]==4) {
			    double min_inv_mass_lep_bjet_from_top_tmp = min( (jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() );
			    double max_inv_mass_lep_bjet_from_top_tmp = max( (jets_lvec[jet_i] + mu_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() );
			    min_inv_mass_lep_bjet_from_top = min(min_inv_mass_lep_bjet_from_top_tmp, min_inv_mass_lep_bjet_from_top);
			    max_inv_mass_lep_bjet_from_top = max(max_inv_mass_lep_bjet_from_top_tmp, max_inv_mass_lep_bjet_from_top_tmp); }
			  else {
			    double min_inv_mass_lep_bjet_not_from_top_tmp = min( (jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() );
			    double max_inv_mass_lep_bjet_not_from_top_tmp = max( (jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() );
			    min_inv_mass_lep_bjet_not_from_top = min(min_inv_mass_lep_bjet_not_from_top_tmp, min_inv_mass_lep_bjet_not_from_top);
			    max_inv_mass_lep_bjet_not_from_top = max(max_inv_mass_lep_bjet_not_from_top_tmp, max_inv_mass_lep_bjet_not_from_top); } }
			// other than bjets and leptons - min and max invariant masses
			else {
			  double min_inv_mass_lep_other_jet_tmp = min( (jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() );
			  double max_inv_mass_lep_other_jet_tmp = max( (jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() );
			  min_inv_mass_lep_other_jet = min( (jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() );
			  max_inv_mass_lep_other_jet = max( (jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M() ); }
		      }
		      
		      // Fill the min/max invariant mass hists
		      if (min_inv_mass_lep_bjet_from_top!=999999) h_min_inv_mass_lep_bjet_from_top->Fill(min_inv_mass_lep_bjet_from_top, weights);
		      if (max_inv_mass_lep_bjet_from_top!=0) h_max_inv_mass_lep_bjet_from_top->Fill(max_inv_mass_lep_bjet_from_top, weights);
		      if (min_inv_mass_lep_bjet_not_from_top!=999999) h_min_inv_mass_lep_bjet_not_from_top->Fill(min_inv_mass_lep_bjet_not_from_top, weights);
		      if (max_inv_mass_lep_bjet_not_from_top!=0) h_max_inv_mass_lep_bjet_not_from_top->Fill(max_inv_mass_lep_bjet_not_from_top, weights);
		      if (min_inv_mass_lep_other_jet!=999999) h_min_inv_mass_lep_other_jet->Fill(min_inv_mass_lep_other_jet, weights);
		      if (max_inv_mass_lep_other_jet!=0) h_max_inv_mass_lep_other_jet->Fill(max_inv_mass_lep_other_jet, weights);

		      
		      

		      // Sort jets wrt DL1r tag weights
		      sort (jet_DL1r->begin(), jet_DL1r->end(), greater<int>());


		      // Fill DL1r tag weight histos for the first three tags, also sort wrt topHFFF
		      h_tag0_DL1r[topHFFF]->Fill((*jet_DL1r)[0], weights);
                      h_tag1_DL1r[topHFFF]->Fill((*jet_DL1r)[1], weights);
                      h_tag2_DL1r[topHFFF]->Fill((*jet_DL1r)[2], weights);

		    } // 2+b, emu, OS cuts 

		} // [entry] - loop over entries
	      
	      // Close ntuple after we're done with it
	      ntuple->Close();
	          
	    } // [ntuple_number] - loop over ntuples of a particular job
	  
	} // [job_number] - loop over jobs (pieces) of a collection.
      
    } // [dir_counter] - loop over directories names with jobs folders: mc16a, mc16d, mc16e, data



  // Save histograms

  TFile *hists_file = new TFile("hists_mc.root", "RECREATE");

  // dR_min between bjets and leptons, 3b channel  
  h_minDeltaR_lep0_bjets_from_top->Write("3b_emu_OS_min_dR_lep0_b_from_top");
  h_minDeltaR_lep1_bjets_from_top->Write("3b_emu_OS_min_dR_lep1_b_from_top");
  h_minDeltaR_lep0_bjets_not_from_top->Write("3b_emu_OS_min_dR_lep0_b_not_from_top");
  h_minDeltaR_lep1_bjets_not_from_top->Write("3b_emu_OS_min_dR_lep1_b_not_from_top");

  // dR_min between btags and leptons, 2b channel
  h_minDeltaR_lep0_btags_from_top->Write("2b_emu_OS_min_dR_lep0_b_from_top");
  h_minDeltaR_lep1_btags_from_top->Write("2b_emu_OS_min_dR_lep1_b_from_top");
  h_minDeltaR_lep0_btags_not_from_top->Write("2b_emu_OS_min_dR_lep0_b_not_from_top");
  h_minDeltaR_lep1_btags_not_from_top->Write("2b_emu_OS_min_dR_lep1_b_not_from_top");

  // dR_min, 2b channel
  h_minDeltaR_b_from_top_to_b->Write("2b_emu_OS_h_minDeltaR_b_from_top_to_b");
  h_minDeltaR_b_not_from_top_to_b->Write("2b_emu_OS_h_minDeltaR_b_not_from_top_to_b");
  h_minDeltaR_not_b_to_b->Write("2b_emu_OS_h_minDeltaR_not_b_to_b");
  h_minDeltaR_b_from_top_to_jet->Write("2b_emu_OS_h_minDeltaR_b_from_top_to_jet");
  h_minDeltaR_b_not_from_top_to_jet->Write("2b_emu_OS_h_minDeltaR_b_not_from_top_to_jet");
  h_minDeltaR_not_b_to_jet->Write("2b_emu_OS_h_minDeltaR_not_b_to_jet");
  h_minDeltaR_b_from_top_to_lep->Write("2b_emu_OS_h_minDeltaR_b_from_top_to_lep");
  h_minDeltaR_b_not_from_top_to_lep->Write("2b_emu_OS_h_minDeltaR_b_not_from_top_to_lep");
  h_minDeltaR_not_b_to_lep->Write("2b_emu_OS_h_minDeltaR_not_b_to_lep");

  // pT of the three leading jets, 2b channel 
  for (int i=0; i<3; i++) { 
    TString title = "2b_emu_OS_jet_pt_" + to_string(i);
    h_jet_pt[i]->Write(title); }
  

  // the first three DL1r tag distributions for 2b1l / 4b / 3b / 2b1c, 2b channel
  for (int topHFFF_i=0; topHFFF_i<4; topHFFF_i++) {
    TString process = "";
    if (topHFFF_i==0) process = "2b1l";
    if (topHFFF_i==1) process = "4b";
    if (topHFFF_i==2) process = "3b";
    if (topHFFF_i==3) process = "2b1c";
    h_tag0_DL1r[topHFFF_i]->Write("DL1r_templates_"+process+"_1st_tag");
    h_tag1_DL1r[topHFFF_i]->Write("DL1r_templates_"+process+"_2nd_tag");
    h_tag2_DL1r[topHFFF_i]->Write("DL1r_templates_"+process+"_3rd_tag"); }

  // MET, 2b channel
  h_met->Write("2b_emu_OS_met");
  h_met_phi->Write("2b_emu_OS_met_phi");
  
  // jets_n, 2b channel
  h_jets_n->Write("2b_emu_OS_jets_n");
  h_bjets_n->Write("2b_emu_OS_bjets_n");
  
  // leptons, 2b channel 
  h_lep0_pt->Write("2b_emu_OS_lep0_pt");
  h_lep1_pt->Write("2b_emu_OS_lep1_pt");
  h_lep_pt->Write("2b_emu_OS_lep_pt");
  h_lep0_eta->Write("2b_emu_OS_lep0_eta");
  h_lep1_eta->Write("2b_emu_OS_lep1_eta");
  h_lep_eta->Write("2b_emu_OS_lep_eta");
  h_lep0_phi->Write("2b_emu_OS_lep0_phi");
  h_lep1_phi->Write("2b_emu_OS_lep1_phi");
  h_lep_phi->Write("2b_emu_OS_lep_phi");
  h_dR_lep0_lep1->Write("2b_emu_OS_dR_lep0_lep1");
  
  // Invariant mass
  h_inv_mass_lep_bjet_from_top_min_dR->Write("2b_emu_OS_h_inv_mass_lep_bjet_from_top_min_dR");
  h_inv_mass_lep_bjet_not_from_top_min_dR->Write("2b_emu_OS_h_inv_mass_lep_bjet_not_from_top_min_dR");
  h_inv_mass_lep_btag_from_top_min_dR->Write("2b_emu_OS_h_inv_mass_lep_btag_from_top_min_dR");
  h_inv_mass_lep_btag_not_from_top_min_dR->Write("2b_emu_OS_h_inv_mass_lep_btag_not_from_top_min_dR");
  h_min_inv_mass_lep_bjet_from_top->Write("2b_emu_OS_h_min_inv_mass_lep_bjet_from_top");
  h_max_inv_mass_lep_bjet_from_top->Write("2b_emu_OS_h_max_inv_mass_lep_bjet_from_top");
  h_min_inv_mass_lep_bjet_not_from_top->Write("2b_emu_OS_h_min_inv_mass_lep_bjet_not_from_top");
  h_max_inv_mass_lep_bjet_not_from_top->Write("2b_emu_OS_h_max_inv_mass_lep_bjet_not_from_top");
  h_min_inv_mass_lep_other_jet->Write("2b_emu_OS_h_min_inv_mass_lep_other_jet");
  h_max_inv_mass_lep_other_jet->Write("2b_emu_OS_h_max_inv_mass_lep_other_jet");

  // NN variables hists 
  h_m_bjet_lep_min_dR->Write("NN__2b_emu_OS_m_bjet_lep_min_dR");
  h_m_bjet_lep_min->Write("NN__2b_emu_OS_m_bjet_lep_min");
  h_m_bjet_lep_max->Write("NN__2b_emu_OS_m_bjet_lep_max");
  h_m_bjet_el->Write("NN__2b_emu_OS_m_bjet_el");
  h_m_bjet_mu->Write("NN__2b_emu_OS_m_bjet_mu");
  h_dR_bjet_lep0->Write("NN__2b_emu_OS_dR_bjet_lep0");
  h_dR_bjet_lep1->Write("NN__2b_emu_OS_dR_bjet_lep1");
  h_min_dR_bjet_lep->Write("NN__2b_emu_OS_min_dR_bjet_lep");
  h_min_dR_jet_bjet->Write("NN__2b_emu_OS_min_dR_jet_bjet");

  // topHFFF overlap removal study
  h_bjets_n_411076->Write("topHFFF_study_2b_emu_OS_bjets_n_411076");
  h_bjets_n_411077->Write("topHFFF_study_2b_emu_OS_bjets_n_411077");
  h_bjets_n_411078->Write("topHFFF_study_2b_emu_OS_bjets_n_411078");
  h_bjets_n_410472->Write("topHFFF_study_2b_emu_OS_bjets_n_410472");
  h_topHFFF_411076->Write("topHFFF_study_2b_emu_OS_topHFFF_411076");
  h_topHFFF_411077->Write("topHFFF_study_2b_emu_OS_topHFFF_411077");
  h_topHFFF_411078->Write("topHFFF_study_2b_emu_OS_topHFFF_411078");
  h_topHFFF_410472->Write("topHFFF_study_2b_emu_OS_topHFFF_410472");

  // Close the hists file
  hists_file->Close();


  // Fill NN ROOT file
  TFile *NN_tfile = new TFile("tt_jets_NN_input.root", "RECREATE");
  TTree *NN_ttree = new TTree("nominal", "NN_input");
  vector<int> *NN_tHOF, *NN_jet_DL1r_77;
  NN_ttree->Branch("topHadronOriginFlag", &NN_tHOF);
  NN_ttree->Branch("jet_isbtagged_DL1r_77", &NN_jet_DL1r_77);
  for (int entry=0; entry<NN_tHOF_v.size(); entry++) {
    NN_tHOF = &NN_tHOF_v[entry];
    NN_jet_DL1r_77 = &NN_jet_DL1r_77_v[entry];
    
    NN_ttree->Fill();
  }
  NN_ttree->Write("nominal", TTree::kOverwrite);
  NN_tfile->Close();
}
