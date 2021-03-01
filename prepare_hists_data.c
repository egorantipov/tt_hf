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
  if (files)
    {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next()))
        {
          fname = file->GetName();
          if (fname != "." && fname != "..")
            {
              if (fname.EndsWith(".root")) { container.push_back(dirname + fname); }
              else { container.push_back(dirname + fname + "/"); }
            }
        }
    }
  return container;
}



// ##############
// ##   MAIN   ##
// ##############
void prepare_hists_data()
{
  // Create a list of directories with ntuples
  TString path_to_ntuples =  "/eos/atlas/atlascerngroupdisk/phys-top/ttjets/v4/data/";
  vector<TString> dir_paths = get_list_of_files(path_to_ntuples);
  
  // MET
  TH1 *h_met = new TH1F("h_met", "h_met", 20, 0, 1000);
  TH1 *h_met_phi = new TH1F("h_met_phi", "h_met_phi", 40, -4, 4);
  TH1 *h_jet_pt[3];
  for (int i =0; i<3; i++) {
      TString h_title = "jet_pt_" + to_string(i);
      h_jet_pt[i] = new TH1F(h_title, h_title, 10, 0, 400); }
  
  // jets
  TH1 *h_jets_n = new TH1F("h_jets_n", "h_jets_n", 6, 2, 8);
  TH1 *h_bjets_n = new TH1F("h_bjets_n", "h_bjets_n", 6, 2, 8);

  // NN variables
  TH1 *h_m_bjet_lep_min_dR = new TH1F("h_m_bjet_lep_min_dR", "h_m_bjet_lep_min_dR", 120, 0, 700);
  TH1 *h_m_bjet_lep_min = new TH1F("h_m_bjet_lep_min", "h_m_bjet_lep_min", 120, 0, 700);
  TH1 *h_m_bjet_lep_max = new TH1F("h_m_bjet_lep_max", "h_m_bjet_lep_max", 120, 0, 700);
  TH1 *h_m_bjet_el = new TH1F("h_m_bjet_el", "h_m_bjet_el", 120, 0, 700);
  TH1 *h_m_bjet_mu = new TH1F("h_m_bjet_mu", "h_m_bjet_mu", 120, 0, 700);
  TH1 *h_dR_bjet_lep0 = new TH1F("h_min_dR_bjet_lep0", "h_min_dR_bjet_lep0", 20, 0, 5);
  TH1 *h_dR_bjet_lep1 = new TH1F("h_min_dR_bjet_lep1", "h_min_dR_bjet_lep1", 20, 0, 5);
  TH1 *h_min_dR_bjet_lep = new TH1F("h_min_dR_bjet_lep", "h_min_dR_bjet_lep", 20, 0, 5);
  TH1 *h_min_dR_jet_bjet = new TH1F("h_min_dR_jet_bjet", "h_min_dR_jet_bjet", 20, 0, 5);
  
  
  // Loop over directories with ntuples collections
  for (int dir_counter=0; dir_counter<dir_paths.size(); dir_counter++) {
    // Announce current directory
    cout << "\n\n\n" << dir_paths[dir_counter] << endl;
    
    
    // Check for the content: data/mc? which campaign?
    vector<TString> dir_path_components = split(dir_paths[dir_counter], '/');
    int last_element_index = dir_path_components.size();
    vector<TString> dir_name_components = split(dir_path_components[last_element_index-1], '.');
    bool is_data = false;
    bool is_2015 = false;
    bool is_2016 = false;
    bool is_2017 = false;
    bool is_2018 = false;
    bool is_mc16a = false;
    bool is_mc16d = false;
    bool is_mc16e = false;
    for (int i=0; i<dir_name_components.size(); i++) {
      if (dir_name_components[i] == "periodAllYear") is_data = true; 
      if (dir_name_components[i] == "grp15_v01_p4030") is_2015 = true;
      if (dir_name_components[i] == "grp16_v01_p4030") is_2016 = true;
      if (dir_name_components[i] == "grp17_v01_p4030") is_2017 = true;
      if (dir_name_components[i] == "grp18_v01_p4030") is_2018 = true;
      if (dir_name_components[i] == "mc16a") is_mc16a = true;
      if (dir_name_components[i] == "mc16d") is_mc16d = true;
      if (dir_name_components[i] == "mc16e") is_mc16e = true; }
    
    
    // We work with MC only:For this code working with data
    if (is_data != true) continue;
    cout << dir_path_components[last_element_index-1] << endl;
    
    // Make a list of paths to jobs/DIDs outputs (pieces of a full ntuple)
    vector<TString> paths_to_jobs = get_list_of_files(dir_paths[dir_counter]);
    
      
    // Loop over files
    for (int job_number=0; job_number<paths_to_jobs.size(); job_number++) {  
      TFile *ntuple = new TFile (paths_to_jobs[job_number]);
      TTree *tree_nominal = (TTree*)ntuple->Get("nominal");  
      
      cout << paths_to_jobs[job_number] << endl << endl;
      
      // Set all the needed branches
      vector<Float_t> *jet_pt, *jet_DL1r, *jet_eta, *jet_phi, *jet_e, *mu_pt, *mu_eta, *mu_phi, *mu_charge, *mu_e, *el_pt, *el_eta, *el_phi, *el_charge, *el_e;
      vector<char> *jet_DL1r_77;
      Float_t met, met_phi;
      jet_pt = jet_DL1r = jet_eta = jet_phi = jet_e = mu_pt = mu_eta = mu_phi = mu_charge = mu_e  =  el_pt = el_eta = el_phi = el_charge = el_e  = 0;
      jet_DL1r_77 = 0;
      met = met_phi = 0;
      
      tree_nominal->SetBranchAddress("jet_pt", &jet_pt);
      tree_nominal->SetBranchAddress("jet_eta", &jet_eta);
      tree_nominal->SetBranchAddress("jet_phi", &jet_phi);
      tree_nominal->SetBranchAddress("jet_DL1r", &jet_DL1r);
      tree_nominal->SetBranchAddress("jet_isbtagged_DL1r_77", &jet_DL1r_77);
      tree_nominal->SetBranchAddress("jet_e", &jet_e);
      tree_nominal->SetBranchAddress("el_pt", &el_pt);
      tree_nominal->SetBranchAddress("el_eta", &el_eta);
      tree_nominal->SetBranchAddress("el_phi", &el_phi);
      tree_nominal->SetBranchAddress("el_charge", &el_charge);
      tree_nominal->SetBranchAddress("mu_pt", &mu_pt);
      tree_nominal->SetBranchAddress("mu_eta", &mu_eta);
      tree_nominal->SetBranchAddress("mu_phi", &mu_phi);
      tree_nominal->SetBranchAddress("mu_charge", &mu_charge);
      tree_nominal->SetBranchAddress("mu_e", &mu_e);
      tree_nominal->SetBranchAddress("el_e", &el_e);
      tree_nominal->SetBranchAddress("met_met", &met);
      tree_nominal->SetBranchAddress("met_phi", &met_phi);
	
      
      // Loop over entries
      Int_t nEntries = tree_nominal->GetEntries();
      cout << "\tEntries = " << nEntries << endl;
      for (int entry=0; entry<nEntries; entry++) {
	
	// Show events counter
	if (entry%1000==0) { cout << "\t" << entry << "\r"; cout.flush(); }
	tree_nominal->GetEntry(entry);

	
	// Declare TLorentzVectors
	TLorentzVector el_lvec, mu_lvec;
	vector<TLorentzVector> jets_lvec;
	el_lvec.SetPtEtaPhiE((*el_pt)[0]*0.001, (*el_eta)[0], (*el_phi)[0], (*el_e)[0]*0.001);
	mu_lvec.SetPtEtaPhiE((*mu_pt)[0]*0.001, (*mu_eta)[0], (*mu_phi)[0], (*mu_e)[0]*0.001);
	for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++) {
	  TLorentzVector lvec;
	  lvec.SetPtEtaPhiE((*jet_pt)[jet_i]*0.001, (*jet_eta)[jet_i], (*jet_phi)[jet_i], (*jet_e)[jet_i]*0.001);
	  jets_lvec.push_back(lvec); }
	
	
	// Initiate cuts names
	bool emu_cut = false;
	bool OS_cut = false;
	bool jets_n_cut = false;
	bool bjets_n2_cut = false;
	bool bjets_n3_cut = false;
	
	// Set cuts
	if ((*el_pt).size()==1 && (*mu_pt).size()==1) emu_cut = true;
	if ((*el_charge)[0]!=(*mu_charge)[0]) OS_cut = true;
	
	int bjets_n = 0;
	for (int i=0; i<(*jet_pt).size(); i++) { if ( int((*jet_DL1r_77)[i]==1) ) bjets_n++; }
	if (bjets_n>=2) bjets_n2_cut = true;
	
	int jets_n = (*jet_pt).size();
	if (jets_n >= 3) jets_n_cut = true;
	
	
	// 2+b, emu, OS channel
	if (emu_cut*OS_cut*bjets_n2_cut*jets_n_cut == true) {
	  
	  h_met->Fill(met*0.001);
	  h_met_phi->Fill(met_phi);
	  
	  //jet pT histograms
	  h_jet_pt[0]->Fill((*jet_pt)[0]*0.001);
	  h_jet_pt[1]->Fill((*jet_pt)[1]*0.001);
	  h_jet_pt[2]->Fill((*jet_pt)[2]*0.001);
	  
	  // jets_n
	  h_jets_n->Fill((*jet_pt).size());
	  h_bjets_n->Fill(bjets_n);

	  
	  // NN variables hists
	  double min_dR0 = 999999.;
	  double min_dR1 = 999999.;
	  
	  for (int jet_i=0; jet_i<(*jet_pt).size(); jet_i++) {
	    
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
	      
	      h_dR_bjet_lep0->Fill(dR0);
	      h_dR_bjet_lep1->Fill(dR1);
	      h_min_dR_bjet_lep->Fill( min(dR0, dR1) );
	      
	      // Compute inv masses for bjet-lep pairs
	      // bjet and the closest lepton
	      double m_bjet_lep = 0;
	      double dr_bjet_el = jets_lvec[jet_i].DeltaR(el_lvec);
	      double dr_bjet_mu = jets_lvec[jet_i].DeltaR(mu_lvec);
	      if (dr_bjet_el <= dr_bjet_mu) { m_bjet_lep = (jets_lvec[jet_i] + el_lvec).M(); }
	      else { m_bjet_lep = (jets_lvec[jet_i] + mu_lvec).M(); }
	      if (m_bjet_lep!=0) h_m_bjet_lep_min_dR->Fill(m_bjet_lep);
	      // bjet and el/mu
	      h_m_bjet_el->Fill( (jets_lvec[jet_i] + el_lvec).M() );
	      h_m_bjet_mu->Fill( (jets_lvec[jet_i] + mu_lvec).M() );
	      // bjet and lepton to min/max inv mass
	      double m_max_bjet_lep = max((jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M());
	      double m_min_bjet_lep = min((jets_lvec[jet_i] + el_lvec).M(), (jets_lvec[jet_i] + mu_lvec).M());
	      h_m_bjet_lep_max->Fill(m_max_bjet_lep);
	      h_m_bjet_lep_min->Fill(m_min_bjet_lep);
	      
	    } // [if] - DL1r tagget jet_i
	    
	    // Compute min inv mass for jet-bjet pairs
	    for (int jet_j=0; jet_j<(*jet_pt).size(); jet_j++) {
	      if (jet_i==jet_j) continue;
	      double min_dR_jet_bjet = 999999.;
	      if ( (*jet_DL1r_77)[jet_j]==1) {
		double dR_jet_bjet = jets_lvec[jet_i].DeltaR(jets_lvec[jet_j]);
		min_dR_jet_bjet = min(min_dR_jet_bjet, dR_jet_bjet);
	      } // [if] - DL1r tagger jet_j
	      
	      h_min_dR_jet_bjet->Fill(min_dR_jet_bjet);
	    
	    } // [jet_j] - loop over jets
	    
	  } // [jet_i] - loop over jets
	  
	  
	} // 2+b, emu, OS channel
	    
      } // [entry] - loop over entries
      
      ntuple->Close();
      
    } // [job_number] - loop over data files of one campaign
  
  } // [dir_counter] = loop over 15-16-17-18 campaigns
  


  //Save histograms
  TFile *hists_file = new TFile("hists_data.root", "RECREATE");

  // Regular hists
  h_met->Write("2b_emu_OS_met");
  h_met_phi->Write("2b_emu_OS_met_phi");
  h_jet_pt[0]->Write("2b_emu_OS_jet_pt_0");
  h_jet_pt[1]->Write("2b_emu_OS_jet_pt_1");
  h_jet_pt[2]->Write("2b_emu_OS_jet_pt_2");
  h_jets_n->Write("2b_emu_OS_jets_n");
  h_bjets_n->Write("2b_emu_OS_bjets_n");
  
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


  hists_file->Close();
}
