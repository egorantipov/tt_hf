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
  TH1 *data_met = new TH1F("data_met", "data_met", 20, 0, 1000);
  TH1 *data_met_phi = new TH1F("data_met_phi", "data_met_phi", 40, -4, 4);
  TH1 *data_jet_pt[3];
  for (int i =0; i<3; i++)
    {
      TString h_title = "jet_pt_" + to_string(i);
      data_jet_pt[i] = new TH1F(h_title, h_title, 10, 0, 500);
    }


// Loop over directories with ntuples collections
  for (int dir_counter=0; dir_counter<dir_paths.size(); dir_counter++)
    {
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
      for (int i=0; i<dir_name_components.size(); i++)
	{
	  if (dir_name_components[i] == "periodAllYear") is_data = true; 
          if (dir_name_components[i] == "grp15_v01_p4030") is_2015 = true;
          if (dir_name_components[i] == "grp16_v01_p4030") is_2016 = true;
          if (dir_name_components[i] == "grp17_v01_p4030") is_2017 = true;
          if (dir_name_components[i] == "grp18_v01_p4030") is_2018 = true;
          if (dir_name_components[i] == "mc16a") is_mc16a = true;
          if (dir_name_components[i] == "mc16d") is_mc16d = true;
          if (dir_name_components[i] == "mc16e") is_mc16e = true;
	}


// We work with MC only:For this code working with data
      if (is_data != true) continue;
      cout << dir_path_components[last_element_index-1] << endl;

      // Make a list of paths to jobs/DIDs outputs (pieces of a full ntuple)
      vector<TString> paths_to_jobs = get_list_of_files(dir_paths[dir_counter]);



      // Loop over files
      for (int job_number=0; job_number<paths_to_jobs.size(); job_number++)
	{  TFile *ntuple = new TFile (paths_to_jobs[job_number]);
	   TTree *tree_nominal = (TTree*)ntuple->Get("nominal");  
        
	   cout << paths_to_jobs[job_number] << endl << endl;
	     // Set all the needed branches
	   vector<Float_t> *jet_pt, *jet_DL1r, *jet_eta, *jet_phi, *mu_pt, *mu_eta, *mu_phi, *mu_charge, *mu_e, *el_pt, *el_eta, *el_phi, *el_charge, *el_e;
	     
	      vector<char> *jet_DL1r_77;
	      jet_pt = jet_DL1r = jet_eta = jet_phi = mu_pt = mu_eta = mu_phi = mu_charge = mu_e  =  el_pt = el_eta = el_phi = el_charge = el_e  = 0;
	      jet_DL1r_77 = 0;
	      Float_t met, met_phi;
	      tree_nominal->SetBranchAddress("jet_pt", &jet_pt);
              tree_nominal->SetBranchAddress("jet_eta", &jet_eta);
              tree_nominal->SetBranchAddress("jet_phi", &jet_phi);
              tree_nominal->SetBranchAddress("jet_DL1r", &jet_DL1r);
              tree_nominal->SetBranchAddress("jet_isbtagged_DL1r_77", &jet_DL1r_77);
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
	      for (int entry=0; entry<nEntries; entry++)
		{
		  // Show events counter
		  if (entry%1000==0) { cout << "\t" << entry << "\r"; cout.flush(); }
                  tree_nominal->GetEntry(entry);

		  
		  // Initiate cuts names
		  bool emu_cut = false;
                  bool OS_cut = false;
                  bool jets_n_cut = false;
                  bool bjets_n2_cut = false;
                  bool bjets_n3_cut = false;
                  bool lep_pt_cut = false;
                  bool jet_pt_cut = true;
                  bool all_cuts = false;
                  

		  //Next is cuts
		  if ((*el_pt).size()==1 && (*mu_pt).size()==1) emu_cut = true;
		  if ((*el_charge)[0]!=(*mu_charge)[0]) OS_cut = true;

		  int bjets_n = 0;
		  for (int i=0; i<(*jet_pt).size(); i++) { if ( int((*jet_DL1r_77)[i]==1) ) bjets_n++; }
		  if (bjets_n>=2) bjets_n2_cut = true;

		   // 2+b, emu, OS channel
		  if (emu_cut*OS_cut*bjets_n2_cut == true)
		    {
		      //MET plots and filled with MeV need GeV that is why .001
		      data_met->Fill(met*0.001);
		      data_met_phi->Fill(met_phi);
		      
		      //jet pT histograms
		      
		      data_jet_pt[0]->Fill((*jet_pt)[0]*0.001);
		      data_jet_pt[1]->Fill((*jet_pt)[1]*0.001);
		      data_jet_pt[2]->Fill((*jet_pt)[2]*0.001);

		      
		    }
		  //When setting up the histo do I just call data for my function like how you called mc16_met etc?
		}
	      ntuple->Close();
	}
    }
  //Save histograms
  TFile *hists_file = new TFile("hists_data.root", "RECREATE");

  data_met->Write("2b_emu_OS_met");
  data_met_phi->Write("2b_emu_OS_met_phi");
  data_jet_pt[0]->Write("2b_emu_OS_jet_pt_0");
  data_jet_pt[1]->Write("2b_emu_OS_jet_pt_1");
  data_jet_pt[2]->Write("2b_emu_OS_jet_pt_2");
  hists_file->Close();
}
