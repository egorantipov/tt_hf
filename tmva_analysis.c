// reference: https://github.com/root-project/root/blob/bec9f2d1f778771703f0006a32436e5e23180c00/tutorials/tmva/TMVAClassification.C

#include <TFile.h>
#include <TTree.h>
#include <map>
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void tmva_analysis()
{
  // Default MVA methods to be trained + tested
  map<string, bool> Use;
  Use["BDT"] = false;
  Use["BDTG"] = false;
  Use["BDTB"] = false;
  Use["BDTD"] = false;
  Use["BDTF"] = false;


  // Create an output file
  TFile *outputFile = new TFile("TMVA.root", "RECREATE");
  cout << "==> Created an output file" << endl;
  

  // Create Factory
  TMVA::Factory *factory = new TMVA::Factory("MVAnalysis", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");
  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");
  cout << "==> Created a Factory and a DataLoader" << endl;
  
  
  // Open the input file
  TFile *input = TFile::Open("tt_jets_NN_input.root");
  cout << "==> Opened an input file" << endl;
  

  // Add variables
  //dataloader->AddVariable("topHadronOriginFlag", 'F');
  //dataloader->AddVariable("jet_isbtagged_DL1r_77", 'F');
  dataloader->AddVariable("NN_dR_bjet_lep0", 'F');
  dataloader->AddVariable("NN_dR_bjet_lep1", 'F');
  dataloader->AddVariable("NN_min_dR_bjet_lep", 'F');
  //dataloader->AddVariable("NN_m_bjet_lep_min_dR", 'F');
  dataloader->AddVariable("NN_m_bjet_el", 'F');
  dataloader->AddVariable("NN_m_njet_mu", 'F');
  dataloader->AddVariable("NN_m_bjet_lep_max", 'F');
  //dataloader->AddVariable("NN_m_bjet_lep_min", 'F');
  dataloader->AddVariable("NN_min_dR_jet_bjet", 'F');
  cout << "==> Added variables" << endl;


  // Initialize trees
  dataloader->AddSignalTree( (TTree*)input->Get("NN_signal"), 1.0 );
  dataloader->AddBackgroundTree( (TTree*)input->Get("NN_background"), 1.0);
  cout << "==> Loaded Signal and Bkgd trees" << endl;
  dataloader->PrepareTrainingAndTestTree("", "SplitMode=random:!V");
  cout << "==> Prepared the trees" << endl;
  

  // Book MVA methods
  // There are different BDTs (Boosted Decision Trees)
  /*  
  if (Use["BDTG"]==true) // Gradient Boost
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTG", 
			"!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2");

  if (Use["BDT"]==true) // Adaptive Boost
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
			"!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
  
  if (Use["BDTB"]==true) // Bagging
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTB",
			"!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20");

  if (Use["BDTD"]==true) // Decorrelation + Adaptive Boost
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTD",
			"!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

  if (Use["BDTF"]==true) // Allow using Fisher discriminant in mode splitting for (strong) linearly correlated variables
    factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTF",
			"!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20");
  */
  // Testing BDT setup
  factory->BookMethod(dataloader, TMVA::Types::kBDT, "KBDT",
		      "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate");

  cout << "==> Booked a specific BDT method" << endl;
  

  // Train, test and evaluate
  factory->TrainAllMethods();
  cout << "==> Trained" << endl;
  factory->TestAllMethods();
  cout << "==> Tested" << endl;
  factory->EvaluateAllMethods();
  cout << "==> Evaluated" << endl;

  cout << "==> Wrote root file: " << outputFile->GetName() << endl;
  cout << "==> TMVAClassification is done!" << endl;
  
  outputFile->Close();
  delete factory;
  delete dataloader;

}
