{
  gSystem->Load("./KLFitter/build/lib/libBAT");
  gSystem->Load("./KLFitter/build/lib/libKLFitter");

  gROOT->ProcessLine(".include ./KLFitter/include");
  gROOT->ProcessLine(".include ./KLFitter/build/include");
  gROOT->ProcessLine(".x prepare_hists_mc.c+");
}
