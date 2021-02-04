#include "KLFitter/DetectorAtlas_8TeV.h"
#include "KLFitter/Fitter.h"
#include "KLFitter/LikelihoodTopDilepton.h"
#include "KLFitter/Permutations.h"

#include <iostream>
using namespace std;

void test_klf()
{
  // initialize KLFitter
  KLFitter::Fitter fitter{};

  KLFitter::DetectorAtlas_8TeV detector{"/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/dev/AnalysisTop/KLFitterTFs/mc15c/akt4_EMtopo_PP6"};
  fitter.SetDetector(&detector);

  KLFitter::LikelihoodTopDilepton likelihood{};
  likelihood.PhysicsConstants()->SetMassTop(172.5);
  likelihood.SetBTagging(KLFitter::LikelihoodBase::BtaggingMethod::kNotag);
  likelihood.SetFlagTopMassFixed(true);
  fitter.SetLikelihood(&likelihood);
}
