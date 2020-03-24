#include <iostream>
#include <cmath>
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TCanvas.h"
#include <TLegend.h>
#include <TPad.h>
#include "TF1.h"
#include "TF2.h"
#include <TStyle.h>
#include "TLine.h"
#include "TProfile.h"
#include "TAttFill.h"
#include <iostream>
#include <cstring>
#include <string>
#include <TGraphErrors.h>
#include <Riostream.h>
#include "TFile.h"
#include <TChain.h>
#include <TClonesArray.h>
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include "../../interface/Analysis_mc.h"
#include <THStack.h>
#include <TPaveText.h>
#include <THStack.h>

int main ();


// ********************************************************************
int main(){
  // Are we running in local or on T2B?
  //TString cwd(gSystem->pwd());
  //bool ist2b = cwd.BeginsWith("/storage_mnt");
  //bool isdtl = (ist2b==false && cwd.Contains("trocino"));
  //if(isdtl==false) ist2b = true;

  // Double_t pigreco= TMath::ACos(-1);

  // std::cout << " >>> dummy: " << inputRootFile.c_str() << std::endl;
  std::cout << "---------------------------" << std::endl;
  unsigned year = 0;  // 2016: 0; 2017: 1; 2018: 2; 

  //==========================================================================================
  //Analysis_mc all(year);
  // --> questa cosa e'molto smart XD
  //
  // 2016: "/2016/"; 2017: "/2017/"; 2018: "/2018/"
  std::string adir  = "/pnfs/iihe/cms/store/user/mvit/samples/FINAL/2016/"; // 2016
  if(year==1) adir  = "/pnfs/iihe/cms/store/user/mvit/samples/FINAL/2017/"; // 2017
  if(year==2) adir  = "/pnfs/iihe/cms/store/user/mvit/samples/FINAL/2018/"; // 2018
  //
  // For signal re-weighting
  std::string alist = "jobs/samples/samples_TEMPLABEL.txt";
  //
  // 2016: "2016_new.txt"; 2017: "2017.txt"; 2018: "2018.txt"
  //std::string alist = "/user/mvit/CMSSW_9_4_4/src/HNL_analysis/sampleLists/2018.txt";

  Analysis_mc all(year, alist, adir);

  std::string basename = "shape_file";
  //                            skipData, skipSignal, skipBackground, skipPlotting, skipLimits, skipTables
  all.analisi(basename.c_str(), true    , false     , true          , true        , false     , true, true      /*, 0, 0*/);
  return 0;
}
