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
#include "interface/Analysis_mc.h"
#include <THStack.h>
#include <TPaveText.h>
#include <THStack.h>

int main ();


// ********************************************************************
int main(){
  // Are we running in local or on T2B?
  TString cwd(gSystem->pwd());
  bool ist2b = cwd.BeginsWith("/storage_mnt");
  //bool isdtl = (ist2b==false && cwd.Contains("trocino"));
  //if(isdtl==false) ist2b = true;

  // Double_t pigreco= TMath::ACos(-1);

  // std::cout << " >>> dummy: " << inputRootFile.c_str() << std::endl;
  std::cout << "---------------------------" << std::endl;
  unsigned year = 2;  // 2016: 0; 2017: 1; 2018: 2; 

  //==========================================================================================
  Analysis_mc all(year);
// --> questa cosa e'molto smart XD
  std::string alist = ist2b ?
    "/user/mvit/CMSSW_9_4_4/src/HNL_analysis/sampleLists/2018.txt" :  // 2016: 2016_new.txt; 2017: 2017.txt; 2018: <none>
    "/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/ANALYSIS/20190419_MartinasCode/HNL_analysis/sampleLists/test_daniel.txt";
  std::string adir = ist2b ?
    "/pnfs/iihe/cms/store/user/mvit/samples/FINAL/2018/" :  // 2016: /2016/; 2017: /2017/; 2018: /2018/
    "/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/ANALYSIS/20190419_MartinasCode/HNL_analysis/samples.noSync/";
  std::string basename = "shape_file";

  all.analisi(alist, adir, basename.c_str(), 0, 0);
  // for(size_t i=1; i<nsysts; ++i) {  // skip i=0
  //   // Skip the following for now:
  //   // pu      qcd,    pdf,    pEle,   pMuo,   jec,    jer
  //   //if(i==1 || i==4 || i==5 || i==8 || i==9) continue;    
  //   all.analisi(alist, adir, basename.c_str(), i, 0);  
  //   all.analisi(alist, adir, basename.c_str(), i, 1);    
  // }
  return 0;
}
