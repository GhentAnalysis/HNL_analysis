#include <iostream>
#include <cmath>
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
  Double_t pigreco= TMath::ACos(-1);
    
    // std::cout << " >>> dummy: " << inputRootFile.c_str() << std::endl;
    std::cout << "---------------------------" << std::endl;
    unsigned year = 0;

//==========================================================================================
    Analysis_mc all(year, "/user/mvit/CMSSW_9_4_4/src/HNL_analysis/sampleLists/2016.txt", "/pnfs/iihe/cms/store/user/mvit/samples/2016/");

    std::string alist = "/user/mvit/CMSSW_9_4_4/src/HNL_analysis/sampleLists/2016.txt";
    std::string adir  = "/pnfs/iihe/cms/store/user/mvit/samples/2016/";
    std::string basename = "shape_file";
    //                                 0   1     2      3      4       5       6        7        8      9      10
    const std::string systNames[] = { "", "pu", "qcd", "pdf", "pEle", "pMuo", "npEle", "npMuo", "jec", "jer", "btag"};
    const size_t nsysts = sizeof(systNames)/sizeof(systNames[0]);

    all.analisi(alist, adir, basename.c_str(), 0, 0);
    /*
    for(size_t i=1; i<nsysts; ++i) {  // skip i=0
      // Skip the following for now:
      // pu      qcd,    pdf,    pEle,   pMuo,   jec,    jer
      if(i==1 || i==4 || i==5 || i==8 || i==9) continue;    
      all.analisi(alist, adir, basename.c_str(), i, 0);  
      all.analisi(alist, adir, basename.c_str(), i, 1);    
      }*/
    return 0;
}
