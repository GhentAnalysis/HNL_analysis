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
#include <Analysis_mc.h>
#include <THStack.h>
#include <TPaveText.h>
#include <THStack.h>

void majo (unsigned year = 0, int selezione = 1, string inputRootFile= "file_mva_gen1.root");


// ********************************************************************
void majo(unsigned year, int selezione, string inputRootFile){
    Double_t pigreco= TMath::ACos(-1);
    
    std::cout << " >>> dummy: " << inputRootFile.c_str() << std::endl;
    std::cout << "---------------------------" << std::endl;
    

//==========================================================================================
    Analysis_mc all(year, selezione, "/Users/Martina/Desktop/CMS/file_bck/zg.root");
    std::string basename = "stica";
    //                                 0     1      2      3       4       5        6        7      8
    const std::string systNames[] = { "pu", "qcd", "pdf", "pEle", "pMuo", "npEle", "npMuo", "jec", "jer"};
    const size_t nsysts = std::sizeof(systNames)/std::sizeof(systNames[0]);

    all.analisi(selezione, 1, (basename+".root").c_str());
    for(size_t i=0; i<nsysts; ++i) {
      // Skip the following for now:
      // qcd,    pdf,    pEle,   pMuo
      if(i==1 || i==2 || i==3 || i==4) continue;
      all.analisi(selezione, 1, (basename+"_"+systNames[i]+"_down.root" ).c_str(), i, 0);
      all.analisi(selezione, 1, (basename+"_"+systNames[i]+"_up.root"   ).c_str(), i, 1);
    }
}
