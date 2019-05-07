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
    Analysis_mc all(year, "/sampleList/2016.txt", "/pnfs/iihe/cms/store/user/mvit/samples/2016/");
    all.analisi(0, "/sampleList/2016.txt", "/pnfs/iihe/cms/store/user/mvit/samples/2016/","pippo",0,0);

    
    
    //all.analisi(selezione, 1, "prova_qcd.root", 1);
    //all.analisi(selezione, 1, "prova_pdf.root", 2);
    
    /* all.analisi(selezione, 1, "prova_jec_down.root" ,  8, 0);
    all.analisi(selezione, 1, "prova_jec_up.root"   ,  8, 1);
    all.analisi(selezione, 1, "prova_jer_down.root" ,  9, 0);
    all.analisi(selezione, 1, "prova_jer_up.root"   ,  9, 1);
    all.analisi(selezione, 1, "prova_btag_down.root", 10, 0);
    all.analisi(selezione, 1, "prova_btag_up.root"  , 10, 1);
    all.analisi(selezione, 1, "prova_stat_down.root", 11, 0);
    all.analisi(selezione, 1, "prova_stat_up.root"  , 11, 1);*/

    
    return 0;
}
