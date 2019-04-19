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
    all.analisi(selezione, 1, "prova.root");
    //all.analisi(selezione, 1, "prova_qcd.root", 1);
    //all.analisi(selezione, 1, "prova_pdf.root", 2);
    all.analisi(selezione, 1, "prova_jec_down.root" ,  8, 0);
    all.analisi(selezione, 1, "prova_jec_up.root"   ,  8, 1);
    all.analisi(selezione, 1, "prova_jer_down.root" ,  9, 0);
    all.analisi(selezione, 1, "prova_jer_up.root"   ,  9, 1);
    all.analisi(selezione, 1, "prova_btag_down.root", 10, 0);
    all.analisi(selezione, 1, "prova_btag_up.root"  , 10, 1);
    all.analisi(selezione, 1, "prova_stat_down.root", 11, 0);
    all.analisi(selezione, 1, "prova_stat_up.root"  , 11, 1);

}
