#ifndef PlotScript
#define PlotScript
 
//import ROOT classes
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TH2D.h"
 //kGreen+3
const double xPad =0.25;
//const Color_t colors[9] ={ 91,  91,91, kRed-10,51, kRed-10, kGreen+3, kGreen, kGreen  };
//const Color_t colors[9] ={ 91,91,91,  2,kViolet+1, kRed-10,kGreen+3, kGreen, 7 };
const Color_t colors[9] ={ kViolet+1,  2,kViolet+1, kRed-10,kGreen+3, kGreen, 7, 8, kViolet+1  };
//const Color_t colors[9] ={ 1,  2,3, 4,5, 6, 7, 8, 9  };

const Color_t sigCols[10] = {1, kAzure +10, kMagenta, 4, kCyan, 4, 93, kRed-3, kBlue-3, kMagenta -4 };
//Set histogram colors and lines
void histcol(TH1D *, const Color_t);
//Return histogram divided by other histogram (both are normalized
TH1D *HistDiv(TH1D *, TH1D *, const bool abs = false);
//Set Histogram labelsizes
void HistLabelSizes(TH1D *h, const double xlabel = 0.045, const double xtitle = 0.05, const double ylabel = 0.045, const double ytitle = 0.045);
void HistLabelSizes(TH2D *h, const double xlabel = 0.045, const double xtitle = 0.05, const double ylabel = 0.045, const double ytitle = 0.045);
//Set Stack colors
void StackCol(TH1D *h, const Color_t);
//Order histograms in terms of yields (biggest first)
void yieldOrder(TH1D**&, unsigned*, const unsigned);
void plotDataVSMC_mu(int categoria,int channel,int istogramma,
		     TH1D* data, TH1D** bkg,
		     const TString* names, const unsigned nHist,
		     const TString& name_cut,const TString& name_channel, const TString& name_histo,
		     const bool ylog,
		     const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm);
void plotDataVSMC_e(int categoria,int channel,int istogramma,
		     TH1D* data, TH1D** bkg,
		     const TString* names, const unsigned nHist,
		     const TString& name_cut,const TString& name_channel, const TString& name_histo,
		     const bool ylog,
		    const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm);
  void plotDataVSMC(int categoria,int channel,int istogramma,
                    TH1D* data, TH1D** bkg,
                    const TString* names, const unsigned nHist,
                    const TString& name_cut,const TString& name_channel, const TString& name_histo,
                    const bool ylog,
                    const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm,
		   const int lumi_year = 0);
  void plotDataVSMC_SR(//int categoria,
		       int channel,
		       TH1D* plot_variation[3],  
		       //const TString& name_cut,
		       const TString& name_channel, const TString& name_histo
		       //, const unsigned widthopt
		       );
#endif
