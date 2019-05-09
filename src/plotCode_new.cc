//include c++ library classes
#include <fstream>
#include <set>
//include Root classes
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
//Include other parts of the code
//#include "MultilepSUSYfunc.h"
#include "plotCode_new.h"
#include "drawLumi.h"

extern const double xPad;
extern const Color_t colors[];
//extern const Color_t sigcolors[];


void histcol(TH1D *h, const Color_t color){
	h->SetLineColor(color);
	h->SetMarkerColor(color);
	h->SetLineWidth(2);
}

TH1D *HistDiv(TH1D *h1, TH1D *h2, const bool abs){
	TH1D *h1c = (TH1D*) h1->Clone();
	TH1D *h2c = (TH1D*) h2->Clone();
	if(!abs){
		h1c->Scale(1/h1c->Integral(), "width");
		h2c->Scale(1/h2c->Integral(), "width");
	}
	h1c->Divide(h2c);
	return h1c;
}

void HistLabelSizes(TH1D *h, const double xlabel, const double xtitle, const double ylabel, const double ytitle){
  	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleFont(42);
	h->GetXaxis()->SetLabelSize(xlabel);
	h->GetXaxis()->SetTitleSize(xtitle);
	h->GetYaxis()->SetLabelSize(ylabel);
	h->GetYaxis()->SetTitleSize(ytitle);
}

void HistLabelSizes(TH2D *h, const double xlabel, const double xtitle, const double ylabel, const double ytitle){
  	h->GetXaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleFont(42);
	h->GetXaxis()->SetLabelSize(xlabel);
	h->GetXaxis()->SetTitleSize(xtitle);
	h->GetYaxis()->SetLabelSize(ylabel);
	h->GetYaxis()->SetTitleSize(ytitle);
}

void StackCol(TH1D *h, const Color_t color){
	histcol(h,color);
	h->SetFillColor(color);
	h->SetLineWidth(1);
	h->SetLineColor(color);
}

void yieldOrder(TH1D**& hists, unsigned* histInd, const unsigned nHist){
	unsigned ordered[nHist];
	for(unsigned h = 0; h < nHist; ++h) ordered[h] = 999;
	for(unsigned h = 0; h < nHist; ++h){
		unsigned maxH = 999;
		double maxYield = -9999.;
		for(unsigned k = 0; k <nHist; ++k){
			bool found = false;
			for(unsigned i = 0; i < nHist; ++i){
				if(ordered[i] == k){
					found = true;
					break;
				}
			}
			if(!found){
				double yield = hists[k]->GetSumOfWeights();
				if(yield > maxYield){
					maxYield = yield;
					maxH = k;
				}
			}
		}
		//ordered[h] = maxH;
		ordered[h] = h;

	}
	TH1D* histC[nHist];
	for(unsigned h = 0; h < nHist; ++h){
		histC[h] = (TH1D*) hists[ordered[h]]->Clone();
		histInd[h] = ordered[h];
	}
	for(unsigned h = 0; h < nHist; ++h){
		hists[h] = (TH1D*) histC[h]->Clone();
	}
}

void plotDataVSMC_mu(int categoria,int channel,int istogramma,
		     TH1D* data, TH1D** bkg,
		     const TString* names, const unsigned nHist,
		     const TString& name_cut,const TString& name_channel, const TString& name_histo,
		     const bool ylog,
		     const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm){
  
	//Order background histograms in terms of yields
	unsigned histI[nHist];
		yieldOrder(bkg, histI, nHist);
	//Calculate total Bkg yields
	TH1D* bkgTot = (TH1D*) bkg[0]->Clone();
	for(unsigned int i = 1; i <  nHist; ++i){
		bkgTot->Add(bkg[i]);
	}
	//Make a stack containing all backgrounds
	THStack* bkgStack;
	bkgStack = new THStack("bkgStack", "bkgStack");	
	for(int effsam = nHist -1; effsam > -1 ; --effsam){
		StackCol(bkg[effsam], colors[effsam]);
		bkgStack->Add(bkg[effsam], "f");
	}

	if(signames==nullptr) {} // dummy, just to avoid warning

	//Make a legend for data and all backgrounds
	TLegend* legend = new TLegend(0.81,0.15,0.99,0.88,NULL,"brNDC");

	legend->SetFillStyle(0);
	
	//Add signal to the legenD
	if(plotsig){
	  for(unsigned sig = 0; sig < nSig; ++sig){		  
	    signal[sig]->SetLineColor(sigCols[sig]);
	    signal[sig]->SetLineColor(sigCols[sig]);
	    signal[sig]->SetMarkerColor(sigCols[sig]);
	    signal[sig]->SetLineWidth(4);
	    legend->AddEntry(signal[sig], signames[sig]);
	  }
	}
	
	for(int effsam = nHist - 1; effsam > -1; --effsam){
	  legend->AddEntry(bkg[effsam], names[histI[effsam] + 1 + nSig]);
	  legend->	 SetNColumns(1);
	}
	//Make canvas and pads for plotting
	double width, height;
	if(widthopt == 0){
		width = 500;
		height = 600;
	} else if(widthopt == 1){
		width = 2000;
		height = 500;
	} else if(widthopt == 2){
		width = 700;
		height = 600;
	} else{
		std::cerr << "Incorrect width option given can't make plot" << std::endl;
		return;
	}
    TCanvas *c =  new TCanvas(file,"",width*(1-xPad),height);   //1000/500
    c->cd();
	

    TPad* p1, *p2;
	//Plot data and MC yields in first pad
    p1 = new TPad(file,"",0,xPad,1,1);
    p1->Draw();
    p1->cd();
    p1->SetTopMargin(0.1);//0.1*(width*(1-xPad)/650)  CHANGE THIS BACK
    p1->SetBottomMargin(0.15);
    p1->SetRightMargin(0.2);

    bkgTot->SetFillStyle(3005);
    bkgTot->SetFillColor(kGray+2);
    bkgTot->SetMarkerStyle(1);
    data->SetMinimum(0.1);
    bkgTot->SetMinimum(0.1);
    bkgStack->SetMinimum(0.1);

    if(ylog) p1->SetLogy();
    HistLabelSizes(data,0.1,0.1,0.07,0.07);
    
    //Determine the maximum range of the histogram, depending on the maximum range of the bkg or data
    if(bkgTot->GetBinContent(bkgTot->GetMaximumBin()) > signal[3]->GetBinContent(signal[3]->GetMaximumBin()) ){
      if(!ylog) signal[3]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*1.5);
      else signal[3]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*10);
    }
    else{
      if(!ylog) signal[3]->SetMaximum(signal[3]->GetBinContent(signal[3]->GetMaximumBin())*1.5);
      else signal[3]->SetMaximum(signal[3]->GetBinContent(signal[3]->GetMaximumBin())*10);
    }

   
    //Draw signal plots
    if(plotsig){
      for(unsigned sig = 0; sig < nSig; ++sig){
	if (sig == 3){
	  if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	  signal[sig]->SetMinimum(0);
	  signal[sig]->SetMinimum(0.1);
	  signal[sig]->Draw("histe ");
	}
      }
    }
   
    bkgStack->Draw("hist same ");
    legend->Draw("same");
    bkgTot->Draw("e2same");
    if(plotsig){
      for(unsigned sig = 0; sig < nSig; ++sig){
	if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	signal[sig]->SetMinimum(0.1);
	signal[sig]->Draw("histe same");
      }
    }
    //redraw axis over histograms
    gPad->RedrawAxis();
    //CMS_lumi(c,"Preliminary", true);
    drawLumi(p1);
    if (channel == 0){
      c->SaveAs("plots_pdf/" +  "mmm" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "mmm" + "/"+ name_cut + "/"+ name_histo ".root");
    }
    if (channel == 1){
      c->SaveAs("plots_pdf/" +  "mmeOS" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "mmeOS" + "/"+ name_cut + "/"+ name_histo ".root");
    }
    if (channel == 2){
      c->SaveAs("plots_pdf/" +  "mmeSS" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "mmeSS" + "/"+ name_cut + "/"+ name_histo ".root");
    }
    if (channel == 6){
      c->SaveAs("plots_pdf/" +  "mu" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "mu" + "/"+ name_cut + "/"+ name_histo ".root");
    }
}

void plotDataVSMC_e(int categoria,int channel,int istogramma,
		     TH1D* data, TH1D** bkg,
		     const TString* names, const unsigned nHist,
		     const TString& name_cut,const TString& name_channel, const TString& name_histo,
		     const bool ylog,
		     const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm){
  
	//Order background histograms in terms of yields
	unsigned histI[nHist];
		yieldOrder(bkg, histI, nHist);
	//Calculate total Bkg yields
	TH1D* bkgTot = (TH1D*) bkg[0]->Clone();
	for(unsigned int i = 1; i <  nHist; ++i){
		bkgTot->Add(bkg[i]);
	}
	//Make a stack containing all backgrounds
	THStack* bkgStack;
	bkgStack = new THStack("bkgStack", "bkgStack");	
	for(int effsam = nHist -1; effsam > -1 ; --effsam){
		StackCol(bkg[effsam], colors[effsam]);
		bkgStack->Add(bkg[effsam], "f");
	}

	if(signames==nullptr) {} // dummy, just to avoid warning

	//Make a legend for data and all backgrounds
	TLegend* legend = new TLegend(0.81,0.15,0.99,0.88,NULL,"brNDC");

	legend->SetFillStyle(0);
	
	//Add signal to the legenD
	if(plotsig){
	  for(unsigned sig = 0; sig < nSig; ++sig){		  
	    signal[sig]->SetLineColor(sigCols[sig]);
	    signal[sig]->SetLineColor(sigCols[sig]);
	    signal[sig]->SetMarkerColor(sigCols[sig]);
	    signal[sig]->SetLineWidth(4);
	    legend->AddEntry(signal[sig], signames[sig]);
	  }
	}
	
	for(int effsam = nHist - 1; effsam > -1; --effsam){
	  legend->AddEntry(bkg[effsam], names[histI[effsam] + 1 + nSig]);
	  legend->	 SetNColumns(1);
	}
	//Make canvas and pads for plotting
	double width, height;
	if(widthopt == 0){
		width = 500;
		height = 600;
	} else if(widthopt == 1){
		width = 2000;
		height = 500;
	} else if(widthopt == 2){
		width = 700;
		height = 600;
	} else{
		std::cerr << "Incorrect width option given can't make plot" << std::endl;
		return;
	}
    TCanvas *c =  new TCanvas(file,"",width*(1-xPad),height);   //1000/500
    c->cd();
	

    TPad* p1, *p2;
	//Plot data and MC yields in first pad
    p1 = new TPad(file,"",0,xPad,1,1);
    p1->Draw();
    p1->cd();
    p1->SetTopMargin(0.1);//0.1*(width*(1-xPad)/650)  CHANGE THIS BACK
    p1->SetBottomMargin(0.15);
    p1->SetRightMargin(0.2);

    bkgTot->SetFillStyle(3005);
    bkgTot->SetFillColor(kGray+2);
    bkgTot->SetMarkerStyle(1);
    data->SetMinimum(0.1);
    bkgTot->SetMinimum(0.1);
    bkgStack->SetMinimum(0.1);

    if(ylog) p1->SetLogy();
    HistLabelSizes(data,0.1,0.1,0.07,0.07);
    
    //Determine the maximum range of the histogram, depending on the maximum range of the bkg or data
    if(bkgTot->GetBinContent(bkgTot->GetMaximumBin()) > signal[3]->GetBinContent(signal[3]->GetMaximumBin()) ){
      if(!ylog) signal[3]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*1.5);
      else signal[3]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*10);
    }
    else{
      if(!ylog) signal[3]->SetMaximum(signal[3]->GetBinContent(signal[3]->GetMaximumBin())*1.5);
      else signal[3]->SetMaximum(signal[3]->GetBinContent(signal[3]->GetMaximumBin())*10);
    }

   
    //Draw signal plots
    if(plotsig){
      for(unsigned sig = 0; sig < nSig; ++sig){
	if (sig == 3){
	  if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	  signal[sig]->SetMinimum(0);
	  signal[sig]->SetMinimum(0.1);
	  signal[sig]->Draw("histe ");
	}
      }
    }
   
    bkgStack->Draw("hist same ");
    legend->Draw("same");
    bkgTot->Draw("e2same");
    if(plotsig){
      for(unsigned sig = 0; sig < nSig; ++sig){
	if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	signal[sig]->SetMinimum(0.1);
	signal[sig]->Draw("histe same");
      }
    }
    //redraw axis over histograms
    gPad->RedrawAxis();
    //CMS_lumi(c,"Preliminary", true);
    drawLumi(p1);
    if (channel == 3){
      c->SaveAs("plots_pdf/" +  "eee" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "eee" + "/"+ name_cut + "/"+ name_histo ".root");
    }
    if (channel == 4){
      c->SaveAs("plots_pdf/" +  "eemOS" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "eeeOS" + "/"+ name_cut + "/"+ name_histo ".root");
    }
    if (channel == 5){
      c->SaveAs("plots_pdf/" +  "eeeSS" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "eeeSS" + "/"+ name_cut + "/"+ name_histo ".root");
    }
    if (channel == 7){
      c->SaveAs("plots_pdf/" +  "e" + "/"+ name_cut + "/"+ name_histo ".pdf");
      c->SaveAs("plots_root/" + "e" + "/"+ name_cut + "/"+ name_histo ".root");
    }
}

