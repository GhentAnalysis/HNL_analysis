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

void plotDataVSMC(int categoria, int istogramma,TH1D* data, TH1D** bkg, const TString* names, const unsigned nHist, const TString& file,const TString& file2, const bool ylog, const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm){
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

	//if (categoria <= 6 ) TLegend* legend = new TLegend(0.5,0.79,0.95,0.88,NULL,"brNDC");
	//	if (categoria > 6 ) TLegend* legend = new TLegend(0.65,0.79,0.95,0.88,NULL,"brNDC");
    legend->SetFillStyle(0);
	//Add data to the legend
    // legend->AddEntry(data,names[0]);
	//Add signal to the legenD
	if(plotsig){
		//unsigned sigstyles[5] = {
		for(unsigned sig = 0; sig < nSig; ++sig){
		  if (sig != 4 && sig != 21 && sig != 59 && sig != 43) continue;
			signal[sig]->SetLineColor(sigCols[sig]);
			signal[sig]->SetLineColor(sigCols[sig]);
			signal[sig]->SetMarkerColor(sigCols[sig]);
			signal[sig]->SetLineWidth(3);
			signal[3] ->SetLineColor(1);
			signal[1] ->SetLineColor(2);

			signal[3]->SetMarkerColor(1);
			signal[1]->SetMarkerColor(2);

			signal[sig]->SetLineWidth(4);
		
			//legend->AddEntry(signal[sig], signames[sig]);
			//if (sig == 4) legend->AddEntry(signal[sig], "m_{N}=2 |V|^{2}=3*10^{-5}");
			//if (sig == 21) legend->AddEntry(signal[sig], "m_{N}=4 |V|^{2}==6*10^{-6}");
			//if (sig == 59) legend->AddEntry(signal[sig], "m_{N}=6 |V|^{2}==4.1*10^{-6}");
			//if (sig == 43) legend->AddEntry(signal[sig], "m_{N}=10 |V|^{2}==6*10^{-6}");


			if (sig == 4)  signal[sig]->SetLineColor(1);
			if (sig == 21) signal[sig]->SetLineColor(kAzure +10);
			if (sig == 59) signal[sig]->SetLineColor(kMagenta);
			if (sig == 43) signal[sig]->SetLineColor(kGreen);
			if (sig == 4)  signal[sig]->SetMarkerColor(1);
			if (sig == 21) signal[sig]->SetMarkerColor(kAzure +10);
			if (sig == 59) signal[sig]->SetMarkerColor(kMagenta);
			if (sig == 43) signal[sig]->SetMarkerColor(kGreen);

		}
	}


	legend->AddEntry(signal[4], "M2 |V|^{2}=3*10^{-5}");
	legend->AddEntry(signal[21], "M4 |V|^{2}=6*10^{-6}");
	legend->AddEntry(signal[59], "M6 |V|^{2}=4.1*10^{-6}");
	legend->AddEntry(signal[43], "M10 |V|^{2}=6*10^{-6}");
	
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
    if(bkgTot->GetBinContent(bkgTot->GetMaximumBin()) > signal[59]->GetBinContent(signal[59]->GetMaximumBin()) ){
      if(!ylog) signal[59]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*1.5);
      else signal[59]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*10);
    }
    else{
      if(!ylog) signal[59]->SetMaximum(signal[59]->GetBinContent(signal[59]->GetMaximumBin())*1.5);
      else signal[59]->SetMaximum(signal[59]->GetBinContent(signal[59]->GetMaximumBin())*10);
    }

    
   

    signal[59]->GetYaxis()->SetTitle("Events/search region");
    signal[59]->GetXaxis()->SetTitle("");

    //Draw signal plots
    if(plotsig){
      for(unsigned sig = 0; sig < nSig; ++sig){
	if (sig == 59){
	  if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	  signal[sig]->SetMinimum(0);
	  signal[sig]->SetMinimum(0.1);
	  signal[sig]->Draw("histe ");
	}
      }
    }
    cout<<"---"<<endl;
    for (int i = 0; i < 22; i ++){

      //cout<<"bgk         bin# "<<i+1<<") bgk: "<<((TH1*)(bkgStack->GetStack()->Last()))->  GetBinContent (i+1)<<endl;
    }
    //data->Draw("pe");	//The range used is now that of the data histogra
        bkgStack->Draw("hist same ");
	
    //data->Draw("pe same");
    legend->Draw("same");
   
    if (categoria == 0 && istogramma == 1){
      //  t_zeta2->Draw();
      // t_zeta1->Draw();
      //  ld2->Draw("lsame");
    }
    bkgTot->Draw("e2same");
    if(plotsig){
      for(unsigned sig = 0; sig < nSig; ++sig){
	if (sig == 4 || sig == 21 || sig == 43 || sig == 59){
	  if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	  signal[sig]->SetMinimum(0.1);
	  signal[sig]->Draw("histe same");
	}
      }
    }

	//redraw axis over histograms
    gPad->RedrawAxis();




    //CMS_lumi(c,"Preliminary", true);
	drawLumi(p1);
 	/*	  c->cd();
	//Make ratio plot in second pad
    p2 = new TPad(file + "2","",0,0.0,1,xPad);
    p2->Draw();
    p2->cd();
    p2->SetTopMargin(0);
    p2->SetBottomMargin(0.4);
	
	TH1D* dataC = (TH1D*) data->Clone();
	TH1D* bkgTotC = (TH1D*) bkgTot->Clone();

    dataC->Divide(bkgTotC);
    dataC->SetMarkerColor(1);
    dataC->SetLineColor(1);
    dataC->GetYaxis()->SetRangeUser(0.,1.999);
    dataC->GetYaxis()->SetTitle("obs/pred");
    dataC->GetYaxis()->SetTitleOffset(0.9/((1.-xPad)/xPad));
    dataC->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.06
    dataC->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.09
    dataC->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.05
    dataC->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.075

    dataC->Draw("pe");
	//Draw line at 1 on ratio plot
    double xmax = dataC->GetBinCenter(dataC->GetNbinsX()) + dataC->GetBinWidth(dataC->GetNbinsX())/2;
    double xmin = dataC->GetBinCenter(0) + dataC->GetBinWidth(0)/2;
    TLine *line = new TLine(xmin, 1, xmax, 1);
    line->SetLineStyle(2);
    line->Draw("same");
	*/
	//c->SaveAs("plots_pdf/" + file2 + "_"+file + ".pdf");
    c->SaveAs("plots_root/" + file2 + "_"+file + ".root");

}


