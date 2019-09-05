//include c++ library classes
#include <fstream>
#include <set>
//include Root classes
#include "TCanvas.h"
#include "TColor.h"
#include "TLatex.h"
#include "TLine.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
//Include other parts of the code
//#include "MultilepSUSYfunc.h"
#include "../interface/plotCode_new.h"
#include "../interface/drawLumi.h"

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
    //unsigned maxH = 999;
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
	  //maxH = k;
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
  // Dummy: to be removed!
  if(categoria==-3672 && istogramma==-3276) std::cout << name_channel.Data() << std::endl;

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
  TLegend* legend = new TLegend(0.2,0.79,0.95,0.88,NULL,"brNDC");
    
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
    legend->     SetNColumns(1);
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
  TCanvas *c =  new TCanvas(name_histo,"",width*(1-xPad),height);   //1000/500
  c->cd();
    
    
  TPad *p1; //, *p2;
  //Plot data and MC yields in first pad
  p1 = new TPad(name_histo,"",0,xPad,1,1);
  p1->Draw();
  p1->cd();
  p1->SetTopMargin(0.1);//0.1*(width*(1-xPad)/650)  CHANGE THIS BACK
  p1->SetBottomMargin(0);
  bkgTot->SetFillStyle(3005);
  bkgTot->SetFillColor(kGray+2);
  bkgTot->SetMarkerStyle(1);
  data->SetMinimum(0.1);
  bkgTot->SetMinimum(0.1);
  bkgStack->SetMinimum(0.1);
  if(!ylog) data->SetMinimum(0.01);
  else if(ylog) p1->SetLogy();
    
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
    c->SaveAs("plots_pdf/mmm/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mmm/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 1){
    c->SaveAs("plots_pdf/mmeOS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mmeOS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 2){
    c->SaveAs("plots_pdf/mmeSS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mmeSS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 6){
    c->SaveAs("plots_pdf/mu/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mu/"+ name_cut + "/"+ name_histo + ".root");
  }
}

void plotDataVSMC_e(int categoria,int channel,int istogramma,
                    TH1D* data, TH1D** bkg,
                    const TString* names, const unsigned nHist,
                    const TString& name_cut,const TString& name_channel, const TString& name_histo,
                    const bool ylog,
                    const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm){
  // Dummy: to be removed!
  if(categoria==-3672 && istogramma==-3276) std::cout << name_channel.Data() << std::endl;
    
  //Order background histograms in terms of yields
  unsigned histI[nHist];
  yieldOrder(bkg, histI, nHist);
  //Calculate total Bkg yields
  TH1D* bkgTot = (TH1D*) bkg[0]->Clone();
  for(unsigned int i = 1; i <  nHist; ++i){
    std::cout<<i<<")  "<< names[i]<< "  "<< name_histo[i]<<std::endl;
    bkgTot->Add(bkg[i]);
  }
  //Make a stack containing all backgrounds
  THStack* bkgStack;
  bkgStack = new THStack("bkgStack", "bkgStack");
  for(int effsam = nHist -1; effsam > -1 ; --effsam){
    StackCol(bkg[effsam], colors[effsam]);
    std::cout<<effsam<<" "<< names[effsam]<<std::endl;
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
    legend->     SetNColumns(1);
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
  TCanvas *c =  new TCanvas(name_histo,"",width*(1-xPad),height);   //1000/500
  c->cd();
    
    
  TPad *p1; //, *p2;
  //Plot data and MC yields in first pad
  p1 = new TPad(name_histo,"",0,xPad,1,1);
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
    c->SaveAs("plots_pdf/eee/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/eee/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 4){
    c->SaveAs("plots_pdf/eemOS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/eemOS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 5){
    c->SaveAs("plots_pdf/eemSS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/eemSS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 7){
    c->SaveAs("plots_pdf/e/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/e/"+ name_cut + "/"+ name_histo + ".root");
  }
}



void plotDataVSMC(int categoria,int channel,int istogramma,
		  TH1D* data, TH1D** bkg,
		  const TString* names, const unsigned nHist,
		  const TString& name_cut,const TString& name_channel, const TString& name_histo,
		  const bool ylog,
		  const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm){
  // Dummy: to be removed!
  if(categoria==-3672) std::cout << name_channel.Data() << std::endl;
    
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
    if (names[histI[effsam] + 1 + nSig] == "nonprompt DF" ) bkg[effsam]->SetFillStyle(3020); 
    bkgStack->Add(bkg[effsam], "f");
  }
    
  if(signames==nullptr) {} // dummy, just to avoid warning
    
  //Make a legend for data and all backgrounds
  TLegend* legend = new TLegend(0.15,0.72,0.95,0.88,NULL,"brNDC");
    
  legend->SetFillStyle(0);
    
  //Add signal to the legenD
  if(plotsig){
    for(unsigned sig = 0; sig < nSig; ++sig){
      if (sig == 1 || sig == 4 || sig == 6 || sig == 8 || sig == 2 || sig == 11 || sig == 14 || sig == 16 || sig == 18 || sig == 12 ) continue;
	  
      if ((channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 )  && sig >= 10) continue; 
      if ((channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ) && sig < 10) continue; 
      signal[sig]->SetLineColor(sigCols[sig]);
      signal[sig]->SetMarkerColor(sigCols[sig]);
      signal[sig]->SetLineWidth(4);
	
      if (sig >= 10) {
	signal[sig]->SetLineColor(sigCols[sig-10]);
	signal[sig]->SetMarkerColor(sigCols[sig-10]);
	signal[sig]->SetLineWidth(4);
      }
      legend->AddEntry(signal[sig], signames[sig]);

      //if (sig == 0 ||sig == 10 ) legend->AddEntry(signal[sig], "m_{N}=1 GeV");
      //if (sig == 2 ||sig == 12 ) legend->AddEntry(signal[sig], "m_{N}=2 GeV");
      //if (sig == 4 ||sig == 14 ) legend->AddEntry(signal[sig], "m_{N}=4 GeV");
      //if (sig == 6 ||sig == 16 ) legend->AddEntry(signal[sig], "m_{N}=6 GeV");
      //if (sig == 8 ||sig == 18 ) legend->AddEntry(signal[sig], "m_{N}=10 GeV");
		 
    }
  }


  /* legend->AddEntry(signal[0], "M = 1 GeV");
     legend->AddEntry(signal[2], "M = 2 GeV");
     legend->AddEntry(signal[4], "M = 4 GeV");
     legend->AddEntry(signal[6], "M = 6 GeV");
     legend->AddEntry(signal[8], "M = 10 GeV");
     legend->AddEntry(signal[10], "M = 1 GeV");
     legend->AddEntry(signal[12], "M = 2 GeV");
     legend->AddEntry(signal[14], "M = 4 GeV");
     legend->AddEntry(signal[16], "M = 6 GeV");
     legend->AddEntry(signal[18], "M = 10 GeV");*/

    
  for(int effsam = nHist - 1; effsam > -1; --effsam){
	  if (names[histI[effsam] + 1 + nSig] == "Xgamma") legend->AddEntry(bkg[effsam], "X+#gamma");
	  else if  (names[histI[effsam] + 1 + nSig] == "XTTX") legend->AddEntry(bkg[effsam], "TT/T+X");
	  else if (names[histI[effsam] + 1 + nSig] == "TTX") legend->AddEntry(bkg[effsam], "TT/T+X");
	  else if (names[histI[effsam] + 1 + nSig] == "ttbar") legend->AddEntry(bkg[effsam], "t#bar{t}");
	  else if (names[histI[effsam] + 1 + nSig] == "WJets") legend->AddEntry(bkg[effsam], "W+jets");
    	  else legend->AddEntry(bkg[effsam], names[histI[effsam] + 1 + nSig]);
          legend->     SetNColumns(5);
  }

  if (istogramma == 1 ){
    signal[0]->SetStats(0);
    signal[0]-> GetXaxis()->LabelsOption("vu");
    signal[0]-> GetXaxis()->SetBinLabel(1, "3 leptons");
    signal[0]-> GetXaxis()->SetBinLabel(2, "#DeltaR (l2-l3)");
    signal[0]-> GetXaxis()->SetBinLabel(3, "bveto");
    signal[0]-> GetXaxis()->SetBinLabel(4, "M_{lll}");
    signal[0]-> GetXaxis()->SetBinLabel(5, "min  #Delta #phi");
    signal[0]-> GetXaxis()->SetBinLabel(6, "cos(SV,HNL)");
    signal[0]-> GetXaxis()->SetBinLabel(7, "M_{ll} (l2-l3)");

    signal[10]->SetStats(0);
    signal[10]-> GetXaxis()->LabelsOption("vu");
    signal[10]-> GetXaxis()->SetBinLabel(1, "3 leptons");
    signal[10]-> GetXaxis()->SetBinLabel(2, "#DeltaR (l2-l3)");
    signal[10]-> GetXaxis()->SetBinLabel(3, "bveto");
    signal[10]-> GetXaxis()->SetBinLabel(4, "M_{lll}");
    signal[10]-> GetXaxis()->SetBinLabel(5, "min  #Delta #phi");
    signal[10]-> GetXaxis()->SetBinLabel(6, "cos(SV,HNL)");
    signal[10]-> GetXaxis()->SetBinLabel(7, "M_{ll} (l2-l3)");

    signal[0]-> GetXaxis()->SetLabelSize(0.045);
    signal[0]-> GetXaxis()->SetLabelOffset(0.01);
    signal[10]-> GetXaxis()->SetLabelSize(0.045);
    signal[10]-> GetXaxis()->SetLabelOffset(0.01);
  }

  // isotgramma delle SR --> linee e roba varia
  if (istogramma == 0 ){
    signal[0]->SetStats(0);
    signal[0]-> GetXaxis()->LabelsOption("hu");
    signal[0]-> GetXaxis()->SetTitle ("M_{l_{2}l_{3}} (GeV)");	  
   signal[0]->GetXaxis()->SetTitleSize(0.06);
    signal[0]->GetXaxis()->SetTitleOffset(1);
    signal[0]->GetXaxis() ->SetTitleFont(132);
    signal[0]-> GetXaxis()->SetBinLabel(1, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(2, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(3, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(4, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(5, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(6, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(7, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(8, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(9, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(10, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(11, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(12, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(13, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(14, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(15, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(16, "> 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(17, "< 4 ");
    signal[0]-> GetXaxis()->SetBinLabel(18, "> 4 ");
 
    signal[10]->SetStats(0);
    signal[10]-> GetXaxis()->LabelsOption("hu");
	signal[10]-> GetXaxis()->SetTitle ("M_{l_{2}l_{3}} (GeV)");	  
   signal[10]->GetXaxis()->SetTitleSize(0.06);
    signal[10]->GetXaxis()->SetTitleOffset(1);
    signal[10]->GetXaxis() ->SetTitleFont(132);
    signal[10]-> GetXaxis()->SetBinLabel(1, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(2, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(3, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(4, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(5, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(6, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(7, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(8, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(9, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(10, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(11, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(12, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(13, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(14, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(15, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(16, "> 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(17, "< 4 ");
    signal[10]-> GetXaxis()->SetBinLabel(18, "> 4 ");

    signal[0]-> GetXaxis()->SetLabelSize(0.02);
    signal[0]-> GetXaxis()->SetLabelOffset(0.01);
    signal[10]-> GetXaxis()->SetLabelSize(0.02);
    signal[10]-> GetXaxis()->SetLabelOffset(0.01);
	    
	    
	    
	    
	    
	    
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
  if (istogramma == 0) width = 800;
  if (istogramma == 0) height = 500;

  TCanvas *c =  new TCanvas(name_histo,"",width*(1-xPad),height);   //1000/500
  c->cd();

  TPad *p1; //, *p2;
  //Plot data and MC yields in first pad
  p1 = new TPad(name_histo,"",0,xPad,1,1);
  p1->Draw();
  p1->cd();
  p1->SetTopMargin(0.1);//0.1*(width*(1-xPad)/650)  CHANGE THIS BACK
  p1->SetBottomMargin(0.15);
  bkgTot->SetFillStyle(3005);
  bkgTot->SetFillColor(kGray+2);
  bkgTot->SetMarkerStyle(1);
  data->SetMinimum(0.1);
  bkgTot->SetMinimum(0.1);
  bkgStack->SetMinimum(0.1);
  if(!ylog) data->SetMinimum(0.01);
  else if(ylog) p1->SetLogy();
    
  if(ylog) p1->SetLogy();
  HistLabelSizes(data,0.1,0.1,0.07,0.07);
    
  //double scaling = 30;
  //Determine the maximum range of the histogram, depending on the maximum range of the bkg or data
  if (channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 ){
    if(bkgTot->GetBinContent(bkgTot->GetMaximumBin()) > signal[0]->GetBinContent(signal[0]->GetMaximumBin()) ){
      if(!ylog) signal[0]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*1.5);
      if(!ylog && istogramma ==0) signal[0]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*20);

      else {
	if (istogramma !=0)signal[0]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*30);
	if (istogramma ==0)signal[0]->SetMaximum(4000000);

      }
    }
    else{
      if(!ylog) signal[0]->SetMaximum(signal[0]->GetBinContent(signal[0]->GetMaximumBin())*1.5);
       if(!ylog && istogramma ==0) signal[0]->SetMaximum(signal[0]->GetBinContent(bkgTot->GetMaximumBin())*20);

      else {
	if (istogramma !=0)signal[0]->SetMaximum(signal[0]->GetBinContent(signal[0]->GetMaximumBin())*30);
	if (istogramma ==0)signal[0]->SetMaximum(4000000);

      } 
    }
  }
  if (channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ){
    if(bkgTot->GetBinContent(bkgTot->GetMaximumBin()) > signal[10]->GetBinContent(signal[10]->GetMaximumBin()) ){
      if(!ylog) signal[10]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*1.5);
      if(!ylog && istogramma ==0) signal[10]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*20);

      else {
	if (istogramma !=0)signal[10]->SetMaximum(bkgTot->GetBinContent(bkgTot->GetMaximumBin())*30);
	if (istogramma ==0)signal[10]->SetMaximum(4000000);

      }
    }
    else{
      if(!ylog) signal[10]->SetMaximum(signal[10]->GetBinContent(signal[10]->GetMaximumBin())*1.5);
      if(!ylog && istogramma ==0) signal[10]->SetMaximum(signal[10]->GetBinContent(bkgTot->GetMaximumBin())*20);
	    
      else {
	if (istogramma !=0)signal[10]->SetMaximum(signal[10]->GetBinContent(signal[10]->GetMaximumBin())*30);
	if (istogramma ==0)signal[10]->SetMaximum(4000000);

      }
    }
  }
    
  //Draw signal plots
  if(plotsig){
    for(unsigned sig = 0; sig < nSig; ++sig){

      if ((channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 ) && sig == 0){
	if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	signal[sig]->SetMinimum(0);
	signal[sig]->SetMinimum(0.1);
	signal[sig]->Draw("histe ");
      }
      if ((channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ) && sig == 10){
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
      if (sig == 1 || sig == 4 || sig == 6 || sig == 8 || sig == 2 || sig == 11 || sig == 14 || sig == 16 || sig == 18 || sig == 12 ) continue;
      if ((channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 )  && sig >= 10) continue; 
      if ((channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ) && sig < 10) continue; 

      if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
      signal[sig]->SetMinimum(0.1);
      signal[sig]->Draw("histe same");
    }
  }
	
  if (istogramma == 0){
    Int_t ci;      // for color index setting
    //TColor *color; // for color definition with alpha
	  //    line = new TLine(12.5,0.02,12.5,50000);

	  
    TLine *line = new TLine(6.5,0.02,6.5,50000);
    line->SetLineWidth(5);
    line->Draw();
    line = new TLine(12.5,0.02,12.5,50000);
    line->SetLineWidth(5);
    line->Draw();
	  //100000
    line = new TLine(2.5,20000,2.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(3);
    line->Draw();
    line = new TLine(4.5,20000,4.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(3);
    line->Draw();
    line = new TLine(8.5,20000,8.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(3);
    line->Draw();
    line = new TLine(10.5,20000,10.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(3);
    line->Draw();
    line = new TLine(14.5,20000,14.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(3);
    line->Draw();
    line = new TLine(16.5,20000,16.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(3);
    line->Draw();
    TLatex *    tex = new TLatex(0.8748578,17546.74,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.031586,21625.48,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.25285,33424.77,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.492553,48185.32,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.66772,54433.15,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.833668,62571.47,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.907423,91788.17,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(2.598872,137012.8,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.584746,8013.302,"");
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(1.12378,4134.047,"");
    tex->SetTextSize(0.003034901);
    tex->SetLineWidth(2);
    tex->Draw();
	//7000.46	
    tex = new TLatex(0.9578318,7000.46,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(2.857013,7000.46,"#DeltaR [2,10]cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(4.912923,7000.46,"#DeltaR > 10cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.9578318,7000.46,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
		
    tex = new TLatex(6.9578318,7000.46,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(8.857013,7000.46,"#DeltaR [2,10]cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(10.912923,7000.46,"#DeltaR > 10cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();	
		
    tex = new TLatex(12.9578318,7000.46,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(14.857013,7000.46,"#DeltaR [2,10]cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(16.912923,7000.46,"#DeltaR > 10cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();	
		
	//25000.83	
    if (channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 ){	
	 //         tex = new TLatex(2.857013,signal[0]->GetBinContent(signal[0]->GetMaximumBin())* 8,"#mu#mu#mu");

      tex = new TLatex(2.857013,25000.83,"#mu#mu#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.07);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(8.857013,25000.83,"#mu^{#pm}#mu^{#mp}e");
      tex->SetTextColor(1);
      tex->SetTextSize(0.07);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(14.857013,25000.83,"#mu^{#pm}#mu^{#pm}e");
      tex->SetTextColor(1);
      tex->SetTextSize(0.07);
      tex->SetLineWidth(2);
      tex->Draw();
    }	
    if (channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ){	
      tex = new TLatex(2.857013,25000.83,"eee");
      tex->SetTextColor(1);
      tex->SetTextSize(0.07);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(8,25000.83,"e^{#pm}e^{#mp}#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.07);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(13.857013,25000.83,"e^{#pm}e^{#pm}#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.07);
      tex->SetLineWidth(2);
      tex->Draw();
    }	
  }
	
	
	
	
	
	
  //redraw axis over histograms
  gPad->RedrawAxis();
  //CMS_lumi(c,"Preliminary", true);
  drawLumi(p1);
  if (channel == 3){
    c->SaveAs("plots_pdf/eee/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/eee/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 4){
    c->SaveAs("plots_pdf/eemOS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/eemOS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 5){
    c->SaveAs("plots_pdf/eemSS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/eemSS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 7){
    c->SaveAs("plots_pdf/e/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/e/"+ name_cut + "/"+ name_histo + ".root");
  }

  if (channel == 0){
    c->SaveAs("plots_pdf/mmm/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mmm/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 1){
    c->SaveAs("plots_pdf/mmeOS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mmeOS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 2){
    c->SaveAs("plots_pdf/mmeSS/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mmeSS/"+ name_cut + "/"+ name_histo + ".root");
  }
  if (channel == 6){
    c->SaveAs("plots_pdf/mu/"+ name_cut + "/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mu/"+ name_cut + "/"+ name_histo + ".root");
  }
}



void plotDataVSMC_SR(int categoria,int channel,
		  TH1D* data_nominal, TH1D* data_down, TH1D* data_up, 
		  const TString& name_cut,const TString& name_channel, const TString& name_histo,
		  const unsigned widthopt){
    
      
  //Make a legend for data and all backgrounds
  TLegend* legend = new TLegend(0.15,0.72,0.95,0.88,NULL,"brNDC");    
  legend->SetFillStyle(0);   
  legend->AddEntry(data_nominal, name_channel+"_nominal");
  legend->AddEntry(data_down, name_channel+"_down");
  legend->AddEntry(data_up, name_channel+"_up");

  
 

  // isotgramma delle SR --> linee e roba varia
  
    data_nominal->SetStats(0);
    data_nominal-> GetXaxis()->LabelsOption("hu");
    data_nominal-> GetXaxis()->SetTitle ("M_{l_{2}l_{3}} (GeV)");	  
    data_nominal->GetXaxis()->SetTitleSize(0.06);
    data_nominal->GetXaxis()->SetTitleOffset(1);
    data_nominal->GetXaxis() ->SetTitleFont(132);
    data_nominal-> GetXaxis()->SetBinLabel(1, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(2, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(3, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(4, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(5, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(6, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(7, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(8, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(9, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(10, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(11, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(12, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(13, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(14, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(15, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(16, "> 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(17, "< 4 ");
    data_nominal-> GetXaxis()->SetBinLabel(18, "> 4 ");   
    data_nominal-> GetXaxis()->SetLabelSize(0.02);
    data_nominal-> GetXaxis()->SetLabelOffset(0.01);
  
  //Make canvas and pads for plotting
  double width, height;  
  width = 800;
  height = 500;

  TCanvas *c =  new TCanvas(name_histo,"",width*(1-xPad),height);   //1000/500
  c->cd();

  TPad *p1; //, *p2;
  //Plot data and MC yields in first pad
  p1 = new TPad(name_histo,"",0,xPad,1,1);
  p1->Draw();
  p1->cd();
  p1->SetTopMargin(0.1);//0.1*(width*(1-xPad)/650)  CHANGE THIS BACK
  p1->SetBottomMargin(0.15);
	
	
	
  for (int i =1; i <= data_nominal->GetNbinsX(); i++){	
	data_down -> SetBinContent(i,data_down -> GetBinContent(i)/data_nominal-> GetBinContent(i));  
	data_up -> SetBinContent(i,data_up -> GetBinContent(i)/data_nominal-> GetBinContent(i));    
	data_nominal-> SetBinContent(i,data_nominal -> GetBinContent(i)/data_nominal-> GetBinContent(i));  	  
  }	  
		
	
  data_nominal -> SetLineColor(kBlack);	
  data_nominal -> SetLineWidth(1);
  data_nominal-> SetMarkerStyle(8);
  data_nominal-> SetMarkerColor(1);
	
  data_down -> SetLineWidth(3);
  data_down -> SetLineColor(kRed);
  data_up -> SetLineWidth(3);
  data_up -> SetLineColor(kBlue);
	
  data_up->SetMinimum(0.6);
  data_up->SetMaximum(1.4);
	
  HistLabelSizes(data_nominal,0.1,0.1,0.07,0.07);
  data_up -> Draw("hist");	
  data_down -> Draw("hist same");	
  data_nominal -> Draw("ep same");	
  legend->Draw("same");		

    
  
  
    Int_t ci;      // for color index setting
   
    TLine *line = new TLine(6.5,0.02,6.5,1.25);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(12.5,0.02,12.5,1.25);
    line->SetLineWidth(2);
    line->Draw();
	
    line = new TLine(2.5,1.1,2.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(1);
    line->Draw();
    line = new TLine(4.5,1.1,4.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(1);
    line->Draw();
    line = new TLine(8.5,1.1,8.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(1);
    line->Draw();
    line = new TLine(10.5,1.1,10.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(1);
    line->Draw();
    line = new TLine(14.5,1.1,14.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(1);
    line->Draw();
    line = new TLine(16.5,1.1,16.5,0.1);

    ci = TColor::GetColor("#ff6600");
    line->SetLineColor(ci);
    line->SetLineWidth(3);
    line->Draw();
    TLatex *    tex = new TLatex(0.8748578,17546.74,"");	
    tex = new TLatex(0.9578318,0.62,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(2.857013,0.62,"#DeltaR [2,10]cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(4.912923,0.62,"#DeltaR > 10cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(0.9578318,0.62,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
		
    tex = new TLatex(6.9578318,0.62,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(8.857013,0.62,"#DeltaR [2,10]cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(10.912923,0.62,"#DeltaR > 10cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();	
		
    tex = new TLatex(12.9578318,0.62,"#DeltaR < 2cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(14.857013,0.62,"#DeltaR [2,10]cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(16.912923,0.62,"#DeltaR > 10cm");

    ci = TColor::GetColor("#ff6600");
    tex->SetTextColor(ci);
    tex->SetTextSize(0.02);
    tex->SetLineWidth(2);
    tex->Draw();	
		
	//25000.83	
    if (channel == 0 ){	
	 //         tex = new TLatex(2.857013,signal[0]->GetBinContent(signal[0]->GetMaximumBin())* 8,"#mu#mu#mu");

      tex = new TLatex(2.857013,0.7,"#mu#mu#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(8.857013,0.7,"#mu^{#pm}#mu^{#mp}e");
      tex->SetTextColor(1);
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(14.857013,0.7,"#mu^{#pm}#mu^{#pm}e");
      tex->SetTextColor(1);
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      tex->Draw();
    }	
    if (channel == 1){	
      tex = new TLatex(2.857013,0.7,"eee");
      tex->SetTextColor(1);
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(8,0.7,"e^{#pm}e^{#mp}#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(13.857013,0.7,"e^{#pm}e^{#pm}#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.04);
      tex->SetLineWidth(2);
      tex->Draw();
    }	
  
	
		
	
  //redraw axis over histograms
  gPad->RedrawAxis();
  //CMS_lumi(c,"Preliminary", true);
  drawLumi(p1);
  if (channel == 1){
    c->SaveAs("plots_pdf/ele_SR/" + name_histo + ".pdf");
    c->SaveAs("plots_root/ele_SR/"+ name_histo + ".root");
  }
	
	
  if (channel == 0){
    c->SaveAs("plots_pdf/mu_SR/"+ name_histo + ".pdf");
    c->SaveAs("plots_root/mu_SR/"+ name_histo + ".root");
  }
}
