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



void plotDataVSMC(int categoria,int channel,int istogramma,
		  TH1D* data, TH1D** bkg,
		  const TString* names, const unsigned nHist,
		  const TString& name_cut,const TString& name_channel, const TString& name_histo,
		  const bool ylog,
		  const unsigned widthopt, const bool plotsig, TH1D** signal , const TString* signames, const unsigned nSig, const bool signorm,
		  const int lumi_year){
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
    if (names[histI[effsam] + 1 + nSig] == "nonprompt DF" ) bkg[effsam]->SetFillColor(kGreen); 
    if (names[histI[effsam] + 1 + nSig] == "nonprompt DF" ) bkg[effsam]->SetLineColor(kGreen);
    if (names[histI[effsam] + 1 + nSig] == "nonprompt SF" ) bkg[effsam]->SetFillColor(kGreen+3); 
    if (names[histI[effsam] + 1 + nSig] == "nonprompt SF" ) bkg[effsam]->SetLineColor(kGreen+3); 

    bkgStack->Add(bkg[effsam], "f");
  }
    
  if(signames==nullptr) {} // dummy, just to avoid warning
    
  //Make a legend for data and all backgrounds
  TLegend* legend = new TLegend(0.16,0.75,0.92,0.87,NULL,"brNDC");
    
  legend->SetFillStyle(0);
  const int signal_out= 14;	 
  unsigned list_signal_out[signal_out] = {1,3,4,5,6,7,8,11,13,14,15,16,17,18};	
	
  //Add signal to the legenD
  if(plotsig){
    for(unsigned sig = 0; sig < nSig; ++sig){
      bool skip_signal = false;   
      for (int j = 0; j < signal_out; j++){
	 if (sig == list_signal_out[j])   skip_signal=true;   	      
      }	  
      if (skip_signal)    continue;
      if ((channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 )  && sig >= 10) continue; 
      if ((channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ) && sig < 10) continue; 
      signal[sig]->SetLineColor(sigCols[sig]);
      signal[sig]->SetMarkerColor(sigCols[sig]);
      if (sig == 0 ||sig == 10 ) signal[sig]->SetLineColor(kBlack);
      if (sig == 0 ||sig == 10 ) signal[sig]->SetMarkerColor(kBlack);
      if (sig == 2 ||sig == 12 ) signal[sig]->SetLineColor(kBlue);
      if (sig == 2 ||sig == 12 ) signal[sig]->SetMarkerColor(kBlue);
      if (sig == 9 ||sig == 19 ) signal[sig]->SetLineColor(kMagenta);
      if (sig == 9 ||sig == 19 ) signal[sig]->SetMarkerColor(kMagenta);

      
      signal[sig]->SetLineWidth(4);
	
      legend->SetTextFont(42);
      //legend->AddEntry(signal[sig], signames[sig]);

	if (sig == 0 ||sig == 10 ) legend->AddEntry(signal[sig], "M = 1GeV, |V|^{2} = 5x10^{-4}");   
	if (sig == 2 ||sig == 12 ) legend->AddEntry(signal[sig], "M = 2GeV, |V|^{2} = 2.4x10^{-4}");  
	if (sig == 19 || sig == 9 ) legend->AddEntry(signal[sig], "M = 10GeV, |V|^{2} = 1x10^{-5}");    
	   		 
    }
  }

    
  for(int effsam = nHist - 1; effsam > -1; --effsam){
    legend->SetTextFont(42);
    //if ( names[histI[effsam] + 1 + nSig] != "TTX"  && names[histI[effsam] + 1 + nSig] != "Xgamma"  && names[histI[effsam] + 1 + nSig] != "nonprompt SF"  && names[histI[effsam] + 1 + nSig] != "nonprompt DF" ) continue;
    //if (names[histI[effsam] + 1 + nSig] == "Xgamma") legend->AddEntry(bkg[effsam], "Conversions");
    //else if (names[histI[effsam] + 1 + nSig] == "TTX") legend->AddEntry(bkg[effsam], "Other");
    //if (effsam == 21 && (names[histI[effsam] + 1 + nSig] != "nonprompt SF" && names[histI[effsam] + 1 + nSig] != "nonprompt DF")) legend->AddEntry(bkg[effsam], "Conversions");  
    if (names[histI[effsam] + 1 + nSig] == "nonprompt SF") legend->AddEntry(bkg[effsam], "nonprompt SF");
    if (names[histI[effsam] + 1 + nSig] == "Xgamma") legend->AddEntry(bkg[effsam], "Conversions");
    if (names[histI[effsam] + 1 + nSig] == "nonprompt DF") legend->AddEntry(bkg[effsam], "nonprompt DF");
    //else if (names[histI[effsam] + 1 + nSig] == "DY") legend->AddEntry(bkg[effsam], "Z#gamma^{*}"); 
    //else 
    //else legend->AddEntry(bkg[effsam], names[histI[effsam] + 1 + nSig]);
    legend->     SetNColumns(3);



    
  }

  if (istogramma == 1 ){
    signal[0]->SetStats(0);
    signal[0]-> GetXaxis()->LabelsOption("vu");
    signal[0]-> GetXaxis()->SetBinLabel(1, "3 leptons");
    signal[0]-> GetXaxis()->SetBinLabel(2, "#DeltaR (l2-l3)");
    signal[0]-> GetXaxis()->SetBinLabel(3, "min  #Delta #phi");
    signal[0]-> GetXaxis()->SetBinLabel(4, "bveto");
    signal[0]-> GetXaxis()->SetBinLabel(5, "resonances");
    signal[0]-> GetXaxis()->SetBinLabel(6, "M_{lll}");
    signal[0]-> GetXaxis()->SetBinLabel(7, "P_{T} (l2-l3)");
    signal[0]-> GetXaxis()->SetBinLabel(8, "cos(SV,HNL)");
    signal[0]-> GetXaxis()->SetBinLabel(9, "SV prob");	
    signal[0]-> GetXaxis()->SetBinLabel(10, "#sigma #Delta (PV-SV)");
    signal[0]-> GetXaxis()->SetBinLabel(11, "M_{ll} (l2-l3)");

    signal[10]->SetStats(0);
    signal[10]-> GetXaxis()->LabelsOption("vu");
    signal[10]-> GetXaxis()->SetBinLabel(1, "3 leptons");
    signal[10]-> GetXaxis()->SetBinLabel(2, "#DeltaR (l2-l3)");
    signal[10]-> GetXaxis()->SetBinLabel(3, "min  #Delta #phi");
    signal[10]-> GetXaxis()->SetBinLabel(4, "bveto");
    signal[10]-> GetXaxis()->SetBinLabel(5, "resonances");
    signal[10]-> GetXaxis()->SetBinLabel(6, "M_{lll}");
    signal[10]-> GetXaxis()->SetBinLabel(7, "P_{T} (l2-l3)");
    signal[10]-> GetXaxis()->SetBinLabel(8, "cos(SV,HNL)");
    signal[10]-> GetXaxis()->SetBinLabel(9, "SV prob");	
    signal[10]-> GetXaxis()->SetBinLabel(10, "#sigma #Delta (PV-SV)");
    signal[10]-> GetXaxis()->SetBinLabel(11, "M_{ll} (l2-l3)");

    signal[0]-> GetXaxis()->SetLabelSize(0.045);
    signal[0]-> GetXaxis()->SetLabelOffset(0.01);
    signal[10]-> GetXaxis()->SetLabelSize(0.045);
    signal[10]-> GetXaxis()->SetLabelOffset(0.01);
  }
  if (istogramma == 2 ){
    signal[0]->SetStats(0);
    signal[0]-> GetXaxis()->LabelsOption("vu");
    signal[0]-> GetXaxis()->SetBinLabel(1, "#DeltaR (l2-l3)");
    signal[0]-> GetXaxis()->SetBinLabel(2, "min  #Delta #phi");
    signal[0]-> GetXaxis()->SetBinLabel(3, "bveto");
    signal[0]-> GetXaxis()->SetBinLabel(4, "resonances");
    signal[0]-> GetXaxis()->SetBinLabel(5, "M_{lll}");
    signal[0]-> GetXaxis()->SetBinLabel(6, "P_{T} (l2-l3)");
    signal[0]-> GetXaxis()->SetBinLabel(7, "cos(SV,HNL)");
    signal[0]-> GetXaxis()->SetBinLabel(8, "SV prob");	
    signal[0]-> GetXaxis()->SetBinLabel(9, "#sigma #Delta (PV-SV)");
    signal[0]-> GetXaxis()->SetBinLabel(10, "M_{ll} (l2-l3)");
    signal[0]-> GetXaxis()->SetBinLabel(11, "full sel");


    signal[10]->SetStats(0);
    signal[10]-> GetXaxis()->LabelsOption("vu");
    signal[10]-> GetXaxis()->SetBinLabel(1, "#DeltaR (l2-l3)");
    signal[10]-> GetXaxis()->SetBinLabel(2, "min  #Delta #phi");
    signal[10]-> GetXaxis()->SetBinLabel(3, "bveto");
    signal[10]-> GetXaxis()->SetBinLabel(4, "resonances");
    signal[10]-> GetXaxis()->SetBinLabel(5, "M_{lll}");
    signal[10]-> GetXaxis()->SetBinLabel(6, "P_{T} (l2-l3)");
    signal[10]-> GetXaxis()->SetBinLabel(7, "cos(SV,HNL)");
    signal[10]-> GetXaxis()->SetBinLabel(8, "SV prob");	
    signal[10]-> GetXaxis()->SetBinLabel(9, "#sigma #Delta (PV-SV)");
    signal[10]-> GetXaxis()->SetBinLabel(10, "M_{ll} (l2-l3)");
    signal[10]-> GetXaxis()->SetBinLabel(11, "full sel");

    signal[0]-> GetXaxis()->SetLabelSize(0.045);
    signal[0]-> GetXaxis()->SetLabelOffset(0.01);
    signal[10]-> GetXaxis()->SetLabelSize(0.045);
    signal[10]-> GetXaxis()->SetLabelOffset(0.01);
  }

  // isotgramma delle SR --> linee e roba varia
	
  TString labels_sr[18]={"0-0.5","0.5-1.5","1.5-4",">4","0-0.5",">0.5","0-0.5","0.5-1.5","1.5-4",">4","0-0.5",">0.5","0-0.5","0.5-1.5","1.5-4",">4","0-0.5",">0.5"};	

  if (istogramma == 0 ){
    signal[0]->SetStats(0);
    signal[0]-> GetXaxis()->LabelsOption("vu");
    signal[0]-> GetXaxis()->SetTitle ("#Delta (PV-SV)_{2D} (cm)");	  
    signal[0]->GetXaxis()->SetTitleSize(0.06);
    signal[0]->GetXaxis()->SetTitleOffset(1.1);
    signal[0]->GetXaxis() ->SetTitleFont(132);
    for (int i =0; i<18; i++){	  
      signal[0]-> GetXaxis()->SetBinLabel(i+1, labels_sr[i]);
    }					     
    signal[10]->SetStats(0);
    signal[10]-> GetXaxis()->LabelsOption("vu");
    signal[10]-> GetXaxis()->SetTitle ("#Delta (PV-SV)_{2D} (cm)");
    signal[10]->GetXaxis()->SetTitleSize(0.06);
    signal[10]->GetXaxis()->SetTitleOffset(1.1);
    signal[10]->GetXaxis() ->SetTitleFont(132);
    for (int i =0; i<18; i++){	  
      signal[10]-> GetXaxis()->SetBinLabel(i+1, labels_sr[i]);
    } 
    signal[0]-> GetXaxis()->SetLabelSize(0.045);
    signal[0]-> GetXaxis()->SetLabelOffset(0.005);
    signal[10]-> GetXaxis()->SetLabelSize(0.045);
    signal[10]-> GetXaxis()->SetLabelOffset(0.005);
	    	    
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

  TPad *p1, *p2;
  //Plot data and MC yields in first pad
  p1 = new TPad(name_histo,"",0,xPad,1,1);
  p1->Draw();
  p1->cd();
  p1->SetTopMargin(0.1);//0.1*(width*(1-xPad)/650)  CHANGE THIS BACK
  //p1->SetBottomMargin(0.15);
  p1->SetBottomMargin(0.0);

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
	if (istogramma == 0) signal[sig] ->GetXaxis()->LabelsOption("v");    
      }
      if ((channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ) && sig == 10){
	if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
	signal[sig]->SetMinimum(0);
	signal[sig]->SetMinimum(0.1);
	signal[sig]->Draw("histe ");
	if (istogramma == 0) signal[sig] ->GetXaxis()->LabelsOption("v");    

      }
    }
  }
  data -> SetMarkerStyle(20);
  data -> SetMarkerColor(1);
  data -> SetLineColor(1); 
  data->Draw("pe same");	//The range used is now that of the data histogra  
  bkgStack->Draw("hist same ");
  data->Draw("pe same");	//The range used is now that of the data histogra  
  legend->AddEntry(data, "data");
  legend->Draw("same");

	
  bkgTot->Draw("e2same");
  if(plotsig){
    for(unsigned sig = 0; sig < nSig; ++sig){
	 bool skip_signal = false;   
         for (int j = 0; j < signal_out; j++){
	    if (sig == list_signal_out[j])   skip_signal=true;   	      
         }	  
      if (skip_signal)    continue;
	    
      if ((channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 )  && sig >= 10) continue; 
      if ((channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ) && sig < 10) continue; 

      if(signorm && signal[sig]->GetSumOfWeights() != 0) signal[sig]->Scale(bkgTot->GetSumOfWeights()/ signal[sig]->GetSumOfWeights());
      signal[sig]->SetMinimum(0.1);
      signal[sig]->Draw("histe same");
    }
  }
	
  if (istogramma == 0){
    //Int_t ci;      // for color index setting
    //TColor *color; // for color definition with alpha
    //    line = new TLine(12.5,0.02,12.5,50000);

    double high_flav=50000;
    double high_mll=6000;
    double left_mll=1.1;
	      
	  
    TLine *line = new TLine(6.5,0.07,6.5, high_flav);
    line->SetLineWidth(2);
    line->Draw();
    line = new TLine(12.5,0.07,12.5, high_flav);
    line->SetLineWidth(2);
    line->Draw();
    //100000
    
    line = new TLine(4.5,high_mll,4.5,0.1);
    //ci = TColor::GetColor("#ff6600");
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw();
    line = new TLine(10.5,high_mll,10.5,0.1);
    //ci = TColor::GetColor("#ff6600");
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw();
    line = new TLine(16.5,high_mll,16.5,0.1);
    //ci = TColor::GetColor("#ff6600");
    line->SetLineStyle(2);
    line->SetLineWidth(1);
    line->Draw();
    
   
    TLatex *    tex = new TLatex(0.8748578,17546.74,"");
    tex = new TLatex(left_mll,high_mll,"M_{ll} < 4 GeV");
    tex->SetTextSize(0.02);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(left_mll+4,high_mll,"M_{ll} > 4 GeV");
    tex->SetTextSize(0.02);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(left_mll+6,high_mll,"M_{ll} < 4 GeV");
    tex->SetTextSize(0.02);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(left_mll+10,high_mll,"M_{ll} > 4 GeV");
    tex->SetTextSize(0.02);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();

    tex = new TLatex(left_mll+12,high_mll,"M_{ll} < 4 GeV");
    tex->SetTextSize(0.02);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
    tex = new TLatex(left_mll+16,high_mll,"M_{ll} > 4 GeV");
    tex->SetTextSize(0.02);
    tex->SetTextFont(42);
    tex->SetLineWidth(2);
    tex->Draw();
   
		
    //25000.83	
    if (channel == 0 ||channel == 1 ||channel == 2 ||channel == 6 ){	
      //         tex = new TLatex(2.857013,signal[0]->GetBinContent(signal[0]->GetMaximumBin())* 8,"#mu#mu#mu");

      tex = new TLatex(2.857013,high_flav/2,"#mu#mu#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.06);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(8.857013,high_flav/2,"#mu^{#pm}#mu^{#mp}e");
      tex->SetTextColor(1);
      tex->SetTextSize(0.06);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(14.857013,high_flav/2,"#mu^{#pm}#mu^{#pm}e");
      tex->SetTextColor(1);
      tex->SetTextSize(0.06);
      tex->SetLineWidth(2);
      tex->Draw();
    }	
    if (channel == 3 ||channel == 4 ||channel == 5 ||channel == 7 ){	
      tex = new TLatex(2.857013,high_flav/2,"eee");
      tex->SetTextColor(1);
      tex->SetTextSize(0.06);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(8,high_flav/2,"e^{#pm}e^{#mp}#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.06);
      tex->SetLineWidth(2);
      tex->Draw();
      tex = new TLatex(14.857013,high_flav/2,"e^{#pm}e^{#pm}#mu");
      tex->SetTextColor(1);
      tex->SetTextSize(0.06);
      tex->SetLineWidth(2);
      tex->Draw();
    }	
  }
	
	
	
	
	
	
  //redraw axis over histograms
  gPad->RedrawAxis();
  //CMS_lumi(c,"Preliminary", true);
  if (channel == 3 || channel == 4 ||channel == 5 ||channel == 7)drawLumi(p1,1, lumi_year);	 
  if (channel == 0 || channel == 1 ||channel == 2 ||channel == 6 )drawLumi(p1,0, lumi_year);




  c->cd();
  
  p2 = new TPad(name_histo + "2","",0,0.0,1,xPad);
  p2->Draw();
  p2->cd();
  p2->SetTopMargin(0);
  p2->SetBottomMargin(0.4);
  // p2->SetLeftMargin(0.1);
	
  TH1D* dataC = (TH1D*) data->Clone();
  TH1D* bkgTotC = (TH1D*) bkgTot->Clone();

  dataC->Divide(bkgTotC);
  dataC->SetMarkerColor(1);
  dataC->SetLineColor(1);
  dataC->GetYaxis()->SetRangeUser(0.,1.999);
  dataC->GetYaxis()->SetTitle("obs/pred");
  dataC->GetYaxis()->SetTitleOffset(0.9/((1.-xPad)/xPad));
  
  dataC->GetYaxis()->SetNdivisions(505);
  
  dataC->GetYaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.06
  dataC->GetXaxis()->SetTitleSize((1.-xPad)/xPad*0.06); //originally 0.09
  
  dataC->GetYaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.05
  dataC->GetXaxis()->SetLabelSize((1.-xPad)/xPad*0.05); //originally 0.075

  if (istogramma == 0 && categoria == 6 && channel == 6){
    std::cout<<"----------------     muon coupling: "<<std::endl;
    for (int i =0; i <dataC->GetNbinsX(); i++ ){
      std::cout<<i<<") "<< dataC->GetBinContent(i+1) << "  ±  "<<dataC->GetBinError(i+1) <<std::endl;
    }
  }
  if (istogramma == 0 && categoria == 6 && channel == 7){
    std::cout<<"----------------    electron coupling: "<<std::endl;
    for (int i =0; i <dataC->GetNbinsX(); i++ ){
      std::cout<<i<<") "<< dataC->GetBinContent(i+1)<<  "  ±  "<<dataC->GetBinError(i+1) <<std::endl;
    }
  }
  dataC->Draw("pe");
  //Draw line at 1 on ratio plot
  double xmax = dataC->GetBinCenter(dataC->GetNbinsX()) + dataC->GetBinWidth(dataC->GetNbinsX())/2;
  double xmin = dataC->GetBinCenter(0) + dataC->GetBinWidth(0)/2;
  TLine *line = new TLine(xmin, 1, xmax, 1);
  line->SetLineStyle(2);
  line->Draw("same");



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


