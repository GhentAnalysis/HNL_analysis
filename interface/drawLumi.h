#ifndef DRAW_LUMI
#define DRAW_LUMI
#include "TPad.h"
#include "TLatex.h"
#include "TLine.h"

void drawLumi(TPad*,int channel = 0, const TString& extraText = "Preliminary", const bool data = true, const int lumi_year =0 );


void drawLumi(TPad* pad,int channel, const TString& extraText, const bool data, const int lumi_year){
	
	TString lumiText;
	TString channelText;
	lumiText = "35.9 fb^{-1} (13 TeV)";
	//if(data && lumi_year==1) lumiText = "41.5 fb^{-1} (13 TeV)";
	//if(data && lumi_year==2) lumiText = "59.7 fb^{-1} (13 TeV)";
	//else lumiText = "(13 TeV)";

	//const float H = pad->GetWh();
  	//const float W = pad->GetWw();
  	const float l = pad->GetLeftMargin();
  	const float t = pad->GetTopMargin();
  	const float r = pad->GetRightMargin();
  	//const float b = pad->GetBottomMargin();

	float CMSTextSize = pad->GetTopMargin()*0.6;
	float lumiTextSize = pad->GetTopMargin()*0.4;
	float channelTextSize = pad->GetTopMargin()*0.6;

	//float CMSTextOffset = pad->GetTopMargin()*0.7;
	float lumiTextOffset = pad->GetTopMargin()*0.7;
	
	pad->cd();
	//Define latex text to draw on plot
	TLatex latex(l,1-t+lumiTextOffset*t,"CMS");
	latex.SetNDC();
	latex.SetTextAngle(0);
	latex.SetTextColor(kBlack); 

	latex.SetTextFont(61);
	latex.SetTextAlign(11); 
	latex.SetTextSize(CMSTextSize);
	const float cmsX = latex.GetXsize();
	latex.DrawLatex(l,1-t+lumiTextOffset*t,"CMS");

	const float extraTextSize = CMSTextSize*0.76;	 
	latex.SetTextFont(52);
	latex.SetTextSize(extraTextSize);
	latex.SetTextAlign(11);
	latex.DrawLatex(l + 1.2*cmsX, 1-t+lumiTextOffset*t, extraText); //relPosX*(1-l-r)

	latex.SetTextFont(42);
	latex.SetTextAlign(31); //31
	latex.SetTextSize(lumiTextSize);  
	latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);
	
	latex.SetTextFont(42);
	latex.SetTextAlign(11); 
	latex.SetTextSize(channelTextSize);	
	if (channel == 0) latex.DrawLatex(l + 4.5*cmsX,1-t+lumiTextOffset*t,"#mu#mu + l");
	if (channel == 1) latex.DrawLatex(l + 4.5*cmsX,1-t+lumiTextOffset*t,"ee + l");

	return;
}
#endif

