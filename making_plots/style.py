import ROOT
from ROOT import gROOT, TStyle

def setStyle():

    style = TStyle( 'style', 'style' ) 
    
    style.SetCanvasBorderMode(0)
    style.SetCanvasColor(ROOT.kWhite)
    style.SetCanvasDefH(600)
    style.SetCanvasDefW(600)
    style.SetCanvasDefX(0)
    style.SetCanvasDefY(0)
    
    style.SetPadBorderMode(0)
    style.SetPadColor(ROOT.kWhite)
    style.SetPadGridX(False)
    style.SetPadGridY(False)
    style.SetGridColor(0)
    style.SetGridStyle(3)
    style.SetGridWidth(1)
    
    style.SetFrameBorderMode(0)
    style.SetFrameBorderSize(1)
    style.SetFrameFillColor(0)
    style.SetFrameFillStyle(0)
    style.SetFrameLineColor(1)
    style.SetFrameLineStyle(1)
    style.SetFrameLineWidth(1)
    
    style.SetHistLineColor(1)
    style.SetHistLineStyle(0)
    style.SetHistLineWidth(1)
    
    style.SetEndErrorSize(2)
    style.SetMarkerStyle(20)
    
    #//For the fit/function:
    #gStyle->SetOptFit(1);
    #gStyle->SetFitFormat("5.4g");
    #gStyle->SetFuncColor(2);
    #gStyle->SetFuncStyle(1);
    #gStyle->SetFuncWidth(1);
    
    #//For the date:
    #gStyle->SetOptDate(0);
    #// gStyle->SetDateX(Float_t x = 0.01);
    #// gStyle->SetDateY(Float_t y = 0.01);
    
    #// For the statistics box:
    style.SetOptFile(0)
    style.SetOptStat(0)
    style.SetStatColor(ROOT.kWhite)
    style.SetStatFont(42)
    style.SetLegendFont(42)
    style.SetStatFontSize(0.08)
    style.SetStatTextColor(1)
    style.SetStatFormat("6.4g")
    style.SetStatBorderSize(1)
    style.SetStatH(0.7)
    style.SetStatW(0.15)
    #//gStyle->SetStatTextSize(2.5);
    
    #gStyle->SetStatX(0.96);
    #//gStyle->SetStatY(0.35);
    #// gStyle->SetStatStyle(Style_t style = 1001);
    #// gStyle->SetStatX(Float_t x = 0);
    #// gStyle->SetStatY(Float_t y = 0);
    
    #Margins:
    style.SetPadTopMargin(0.05)
    style.SetPadBottomMargin(0.13)
    style.SetPadLeftMargin(0.16)
    style.SetPadRightMargin(0.04)
    
    # For the Global title:
    style.SetOptTitle(0)
    style.SetTitleFont(42)
    style.SetTitleColor(1)
    style.SetTitleTextColor(1)
    style.SetTitleFillColor(10)
    style.SetTitleFontSize(0.05)
    
    
    #For the axis titles:
    style.SetTitleColor(1, "XYZ")
    style.SetTitleFont(42, "XYZ")
    style.SetTitleSize(0.06, "XYZ")
    style.SetTitleXOffset(0.9)
    style.SetTitleYOffset(1.25)
    
    #For the axis labels:
    style.SetLabelColor(1, "XYZ")
    style.SetLabelFont(42, "XYZ")
    style.SetLabelOffset(0.007, "XYZ")
    style.SetLabelSize(0.05, "XYZ")
    
    #For the axis:
    style.SetAxisColor(1, "XYZ")
    style.SetStripDecimals(ROOT.kTRUE)
    style.SetTickLength(0.03, "XYZ")
    style.SetNdivisions(505, "XYZ")
    style.SetPadTickX(1)
    style.SetPadTickY(1)
    
    ##Change for log plots:
    #gStyle->SetOptLogx(0);
    #gStyle->SetOptLogy(0);
    #gStyle->SetOptLogz(0);
    
    #// Postscript options:
    #gStyle->SetPaperSize(20.,20.);
    #// gStyle->SetLineScalePS(Float_t scale = 3);
    #// gStyle->SetLineStyleString(Int_t i, const char* text);
    #// gStyle->SetHeaderPS(const char* header);
    #// gStyle->SetTitlePS(const char* pstitle);
    
    #// gStyle->SetBarOffset(Float_t baroff = 0.5);
    #// gStyle->SetBarWidth(Float_t barwidth = 0.5);
    #// gStyle->SetPaintTextFormat(const char* format = "g");
    #// gStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
    #// gStyle->SetTimeOffset(Double_t toffset);
    #// gStyle->SetHistMinimumZero(kTRUE);
    
    #//For error in SF maps
    #gStyle->SetMarkerSize(0.9); //TEMPORARY SIZE FOR DILEPTON PLOTS. SET BACK TO DEFAULT FOR TRILEPTON
    #gStyle->SetPaintTextFormat("4.2f");  //4.2
    #//gStyle->SetHatchesLineWidth(5);
    #//gStyle->SetHatchesSpacing(0.05);
    
    style.SetLegendBorderSize(0)
    gROOT.SetStyle('style')
    gROOT.ForceStyle()

    return style
