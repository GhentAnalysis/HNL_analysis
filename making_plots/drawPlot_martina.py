
import os
import math
import sys


import ROOT
from ROOT import TCanvas, TPad, TGraphAsymmErrors, TLegend, TLine, TLatex


from drawCMSHeader import drawCMSHeader


def colorHistogram( hist, color ):
	
    #black lines between processes
    #hist.SetLineColor( ROOT.kBlack )
    hist.SetLineColor( color )
    hist.SetLineWidth( 1 )
    
    #set the requested histogram color
    hist.SetFillColor( color )
    hist.SetMarkerColor( color )


def colorProcesses( processCollection, color_dict ):
    for p in processCollection:
        if p.isSignal() : continue
        colorHistogram( p.nominal(), color_dict[ p.name() ] )



def prepareData( data ):

    #set poison errors
    data.SetBinErrorOption( ROOT.TH1.kPoisson )

    #build TGraphAsymmErrors for plotting
    data_graph = TGraphAsymmErrors( data )
    for b in range( 1, data.GetNbinsX() + 1 ):
        bin_content = data.GetBinContent( b )
        data_graph.SetPointError( b - 1, 0, 0,  data.GetBinErrorLow( b ) if bin_content >= 0. else 0., data.GetBinErrorUp( b ) if bin_content >= 0. else 0. )

    #hack for not displaying points at zero 
    maximum = ( data.GetBinContent( data.GetMaximumBin() ) + data.GetBinErrorUp( data.GetMaximumBin() ) )
    for b in range( 1, data.GetNbinsX() + 1 ):
        if data.GetBinContent( b ) < 1e-8:
            data_graph.GetY()[b - 1] += ( maximum * 1e8 )

    return data, data_graph
    

def makeAndDivideCanvas( width, height, lower_pad_fraction ):
    c = TCanvas( "", "", width, height )

    upper_pad = TPad( "", "", 0, lower_pad_fraction, 1, 1 )
    upper_pad.SetBottomMargin( 0.03 )
    upper_pad.SetTopMargin( 0.08 )
    upper_pad.SetLeftMargin( 0.1 )



    lower_pad = TPad( "", "", 0, 0, 1, lower_pad_fraction )
    lower_pad.SetTopMargin( 0.01 )
    lower_pad.SetBottomMargin( 0.4 )
    lower_pad.SetLeftMargin( 0.1 )

    return c, upper_pad, lower_pad


def makeUpperLegend( data, processCollection, bkg_total, plot_axis= 'SR', legend_names = None ):

    if plot_axis == 'SR': legend = TLegend(0.12, 0.73, 0.92, 0.90, '', 'brNDC');
    else : legend = TLegend(0.12, 0.70, 0.92, 0.90, '', 'brNDC');
        
    legend.SetNColumns( 3 )
    legend.SetFillStyle( 0 ) #avoid box
    legend.SetBorderSize(0)
    
    legend.AddEntry( data, 'Data', 'pe1' )
    for p in processCollection:
        if p.isSignal(): continue
        name = p.name()
        if legend_names is not None:
            try:
                name = legend_names[ p.name () ]
            except KeyError:
                pass
        legend.AddEntry( p.nominal(), name, 'f' )
    legend.AddEntry( bkg_total, 'Total bkg. unc.', 'f' )
    legend.AddEntry( None, "", "");
    return legend


#compute maximum entry to be drawn on plot
def maximum( data, bkg_total ):
	data_max = ( data.GetBinContent( data.GetMaximumBin() ) + data.GetBinErrorUp( data.GetMaximumBin() ) )
	bkg_max = ( bkg_total.GetBinContent( bkg_total.GetMaximumBin() ) + bkg_total.GetBinErrorUp( bkg_total.GetMaximumBin() ) )
	return max( data_max, bkg_max )


#compute minumum entry to be drawn on plot, but ignore zeroes
def minimum( data, bkg_total ):
	total_min = maximum( data, bkg_total )

	assert ( data.GetNbinsX() == bkg_total.GetNbinsX() )
	for b in range( 1, data.GetNbinsX() + 1 ):
		data_bin = data.GetBinContent( b )
		bkg_bin = bkg_total.GetBinContent( b )
		if data_bin > 0 and data_bin < total_min :
			total_min = data_bin
		if bkg_bin > 0 and bkg_bin < total_min :
			total_min = bkg_bin
	return total_min


#compute range of upper plot 
def rangeLinear( data, background, signal):
	total_max = maximum( signal, background )
	total_max = (total_max * 1.3)
	return 0, total_max


def rangeLog( data, background,signal, flav_name, plot_axis, lumi_text  ):
	
    total_min = minimum( data, background )
    #total_max = maximum( data, background )
    total_max = maximum( signal, background )

    #set minimum to be 5 times smaller than the smallest background yield 
    #total_min /= 2
    #total_min = 0.02
    if total_min < 0.08:
       total_min = 0.08
    #total_min = 0.0001    
    print( total_min )
    
    #compute the number of axis divisions ( powers of 10 ) between minimum and maximum
    number_of_orders = math.log10( total_max / total_min )
    
    #the plot maximum should be 50% higher in terms of relative canvas size than total_max 
    if flav_name == 'muon' and lumi_text != '35.9 fb^{-1} (13 TeV)' :
        total_max = total_max * 10**( 0.8 * number_of_orders )
    if flav_name == 'muon' and lumi_text == '35.9 fb^{-1} (13 TeV)' :
        total_max = total_max * 10**( 1.0 * number_of_orders )    
    if flav_name == 'ele' and lumi_text == '41.5 fb^{-1} (13 TeV)':
        total_max = total_max * 10**( 1.3 * number_of_orders )
    if flav_name == 'ele' and lumi_text != '41.5 fb^{-1} (13 TeV)':
        total_max = total_max * 10**( 0.8 * number_of_orders )    
    return total_min, total_max


def setRelativeUncStyle( relative_unc, color ):
	relative_unc.SetFillStyle( 1001 )
	relative_unc.SetMarkerStyle( 1 )
	relative_unc.SetFillColor( color )
	relative_unc.SetLineColor( color )


def setRatioPlotStyle( first_hist, lower_pad_fraction,plot_axis ):
    scale_factor = ( 1. - lower_pad_fraction ) / lower_pad_fraction 
    first_hist.GetYaxis().SetTitleOffset( 0.8 / scale_factor )
    first_hist.GetXaxis().SetTitleOffset( 1 )
    first_hist.GetYaxis().SetTitleSize( scale_factor * .06 )
    first_hist.GetXaxis().SetTitleSize( scale_factor * .06 )
    first_hist.GetYaxis().SetLabelSize( scale_factor * .05 )
    first_hist.GetYaxis().SetLabelOffset( 0.013 )    
    first_hist.GetXaxis().SetLabelSize( scale_factor * .06 )
    first_hist.GetXaxis().LabelsOption("vu");
    first_hist.GetXaxis().SetLabelOffset( scale_factor * 0.01 )
    if plot_axis == 'massl2l3':
        first_hist.GetXaxis().SetNdivisions(206)
        first_hist.GetXaxis().SetTickLength(0.08)
    if plot_axis == 'displacement':
        first_hist.GetXaxis().SetNdivisions(506)
        first_hist.GetXaxis().SetTickLength(0.08)    
    labels_sr=["0-0.5","0.5-1.5","1.5-4",">4","0-0.5",">0.5","0-0.5","0.5-1.5","1.5-4",">4","0-0.5",">0.5","0-0.5","0.5-1.5","1.5-4",">4","0-0.5",">0.5"]
    if plot_axis == 'SR':
        for b in range( 1, first_hist.GetNbinsX() + 1 ):
            name = labels_sr[b-1]
            first_hist.GetXaxis().SetBinLabel(b, name)
            first_hist.GetXaxis().LabelsOption("vu");

#make uncertainty band for ratio plot
def relativeUncBand( bkg_hist ):

    for b in range( 1, bkg_hist.GetNbinsX() + 1 ):
    
        #transform total uncertainties to relative uncertainties
        bin_content = bkg_hist.GetBinContent( b )
        if bin_content > 0:
            bkg_hist.SetBinError( b, bkg_hist.GetBinError( b ) / bin_content )
        else:
            bkg_hist.SetBinError( b, 0. )
        
        #center band at 1
        bkg_hist.SetBinContent( b, 1. )
    
    return bkg_hist


#Make TGraphAsymmErrors representing observed/predicted with statistical uncertainties from data 
def obsOverPredWithDataStat( data_hist, bkg_total ):
    obs_over_pred = TGraphAsymmErrors( data_hist )
    for b in range( 1, data_hist.GetNbinsX() + 1 ):
    
        #divide data by background
        bkg_bin = bkg_total.GetBinContent( b )
        if abs( bkg_bin ) > 1e-8:
            obs_over_pred.GetY()[ b - 1] /= bkg_bin
            obs_over_pred.SetPointError( b - 1, 0., 0., data_hist.GetBinErrorLow( b ) / bkg_bin, data_hist.GetBinErrorUp( b ) / bkg_bin )
        else:
            obs_over_pred.GetY()[ b - 1 ] = 0
            obs_over_pred.SetPointError( b - 1, 0., 0., 0, 0 )
    
    return obs_over_pred
 
	
#draw horizontal line with the same range as the plot
def horizontalLine( data, center ):
    line_begin = data.GetBinLowEdge( 1 )
    line_end = data.GetBinLowEdge( data.GetNbinsX() ) + data.GetBinWidth( data.GetNbinsX() )
    line = TLine( line_begin, center, line_end, center )
    line.SetLineStyle(2);
    return line

#draw vertical line with the same range as the plot
def verticalLine_1( data, center,line_begin,line_end  ):
    line = TLine(center, line_begin, center, line_end)
    line.SetLineWidth(2);
    line.SetLineColor(ROOT.kCyan+3);
    return line
def verticalLine_2( data, center,line_begin,line_end  ):
    line = TLine(center, line_begin, center, line_end)
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    return line


#make ratio plot legend
def makeLowerLegend( obs_over_pred, relative_bkg_stat_unc, relative_bkg_total_unc ):
    legend = TLegend( 0.18, 0.83, 0.94, 0.95, '', 'brNDC' )
    legend.SetNColumns( 3 )
    legend.SetFillStyle( 0 )
    legend.SetBorderSize( 0 )
    legend.AddEntry( relative_bkg_stat_unc, 'Stat. pred. unc.', 'f' )
    legend.AddEntry( relative_bkg_total_unc, 'Total pred. unc.', 'f' )
    legend.AddEntry( obs_over_pred, 'Obs./Bgk.', 'pe12' )
    return legend



def makeLabel1( label_text, xposition, yposition, taglias ):	
    label = TLatex( xposition, yposition, label_text )
    #label.SetNDC()
    label.SetTextAlign( 21 )
    label.SetTextFont( 42 )
    label.SetTextColor( ROOT.kCyan+3 )
    label.SetTextSize(  taglias )
    label.SetTextAngle( 0 )
    return label

def makeLabel( label_text, xposition, yposition, taglias ):
    label = TLatex( xposition, yposition, label_text )
    #label.SetNDC()
    label.SetTextAlign( 21 )
    label.SetTextFont( 42 )
    label.SetTextColor( ROOT.kBlack )
    label.SetTextSize( taglias )
    label.SetTextAngle( 0 )
    return label




    
        
def drawPlot( data, processCollection, processCollection2, processCollection3,processCollection4, plot_name, color_dict = None, log = True, lower_pad_fraction = 0.25, legend_names = None, lumi_text = '137 fb^{-1} (13 TeV)',plot_axis = 'SR',flav_name = 'muon',mass1_name = '1gev',mass2_name = '2gev',mass3_name = '3gev',v1_name = 'v1',v2_name = 'v2',v3_name = 'v3',width = 800, height = 500, additional_label = 'pippo', additional_label2 = 'pippo'):
    
    if plot_axis == 'SR':
        width = 1000
        height = 600
        log = True
    if plot_axis == 'massl2l3' or plot_axis == 'displacement' or plot_axis == 'mass3':
        width = 700
        height = 600
        log = True

  
    v1_name = v1_name.replace('_', '.' )
    v1_name = v1_name.replace('=', '#times10^{#minus' )
    v2_name = v2_name.replace('_', '.' )
    v2_name = v2_name.replace('=', '#times10^{#minus' )
    v3_name = v3_name.replace('_', '.' )
    v3_name = v3_name.replace('=', '#times10^{#minus' )
    legend_label_signal1 = 'M = ' + mass1_name + ", |V|^{2} = " + v1_name + '}'
    legend_label_signal2 = 'M = ' + mass2_name + ", |V|^{2} = " + v2_name + '}'
    legend_label_signal3 = 'M = ' + mass3_name + ", |V|^{2} = " + v3_name + '}'
    print   legend_label_signal1   
    print   legend_label_signal2  
    print   legend_label_signal3   
    

        
      
    #order background processes by yield
    processCollection.orderByYield()
    
    #color the histograms
    if color_dict is not None:
        colorProcesses( processCollection, color_dict )
    
    
    #Replace data by TGRaphAsymmErrors for plotting
    data, data_graph = prepareData( data )
    
    
    stack = processCollection.backgroundStack()
    bkg_total_syst = processCollection.totalBkgWithAllErrors()
    signal_histo =  processCollection.signal_get()
    signal_histo2 =  processCollection2.signal_get()
    signal_histo3 =  processCollection3.signal_get()
    signal_histo4 =  processCollection4.signal_get()

    #nasty hack for martina
    if log:
        print ("we enter in the log")
        for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
            print b
            print bkg_total_syst.GetBinContent( b )
            if bkg_total_syst.GetBinContent( b ) < 0.08:
                print ("the content is below 1e-5")
                print bkg_total_syst.GetBinContent( b )
                error = bkg_total_syst.GetBinError( b )
                print ("error of the bin")
                print bkg_total_syst.GetBinError( b )
                bkg_total_syst.SetBinError( b, error / 2 )
                bkg_total_syst.SetBinContent( b, error / 2 )

    bkg_total_syst.SetFillStyle( 3005 )#3244    
    bkg_total_syst.SetFillColor( ROOT.kGray + 2 )
    bkg_total_syst.SetMarkerStyle( 0 )
    signal_histo.SetMarkerStyle( 1 )
    signal_histo.SetLineColor(ROOT.kBlack )
    signal_histo.SetLineWidth( 3 )
    signal_histo.SetMarkerColor( ROOT.kBlack )

    signal_histo2.SetMarkerStyle( 1 )
    signal_histo2.SetLineColor(ROOT.kBlue )
    signal_histo2.SetLineWidth( 3 )
    signal_histo2.SetMarkerColor( ROOT.kBlue  )

    signal_histo3.SetMarkerStyle( 1 )
    signal_histo3.SetLineColor(ROOT.kRed   )
    signal_histo3.SetLineWidth( 3 )
    signal_histo3.SetMarkerColor( ROOT.kRed  )

    signal_histo4.SetMarkerStyle( 1 )
    signal_histo4.SetLineColor(ROOT.kMagenta )
    signal_histo4.SetLineWidth( 2 )
    signal_histo4.SetMarkerColor( ROOT.kMagenta )

    


    c, upper_pad, lower_pad = makeAndDivideCanvas( width, height, lower_pad_fraction )
    
    
    #draw data and backgrounds on upper pad
    upper_pad.Draw()
    upper_pad.cd()
    if log:
        upper_pad.SetLogy()
    
    legend = makeUpperLegend( data_graph, processCollection, bkg_total_syst, plot_axis, legend_names)
    
    #determine plot range 
    if log:
        range_min, range_max = rangeLog( data, bkg_total_syst,signal_histo, flav_name, plot_axis, lumi_text )
    else:
        range_min, range_max = rangeLinear( data, bkg_total_syst, signal_histo )
    range_min = 0.08    
    bkg_total_syst.SetMinimum( range_min )
    bkg_total_syst.SetMaximum( range_max )

    #massimo = maximum( signal_histo, bkg_total_syst )
    massimo = maximum( signal_histo, data )

    
    #only draw labels in bottom pad
    bkg_total_syst.GetXaxis().SetLabelSize(0);
    bkg_total_syst.GetYaxis().SetTitle( 'Events' )
    bkg_total_syst.GetYaxis().SetTitleOffset( 0.8 )
    
    bkg_total_syst.Draw( 'e2' )
    stack.Draw( 'histsame' )
    #legend.AddEntry (signal_histo,  'M = 2GeV, |V|^{2} = 1x10^{-4}')
    #legend.AddEntry (signal_histo2,  'M = 4GeV, |V|^{2} = 8x10^{-6}')
    #legend.AddEntry (signal_histo3,  'M = 8GeV, |V|^{2} = 2x10^{-6}')
    legend.AddEntry (signal_histo,   legend_label_signal1  )
    legend.AddEntry (signal_histo2,  legend_label_signal2  )
    legend.AddEntry (signal_histo3,  legend_label_signal3  )
    legend.Draw( 'same' )

   
    
    max_flav_line = massimo*4.5
    max_flav_name = massimo*4
    max_mass_line = massimo*2
    max_mass_name = massimo*1.7
        
    if flav_name == 'ele':
        max_flav_line = massimo*3
        max_flav_name = massimo*2.
        max_mass_line = massimo*1.5
        max_mass_name = massimo*1.3

     
        

    #lines for cosmesi
    line1 = verticalLine_1( data, 6.5, range_min, max_flav_line)
    line2 = verticalLine_1( data, 12.5, range_min, max_flav_line)
    line3 = verticalLine_2( data, 4.5, range_min, max_mass_line)
    line4 = verticalLine_2( data, 10.5, range_min, max_mass_line)
    line5 = verticalLine_2( data, 16.5, range_min, max_mass_line)
    if plot_axis == 'SR':
        line1.Draw('same')    
        line2.Draw('same')
        line3.Draw('same')    
        line4.Draw('same')
        line5.Draw('same')

    #labales channels
    if flav_name == 'muon':
        label1 = makeLabel1( "#mu#mu#mu", 3, max_flav_name, 0.08 )
        label2 = makeLabel1( "#mu^{#pm}#mu^{#mp}e",9, max_flav_name, 0.08 )
        label3 = makeLabel1( "#mu^{#pm}#mu^{#pm}e", 15, max_flav_name, 0.08 )
    if flav_name == 'ele':     
        label1 = makeLabel1( "eee", 3,  max_flav_name, 0.08 )
        label2 = makeLabel1( "e^{#pm}e^{#mp}#mu", 9,  max_flav_name, 0.08 )
        label3 = makeLabel1( "e^{#pm}e^{#pm}#mu", 15,  max_flav_name, 0.08 )
    if plot_axis == 'SR':
        label1.Draw( 'same' )
        label2.Draw( 'same' )
        label3.Draw( 'same' )

    label11 = makeLabel( "M_{ll} < 4 GeV", 3.4, max_mass_name, 0.04 )
    label22 = makeLabel( "M_{ll} > 4 GeV", 5.5, max_mass_name, 0.04 )
    label33 = makeLabel( "M_{ll} < 4 GeV", 9.4, max_mass_name, 0.04 )
    label44 = makeLabel( "M_{ll} > 4 GeV", 11.5, max_mass_name, 0.04 )
    label55 = makeLabel( "M_{ll} < 4 GeV", 15.4, max_mass_name, 0.04 )
    label66 = makeLabel( "M_{ll} > 4 GeV", 17.5, max_mass_name, 0.04 )
    if plot_axis == 'SR':
        label11.Draw( 'same' )
        label22.Draw( 'same' )
        label33.Draw( 'same' )
        label44.Draw( 'same' )
        label55.Draw( 'same' )
        label66.Draw( 'same' )
    

    #redraw total background uncertainty so it overlays the stack
    bkg_total_syst.SetLineColor( ROOT.kGray )
    bkg_total_syst.Draw( 'e2same' )
    data_graph.Draw( 'pe1same' )
    signal_histo.Draw('histe same')
    signal_histo2.Draw('histe same')
    signal_histo3.Draw('histe same')

    upper_pad.RedrawAxis()
    if plot_axis == 'SR':
        if flav_name == 'muon': drawCMSHeader( upper_pad, lumi_text, 'Preliminary', '' )
        if flav_name == 'ele': drawCMSHeader( upper_pad, lumi_text, 'Preliminary', '' )
    else:
        if flav_name == 'muon': drawCMSHeader( upper_pad, lumi_text, 'Preliminary', '#mu#mu + l' )
        if flav_name == 'ele': drawCMSHeader( upper_pad, lumi_text, 'Preliminary', 'ee + l' )

    upper_pad.RedrawAxis()
    
   #  if additional_label is not None:
   #  	label = makeLabel( additional_label )
    # 	label.Draw( 'same' )
    
    c.cd()
    
    #make ratio plot in lower pad 
    lower_pad.Draw()
    lower_pad.cd()
    
    #make separate histograms containing total and statistical background uncertainty which will be used to plot uncertainty band
    relative_bkg_stat_unc = relativeUncBand( processCollection.totalBkgWithStatErrors() )
    setRelativeUncStyle( relative_bkg_stat_unc, ROOT.kCyan - 4 )
    relative_bkg_total_unc = relativeUncBand( processCollection.totalBkgWithAllErrors() )
    setRelativeUncStyle( relative_bkg_total_unc, ROOT.kOrange - 4 )
    
    
    #make TGraphAsymmErrors representing observed/predicted yields with statistical uncertainties from data
    obs_over_pred = obsOverPredWithDataStat( data, bkg_total_syst )
    
    #set name and range of ratio plot 
    relative_bkg_total_unc.GetYaxis().SetRangeUser( 0.,2.45 )
    relative_bkg_total_unc.GetYaxis().SetTitle( 'Obs./Bgk.' )
    relative_bkg_total_unc.GetYaxis().SetLabelSize( 0.7 )    
    
    if plot_axis == 'SR':
        relative_bkg_total_unc.GetXaxis().SetTitle( '#Delta (PV-SV)_{2D} (cm)' )
    elif plot_axis == 'massl2l3':
        relative_bkg_total_unc.GetXaxis().SetTitle( 'M_{ll}#left(l_{2}+l_{3} #right) (GeV)' )
    elif plot_axis == 'displacement':
        relative_bkg_total_unc.GetXaxis().SetTitle( '#Delta (PV-SV)_{2D} (cm)' )
    elif plot_axis == 'mass3':
        relative_bkg_total_unc.GetXaxis().SetTitle( 'M_{lll} (GeV)' )      
    #set label sizes for ratio plot 
    setRatioPlotStyle( relative_bkg_total_unc, lower_pad_fraction, plot_axis )
    
    #make legend for ratio plot
    lower_legend = makeLowerLegend( obs_over_pred, relative_bkg_stat_unc, relative_bkg_total_unc )
    relative_bkg_total_unc.Draw( 'e2' )
    relative_bkg_stat_unc.Draw( 'e2same' )
    obs_over_pred.Draw( 'pe01same' )
    lower_legend.Draw( 'same' )
    if plot_axis == 'SR': relative_bkg_total_unc.GetXaxis().LabelsOption("vu");  
    #relative_bkg_total_unc.GetXaxis().LabelsOption("v");

    #draw horizontal line at 1 on ratio plot
    line = horizontalLine( data, 1. )
    line.Draw('same')


    line11 = verticalLine_1( relative_bkg_total_unc, 6.5, 0.,2.3 )
    line22 = verticalLine_1( relative_bkg_total_unc, 12.5,0.,2.3 )
    line33 = verticalLine_2( relative_bkg_total_unc, 4.5, 0.,2.3 )
    line44 = verticalLine_2( relative_bkg_total_unc, 10.5,0.,2.3 )
    line55 = verticalLine_2( relative_bkg_total_unc, 16.5, 0.,2.3 )
    if plot_axis == 'SR':
        line11.Draw('same')    
        line22.Draw('same')
        line33.Draw('same')    
        line44.Draw('same')
        line55.Draw('same')
    
    lower_pad.RedrawAxis();
    
    #remove possible file extension from plot name and print it as pdf and png
    print ("qui dovrebbe printate")
    plot_name = os.path.splitext( plot_name )[0]
    c.SaveAs( 'pPdf/'+ plot_name + '.pdf' )
    c.SaveAs( 'pPng/'+ plot_name + '.png' )
    c.SaveAs('pRoot/' + plot_name + '.root' )










    
def writeTables_sumErrors(data, processCollection, processCollection2, processCollection3,processCollection4, plot_name,plot_axis = 'SR',flav_name = 'muon', year_label = '16',mass1_name = '1gev',mass2_name = '2gev',mass3_name = '3gev',v1_name = 'v1',v2_name = 'v2',v3_name = 'v3'):



    
	#Replace data by TGRaphAsymmErrors for plotting
    data, data_graph = prepareData( data )
 
    stack = processCollection.backgroundStack()
    bkg_total_syst = processCollection.totalBkgWithAllErrors()
    bkg_total_onlysyst = processCollection.totalBkgWithSystErrors()
    bkg_total_onlystat = processCollection.totalBkgWithStatErrors()

    bkg_np_syst =     processCollection.npBkgWithAllErrors()
    bkg_np_onlysyst = processCollection.npBkgWithSystErrors()
    bkg_np_onlystat = processCollection.npBkgWithStatErrors()

    bkg_p_syst =     processCollection.pBkgWithAllErrors()
    bkg_p_onlysyst = processCollection.pBkgWithSystErrors()
    bkg_p_onlystat = processCollection.pBkgWithStatErrors()


    signal_histo1 =  processCollection.signal_get()
    signal_histo1_onlysyst = processCollection.signalWithSystErrors()
    signal_histo1_onlystat = processCollection.signalWithStatErrors()

    signal_histo2 =  processCollection2.signal_get()
    signal_histo2_onlysyst = processCollection2.signalWithSystErrors()
    signal_histo2_onlystat = processCollection2.signalWithStatErrors()
    
    signal_histo3 =  processCollection3.signal_get()
    signal_histo3_onlysyst = processCollection3.signalWithSystErrors()
    signal_histo3_onlystat = processCollection3.signalWithStatErrors()
    

    name_txt = 'table_bgk'
    if year_label == 'combined':
        if flav_name == 'ele': name_txt = 'bgk_sum_error_ele_combined.txt'
        if flav_name == 'muon': name_txt = 'bgk_sum_error_muon_combined.txt'
    else:
        if flav_name == 'ele': name_txt = 'bgk_sum_error_ele'  +'_' + year_label  +'.txt'
        if flav_name == 'muon': name_txt = 'bgk_sum_error_muon' + '_' + year_label  +'.txt'

    name_txt_signal = 'table_signal'
    if year_label == 'combined':
        if flav_name == 'ele': name_txt_signal = mass1_name+'_signal_sum_error_ele_combined.txt'
        if flav_name == 'muon': name_txt_signal = mass1_name+'_signal_sum_error_muon_combined.txt'
    else:
        if flav_name == 'ele': name_txt_signal = mass1_name+'_signal_sum_error_ele'  +'_' + year_label  +'.txt'
        if flav_name == 'muon': name_txt_signal = mass1_name+'_signal_sum_error_muon' + '_' + year_label  +'.txt'


    mass1_name = mass1_name.replace('GeV', ' \GeV')
    mass2_name = mass2_name.replace('GeV', ' \GeV')
    mass3_name = mass3_name.replace('GeV', ' \GeV')
    

    v1_name = v1_name.replace('_', '.' )
    v1_name = v1_name.replace('=', ' \cdot10^{-' )
    v2_name = v2_name.replace('_', '.' )
    v2_name = v2_name.replace('=', ' \cdot10^{-' )
    v3_name = v3_name.replace('_', '.' )
    v3_name = v3_name.replace('=', ' \cdot10^{-' )
    legend_label_signal1 = '$M = ' + mass1_name + ", \lvert V^{2} \\rvert = " + v1_name + '}$       '
    legend_label_signal2 = '$M = ' + mass2_name + ", \lvert V^{2} \\rvert = " + v2_name + '}$       '
    legend_label_signal3 = '$M = ' + mass3_name + ", \lvert V^{2} \\rvert = " + v3_name + '}$       '
    print legend_label_signal1
            
    simb = '\\\\'
    begin_t = '\\begin'       
            
    sys.stdout=open(name_txt,'w')
    #eee/mmm case
    if flav_name == 'muon': print '-------------------------------- mmm -------------------------------- '
    if flav_name == 'ele': print '-------------------------------- eee -------------------------------- '
    print 'first 4 bins M< 4, last 2 > 4GeV. delta 2D:  0-0.5/ 0.5-1.5/1.5-4/>4 '
   
    print ('\n\n\n')
    print '\\begin{table}'
    print '\footnotesize'
    print '\\centering'
    print ' \\caption{\label{tab: TO FILL} blablabla.}'
    print begin_t +'{tabular}{l|c|c|c|c|c|c}'
    print'\hline'
    if  year_label == '16': print' & \multicolumn{6}{c}{\\textbf{Yields (2016)}}' + simb
    if  year_label == '17': print' & \multicolumn{6}{c}{\\textbf{Yields (2017)}}' + simb
    if  year_label == '18': print' & \multicolumn{6}{c}{\\textbf{Yields (2018)}}' + simb
    if  year_label == 'combined': print'& \multicolumn{6}{c}{\\textbf{Yields (Full Run2)}}' + simb
    print'\cline{2-7}'
    print'\multicolumn{1}{r}{\mtwol}    & \multicolumn{4}{|c|}{$< 4 \GeV$} & \multicolumn{2}{c}{$> 4 \GeV$}' + simb
    print'\cline{2-7}'
    print'\multicolumn{1}{r|}{\Deltwod (cm)} & $[0,0.5]$ & $[0.5,1.5]$ & $[1.5,4]$ &  $>4$ &  $[0,0.5]$ &  $>0.5$' +simb
    print'\hline'
    if flav_name == 'ele': print'$\Pe\Pe\Pe$ & & & & & & ' +simb
    if flav_name == 'muon': print'$\PGm\PGm\PGm$   & & & & & &  ' +simb
    print'Total bkg        ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        tt =  str("%.2g" %bkg_total_syst.GetBinContent( b ) if bkg_total_syst.GetBinContent( b ) > 0.01 else 0)
        ttt = str("%.1g" %math.sqrt(bkg_total_onlystat.GetBinError(b)* bkg_total_onlystat.GetBinError(b)  + bkg_total_onlysyst.GetBinError( b )*bkg_total_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        if tt == '1e-06' : tt = 0        
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $     '),
    print simb
    print 'Observed        ',
    for g in range( 1, data.GetNbinsX() + 1 ):
        if g > 6: continue
        tt =  str("%.2g" %data.GetBinContent( g ))
        ttt = str("%.1g" %data.GetBinError(g))      
        print ('&  $  '),
        if tt == '1e-06' : tt = 0        
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $     '),
    print simb
    print'\hline'
    if flav_name == 'ele': print ' $\Pe^\pm\Pe^\mp\PGm^\pm$  & & & & & & ' +simb
    if flav_name == 'muon': print ' $\PGm^\pm\PGm^\mp\Pe^\pm$  & & & & & & ' +simb      
    print'Total bkg       ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        tt =  str("%.2g" %bkg_total_syst.GetBinContent( b ) if bkg_total_syst.GetBinContent( b ) > 0.01 else 0)
        ttt = str("%.1g" %math.sqrt(bkg_total_onlystat.GetBinError(b)* bkg_total_onlystat.GetBinError(b)  + bkg_total_onlysyst.GetBinError( b )*bkg_total_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        if tt == '1e-06' : tt = 0                
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $    '),
    print simb
    print 'Observed        ',
    for b in range( 1, data.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        tt =  str("%.2g" %data.GetBinContent( b ))
        ttt = str("%.1g" %data.GetBinError(b))      
        print ('&  $  '),
        if tt == '1e-06' : tt = 0        
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $     '),
    print simb
    print'\hline'
    if flav_name == 'ele': print '$\Pe^\pm\Pe^\pm\PGm^\mp$  & & & & & & ' +simb
    if flav_name == 'muon': print '$\PGm^\pm\PGm^\pm\Pe^\mp$  & & & & & & ' +simb       
    print'Total bkg        ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        tt =  str("%.2g" %bkg_total_syst.GetBinContent( b ) if bkg_total_syst.GetBinContent( b ) > 0.01 else 0)
        ttt = str("%.1g" %math.sqrt(bkg_total_onlystat.GetBinError(b)* bkg_total_onlystat.GetBinError(b)  + bkg_total_onlysyst.GetBinError( b )*bkg_total_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        if tt == '1e-06' : tt = 0                
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $     '),
    print simb
    print 'Observed        ',
    for b in range( 1, data.GetNbinsX() + 1 ):
        if b <= 12: continue
        tt =  str("%.2g" %data.GetBinContent( b ))
        ttt = str("%.1g" %data.GetBinError(b))      
        print ('&  $  '),
        if tt == '1e-06' : tt = 0        
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $     '),
    print simb
    print'\hline'
    print'\end{tabular}'
    print'\end{table}'
    sys.stdout.close()


    if flav_name == 'ele':   label_xxx = '$\Pe\Pe\Pe$   '
    if flav_name == 'muon':  label_xxx = '$\PGm\PGm\PGm$    '
    if flav_name == 'ele':   label_xxOS = '$\Pe^\pm\Pe^\mp\PGm^\pm$    '
    if flav_name == 'muon':  label_xxOS = '$\PGm^\pm\PGm^\mp\Pe^\pm$    '
    if flav_name == 'ele':   label_xxSS = '$\Pe^\pm\Pe^\pm\PGm^\mp$    '
    if flav_name == 'muon':  label_xxSS = '$\PGm^\pm\PGm^\pm\Pe^\mp$   '

            
            
    sys.stdout=open(name_txt_signal,'w')
    #eee/mmm case
    if flav_name == 'muon': print '-------------------------------- mmm -------------------------------- '
    if flav_name == 'ele': print '-------------------------------- eee -------------------------------- '
    print 'first 4 bins M< 4, last 2 > 4GeV. delta 2D:  0-0.5/ 0.5-1.5/1.5-4/>4 '
   
    print ('\n\n\n')
    print '\ begin{table}'
    print '\footnotesize'
    print '\ centering'
    print ' \ caption{\label{tab: TO FILL} blablabla.}'
    print begin_t +'{tabular}{l|c|c|c|c|c|c}'
    print'\hline'
    if  year_label == '16': print' & \multicolumn{6}{c}{\\textbf{Yields (2016)}}' + simb
    if  year_label == '17': print' & \multicolumn{6}{c}{\\textbf{Yields (2017)}}' + simb
    if  year_label == '18': print' & \multicolumn{6}{c}{\\textbf{Yields (2018)}}' + simb
    if  year_label == 'combined': print'& \multicolumn{6}{c}{\\textbf{Yields (Full Run2)}}' + simb
    print'\cline{2-7}'
    print'\multicolumn{1}{r}{\mtwol}    & \multicolumn{4}{|c|}{$< 4 \GeV$} & \multicolumn{2}{c}{$> 4 \GeV$}' +simb
    print'\cline{2-7}'
    print'\multicolumn{1}{r|}{\Deltwod (cm)} & $[0,0.5]$ & $[0.5,1.5]$ & $[1.5,4]$ &  $>4$ &  $[0,0.5]$ &  $>0.5$' +simb
    print'\hline'
    print legend_label_signal1 +'  & & & & & & ' +simb   
    print label_xxx ,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        ss =  str("%.2g" %signal_histo1.GetBinContent( b ) if signal_histo1.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo1_onlystat.GetBinError(b)* signal_histo1_onlystat.GetBinError(b)  + signal_histo1_onlysyst.GetBinError( b )*signal_histo1_onlysyst.GetBinError( b )) if math.sqrt(signal_histo1_onlystat.GetBinError(b)* signal_histo1_onlystat.GetBinError(b)  + signal_histo1_onlysyst.GetBinError( b )*signal_histo1_onlysyst.GetBinError( b )) > 0.001 else 0)    
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print ('&  $  '),
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $     '),    
    print simb
    print label_xxOS,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        ss =  str("%.2g" %signal_histo1.GetBinContent( b ) if signal_histo1.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo1_onlystat.GetBinError(b)* signal_histo1_onlystat.GetBinError(b)  + signal_histo1_onlysyst.GetBinError( b )*signal_histo1_onlysyst.GetBinError( b )) if math.sqrt(signal_histo1_onlystat.GetBinError(b)* signal_histo1_onlystat.GetBinError(b)  + signal_histo1_onlysyst.GetBinError( b )*signal_histo1_onlysyst.GetBinError( b )) > 0.001 else 0)    
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print simb
    print label_xxSS,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        ss =  str("%.2g" %signal_histo1.GetBinContent( b ) if signal_histo1.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo1_onlystat.GetBinError(b)* signal_histo1_onlystat.GetBinError(b)  + signal_histo1_onlysyst.GetBinError( b )*signal_histo1_onlysyst.GetBinError( b )) if math.sqrt(signal_histo1_onlystat.GetBinError(b)* signal_histo1_onlystat.GetBinError(b)  + signal_histo1_onlysyst.GetBinError( b )*signal_histo1_onlysyst.GetBinError( b )) > 0.001 else 0)    
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print simb
    print'\hline'
    print legend_label_signal2 +'  & & & & & &' +simb  
    print label_xxx ,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        ss =  str("%.2g" %signal_histo2.GetBinContent( b ) if signal_histo2.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo2_onlystat.GetBinError(b)* signal_histo2_onlystat.GetBinError(b)  + signal_histo2_onlysyst.GetBinError( b )*signal_histo2_onlysyst.GetBinError( b )) if math.sqrt(signal_histo2_onlystat.GetBinError(b)* signal_histo2_onlystat.GetBinError(b)  + signal_histo2_onlysyst.GetBinError( b )*signal_histo2_onlysyst.GetBinError( b )) > 0.001 else 0)
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $     '),    
    print simb
    print label_xxOS,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        ss =  str("%.2g" %signal_histo2.GetBinContent( b ) if signal_histo2.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo2_onlystat.GetBinError(b)* signal_histo2_onlystat.GetBinError(b)  + signal_histo2_onlysyst.GetBinError( b )*signal_histo2_onlysyst.GetBinError( b )) if math.sqrt(signal_histo2_onlystat.GetBinError(b)* signal_histo2_onlystat.GetBinError(b)  + signal_histo2_onlysyst.GetBinError( b )*signal_histo2_onlysyst.GetBinError( b )) > 0.001 else 0)
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print simb
    print label_xxSS,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        ss =  str("%.2g" %signal_histo2.GetBinContent( b ) if signal_histo2.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo2_onlystat.GetBinError(b)* signal_histo2_onlystat.GetBinError(b)  + signal_histo2_onlysyst.GetBinError( b )*signal_histo2_onlysyst.GetBinError( b )) if math.sqrt(signal_histo2_onlystat.GetBinError(b)* signal_histo2_onlystat.GetBinError(b)  + signal_histo2_onlysyst.GetBinError( b )*signal_histo2_onlysyst.GetBinError( b )) > 0.001 else 0)
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print simb
    print'\hline'
    print legend_label_signal3 +'  & & & & & & ' +simb   
    print label_xxx ,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        ss =  str("%.2g" %signal_histo3.GetBinContent( b ) if signal_histo3.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo3_onlystat.GetBinError(b)* signal_histo3_onlystat.GetBinError(b)  + signal_histo3_onlysyst.GetBinError( b )*signal_histo3_onlysyst.GetBinError( b )) if math.sqrt(signal_histo3_onlystat.GetBinError(b)* signal_histo3_onlystat.GetBinError(b)  + signal_histo3_onlysyst.GetBinError( b )*signal_histo3_onlysyst.GetBinError( b )) > 0.001 else 0)
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $     '),    
    print simb
    print label_xxOS,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        ss =  str("%.2g" %signal_histo3.GetBinContent( b ) if signal_histo3.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo3_onlystat.GetBinError(b)* signal_histo3_onlystat.GetBinError(b)  + signal_histo3_onlysyst.GetBinError( b )*signal_histo3_onlysyst.GetBinError( b )) if math.sqrt(signal_histo3_onlystat.GetBinError(b)* signal_histo3_onlystat.GetBinError(b)  + signal_histo3_onlysyst.GetBinError( b )*signal_histo3_onlysyst.GetBinError( b )) > 0.001 else 0)
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print simb
    print label_xxSS,
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        ss =  str("%.2g" %signal_histo3.GetBinContent( b ) if signal_histo3.GetBinContent( b ) > 0.01 else 0)
        sss = str("%.1g" %math.sqrt(signal_histo3_onlystat.GetBinError(b)* signal_histo3_onlystat.GetBinError(b)  + signal_histo3_onlysyst.GetBinError( b )*signal_histo3_onlysyst.GetBinError( b )) if math.sqrt(signal_histo3_onlystat.GetBinError(b)* signal_histo3_onlystat.GetBinError(b)  + signal_histo3_onlysyst.GetBinError( b )*signal_histo3_onlysyst.GetBinError( b )) > 0.001 else 0)
        print ('&  $  '),
        if ss == '1e-06' : ss = 0
        if sss == '1e-06' : sss = 0
        if ss == '3e-06' : ss = 0
        if sss == '1.7e-06' : sss = 0            
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print simb
    print'\hline'
    print'\end{tabular}'
    print'\end{table}'
    sys.stdout.close()





