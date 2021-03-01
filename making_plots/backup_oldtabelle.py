
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


    lower_pad = TPad( "", "", 0, 0, 1, lower_pad_fraction )
    lower_pad.SetTopMargin( 0.01 )
    lower_pad.SetBottomMargin( 0.4 )

    return c, upper_pad, lower_pad


def makeUpperLegend( data, processCollection, bkg_total, legend_names = None ):

    legend = TLegend(0.20, 0.73, 0.92, 0.90, '', 'brNDC');
    legend.SetNColumns( 3 )
    legend.SetFillStyle( 0 ) #avoid box
    legend.SetBorderSize(0)
    
    legend.AddEntry( data, 'Total bkg', 'pe1' )
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


def rangeLog( data, background,signal, flav_name  ):
	
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
    if flav_name == 'muon' :
        total_max = total_max * 10**( 0.8 * number_of_orders )
    if flav_name == 'ele' :
        total_max = total_max * 10**( 0.8 * number_of_orders )
        
    return total_min, total_max


def setRelativeUncStyle( relative_unc, color ):
	relative_unc.SetFillStyle( 1001 )
	relative_unc.SetMarkerStyle( 1 )
	relative_unc.SetFillColor( color )


def setRatioPlotStyle( first_hist, lower_pad_fraction,plot_axis ):
    scale_factor = ( 1. - lower_pad_fraction ) / lower_pad_fraction 
    first_hist.GetYaxis().SetTitleOffset( 0.8 / scale_factor )
    first_hist.GetXaxis().SetTitleOffset( 3.0 / scale_factor )
    first_hist.GetYaxis().SetTitleSize( scale_factor * .06 )
    first_hist.GetXaxis().SetTitleSize( scale_factor * .06 )
    first_hist.GetYaxis().SetLabelSize( scale_factor * .04 )
    first_hist.GetXaxis().SetLabelSize( scale_factor * .05 )
    first_hist.GetXaxis().LabelsOption("v");
    first_hist.GetXaxis().SetLabelOffset( scale_factor * 0.009 )
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
    legend = TLegend( 0.18, 0.85, 0.94, 0.98, '', 'brNDC' )
    legend.SetNColumns( 3 )
    legend.SetFillStyle( 0 )
    legend.SetBorderSize( 0 )
    legend.AddEntry( relative_bkg_stat_unc, 'Stat. pred. unc.', 'f' )
    legend.AddEntry( relative_bkg_total_unc, 'Total pred. unc.', 'f' )
    legend.AddEntry( obs_over_pred, 'Obs./Pred.', 'pe12' )
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


def writeTables_sumErrors(data, processCollection, processCollection2, processCollection3,processCollection4, plot_name,plot_axis = 'SR',flav_name = 'muon', year_label = '16'):

 
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


    signal_histo =  processCollection.signal_get()
    signal_histo_onlysyst = processCollection.signalWithSystErrors()
    signal_histo_onlystat = processCollection.signalWithStatErrors() 

    name_txt = 'table'
    if year_label == 'combined':
        if flav_name == 'ele': name_txt = 'sum_error_ele_combined.txt'
        if flav_name == 'muon': name_txt = 'sum_error_muon_combined.txt'
    else:
        if flav_name == 'ele': name_txt = 'sum_error_ele'  +'_' + year_label  +'.txt'
        if flav_name == 'muon': name_txt = 'sum_error_muon' + '_' + year_label  +'.txt'
            
    sys.stdout=open(name_txt,'w')
    #eee/mmm case
    if flav_name == 'muon': print '-------------------------------- mmm -------------------------------- '
    if flav_name == 'ele': print '-------------------------------- eee -------------------------------- '
    print 'first 4 bins M< 4, last 2 > 4GeV. delta 2D:  0-0.5/ 0.5-1.5/1.5-4/>4 '
   
    print ('\n\n\n')
    print '\ begin{table}'
    print '\footnotesize'
    print '\ centering'
    print ' \ caption{\label{tab: TO FILL} blablabla.}'
    print ' \ begin{tabular}{l|c|c|c|c|c|c}'
    print'\hline'
    if  year_label == '16': print'\multirow{3}{*}{Process} & \multicolumn{6}{c}{\ textbf{Yields (2016)}}' + '\\'
    if  year_label == '17': print'\multirow{3}{*}{Process} & \multicolumn{6}{c}{\ textbf{Yields (2017)}}' + '\\'
    if  year_label == '18': print'\multirow{3}{*}{Process} & \multicolumn{6}{c}{\ textbf{Yields (2018)}}' + '\\'
    if  year_label == 'combined': print'\multirow{3}{*}{Process} & \multicolumn{6}{c}{\ textbf{Yields (Full Run2)}}' + '\\'
    print'\cline{2-7}'
    print'& \multicolumn{4}{|c|}{\mtwol $< 4$} & \multicolumn{2}{|c|}{\mtwol $> 4$}  \\ '
    print'\cline{2-7}'
    print'&\Deltwod (cm) $\ rightarrow$ $[0,0.5]$ & $[0.5,1.5]$ & $[1.5,4]$ &  $>4$ &  $[0,0.5]$ &  $>0.5$  \\ '
    print'\hline'
    print'\hline'
    if flav_name == 'ele': print'$\Pe\Pe\Pe$ & \multicolumn{6}{c}{} \\ '
    if flav_name == 'muon': print'$\PGm\PGm\PGm$  & \multicolumn{6}{c}{} \\ '
    print'\hline'
    print'Total bkg        ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        tt =  str("%.2g" %bkg_total_syst.GetBinContent( b ))
        ttt = str("%.2g" %math.sqrt(bkg_total_onlystat.GetBinError(b)* bkg_total_onlystat.GetBinError(b)  + bkg_total_onlysyst.GetBinError( b )*bkg_total_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $     '),
    print '\\ '
    print'\hline'
    print 'Observed        & 0.&0.&0.&0.&0.&0. \\ '
    print'\hline'
    print '$M = 2 GeV, \lvert V^{2} \ rvert = 1*10^{-4}$       ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        ss =  str("%.2g" %signal_histo.GetBinContent( b ))
        sss = str("%.2g" %math.sqrt(signal_histo_onlystat.GetBinError(b)* signal_histo_onlystat.GetBinError(b)  + signal_histo_onlysyst.GetBinError( b )*signal_histo_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $     '),    
    print '\\'
    print'\hline'
    print'\hline'
    if flav_name == 'ele': print ' $\Pe^\pm\Pe^\mp\PGm^\pm$ & \multicolumn{6}{c}{} \\'
    if flav_name == 'muon': print ' $\PGm^\pm\PGm^\mp\Pe^\pm$ & \multicolumn{6}{c}{} \\'       
    print'\hline'
    print'Total bkg       ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        tt =  str("%.2g" %bkg_total_syst.GetBinContent( b ))
        ttt = str("%.2g" %math.sqrt(bkg_total_onlystat.GetBinError(b)* bkg_total_onlystat.GetBinError(b)  + bkg_total_onlysyst.GetBinError( b )*bkg_total_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $    '),
    print '\\ '
    print'\hline'
    print 'Observed        & 0.&0.&0.&0.&0.&0. \\ '
    print'\hline'
    print '$M = 2 GeV, \lvert V^{2} \ rvert = 1*10^{-4}$       ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        ss =  str("%.2g" %signal_histo.GetBinContent( b ))
        sss = str("%.2g" %math.sqrt(signal_histo_onlystat.GetBinError(b)* signal_histo_onlystat.GetBinError(b)  + signal_histo_onlysyst.GetBinError( b )*signal_histo_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print '\\'
    print'\hline'
    print'\hline'
    if flav_name == 'ele': print '$\Pe^\pm\Pe^\pm\PGm^\mp$ & \multicolumn{6}{c}{} \\'
    if flav_name == 'muon': print '$\PGm^\pm\PGm^\pm\Pe^\mp$ & \multicolumn{6}{c}{} \\'        
    print'\hline'
    print'Total bkg        ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        tt =  str("%.2g" %bkg_total_syst.GetBinContent( b ))
        ttt = str("%.2g" %math.sqrt(bkg_total_onlystat.GetBinError(b)* bkg_total_onlystat.GetBinError(b)  + bkg_total_onlysyst.GetBinError( b )*bkg_total_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        print (tt),
        print ('  \pm   '),
        print (ttt),
        print (' $     '),
    print '\\ '
    print'\hline'
    print 'Observed        & 0.&0.&0.&0.&0.&0. \\ '
    print'\hline'
    print '$M = 2 GeV, \lvert V^{2} \ rvert = 1*10^{-4}$      ',
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        ss =  str("%.2g" %signal_histo.GetBinContent( b ))
        sss = str("%.2g" %math.sqrt(signal_histo_onlystat.GetBinError(b)* signal_histo_onlystat.GetBinError(b)  + signal_histo_onlysyst.GetBinError( b )*signal_histo_onlysyst.GetBinError( b )))      
        print ('&  $  '),
        print (ss),
        print ('  \pm   '),
        print (sss),
        print (' $      '),    
    print '\\'
    print'\hline'
    print'\hline'
    print'\end{tabular}'
    print'\end{table}'
    sys.stdout.close()


def writeTables(data, processCollection, processCollection2, processCollection3,processCollection4, plot_name,plot_axis = 'SR',flav_name = 'muon', year_label = '16'):
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


    signal_histo =  processCollection.signal_get()
   

    name_txt = 'table'
    if year_label == 'combined':
        if flav_name == 'ele': name_txt = 'ele_combined.txt'
        if flav_name == 'muon': name_txt = 'muon_combined.txt'
    else:
        if flav_name == 'ele': name_txt = 'ele'  +'_' + year_label  +'.txt'
        if flav_name == 'muon': name_txt = 'muon' + '_' + year_label  +'.txt'
            
    sys.stdout=open("name_txt",'w')
    #eee/mmm case
    if flav_name == 'muon': print '-------------------------------- mmm -------------------------------- '
    if flav_name == 'ele': print '-------------------------------- eee -------------------------------- '
    print 'first 4 bins M< 4, last 2 > 4GeV. delta 2D:  0-0.5/ 0.5-1.5/1.5-4/>4 '
    print 'prompt bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        s = repr("%.3f" %bkg_p_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_p_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_p_onlysyst.GetBinError( b )) +' $  &  $'  
        print s
    print ('\n')
    print 'non prompt bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        t = repr("%.3f" %bkg_np_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_np_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_np_onlysyst.GetBinError( b )) +' $  &  $'  
        print t
    print ('\n')
    print 'tot bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        tt = repr("%.3f" %bkg_total_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_total_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_total_onlysyst.GetBinError( b )) +' $  &  $'  
        print tt
    print ('\n')
    print 'obs'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        g = repr(0) + '  \pm  ' + repr(0) + '  \pm  ' + repr(0) +' $  &  $'  
        print g
    print ('\n')
    print 'signal'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b > 6: continue
        ss = repr("%.3f" %signal_histo.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %signal_histo.GetBinError( b ))  +' $  &  $'
        print ss
    print ('\n')
     #sys.stdout.close()
    print ('\n\n\n')
    
    if flav_name == 'muon': print '-------------------------------- mme OS -------------------------------- '
    if flav_name == 'ele': print '-------------------------------- eem OS -------------------------------- '
    print 'first 4 bins M< 4, last 2 > 4GeV. delta 2D:  0-0.5/ 0.5-1.5/1.5-4/>4 '
    print 'prompt bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        s = repr("%.3f" %bkg_p_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_p_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_p_onlysyst.GetBinError( b )) +' $  &  $'  
        print s
    print ('\n')
    print 'non prompt bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        t = repr("%.3f" %bkg_np_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_np_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_np_onlysyst.GetBinError( b )) +' $  &  $'  
        print t
    print ('\n')
    print 'tot bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        tt = repr("%.3f" %bkg_total_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_total_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_total_onlysyst.GetBinError( b )) +' $  &  $'  
        print tt
    print ('\n')
    print 'obs'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        g = repr(0) + '  \pm  ' + repr(0) + '  \pm  ' + repr(0) +' $  &  $'  
        print g
    print ('\n')
    print 'signal'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 6 or b>12: continue
        ss = repr("%.3f" %signal_histo.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %signal_histo.GetBinError( b ))  +' $  &  $'
        print ss
    print ('\n')
     #sys.stdout.close()
    print ('\n\n\n')

    if flav_name == 'muon': print '-------------------------------- mme SS -------------------------------- '
    if flav_name == 'ele': print '-------------------------------- eem SS -------------------------------- '
    print 'first 4 bins M< 4, last 2 > 4GeV. delta 2D:  0-0.5/ 0.5-1.5/1.5-4/>4 '
    print 'prompt bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        s = repr("%.3f" %bkg_p_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_p_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_p_onlysyst.GetBinError( b )) +' $  &  $'  
        print s
    print ('\n')
    print 'non prompt bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        t = repr("%.3f" %bkg_np_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_np_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_np_onlysyst.GetBinError( b )) +' $  &  $'  
        print t
    print ('\n')
    print 'tot bkg'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        tt = repr("%.3f" %bkg_total_syst.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %bkg_total_onlystat.GetBinError( b )) + '  \pm  ' + repr("%.3f" %bkg_total_onlysyst.GetBinError( b )) +' $  &  $'  
        print tt
    print ('\n')
    print 'obs'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        g = repr(0) + '  \pm  ' + repr(0) + '  \pm  ' + repr(0) +' $  &  $'  
        print g
    print ('\n')
    print 'signal'
    for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
        if b <= 12: continue
        ss = repr("%.3f" %signal_histo.GetBinContent( b )) + '  \pm  ' + repr("%.3f" %signal_histo.GetBinError( b ))  +' $  &  $'
        print ss
    print ('\n')
     #sys.stdout.close()
    print ('\n\n\n')
    sys.stdout.close()








    
        
def drawPlot( data, processCollection, processCollection2, processCollection3,processCollection4, plot_name, color_dict = None, log = True, lower_pad_fraction = 0.25, legend_names = None, lumi_text = '137 fb^{-1} (13 TeV)',plot_axis = 'SR',flav_name = 'muon',mass1_name = '1gev',mass2_name = '2gev',mass3_name = '3gev',v1_name = 'v1',v2_name = 'v2',v3_name = 'v3',width = 800, height = 500, additional_label = 'pippo', additional_label2 = 'pippo'):
    
    if plot_axis == 'SR':
        width = 1000
        height = 600
        log = True
    if plot_axis == 'massl2l3' or plot_axis == 'displacement' or plot_axis == 'mass3':
        width = 700
        height = 600
        log = False

  
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
        for b in range( 1, bkg_total_syst.GetNbinsX() + 1 ):
            if bkg_total_syst.GetBinContent( b ) < 1e-8:
                error = bkg_total_syst.GetBinError( b )
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
    
    legend = makeUpperLegend( data_graph, processCollection, bkg_total_syst, legend_names )
    
    #determine plot range 
    if log:
        range_min, range_max = rangeLog( data, bkg_total_syst,signal_histo, flav_name )
    else:
        range_min, range_max = rangeLinear( data, bkg_total_syst, signal_histo )
    bkg_total_syst.SetMinimum( range_min )
    bkg_total_syst.SetMaximum( range_max )
    
    #only draw labels in bottom pad
    bkg_total_syst.GetXaxis().SetLabelSize(0);
    
    bkg_total_syst.Draw( 'e2' )
    stack.Draw( 'histsame' )
    #legend.AddEntry (signal_histo,  'M = 2GeV, |V|^{2} = 1x10^{-4}')
    #legend.AddEntry (signal_histo2,  'M = 4GeV, |V|^{2} = 8x10^{-6}')
    #legend.AddEntry (signal_histo3,  'M = 8GeV, |V|^{2} = 2x10^{-6}')
    legend.AddEntry (signal_histo,  'M = 3GeV, |V|^{2} = 5x10^{-5}')
    legend.AddEntry (signal_histo2,  'M = 4GeV, |V|^{2} = 8x10^{-6}')
    legend.AddEntry (signal_histo3,  'M = 8GeV, |V|^{2} = 2x10^{-6}')
    legend.Draw( 'same' )

    max_flav_line = range_max/40
    max_flav_name = range_max/40
    max_mass_line = range_max/80
    max_mass_name = range_max/80

    if flav_name == 'ele':
        max_flav_line = range_max/15
        max_flav_name = range_max/15
        max_mass_line = range_max/40
        max_mass_name = range_max/40

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

    label11 = makeLabel( "M_{ll} < 4 GeV", 3.6, max_mass_name, 0.03 )
    label22 = makeLabel( "M_{ll} > 4 GeV", 5.4, max_mass_name, 0.03 )
    label33 = makeLabel( "M_{ll} < 4 GeV", 9.6, max_mass_name, 0.03 )
    label44 = makeLabel( "M_{ll} > 4 GeV", 11.4, max_mass_name, 0.03 )
    label55 = makeLabel( "M_{ll} < 4 GeV", 15.6, max_mass_name, 0.03 )
    label66 = makeLabel( "M_{ll} > 4 GeV", 17.4, max_mass_name, 0.03 )
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
    relative_bkg_total_unc.GetYaxis().SetRangeUser( 0.,2.3 )
    relative_bkg_total_unc.GetYaxis().SetTitle( 'Obs./Pred.' )
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
