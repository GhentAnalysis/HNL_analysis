
import sys
import os
import time

from ROOT import TFile, TH1, TDirectory
from ROOT import TROOT
import ROOT

from Uncertainty import Uncertainty
from Process import Process
from Datacard import Datacard
from ProcessCollection import ProcessCollection
from style import setStyle
from drawPlot import drawPlot
from drawPlot import writeTables_sumErrors

def buildHistogramDictionary( root_file_path ):
    root_file = TFile( root_file_path )
    histogram_dict = {}
    for entry in root_file.GetListOfKeys():
        name = entry.GetName()
        hist = root_file.Get( name )
        hist.SetDirectory( 0 )
        if isinstance( hist, TH1 ) :
            histogram_dict[ name ] = hist
    return histogram_dict

if __name__ == '__main__':

    #set plotting style
    style = setStyle()

    if len( sys.argv ) > 1: base_card_path = sys.argv[1]

    if len( sys.argv ) > 2: base_card_path2 = sys.argv[2]     

    if len( sys.argv ) > 3: base_card_path3 = sys.argv[3]
        
    if len( sys.argv ) > 4: base_card_path4 = sys.argv[4]     

    if len( sys.argv ) > 5: flav = [ sys.argv[5] ]
    else: flav = ['muon', 'ele']

    if len( sys.argv ) > 6: plot = [ sys.argv[6] ]
    else: plot = ['SR', 'M', 'D','M3' ]

    if len( sys.argv ) > 7: mass1 = [ sys.argv[7] ]
    else: mass1 = ['5GeV']
    
    if len( sys.argv ) > 8: mass2 = [ sys.argv[8] ]
    else: mass2 = ['5GeV']
    if len( sys.argv ) > 9: mass3 = [ sys.argv[9] ]
    else: mass3 = ['5GeV']        

    if len( sys.argv ) > 10: v1 = [ sys.argv[10] ]
    else: v1 = ['01']
    if len( sys.argv ) > 11: v2 = [ sys.argv[11] ]
    else: v2 = ['01']
    if len( sys.argv ) > 12: v3 = [ sys.argv[12] ]
    else: v3 = ['01']

    if len( sys.argv ) > 13: years = [ sys.argv[13] ]
    if sys.argv[13] == 'Run2': years = ['16', '17', '18' ]
    
    if len( sys.argv ) > 14: cr_card_path = sys.argv[14]

    if flav[0] == 'muon':
        base_card_path = base_card_path.replace( '_mu_muo', '_mu_muo' )
        base_card_path2 = base_card_path2.replace( '_mu_muo', '_mu_muo' )
        base_card_path3 = base_card_path3.replace( '_mu_muo', '_mu_muo' )
        base_card_path4 = base_card_path4.replace( '_mu_muo', '_mu_muo' )
        cr_card_path = cr_card_path.replace( '_mu_muo', '_mu_muo' )
    if flav[0] == 'ele':
        base_card_path = base_card_path.replace( '_mu_muo', '_e_ele' )
        base_card_path2 = base_card_path2.replace( '_mu_muo', '_e_ele' )
        base_card_path3 = base_card_path3.replace( '_mu_muo', '_e_ele' )
        base_card_path4 = base_card_path4.replace( '_mu_muo', '_e_ele' )
        cr_card_path = cr_card_path.replace( '_mu_muo', '_e_ele' )
        
    if plot[0] == 'M':
        base_card_path = base_card_path.replace( '_datacard', '_mass_datacard' )
        base_card_path2 = base_card_path2.replace( '_datacard', '_mass_datacard' )
        base_card_path3 = base_card_path3.replace( '_datacard', '_mass_datacard' )
        base_card_path4 = base_card_path4.replace( '_datacard', '_mass_datacard' )
        cr_card_path = cr_card_path.replace( '_datacard', '_mass_datacard' )
    elif plot[0] == 'D':
        base_card_path = base_card_path.replace( '_datacard', '_disp_datacard' )
        base_card_path2 = base_card_path2.replace( '_datacard', '_disp_datacard' )
        base_card_path3 = base_card_path3.replace( '_datacard', '_disp_datacard' )
        base_card_path4 = base_card_path4.replace( '_datacard', '_disp_datacard' )
        cr_card_path = cr_card_path.replace( '_datacard', '_disp_datacard' )
    elif plot[0] == 'M3':
        base_card_path = base_card_path.replace( '_datacard', '_mass3_datacard' )
        base_card_path2 = base_card_path2.replace( '_datacard', '_mass3_datacard' )
        base_card_path3 = base_card_path3.replace( '_datacard', '_mass3_datacard' )
        base_card_path4 = base_card_path4.replace( '_datacard', '_mass3_datacard' )
        cr_card_path = cr_card_path.replace( '_datacard', '_mass3_datacard' )
    elif plot[0] == 'SR':
        base_card_path = base_card_path.replace( '_datacard', '_datacard' )
        base_card_path2 = base_card_path2.replace( '_datacard', '_datacard' )    
        base_card_path3 = base_card_path3.replace( '_datacard', '_datacard' )    
        base_card_path4 = base_card_path4.replace( '_datacard', '_datacard' )
        cr_card_path = cr_card_path.replace( '_datacard', '_datacard' )
 
    datacard_paths = [ '{}_{}.txt'.format( base_card_path, y ) for y in years ]
    datacards = [ Datacard( path ) for path in datacard_paths ]
    #additional signal sample to add to the plot
    datacard_paths2 = [ '{}_{}.txt'.format( base_card_path2, y ) for y in years ]
    datacards2 = [ Datacard( path ) for path in datacard_paths2 ]
    datacard_paths3 = [ '{}_{}.txt'.format( base_card_path3, y ) for y in years ]
    datacards3 = [ Datacard( path ) for path in datacard_paths3 ]
    datacard_paths4 = [ '{}_{}.txt'.format( base_card_path4, y ) for y in years ]
    datacards4 = [ Datacard( path ) for path in datacard_paths4 ]

    crf, sf, df = [], [], []
    for y in years:
        if flav[0] == 'muon': crf.append(TFile( cr_card_path+'/muon_'+plot[0]+'_'+y+'.root', 'READ' ))
        else: crf.append(TFile( cr_card_path+'/ele_'+plot[0]+'_'+y+'.root', 'READ' ))
        sf.append(crf[-1].Get('nonpromptSFpos').Clone('sf'+y))
        df.append(crf[-1].Get('nonpromptDF').Clone('df'+y))
    
    processCollections = []
    total_data = None
    for i, card in enumerate( datacards ):  
        #build the backgrounds
        histogram_dict = buildHistogramDictionary( card.shapeFilePath() )
        process_list = []
        for p in card.processNames():
            process_list.append( Process( p, histogram_dict, card, 'signal' in p ) )
        processCollections.append( ProcessCollection( process_list, sf = sf[i], df = df[i] ) )
        #merge the data
        if total_data is None:
            total_data = histogram_dict[ 'data_obs' ]
        else :
            total_data.Add( histogram_dict[ 'data_obs' ] )   
    total_process_collection = processCollections[0]
    for col in processCollections[1:]:
        total_process_collection += col

    processCollections2 = []
    for i, card2 in enumerate( datacards2 ):  
        #build the backgrounds
        histogram_dict2 = buildHistogramDictionary( card2.shapeFilePath() )
        process_list2 = []
        for p in card2.processNames():
            process_list2.append( Process( p, histogram_dict2, card2, 'signal' in p ) )
        processCollections2.append( ProcessCollection( process_list2, sf = sf[i], df = df[i] ) )
    total_process_collection2 = processCollections2[0]
    for col in processCollections2[1:]:
        total_process_collection2 += col

    processCollections3 = []
    for i, card3 in enumerate( datacards3 ):  
        #build the backgrounds
        histogram_dict3 = buildHistogramDictionary( card3.shapeFilePath() )
        process_list3 = []
        for p in card3.processNames():
            process_list3.append( Process( p, histogram_dict3, card3, 'signal' in p ) )
        processCollections3.append( ProcessCollection( process_list3, sf = sf[i], df = df[i] ) )
    total_process_collection3 = processCollections3[0]
    for col in processCollections3[1:]:
        total_process_collection3 += col

    processCollections4 = []
    for i, card4 in enumerate( datacards4 ):  
        #build the backgrounds
        histogram_dict4 = buildHistogramDictionary( card4.shapeFilePath() )
        process_list4 = []
        for p in card4.processNames():
            process_list4.append( Process( p, histogram_dict4, card4, 'signal' in p ) )
        processCollections4.append( ProcessCollection( process_list4, sf = sf[i], df = df[i] ) )
    total_process_collection4 = processCollections4[0]
    for col in processCollections4[1:]:
        total_process_collection4 += col
        
   #correct plot label for each variable:
    if len( plot ) == 1:
        if plot[0] == 'SR': plot_label = 'SR'
        elif plot[0] == 'M': plot_label = 'massl2l3'
        elif plot[0] == 'D': plot_label = 'displacement'
        elif plot[0] == 'M3': plot_label = 'mass3'    

    #correct plot label for each variable:
    if len( flav ) == 1:
        if flav[0] == 'muon': flav_label = 'muon'
        elif flav[0] == 'ele': flav_label = 'ele'
        
    #correct lumi label for each year:
    if len( years ) == 1:
        if years[0] == '16': lumi_label = '35.9 fb^{-1} (13 TeV)'
        elif years[0] == '17': lumi_label = '41.5 fb^{-1} (13 TeV)'
        elif years[0] == '18': lumi_label = '59.7 fb^{-1} (13 TeV)'
    elif len( years ) == 3: lumi_label = '137 fb^{-1} (13 TeV)'

    #correct signal legend:
    if len( mass1 ) == 1: mass1_label = mass1[0]
    else: print( 'no mass1' )
    if len( mass2 ) == 1: mass2_label = mass2[0]
    else: print( 'no mass2' )
    if len( mass3 ) == 1: mass3_label = mass3[0]
    else: print( 'no mass3' )
    if len( v1 ) == 1: v1_label = v1[0]
    else: print( 'no v1' )
    if len( v2 ) == 1: v2_label = v2[0]
    else: print( 'no v2' )
    if len( v3 ) == 1: v3_label = v3[0]
    else: print( 'no v3' )        
               
    year_name = years[0] if len( years ) == 1 else 'combined'
   
    base_card_path = base_card_path.replace( '.', 'p' )
    plot_name = os.path.basename( base_card_path ) + '_' + year_name + '_' + plot_label
    
#    colors = { 'nonpromptDF' : ROOT.kGreen, 'nonpromptSF' : ROOT.kGreen + 2, 'Xgamma' : ROOT.kViolet+1 }
    colors = { 'nonpromptDF' : ROOT.kGreen-9, 'nonpromptSF' : ROOT.kGreen + 1, 'Xgamma' : ROOT.kYellow-8 }
    legend_names = { 'nonpromptDF' : 'Nonprompt (DF)', 'nonpromptSF' : 'Nonprompt (SF)', 'Xgamma' : 'Conversions' }
    
    drawPlot( total_data, total_process_collection, \
    total_process_collection2, total_process_collection3, \
    total_process_collection4, plot_name, colors, \
    legend_names = legend_names, lumi_text = lumi_label, \
    plot_axis= plot_label, flav_name = flav_label, mass1_name = mass1_label, \
    mass2_name = mass2_label,mass3_name = mass3_label, v1_name = v1_label, v2_name = v2_label,v3_name = v3_label)
    
    if plot_label == 'SR': writeTables_sumErrors( total_data, total_process_collection,total_process_collection2,total_process_collection3, total_process_collection4, plot_name, plot_axis= plot_label, flav_name = flav_label, year_label = year_name, mass1_name = mass1_label, mass2_name = mass2_label,mass3_name = mass3_label, v1_name = v1_label, v2_name = v2_label,v3_name = v3_label)
