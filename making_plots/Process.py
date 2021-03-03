

from Uncertainty import Uncertainty, addUncLinear
from Datacard import Datacard


def buildNormUncertainty( nominal_histogram, uncertainty_size ):

    #datacard represents uncertainties as 1 + the size of the uncertainty
    uncertainty_size -= 1.
        
    histogram_down = nominal_histogram.Clone()
    histogram_down.Scale( ( 1. - uncertainty_size ) )

    histogram_up = nominal_histogram.Clone()
    histogram_up.Scale( ( 1. + uncertainty_size ) )
    return Uncertainty( histogram_down, histogram_up )


#combine returns post-fit histograms with bins of width 1, going from 0 to Nbins, so they can't be added to the original histograms in all cases
def setSameRange( old_hist, range_hist ):
    assert( old_hist.GetNbinsX() == range_hist.GetNbinsX() )
    new_hist = range_hist.Clone()
    for b in range( 1, range_hist.GetNbinsX() ):
        new_hist.SetBinContent( b, old_hist.GetBinContent( b ) )
        new_hist.SetBinError( b, old_hist.GetBinError( b ) )
    return new_hist

    

class Process:

    def __init__( self, name, histogram_dict, datacard, is_signal = False ):
        self.__name = name
        self.__is_signal = is_signal
        
        process_variations_down = {}
        process_variations_up = {}
        self.__nominal = None
        for key, hist in histogram_dict.items():
        
            #ignore histograms corresponding to other processes 
            if not key.startswith( self.__name ):
                continue
            
            #the nominal histogram should just have the process name
            if key == self.__name:
                self.__nominal = hist
            else :
            
                #extract name of uncertainty
                unc_name = key.replace( self.__name + '_', '', 1 )
                unc_name = unc_name.split('Down')[0]
                unc_name = unc_name.split('Up')[0]

                #verify that this uncertainty is present for this process, if not skip it
                if datacard.shapeUncertaintySize( unc_name, self.__name ) is None:
                    continue

                if 'Down' in key :
                    process_variations_down[ unc_name ] = hist
                
                elif 'Up' in key:
                    process_variations_up[ unc_name ] = hist 
            
                else:
                    raise KeyError( 'Histogram name "{}" corresponds to neither an up or down variation.'.format( key ) )
        
        #make sure each uncertainty has both a down and up variation
        if sorted( list( process_variations_down.keys() ) ) != sorted( list( process_variations_up.keys() ) ):
            raise KeyError( 'Different uncertainty sources present for up and down variations in histogram dictionary.' )
        
        #convert the separate dictionaries to Uncertainty objects 
        self.__uncertainty_dict = {}
        for unc in process_variations_down.keys():
            self.__uncertainty_dict[ unc ] = Uncertainty( process_variations_down[ unc ], process_variations_up[ unc ] )


        #add normalization uncertainties 
        for unc in datacard.normalizationUncertainties():
            size = datacard.normalizationUncertaintySize( unc, self.__name )
            if size is not None:
                self.__uncertainty_dict[ unc ] = buildNormUncertainty( self.__nominal, size )

        #transform all variations to represent differences from the nominal yield rather than varied yields
        for _, unc in self.__uncertainty_dict.items():
            unc.varDown().Add( self.__nominal, -1. )
            unc.varUp().Add( self.__nominal, -1. )


    def isSignal( self ):
        return self.__is_signal


    def uncertaintySources( self ):
        return self.__uncertainty_dict.keys()


    def hasUncertainty( self, unc_name ):
        return ( unc_name in self.__uncertainty_dict.keys() )


    def nominal( self ):
        return self.__nominal


    def setPostFitNominal( self, post_fit_shape ):
        post_fit_shifted = setSameRange( post_fit_shape, self.nominal() )
        self.__nominal = post_fit_shifted


    def name( self ):
        return self.__name


    def uncertainty( self, source ):
        try:
            return self.__uncertainty_dict[ source ]
        except KeyError:
            raise KeyError( 'Uncertainty source {} is not known for process {}.'.format( source, self.__name ) )


    def scalePostFitUncertainty( self, source, scale_factor ):
        self.uncertainty( source ).scale( scale_factor )
    
    
    def __findVariation( self, source, func ):
        return func( self.uncertainty( source ) )
    

    def varDown( self, source ):
        return self.__findVariation( source, Uncertainty.varDown )
    
    
    def varUp( self, source ):
        return self.__findVariation( source, Uncertainty.varUp )


    def __iadd__( self, rhs ):
        if self.name() != rhs.name():
            raise NameError( 'Processes "{}" and "{}" can not be added. Only processes of the same name can be added.'.format( self.name(), rhs.name() ) )

        #add nominal yields
        self.__nominal.Add( rhs.nominal() )

        for rhs_unc_key, rhs_unc in rhs.__uncertainty_dict.items():

            #add shared uncertainties together ( linear because uncertainties of the same name are correlated )
            try:
                lhs_unc = self.__uncertainty_dict[ rhs_unc_key ]
                self.__uncertainty_dict[ rhs_unc_key ] = addUncLinear( lhs_unc, rhs_unc )
            except KeyError:
                self.__uncertainty_dict[ rhs_unc_key ] = rhs_unc

        return self 


if __name__ == '__main__':
    pass

    
