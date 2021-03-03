
from ROOT import THStack

from Uncertainty import addUncLinear, addUncQuadratic, addHistQuadratic


class ProcessCollection:

    def __init__( self, process_list ):
        self.__processes = process_list
        self.__uncertainties = set()
        for p in self.__processes:
            self.__uncertainties = self.__uncertainties.union( set( p.uncertaintySources() ) )



    def __iter__( self ):
        for p in self.__processes:
            yield p


    def backgroundStack( self ):

        #implement reordering here by total yield 
        bkgStack = THStack()
        for p in reversed( self.__processes ):
            if not p.isSignal():
                bkgStack.Add( p.nominal() )
        return bkgStack


    def signalUncertainty( self, unc_name ):
        if unc_name not in self.__uncertainties:
            raise KeyError( 'Uncertainty "{}" does not affect any process in the collection.'.format( unc_name ) )
    
        #sum it across all affected processed 
        #note that since an uncertainty is correlated between processed the variations must be summed linearly rather than quadratically
        total_var = None
        for p in self.__processes:
            if not p.isSignal():
                continue
            if ( not p.hasUncertainty( unc_name ) ) :
                continue          
            if total_var is None:
                total_var = p.uncertainty( unc_name )
            else:
                total_var = addUncLinear( total_var, p.uncertainty( unc_name ) )

        return total_var

    def bkgUncertainty( self, unc_name ):
        if unc_name not in self.__uncertainties:
            raise KeyError( 'Uncertainty "{}" does not affect any process in the collection.'.format( unc_name ) )
    
        #sum it across all affected processed 
        #note that since an uncertainty is correlated between processed the variations must be summed linearly rather than quadratically
        total_var = None
        for p in self.__processes:
            if ( not p.hasUncertainty( unc_name ) ) or p.isSignal():
                continue

            if total_var is None:
                total_var = p.uncertainty( unc_name )
            else:
                total_var = addUncLinear( total_var, p.uncertainty( unc_name ) )

        return total_var


    def orderByYield( self ):
        def getYield( process ):
            return process.nominal().GetSumOfWeights()
        self.__processes.sort( key = getYield, reverse = True )


    #make the total uncertainty from one particular source symmetric for later summation of several uncertainty sources
    def bkgSymmetricEnvelope( self, unc_name ):
        try:
            return self.bkgUncertainty( unc_name ).envelope() 
        except AttributeError:
            return None

    #make the total uncertainty from one particular source symmetric for later summation of several uncertainty sources
    def signalSymmetricEnvelope( self, unc_name ):
        try:
            return self.signalUncertainty( unc_name ).envelope() 
        except AttributeError:
            return None    


    #make total systematic uncertainty band 
    def totalBkgSystUncertainty( self ):
        total_syst_hist = None
        for unc in self.__uncertainties:
            envelope = self.bkgSymmetricEnvelope( unc )
            if envelope is None:
                continue

            if total_syst_hist is None:
                total_syst_hist = envelope
            else:
                total_syst_hist = addHistQuadratic( total_syst_hist, envelope )
        return total_syst_hist
    
    def totalBkgWithStatErrors( self ):
        total_bkg_hist = None
        for p in self.__processes:
            if p.isSignal():
                continue
            if total_bkg_hist is None:
                total_bkg_hist = p.nominal().Clone()
            else:
                total_bkg_hist.Add( p.nominal() )
        return total_bkg_hist
    
    def totalBkgWithSystErrors( self ):
        total_bkg_hist = self.totalBkgWithStatErrors()
        total_bkg_syst = self.totalBkgSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            stat_error = 0.
            syst_error = total_bkg_syst.GetBinContent( b )
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist

    def totalBkgWithAllErrors( self ):
        total_bkg_hist = self.totalBkgWithStatErrors()
        total_bkg_syst = self.totalBkgSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            stat_error = total_bkg_hist.GetBinError( b )
            syst_error = total_bkg_syst.GetBinContent( b )
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist


    ########### signal
    def signal_get( self ):
       signal_hist = None
       for p in self.__processes:
            if p.isSignal():
                signal_hist = p.nominal().Clone()
       return signal_hist

    def signalSystUncertainty( self ):
        total_syst_hist = None
        for unc in self.__uncertainties:
            if unc == "dfShape" or unc == "dfLowStat" or unc == "dfmm" or unc == "dfem" or unc == "dfee": continue
            envelope = self.signalSymmetricEnvelope( unc )
            if envelope is None:
                continue

            if total_syst_hist is None:
                total_syst_hist = envelope
            else:
                total_syst_hist = addHistQuadratic( total_syst_hist, envelope )
        return total_syst_hist

    def signalWithStatErrors( self ):
        total_bkg_hist = None
        for p in self.__processes:
            if not p.isSignal():
                continue
            if total_bkg_hist is None:
                total_bkg_hist = p.nominal().Clone()
            else:
                total_bkg_hist.Add( p.nominal() )
        return total_bkg_hist
    
    def signalWithSystErrors( self ):
        total_bkg_hist = self.signalWithStatErrors()
        total_bkg_syst = self.signalSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            if total_bkg_hist.GetBinContent( b ) !=0:
                stat_error = 0.
                syst_error = total_bkg_syst.GetBinContent( b )
            else:
                stat_error = 0.
                syst_error = 0.
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist

    def signalWithAllErrors( self ):
        total_bkg_hist = self.signalWithStatErrors()
        total_bkg_syst = self.signalSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            stat_error = total_bkg_hist.GetBinError( b )
            syst_error = total_bkg_syst.GetBinContent( b )
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist


   

    #only nonprompt
    def npBkgSystUncertainty( self ):
        total_syst_hist = None
        for unc in self.__uncertainties:
            if unc != "dfShape" and unc != "dfLowStat" and unc != "dfmm" and unc != "dfem" and unc != "dfee": continue
            envelope = self.bkgSymmetricEnvelope( unc )
            if envelope is None:
                continue

            if total_syst_hist is None:
                total_syst_hist = envelope
            else:
                total_syst_hist = addHistQuadratic( total_syst_hist, envelope )
        return total_syst_hist

    def npBkgWithStatErrors( self ):
        total_bkg_hist = None
        for p in self.__processes:
            if p.isSignal() or p.name() == "Xgamma":
                continue
            if total_bkg_hist is None:
                total_bkg_hist = p.nominal().Clone()
            else:
                total_bkg_hist.Add( p.nominal() )
        return total_bkg_hist
    
    def npBkgWithSystErrors( self ):
        total_bkg_hist = self.npBkgWithStatErrors()
        total_bkg_syst = self.npBkgSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            stat_error = 0.
            syst_error = total_bkg_syst.GetBinContent( b )
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist

    def npBkgWithAllErrors( self ):
        total_bkg_hist = self.npBkgWithStatErrors()
        total_bkg_syst = self.npBkgSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            stat_error = total_bkg_hist.GetBinError( b )
            syst_error = total_bkg_syst.GetBinContent( b )
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist


     #only prompt
    def pBkgSystUncertainty( self ):
        total_syst_hist = None
        for unc in self.__uncertainties:
            if unc == "dfShape" or unc == "dfLowStat" or unc == "dfmm" or unc == "dfem" or unc == "dfee": continue
            envelope = self.bkgSymmetricEnvelope( unc )
            if envelope is None:
                continue

            if total_syst_hist is None:
                total_syst_hist = envelope
            else:
                total_syst_hist = addHistQuadratic( total_syst_hist, envelope )
        return total_syst_hist

    
    def pBkgWithStatErrors( self ):
        total_bkg_hist = None
        for p in self.__processes:
            if p.isSignal() or p.name() == "nonpromptSF" or p.name() == "nonpromptDF":
                continue
            if total_bkg_hist is None:
                total_bkg_hist = p.nominal().Clone()
            else:
                total_bkg_hist.Add( p.nominal() )
        return total_bkg_hist
    
    def pBkgWithSystErrors( self ):
        total_bkg_hist = self.pBkgWithStatErrors()
        total_bkg_syst = self.pBkgSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            if total_bkg_hist.GetBinContent( b ) !=0:
                stat_error = 0.
                syst_error = total_bkg_syst.GetBinContent( b )
            else:
                stat_error = 0.
                syst_error = 0.
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist

    def pBkgWithAllErrors( self ):
        total_bkg_hist = self.pBkgWithStatErrors()
        total_bkg_syst = self.pBkgSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            stat_error = total_bkg_hist.GetBinError( b )
            syst_error = total_bkg_syst.GetBinContent( b )
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist






    
   
    def __getitem__( self, process_name ):
        for p in self.__processes:
            if p.name() == process_name:
                return p
        raise KeyError( 'Process "{}" is not part of the collection.'.format( process_name ) )


    def __contains__( self, process_name ):
        try:
            self.__getitem__( process_name )
            return True
        except KeyError:
            return False


    def __iadd__( self, rhs ):
        for rhs_p in rhs.__processes:
            lhs_p = self.__getitem__( rhs_p.name() )
            lhs_p += rhs_p
        return self
            




if __name__ == '__main__':
    pass
