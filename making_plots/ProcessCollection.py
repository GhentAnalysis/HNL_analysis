from ROOT import Math
from ROOT import THStack
from ROOT import TROOT, TH1, TGraphAsymmErrors
import math

from Uncertainty import addUncLinear, addUncQuadratic, addHistQuadratic

class ProcessCollection:

    def __init__( self, process_list, sf = None, df = None ):
        
        self.__processes = process_list
        
        if sf and df:
            
            self.__sf = sf
            self.__df = df
            
            alpha = 1 - 0.6827            

            self.__sfCRGraph = TGraphAsymmErrors(self.__sf)
            self.__dfCRGraph = TGraphAsymmErrors(self.__df)
            
            for g in [self.__sfCRGraph, self.__dfCRGraph]:
                for i in range(g.GetN()):
                    N = g.GetY()[i]
                    LowErr = 0 if N == 0 else Math.gamma_quantile(alpha/2, N, 1.)
                    UpErr = Math.gamma_quantile_c(alpha/2, N+1, 1.) # recommended by StatComm (1.8)
                    # Math.gamma_quantile_c(alpha, N+1, 1.) # gives 1.2, strictly 68% one-sided
                    g.SetPointEYlow(i, N-LowErr)
                    g.SetPointEYhigh(i, UpErr-N)
            
            for p in self.__processes:
                
                if p.name() in ['nonpromptSF', 'nonpromptDF']:
                    
                    nom, low, high = [], [], []
                    
                    for ib in range(p.nominal().GetXaxis().GetNbins()):
                        if p.name() in ['nonpromptDF']: hnp = self.__dfCRGraph
                        else: hnp = self.__sfCRGraph
                        npcr = hnp.GetY()[ib]
                        npcrErrLow = hnp.GetEYlow()[ib]
                        npcrErrHigh = hnp.GetEYhigh()[ib]
                        np = p.nominal().GetBinContent(ib+1)
                        npErr = p.nominal().GetBinError(ib+1)
                        alpha = np/npcr if npcr > 0 else npErr
                        nom.append(npcr*alpha)
                        low.append(npcrErrLow*alpha)
                        high.append(npcrErrHigh*alpha)
                            
                    hnom = p.nominal().Clone()
                    for ib in range(hnom.GetXaxis().GetNbins()):
                        hnom.SetBinContent(ib+1, nom[ib])
                    npGraph = TGraphAsymmErrors(hnom)
                    for ib in range(hnom.GetXaxis().GetNbins()):
                        npGraph.SetPointEYlow(ib, low[ib])
                        npGraph.SetPointEYhigh(ib, high[ib])
                    if p.name() in ['nonpromptDF']:
                        self.__dfGraph = npGraph
                    else:
                        self.__sfGraph = npGraph
                        
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
        low, high = None, None
        for p in self.__processes:
            if p.isSignal():
                continue
            nbins = p.nominal().GetXaxis().GetNbins()
            if low is None: low = [0]*nbins
            if high is None: high = [0]*nbins
            
            if total_bkg_hist is None:
                total_bkg_hist = p.nominal().Clone()
            else:
                total_bkg_hist.Add( p.nominal() )
            gr = None
            if p.name() in ['nonpromptDF']: gr = self.__dfGraph
            elif p.name() in ['nonpromptSF']: gr = self.__sfGraph
            if gr:
                for ib in range(nbins):
                    low[ib] += math.pow(gr.GetEYlow()[ib], 2)
                    high[ib] += math.pow(gr.GetEYhigh()[ib], 2)
            else:
                for ib in range(nbins):
                    low[ib] += math.pow(p.nominal().GetBinError(ib), 2)
                    high[ib] += math.pow(p.nominal().GetBinError(ib), 2)
                    
        for i in range(len(low)):
            low[ib] = math.sqrt(low[ib])
            high[ib] = math.sqrt(high[ib])
                    
        total_bkg_graph = TGraphAsymmErrors(total_bkg_hist)
        for ib in range(total_bkg_graph.GetN()):
            total_bkg_graph.SetPointEYlow(ib, low[ib])
            total_bkg_graph.SetPointEYhigh(ib, high[ib])
                    
        return [total_bkg_hist, total_bkg_graph]
    
    def totalBkgWithSystErrors( self ):
        total_bkg_hist, total_bkg_graph = self.totalBkgWithStatErrors()
        total_bkg_syst = self.totalBkgSystUncertainty()
        for b in range( 1, total_bkg_hist.GetNbinsX() + 1 ):
            stat_error = 0.
            syst_error = total_bkg_syst.GetBinContent( b )
            total_bkg_hist.SetBinError( b, ( stat_error*stat_error + syst_error*syst_error )**0.5 )
        return total_bkg_hist

    def totalBkgWithAllErrors( self ):
        total_bkg_hist, total_bkg_graph = self.totalBkgWithStatErrors()
        total_bkg_syst = self.totalBkgSystUncertainty()
        for b in range( total_bkg_hist.GetNbinsX() ):
            stat_errorLow = total_bkg_graph.GetEYlow()[b]
            stat_errorHigh = total_bkg_graph.GetEYhigh()[b]
            syst_error = total_bkg_syst.GetBinContent( b+1 )
            total_bkg_graph.SetPointEYlow( b, ( stat_errorLow*stat_errorLow + syst_error*syst_error )**0.5 )
            total_bkg_graph.SetPointEYhigh( b, ( stat_errorHigh*stat_errorHigh + syst_error*syst_error )**0.5 )
        return total_bkg_graph


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
