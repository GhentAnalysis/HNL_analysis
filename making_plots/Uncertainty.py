
def checkBinConsistency( lhs_hist, rhs_hist ):
    if lhs_hist.GetNbinsX() != rhs_hist.GetNbinsX():
        raise IndexError( 'left histograms has {} bins, while right histogram has {} bins.'.format( lhs_hist.GetNbinsX(), rhs_hist.GetNbinsX() ) )



class Uncertainty:

    def __init__( self, varDown, varUp ):
        checkBinConsistency( varDown, varUp )
            
        self.__varDown = varDown
        self.__varUp = varUp


    def varDown( self ):
        return self.__varDown


    def varUp( self ):
        return self.__varUp


    def envelope( self ):
        envelope_hist = self.__varDown.Clone()
        for b in range( 1, envelope_hist.GetNbinsX() + 1 ):
            max_var = max( abs( self.varDown().GetBinContent( b ) ), abs( self.varUp().GetBinContent( b ) ) )
            envelope_hist.SetBinContent( b, max_var )
        return envelope_hist


    def scale( self, scale_factor ):
        self.__varDown.Scale( scale_factor )
        self.__varUp.Scale( scale_factor )



def addUncLinear( lhs, rhs ):
    summedHistDown = lhs.varDown().Clone()
    summedHistDown.Add( rhs.varDown() )

    summedHistUp = lhs.varUp().Clone()
    summedHistUp.Add( rhs.varUp() )

    return Uncertainty( summedHistDown, summedHistUp )


def addHistQuadratic( lhs, rhs ):
    checkBinConsistency( lhs, rhs )    

    summed_hist = lhs.Clone()
    for b in range( 1, summed_hist.GetNbinsX() + 1 ):
        lhs_content = lhs.GetBinContent( b )
        rhs_content = rhs.GetBinContent( b )
        
        summed_hist.SetBinContent( b, ( lhs_content*lhs_content + rhs_content*rhs_content )**0.5 )

    return summed_hist
        
        
def addUncQuadratic( lhs, rhs ):
    summedHistDown = addHistQuadratic( lhs.varDown(), rhs.varDown() )
    summedHistUp = addHistQuadratic( lhs.varUp(), rhs.varUp() )
    return Uncertainty( summedHistDown, summedHistUp )


if __name__ == '__main__':
    a = Uncertainty( "up", "down" )
