import os


def readProcessNamesAndIndices( datacard_lines ):
    processes = []
    process_indices = []
    counter = 0
    for line in datacard_lines:

        #datacards typically have two lines starting with process, the first one giving the names of each process
        if line.startswith( 'process' ) and counter == 0:
            processes = line.split()[1:]
            counter += 1
        elif line.startswith( 'process' ):
            process_indices = [ int(i) for i in line.split()[1:] ]

    if ( not processes ) or ( not process_indices ):
        raise LookupError( 'Either no list of processes or process indices was found in the datacard.' )

    if ( len( processes ) != len( process_indices ) ):
        raise LookupError( 'There are {} processes, while there are {} indices'.format( len( processes ), len( process_indices ) ) )

    processes_with_indices = {}
    for p, i in zip( processes, process_indices ):
        processes_with_indices[ i ] = p

    return processes_with_indices



def lineIsNuisance( line ):
    return ( 'lnN' in line or ( 'shape' in line and not 'shapes' in line ) )


def lineIsNormNuisance( line ):
    return ( lineIsNuisance( line ) and ( not 'shape' in line ) )


def lineIsShapeNuisance( line ):
    return ( lineIsNuisance( line ) and ( 'shape' in line ) )


def readNuisanceLine( line ):
    entries = line.split()
    name = entries[ 0 ]
    nuisance_size_str = entries[ 2: ]

    nuisance_sizes = {}
    for i, nuisance in enumerate( nuisance_size_str ):
        try:
            nuisance_sizes[ i ] = float( nuisance )
        except ValueError:
            pass
    return name, nuisance_sizes



def indexToProcessNameNuisances( processIndex_to_nuisance, processIndex_to_name ):
    return { processIndex_to_name[ key ] : processIndex_to_nuisance[ key ] for key in processIndex_to_nuisance }


def buildProcessToUncDict( nuisance_lines, processIndex_to_name ):
    nuisance_dict = {}
    for l in nuisance_lines:
        nuisance_name, processIndex_to_nuisance = readNuisanceLine( l )
        nuisance_dict[ nuisance_name ] = indexToProcessNameNuisances( processIndex_to_nuisance, processIndex_to_name )

    return nuisance_dict


#functions to retrieve the shapes to which the datacard points 
def lineIsShapePath( line ):
    return ( '$PROCESS_$SYSTEMATIC' in line )


def extractShapePathLine( datacard_lines ):
    for l in datacard_lines:
        if lineIsShapePath( l ):
            return l
    raise KeyError( 'No line specifying the shape is found in the datacard.' )



def extractShapeFilePath( datacard_lines ):
    shape_path_line = extractShapePathLine( datacard_lines )
    shape_path = None
    for part in shape_path_line.split():
        if '.root' in part:
            return part
    raise KeyError( 'No .root file is specified in line "{}".'.format( shape_path_line ) )


class Datacard:

    def __init__( self, datacard_path ):
        
        lines = None
        with open( datacard_path ) as f:
            lines = f.readlines()
        
        processIndex_to_name = readProcessNamesAndIndices( lines )
        self.__process_names = [ p for _, p in processIndex_to_name.items() ]
        
        normalization_uncertainty_lines = [ l for l in lines if lineIsNormNuisance( l ) ]
        self.__normalization_uncertainties = buildProcessToUncDict( normalization_uncertainty_lines, processIndex_to_name )
        
        shape_uncertainty_lines = [ l for l in lines if lineIsShapeNuisance( l ) ]
        self.__shape_uncertainties = buildProcessToUncDict(  shape_uncertainty_lines, processIndex_to_name )

        relative_shape_path = extractShapeFilePath( lines )
        datacard_directory = os.path.dirname(  datacard_path )
        self.__shape_file_path = os.path.join( datacard_directory, relative_shape_path )


    def processNames( self ):
        return self.__process_names


    def __findUncertaintySize( self, uncertainty_name, process_name, uncertainty_dictionary ):
        try:
            uncertainty_sizes = uncertainty_dictionary[ uncertainty_name ]
            try:
                return uncertainty_sizes[ process_name ]
            except KeyError:
                return None
                #raise KeyError( 'Uncertainty source "{}" is not available for process "{}"'.format( uncertainty_name, process_name ) )
            
        except KeyError:
            return None
            #print( uncertainty_dictionary.keys() )
            #raise KeyError( 'Uncertainty source "{}" is not present in datacard.'.format( uncertainty_name ) )


    def normalizationUncertainties( self ):
        return [ key for key in self.__normalization_uncertainties ]


    def shapeUncertainties( self ):
        return [ key for key in self.__shape_uncertainties ]
        

    def normalizationUncertaintySize( self, uncertainty_name, process_name ):
        return self.__findUncertaintySize( uncertainty_name, process_name, self.__normalization_uncertainties )


    def shapeUncertaintySize( self, uncertainty_name, process_name ):
        return self.__findUncertaintySize( uncertainty_name, process_name, self.__shape_uncertainties )

    
    def shapeFilePath( self ):
        return self.__shape_file_path




if __name__ == '__main__':
    pass
