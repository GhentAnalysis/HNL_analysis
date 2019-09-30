import sys
import os


if __name__ == '__main__' :

    source_list = []
    with open( sys.argv[1] ) as makefile:
        for line in makefile.readlines():
            if line.startswith( 'SOURCES' ):
                print( line )
                for entry in line.split():
                    if entry.endswith( '.cc' ) or entry.endswith('.cpp' ) or entry.endswith('.cxx') or entry.endswith('.C') :
                        source_list.append( entry )

    command = 'clang++ --analyze'
    for source_file in source_list:
        command += ' {}'.format( source_file )
    command += ' `root-config --glibs --cflags`'
    os.system( command )
