#include "Bug.h"

//include c++ library classes
#include <fstream>
#include <sstream>

Bug::Bug(){
    std::ifstream textStream( "Bug.txt" );
    std::string line;
    while( std::getline( textStream, line ) ){
        lines.push_back( line );
    }
    textStream.close();
}


void Bug::printBug( std::ostream& os ) const{
    for( const auto& line : lines ){
        os << line << "\n";
    }
    os << std::flush;
}
