#include "../interface/stringTools.h"

//Include c++ library classes
#include <algorithm>
#include <sstream>
#include <iomanip>


std::string stringTools::removeBackSpaces( const std::string& s ){
    std::string formated(s);
    while(formated.find_last_of("\t") == formated.size() - 1 || formated.find_last_of(" ") == formated.size() - 1){
        formated.erase(formated.size() - 1);
    }
    return formated;
}


std::string stringTools::removeFrontSpaces( const std::string& s ){
    std::string formated(s);
    while(formated.find("\t") == 0 || formated.find(" ") == 0){
        formated.erase(0, 1);
    }
    return formated;
}


std::string stringTools::cleanSpaces( const std::string& s ){
    std::string formated(s);

    //remove tabs and spaces in the front of the string
    formated = removeFrontSpaces(formated);

    //remove tabs and spaces in the end of the string
    formated = removeBackSpaces(formated);

    return formated;
}


std::string stringTools::extractFirst( std::string& s ){
    
    //clean all spaces 
    s = cleanSpaces(s);

    //find the location of the first tab or space ( the word ends here )
    auto space = std::min(s.find(" "), s.find("\t") );
    std::string substr = s.substr(0, space);
    s.erase(0, space);	
    s = removeFrontSpaces(s);
    return substr;
}


//add training / to directoryName if needed
std::string stringTools::directoryName( const std::string& directory ){
    std::string formated(directory);
    if(formated.back() != '/') formated.append("/");
    return formated;
}


bool stringTools::stringContains( const std::string& s, const std::string& substring ){
    return (s.find(substring) != std::string::npos);
}

bool stringTools::stringEndsWith( const std::string& s, const std::string& ending ){
    return std::equal(ending.crbegin(), ending.crend(), s.crbegin() );
}

std::pair< std::string, std::string > stringTools::splitFileExtension( const std::string& fileName ){
    size_t dotPosition = fileName.find_last_of(".");
    return { fileName.substr(0, dotPosition),  fileName.substr(dotPosition, fileName.size()) };
}

std::string stringTools::fileWithoutExtension( const std::string& fileName ){
    return stringTools::splitFileExtension(fileName).first;
}

std::string stringTools::doubleToString( const double input, const unsigned precision ){
    std::ostringstream stringStream;
    if( precision == 0 ){
        stringStream << input;
    } else{
        stringStream << std::setprecision( precision ) << input;
    }
    return stringStream.str();
}
