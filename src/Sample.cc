#include "../interface/Sample.h"

//include c++ library classes 
#include <sstream>
#include <fstream>

//include other parts of code 
#include "../interface/stringTools.h"


Sample::Sample( const std::string& line, const std::string& sampleDirectory ) :
    directory( sampleDirectory)
{
    /*
    only works if input line is formatted as:
    processName    fileName    xSec
    */
    std::string xSecString;   //temporary string to read xSection
    std::string signalString; //temporary string to fill signal boolean

    //first 3 words on the line are the process name, filename and cross section
    std::istringstream stream(line);
    stream >> process >> fileName >> xSecString;

    //if not Xsection is specified it is zero
    xSec = ( xSecString == "" ? 0 : std::stod(xSecString) );

    setData();
    set2017();

    //unique name is equal to fileName without file extension
    uniqueName = stringTools::fileWithoutExtension( fileName );

    //data has no xSection
    if(isData() && xSecString != ""){
        std::cerr << "xSection specified for data: are you sure this was intended?" << std::endl;
    }

    //extract all optional strings at the end of the line
    std::string optionString;
    std::string tempString;
    while(stream){
        stream >> tempString;
        optionString.append(tempString);
    }

    //read options
    //This might modify uniqueName. uniqueName has to be set before calling this function!
    setOptions(optionString);
}


Sample::Sample( std::istream& is, const std::string& directory ){

    //read sample info from txt file
    std::string line;

    //jump to next line if current line is a comment
    bool nextLineIsComment;
    do{
        nextLineIsComment = false;
        if(std::getline(is, line)){
            nextLineIsComment =  (line[line.find_first_not_of(" \t")] == '#');
            if(!nextLineIsComment){
                *this = Sample(line, directory); 
            }
        }
    } while(nextLineIsComment);
}


void Sample::setData(){
    isDataSample = false;
    static std::vector<std::string> dataNames = {"data", "SingleMuon", "SingleElectron", "SingleMuon", "DoubleMuon", "DoubleEG"};
    for(auto it = dataNames.cbegin(); it != dataNames.cend(); ++it){
        if(fileName.find(*it) != std::string::npos){
            isDataSample = true;
        }
    }
}


void Sample::set2017(){
    is2017Sample = (fileName.find("Fall17") != std::string::npos) || (fileName.find("2017") != std::string::npos);
}


void Sample::setOptions( const std::string& optionString ){
    if(optionString == ""){
        smSignal = false;
        newPhysicsSignal = false;
        return;
    } 

    //signal flags
    //determine whether process is some kind of signal
    smSignal = ( optionString.find("SMSignal") != std::string::npos );
    newPhysicsSignal = ( optionString.find("newPhysicsSignal") != std::string::npos );

    //signal can not be both SM and BSM sigal
    if(smSignal && newPhysicsSignal){
        std::cerr << "Error in sample construction: sample is both SM and BSM signal" << std::endl;
    }
    
    //check if sample needs to be used in different era it was intended for (i.e. 2016 sample when comparing to 2017 data)
    bool flag2017 = ( optionString.find("forceIs2017") != std::string::npos );
    bool flag2016 = ( optionString.find("forceIs2016") != std::string::npos );
    if(flag2016 && flag2017){
        std::cerr << "Error in sample construction: both forceIs2016 and forceIs2017 flags were set, can not set both!" << std::endl;
    }
    if(flag2017){
        is2017Sample = true;
        uniqueName += "_forcedIs2017";
    } else if(flag2016){
        is2017Sample = false;
        uniqueName += "_forcedIs2016";
    }
}


std::shared_ptr<TFile> Sample::getFile() const{
    return std::make_shared<TFile>( ( stringTools::directoryName(directory) + fileName).c_str() , "read");
}


//print Sample info
std::ostream& operator<<( std::ostream& os, const Sample& sam ){
    os << sam.process << "\t" << 
        sam.fileName << "\t" << 
        sam.xSec << "\t" << 
        ( sam.isData() ? "data" : "MC" ) << "\t" << 
        ( sam.is2017() ? "Fall17" : "Summer16" ) << 
        ( sam.smSignal ? "\tSM signal" : "" ) << 
        ( sam.newPhysicsSignal ? "\tBSM signal" : "" );
    return os;
}

//read a list of samples into a vector 
std::vector< Sample > readSampleList( const std::string& listFile, const std::string& directory ){
	std::vector< Sample> sampleList;

    //read sample info from txt file
    std::ifstream inFile(listFile);
   

    while( !inFile.eof() ){

        sampleList.push_back( Sample( inFile, directory ) );

    }
    sampleList.pop_back();

    //close file after usage
    inFile.close();

    return sampleList;
}
