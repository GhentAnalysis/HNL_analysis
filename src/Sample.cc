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
    processName    fileName    xSec    (ctau[mm]   (new-V2))
    */
    std::string xSecString;     // temporary string to read xSection
    std::string ctauHnlString;  // temporary string to read lifetime ctau (if any)
    std::string v2HnlNewString; // temporary string to read new V^2 value (if any)
    //std::string signalString; // temporary string to fill signal boolean

    //first 3 words on the line are the process name, filename and cross section
    std::istringstream stream(line);
    stream >> process >> fileName >> xSecString >> ctauHnlString >> v2HnlNewString;

    // if no xSec is specified, set it to zero
    xSec = ( xSecString == "" ? 0. : std::stod(xSecString) );

    // if no ctau, it is irrelevant: set it to -1
    ctauHnl = ( ctauHnlString == "" ? -1. : std::stod(ctauHnlString) );

    // if no v2 new, it is irrelevant: set it to -1
    v2HnlNew = ( v2HnlNewString == "" ? -1. : std::stod(v2HnlNewString) );

    // Initialize HNL parameters
    couplHnl   = "";
    isDiracHnl = false;
    massHnl    = -1.;
    xSecNew    = -1.;
    ctauHnlNew = -1.;

    setHNL();
    setData();
    set2017();
    set2018();

   std::cout<<"is2017Sample: "<<is2017Sample<<std::endl;	    
    std::cout<<"is2016Sample: "<<is2016Sample<<std::endl;	    
   std::cout<<"is2018Sample: "<<is2018Sample<<std::endl;	    

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
    //setOptions(optionString);

    // Tmp
    std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << std::endl;
    std::cout << " File       : " << fileName         << std::endl;
    std::cout << " couplHnl   : " << getHNLcoupling() << std::endl;
    std::cout << " isDiracHnl : " << isHNLdirac()     << std::endl;
    std::cout << " massHnl    : " << getHNLmass()     << std::endl;
    std::cout << " v2Hnl      : " << getHNLV2()       << std::endl;
    std::cout << " v2HnlNew   : " << getHNLV2New()    << std::endl;
    std::cout << " ctauHnl    : " << getHNLctau()     << std::endl;
    std::cout << " ctauHnlNew : " << getHNLctauNew()  << std::endl;
    std::cout << " xSec       : " << getXSec()        << std::endl;
    std::cout << " xSecOrig   : " << getXSecOrig()    << std::endl;
    std::cout << " xSecNew    : " << getXSecNew()     << std::endl;
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

void Sample::setHNL(){
    if(fileName.find("HeavyNeutrino") != std::string::npos) {
        newPhysicsSignal = true;
	size_t pos = 26; // length of "HeavyNeutrino_trilepton_M-"
	std::string tmpstr = fileName.substr(pos);
	pos = tmpstr.find("_");
	massHnl = std::stod(tmpstr.substr(0, pos));
	tmpstr = tmpstr.substr(pos+3); // length of "_V-"
	pos = tmpstr.find("_");
        v2Hnl = std::stod(tmpstr.substr(0, pos));
	v2Hnl *= v2Hnl;
	tmpstr = tmpstr.substr(pos+1); // length of "_"
	pos = tmpstr.find("_");
	couplHnl = tmpstr.substr(0, pos);
	isDiracHnl = (tmpstr.find("_Dirac_" ) != std::string::npos);
	if(v2HnlNew>0.) {
	  xSecNew = xSec * (v2HnlNew/v2Hnl);
	  ctauHnlNew = ctauHnl * (v2Hnl/v2HnlNew);
	}
    }
}


void Sample::setData(){
    isDataSample = false;
    static std::vector<std::string> dataNames = {"data", "SingleMuon", "SingleElectron", "DoubleMuon", "DoubleEG", "EGamma"};
    for(auto it = dataNames.cbegin(); it != dataNames.cend(); ++it){
        if(fileName.find(*it) != std::string::npos){
            isDataSample = true;
        }
    }
}


void Sample::set2017(){
    is2017Sample = (fileName.find("Fall17"  ) != std::string::npos) || (fileName.find("2017") != std::string::npos);
}

void Sample::set2018(){
    is2018Sample = (fileName.find("Autumn18") != std::string::npos) || (fileName.find("2018") != std::string::npos);
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
    bool flag2016 = ( optionString.find("forceIs2016") != std::string::npos );
    bool flag2017 = ( optionString.find("forceIs2017") != std::string::npos );
    bool flag2018 = ( optionString.find("forceIs2018") != std::string::npos );
    unsigned ntrue = ((unsigned)flag2016) + ((unsigned)flag2017) + ((unsigned)flag2018);
    if(ntrue>1) {
      std::cerr << "Error in sample construction: forceIs2016=" << flag2016 << ", forceIs2017=" << flag2017 << ", forceIs2018=" << flag2018 << ", they cannot be all set to true!" << std::endl;
    }
    if(flag2018){
        is2017Sample = false;
        is2018Sample = true;
        uniqueName += "_forcedIs2018";
    } else if(flag2017){
        is2017Sample = true;
        is2018Sample = false;
        uniqueName += "_forcedIs2017";
    } else if(flag2016){
        is2017Sample = false;
        is2018Sample = false;
        uniqueName += "_forcedIs2016";
	    
    }
	std::cout<<"=======================.........  "<< is2018Sample<< "  "<< is2017Sample<<" . "<<flag2016<<std::endl;
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
        ( sam.is2018() ? "Autumn18" : (sam.is2017() ? "Fall17" : "Summer16") ) << 
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
