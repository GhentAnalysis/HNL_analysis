/*
Class used for storing the info related to a particular sample: xSec, process name, file name
*/
#ifndef Sample_H
#define Sample_H

//include c++ library classes
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

//include ROOT classes
#include "TFile.h"

class Sample{
    friend std::ostream& operator<<( std::ostream& os, const Sample& ); 
    public:
        Sample() = default;
        //Sample(const std::string& file, const std::string& proc, double cross = 0);
        Sample( const std::string& line, const std::string& directory ); 
        Sample( std::istream&, const std::string& directory ); 

        std::string getFileName() const { return fileName; }

        std::string getProcessName() const { return process; } 
    
        //to prevent overlapping file names when re-using a sample in both the 2016 and 2017 data lists 
        std::string getUniqueName() const { return uniqueName; }

        double getXSecOrig() const { return xSec; }

        double getXSecNew() const { return xSecNew; }

        double getXSec() const { return (v2HnlNew<0. ? xSec : xSecNew); }

        double getHNLmass() const { return massHnl; }

        double getHNLV2Orig() const { return v2Hnl; }

        double getHNLV2New() const { return v2HnlNew; }

        double getHNLV2() const { return (v2HnlNew<0. ? v2Hnl : v2HnlNew); }

        double getHNLctauOrig() const { return ctauHnl; }

        double getHNLctauNew() const { return ctauHnlNew; }

        double getHNLctau() const { return (ctauHnlNew<0. ? ctauHnl : ctauHnlNew); }

        std::string getHNLcoupling() const { return couplHnl; }

        bool isHNLdirac() const { return isDiracHnl; }

        bool isMajoranaToDiracSimul() const { return convMajToDir; }

        bool useLNCeventsOnly() const { return useLNConly; }

        bool isAMergedSample() const { return isMergedSample; }

        bool isData() const { return isDataSample; }

        bool isMC() const { return !isDataSample; }

        bool is2016() const { return !is2017Sample && !is2018Sample; }

        bool is2017() const { return  is2017Sample && !is2018Sample; }

        bool is2018() const { return !is2017Sample &&  is2018Sample; }

        bool isSMSignal() const { return smSignal; }

        bool isNewPhysicsSignal() const { return newPhysicsSignal; }

        std::shared_ptr<TFile> getFile() const;

    private:
        void setHNL(); 

        void setData(); 

        void set2017();
  
        void set2018();

        void setOptions(const std::string&);

        std::string fileName;
        std::string process;
        std::string uniqueName;
        std::string directory;

        double xSec;
        bool isDataSample;
        bool is2017Sample;
        bool is2018Sample;
        bool smSignal;
        bool newPhysicsSignal;
        // HNL parameters
        std::string couplHnl;
        bool isDiracHnl;
        bool convMajToDir;
        bool useLNConly;
        bool isMergedSample;
        double massHnl;
        double v2Hnl;
        double ctauHnl;
        // For re-weighting
        double xSecNew;
        double v2HnlNew;
        double ctauHnlNew;
};

//read a txt file containing a list of samples
std::vector< Sample > readSampleList( const std::string& list, const std::string& directory );
        
#endif
