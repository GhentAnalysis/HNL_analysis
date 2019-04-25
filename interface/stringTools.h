/*
Some simple tools which are useful in most physics analysis tasks
*/

#ifndef analysisTools_h
#define analysisTools_h

//include c++ library classes
#include<string>

//include ROOT classes
#include "TH1D.h"

namespace analysisTools{

    //show progress bar
    void printProgress(double);

    //Set negative histogram bins to zero
    void setNegativeBinsToZero(TH1D*);
    void setNegativeBinsToZero(const std::shared_ptr<TH1D>&);

    //print a combine datacard
    void printDataCard(const double obsYield, const double sigYield, const std::string& sigName, const double* bkgYield, const unsigned nBkg, const std::string* bkgNames, const std::vector<std::vector<double> >& systUnc, const unsigned nSyst, const std::string* systNames, const std::string* systDist, const std::string& cardName, const bool shapeCard = false, const std::string& shapeFileName = "");
}
#endif
