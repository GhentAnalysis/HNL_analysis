#ifndef Reweighter_H
#define Reweighter_H

//include c++ library classes

//include ROOT classes
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

//include other parts of code 
#include "../bTagSFCode/BTagCalibrationStandalone.h"
#include "../interface/Sample.h"

//Class storing scale-factor weights to be used in events
class Reweighter{
    public:
        Reweighter(const std::vector<Sample>&, const bool is2016, const std::string& bTagWorkingPoint = "medium");
        ~Reweighter();

        //pileup weight
        double puWeight(const double nTrueInt, const Sample&, const unsigned unc = 0) const;

    private:
        //boolean flagging weights as 2016 or 2017
        bool is2016;

        //pu weights (one for every sample)
        std::map< std::string, std::vector< std::shared_ptr<TH1D> > > puWeights;

        //initialize all weight histograms
        void initialize2016Weights();
        void initialize2017Weights();

        //return jet flavor index 1 -> 0, 4 -> 1, 5 -> 2
        unsigned flavorInd(const unsigned jetFlavor) const{ 
            return 0 + (jetFlavor == 4) + 2*(jetFlavor == 5);
        }

   
        //read pu weights for a given list of samples
        void initializePuWeights(const std::vector< Sample >&); 
    
        //initialize all weights 
        void initializeAllWeights(const std::vector< Sample>&);

   
};
#endif
