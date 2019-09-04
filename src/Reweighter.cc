#include "../interface/Reweighter.h"

//include c++ library classes
#include <map>
//include ROOT classes
#include "TFile.h"
#include "TROOT.h"

//include other parts of code 

Reweighter::Reweighter(const std::vector<Sample>& samples, const bool sampleIs2016) {
    initializePuWeights(samples);
}
void Reweighter::initializePuWeights(const std::vector< Sample >& sampleList){
std::cout<<"==============================================="<<std::endl;
    std::cout<<"==============================================="<<std::endl;
std::cout<<"==============================================="<<std::endl;
std::cout<<"==============================================="<<std::endl;
std::cout<<"==============================================="<<std::endl;

    static const std::string minBiasVariation[3] = {"central", "down", "up"};
    for(auto& sample : sampleList){

        //no pu weights for data 
        if( sample.isData() ) continue;

        //open root file corresponding to sample
        //TFile* puFile = TFile::Open( (const TString&) "PU/puWeights_" + sample.getFileName() );
        TFile* puFile = TFile::Open("/user/mvit/CMSSW_9_4_4/src/HNL_analysis/PU/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root");
        //extract pu weights 
        for(unsigned var = 0; var < 3; ++var){
            std::string histName = "puw_Run";
            if (sample.is2016()) histName += "2016" ;
            if (sample.is2017()) histName += "2017" ;
            if (sample.is2018()) histName += "2018" ;
            histName += "Inclusive_" + minBiasVariation[var];
            puWeights[sample.getUniqueName()].push_back( std::shared_ptr<TH1D> ( (TH1D*)puFile->Get( (const TString&) histName ) ) );

            //make sure histogram does not get deleted when closing file
            puWeights[sample.getUniqueName()].back()->SetDirectory(gROOT);
        }
        puFile->Close();
    }
    for( auto entry : puWeights ){
     std::cout << entry.first << std::endl;   
    }
}

Reweighter::~Reweighter(){}

double Reweighter::puWeight(const double nTrueInt, const Sample& sample, const unsigned unc) const{
    if(unc < 3){

        //find weights for given sample
        const auto& weightVectorIter = puWeights.find(sample.getUniqueName() );

        //check if pu weights are available for this sample
        if( weightVectorIter == puWeights.cend() ){
            std::cerr << "Error: no pu weights found for sample : " << sample.getUniqueName() << " returning weight 0  " << std::endl;
            return 0.;
        }

        //WARNING: 2016 samples might be used when making 2017 plots when no 2017 equivalent is available
        // In this case is2016() will return false, but the pileup weights for the sample will only exist
        // up to 50 interactions. So one needs to be sure to check for what campaign the sample was generated.
        // When using a sample (potentially simulated for 2017) for 2016 data analysis, always limit the bin
        // to 50 interactions which is the limit to which the data pileup profile was evaluated.
        bool sampleIsSimulatedFor2016 =  ( sample.getFileName().find("Summer16") != std::string::npos );
        double maxBin = ( (sampleIsSimulatedFor2016 || is2016) ?  49.5 : 99.5 );

        TH1D* weights = ( (*weightVectorIter).second)[unc].get();

        //buggy events in 2017 MC can have negative pileup! catch them 
        if( nTrueInt >= 0){
            return weights->GetBinContent(weights->FindBin(std::min(nTrueInt, maxBin) ) );
        } else {
            std::cerr << "Error: event with nTrueInt = " << nTrueInt << " , returning weight 9999." << std::endl;
            return 9999.;
        }
    }
    else {
        std::cerr << "Error: invalid pu uncertainty requested: returning weight 0" << std::endl;
        return 0.;
    }
}

