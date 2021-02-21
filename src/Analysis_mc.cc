#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <Riostream.h>
#include <cmath>
#include <cstring>
#include <string>
#include <tuple>
#include <set>

#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLorentzRotation.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TStyle.h"
#include "TLine.h"
#include "TProfile.h"
#include "TAttFill.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TLegendEntry.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"
#include "TPaveText.h"
#include "../interface/Analysis_mc.h"
#include "TApplication.h"
#include "TColor.h"


//include C++ library classes
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

//include other parts of the code
#include "../interface/tdrstyle.h"
#include "../interface/plotCode_new.h"
#include "../interface/kinematicTools.h"
#include "../interface/Reweighter.h"

// For b-tagging SFs and variations thereof
#include "../interface/BTagCalibrationStandalone.h"
//#include "BTagCalibrationStandalone.cpp"



using namespace std;

//ClassImp(Analysis_mc)

//_______________________________________________________default constructor_____
Analysis_mc::Analysis_mc():TObject()
{
}

//_______________________________________________________ constructor_____
Analysis_mc::Analysis_mc(unsigned jaar, const std::string& list, const std::string& directory) : TObject() {
  if(jaar>2) {
    std::cout << " --- WARNING: invalid value for 'year' variable (" << jaar
	      << "), setting it to 0 (i.e. 2016) ---" << std::endl;
    year = 0;
  }
  else {
    std::cout << " >>> Applying selection for "
	      << ( jaar==0 ? "2016" : ( jaar==1 ? "2017" : "2018" ) )
	      << " data set." << std::endl;
    year = jaar;
  }


  // >>>>>>>>>>>>>>>>>>> samples info >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  samples = readSampleList(list, directory);
  //std::cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> " << samples.size() << std::endl;
  nSamples_eff       = 0;
  nSamples_signal    = 0;
  nSamples_signal_e  = 0;
  nSamples_signal_mu = 0;
  TString *bpro = std::begin(eff_names);
  TString *epro = std::end(eff_names);

  for(size_t is=0; is<samples.size(); ++is) {
    TString procname(samples[is].getProcessName().c_str());
    TString *ipro = std::find(bpro, epro, procname);
    if(ipro!=epro) continue;
    if(nSamples_eff==max_nSamples_eff+1)
      throw std::runtime_error("nSamples_eff == max_nSamples_eff + 1");
    eff_names[nSamples_eff] = procname;
    if(samples[is].isNewPhysicsSignal()) {
      if(nSamples_signal==max_nSamples_signal) throw std::runtime_error("nSamples_signal == max_nSamples_signal");
      sigNames[nSamples_signal] = eff_names[nSamples_eff];
      sigNames_short[nSamples_signal] = "M = ";
      sigNames_short[nSamples_signal] += samples[is].getHNLmass();
      sigNames_short[nSamples_signal] += "GeV, |V_{NYY}|^{2} = ";
      sigNames_short[nSamples_signal] += (samples[is].getHNLV2New()>0. ? samples[is].getHNLV2New() : samples[is].getHNLV2());
      //if(eff_names[nSamples_eff].EndsWith("_e")) {
      if(eff_names[nSamples_eff].Contains("_e")) {
	if(nSamples_signal_e==max_nSamples_signal_e) throw std::runtime_error("nSamples_signal_e == max_nSamples_signal_e");
	sigNames_short[nSamples_signal].ReplaceAll("YY", "e");
	sigNames_e[nSamples_signal_e] = eff_names[nSamples_eff];
	string_sigNames_e[nSamples_signal_e] = samples[is].getProcessName();
	++nSamples_signal_e;
      }
      //else if(eff_names[nSamples_eff].EndsWith("_mu")) {
      else if(eff_names[nSamples_eff].Contains("_mu")) {
	if(nSamples_signal_mu==max_nSamples_signal_mu) throw std::runtime_error("nSamples_signal_mu == max_nSamples_signal_mu");
	sigNames_short[nSamples_signal].ReplaceAll("YY", "#mu");
	sigNames_mu[nSamples_signal_mu] = eff_names[nSamples_eff];
	string_sigNames_mu[nSamples_signal_mu] = samples[is].getProcessName();
	++nSamples_signal_mu;
      }
      else {
	throw std::runtime_error("Signal sample name: "+eff_names[nSamples_eff]);
      }
      ++nSamples_signal;
    } // end if(samples[is].isNewPhysicsSignal())
    ++nSamples_eff;
  } // end for(size_t is=0; is<samples.size(); ++is)

  // Add single- and double-fake backgrounds
  eff_names[nSamples_eff] = "nonprompt SF neg";
  ++nSamples_eff;
  eff_names[nSamples_eff] = "nonprompt SF pos";
  ++nSamples_eff;
  eff_names[nSamples_eff] = "nonprompt DF";
  ++nSamples_eff;

  // Remove "Obs" from count
  --nSamples_eff;

  std::cout << " >>>> nSamples_eff = " << nSamples_eff << " - nSamples_signal = " << nSamples_signal
	    << " - nSamples_signal_e = " << nSamples_signal_e << " - nSamples_signal_mu = " << nSamples_signal_mu
	    << std::endl;
  for(size_t ii=0; ii<nSamples_eff; ++ii) std::cout << eff_names[ii] << std::endl;

  if(nSamples_signal!=(nSamples_signal_e+nSamples_signal_mu))
    throw std::runtime_error("nSamples_signal != nSamples_signal_e + max_nSamples_signal_mu");


  const int nBin_probvertex = 18;	
  Double_t bin_probvertex[nBin_probvertex+1] = {0,0.0005,0.001,0.005,0.01,0.03,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.75,1.};
	
  const int nBin_cos = 28;	
  Double_t bin_cos[nBin_cos+1] = {-1.,-0.99,-0.98,-0.9,-0.8,-0.7, -0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9, 0.92,0.94,0.96,0.97,0.98,0.99,1.};

  const int nBin_2d = 13;	
  Double_t bin_2d[nBin_2d+1] = {0,1,2,3,4,5,6,7,8,9,10,15,25,50};

  const int nBin_2d_zoom = 8;	
  Double_t bin_2d_zoom[nBin_2d_zoom+1] = {0.,0.25,0.5,0.75,1.,1.50,2.,3.,4.};

  const int nBin_mass = 7;	
  Double_t bin_mass[nBin_mass+1] = {0.,1.,2.,3.,4.,5.,8.,12.};

  
	
  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> histogramms creation >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  for(int i = 0; i < nDist; ++i){
    float BinWidth = (HistMax[i] - HistMin[i])/nBins[i];
    std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
    for(size_t effsam = 0; effsam < nSamples_eff + 1; ++effsam){
      for(int cat = 0; cat < nCat; ++cat){
	//if (cat !=0 && cat !=6) continue;
	for(int cha = 0; cha < nChannel; ++cha){  
	  if (i != 34 && i!= 37 && i!=40 && i!=41){
	    Histos[i][cha][cat][effsam] =  new TH1D(eff_names[effsam] +"_"+ channelNames[cha] +"_"+ catNames[cat] +"_"+ Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + ";events/" + Yaxis + Units[i], nBins[i], HistMin[i], HistMax[i]);  
	  }	
	  if (i == 34) Histos[i][cha][cat][effsam] =  new TH1D(eff_names[effsam] +"_"+ channelNames[cha] +"_"+ catNames[cat] +"_"+ Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + ";events/" + Yaxis + Units[i], nBin_probvertex, bin_probvertex);  
	  if (i == 37) Histos[i][cha][cat][effsam] =  new TH1D(eff_names[effsam] +"_"+ channelNames[cha] +"_"+ catNames[cat] +"_"+ Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + ";events/" + Yaxis + Units[i], nBin_cos, bin_cos);  
	  if (i == 40) Histos[i][cha][cat][effsam] =  new TH1D(eff_names[effsam] +"_"+ channelNames[cha] +"_"+ catNames[cat] +"_"+ Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + ";events/" + Yaxis + Units[i], nBin_2d, bin_2d);
	  if (i == 41) Histos[i][cha][cat][effsam] =  new TH1D(eff_names[effsam] +"_"+ channelNames[cha] +"_"+ catNames[cat] +"_"+ Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + ";events/" + Yaxis + Units[i], nBin_2d_zoom, bin_2d_zoom);  			  
	  Histos[i][cha][cat][effsam]->Sumw2();	
	}
      }
    }
  }
  //Calculate the center of the maximum bin of each histogram
  for(int i = 0; i < nDist; ++i){
    maxBinC[i] = Histos[i][0][0][0]->GetBinCenter(Histos[i][0][0][0]->GetNbinsX());
  }

  /*for(int i = 0; i < nhist2d; ++i){
    for(size_t effsam = 0; effsam < number_2d; ++effsam){
      for(int cat = 0; cat < nCat2D; ++cat){
	for(int flav = 0; flav < nFlavor; ++flav){

	  Histos2d[i][cat][flav][effsam] =  new TH2D(names2d[effsam]+"_"+ catNames2D[cat] +"_"+ flavorNames[flav]+ "_"+ Hist2d_ossf[i] ,names2d[effsam]+"_"+ catNames2D[cat] +"_"+ flavorNames[flav] +"_"+ Hist2d_ossf[i], XnBins2d[i], XHistMin2d[i], XHistMax2d[i], YnBins2d[i], YHistMin2d[i], YHistMax2d[i]);
	
	 
	  std::cout<<"dist: "<<i<<"  2d  "<<effsam<<"  "<<cat<<"  "<<flav<<"    --> "<< names2d[effsam]+"_"+ catNames2D[cat] +"_"+ flavorNames[flav]+ "_"+ Hist2d_ossf[i] <<std::endl;
	  Histos2d[i][cat][flav][effsam] ->Sumw2();
	}
      }
    }
  }
  */


  
  for(size_t effsam = 0; effsam < nSamples_eff + 1; ++effsam){
    std::cout<<"-----------------> pippo: "<< eff_names[effsam] <<std::endl;
  }	
	
 
  // plot for limits 
  // weights for limits	
  for(size_t effsam = 0; effsam < nSamples_eff + 1; ++effsam){
    for(int var = 0; var < nVariation; ++var){
      for (int syst = 0; syst < nSystematic; ++syst)	{
	for(int cha = 0; cha < nCoupling; ++cha){
	  
	  for (int nplo = 0; nplo < nPlots; ++nplo){
	    if (nplo == 0) plots_SR2[nplo][cha][syst][var][effsam] = new TH1D(std::to_string(jaar)+"_"+plots_names[nplo]+"_"+eff_names[effsam]+"_"+ chaNames[cha] +"_"+systNamesT[syst]+"_"+varNames[var], eff_names[effsam]+"_"+ chaNames[cha] +"_"+systNamesT[syst]+"_"+varNames[var], nBin_mass, bin_mass);
	    if (nplo == 1) plots_SR2[nplo][cha][syst][var][effsam] = new TH1D(std::to_string(jaar)+"_"+plots_names[nplo]+"_"+eff_names[effsam]+"_"+ chaNames[cha] +"_"+systNamesT[syst]+"_"+varNames[var], eff_names[effsam]+"_"+ chaNames[cha] +"_"+systNamesT[syst]+"_"+varNames[var], nBin_2d_zoom, bin_2d_zoom);	    
	    plots_SR2[nplo][cha][syst][var][effsam]-> Sumw2();
	  }
	  
	  plots_SR[cha][syst][var][effsam] = new TH1D(std::to_string(jaar)+"_"+eff_names[effsam]+"_"+ chaNames[cha] +"_"+systNamesT[syst]+"_"+varNames[var], eff_names[effsam]+"_"+ chaNames[cha] +"_"+systNamesT[syst]+"_"+varNames[var], 18, 0.5, 18.5);
	  plots_SR[cha][syst][var][effsam]-> Sumw2();
	  weight_SR[cha][syst][var][effsam]=1.;
	}
      }
    }
  }
  for(int i = 0; i < nPlots; ++i){
    maxBinC_dc[i] = plots_SR2[i][0][0][0][0]->GetBinCenter(plots_SR2[i][0][0][0][0]->GetNbinsX());
  }


  for (int ibin =0; ibin<18; ++ibin){
    for(size_t effsam = 0; effsam < nSamples_eff + 1; ++effsam){
      for(int cha = 0; cha < nCoupling; ++cha){
	syst_error[ibin][cha][effsam] = 0.;
      }
    }
  }
  


  
}
//________________________________________________________________distruttore_____
Analysis_mc::~Analysis_mc()	 {
}
//          ================= ================= ================= ================= ================= =================          // 
//          ================= ================= ================= ================= ================= =================          // 
//_______________________________________________________ print status _____
void Analysis_mc::printProgress(double progress){
  const unsigned barWidth = 100;
  std::cout << "[";
  unsigned pos = barWidth * progress;
  for (unsigned i = 0; i < barWidth; ++i) {
    if (i < pos) std::cout << "=";
    else if (i == pos) std::cout << ">";
    else std::cout << " ";
  }
  std::cout << "] " << unsigned(progress * 100.0) << " %\r" << std::flush;
}
//_______________________________________________________ read sample _____
void Analysis_mc::readSamples(const std::string& list, const std::string& directory, std::vector<Sample>& sampleVector){

  //clean current sample list 
  sampleVector.clear();
  //read list of samples from file
  sampleVector = readSampleList(list, directory);
  //print sample information
  for(auto& sample : sampleVector){
    std::cout << sample << std::endl;
  }
}
//_______________________________________________________ read sample _____
void Analysis_mc::readSamples(const std::string& list, const std::string& directory){
  readSamples(list, directory, this->samples);
}
//_______________________________________________________ initialize sample _____
void Analysis_mc::initSample(const Sample& samp){ 
  std::cout << "In initSample(sample)"<<std::endl;
  //update current sample
  currentSample = samp;
  sampleFile = samp.getFile();
  sampleFile->cd("blackJackAndHookers");
  fChain = (TTree*) sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree");
  initTree(fChain, samp.isData(), samp.isNewPhysicsSignal());
  nEntries = fChain->GetEntries();
  if(!samp.isData()){

    //read sum of simulated event weights
    TH1D* hCounter = new TH1D("hCounter", "Events counter", 1, 0, 1);
    hCounter->Read("hCounter");
    sumSimulatedEventWeights = hCounter->GetBinContent(1);
    delete hCounter;

    //event weights set with lumi depending on sample's era 
    double dataLumi;
    if( year == 0 ){
      dataLumi = lumi2016;
    }
    else if ( year == 1 ){
      dataLumi = lumi2017;
    }
    else dataLumi = lumi2018;
 
    std::cout<<"''''''''''''''''''''' lumi ''''''''''"<<dataLumi<<std::endl;	  
    // N.B.: getXSec() returns the cross section, or the *re-weighted* cross section in case of V^2 (or ctau) re-weighting
    scale = samp.getXSec()*dataLumi*1000/sumSimulatedEventWeights;       //xSec*lumi divided by total sum of simulated event weights
    //
    // Add a factor 2 if it's a Majorana sample used to simulate a Dirac one, and
    // we are only using LNC events (i.e. half the statistics in the Majorana sample)
    //if(samp.isNewPhysicsSignal() && samp.isMajoranaToDiracSimul() && samp.useLNCeventsOnly()) scale *= 2;
    if(samp.isNewPhysicsSignal() && samp.useLNCeventsOnly()) scale *= 2;
  }
}
//_______________________________________________________ initialize weight for PU ____
void Analysis_mc::initializeWeights(){
  static bool weightsAre2016 = is2016();
  bool firstTime = ( reweighter.use_count() == 0 );
  bool changedEra = ( weightsAre2016 != is2016() );
  if( firstTime || changedEra){
    weightsAre2016 = is2016();
    //automatically use b-tag reshaping for nows
    reweighter.reset(new Reweighter(samples /*, is2016()*/) );
  } 
}
//_______________________________________________________  weight for PU ____
double Analysis_mc::PUWeight(){
  //check if weights are initialized, and initialize if needed 
  initializeWeights();
  //pileup reweighting
  double sf = puWeight();
  if (sf == 0) 	std::cout<<_nTrueInt<<std::endl;
  if( _nTrueInt < 0){
    std::cerr << "Error: event with negative pileup, returning SF weight 0." << std::endl;
    return 0.;
  }
  return sf;	
}
double Analysis_mc::puWeight(const unsigned unc) const{
  return reweighter->puWeight(_nTrueInt, currentSample, unc);
}

//_______________________________________________________ initialize sample ____
void Analysis_mc::initSample(){ //initialize the next sample in the list 
  initSample(samples[++currentSampleIndex]);
}
//_______________________________________________________ initialize sample ____
void Analysis_mc::GetEntry(const Sample& samp, long unsigned entry)
{
  if (!fChain) return;
  fChain->GetEntry(entry);
  //Set up correct weights
  if(!samp.isData() ) weight = _weight*scale; //MC
  else weight = 1;                            //data
}
//_______________________________________________________ initialize sample ____
void Analysis_mc::GetEntry(long unsigned entry){    //currently initialized sample when running serial
  GetEntry(samples[currentSampleIndex], entry);
}
//_______________________________________________________ initialize tree ____
void Analysis_mc::initTree(TTree *tree, const bool isData, const bool isNewPhys)
{
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("_runNb", &_runNb, &b__runNb);
  fChain->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
  fChain->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
  fChain->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);   
  fChain->SetBranchAddress("_passMETFilters", &_passMETFilters, &b__passMETFilters); 
  fChain->SetBranchAddress("_passTrigger_1l", &_passTrigger_1l, &b__passTrigger_1l);   
  fChain->SetBranchAddress("_HLT_IsoMu24", &_HLT_IsoMu24, &b__HLT_IsoMu24);
  if (year==0) fChain->SetBranchAddress("_HLT_IsoTkMu24", &_HLT_IsoTkMu24, &b__HLT_IsoTkMu24);
  if (year==0) fChain->SetBranchAddress("_HLT_Ele27_WPTight_Gsf", &_HLT_Ele27_WPTight_Gsf, &b__HLT_Ele27_WPTight_Gsf);   
  if (year!=0)fChain->SetBranchAddress("_HLT_IsoMu27", &_HLT_IsoMu27, &b__HLT_IsoMu27);  
  if (year!=0)fChain->SetBranchAddress("_HLT_Ele32_WPTight_Gsf", &_HLT_Ele32_WPTight_Gsf, &b__HLT_Ele32_WPTight_Gsf);
  // if (year!=0)fChain->SetBranchAddress("_HLT_Ele35_WPTight_Gsf", &_HLT_Ele35_WPTight_Gsf, &b__HLT_Ele35_WPTight_Gsf);
  // if (year!=0)fChain->SetBranchAddress("_HLT_Ele32_WPTight_Gsf_L1DoubleEG", &_HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b__HLT_Ele32_WPTight_Gsf_L1DoubleEG);   
  if (year==1)fChain->SetBranchAddress("_HLT_Ele35_WPTight_Gsf", &_HLT_Ele35_WPTight_Gsf, &b__HLT_Ele35_WPTight_Gsf);
  if (year==1)fChain->SetBranchAddress("_HLT_Ele32_WPTight_Gsf_L1DoubleEG", &_HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b__HLT_Ele32_WPTight_Gsf_L1DoubleEG);   
  fChain->SetBranchAddress("_nL", &_nL, &b__nL);
  fChain->SetBranchAddress("_nMu", &_nMu, &b__nMu);
  fChain->SetBranchAddress("_nEle", &_nEle, &b__nEle);
  fChain->SetBranchAddress("_nLight", &_nLight, &b__nLight);
  fChain->SetBranchAddress("_nTau", &_nTau, &b__nTau);
  fChain->SetBranchAddress("_pvX", &_pvX, &b__pvX);
  fChain->SetBranchAddress("_pvY", &_pvY, &b__pvY);
  fChain->SetBranchAddress("_pvZ", &_pvZ, &b__pvZ);
  fChain->SetBranchAddress("_pvXErr", &_pvXErr, &b__pvXErr);
  fChain->SetBranchAddress("_pvYErr", &_pvYErr, &b__pvYErr);
  fChain->SetBranchAddress("_pvZErr", &_pvZErr, &b__pvZErr);
  fChain->SetBranchAddress("_nVFit_os", &_nVFit_os, &b__nVFit_os);
  fChain->SetBranchAddress("_nVFit", &_nVFit, &b__nVFit);
  fChain->SetBranchAddress("_nGoodLeading", &_nGoodLeading, &b__nGoodLeading);
  fChain->SetBranchAddress("_nGoodDisplaced", &_nGoodDisplaced, &b__nGoodDisplaced);
  fChain->SetBranchAddress("_vertices_os", _vertices_os, &b__vertices_os);
  fChain->SetBranchAddress("_lDisplaced_os", _lDisplaced_os, &b__lDisplaced_os);
  fChain->SetBranchAddress("_vertices", _vertices, &b__vertices);
  fChain->SetBranchAddress("_lDisplaced", _lDisplaced, &b__lDisplaced);
  fChain->SetBranchAddress("_lHasTrigger", _lHasTrigger, &b__lHasTrigger);
  fChain->SetBranchAddress("_lPt", _lPt, &b__lPt);
  fChain->SetBranchAddress("_lEta", _lEta, &b__lEta);
  fChain->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
  fChain->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
  fChain->SetBranchAddress("_lE", _lE, &b__lE);
  fChain->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
  fChain->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
  fChain->SetBranchAddress("_dxy", _dxy, &b__dxy);
  fChain->SetBranchAddress("_dz", _dz, &b__dz);
  fChain->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
  fChain->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
  fChain->SetBranchAddress("_2dIP", _2dIP, &b__2dIP);
  fChain->SetBranchAddress("_2dIPSig", _2dIPSig, &b__2dIPSig);
  fChain->SetBranchAddress("_lElectronPassEmu", _lElectronPassEmu, &b__lElectronPassEmu);
  fChain->SetBranchAddress("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, &b__lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto);
  fChain->SetBranchAddress("_lElectronPassConvVeto", _lElectronPassConvVeto, &b__lElectronPassConvVeto);
  fChain->SetBranchAddress("_lElectronChargeConst", _lElectronChargeConst, &b__lElectronChargeConst);
  fChain->SetBranchAddress("_lElectronMissingHits", _lElectronMissingHits, &b__lElectronMissingHits);
  fChain->SetBranchAddress("_lPOGVeto", _lPOGVeto, &b__lPOGVeto);
  fChain->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
  fChain->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
  fChain->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);
  fChain->SetBranchAddress("_lGlobalMuon", _lGlobalMuon, &b__lGlobalMuon);
  fChain->SetBranchAddress("_lTrackerMuon", _lTrackerMuon, &b__lTrackerMuon);
  fChain->SetBranchAddress("_lInnerTrackValidFraction", _lInnerTrackValidFraction, &b__lInnerTrackValidFraction);
  fChain->SetBranchAddress("_lGlobalTrackNormalizeChi2", _lGlobalTrackNormalizeChi2, &b__lGlobalTrackNormalizeChi2);
  fChain->SetBranchAddress("_lCQChi2Position", _lCQChi2Position, &b__lCQChi2Position);
  fChain->SetBranchAddress("_lCQTrackKink", _lCQTrackKink, &b__lCQTrackKink);
  fChain->SetBranchAddress("_lNumberOfMatchedStation", _lNumberOfMatchedStation, &b__lNumberOfMatchedStation);
  fChain->SetBranchAddress("_lNumberOfValidPixelHits", _lNumberOfValidPixelHits, &b__lNumberOfValidPixelHits);
  fChain->SetBranchAddress("_lTrackerLayersWithMeasurement", _lTrackerLayersWithMeasurement, &b__lTrackerLayersWithMeasurement);
  fChain->SetBranchAddress("_lSimType", _lSimType, &b__lSimType);
  fChain->SetBranchAddress("_lSimExtType", _lSimExtType, &b__lSimExtType);
  fChain->SetBranchAddress("_lSimFlavour", _lSimFlavour, &b__lSimFlavour);
  fChain->SetBranchAddress("_muDTStationsWithValidHits", _muDTStationsWithValidHits, &b__muDTStationsWithValidHits);
  fChain->SetBranchAddress("_muCSCStationsWithValidHits", _muCSCStationsWithValidHits, &b__muCSCStationsWithValidHits);
  fChain->SetBranchAddress("_muRPCStationsWithValidHits", _muRPCStationsWithValidHits, &b__muRPCStationsWithValidHits);
  fChain->SetBranchAddress("_muMuonStationsWithValidHits", _muMuonStationsWithValidHits, &b__muMuonStationsWithValidHits);
  fChain->SetBranchAddress("_lMuRPCTimenDof", _lMuRPCTimenDof, &b__lMuRPCTimenDof);
  fChain->SetBranchAddress("_lMuTimenDof", _lMuTimenDof, &b__lMuTimenDof);
  fChain->SetBranchAddress("_lMuRPCTime", _lMuRPCTime, &b__lMuRPCTime);
  fChain->SetBranchAddress("_lMuRPCTimeErr", _lMuRPCTimeErr, &b__lMuRPCTimeErr);
  fChain->SetBranchAddress("_lMuTime", _lMuTime, &b__lMuTime);
  fChain->SetBranchAddress("_lMuTimeErr", _lMuTimeErr, &b__lMuTimeErr);
  fChain->SetBranchAddress("_muNumberInnerHits", _muNumberInnerHits, &b__muNumberInnerHits);
  fChain->SetBranchAddress("_lEleIsEB", _lEleIsEB, &b__lEleIsEB);
  fChain->SetBranchAddress("_lEleIsEE", _lEleIsEE, &b__lEleIsEE);
  fChain->SetBranchAddress("_lEleSuperClusterOverP", _lEleSuperClusterOverP, &b__lEleSuperClusterOverP);
  fChain->SetBranchAddress("_lEleEcalEnergy", _lEleEcalEnergy, &b__lEleEcalEnergy);
  fChain->SetBranchAddress("_lElefull5x5SigmaIetaIeta", _lElefull5x5SigmaIetaIeta, &b__lElefull5x5SigmaIetaIeta);
  fChain->SetBranchAddress("_lEleDEtaInSeed", _lEleDEtaInSeed, &b__lEleDEtaInSeed);
  fChain->SetBranchAddress("_lEleDeltaPhiSuperClusterTrackAtVtx", _lEleDeltaPhiSuperClusterTrackAtVtx, &b__lEleDeltaPhiSuperClusterTrackAtVtx);
  fChain->SetBranchAddress("_lElehadronicOverEm", _lElehadronicOverEm, &b__lElehadronicOverEm);
  fChain->SetBranchAddress("_lEleInvMinusPInv", _lEleInvMinusPInv, &b__lEleInvMinusPInv);
  fChain->SetBranchAddress("_puCorr", _puCorr, &b__puCorr);
  fChain->SetBranchAddress("_absIso03", _absIso03, &b__absIso03);
  fChain->SetBranchAddress("_absIso04", _absIso04, &b__absIso04);
  fChain->SetBranchAddress("_sumNeutralHadronEt04", _sumNeutralHadronEt04, &b__sumNeutralHadronEt04);
  fChain->SetBranchAddress("_sumChargedHadronPt04", _sumChargedHadronPt04, &b__sumChargedHadronPt04);
  fChain->SetBranchAddress("_sumPhotonEt04", _sumPhotonEt04, &b__sumPhotonEt04);
  fChain->SetBranchAddress("_sumNeutralHadronEt03", _sumNeutralHadronEt03, &b__sumNeutralHadronEt03);
  fChain->SetBranchAddress("_sumChargedHadronPt03", _sumChargedHadronPt03, &b__sumChargedHadronPt03);
  fChain->SetBranchAddress("_sumPhotonEt03", _sumPhotonEt03, &b__sumPhotonEt03);
  fChain->SetBranchAddress("_trackIso", _trackIso, &b__trackIso);
  fChain->SetBranchAddress("_ecalIso", _ecalIso, &b__ecalIso);
  fChain->SetBranchAddress("_hcalIso", _hcalIso, &b__hcalIso);
  fChain->SetBranchAddress("_ecalPFClusterIso", _ecalPFClusterIso, &b__ecalPFClusterIso);
  fChain->SetBranchAddress("_hcalPFClusterIso", _hcalPFClusterIso, &b__hcalPFClusterIso);
  fChain->SetBranchAddress("_relIso", _relIso, &b__relIso);
  fChain->SetBranchAddress("_relIso0p4", _relIso0p4, &b__relIso0p4);
  fChain->SetBranchAddress("_relIso0p4MuDeltaBeta", _relIso0p4MuDeltaBeta, &b__relIso0p4MuDeltaBeta);
  fChain->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
  fChain->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
  fChain->SetBranchAddress("_closestJetCsvV2", _closestJetCsvV2, &b__closestJetCsvV2);
  fChain->SetBranchAddress("_closestJetDeepCsv_b", _closestJetDeepCsv_b, &b__closestJetDeepCsv_b);
  fChain->SetBranchAddress("_closestJEC", _closestJEC, &b__closestJEC);
  fChain->SetBranchAddress("_closest_lepAwareJetE", _closest_lepAwareJetE, &b__closest_lepAwareJetE);
  fChain->SetBranchAddress("_closest_lepAwareJetPx", _closest_lepAwareJetPx, &b__closest_lepAwareJetPx);
  fChain->SetBranchAddress("_closest_lepAwareJetPy", _closest_lepAwareJetPy, &b__closest_lepAwareJetPy);
  fChain->SetBranchAddress("_closest_lepAwareJetPz", _closest_lepAwareJetPz, &b__closest_lepAwareJetPz);
  fChain->SetBranchAddress("_closest_l1JetE", _closest_l1JetE, &b__closest_l1JetE);
  fChain->SetBranchAddress("_closest_l1JetPx", _closest_l1JetPx, &b__closest_l1JetPx);
  fChain->SetBranchAddress("_closest_l1JetPy", _closest_l1JetPy, &b__closest_l1JetPy);
  fChain->SetBranchAddress("_closest_l1JetPz", _closest_l1JetPz, &b__closest_l1JetPz);
  fChain->SetBranchAddress("_closest_lJetE", _closest_lJetE, &b__closest_lJetE);
  fChain->SetBranchAddress("_closest_lJetPx", _closest_lJetPx, &b__closest_lJetPx);
  fChain->SetBranchAddress("_closest_lJetPy", _closest_lJetPy, &b__closest_lJetPy);
  fChain->SetBranchAddress("_closest_lJetPz", _closest_lJetPz, &b__closest_lJetPz);
  fChain->SetBranchAddress("_closestJetDeepCsv_bb", _closestJetDeepCsv_bb, &b__closestJetDeepCsv_bb);    
  fChain->SetBranchAddress("_lElectronSummer16MvaGP", _lElectronSummer16MvaGP, &b__lElectronSummer16MvaGP);
  fChain->SetBranchAddress("_lElectronSummer16MvaHZZ", _lElectronSummer16MvaHZZ, &b__lElectronSummer16MvaHZZ);
  fChain->SetBranchAddress("_lElectronMvaFall17v1NoIso", _lElectronMvaFall17v1NoIso, &b__lElectronMvaFall17v1NoIso);
  fChain->SetBranchAddress("_lElectronMvaFall17Iso", _lElectronMvaFall17Iso, &b__lElectronMvaFall17Iso);
  fChain->SetBranchAddress("_lElectronMvaFall17NoIso", _lElectronMvaFall17NoIso, &b__lElectronMvaFall17NoIso);  
  fChain->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
  fChain->SetBranchAddress("_lMuonSegComp", _lMuonSegComp, &b__lMuonSegComp);
  fChain->SetBranchAddress("_lMuonTrackPt", _lMuonTrackPt, &b__lMuonTrackPt);
  fChain->SetBranchAddress("_lMuonTrackPtErr", _lMuonTrackPtErr, &b__lMuonTrackPtErr);   
  // fChain->SetBranchAddress("_lPtCorr", _lPtCorr, &b__lPtCorr);
  // fChain->SetBranchAddress("_lPtScaleUp", _lPtScaleUp, &b__lPtScaleUp);
  // fChain->SetBranchAddress("_lPtScaleDown", _lPtScaleDown, &b__lPtScaleDown);
  // fChain->SetBranchAddress("_lPtResUp", _lPtResUp, &b__lPtResUp);
  // fChain->SetBranchAddress("_lPtResDown", _lPtResDown, &b__lPtResDown);
  // fChain->SetBranchAddress("_lECorr", _lECorr, &b__lECorr);
  // fChain->SetBranchAddress("_lEScaleUp", _lEScaleUp, &b__lEScaleUp);
  // fChain->SetBranchAddress("_lEScaleDown", _lEScaleDown, &b__lEScaleDown);
  // fChain->SetBranchAddress("_lEResUp", _lEResUp, &b__lEResUp);
  // fChain->SetBranchAddress("_lEResDown", _lEResDown, &b__lEResDown);
  //if(!isData) fChain->SetBranchAddress("_nJets", &_nJets, &b__nJets);
  //if(isData) fChain->SetBranchAddress("_nJets", &_nJets_data, &b__nJets);
  fChain->SetBranchAddress("_nJets", &_nJets, &b__nJets);
  fChain->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
  fChain->SetBranchAddress("_jetPt_JECDown", _jetPt_JECDown, &b__jetPt_JECDown);
  fChain->SetBranchAddress("_jetPt_JECUp", _jetPt_JECUp, &b__jetPt_JECUp);
  fChain->SetBranchAddress("_jetSmearedPt", _jetSmearedPt, &b__jetSmearedPt);
  fChain->SetBranchAddress("_jetSmearedPt_JECDown", _jetSmearedPt_JECDown, &b__jetSmearedPt_JECDown);
  fChain->SetBranchAddress("_jetSmearedPt_JECUp", _jetSmearedPt_JECUp, &b__jetSmearedPt_JECUp);
  fChain->SetBranchAddress("_jetSmearedPt_JERDown", _jetSmearedPt_JERDown, &b__jetSmearedPt_JERDown);
  fChain->SetBranchAddress("_jetSmearedPt_JERUp", _jetSmearedPt_JERUp, &b__jetSmearedPt_JERUp);
  fChain->SetBranchAddress("_jetPt_Uncorrected", _jetPt_Uncorrected, &b__jetPt_Uncorrected);
  fChain->SetBranchAddress("_jetPt_L1", _jetPt_L1, &b__jetPt_L1);
  fChain->SetBranchAddress("_jetPt_L2", _jetPt_L2, &b__jetPt_L2);
  fChain->SetBranchAddress("_jetPt_L3", _jetPt_L3, &b__jetPt_L3);
  fChain->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
  fChain->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
  fChain->SetBranchAddress("_jetE", _jetE, &b__jetE);
  fChain->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
  fChain->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
  fChain->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
  fChain->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
  fChain->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
  fChain->SetBranchAddress("_jetHadronFlavor", _jetHadronFlavor, &b__jetHadronFlavor);
  //fChain->SetBranchAddress("_jetIsLoose", _jetIsLoose, &b__jetIsLoose);
  fChain->SetBranchAddress("_jetIsTight", _jetIsTight, &b__jetIsTight);
  fChain->SetBranchAddress("_jetIsTightLepVeto", _jetIsTightLepVeto, &b__jetIsTightLepVeto);
  fChain->SetBranchAddress("_jetNeutralHadronFraction", _jetNeutralHadronFraction, &b__jetNeutralHadronFraction);
  fChain->SetBranchAddress("_jetChargedHadronFraction", _jetChargedHadronFraction, &b__jetChargedHadronFraction);
  fChain->SetBranchAddress("_jetNeutralEmFraction", _jetNeutralEmFraction, &b__jetNeutralEmFraction);
  fChain->SetBranchAddress("_jetChargedEmFraction", _jetChargedEmFraction, &b__jetChargedEmFraction);
  fChain->SetBranchAddress("_jetHFHadronFraction", _jetHFHadronFraction, &b__jetHFHadronFraction);
  fChain->SetBranchAddress("_jetHFEmFraction", _jetHFEmFraction, &b__jetHFEmFraction);
  fChain->SetBranchAddress("_met", &_met, &b__met);
  // fChain->SetBranchAddress("_metRaw", &_metRaw, &b__metRaw);
  fChain->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
  fChain->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
  fChain->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
  fChain->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
  fChain->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
  // fChain->SetBranchAddress("_metRawPhi", &_metRawPhi, &b__metRawPhi);
  fChain->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
  fChain->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
  fChain->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
  fChain->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
  fChain->SetBranchAddress("_metSignificance", &_metSignificance, &b__metSignificance);

  if(!isData){
    fChain->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
    fChain->SetBranchAddress("_weight", &_weight, &b__weight);
    // fChain->SetBranchAddress("_lheHTIncoming", &_lheHTIncoming, &b__lheHTIncoming);
    if(isNewPhys) fChain->SetBranchAddress("_ctauHN", &_ctauHN, &b__ctauHN);
    // fChain->SetBranchAddress("_nLheTau", &_nLheTau, &b__nLheTau);
    fChain->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
    fChain->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
    // fChain->SetBranchAddress("_nPsWeights", &_nPsWeights, &b__nPsWeights);
    // fChain->SetBranchAddress("_psWeight", _psWeight, &b__psWeight);
    fChain->SetBranchAddress("_leptonNumber", &_leptonNumber, &b__leptonNumber);
    fChain->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
    // fChain->SetBranchAddress("_gen_pdgID", _gen_pdgID, &b__gen_pdgID);
    fChain->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
    fChain->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
    fChain->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
    fChain->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
    fChain->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
    fChain->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
    fChain->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
    // fChain->SetBranchAddress("_gen_vertex_x", _gen_vertex_x, &b__gen_vertex_x);
    // fChain->SetBranchAddress("_gen_vertex_y", _gen_vertex_y, &b__gen_vertex_y);
    // fChain->SetBranchAddress("_gen_vertex_z", _gen_vertex_z, &b__gen_vertex_z);
    // fChain->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
    // fChain->SetBranchAddress("_gen_lMinDeltaR", _gen_lMinDeltaR, &b__gen_lMinDeltaR);
    // fChain->SetBranchAddress("_gen_lPassParentage", _gen_lPassParentage, &b__gen_lPassParentage);    
    fChain->SetBranchAddress("_lGenIndex", _lGenIndex, &b__lGenIndex);
    fChain->SetBranchAddress("_lMatchType", _lMatchType, &b__lMatchType);
    fChain->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
    fChain->SetBranchAddress("_lIsPromptFinalState", _lIsPromptFinalState, &b__lIsPromptFinalState);
    fChain->SetBranchAddress("_lIsPromptDecayed", _lIsPromptDecayed, &b__lIsPromptDecayed);
    fChain->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
    fChain->SetBranchAddress("_lMomPdgId", _lMomPdgId, &b__lMomPdgId);
    fChain->SetBranchAddress("_lProvenance", _lProvenance, &b__lProvenance);
    fChain->SetBranchAddress("_lProvenanceCompressed", _lProvenanceCompressed, &b__lProvenanceCompressed);
    fChain->SetBranchAddress("_lProvenanceConversion", _lProvenanceConversion, &b__lProvenanceConversion);
    fChain->SetBranchAddress("_lMatchPt", _lMatchPt, &b__lMatchPt);
    fChain->SetBranchAddress("_lMatchEta", _lMatchEta, &b__lMatchEta);
    fChain->SetBranchAddress("_lMatchPhi", _lMatchPhi, &b__lMatchPhi);
    fChain->SetBranchAddress("_lMatchVertexX", _lMatchVertexX, &b__lMatchVertexX);
    fChain->SetBranchAddress("_lMatchVertexY", _lMatchVertexY, &b__lMatchVertexY);
    fChain->SetBranchAddress("_lMatchVertexZ", _lMatchVertexZ, &b__lMatchVertexZ);

    fChain->SetBranchAddress("_hasInternalConversion", &_hasInternalConversion, &b__hasInternalConversion);
    fChain->SetBranchAddress("_zgEventType", &_zgEventType, &b__zgEventType);
  }
}

//          ================= ================= ================= ================= ================= =================          // 
//          ================= ================= ================= ================= ================= =================          //

//_______________________________________________________ analysis function ____
void Analysis_mc::analisi( //const std::string& list, const std::string& directory,
			  std::string outfilename,
			  bool skipData, bool skipSignal, bool skipBackground,
			  bool skipPlotting, bool skipLimits, bool skipTables, bool skipXCheck 
			  /*int systcat, int systdir*/) {

  // std::ofstream zero("zero.txt"); 
  // std::ofstream one("one.txt");  
  // std::ofstream two("two.txt"); 
  // std::ofstream three("three.txt"); 
  // std::ofstream four("four.txt"); 
  // std::ofstream ratios_n_1("ratios_n_1.txt"); 
  // std::ofstream yields_check("yields_check.txt");
  std::ofstream single_fake_txt("single_fake.txt");


  cout<<"in analisi"<<endl;
  cout<<"---------------------------"<<endl;   
  setTDRStyle();

  // TEMPORARY --> this code SUCKS!!!!!!
  TString signProcName;

  // Are we running in local or on T2B?
  TString cwd(gSystem->pwd());
  bool ist2b = true;
  //bool ist2b = cwd.BeginsWith("/storage_mnt");
  //bool isdtl = (ist2b==false && cwd.Contains("trocino"));
  //if(isdtl==false) ist2b = true;

  // ------------ pile up -----------------------------------------------//
  //TH1D *pileUpWeight[1];
  // TFile *hfile_pu = ist2b ?
  //   TFile::Open("/user/mvit/CMSSW_9_4_4/src/HNL_analysis/PU/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root") :
  //   TFile::Open("/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/ANALYSIS/20190419_MartinasCode/HNL_analysis/PU/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root");
  //pileUpWeight[0] = (TH1D*)hfile_pu->Get("puw_Run2016Inclusive_central");

  const TString names_FR_files[5]= {"FR/" +std::to_string(year)+ "/fake_rate_mu.root",
				    "FR/" +std::to_string(year)+ "/fake_rate_e.root",
				    "FR/" +std::to_string(year)+ "/fake_rate_mumu_1eta.root",
				    "FR/" +std::to_string(year)+ "/fake_rate_ee_1eta.root",
				    "FR/" +std::to_string(year)+ "/fake_rate_emu_1eta.root"};
	
	
	


  // FR histograms
  TGraphAsymmErrors *fakeRate_mu[3];
  TGraphAsymmErrors *fakeRate_e[3];
  TGraphAsymmErrors *fakeRate_mumu[3];
  TGraphAsymmErrors *fakeRate_ee[3];
  TGraphAsymmErrors *fakeRate_mue[3];
  TFile *hfile1 =  TFile::Open(names_FR_files[0]) ;
  fakeRate_mu[0] = (TGraphAsymmErrors*)hfile1->Get("fakeRate_mu_eta1");
  fakeRate_mu[1] = (TGraphAsymmErrors*)hfile1->Get("fakeRate_mu_eta2");
  fakeRate_mu[2] = (TGraphAsymmErrors*)hfile1->Get("fakeRate_mu_eta3");
  TFile *hfile2 =  TFile::Open(names_FR_files[1]) ;
  fakeRate_e[0] = (TGraphAsymmErrors*)hfile2->Get("fakeRate_e_eta1");
  fakeRate_e[1] = (TGraphAsymmErrors*)hfile2->Get("fakeRate_e_eta2");
  fakeRate_e[2] = (TGraphAsymmErrors*)hfile2->Get("fakeRate_e_eta3");
  TFile *hfile_dfr1=  TFile::Open(names_FR_files[2]) ;
  fakeRate_mumu[0]= (TGraphAsymmErrors*)hfile_dfr1->Get("fakeRate_mu_eta1");
  fakeRate_mumu[1]= (TGraphAsymmErrors*)hfile_dfr1->Get("fakeRate_mu_eta1");
  fakeRate_mumu[2]= (TGraphAsymmErrors*)hfile_dfr1->Get("fakeRate_mu_eta1");
  TFile *hfile_dfr2 = TFile::Open(names_FR_files[3]) ;
  fakeRate_ee[0]= (TGraphAsymmErrors*)hfile_dfr2->Get("fakeRate_e_eta1");
  fakeRate_ee[1]= (TGraphAsymmErrors*)hfile_dfr2->Get("fakeRate_e_eta1");
  fakeRate_ee[2]= (TGraphAsymmErrors*)hfile_dfr2->Get("fakeRate_e_eta1");
  TFile *hfile_dfr3 =  TFile::Open(names_FR_files[4]) ;
  fakeRate_mue[0]= (TGraphAsymmErrors*)hfile_dfr3->Get("fakeRate_emu_eta1");
  fakeRate_mue[1]= (TGraphAsymmErrors*)hfile_dfr3->Get("fakeRate_emu_eta1");
  fakeRate_mue[2]= (TGraphAsymmErrors*)hfile_dfr3->Get("fakeRate_emu_eta1");
 
	
  //   SF leptons histograms	
  TH2D *sf_prompt_muon[1];
  TH2D *sf_prompt_muon_syst[1];
  TH2F *sf_iso_ip_prompt_muon[1];
  TH2F *sf_iso_ip_prompt_muon_syst[1];
  TH2F *sf_trigger_muon[1];
  TH2F *sf_trigger_ele[1]; 
  TH2F *sf_prompt_ele[1];
  //non prompt
  TH2F *sf_isoID_nPMuon[1];
  TH2F *sf_isoID_nPMuon_syst[1];
  TH2F *sf_btag_eff[3];	//0 uds, 1 charm, 2 b

  //by luka
  TH2F *sf_sv_effcy_num[1];
  TH2F *sf_sv_effcy_den[1];
  
  for (int i = 0; i < 3; i++){   
    if (year == 0 && i!= 0) continue;
    if (year == 1 && i!= 1) continue;
    if (year == 2 && i!= 2) continue;

    //ID prompt muon
    TFile *hfile1_names_SF_muon_files = TFile::Open(names_SF_muon_files[i]);
    if (i == 0 ) sf_prompt_muon[0] = (TH2D*)hfile1_names_SF_muon_files->Get("NUM_MediumID_DEN_genTracks_eta_pt");
    if (i == 1 ) sf_prompt_muon[0] = (TH2D*)hfile1_names_SF_muon_files->Get("NUM_MediumID_DEN_genTracks_pt_abseta");
    if (i == 2 ) sf_prompt_muon[0] = (TH2D*)hfile1_names_SF_muon_files->Get("NUM_MediumID_DEN_TrackerMuons_pt_abseta");
    TFile *hfile1_names_SFSY_muon_files = TFile::Open(names_SFSY_muon_files[i]);
    if (i == 0 ) sf_prompt_muon_syst[0] = (TH2D*)hfile1_names_SFSY_muon_files->Get("NUM_MediumID_DEN_genTracks_eta_pt");
    if (i == 1 ) sf_prompt_muon_syst[0] = (TH2D*)hfile1_names_SFSY_muon_files->Get("NUM_MediumID_DEN_genTracks_pt_abseta_syst");
    if (i == 2 ) sf_prompt_muon_syst[0] = (TH2D*)hfile1_names_SFSY_muon_files->Get("NUM_MediumID_DEN_TrackerMuons_pt_abseta_syst");
    //trigger muon
    TFile *hfile1_names_trigger_muon_files = TFile::Open(names_trigger_muon_files[i]);
    if (i == 0){
      hfile1_names_trigger_muon_files->cd("IsoMu24_OR_IsoTkMu24_PtEtaBins");
      sf_trigger_muon[0] = (TH2F*)hfile1_names_trigger_muon_files->Get("IsoMu24_OR_IsoTkMu24_PtEtaBins/abseta_pt_ratio");
    }
    if (i == 1){
      hfile1_names_trigger_muon_files->cd("IsoMu27_PtEtaBins");
      sf_trigger_muon[0] = (TH2F*)hfile1_names_trigger_muon_files->Get("IsoMu27_PtEtaBins/abseta_pt_ratio");
    }
    if (i == 2){
      hfile1_names_trigger_muon_files->cd("IsoMu24_PtEtaBins");
      sf_trigger_muon[0] = (TH2F*)hfile1_names_trigger_muon_files->Get("IsoMu24_PtEtaBins/abseta_pt_ratio");
    }
    //ISO prompt muon
    TFile *hfile1_names_SF_isoIP_muon_files = TFile::Open(names_SF_isoIP_muon_files[i]);
    sf_iso_ip_prompt_muon[0]      = (TH2F*)hfile1_names_SF_isoIP_muon_files->Get("NUM_displacedIsoIP_DEN_displaced_abseta_pt");
    sf_iso_ip_prompt_muon_syst[0] = (TH2F*)hfile1_names_SF_isoIP_muon_files->Get("NUM_displacedIsoIP_DEN_displaced_abseta_pt_combined_syst");
    //trigger electron
    TFile *hfile1_names_trigger_ele_files = TFile::Open(names_trigger_ele_files[i]);
    sf_trigger_ele[0] = (TH2F*)hfile1_names_trigger_ele_files->Get("EGamma_SF2D");
    //ID prompt electron 
    TFile *hfile1_names_SF_ele_files = TFile::Open(names_SF_ele_files[i]);
    sf_prompt_ele[0] = (TH2F*)hfile1_names_SF_ele_files->Get("EGamma_SF2D");
    //ISO displaced muons
    TFile *hfile1_names_SF_isoID_nPMuon_files = TFile::Open(names_SF_isoID_nPMuon_files[i]);
    sf_isoID_nPMuon[0]      = (TH2F*)hfile1_names_SF_isoID_nPMuon_files->Get("NUM_displacedIso_DEN_trackerMuons_abseta_pt");
    sf_isoID_nPMuon_syst[0] = (TH2F*)hfile1_names_SF_isoID_nPMuon_files->Get("NUM_displacedIso_DEN_trackerMuons_abseta_pt_combined_syst");
    //ID displaced muons missing

    //--------------- btagging
    TFile *hfile1_names_btagging_eff_files = TFile::Open(names_btagging_eff_files[i]);
    sf_btag_eff[0] = (TH2F*)hfile1_names_btagging_eff_files->Get("bTagEff_loose_udsg");
    sf_btag_eff[1] = (TH2F*)hfile1_names_btagging_eff_files->Get("bTagEff_loose_charm");
    sf_btag_eff[2] = (TH2F*)hfile1_names_btagging_eff_files->Get("bTagEff_loose_beauty");
    //--------------- SV study
    TFile *hfile1_names_SV_eff_files = TFile::Open(names_SV_eff_files[i]);
    sf_sv_effcy_num[0] =  (TH2F*)hfile1_names_SV_eff_files->Get("datahistn0");
    sf_sv_effcy_den[0] =  (TH2F*)hfile1_names_SV_eff_files->Get("mchistn0");
  }//end loop
  
 
  
  // by tom
  double displEleVars[8] = {0.,0.,0.,0.,0.,0.,0.}; // 2016
  if (year == 0){
    displEleVars[0] = 1.0124869154860976;
    displEleVars[1] = 0.9663328195490085;
    displEleVars[2] = 0.8775724113625718;
    displEleVars[3] = 0.8849044007496563;
    displEleVars[4] = 0.87900667647424;
    displEleVars[5] = 0.8826087597277812;
    displEleVars[6] = 0.9181381563802627;
    displEleVars[7] = 1.049691544994153;

  }  
  if (year == 1){
    displEleVars[0] = 1.1278756580420497;
    displEleVars[1] = 1.043668938242718;
    displEleVars[2] = 0.9777185186569067;
    displEleVars[3] = 0.9874760079682131;
    displEleVars[4] = 0.9517916002237199;
    displEleVars[5] = 0.9863389679172274;
    displEleVars[6] = 0.9797529808714778;
    displEleVars[7] = 0.9842090048414325;

  }
  if (year == 2){
    displEleVars[0] = 1.1190865037865663;
    displEleVars[1] = 1.0166429650926398;
    displEleVars[2] = 0.92393001102722;
    displEleVars[3] = 0.8820411147655539;
    displEleVars[4] = 0.8824176687641517;
    displEleVars[5] = 0.9403542112673093;
    displEleVars[6] = 1.0016620225340234;
    displEleVars[7] = 0.9688345234979848;
  }

  

  double stat_DF_factors[3]= {0.,0.,0.};
  if (year == 0){
    stat_DF_factors[0]= 0.17;
    stat_DF_factors[1]= 0.2;
    stat_DF_factors[2]= 0.4;
  }
  if (year == 1){
    stat_DF_factors[0]= 0.14;
    stat_DF_factors[1]= 0.14;
    stat_DF_factors[2]= 0.33;
  }
  if (year == 2){
    stat_DF_factors[0]= 0.1;
    stat_DF_factors[1]= 0.11;
    stat_DF_factors[2]= 0.25;
  }
  // ------------ b tagging -----------------------------------------------//
  // b-tagging working points (DeepCsv_b + DeepCsv_bb)
  BTagEntry::OperatingPoint bwp = BTagEntry::OP_LOOSE;    // = 0

  // B-tagging calibration + reader
  BTagCalibration calib("DeepCSV", (year==0 ? "DeepCSV_2016LegacySF_WP_V1.csv" : (year==1 ? "DeepCSV_94XSF_WP_V4_B_F.csv" : "DeepCSV_102XSF_WP_V1.csv")));
  BTagCalibrationReader reader(bwp,             // working point
			       "central",       // central sys type
			       {"up", "down"}); // other sys types

  reader.load(calib,             // calibration instance
	      BTagEntry::FLAV_B, // b-tag flavor
	      "comb");           // measurement type
  reader.load(calib,             // calibration instance
	      BTagEntry::FLAV_C, // b-tag flavor
	      "comb");           // measurement type
  reader.load(calib,             // calibration instance
	      BTagEntry::FLAV_UDSG, // b-tag flavor
	      "incl");           // measurement type



  // // ------------   samples info -----------------------------------------------//
  // samples = readSampleList(list, directory);
  // std::cout << "I am in the analysis number:  " << systcat << std::endl;
  // pdf!
  std::vector<unsigned> qcdSystVars, pdfSystVars;
  //bool runtheosyst = (systcat==2 || systcat==3);
  //if(systcat==2) {
  qcdSystVars.push_back(2);
  qcdSystVars.push_back(3);
  qcdSystVars.push_back(4);
  qcdSystVars.push_back(5);
  qcdSystVars.push_back(7);
  qcdSystVars.push_back(9);
  //}
  //else if(systcat==3) {
  for(unsigned l=10; l<110; ++l)
    pdfSystVars.push_back(l);
  //}
  const unsigned nQcdVars = qcdSystVars.size();
  const unsigned nPdfVars = pdfSystVars.size();

  /*TH1D* qcdHistos[nQcdVars][nCoupling][nSamples_eff+1];
    TH1D* qcdHistosNorm[nQcdVars][nCoupling][nSamples_eff+1];
    TH1D* pdfHistos[nPdfVars][nCoupling][nSamples_eff+1];
    TH1D* pdfHistosNorm[nPdfVars][nCoupling][nSamples_eff+1];
    //if(runtheosyst) {
    float binWidth = (HistMax[0] - HistMin[0])/nBins[0];
    std::ostringstream strs; strs << binWidth; std::string Yaxis = strs.str();
    for(size_t effsam=0; effsam<nSamples_eff+1; ++effsam) {
    for(int cha=0; cha<nCoupling; ++cha) {
    //if(cha!=6 && cha!=7) continue;
    // Only for theory systs
    for(unsigned sidx=0; sidx<nQcdVars; ++sidx) {
    qcdHistos    [sidx][cha][effsam] = new TH1D(eff_names[effsam] + "_qcdSyst_" + std::to_string(qcdSystVars[sidx]) + "_" + chaNames[cha] + "_" + Histnames_ossf[0] , eff_names[effsam] + "_qcdSyst_" + std::to_string(qcdSystVars[sidx]) + "_" + Histnames_ossf[0] + ";" + Xaxes[0] + ";events/" + Yaxis + Units[0], nBins[0], HistMin[0], HistMax[0]);
    qcdHistosNorm[sidx][cha][effsam] = new TH1D(eff_names[effsam] + "_qcdSystNorm_" + std::to_string(qcdSystVars[sidx]) + "_" + chaNames[cha] + "_" + Histnames_ossf[0] , eff_names[effsam] + "_qcdSystNorm_" + std::to_string(qcdSystVars[sidx]) + "_" + Histnames_ossf[0] + ";" + Xaxes[0] + ";events", 1, 0., 2.);
    qcdHistos    [sidx][cha][effsam]->Sumw2();
    qcdHistosNorm[sidx][cha][effsam]->Sumw2();
    }
    for(unsigned sidx=0; sidx<nPdfVars; ++sidx) {
    pdfHistos    [sidx][cha][effsam] = new TH1D(eff_names[effsam] + "_pdfSyst_" + std::to_string(pdfSystVars[sidx]) + "_" + chaNames[cha] + "_" + Histnames_ossf[0] , eff_names[effsam] + "_pdfSyst_" + std::to_string(pdfSystVars[sidx]) + "_" + Histnames_ossf[0] + ";" + Xaxes[0] + ";events/" + Yaxis + Units[0], nBins[0], HistMin[0], HistMax[0]);
    pdfHistosNorm[sidx][cha][effsam] = new TH1D(eff_names[effsam] + "_pdfSystNorm_" + std::to_string(pdfSystVars[sidx]) + "_" + chaNames[cha] + "_" + Histnames_ossf[0] , eff_names[effsam] + "_pdfSystNorm_" + std::to_string(pdfSystVars[sidx]) + "_" + Histnames_ossf[0] + ";" + Xaxes[0] + ";events", 1, 0., 2.);
    pdfHistos    [sidx][cha][effsam]->Sumw2();
    pdfHistosNorm[sidx][cha][effsam]->Sumw2();
    }
    }
    }*/
  //}
	
  // ------------   run over samples -----------------------------------------------//
  std::set<std::tuple<long, long, long> > usedEvents;
  for(size_t sam=0, effsam=0; sam<samples.size(); ++sam, ++effsam) {
    // >>>>>> FILTER OUT UNWANTED SAMPLES!!!!!!! <<<<<<
    // NOTE: this assumes that the list of samples in the txt follows a strict order
    //       (data, signal, background). Otherwise you need to check the file name,
    //       as is done in initSample(...)
    //if(sam==0 || sam>nSamples_signal) continue; // only signal, sam = 1-20 (10 mu + 10 e samples)
    if(skipData       && sam==0) continue;
    if(skipSignal     && sam>0 && sam<=nSamples_signal) continue;
    if(skipBackground && sam>nSamples_signal) continue;


    initSample(samples[sam]);
    TH1D *hLheCounter = nullptr;
    if(!samples[sam].isData()){
      //read LHE weights of simulated events
      hLheCounter = new TH1D("lheCounter", "Events counter", 110, 0., 110.);
      hLheCounter->Read("lheCounter");
    }

    //check consistency
    std::cout << "sample initialized: --> " << std::endl;
    std::cout << "fileName: " << samples[sam].getFileName() << "  process name: " << samples[sam].getProcessName() << "   xsec: " << samples[sam].getXSec() << std::endl;
    if(samples[sam].isData()) std::cout << " is Data" << std::endl;
    if(samples[sam].isMC()  ) std::cout << " is MC"   << std::endl;
    if(sam != 0){
      if(samples[sam].getProcessName() == samples[sam-1].getProcessName()) --effsam;     
    }
    //if (samples[sam].isData()) continue; 
    std::cout<<	"sam.is2018() "<< samples[sam].is2018()  <<std::endl;
    std::cout<<	"sam.is2017() "<< samples[sam].is2017()  <<std::endl;

      
    bool isSignal = samples[sam].isNewPhysicsSignal();   
    if (eff_names[effsam] == "ttbar") continue;  
    if (eff_names[effsam] == "WJets") continue;

    // TEMPORARY ---> this code SUCKS!!!!!!!!!!!!
    if(signProcName.Length()==0) signProcName = samples[sam].getProcessName();

    bool transformCtau = false;
    bool skipSignalLNV = false;
    bool useAllMajEvts = false;

    
    // For lifetime re-weighting (hip hip hip hurray)
    double ctauOld(0.), ctauNew(0.); //, ctWeight(1.);
    if(isSignal) {
      std::cout << " is signal" << std::endl;
      //if(samples[sam].getHNLV2New()>0.) {
      if(samples[sam].getHNLctauNew()>0.) {
	transformCtau = true;
	skipSignalLNV = samples[sam].useLNCeventsOnly();
	useAllMajEvts = samples[sam].isMajoranaToDiracSimul() && !skipSignalLNV;
	ctauOld = samples[sam].getHNLctau();
	ctauNew = samples[sam].getHNLctauNew();
	std::cout << "  ==> HNL lifetime re-weighting: " << std::endl;
	std::cout << "      (" << samples[sam].getHNLV2()    << ", " << ctauOld
		  << ") --> (" << samples[sam].getHNLV2New() << ", " << ctauNew
		  << ")" << std::endl;

	//ctWeight = (ctauOld/ctauNew) * TMath::Exp(((1./ctauOld)-(1./ctauNew))*_ctauHN);
      }
    }
    // ------------   run over entries -----------------------------------------------//  
   	  
    for(ULong64_t it=0; it<nEntries; ++it) {
      GetEntry(samples[sam], it);
      // Skip event if
      //  (1) it's a Majorana sample used to simulate a Dirac one,
      //  (2) it's a LNV event, and
      //  (3) we only want to use LNC events (half statistics)
      if(skipSignalLNV && _leptonNumber!=0) continue;

      if (samples[sam].isData()){
	auto event = usedEvents.find(std::make_tuple(_eventNb, _lumiBlock, _runNb));
	if(event != usedEvents.end()) continue;
	usedEvents.insert(std::make_tuple(_eventNb, _lumiBlock, _runNb));
      }	    
	
      for(int var = 0; var < nVariation; ++var){
	for (int syst = 0; syst < nSystematic; ++syst)	{
	  for(int cha = 0; cha < nCoupling; ++cha){
	    weight_SR[cha][syst][var][effsam]=1.;
	  }
	}
      }	    

      double ctWeight(1.);
      //if(isSignal && samples[sam].getHNLV2New()>0.) {
      if(transformCtau) {
	ctWeight = (ctauOld/ctauNew) * TMath::Exp(((1./ctauOld)-(1./ctauNew))*_ctauHN);
      }

      // N.B.: ctWeight = 1 unless it is a ctau-reweighted signal sample
      //ctWeight = 1;
      double scal = 0;
      //scal = scale * _weight * ctWeight * pu_weight(*&pileUpWeight[0],_nTrueInt);
      scal = scale * _weight * ctWeight; 
      bwght = 1.;
      if (samples[sam].isData()) scal =1.;
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VARIABLES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      std::vector<unsigned> ind; 
      unsigned         ind_new_leading=0;
      unsigned        l1,l2,l3=0;
      TLorentzVector  v4l1,v4l2,v4l3;
      TLorentzVector  v4l2_propagated,v4l3_propagated;
      double            _vertex_X,_vertex_Y,_vertex_Z,_vertex_r, _vertex_phi,_vertex_eta, _vertex_theta, _vertex_R,_vertex_errX,_vertex_errY,_vertex_errZ,_vertex_chi2,_vertex_normchi2,_vertex_ndf=-1;
 
      // initializations!!!
      goodjet=0;
      bjet=0;
      promptC = 0;
      iV_ls=0;
      iV_lt=0;
      iV_st=0;  
      _mll_min=50000;
      METvec.SetPtEtaPhiE(0.,0.,0.,0.);
      sum_3l_rec.SetPtEtaPhiE(0.,0.,0.,0.);
      sum_2l_rec_pair.SetPtEtaPhiE(0.,0.,0.,0.);
      other[0].SetPtEtaPhiE(0.,0.,0.,0.);
      kind[0]=-1;   
      skip_event[0]= -1;
      for (int i =0; i < 3; i++){
	lepton_reco[i].SetPtEtaPhiE(0.,0.,0.,0.);
	lepton_transv[i].SetPtEtaPhiE(0.,0.,0.,0.);
	if (i !=2)pair[i].SetPtEtaPhiE(0.,0.,0.,0.);
	flavors_3l[i]=0;
	charge_3l[i]=0;	
      }
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //------------------------------------------------------------ lepton selection for FO
      const unsigned lCount = selectLepConeCorr(ind);
      if (lCount < 3) continue;
      //------------------------------------------------------------ jet pt variation and nJet and bjet
      //counting bjet and njet
      for (unsigned j =0; j < _nJets ; j++){
	_jetPt[j]=_jetSmearedPt[j];      
	if (jetIsGood(j, _jetPt[j])) ++goodjet;
	if (jetIsBJet(j, _jetPt[j])) ++bjet;   
      }
      // ------------ ==================== -----------------------------------------------//
      // ------------   event selection   -----------------------------------------------//
      //assign the l1 index
      ind_new_leading = l1Index(ind);
      if (l1Index(ind) == -1) continue; //in case there are not l1 at all
      //check how many displaced there are (displaced --> dxy, common vertex, FO, no l1)
      unsigned displacedC = 0;
      std::vector<TLorentzVector> lepV_displaced;
      std::vector<int> charge_displaced;
      std::vector<unsigned> temp_index;
          
      int index_to_use_for_l2_l3[2]={0,0};
      //find the right OS pair with min invariant mass
      //int min_test= 9999;
      double min_mass=0; // 999 if minMass idea is used
      displacedC=0;
      for(unsigned l = 0; l < lCount-1; ++l){
	for(unsigned j = l+1; j < lCount; ++j){	  	
	  if (!IsDisplacedPair(ind[l], ind[j], ind_new_leading, ind)) continue; 
	  ++displacedC;
	  TLorentzVector temp_displaced1;
	  TLorentzVector temp_displaced2;
	  temp_displaced1.SetPtEtaPhiE(_lPt[ind[l]],_lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
	  temp_displaced2.SetPtEtaPhiE(_lPt[ind[j]],_lEta[ind[j]], _lPhi[ind[j]], _lE[ind[j]]);
	  if ( (temp_displaced1+temp_displaced2).Pt()  > min_mass) {
	    min_mass= (temp_displaced1+temp_displaced2).Pt();
	    if (_lPt[ind[l]]> _lPt[ind[j]]){
	      index_to_use_for_l2_l3[0] = ind[l];
	      index_to_use_for_l2_l3[1] = ind[j];
	    }		    
	    else{
	      index_to_use_for_l2_l3[0] = ind[j];
	      index_to_use_for_l2_l3[1] = ind[l];
	    }	    
	  }
	}//end loop2
      }//end loop1

      if (displacedC < 1) continue;     
      //trigger NOT trigger matching!!!!!!
      if (!_passTrigger_1l) continue;

      // ------------ changing all the lep info and vertex-----------------------------------------------//
      l1 =-1;
      l2 =-1;
      l3 =-1; 
      l1=ind_new_leading;
      l2=index_to_use_for_l2_l3[0];
      l3=index_to_use_for_l2_l3[1];
      v4l1.SetPtEtaPhiE(_lPt[l1],_lEta[l1], _lPhi[l1], _lE[l1]);
      // first definition 
      v4l2.SetPtEtaPhiE(_lPt[l2],_lEta[l2], _lPhi[l2], _lE[l2]);
      v4l3.SetPtEtaPhiE(_lPt[l3],_lEta[l3], _lPhi[l3], _lE[l3]);      
      //-----------------------------------------------
      bool single_fake=false;
      bool Double_fake=false;     
     
      TLorentzVector  lJet[2] ;
      float JECl2       ;
      float JECl3      ;
      TLorentzVector  lepAwareJetl2[1] ;
      TLorentzVector  lepAwareJetl3[1] ;
      lJet[0].SetPxPyPzE(_closest_l1JetPx[l2],_closest_l1JetPy[l2],_closest_l1JetPz[l2],_closest_l1JetE[l2]);
      lJet[1].SetPxPyPzE(_closest_l1JetPx[l3],_closest_l1JetPy[l3],_closest_l1JetPz[l3],_closest_l1JetE[l3]);
      JECl2             = _closestJEC[l2];
      JECl3             = _closestJEC[l3];
      lepAwareJetl2[0] = (lJet[0] - v4l2 )*JECl2 + v4l2;
      lepAwareJetl3[0] = (lJet[1] - v4l3 )*JECl3 + v4l3;
       // ------------ closest jet info --------------------------------------//
      TLorentzVector  l1Jet[1] ;
      float JEC       ;
      TLorentzVector  lepAwareJet[1] ;
      l1Jet[0].SetPxPyPzE(_closest_l1JetPx[l2],_closest_l1JetPy[l2],_closest_l1JetPz[l2],_closest_l1JetE[l2]);
      JEC             = _closestJEC[l2];
      lepAwareJet[0] = (l1Jet[0] - v4l2 - v4l3)*JEC + v4l3 + v4l2;  
     
      if ( _closest_l1JetE[l2] ==  _closest_l1JetE[l3] ) Double_fake = true;
      if (v4l2.DeltaR(v4l3) < 0.45) {
	Double_fake = true;
	lepAwareJet[0] = lepAwareJetl2[0] + lepAwareJetl3[0];
      }
      if(v4l2.DeltaR(v4l3) > 0.45 && (_closest_l1JetE[l2] ==0 || _closest_l1JetE[l3] ==0)) {
	Double_fake = true;
	lepAwareJet[0] = lepAwareJetl2[0] + lepAwareJetl3[0];
      }
      if (!Double_fake) single_fake = true;
      
      double momentum_jet=0.;
      momentum_jet = lepAwareJet[0].Pt();
      if (momentum_jet<10) momentum_jet=12;
      
      //-----------------------------------------------
      if (single_fake)  v4l2.SetPtEtaPhiE(_lPt[l2]*(1+std::max(_relIso[l2]-0.2,0.)),_lEta[l2], _lPhi[l2], _lE[l2]*(1+std::max(_relIso[l2]-0.2,0.)));
      if (single_fake)  v4l3.SetPtEtaPhiE(_lPt[l3]*(1+std::max(_relIso[l3]-0.2,0.)),_lEta[l3], _lPhi[l3], _lE[l3]*(1+std::max(_relIso[l3]-0.2,0.)));
      flavors_3l[0]=_lFlavor[l1];
      flavors_3l[1]=_lFlavor[l2];
      flavors_3l[2]=_lFlavor[l3];
      charge_3l[0]=_lCharge[l1];
      charge_3l[1]=_lCharge[l2];
      charge_3l[2]=_lCharge[l3];
      // Flip l1 charge if (1) Majorana-->Dirac signal and (2) it's a LNV event
      if(useAllMajEvts && _leptonNumber!=0)
	charge_3l[0] *= (-1);

      //vertex l2l3 info
      int index_l2l3= l2l3_vertex_variable (l2,l3);      
      _vertex_X=_vertices[index_l2l3][1];
      _vertex_Y=_vertices[index_l2l3][2];
      _vertex_Z=_vertices[index_l2l3][3];
      _vertex_errX=TMath::Sqrt(_vertices[index_l2l3][4]);
      _vertex_errY=TMath::Sqrt(_vertices[index_l2l3][5]);
      _vertex_errZ=TMath::Sqrt(_vertices[index_l2l3][6]);
      _vertex_chi2=_vertices[index_l2l3][11];
      _vertex_normchi2= _vertices[index_l2l3][11]/_vertices[index_l2l3][10];
      _vertex_ndf =_vertices[index_l2l3][10];

      _vertex_r = TMath::Sqrt(_vertex_X*_vertex_X + _vertex_Y*_vertex_Y);
      _vertex_R = TMath::Sqrt(_vertex_X*_vertex_X + _vertex_Y*_vertex_Y + _vertex_Z*_vertex_Z);

      if (_vertex_r == 0)  _vertex_phi = 0.;
      else _vertex_phi= (_vertex_Y >= 0) ? TMath::ACos(_vertex_X/_vertex_r) : -(TMath::ACos(_vertex_X/_vertex_r));
      if (_vertex_R == 0) _vertex_theta = -10;
      else _vertex_theta = TMath::ACos (_vertex_Z/ _vertex_R);
      if (_vertex_R == 0) _vertex_eta = -10;
      else _vertex_eta = -TMath::Log(TMath::Tan(_vertex_theta*0.5));

      
      // ------------ ==================== -----------------------------------------------//
      // ------------   tight selection   -----------------------------------------------//
      unsigned* _isT= new unsigned[_nL];
      unsigned* _isT_prompt= new unsigned[_nL];  // only for CR--> so higher pT threshold and the same for mu and e
      for(unsigned l = 0; l < lCount; ++l){
	_isT[ind[l]] = false;
	_isT_prompt[ind[l]] = false;
      }
      int tightC=0;
      if (lepIsTightDisplaced(l2)) _isT[l2] = true;
      if (lepIsTightDisplaced(l3)) _isT[l3] = true;
      if (_isT[l2]) tightC++;
      if (_isT[l3]) tightC++;
      // -----------------------------------------------------------//
      if (single_fake && flavors_3l[1] == 1 && v4l2.Pt() < 5) continue;
      if (single_fake && flavors_3l[2] == 1 && v4l3.Pt() < 5) continue;
      if (!isSignal && single_fake && flavors_3l[1] == 0 && v4l2.Pt() < 10) continue;
      if (!isSignal && single_fake && flavors_3l[2] == 0 && v4l3.Pt() < 10) continue;
     
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     analysis   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      // ----------------- conversion overlap removal    --------------------------------//
      if ((samples[sam].getFileName() == "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root"  || samples[sam].getFileName() == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root"  ||
	   samples[sam].getFileName() == "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v10_Fall17.root"  || samples[sam].getFileName() == "DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_realistic_v14_Fall17.root"  ||
	   samples[sam].getFileName() == "DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_Autumn18.root"  || samples[sam].getFileName() == "DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_Autumn18.root") &&   _zgEventType>=3 ) continue;
      
      if ((samples[sam].getFileName()== "ZGToLLG_01J_5f_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" || samples[sam].getFileName()== "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_realistic_v14_Fall17.root" ||  samples[sam].getFileName()== "ZGToLLG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8_Autumn18.root")   &&   _hasInternalConversion == 1 ) continue;

      
      // -----------------   function useful    --------------------------------//
      zCandidate( pair,other, v4l1, v4l2, v4l3, flavors_3l, charge_3l);
      // -----------------   variables useful    --------------------------------//
      double min_delta_phi = 0;
      min_delta_phi = fabs(v4l1.DeltaPhi(v4l2));
      if (fabs(v4l1.DeltaPhi(v4l3)) < min_delta_phi)  min_delta_phi = fabs(v4l1.DeltaPhi(v4l3));
      //vertex
      TVector3 primary_vertex[1];
      TVector3 secondary_vertex[1];
      primary_vertex[0].SetXYZ(_pvX,_pvY,_pvZ);
      secondary_vertex[0].SetXYZ(_vertex_X,_vertex_Y,_vertex_Z);
      double D3_delta_pv_sv=  (primary_vertex[0] - secondary_vertex[0]).Mag();
      double D2_delta_pv_sv= sqrt(  (primary_vertex[0].X()-secondary_vertex[0].X())*(primary_vertex[0].X()-secondary_vertex[0].X())   +    (primary_vertex[0].Y()-secondary_vertex[0].Y())*(primary_vertex[0].Y()-secondary_vertex[0].Y()) );
      double prob_vertex= TMath::Prob(_vertex_chi2,_vertex_ndf);
      TVector3 l2plusl3=  (v4l2 + v4l3).Vect().Unit();
      //l2plusl3.SetXYZ(l2plusl3.X(), l2plusl3.Y(), 0);
      TVector3 svMpv =secondary_vertex[0]- primary_vertex[0];
      // svMpv.SetXYZ(svMpv.X(), svMpv.Y(), 0);
      double vtxR     = svMpv.Mag();
      //double vtxRvtxPcosAlpha = svMpv.Dot(l2plusl3)/vtxR;	// previous UGent definition    
      double D2_delta_pv_svSig= D2_delta_pv_sv*D2_delta_pv_sv/(TMath::Sqrt(_vertex_X*_vertex_X*_vertex_errX*_vertex_errX  +   _vertex_Y*_vertex_Y*_vertex_errY*_vertex_errY));
      double D3_delta_pv_svSig= D3_delta_pv_sv*D3_delta_pv_sv/(TMath::Sqrt(_vertex_X*_vertex_X*_vertex_errX*_vertex_errX  +   _vertex_Y*_vertex_Y*_vertex_errY*_vertex_errY +   _vertex_Z*_vertex_Z*_vertex_errZ*_vertex_errZ));
      TVector3 perp;
      perp.SetXYZ(   (v4l2+v4l3).X(),   (v4l2+v4l3).Y(),   0.);
      TVector3 vperp;
      vperp.SetXYZ(   svMpv.X(),   svMpv.Y(),   0.);
      double vtxRvtxPcosAlpha = vperp.Dot(perp)/(vperp.Mag()*perp.Mag());    // new ETH definition -> going 2D


      //conversion variable
      double R1_l2 = v4l2.Pt()/(0.3*3.8)*100.0; //cm
      double displ2 = TMath::Sqrt(2.*fabs(_dxy[l2])*R1_l2+_dxy[l2]*_dxy[l2]);
      double R1_l3 = v4l3.Pt()/(0.3*3.8)*100.0; //cm
      double displ3 = TMath::Sqrt(2.*fabs(_dxy[l3])*R1_l3+_dxy[l3]*_dxy[l3]);

      
      // -----------------   masses, mT and pT
      double M_ZPair = (pair[0]+pair[1]).M();
      double M_l2l3 = (v4l2 + v4l3).M();
      double M_3L_combined = (v4l2 + v4l3 + v4l1).M();
      double M_l2l3_combined = (v4l2 + v4l3).M();
      double M_l1l2_combined = (v4l2 + v4l1).M();
      double M_l1l3_combined = (v4l3 + v4l1).M();
      
      METvec.SetPtEtaPhiE(_met, 0, _metPhi,_met);    
      TLorentzVector to_use_mT;
      to_use_mT.SetPtEtaPhiE(other[0].Pt(),0, other[0].Phi(), other[0].Pt());
      double mT=(to_use_mT+METvec).M();
      double sumjet=0.;
      for (unsigned j =0; j < _nJets ; j++){
	if (jetIsGood(j, _jetPt[j])){
	  sumjet= sumjet+ _jetPt[j];
	}	
      }
      double met_sumjet=_met+sumjet;    
      double sum_vec_l2l3=0.;
      sum_vec_l2l3 = (v4l2+v4l3).Pt();
      // -----------------   function useful  2 --> SR also    --------------------------------//
      // 0 = mmm
      // 1 = mme OS
      // 2 = mme SS
      // 3 = eee
      // 4 = eem OS
      // 5 = eem SS   
      int SR_channel=0;
      SR_channel=channel(flavors_3l, charge_3l);
      if (SR_channel == -1 ) continue;	    
      //avoid +++ or ---
      if (SR_channel == 0 && charge_3l[0] == charge_3l[1] && charge_3l[0] == charge_3l[2]) continue;
      if (SR_channel == 3 && charge_3l[0] == charge_3l[1] && charge_3l[0] == charge_3l[2]) continue;  
      //bin histogram SR
      int bin_SR_muonCoupling =0;
      int bin_SR_eleCoupling =0; 
      bin_SR_muonCoupling = SR_bin_muon( SR_channel, D2_delta_pv_sv,  M_l2l3_combined );
      bin_SR_eleCoupling =  SR_bin_ele( SR_channel, D2_delta_pv_sv,  M_l2l3_combined  );
      
      //----------- trigger checks according to the channel -->
      if (SR_channel <= 2 && !(_HLT_IsoMu27 || _HLT_IsoMu24 || _HLT_IsoTkMu24)) continue;
      if (SR_channel > 2 && !(_HLT_Ele27_WPTight_Gsf || _HLT_Ele32_WPTight_Gsf || _HLT_Ele35_WPTight_Gsf  || _HLT_Ele32_WPTight_Gsf_L1DoubleEG)) continue;

      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     sFR and  dRF   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
      // ------------ closest jet info --------------------------------------//
      int flav_dRF = -1;
      if (_lFlavor[l2]==1 && _lFlavor[l3]==1) flav_dRF=1;
      if (_lFlavor[l2]==0 && _lFlavor[l3]==0) flav_dRF=0;
      if ((_lFlavor[l2]==1 && _lFlavor[l3]==0) || (_lFlavor[l2]==0 && _lFlavor[l3]==1))  flav_dRF=2;
      int index_eta = 0;
      if(TMath::Abs(lepAwareJet[0].Eta()) < 0.8 ) index_eta = 1;
      else if(TMath::Abs(lepAwareJet[0].Eta()) < 1.479 )index_eta = 2;
      else index_eta = 3;
      // -----------------    variables for sFR and dFR    --------------------------------//
      bool tight_lepton_dFR = false;
      bool loose_lepton_dFR = false;
      if (_isT[l2] && _isT[l3]) tight_lepton_dFR = true;
      if (!tight_lepton_dFR) loose_lepton_dFR = true;
      bool tightFail_sFR=false;
      tightFail_sFR = (tightC < 2);
      //where the FR has to be applied
      bool sideBandRegion= false;
      if ( tightFail_sFR     && single_fake)     sideBandRegion= true;
      if ( loose_lepton_dFR  && Double_fake)     sideBandRegion= true;
      // ------------------ prompt check for MC ------------------------//
      promptC=0;
      if (_lIsPrompt[l1] || _lProvenanceCompressed[l1]==0) promptC++;
      if (_lIsPrompt[l2] || _lProvenanceCompressed[l2]==0) promptC++;
      if (_lIsPrompt[l3] || _lProvenanceCompressed[l3]==0) promptC++;
      if (!samples[sam].isData() && promptC!=3) continue;
      // -----------------    applying the FRs    --------------------------------//
      bool singleSF_negative= false;
      bool singleSF_positive= false;
      if (single_fake &&  tightC == 0) singleSF_negative = true;
      if (single_fake &&  tightC == 1) singleSF_positive = true;

      
      if (sideBandRegion){
	if ( samples[sam].isData() )   scal *= -1;
	if (!samples[sam].isData() )   scal  = 1 * scal;
	if (single_fake){
	  if (!_isT[l2]) {
	    double fr = FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
				   _lEta[l2], _lFlavor[l2], v4l2.Pt(), index_eta,flav_dRF, momentum_jet);
	    scal *= -fr/(1-fr);
	  }
	  if (!_isT[l3]) {
	    double fr = FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
				   _lEta[l3], _lFlavor[l3], v4l3.Pt(), index_eta,flav_dRF, momentum_jet);
	    scal *= -fr/(1-fr);
	  }	  
	}//sFR
	if (loose_lepton_dFR &&  Double_fake) {
	  double fr = FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
				 _lEta[l2], _lFlavor[l2], _lPt[l2], index_eta,flav_dRF, momentum_jet);
	  scal *= -fr/(1-fr);
	  if (SR_channel == 0) weight_SR[0][dfmm_index][1][effsam] = ((fr*(1+stat_DF_factors[0]))/(1-fr*(1+stat_DF_factors[0])))/scal;
	  if (SR_channel == 0) weight_SR[0][dfmm_index][2][effsam] = ((fr*(1-stat_DF_factors[0]))/(1-fr*(1-stat_DF_factors[0])))/scal;
	  if (SR_channel == 1 || SR_channel == 2 ) weight_SR[0][dfem_index][1][effsam] = ((fr*(1+stat_DF_factors[1]))/(1-fr*(1+stat_DF_factors[1])))/scal;
	  if (SR_channel == 1 || SR_channel == 2 ) weight_SR[0][dfem_index][2][effsam] = ((fr*(1-stat_DF_factors[1]))/(1-fr*(1-stat_DF_factors[1])))/scal;
	  if (SR_channel == 4 || SR_channel == 5 ) weight_SR[1][dfem_index][1][effsam] = ((fr*(1+stat_DF_factors[1]))/(1-fr*(1+stat_DF_factors[1])))/scal;
	  if (SR_channel == 4 || SR_channel == 5 ) weight_SR[1][dfem_index][2][effsam] =( (fr*(1-stat_DF_factors[1]))/(1-fr*(1-stat_DF_factors[1])))/scal;
	  if (SR_channel == 3) weight_SR[1][dfee_index][1][effsam] = ((fr*(1+stat_DF_factors[2]))/(1-fr*(1+stat_DF_factors[2])))/scal;
	  if (SR_channel == 3) weight_SR[1][dfee_index][2][effsam] = ((fr*(1-stat_DF_factors[2]))/(1-fr*(1-stat_DF_factors[2])))/scal;
	}
      }//FR

      if (single_fake && tightFail_sFR && !_isT[l2] && _relIso[l2] < isolation_tight) continue;
      if (single_fake && tightFail_sFR && !_isT[l3] && _relIso[l3] < isolation_tight) continue;  
      // 0 = mmm
      // 1 = mme OS
      // 2 = mme SS
      // 3 = eee
      // 4 = eem OS
      // 5 = eem SS   


      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      ANALYSIS SELECTIONS     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //-----------
      if (charge_3l[2] == charge_3l[1]) continue;     
      bool selection_0 = true;
      bool selection_final = false;
      bool vetoes=false;
      vetoes = resonanceVeto( SR_channel,  D2_delta_pv_sv, flavors_3l,charge_3l,  M_l2l3_combined , M_l1l2_combined,  M_l1l3_combined);
      
      bool SR_selection = false;  // bveto is not there because we want btagging SF  
      SR_selection = 	v4l2.DeltaR(v4l3) < 1 &&   
					    M_3L_combined > 50 && 
	M_3L_combined < 80 && 
			min_delta_phi > 1 &&
	vtxRvtxPcosAlpha > 0.99  &&
	M_l2l3_combined < 20 &&
			  sum_vec_l2l3 > 15 &&
	D2_delta_pv_svSig > 20 &&	
	prob_vertex > 0.001 &&
	//	bjet == 0 &&
	vetoes;
      if (SR_selection && bjet == 0) selection_final = true; // then it would be == SR_plot[central]
      if (SR_channel == 3 && M_l2l3_combined < 0.5) continue;
      //cut flow bools
      bool _0,_1,_2,_3,_4,_5,_6,_7,_8,_9, _10=false;
      bool sel_0,sel_1,sel_2,sel_3,sel_4,sel_5,sel_6,sel_7,sel_8,sel_9, sel_10=false;
      sel_0 = selection_0;
      sel_1 = selection_0 && v4l2.DeltaR(v4l3) < 1;
      sel_2 = sel_1 && min_delta_phi > 1;
      sel_3 = sel_2 && bjet == 0;
      sel_4 = sel_3 && vetoes;
      sel_5 = sel_4 && M_3L_combined > 50 && M_3L_combined < 80;
      sel_6 = sel_5 && sum_vec_l2l3 > 15;
      sel_7 = sel_6 && vtxRvtxPcosAlpha > 0.99;
      sel_8 = sel_7 && prob_vertex > 0.001;
      sel_9 = sel_8 && D2_delta_pv_svSig > 20;
      sel_10= sel_9 && M_l2l3_combined < 20;
      _0 = SR_selection;
      _1 = M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&bjet == 0 &&vetoes;
      _2 = v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&bjet == 0 &&vetoes;
      _3 = v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&vetoes;
      _4 = v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&bjet == 0;
      _5 = v4l2.DeltaR(v4l3)<1 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&bjet == 0 &&vetoes;
      _6 = v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&bjet == 0 &&vetoes;
      _7 = v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && M_l2l3_combined<20 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&bjet == 0 &&vetoes;
      _8 = v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && bjet == 0 &&vetoes;
      _9 = v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && M_l2l3_combined<20 && sum_vec_l2l3>15 && prob_vertex>0.001 &&bjet == 0 &&vetoes;
      _10= v4l2.DeltaR(v4l3)<1 && M_3L_combined>50 && M_3L_combined<80 && min_delta_phi>1 && vtxRvtxPcosAlpha>0.99 && sum_vec_l2l3>15 && D2_delta_pv_svSig>20 && prob_vertex>0.001 &&bjet == 0 &&vetoes;
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     histogramm   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      double values[nDist] ={static_cast<double>(0) ,static_cast<double>(0) ,
			     static_cast<double>(0) ,
			     v4l2.Pt(),
			     v4l3.Pt(),
			     M_3L_combined,
			     M_l2l3,
			     M_l2l3,
			     M_l2l3_combined,
			     M_l2l3_combined,
			     M_3L_combined,
			     mT,
			     _met,		    
			     static_cast<double>( goodjet),
			     static_cast<double>(bjet),
			     met_sumjet,
			     fabs(_dxy[l1]),fabs(_dz[l1]),fabs(_3dIPSig[l1]), fabs(_2dIPSig[l1]), 
			     fabs(_dxy[l2]),fabs(_dz[l2]),fabs(_3dIPSig[l2]), fabs(_2dIPSig[l2]), 
			     fabs(_dxy[l3]),fabs(_dz[l3]),fabs(_3dIPSig[l3]), fabs(_2dIPSig[l3]), 
			     _relIso[l1],
			     _relIso[l2],
			     _relIso[l3],
			     v4l1.DeltaR(v4l3),
			     v4l2.DeltaR(v4l3),
			     min_delta_phi,
			     prob_vertex,
			     _vertex_normchi2,
			     _vertex_chi2,
			     vtxRvtxPcosAlpha,
			     D3_delta_pv_sv,
			     D3_delta_pv_svSig,
			     D2_delta_pv_sv,
			     D2_delta_pv_sv,
			     D2_delta_pv_svSig,
			     momentum_jet, sum_vec_l2l3};
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  filling   histogramm   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      unsigned fill = effsam;
      bool isDataDrivenBgk= false;
      if (samples[sam].isData() && tightFail_sFR     && single_fake)     isDataDrivenBgk= true;
      if (samples[sam].isData() && loose_lepton_dFR  && Double_fake)     isDataDrivenBgk= true;
      bool isDataYield= false;
      if (samples[sam].isData() && !tightFail_sFR    && single_fake)     isDataYield= true;
      if (samples[sam].isData() && tight_lepton_dFR  && Double_fake)     isDataYield= true;
      if (isDataDrivenBgk &&  single_fake && singleSF_negative) fill = nSamples_eff - 2;
      if (isDataDrivenBgk &&  single_fake && singleSF_positive) fill = nSamples_eff - 1;    
      if (isDataDrivenBgk &&  Double_fake) fill = nSamples_eff;
      if (isSignal &&  tightFail_sFR     && single_fake) continue;
      if (isSignal &&  loose_lepton_dFR  && Double_fake) continue;
 
      if (isDataYield)     fill = 0;
      if (isDataYield)     scal = 1;
      // it is blinded 
      //if (blinded && isDataYield)     continue;
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //<<<<<<<<
      //<<<<<<<<
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      // Pile UP!
      if (!samples[sam].isData()){	    
	for (int w_loop =0; w_loop < nCoupling; w_loop++){
	  weight_SR[w_loop][pu_index][0][effsam] = PUWeight();
	  weight_SR[w_loop][pu_index][1][effsam] = puWeight(1);	
	  weight_SR[w_loop][pu_index][2][effsam] = puWeight(2);	
	  if (PUWeight() == 0) std::cout<<"-----------> pileup e' zero"<<std::endl;	
	}     
      }
      //<<<<<<<<<<<<<<<<<<<<<<<< signal xsection
      if (isSignal){
	for (int w_loop =0; w_loop < nCoupling; w_loop++){
	  weight_SR[w_loop][xsecNorm_index][0][effsam] = 1.089;
	  weight_SR[w_loop][xsecNorm_index][1][effsam] = 1.089-0.041;	
	  weight_SR[w_loop][xsecNorm_index][2][effsam] = 1.089+0.042;	
	}
      }
      
      //-------------------- central values SF calculations and systematics -------------------------
      // l1   
      // Systematics on prompt muons
      if(SR_channel <= 2) {
	//ID MULTIPLIED with ISO and IP 
	weight_SR[muon_case][pMuo_index][0][effsam] = SF_prompt_muon (*&sf_prompt_muon,*&sf_iso_ip_prompt_muon, l1);
	//trigger
	weight_SR[muon_case][trigger_index][0][effsam] = SF_trigger_muon(*&sf_trigger_muon, l1);		
	//
	weight_SR[muon_case][pMuo_index][1][effsam] = SF_prompt_muon (*&sf_prompt_muon,*&sf_iso_ip_prompt_muon, l1) - SF_prompt_muon_error(*&sf_prompt_muon_syst, *&sf_prompt_muon,*&sf_iso_ip_prompt_muon,*&sf_iso_ip_prompt_muon_syst, l1);
	weight_SR[muon_case][pMuo_index][2][effsam] = SF_prompt_muon (*&sf_prompt_muon,*&sf_iso_ip_prompt_muon, l1) + SF_prompt_muon_error(*&sf_prompt_muon_syst, *&sf_prompt_muon,*&sf_iso_ip_prompt_muon,*&sf_iso_ip_prompt_muon_syst, l1);

	
	weight_SR[muon_case][trigger_index][1][effsam] = SF_trigger_muon(*&sf_trigger_muon, l1)-SF_trigger_muon_error(*&sf_trigger_muon, l1);
	weight_SR[muon_case][trigger_index][2][effsam] = SF_trigger_muon(*&sf_trigger_muon, l1)+SF_trigger_muon_error(*&sf_trigger_muon, l1);

      }
      // Systematics on prompt ele
      if(SR_channel > 2) {
	weight_SR[ele_case][pEle_index][0][effsam]  = SF_prompt_ele  (*&sf_prompt_ele, l1);
	weight_SR[ele_case][trigger_index][0][effsam] = SF_trigger_ele(*&sf_trigger_ele, l1);
	//
	weight_SR[ele_case][pEle_index][1][effsam] = SF_prompt_ele(*&sf_prompt_ele, l1)-SF_prompt_ele_error(*&sf_prompt_ele, l1);	  
	weight_SR[ele_case][pEle_index][2][effsam] = SF_prompt_ele(*&sf_prompt_ele, l1)+SF_prompt_ele_error(*&sf_prompt_ele, l1);
	weight_SR[ele_case][trigger_index][1][effsam] = SF_trigger_ele(*&sf_trigger_ele, l1)-SF_trigger_ele_error(*&sf_trigger_ele, l1);	  
	weight_SR[ele_case][trigger_index][2][effsam] = SF_trigger_ele(*&sf_trigger_ele, l1)+SF_trigger_ele_error(*&sf_trigger_ele, l1);
      }
      // Systematics on DISPALCED SIGNATURE

      //std::cout<<"===========  before calling    "<<D2_delta_pv_sv<<"   "<< (v4l2+v4l3).Pt()<<std::endl;

      double central_displaced_signature =   displaced_weight       (flavors_3l,SR_channel,_lElectronMissingHits[l2], _lElectronMissingHits[l3],displ2,displ3, (v4l2+v4l3).Pt(), D2_delta_pv_sv, displEleVars, *&sf_sv_effcy_num, *&sf_sv_effcy_den, *&sf_isoID_nPMuon, *&sf_isoID_nPMuon_syst,l2,l3);

      double variation_displaced_signature = displaced_weight_error (flavors_3l,SR_channel,_lElectronMissingHits[l2], _lElectronMissingHits[l3],displ2,displ3,  (v4l2+v4l3).Pt(), D2_delta_pv_sv, displEleVars, *&sf_sv_effcy_num, *&sf_sv_effcy_den,*&sf_isoID_nPMuon, *&sf_isoID_nPMuon_syst,l2,l3 );


      if (SR_channel== 0){
	weight_SR[muon_case][npLeptons_mm_index][0][effsam] =  central_displaced_signature;
	weight_SR[muon_case][npLeptons_mm_index][1][effsam] =  central_displaced_signature - variation_displaced_signature;
	weight_SR[muon_case][npLeptons_mm_index][2][effsam] =  central_displaced_signature + variation_displaced_signature;
      }
      if (SR_channel == 1 || SR_channel == 2){
	weight_SR[muon_case][npLeptons_em_index][0][effsam] =  central_displaced_signature;
	weight_SR[muon_case][npLeptons_em_index][1][effsam] =  central_displaced_signature - variation_displaced_signature;
	weight_SR[muon_case][npLeptons_em_index][2][effsam] =  central_displaced_signature + variation_displaced_signature;
      }
      if (SR_channel== 3){
	weight_SR[ele_case][npLeptons_ee_index][0][effsam] =  central_displaced_signature;
	weight_SR[ele_case][npLeptons_ee_index][1][effsam] =  central_displaced_signature - variation_displaced_signature;
	weight_SR[ele_case][npLeptons_ee_index][2][effsam] =  central_displaced_signature + variation_displaced_signature;
      }
      if (SR_channel == 4 || SR_channel == 5){
	weight_SR[ele_case][npLeptons_em_index][0][effsam] =  central_displaced_signature;
	weight_SR[ele_case][npLeptons_em_index][1][effsam] =  central_displaced_signature - variation_displaced_signature;
	weight_SR[ele_case][npLeptons_em_index][2][effsam] =  central_displaced_signature + variation_displaced_signature;
      }
	
      // ------------------------- DFR ------------------------- //    
      for (int w_loop =0; w_loop < nCoupling; w_loop++){
	weight_SR[w_loop][dfmm_index][0][effsam] = 1.;
	weight_SR[w_loop][dfee_index][0][effsam] = 1.;
	weight_SR[w_loop][dfem_index][0][effsam] = 1.;
	//flat 30perct on everything 
	weight_SR[w_loop][dfShape_index][0][effsam] = 1.;
	weight_SR[w_loop][dfShape_index][1][effsam] = 0.7;
	weight_SR[w_loop][dfShape_index][2][effsam] = 1.3;

	weight_SR[w_loop][dfLowStat_index][0][effsam] = 1.;
	weight_SR[w_loop][dfLowStat_index][1][effsam] = 1.;
	weight_SR[w_loop][dfLowStat_index][2][effsam] = 1.;
	
	if (M_l2l3_combined < 4 &&  D2_delta_pv_sv > 4) {
	  weight_SR[w_loop][dfLowStat_index][1][effsam] = 0.5;
	  weight_SR[w_loop][dfLowStat_index][2][effsam] = 1.5;
	}
	if (M_l2l3_combined > 4) {
	  weight_SR[w_loop][dfLowStat_index][1][effsam] = 0.5;
	  weight_SR[w_loop][dfLowStat_index][2][effsam] = 1.5;
	}
      }   
      // ------------------------- bjet SF + JEC/JER number of jets ------------------------- //        
      double btag_weight_central=1;
      double btag_weight_down=1; 	    
      double btag_weight_up=1; 	 
	  
      bjet_down_jec   = 0; 
      bjet_up_jec   = 0; 
      bjet_down_jer   = 0; 
      bjet_up_jer   = 0; 
      for (unsigned j =0; j < _nJets ; j++){
	if (jetIsBJet(j, _jetSmearedPt_JECDown[j])) ++bjet_down_jec;    
	if (jetIsBJet(j, _jetSmearedPt_JECUp[j]))   ++bjet_up_jec;     
	if (jetIsBJet(j, _jetSmearedPt_JERDown[j])) ++bjet_down_jer;     
	if (jetIsBJet(j, _jetSmearedPt_JERUp[j]))   ++bjet_up_jer;

	if(jetIsGood(j, _jetPt[j]) && _jetPt[j]<1000. && (_jetHadronFlavor[j] == 5 ||  _jetHadronFlavor[j] == 4 ||  _jetHadronFlavor[j] == 0)) {
	  double eff_cy = 0.;
	  eff_cy = SF_btag_eff(*&sf_btag_eff, _jetEta[j], _jetPt[j], _jetHadronFlavor[j]);
	  if (_jetHadronFlavor[j] == 0){	
	    btag_weight_central *=  (1. - eff_cy* reader.eval_auto_bounds("central", BTagEntry::FLAV_UDSG, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	    btag_weight_down    *=  (1. - eff_cy* reader.eval_auto_bounds("down", BTagEntry::FLAV_UDSG, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	    btag_weight_up      *=  (1. - eff_cy* reader.eval_auto_bounds("up", BTagEntry::FLAV_UDSG, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	  }
	  if (_jetHadronFlavor[j] == 4){	
	    btag_weight_central *=  (1. - eff_cy* reader.eval_auto_bounds("central", BTagEntry::FLAV_C, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	    btag_weight_down    *=  (1. - eff_cy* reader.eval_auto_bounds("down", BTagEntry::FLAV_C, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	    btag_weight_up      *=  (1. - eff_cy* reader.eval_auto_bounds("up", BTagEntry::FLAV_C, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	  }
	  if (_jetHadronFlavor[j] == 5){	
	    btag_weight_central *=  (1. - eff_cy* reader.eval_auto_bounds("central", BTagEntry::FLAV_B, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	    btag_weight_down    *=  (1. - eff_cy* reader.eval_auto_bounds("down", BTagEntry::FLAV_B, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	    btag_weight_up      *=  (1. - eff_cy* reader.eval_auto_bounds("up", BTagEntry::FLAV_B, std::abs(_jetEta[j]), _jetPt[j])) / (1. - eff_cy);
	  }		  
	}	//bjet
      }    //njet
      for (int w_loop =0; w_loop < nCoupling; w_loop++){ 
	weight_SR[w_loop][btag_index][0][effsam] = btag_weight_central;
	weight_SR[w_loop][btag_index][1][effsam] = btag_weight_down;
	weight_SR[w_loop][btag_index][2][effsam] = btag_weight_up;			
      }
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      // !!!!!!!!!!!!!!    filling all the histograms for data cards !!!!!!!!!!!!!!    	    
      double central_total_weight_mu = 1.;
      double central_total_weight_ele = 1.;

      if (!isDataDrivenBgk && !isDataYield){
	for (int w_loop =0; w_loop < nSystematic; w_loop++){
	  if (SR_channel <= 2 ) central_total_weight_mu *= weight_SR[muon_case][w_loop][0][effsam];	 
	  if (SR_channel > 2 ) central_total_weight_ele *= weight_SR[ele_case][w_loop][0][effsam];
	} 	        	
      } 
      // electron case --> eee eeµ eeµ	  
      if (SR_selection){ // only final final step 
	// central distribution --> on_index ==> 0 and  "central" => 0
	if (SR_channel <= 2 && bjet == 0)  {
	  plots_SR[muon_case][on_index][0][fill] ->  Fill(static_cast<double>(bin_SR_muonCoupling), scal*central_total_weight_mu);
	  plots_SR2[0][muon_case][on_index][0][fill] ->  Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), scal*central_total_weight_mu);
	  plots_SR2[1][muon_case][on_index][0][fill] ->  Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), scal*central_total_weight_mu);
	}
	if (SR_channel > 2  && bjet == 0) {
	  plots_SR[ele_case][on_index][0][fill]  ->  Fill(static_cast<double>(bin_SR_eleCoupling), scal*central_total_weight_ele);
	  plots_SR2[0][ele_case][on_index][0][fill]  ->  Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), scal*central_total_weight_ele);
	  plots_SR2[1][ele_case][on_index][0][fill]  ->  Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), scal*central_total_weight_ele);
	}
	// plots for systematics
	//filling the shape histogram for DF background
	if (isDataDrivenBgk && Double_fake){
	  for (int iSystematics = 1; iSystematics <  nSystematic; iSystematics++) { // loop on sys
	    if (iSystematics != dfShape_index && iSystematics != dfLowStat_index && iSystematics != dfmm_index && iSystematics != dfem_index && iSystematics != dfee_index) continue;
	    for (int iVariation = 1; iVariation < nVariation; iVariation++){//loop on up-down
	      if (SR_channel <= 2 && bjet == 0){
		plots_SR[muon_case][iSystematics][iVariation][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		plots_SR2[0][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		plots_SR2[1][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
	      }
	      if (SR_channel > 2  && bjet == 0)  {
		plots_SR[ele_case][iSystematics][iVariation][fill]  -> Fill(static_cast<double>(bin_SR_eleCoupling),  weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		plots_SR2[0][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]),  weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		plots_SR2[1][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]),  weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
	      }
	    }//end variation
	  }//end syst
	}// end isDataDrivenBgk && Double_fake
	if (!samples[sam].isData()){ // only for MC
	  for (int iSystematics = 1; iSystematics <  nSystematic; iSystematics++) { // loop on sys
	    if (iSystematics == dfShape_index) continue;
	    if (iSystematics == dfLowStat_index) continue;
	    if (iSystematics == dfmm_index) continue;
	    if (iSystematics == dfem_index) continue;
	    if (iSystematics == dfee_index) continue;
	    for (int iVariation = 1; iVariation < nVariation; iVariation++){//loop on up-down	      
	      double central_divided_by_sys_ele= 1.;
	      double central_divided_by_sys_muon= 1.;
	      if (weight_SR[muon_case][iSystematics][0][effsam]!=0 && SR_channel <= 2 )central_divided_by_sys_muon  =  central_total_weight_mu/weight_SR[muon_case][iSystematics][0][effsam];
	      if (weight_SR[ele_case][iSystematics][0][effsam] !=0 && SR_channel > 2 ) central_divided_by_sys_ele   =  central_total_weight_ele/weight_SR[ele_case][iSystematics][0][effsam];
	      if (weight_SR[muon_case][iSystematics][0][effsam]==0 && SR_channel <= 2 )central_divided_by_sys_muon  =  0.;
	      if (weight_SR[ele_case][iSystematics][0][effsam] ==0 && SR_channel > 2 ) central_divided_by_sys_ele   =  0.;
	      //
	      if (SR_channel <= 2 && weight_SR[muon_case][iSystematics][0][effsam] == 0) std::cout<<" Warning!!!! divisione per zero muon -----  "<<systNamesT[iSystematics]<< " var: "<< iVariation<<std::endl;
	      if (SR_channel > 2 && weight_SR[ele_case][iSystematics][0][effsam] == 0) std::cout<<" Warning!!!! divisione per zero ele ------  "<<systNamesT[iSystematics]<< " var: "<< iVariation<<std::endl;
	      //!JEC and !JER
	      if (iSystematics!=jec_index && iSystematics!=jer_index){
		if (SR_channel <= 2 && bjet == 0){
		  plots_SR[muon_case][iSystematics][iVariation][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		}
		if (SR_channel > 2  && bjet == 0) {
		  plots_SR[ele_case][iSystematics][iVariation][fill]  -> Fill(static_cast<double>(bin_SR_eleCoupling),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		}
	      }
	      //JEC down
	      if (iSystematics==jec_index && iVariation==1){
		if (SR_channel > 2  && bjet_down_jec == 0) {
		  plots_SR[ele_case][iSystematics][iVariation][fill]  -> Fill(static_cast<double>(bin_SR_eleCoupling),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		}
		if (SR_channel <= 2 && bjet_down_jec == 0) {
		  plots_SR[muon_case][iSystematics][iVariation][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		}
	      }	
	      if (iSystematics==jec_index && iVariation==2){
		if (SR_channel > 2  && bjet_up_jec == 0) {
		  plots_SR[ele_case][iSystematics][iVariation][fill]  -> Fill(static_cast<double>(bin_SR_eleCoupling),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		}
		if (SR_channel <= 2 && bjet_up_jec == 0) {
		  plots_SR[muon_case][iSystematics][iVariation][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		}
	      } 	    
	      if (iSystematics==jer_index && iVariation==1){
		if (SR_channel > 2  && bjet_down_jer == 0) {
		  plots_SR[ele_case][iSystematics][iVariation][fill]  -> Fill(static_cast<double>(bin_SR_eleCoupling),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		}
		if (SR_channel <= 2 && bjet_down_jer == 0) {
		  plots_SR[muon_case][iSystematics][iVariation][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		}
	      }	
	      if (iSystematics==jer_index && iVariation==2){
		if (SR_channel > 2  && bjet_up_jer == 0) {
		  plots_SR[ele_case][iSystematics][iVariation][fill]  -> Fill(static_cast<double>(bin_SR_eleCoupling),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][ele_case][iSystematics][iVariation][fill]  -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]),  central_divided_by_sys_ele* weight_SR[ele_case][iSystematics][iVariation][effsam]*scal);
		}
		if (SR_channel <= 2 && bjet_up_jer == 0) {
		  plots_SR[muon_case][iSystematics][iVariation][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[0][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(M_l2l3_combined, maxBinC_dc[0]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		  plots_SR2[1][muon_case][iSystematics][iVariation][fill] -> Fill(TMath::Min(D2_delta_pv_sv,  maxBinC_dc[1]), central_divided_by_sys_muon*weight_SR[muon_case][iSystematics][iVariation][effsam]*scal);
		}
	      }    
	    } //end loop up-down
	  } // end loop on sys
	    //
		
	    // For QCD scale uncertainties
	    /*if(bjet == 0) {
	      for(unsigned sidx=0; sidx<nQcdVars; ++sidx) {
	      double wghtCorr     = _lheWeight[qcdSystVars[sidx]-1] * (sumSimulatedEventWeights/hLheCounter->GetBinContent(qcdSystVars[sidx]));
	      double wghtCorrNorm = _lheWeight[qcdSystVars[sidx]-1];
	      if(SR_channel> 2) {
	      qcdHistos    [sidx][ ele_case][fill]->Fill(static_cast<double>(bin_SR_eleCoupling) , scal*central_total_weight_ele*wghtCorr    );
	      qcdHistosNorm[sidx][ ele_case][fill]->Fill(1.                                      , scal*central_total_weight_ele*wghtCorrNorm);
	      }
	      if(SR_channel<=2) {
	      qcdHistos    [sidx][muon_case][fill]->Fill(static_cast<double>(bin_SR_muonCoupling), scal*central_total_weight_mu *wghtCorr    );
	      qcdHistosNorm[sidx][muon_case][fill]->Fill(1.                                      , scal*central_total_weight_mu *wghtCorrNorm);
	      }
	      }
	      //
	      // For PDF uncertainties
	      for(unsigned sidx=0; sidx<nPdfVars; ++sidx) {
	      double wghtCorr     = _lheWeight[pdfSystVars[sidx]-1] * (sumSimulatedEventWeights/hLheCounter->GetBinContent(pdfSystVars[sidx]));
	      double wghtCorrNorm = _lheWeight[pdfSystVars[sidx]-1];
	      if(SR_channel> 2) {
	      pdfHistos    [sidx][ ele_case][fill]->Fill(static_cast<double>(bin_SR_eleCoupling) , scal*central_total_weight_ele*wghtCorr    );
	      pdfHistosNorm[sidx][ ele_case][fill]->Fill(1.                                      , scal*central_total_weight_ele*wghtCorrNorm);
	      }
	      if(SR_channel<=2) {
	      pdfHistos    [sidx][muon_case][fill]->Fill(static_cast<double>(bin_SR_muonCoupling), scal*central_total_weight_mu *wghtCorr    );
	      pdfHistosNorm[sidx][muon_case][fill]->Fill(1.                                      , scal*central_total_weight_mu *wghtCorrNorm);
	      }
	      }
	      } // end if(bjet == 0)
	    */
	} // end MC
      } // end SR_selection

      /*if (isDataYield) central_total_weight_mu = 1.;
      if (isDataYield) central_total_weight_ele = 1.;

      if (isDataYield && SR_selection  && bjet == 0){
	std::cout<<"event: "<<"SR_channel: "<< SR_channel<<"   --->>>   _runNb: "<<_runNb<<"   _lumiBlock:  "<<_lumiBlock<<" _eventNb:  "<<_eventNb<<"               deltaR: "<<v4l2.DeltaR(v4l3)<<"   displaced: "<<D2_delta_pv_sv<<"   mass: "<< M_l2l3_combined<<"   sumPT: "<<v4l2.Pt()+v4l3.Pt()<<" conv veto l1: "<< _lElectronPassConvVeto[l1]<< " conv veto l2: "<< _lElectronPassConvVeto[l2]<<" conv veto l3: "<< _lElectronPassConvVeto[l3]<<std::endl;

      }
      */
      

      
      if (isDataDrivenBgk && SR_selection  && bjet == 0 && single_fake){
	if (single_fake) single_fake_txt<<"event: "<<"SR_channel: "<< SR_channel<<"   "<< bin_SR_eleCoupling<<"  "<<bin_SR_muonCoupling<<" --->>>   _runNb: "<<_runNb<<"   _lumiBlock:  "<<_lumiBlock<<" _eventNb:  "<<_eventNb<<std::endl;
	if (single_fake) single_fake_txt<< " deltaR: "<<v4l2.DeltaR(v4l3)<<"   displaced: "<<D2_delta_pv_sv<<"   mass: "<< M_l2l3_combined<<"   sumPT: "<<v4l2.Pt()+v4l3.Pt()<<" conv veto l1: "<< _lElectronPassConvVeto[l1]<< " conv veto l2: "<< _lElectronPassConvVeto[l2]<<" conv veto l3: "<< _lElectronPassConvVeto[l3]<<std::endl;
	single_fake_txt<< "l2: "<< _lFlavor[l2]<<"  "<< v4l2.Pt()<< " relIso: "<< _relIso[l2]<< std::endl;
	single_fake_txt<<  "closest jet ---> px,py,pz,E  "<<_closest_l1JetPx[l2]<<" "<<_closest_l1JetPy[l2]<<" "<<_closest_l1JetPz[l2]<<" "<<_closest_l1JetE[l2]<<"   pt: "<<lepAwareJetl2[0].Pt()<< std::endl;
	single_fake_txt<< "l3: "<< _lFlavor[l3]<<"  "<< v4l3.Pt()<< " relIso: "<< _relIso[l3]<< std::endl;
	single_fake_txt<<  "closest jet ---> px,py,pz,E  "<<_closest_l1JetPx[l3]<<" "<<_closest_l1JetPy[l3]<<" "<<_closest_l1JetPz[l3]<<" "<<_closest_l1JetE[l3]<<"   pt: "<<lepAwareJetl3[0].Pt()<< std::endl;
	single_fake_txt<< "dilepton pT: "<< (v4l2+v4l3).Pt()<<std::endl;
	//single_fake_txt<<"deltaR===>  jetl2-jetl3: "<< lepAwareJetl2[0].DeltaR(lepAwareJetl3[0])<<"  jetl2-l3: "<<lepAwareJetl2[0].DeltaR(v4l3)<<"  jetl3-l2: "<<lepAwareJetl3[0].DeltaR(v4l2)<<"  jetl2-l2: "<<lepAwareJetl2[0].DeltaR(v4l2)<<"  jetl3-l3: "<<lepAwareJetl3[0].DeltaR(v4l3)<<"         l3-l2: "<<v4l3.DeltaR(v4l2)<<std::endl;


	
	if (single_fake) single_fake_txt<<"------- scal "<<scal<< " Tl2: "<<_isT[l2]<<" Tl3: "<<_isT[l3]<<std::endl;
	if (single_fake && !_isT[l2] ) single_fake_txt<<"------- scal "<<scal<< " Tl2: "<<_isT[l2]<<" Tl3: "<<_isT[l3]<<"   fr: "<< FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
																	   _lEta[l2], _lFlavor[l2], v4l2.Pt(), index_eta,flav_dRF, momentum_jet)<<std::endl;
	if (single_fake && !_isT[l3] ) single_fake_txt<<"------- scal "<<scal<< " Tl2: "<<_isT[l2]<<" Tl3: "<<_isT[l3]<<"   fr: "<< FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
																	   _lEta[l3], _lFlavor[l3], v4l3.Pt(), index_eta,flav_dRF, momentum_jet)<<std::endl;

	if (!_isT[l2] && !_isT[l3])single_fake_txt<<" both loose not tight"<<std::endl;
	if ((_isT[l2] && !_isT[l3])    ||    (!_isT[l2] && _isT[l3])) single_fake_txt<<" only 1 loose not tight"<<std::endl;
	single_fake_txt<<""<<std::endl;
	single_fake_txt<<""<<std::endl;
      }
      
     


      /*


      
      unsigned fill_2d = effsam;

      if (isDataYield) fill_2d = 0;
      if (isSignal) fill_2d = 1;
      if (isDataDrivenBgk) fill_2d=2;
      
      //D2_delta_pv_sv
      if (SR_channel <= 2) {
        if (SR_selection   && bjet == 0){
	  if (SR_channel == 0 ) Histos2d[0][0][0][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 0 ) Histos2d[1][0][0][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 0 ) Histos2d[2][0][0][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	  if (SR_channel == 0 ) Histos2d[4][0][0][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	  if (SR_channel == 0 ) Histos2d[6][0][0][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	  if (SR_channel == 0 ) Histos2d[8][0][0][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  //
	  if (SR_channel == 1 ) Histos2d[0][0][1][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 1 ) Histos2d[1][0][1][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 1 ) Histos2d[2][0][1][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	  if (SR_channel == 1 ) Histos2d[4][0][1][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	  if (SR_channel == 1 ) Histos2d[6][0][1][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	  if (SR_channel == 1 ) Histos2d[8][0][1][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  //
	  if (SR_channel == 2 ) Histos2d[0][0][2][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 2 ) Histos2d[1][0][2][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 2 ) Histos2d[2][0][2][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	  if (SR_channel == 2 ) Histos2d[4][0][2][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	  if (SR_channel == 2 ) Histos2d[6][0][2][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	  if (SR_channel == 2 ) Histos2d[8][0][2][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  if (D2_delta_pv_sv > 3.5) {
	    if (SR_channel == 0 ) Histos2d[3][0][0][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	    if (SR_channel == 0 ) Histos2d[5][0][0][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	    if (SR_channel == 0 ) Histos2d[7][0][0][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	    if (SR_channel == 0 ) Histos2d[9][0][0][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	    if (SR_channel == 1 ) Histos2d[3][0][1][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	    if (SR_channel == 1 ) Histos2d[5][0][1][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	    if (SR_channel == 1 ) Histos2d[7][0][1][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	    if (SR_channel == 1 ) Histos2d[9][0][1][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	    if (SR_channel == 2 ) Histos2d[3][0][2][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	    if (SR_channel == 2 ) Histos2d[5][0][2][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	    if (SR_channel == 2 ) Histos2d[7][0][2][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	    if (SR_channel == 2 ) Histos2d[9][0][2][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  }
	}
    

      }


      if (SR_channel > 2) {
	if (SR_selection   && bjet == 0){
	  if (SR_channel == 3 ) Histos2d[0][0][3][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 3 ) Histos2d[1][0][3][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 3 ) Histos2d[2][0][3][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	  if (SR_channel == 3 ) Histos2d[4][0][3][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	  if (SR_channel == 3 ) Histos2d[6][0][3][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	  if (SR_channel == 3 ) Histos2d[8][0][3][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  //
	  if (SR_channel == 4 ) Histos2d[0][0][4][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 4 ) Histos2d[1][0][4][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 4 ) Histos2d[2][0][4][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	  if (SR_channel == 4 ) Histos2d[4][0][4][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	  if (SR_channel == 4 ) Histos2d[6][0][4][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	  if (SR_channel == 4 ) Histos2d[8][0][4][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  //
	  if (SR_channel == 5 ) Histos2d[0][0][5][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 5 ) Histos2d[1][0][5][fill_2d] -> Fill(D2_delta_pv_sv, M_l2l3_combined,1);
	  if (SR_channel == 5 ) Histos2d[2][0][5][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	  if (SR_channel == 5 ) Histos2d[4][0][5][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	  if (SR_channel == 5 ) Histos2d[6][0][5][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	  if (SR_channel == 5 ) Histos2d[8][0][5][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  if (D2_delta_pv_sv > 3.5) {
	    if (SR_channel == 3 ) Histos2d[3][0][3][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	    if (SR_channel == 3 ) Histos2d[5][0][3][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	    if (SR_channel == 3 ) Histos2d[7][0][3][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	    if (SR_channel == 3 ) Histos2d[9][0][3][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	    if (SR_channel == 4 ) Histos2d[3][0][4][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	    if (SR_channel == 4 ) Histos2d[5][0][4][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	    if (SR_channel == 4 ) Histos2d[7][0][4][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	    if (SR_channel == 4 ) Histos2d[9][0][4][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	    if (SR_channel == 5 ) Histos2d[3][0][5][fill_2d] -> Fill(_vertex_X, _vertex_Y,1);
	    if (SR_channel == 5 ) Histos2d[5][0][5][fill_2d] -> Fill(_vertex_phi, _vertex_r,1);
	    if (SR_channel == 5 ) Histos2d[7][0][5][fill_2d] -> Fill(_vertex_eta, _vertex_r,1);
	    if (SR_channel == 5 ) Histos2d[9][0][5][fill_2d] -> Fill(_vertex_theta, _vertex_R,1);
	  }
        }

      }

*/


      

      // Fin.state  SR_channel
      // ---------------------
      // mmm        0
      // mme OS     1
      // mme SS     2
      // eee        3
      // eem OS     4
      // eem SS     5
      // ------------------- Histo SR
      if (SR_channel <= 2) {
	if (selection_0)      Histos[0][SR_channel][0][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal*central_total_weight_mu);
	if (selection_final && bjet == 0)  Histos[0][SR_channel][1][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal*central_total_weight_mu);
	if (selection_0)      Histos[0][6][0][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal*central_total_weight_mu);
	if (selection_final && bjet == 0)  Histos[0][6][1][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal*central_total_weight_mu);
      }	   	    
      if (SR_channel > 2) {
	if (selection_0)      Histos[0][SR_channel][0][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal*central_total_weight_ele);
	if (selection_final  && bjet == 0)  Histos[0][SR_channel][1][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal*central_total_weight_ele);
	if (selection_0)      Histos[0][7][0][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal*central_total_weight_ele);
	if (selection_final  && bjet == 0)  Histos[0][7][1][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal*central_total_weight_ele);
      }
      // ------------------- Histo cut flow  
      if (SR_channel <= 2){
	if (sel_0)      Histos[1][6][0][fill] -> Fill(static_cast<double>(1), scal*central_total_weight_mu);		
	if (sel_1)      Histos[1][6][0][fill] -> Fill(static_cast<double>(2), scal*central_total_weight_mu);
	if (sel_2)      Histos[1][6][0][fill] -> Fill(static_cast<double>(3), scal*central_total_weight_mu);
	if (sel_3)      Histos[1][6][0][fill] -> Fill(static_cast<double>(4), scal*central_total_weight_mu);
	if (sel_4)      Histos[1][6][0][fill] -> Fill(static_cast<double>(5), scal*central_total_weight_mu);
	if (sel_5)      Histos[1][6][0][fill] -> Fill(static_cast<double>(6), scal*central_total_weight_mu);
	if (sel_6)      Histos[1][6][0][fill] -> Fill(static_cast<double>(7), scal*central_total_weight_mu);
	if (sel_7)      Histos[1][6][0][fill] -> Fill(static_cast<double>(8), scal*central_total_weight_mu);
	if (sel_8)      Histos[1][6][0][fill] -> Fill(static_cast<double>(9), scal*central_total_weight_mu);
	if (sel_9)      Histos[1][6][0][fill] -> Fill(static_cast<double>(10), scal*central_total_weight_mu);
	if (sel_10)      Histos[1][6][0][fill] -> Fill(static_cast<double>(11), scal*central_total_weight_mu);	
      }
      if (SR_channel > 2){
	if (sel_0)      Histos[1][7][0][fill] -> Fill(static_cast<double>(1), scal*central_total_weight_ele);		
	if (sel_1)      Histos[1][7][0][fill] -> Fill(static_cast<double>(2), scal*central_total_weight_ele);
	if (sel_2)      Histos[1][7][0][fill] -> Fill(static_cast<double>(3), scal*central_total_weight_ele);
	if (sel_3)      Histos[1][7][0][fill] -> Fill(static_cast<double>(4), scal*central_total_weight_ele);
	if (sel_4)      Histos[1][7][0][fill] -> Fill(static_cast<double>(5), scal*central_total_weight_ele);
	if (sel_5)      Histos[1][7][0][fill] -> Fill(static_cast<double>(6), scal*central_total_weight_ele);
	if (sel_6)      Histos[1][7][0][fill] -> Fill(static_cast<double>(7), scal*central_total_weight_ele);
	if (sel_7)      Histos[1][7][0][fill] -> Fill(static_cast<double>(8), scal*central_total_weight_ele);
	if (sel_8)      Histos[1][7][0][fill] -> Fill(static_cast<double>(9), scal*central_total_weight_ele);
	if (sel_9)      Histos[1][7][0][fill] -> Fill(static_cast<double>(10), scal*central_total_weight_ele);
	if (sel_10)     Histos[1][7][0][fill] -> Fill(static_cast<double>(11), scal*central_total_weight_ele);
      }
      // ------------------- Histo cut flow  N-1 
      if (SR_channel <= 2){
	//Histos[1][6][0][fill]->Fill(static_cast<double>(cut_bin+1), scal);
	if (_1)      Histos[2][6][0][fill] -> Fill(static_cast<double>(1), scal*central_total_weight_mu);
	if (_2)      Histos[2][6][0][fill] -> Fill(static_cast<double>(2), scal*central_total_weight_mu);
	if (_3)      Histos[2][6][0][fill] -> Fill(static_cast<double>(3), scal*central_total_weight_mu);
	if (_4)      Histos[2][6][0][fill] -> Fill(static_cast<double>(4), scal*central_total_weight_mu);
	if (_5)      Histos[2][6][0][fill] -> Fill(static_cast<double>(5), scal*central_total_weight_mu);
	if (_6)      Histos[2][6][0][fill] -> Fill(static_cast<double>(6), scal*central_total_weight_mu);
	if (_7)      Histos[2][6][0][fill] -> Fill(static_cast<double>(7), scal*central_total_weight_mu);
	if (_8)      Histos[2][6][0][fill] -> Fill(static_cast<double>(8), scal*central_total_weight_mu);
	if (_9)      Histos[2][6][0][fill] -> Fill(static_cast<double>(9), scal*central_total_weight_mu);
	if (_10)     Histos[2][6][0][fill] -> Fill(static_cast<double>(10), scal*central_total_weight_mu);
	if (_0)      Histos[2][6][0][fill] -> Fill(static_cast<double>(11), scal*central_total_weight_mu);
      }
      if (SR_channel > 2){
	if (_1)      Histos[2][7][0][fill] -> Fill(static_cast<double>(1), scal*central_total_weight_ele);
	if (_2)      Histos[2][7][0][fill] -> Fill(static_cast<double>(2), scal*central_total_weight_ele);
	if (_3)      Histos[2][7][0][fill] -> Fill(static_cast<double>(3), scal*central_total_weight_ele);
	if (_4)      Histos[2][7][0][fill] -> Fill(static_cast<double>(4), scal*central_total_weight_ele);
	if (_5)      Histos[2][7][0][fill] -> Fill(static_cast<double>(5), scal*central_total_weight_ele);
	if (_6)      Histos[2][7][0][fill] -> Fill(static_cast<double>(6), scal*central_total_weight_ele);
	if (_7)      Histos[2][7][0][fill] -> Fill(static_cast<double>(7), scal*central_total_weight_ele);
	if (_8)      Histos[2][7][0][fill] -> Fill(static_cast<double>(8), scal*central_total_weight_ele);
	if (_9)      Histos[2][7][0][fill] -> Fill(static_cast<double>(9), scal*central_total_weight_ele);
	if (_10)     Histos[2][7][0][fill] -> Fill(static_cast<double>(10), scal*central_total_weight_ele);
	if (_0)     Histos[2][7][0][fill] -> Fill(static_cast<double>(11), scal*central_total_weight_ele);
      }         
      // ------------------- all the other histograms
      for(int numero_histo = 3; numero_histo < nDist; ++numero_histo){
	if (SR_channel <= 2){
	  if (selection_0) {
	    Histos[numero_histo][SR_channel][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_mu);
	    Histos[numero_histo][6][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_mu);
	  }
	  if (selection_final  && bjet == 0) {
	    Histos[numero_histo][SR_channel][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_mu);
	    Histos[numero_histo][6][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_mu);
	  }
	}
	if (SR_channel > 2)  {
	  if (selection_0) {
	    Histos[numero_histo][7][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_ele);
	    Histos[numero_histo][SR_channel][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_ele);
	  }
	  if (selection_final  && bjet == 0) {
	    Histos[numero_histo][7][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_ele);
	    Histos[numero_histo][SR_channel][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal*central_total_weight_ele);
	  }
	}
      }//<<<<<end histo
	
    }//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<end loop over the entries
    std::cout<<"after end loop all entries"<<std::endl;

    delete hLheCounter;
  }//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< loop over samples

  // Theory uncertainties
  //
  /*for(size_t ss=0; ss<nSamples_eff+1; ++ss) {
  // QCD uncertainties
  for(size_t ib=0; ib<nBins[0]; ++ib) {
  //for(size_t ic=0; ic<nCoupling; ++ic) {
  for(size_t ic=0; ic<2; ++ic) { // skip tau for now
  double errorByBin = 0.;
  double iniCont    = plots_SR[ic][on_index][0][ss]->GetBinContent(ib+1);
  double errorNorm  = 0.;
  double iniNorm    = plots_SR[ic][on_index][0][ss]->Integral();
  for(size_t is=0; is<nQcdVars; ++is) {
  double deltabin = iniCont>0. ? std::abs(qcdHistos[is][ic][ss]->GetBinContent(ib+1) - iniCont)/iniCont : 0.;
  if(deltabin>errorByBin) errorByBin = deltabin;
  if(ib==0) {
  double deltanorm = iniNorm>0. ? std::abs(qcdHistosNorm[is][ic][ss]->GetBinContent(ib+1) - iniNorm)/iniNorm : 0.;
  if(deltanorm>errorNorm) errorNorm = deltanorm;
  }
  }
  // Shape, down variation
  plots_SR[ic][qcdShape_index][1][ss]->SetBinContent(ib+1, iniCont/(1.+errorByBin));
  // Shape, up variation
  plots_SR[ic][qcdShape_index][2][ss]->SetBinContent(ib+1, iniCont*(1.+errorByBin));
  if(ib==0) {
  // Normalization, down variation
  plots_SR[ic][qcdNorm_index][1][ss]->Scale(1./(1.+errorNorm));
  // Normalization, up variation
  plots_SR[ic][qcdNorm_index][2][ss]->Scale(1.+errorNorm);
  }
  }
  }
  //
  // PDF uncertainties
  for(size_t ib=0; ib<nBins[0]; ++ib) {
  //for(size_t ic=0; ic<nCoupling; ++ic) {
  for(size_t ic=0; ic<2; ++ic) { // skip tau for now
  double meanByBin  = 0.;
  double errorByBin = 0.;
  double iniCont    = plots_SR[ic][on_index][0][ss]->GetBinContent(ib+1);
  double meanNorm   = 0.;
  double errorNorm  = 0.;
  //double iniNorm  = plots_SR[ic][on_index][0][ss]->Integral(1, plots_SR[ic][on_index][0][ss]->GetNbinsX());
  double iniNorm    = plots_SR[ic][on_index][0][ss]->Integral();
  for(size_t is=0; is<nPdfVars; ++is) {
  double iadd = iniCont>0. ? pdfHistos[is][ic][ss]->GetBinContent(ib+1)/iniCont : 0.;
  meanByBin += iadd;
  errorByBin += iadd*iadd;
  if(ib==0) {
  double iaddnorm = iniNorm>0. ? pdfHistosNorm[is][ic][ss]->GetBinContent(ib+1)/iniNorm : 0.;
  meanNorm += iaddnorm;
  errorNorm += iaddnorm*iaddnorm;
  }
  } // end for(size_t is=6; is<106; ++is)
    //
    // Var[x] = [1/(N-1)] * [Sum(xi^2) - (Sum(xi))^2/N]
    errorByBin = (errorByBin - (meanByBin*meanByBin/100.))/99.;
    errorByBin = std::sqrt(errorByBin);
    // Shape, down variation
    plots_SR[ic][pdfShape_index][1][ss]->SetBinContent(ib+1, iniCont/(1.+errorByBin));
    // Shape, up variation
    plots_SR[ic][pdfShape_index][2][ss]->SetBinContent(ib+1, iniCont*(1.+errorByBin));
    if(ib==0) {
    errorNorm = (errorNorm - (meanNorm*meanNorm/100.))/99.;
    errorNorm = std::sqrt(errorNorm);
    // Normalization, down variation
    plots_SR[ic][pdfNorm_index][1][ss]->Scale(1./(1.+errorNorm));
    // Normalization, up variation
    plots_SR[ic][pdfNorm_index][2][ss]->Scale(1.+errorNorm);
    }
    } // end for(size_t ic=0; ic<nCoupl; ++ic)
    } // end for(size_t ib=0; ib<nBins[0]; ++ib)
    
    }*/ // end for(size_t ss=0; ss<nSamples_eff+1; ++ss)

  // THIS IS THE UNBLIND PLOT ===>IT HAS TO SILENT IN THE PLOTTING!!!!!!!!!!!!!!!!!!
  //                |                         |
  //                V                         V
  //            unblindED           "silent" is not a verb...
  // TH1D* dataYields[nDist][nChannel][nCat];
	
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // all stack etc etc for the right plots to put in the data cards  
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<	
  for(int cha = 0; cha < nCoupling; ++cha) {	
    if (cha == 2) continue; // no taus for the moment
    for (int iSystematics = 0; iSystematics <  nSystematic; iSystematics++){// loop on sys
      for (int iVariation = 0; iVariation < nVariation; iVariation++){//loop on up-down
	for(unsigned effsam1 = 1; effsam1 < nSamples_eff +1 ; ++effsam1){
	  std::cout<<"   nSamples_eff  "<<nSamples_eff<<"  "<<eff_names[nSamples_eff]<<std::endl;
	  if (effsam1 == nSamples_eff){
	    put_at_zero(year,iSystematics,iVariation,cha, 1, *&plots_SR[cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 1, *&plots_SR2[0][cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 1, *&plots_SR2[1][cha][iSystematics][iVariation][effsam1]);
	  }
	  else if (effsam1 == nSamples_eff-1){
	    put_at_zero(year,iSystematics,iVariation,cha, 2, *&plots_SR[cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 2, *&plots_SR2[0][cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 2, *&plots_SR2[1][cha][iSystematics][iVariation][effsam1]);	  
	  }
	  else if (effsam1 == nSamples_eff-2){
	    put_at_zero(year,iSystematics,iVariation,cha, 2, *&plots_SR[cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 2, *&plots_SR2[0][cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 2, *&plots_SR2[1][cha][iSystematics][iVariation][effsam1]);	  
	  }
	  else {
	    put_at_zero(year,iSystematics,iVariation,cha, 0, *&plots_SR[cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 0, *&plots_SR2[0][cha][iSystematics][iVariation][effsam1]);
	    put_at_zero(year,iSystematics,iVariation,cha, 0, *&plots_SR2[1][cha][iSystematics][iVariation][effsam1]);
	  }	  
	}
      }	    
    }
  }
  std::cout<<"after putting at zero"<<std::endl;
  for(int cha = 0; cha < nCoupling; ++cha) {
    if (cha == 2) continue; // no taus for the moment
    for (int iSystematics = 0; iSystematics <  nSystematic; iSystematics++){// loop on sys
      for (int iVariation = 0; iVariation < nVariation; iVariation++){//loop on up-down
	sum_expected_SR[cha][iSystematics][iVariation] = (TH1D*)plots_SR[cha][iSystematics][iVariation][0]->Clone();
	sum_expected_SR2[0][cha][iSystematics][iVariation] = (TH1D*)plots_SR2[0][cha][iSystematics][iVariation][0]->Clone();
	sum_expected_SR2[1][cha][iSystematics][iVariation] = (TH1D*)plots_SR2[1][cha][iSystematics][iVariation][0]->Clone();


	/*sum_expected_SR[cha][iSystematics][iVariation] = (TH1D*)plots_SR[cha][iSystematics][iVariation][nSamples_signal+1]->Clone();
	  sum_expected_SR2[0][cha][iSystematics][iVariation] = (TH1D*)plots_SR2[0][cha][iSystematics][iVariation][nSamples_signal+1]->Clone();
	  sum_expected_SR2[1][cha][iSystematics][iVariation] = (TH1D*)plots_SR2[1][cha][iSystematics][iVariation][nSamples_signal+1]->Clone();*/

	

      }//end loop up-down
    }// end loop on sys
  }// coupling
  std::cout<<"after data"<<std::endl;

  for(int cha = 0; cha < nCoupling; ++cha) {//--
    if (cha == 2) continue; // no taus for the moment
    for (int iSystematics = 0; iSystematics <  nSystematic; iSystematics++){// loop on sys
      for (int iVariation = 0; iVariation < nVariation; iVariation++){//loop on up-down
	for(unsigned effsam1 = nSamples_signal+1; effsam1 < nSamples_eff +1 ; ++effsam1){	  
	  bkgYields_SR[cha][iSystematics][iVariation][effsam1 -nSamples_signal-1] = (TH1D*) plots_SR[cha][iSystematics][iVariation][effsam1]->Clone();
	  bkgYields_SR2[0][cha][iSystematics][iVariation][effsam1 -nSamples_signal-1] = (TH1D*) plots_SR2[0][cha][iSystematics][iVariation][effsam1]->Clone();	  
	  bkgYields_SR2[1][cha][iSystematics][iVariation][effsam1 -nSamples_signal-1] = (TH1D*) plots_SR2[1][cha][iSystematics][iVariation][effsam1]->Clone();	  	  
	 

	  /*bkgYields_SR[cha][iSystematics][iVariation][effsam1 -nSamples_signal-1] = (TH1D*) plots_SR[cha][iSystematics][iVariation][effsam1]->Clone();
	    bkgYields_SR2[0][cha][iSystematics][iVariation][effsam1 -nSamples_signal-1] = (TH1D*) plots_SR2[0][cha][iSystematics][iVariation][effsam1]->Clone();	  
	    bkgYields_SR2[1][cha][iSystematics][iVariation][effsam1 -nSamples_signal-1] = (TH1D*) plots_SR2[1][cha][iSystematics][iVariation][effsam1]->Clone();	  	  
	    if(effsam1 > nSamples_signal+1 && effsam1 <= nSamples_eff){	  
	    sum_expected_SR[cha][iSystematics][iVariation]->Add(bkgYields_SR[cha][iSystematics][iVariation][effsam1 -nSamples_signal-1]);
	    sum_expected_SR2[0][cha][iSystematics][iVariation]->Add(bkgYields_SR2[0][cha][iSystematics][iVariation][effsam1 -nSamples_signal-1]);
	    sum_expected_SR2[1][cha][iSystematics][iVariation]->Add(bkgYields_SR2[1][cha][iSystematics][iVariation][effsam1 -nSamples_signal-1]);

	    }*/

	  
	}
      }	    
    }
  }//--
  std::cout<<"after bgk"<<std::endl;

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  /* TCanvas *c_2dhisto[nhist2d][nCat2D][nFlavor][number_2d];
  std::cout<<"cavas"<<std::endl;

  for(unsigned dist = 0; dist < nhist2d; ++dist){
    for(int cat = 0; cat < nCat2D; ++cat){
      for (int flav = 0; flav < nFlavor; ++flav){
	for(unsigned effsam1 = 0; effsam1 < number_2d; ++effsam1){
	  // std::cout<<"in the loop   "<<dist<<"  "<<cat<<"  "<<flav<<" "<<effsam1<<std::endl;
	  //									  std::cout<< names2d[effsam1]<<"  "<<Hist2d_ossf[dist]<<"  "<< flavorNames[flav]<<"  "<<catNames2D[cat]<<std::endl;
	  c_2dhisto[dist][cat][flav][effsam1] = new TCanvas("c_"+ names2d[effsam1]+"_"+Hist2d_ossf[dist]+"_"+  flavorNames[flav] +"_"+ catNames2D[cat],names2d[effsam1]+ Hist2d_ossf[dist]+"_"+  flavorNames[flav] + "_"+ catNames[cat], 10,10,600,600);
	    
	  c_2dhisto[dist][cat][flav][effsam1]->cd();
	  gPad->SetRightMargin (0.15);
	  Histos2d[dist][cat][flav][effsam1] ->GetXaxis()->SetTitle("#Delta (PV-SV)_{2D}");
	  
	  if (dist == 0) Histos2d[dist][cat][flav][effsam1] ->GetYaxis()->SetTitle("M_{ll}#left(l_{2}+l_{3} #right) (GeV)");
	  if (dist == 1) Histos2d[dist][cat][flav][effsam1] ->GetYaxis()->SetTitle("M_{ll}#left(l_{2}+l_{3} #right) (GeV)");
	  if (dist == 2 || dist == 3) Histos2d[dist][cat][flav][effsam1] ->GetYaxis()->SetTitle("y (cm)");
	  if (dist == 2 || dist == 3) Histos2d[dist][cat][flav][effsam1] ->GetXaxis()->SetTitle("x (cm)");
	  if (dist == 2 || dist == 3) Histos2d[dist][cat][flav][effsam1] ->GetZaxis()->SetTitle("Events/(0.2x0.2 mm^{2})");
	  if (dist == 4 || dist == 5) Histos2d[dist][cat][flav][effsam1] ->GetYaxis()->SetTitle("r, x-y (cm)");
	  if (dist == 4 || dist == 5) Histos2d[dist][cat][flav][effsam1] ->GetXaxis()->SetTitle("#phi (rad)");
	  if (dist == 4 || dist == 5) Histos2d[dist][cat][flav][effsam1] ->GetZaxis()->SetTitle("Events/(0.008 rad x 0.1 mm)");

	  if (dist == 6 || dist == 7 || dist == 8 || dist == 9 ) Histos2d[dist][cat][flav][effsam1] ->GetYaxis()->SetTitle("r, x-y (cm)");
	  if ( dist == 8 || dist == 9 ) Histos2d[dist][cat][flav][effsam1] ->GetYaxis()->SetTitle("r, x-y-z (cm)");
	  if (dist == 6 || dist == 7) Histos2d[dist][cat][flav][effsam1] ->GetXaxis()->SetTitle("#eta");
	  if (dist == 8 || dist == 9) Histos2d[dist][cat][flav][effsam1] ->GetXaxis()->SetTitle("#theta (rad)");
	  if (dist == 6 || dist == 7) Histos2d[dist][cat][flav][effsam1] ->GetZaxis()->SetTitle("Events/(0.01 x 0.1 mm)");
	  if (dist == 8 || dist == 9) Histos2d[dist][cat][flav][effsam1] ->GetZaxis()->SetTitle("Events/(0.004 rad x 0.1 mm)");	  
	  Histos2d[dist][cat][flav][effsam1]->SetStats(0);
	  //Histos2d[dist][cat][flav][effsam1] -> SetMinimum (0.);
	  gStyle->SetPalette(kCool);
	  Histos2d[dist][cat][flav][effsam1] -> Draw("COLZ");
	  //correlation << eff_names[effsam1]+"_"+ Hist2d_ossf[dist]+"_"+ flavorNames[flav] +"_""_"+ "_"+catNames2D[cat]+ "   ----->   "<< Histos2d[dist][cat][flav][effsam1] -> GetCorrelationFactor()<< std::endl;
	  c_2dhisto[dist][cat][flav][effsam1]->Print("plot/c_"+ names2d[effsam1]+"_"+Hist2d_ossf[dist]+"_"+  flavorNames[flav] +"_"+ catNames2D[cat]+".pdf");
	  c_2dhisto[dist][cat][flav][effsam1]->Print("plotroot/c_"+ names2d[effsam1]+"_"+Hist2d_ossf[dist]+"_"+  flavorNames[flav] +"_"+ catNames2D[cat]+".root");
	}
      }
    }
    }*/
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    int numer_plot_class =0;
    numer_plot_class = nSamples_eff -  nSamples_signal;
    for(unsigned dist = 0; dist < nDist; ++dist){
      for(unsigned cat = 0; cat < nCat; ++cat){
	for(int cha = 0; cha < nChannel; ++cha){	
	  for(unsigned effsam1 = 1; effsam1 < nSamples_eff +1 ; ++effsam1){
	    // std::cout<<"in plotting one "<<std::endl;
	    if (effsam1 == nSamples_eff) put_at_zero(year,0,0,cha, 1, *&Histos[dist][cha][cat][effsam1]);
	    if (effsam1 != nSamples_eff) put_at_zero(year,0,0,cha, 0, *&Histos[dist][cha][cat][effsam1]);

	  }
	}
      }
    } 

    for(unsigned dist = 0; dist < nDist; ++dist){
      for(unsigned cat = 0; cat < nCat; ++cat){
	for(int cha = 0; cha < nChannel; ++cha){               
	  dataYields[dist][cha][cat] = (TH1D*)Histos[dist][cha][cat][0]->Clone();
	  //no due unblindng--> 	dataYields[dist][cha][cat] = (TH1D*)Histos[dist][cha][cat][nSamples_signal+1]->Clone();

	}
      }
    }
  
    for(unsigned dist = 0; dist < nDist; ++dist){//--
      for(unsigned cat = 0; cat < nCat; ++cat){
	for(int cha = 0; cha < nChannel; ++cha){	
	  for(unsigned effsam1 = nSamples_signal+1; effsam1 < nSamples_eff +1 ; ++effsam1){	  
	    bkgYields[dist][cha][cat][effsam1 -nSamples_signal-1] = (TH1D*) Histos[dist][cha][cat][effsam1]->Clone();

	    // no due unblindng 
	    //bkgYields[dist][cha][cat][effsam1 -nSamples_signal-1] = (TH1D*) Histos[dist][cha][cat][effsam1]->Clone();	  
	    // if(effsam1 > nSamples_signal+1 && effsam1 <= nSamples_eff){	  
	    //  dataYields[dist][cha][cat]->Add(bkgYields[dist][cha][cat][effsam1 -nSamples_signal-1]);
	    //   }
	  
	  }
	}
      }
    }//--
    const size_t ntab = 14;
  
    if(skipPlotting == false){
      for(unsigned dist = 0; dist < nDist; ++dist){
	for(unsigned cat = 0; cat < nCat; ++cat){
	  for(int cha = 0; cha < nChannel; ++cha){               
	    for (unsigned signal_sample = 0; signal_sample< nSamples_signal; signal_sample++){
	      signals[signal_sample] =(TH1D*)Histos[dist][cha][cat][signal_sample+1]->Clone() ;     
	    }
	    plotDataVSMC(cat,cha,dist,
			 dataYields[dist][cha][cat], bkgYields[dist][cha][cat],
			 eff_names,numer_plot_class ,
			 catNames[cat], channelNames[cha],  channelNames[cha]+"_"+ Histnames_ossf[dist]+"_"+catNames[cat],
			 true,
			 2, true, signals,  sigNames_short, nSamples_signal, false, year);
	  }//end cha
	} // end cat
      } // end histo  
    } // end if skipPlotting
   




  

    if((skipLimits && skipTables)==false) {
      ////////////////                                    ////////////////
      //// List of stuff for data cards, shape ROOT files, and tables ////
      ////////////////                                    ////////////////
      //
      // List of couplings
      //const std::string couplings[] = {"ele", "muo", "tau"};
      const std::string couplings[] = {"muo", "ele"};
      const size_t nCoupl = sizeof(couplings)/sizeof(couplings[0]);

      // List of backgrounds
      //const std::string bkgNames[] = {"DY", "ttbar", "WJets", "multiboson", "Xgamma", "TTTX", "nonpromptSF", "nonpromptDF"};
      //const std::string bkgNames[] = {"DY", "multiboson", "Xgamma", "TTTX", "nonpromptSF", "nonpromptDF"};

      const std::string bkgNames[] = {"Xgamma", "nonpromptSFneg", "nonpromptSFpos","nonpromptDF"};
      const size_t nBkg = sizeof(bkgNames)/sizeof(bkgNames[0]);

      // Output directory for datacards and shape ROOT files
      std::string datacarddir = "dataCards_shapeRoot";

      // List of signal and background labels (for tables)
      std::map<std::string, std::string> labelPerProc;
      labelPerProc["signal"     ] = "signal"; // to be changed...
      //labelPerProc["DY"         ] = "$\\PZ\\rarr\\lept\\lept$";
      //labelPerProc["ttbar"      ] = "Top";
      // labelPerProc["WJets"      ] = "$\\PW +$ jets";
      //labelPerProc["multiboson" ] = "Multiboson";
      labelPerProc["Xgamma"     ] = "X $+ \\gamma$";
      //labelPerProc["TTTX"       ] = "Top $+$ X";
      labelPerProc["nonpromptSFneg"] = "Nonprompt SF neg";
      labelPerProc["nonpromptSFpos"] = "Nonprompt SF pos";
      labelPerProc["nonpromptDF"] = "Nonprompt DF";

      // List of systematics
      //const TString systNamesT[nSystematic] 	= { "on", "pu", "qcdNorm", "qcdShape", "pdfNorm", "pdfShape", "pEle", "pMuo", "npLeptons", "jec", "jer", "btag", "trigger","dfShape"};
 
      const std::string systNames[] = {"n", "pu", "xsecNorm", "pEle", "pMuo", "npLeptons_mm","npLeptons_em","npLeptons_ee", "jec", "jer", "btag", "trigger","dfShape","dfLowStat","dfmm","dfem","dfee","lumi", "npsfnorm"};
      const size_t nSyst = sizeof(systNames)/sizeof(systNames[0]) - 1;

      // List of systematics applicable to each process (signal + backgrounds)
      // (delete sample name from list if the systematic source does not aplpy to it)
      std::map<std::string, std::string> procPerSyst;
      //                       Type     Correl.   Processes
      //                       -------  --------  -------------------------------------------------------------
      procPerSyst["pu"      ] = "shapeN; not_corr; signal,  Xgamma                           ";
      procPerSyst["xsecNorm" ] = "shapeN;  is_corr; signal                           ";
      // procPerSyst["qcdShape"] = "shapeN;  is_corr; signal,  Xgamma                           ";
      //procPerSyst["pdfNorm" ] = "shapeN;  is_corr; signal,  Xgamma                           ";
      // procPerSyst["pdfShape"] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["pEle"    ] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["pMuo"    ] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["npLeptons_mm"] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["npLeptons_em"] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["npLeptons_ee"] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["jec"     ] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["jer"     ] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["btag"    ] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["trigger" ] = "shapeN;  is_corr; signal,  Xgamma                           ";
      procPerSyst["dfShape" ] = "shapeN;  is_corr;                                                                  nonpromptDF";
      procPerSyst["dfLowStat" ] = "shapeN;  is_corr;                                                                  nonpromptDF";
      procPerSyst["dfmm" ] = "shapeN;  not_corr;                                                                  nonpromptDF";
      procPerSyst["dfem" ] = "shapeN;  not_corr;                                                                  nonpromptDF";
      procPerSyst["dfee" ] = "shapeN;  not_corr;                                                                  nonpromptDF";
      procPerSyst["lumi"    ] = "lnN   ; not_corr; signal,  Xgamma                          ";
      procPerSyst["npsfnorm"] = "lnN   ;  is_corr;                                            nonpromptSFneg,  nonpromptSFpos           ";


      /* procPerSyst["pu"      ] = "shapeN; not_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["qcdNorm" ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["qcdShape"] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["pdfNorm" ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["pdfShape"] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["pEle"    ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["pMuo"    ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["npLeptons"] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["jec"     ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["jer"     ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["btag"    ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["trigger" ] = "shapeN;  is_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["dfShape" ] = "shapeN;  is_corr;                                                                  nonpromptDF";
	 procPerSyst["dfLowStat" ] = "shapeN;  is_corr;                                                                  nonpromptDF";
	 procPerSyst["dfmm" ] = "shapeN;  not_corr;                                                                  nonpromptDF";
	 procPerSyst["dfem" ] = "shapeN;  not_corr;                                                                  nonpromptDF";
	 procPerSyst["dfee" ] = "shapeN;  not_corr;                                                                  nonpromptDF";
	 procPerSyst["lumi"    ] = "lnN   ; not_corr; signal, DY,  multiboson, Xgamma, TTTX                          ";
	 procPerSyst["npsfnorm"] = "lnN   ;  is_corr;                                                     nonpromptSF             ";
      */


      std::map<std::string, std::vector<std::string> > normSystsPerYear;
      normSystsPerYear["lumi"    ] = {"1.025", "1.023", "1.025"};
      normSystsPerYear["npsfnorm"] = {"1.400", "1.400", "1.400"};
      //normSystsPerYear["npdfnorm"] = {"1.400", "1.400", "1.400"};

      //if(systcat==0) { // print data card only if systcat==0
      // Size of tab

      for(size_t isign=0; isign<nSamples_signal; ++isign) {
	std::string sgn = sigNames[isign].Data();
	for(size_t icoup=0; icoup<nCoupling; ++icoup) {
	  if(icoup==2) continue;     
	  if(icoup==0 && sgn.find("_mu")==std::string::npos) continue;
	  if(icoup==1 && sgn.find("_e" )==std::string::npos) continue;
	  std::string cpl = couplings[icoup];
	  //
	  // ========================================================
	  //   Write shape ROOT files and data cards (NEED TO CHECK AFTER LAST CHANGES!!!!)
	  // ========================================================
	  //
	  if(skipLimits==false) {
	    std::string appx1 = "";  
	    if(year==0) appx1 += "_16";
	    if(year==1) appx1 += "_17";
	    if(year==2) appx1 += "_18";
      	


	  
	    // ROOT file with shapes
	    std::string rootfilename = outfilename+"_"+sgn+"_"+cpl+appx1+".root";
	    TFile *rootfile = new TFile((datacarddir+"/"+rootfilename).c_str(), "RECREATE");
	    rootfile->cd();
	    sum_expected_SR[icoup][0][0]-> Write ("data_obs");
	    //sum_observed_SR[icoup][0][0]-> Write ("data_obs");      
	    //dataYields[0][couplidx[icoup]][6]->Write("data_obs"); 
		      
	    plots_SR[icoup][0][0][1+isign] ->Write("signal");   
	    //Histos[0][couplidx[icoup]][6][1+isign]->Write("signal");

	    // Stream for writing card and tables
	    std::ofstream card;

	    // Add .txt to name if no file extension is given
	    std::string cardName = datacarddir+"/"+sgn+"_"+cpl+"_datacard" +appx1+".txt";
	    card.open(cardName + ((cardName.find(".txt") == std::string::npos) ? ".txt" : ""));
	    // Define number of channels, background sources and systematics
	    card << "imax 1 number of channels\n";
	    card << "jmax " << nBkg  << " number of backgrounds\n";
	    card << "kmax " << nSyst << " number of nuisance parameters (sources of systematical uncertainties)\n";
	    card << "----------------------------------------------------------------------------------------\n";

	    // Shape file
	    card << "shapes * * " << rootfilename.c_str() << " $PROCESS $PROCESS_$SYSTEMATIC\n";
	    card << "----------------------------------------------------------------------------------------\n";

	    // Define the channels and the number of observed events
	    card << "bin bin1\n";
	    // While we are blinded, dataYields[0][couplidx[icoup]][6] is filled with sum of backgrounds
	    card << "observation " << std::fixed << std::setprecision(7) << sum_expected_SR[icoup][0][0]->Integral(0, -1) << "\n";
	    // Define all backgrounds and their yields
	    card << left << std::setw(2*ntab) << "bin";
	    for(unsigned proc=0; proc<nBkg+1; ++proc) {
	      card << left << std::setw(ntab) << "bin1";
	    }
	    card << "\n";
	    card << left << std::setw(2*ntab) << "process";
	    card << left << std::setw(ntab)   << "signal";
	    for(unsigned bkg=0; bkg<nBkg; ++bkg) {
	      card << left << std::setw(ntab) << bkgNames[bkg];
	    }
	    card << "\n";
	    card << left << std::setw(2*ntab) << "process";
	    for(unsigned bkg=0; bkg<nBkg+1; ++bkg){
	      card << left << std::setw(ntab) << bkg;
	    }
	    card << "\n";
	    card << left << std::setw(2*ntab) << "rate";
	    card << left << std::setw(ntab)   << std::setprecision(5) << plots_SR[icoup][0][0][1+isign]->Integral(0, -1);


	    std::cout<<"check: cpupling "<<icoup<<" .  signal: "<<sgn<<std::endl;
	    std::cout<<"sum: "<< sum_expected_SR[icoup][0][0]->Integral(0, -1) <<" .  "<< sum_expected_SR[icoup][0][0]->GetNbinsX()<<std::endl;

	    for(unsigned bkg=0; bkg<nBkg; ++bkg) {
	      rootfile->cd();
	      plots_SR[icoup][0][0][1+nSamples_signal+bkg] -> Write(bkgNames[bkg].c_str());
	      float iyield = plots_SR[icoup][0][0][1+nSamples_signal+bkg]->Integral(0, -1);
	  
	      if(iyield<=0) card << left << std::setw(ntab) << "0.0000000";
	      else          card << left << std::setw(ntab) << std::setprecision(7) << iyield;
	    }
	    card << "\n";
	    card << "----------------------------------------------------------------------------------------\n";

	    // Define sources of systematic uncertainty, what distibution they follow and how large their effect is
	    for(unsigned syst=1; syst<=nSyst; ++syst) {
	      std::string asyst = systNames[syst];
	      if(procPerSyst.count(asyst)==0) {
		std::cout << " >>> WARNING: systematic source " << asyst << " not found in the list procPerSyst! <<<" << std::endl;
		continue;
	      }
	      // Correlated or uncorrelated
	      if(procPerSyst[asyst].find("not_corr")!=std::string::npos) {
		asyst += (year==0 ? "_16" : (year==1 ? "_17" : "_18"));
	      }

	      card << left << std::setw(ntab) << asyst;
	      // If shape error, set it to 1.000
	      std::string errStr = "1.000";
	      // If normalization error, change it accordingly
	      if(procPerSyst[systNames[syst]].find("lnN")!=std::string::npos) { // normalization error: lnN
		if(normSystsPerYear.count(systNames[syst])==0) {
		  std::cout << " >>> WARNING: normalization systematic uncertainty " << asyst << " not found in the list normSystsPerYear! Set it to 100%! <<<" << std::endl;
		  errStr = "2.000";
		}
		else {
		  errStr = normSystsPerYear[systNames[syst]][year];
		}
		card << left << std::setw(ntab) << "lnN";
	      }
	      else { // all the other systematics: shapeN
		card << left << std::setw(ntab) << "shapeN";
	      }
	      //
	      // Fill in systs for all processes:
	      //
	      //  - signal
	      if(procPerSyst[systNames[syst]].find("signal")==std::string::npos)
		card << left << std::setw(ntab) << "-";
	      else
		card << left << std::setw(ntab) << errStr.c_str();
	      //
	      //  - backgrounds
	      for(unsigned bkg=0; bkg<nBkg; ++bkg) {
		if(procPerSyst[systNames[syst]].find(bkgNames[bkg])==std::string::npos)
		  card << left << std::setw(ntab) << "-";
		else
		  card << left << std::setw(ntab) << errStr;
	      }
	      card << "\n";
	    } // end systs

	    card << "* autoMCStats 0.5\n";
	    card.close();
	    rootfile->Close();


	    //++++++++++++++++++++++++++++++++++++++++ for ml2l3

	    // ROOT file with shapes
	    std::string rootfilename1 = outfilename+"_"+sgn+"_"+cpl+"_mass"+appx1+".root";
	    TFile *rootfile1 = new TFile((datacarddir+"/"+rootfilename1).c_str(), "RECREATE");
	    rootfile1->cd();
	    sum_expected_SR2[0][icoup][0][0]-> Write ("data_obs");
	    //sum_observed_SR[icoup][0][0]-> Write ("data_obs");
	    //dataYields[0][couplidx[icoup]][6]->Write("data_obs");

	    plots_SR2[0][icoup][0][0][1+isign] ->Write("signal");
	    //Histos[0][couplidx[icoup]][6][1+isign]->Write("signal");

	    // Stream for writing card1 and tables
	    std::ofstream card1;

	    // Add .txt to name if no file extension is given
	    std::string cardName1 = datacarddir+"/"+sgn+"_"+cpl+"_mass"+"_datacard" +appx1+".txt";
	    card1.open(cardName1 + ((cardName1.find(".txt") == std::string::npos) ? ".txt" : ""));
	    // Define number of channels, background sources and systematics
	    card1 << "imax 1 number of channels\n";
	    card1 << "jmax " << nBkg  << " number of backgrounds\n";
	    card1 << "kmax " << nSyst << " number of nuisance parameters (sources of systematical uncertainties)\n";
	    card1 << "----------------------------------------------------------------------------------------\n";

	    // Shape file
	    card1 << "shapes * * " << rootfilename1.c_str() << " $PROCESS $PROCESS_$SYSTEMATIC\n";
	    card1 << "----------------------------------------------------------------------------------------\n";

	    // Define the channels and the number of observed events
	    card1 << "bin bin1\n";
	    // While we are blinded, dataYields[0][couplidx[icoup]][6] is filled with sum of backgrounds
	    card1 << "observation " << std::fixed << std::setprecision(7) << sum_expected_SR2[0][icoup][0][0]->Integral(0, -1) << "\n";
	    // Define all backgrounds and their yields
	    card1 << left << std::setw(2*ntab) << "bin";
	    for(unsigned proc=0; proc<nBkg+1; ++proc) {
	      card1 << left << std::setw(ntab) << "bin1";
	    }
	    card1 << "\n";
	    card1 << left << std::setw(2*ntab) << "process";
	    card1 << left << std::setw(ntab)   << "signal";
	    for(unsigned bkg=0; bkg<nBkg; ++bkg) {
	      card1 << left << std::setw(ntab) << bkgNames[bkg];
	    }
	    card1 << "\n";
	    card1 << left << std::setw(2*ntab) << "process";
	    for(unsigned bkg=0; bkg<nBkg+1; ++bkg){
	      card1 << left << std::setw(ntab) << bkg;
	    }
	    card1 << "\n";
	    card1 << left << std::setw(2*ntab) << "rate";
	    card1 << left << std::setw(ntab)   << std::setprecision(5) << plots_SR2[0][icoup][0][0][1+isign]->Integral(0, -1);


	    std::cout<<"check: cpupling "<<icoup<<" .  signal: "<<sgn<<std::endl;
	    std::cout<<"sum: "<< sum_expected_SR2[0][icoup][0][0]->Integral(0, -1) <<" .  "<< sum_expected_SR2[0][icoup][0][0]->GetNbinsX()<<std::endl;

	    for(unsigned bkg=0; bkg<nBkg; ++bkg) {
	      rootfile1->cd();
	      plots_SR2[0][icoup][0][0][1+nSamples_signal+bkg] -> Write(bkgNames[bkg].c_str());
	      float iyield = plots_SR2[0][icoup][0][0][1+nSamples_signal+bkg]->Integral(0, -1);
    
	      if(iyield<=0) card1 << left << std::setw(ntab) << "0.0000000";
	      else          card1 << left << std::setw(ntab) << std::setprecision(7) << iyield;
	    }
	    card1 << "\n";
	    card1 << "----------------------------------------------------------------------------------------\n";

	    // Define sources of systematic uncertainty, what distibution they follow and how large their effect is
	    for(unsigned syst=1; syst<=nSyst; ++syst) {
	      std::string asyst = systNames[syst];
	      if(procPerSyst.count(asyst)==0) {
		std::cout << " >>> WARNING: systematic source " << asyst << " not found in the list procPerSyst! <<<" << std::endl;
		continue;
	      }
	      // Correlated or uncorrelated
	      if(procPerSyst[asyst].find("not_corr")!=std::string::npos) {
		asyst += (year==0 ? "_16" : (year==1 ? "_17" : "_18"));
	      }
    
	      card1 << left << std::setw(ntab) << asyst;
	      // If shape error, set it to 1.000
	      std::string errStr = "1.000";
	      // If normalization error, change it accordingly
	      if(procPerSyst[systNames[syst]].find("lnN")!=std::string::npos) { // normalization error: lnN
		if(normSystsPerYear.count(systNames[syst])==0) {
		  std::cout << " >>> WARNING: normalization systematic uncertainty " << asyst << " not found in the list normSystsPerYear! Set it to 100%! <<<" << std::endl;
		  errStr = "2.000";
		}
		else {
		  errStr = normSystsPerYear[systNames[syst]][year];
		}
		card1 << left << std::setw(ntab) << "lnN";
	      }
	      else { // all the other systematics: shapeN
		card1 << left << std::setw(ntab) << "shapeN";
	      }
	      //
	      // Fill in systs for all processes:
	      //
	      //  - signal
	      if(procPerSyst[systNames[syst]].find("signal")==std::string::npos)
		card1 << left << std::setw(ntab) << "-";
	      else
		card1 << left << std::setw(ntab) << errStr.c_str();
	      //
	      //  - backgrounds
	      for(unsigned bkg=0; bkg<nBkg; ++bkg) {
		if(procPerSyst[systNames[syst]].find(bkgNames[bkg])==std::string::npos)
		  card1 << left << std::setw(ntab) << "-";
		else
		  card1 << left << std::setw(ntab) << errStr;
	      }
	      card1 << "\n";
	    } // end systs

	    card1 << "* autoMCStats 0.5\n";
	    card1.close();
	    rootfile1->Close();

	    //++++++++++++++++++++++++++++++++++++++++ for delta2d

	    // ROOT file with shapes
	    std::string rootfilename2 = outfilename+"_"+sgn+"_"+cpl+"_disp"+appx1+".root";
	    TFile *rootfile2 = new TFile((datacarddir+"/"+rootfilename2).c_str(), "RECREATE");
	    rootfile2->cd();
	    sum_expected_SR2[1][icoup][0][0]-> Write ("data_obs");
	    //sum_observed_SR[icoup][0][0]-> Write ("data_obs");
	    //dataYields[0][couplidx[icoup]][6]->Write("data_obs");

	    plots_SR2[1][icoup][0][0][1+isign] ->Write("signal");
	    //Histos[0][couplidx[icoup]][6][1+isign]->Write("signal");

	    // Stream for writing card2 and tables
	    std::ofstream card2;

	    // Add .txt to name if no file extension is given
	    std::string cardName2 = datacarddir+"/"+sgn+"_"+cpl+"_disp"+"_datacard" +appx1+".txt";
	    card2.open(cardName2 + ((cardName2.find(".txt") == std::string::npos) ? ".txt" : ""));
	    // Define number of channels, background sources and systematics
	    card2 << "imax 1 number of channels\n";
	    card2 << "jmax " << nBkg  << " number of backgrounds\n";
	    card2 << "kmax " << nSyst << " number of nuisance parameters (sources of systematical uncertainties)\n";
	    card2 << "----------------------------------------------------------------------------------------\n";

	    // Shape file
	    card2 << "shapes * * " << rootfilename2.c_str() << " $PROCESS $PROCESS_$SYSTEMATIC\n";
	    card2 << "----------------------------------------------------------------------------------------\n";

	    // Define the channels and the number of observed events
	    card2 << "bin bin1\n";
	    // While we are blinded, dataYields[0][couplidx[icoup]][6] is filled with sum of backgrounds
	    card2 << "observation " << std::fixed << std::setprecision(7) << sum_expected_SR2[1][icoup][0][0]->Integral(0, -1) << "\n";
	    // Define all backgrounds and their yields
	    card2 << left << std::setw(2*ntab) << "bin";
	    for(unsigned proc=0; proc<nBkg+1; ++proc) {
	      card2 << left << std::setw(ntab) << "bin1";
	    }
	    card2 << "\n";
	    card2 << left << std::setw(2*ntab) << "process";
	    card2 << left << std::setw(ntab)   << "signal";
	    for(unsigned bkg=0; bkg<nBkg; ++bkg) {
	      card2 << left << std::setw(ntab) << bkgNames[bkg];
	    }
	    card2 << "\n";
	    card2 << left << std::setw(2*ntab) << "process";
	    for(unsigned bkg=0; bkg<nBkg+1; ++bkg){
	      card2 << left << std::setw(ntab) << bkg;
	    }
	    card2 << "\n";
	    card2 << left << std::setw(2*ntab) << "rate";
	    card2 << left << std::setw(ntab)   << std::setprecision(5) << plots_SR2[1][icoup][0][0][1+isign]->Integral(0, -1);


	    std::cout<<"check: cpupling "<<icoup<<" .  signal: "<<sgn<<std::endl;
	    std::cout<<"sum: "<< sum_expected_SR2[1][icoup][0][0]->Integral(0, -1) <<" .  "<< sum_expected_SR2[1][icoup][0][0]->GetNbinsX()<<std::endl;

	    for(unsigned bkg=0; bkg<nBkg; ++bkg) {
	      rootfile2->cd();
	      plots_SR2[1][icoup][0][0][1+nSamples_signal+bkg] -> Write(bkgNames[bkg].c_str());
	      float iyield = plots_SR2[1][icoup][0][0][1+nSamples_signal+bkg]->Integral(0, -1);
    
	      if(iyield<=0) card2 << left << std::setw(ntab) << "0.0000000";
	      else          card2 << left << std::setw(ntab) << std::setprecision(7) << iyield;
	    }
	    card2 << "\n";
	    card2 << "----------------------------------------------------------------------------------------\n";

	    // Define sources of systematic uncertainty, what distibution they follow and how large their effect is
	    for(unsigned syst=1; syst<=nSyst; ++syst) {
	      std::string asyst = systNames[syst];
	      if(procPerSyst.count(asyst)==0) {
		std::cout << " >>> WARNING: systematic source " << asyst << " not found in the list procPerSyst! <<<" << std::endl;
		continue;
	      }
	      // Correlated or uncorrelated
	      if(procPerSyst[asyst].find("not_corr")!=std::string::npos) {
		asyst += (year==0 ? "_16" : (year==1 ? "_17" : "_18"));
	      }
    
	      card2 << left << std::setw(ntab) << asyst;
	      // If shape error, set it to 1.000
	      std::string errStr = "1.000";
	      // If normalization error, change it accordingly
	      if(procPerSyst[systNames[syst]].find("lnN")!=std::string::npos) { // normalization error: lnN
		if(normSystsPerYear.count(systNames[syst])==0) {
		  std::cout << " >>> WARNING: normalization systematic uncertainty " << asyst << " not found in the list normSystsPerYear! Set it to 100%! <<<" << std::endl;
		  errStr = "2.000";
		}
		else {
		  errStr = normSystsPerYear[systNames[syst]][year];
		}
		card2 << left << std::setw(ntab) << "lnN";
	      }
	      else { // all the other systematics: shapeN
		card2 << left << std::setw(ntab) << "shapeN";
	      }
	      //
	      // Fill in systs for all processes:
	      //
	      //  - signal
	      if(procPerSyst[systNames[syst]].find("signal")==std::string::npos)
		card2 << left << std::setw(ntab) << "-";
	      else
		card2 << left << std::setw(ntab) << errStr.c_str();
	      //
	      //  - backgrounds
	      for(unsigned bkg=0; bkg<nBkg; ++bkg) {
		if(procPerSyst[systNames[syst]].find(bkgNames[bkg])==std::string::npos)
		  card2 << left << std::setw(ntab) << "-";
		else
		  card2 << left << std::setw(ntab) << errStr;
	      }
	      card2 << "\n";
	    } // end systs

	    card2 << "* autoMCStats 0.5\n";
	    card2.close();
	    rootfile2->Close();

	  
	  } // end if(skipLimits)

	  // if(skipTables==false) {
	  //   tabletexS.close();
	  //   tabletexL.close();
	  // }
	} // end for(size_t icoup=0; icoup<nCoupling; ++icoup)
      } // end for(size_t isign=0; isign<nSamples_signal; ++isign)
      //  } // end if(systcat==0)

      //else { // if(systcat!=0)
      if(skipLimits==false) {
	std::string appx1 = "";  
	if(year==0) appx1 += "_16";
	if(year==1) appx1 += "_17";
	if(year==2) appx1 += "_18";
      	

      
	for(unsigned syst=1; syst<=nSyst; ++syst) {
	  if(procPerSyst[systNames[syst]].find("lnN")!=std::string::npos) continue;
	  for (unsigned iVariation = 1; iVariation < nVariation; iVariation++) { //loop on up-down
	    std::string appx = "_" + systNames[syst];
	    if(procPerSyst[systNames[syst]].find("not_corr")!=std::string::npos) {
	      if(year==0) appx += "_16";
	      if(year==1) appx += "_17";
	      if(year==2) appx += "_18";
	    }	     	
	    appx += (iVariation==1 ? "Down" : "Up");    
	    
	    for(size_t isign=0; isign<nSamples_signal; ++isign) {
	      std::string sgn = sigNames[isign].Data();
	      for(size_t icoup=0; icoup<nCoupl; ++icoup) {
		if(icoup==0 && sgn.find("_mu" )==std::string::npos) continue;
		if(icoup==1 && sgn.find("_e")==std::string::npos) continue;
		std::string cpl = couplings[icoup];

		// ROOT file with shapes
		std::string rootfilename = outfilename+"_"+sgn+"_"+cpl+appx1+".root";
		TFile *rootfile = TFile::Open((datacarddir+"/"+rootfilename).c_str(), "UPDATE");
		rootfile->cd();	    
		plots_SR[icoup][syst][iVariation][1+isign] ->Write(("signal"+appx).c_str());

	  
		for(unsigned bkg=0; bkg<nBkg; ++bkg) {
		  if (bkg == nBkg -2) continue;
		  rootfile->cd(); 	
		  plots_SR[icoup][syst][iVariation][1+nSamples_signal+bkg]->Write((bkgNames[bkg]+appx).c_str());
		}
		rootfile->Close();
	      } // end couplings
	    } // end signal samples
	  } // variation up down
	} // loop sty


	//++++++++++++++++++++++++++++++++++++++++ for ml2l3

	for(unsigned syst=1; syst<=nSyst; ++syst) {
	  if(procPerSyst[systNames[syst]].find("lnN")!=std::string::npos) continue;
	  for (unsigned iVariation = 1; iVariation < nVariation; iVariation++) { //loop on up-down
	    std::string appx = "_" + systNames[syst];
	    if(procPerSyst[systNames[syst]].find("not_corr")!=std::string::npos) {
	      if(year==0) appx += "_16";
	      if(year==1) appx += "_17";
	      if(year==2) appx += "_18";
	    }
	    appx += (iVariation==1 ? "Down" : "Up");
            
	    for(size_t isign=0; isign<nSamples_signal; ++isign) {
	      std::string sgn = sigNames[isign].Data();
	      for(size_t icoup=0; icoup<nCoupl; ++icoup) {
		if(icoup==0 && sgn.find("_mu" )==std::string::npos) continue;
		if(icoup==1 && sgn.find("_e")==std::string::npos) continue;
		std::string cpl = couplings[icoup];
                    
		// ROOT file with shapes
		std::string rootfilename1 = outfilename+"_"+sgn+"_"+cpl+"_mass"+appx1+".root";
		TFile *rootfile1 = TFile::Open((datacarddir+"/"+rootfilename1).c_str(), "UPDATE");
		rootfile1->cd();
		plots_SR2[0][icoup][syst][iVariation][1+isign] ->Write(("signal"+appx).c_str());
                    
                    
		for(unsigned bkg=0; bkg<nBkg; ++bkg) {
		  if (bkg == nBkg -2) continue;
		  rootfile1->cd();
		  plots_SR2[0][icoup][syst][iVariation][1+nSamples_signal+bkg]->Write((bkgNames[bkg]+appx).c_str());
		}
		rootfile1->Close();
	      } // end couplings
	    } // end signal samples
	  } // variation up down
	} // loop sty

	//++++++++++++++++++++++++++++++++++++++++ for delta2d
	for(unsigned syst=1; syst<=nSyst; ++syst) {
	  if(procPerSyst[systNames[syst]].find("lnN")!=std::string::npos) continue;
	  for (unsigned iVariation = 1; iVariation < nVariation; iVariation++) { //loop on up-down
	    std::string appx = "_" + systNames[syst];
	    if(procPerSyst[systNames[syst]].find("not_corr")!=std::string::npos) {
	      if(year==0) appx += "_16";
	      if(year==1) appx += "_17";
	      if(year==2) appx += "_18";
	    }
	    appx += (iVariation==1 ? "Down" : "Up");
            
	    for(size_t isign=0; isign<nSamples_signal; ++isign) {
	      std::string sgn = sigNames[isign].Data();
	      for(size_t icoup=0; icoup<nCoupl; ++icoup) {
		if(icoup==0 && sgn.find("_mu" )==std::string::npos) continue;
		if(icoup==1 && sgn.find("_e")==std::string::npos) continue;
		std::string cpl = couplings[icoup];
                    
		// ROOT file with shapes
		std::string rootfilename2 = outfilename+"_"+sgn+"_"+cpl+"_disp"+appx1+".root";
		TFile *rootfile2 = TFile::Open((datacarddir+"/"+rootfilename2).c_str(), "UPDATE");
		rootfile2->cd();
		plots_SR2[1][icoup][syst][iVariation][1+isign] ->Write(("signal"+appx).c_str());
                    
                    
		for(unsigned bkg=0; bkg<nBkg; ++bkg) {
		  if (bkg == nBkg -2) continue;
		  rootfile2->cd();
		  plots_SR2[1][icoup][syst][iVariation][1+nSamples_signal+bkg]->Write((bkgNames[bkg]+appx).c_str());
		}
		rootfile2->Close();
	      } // end couplings
	    } // end signal samples
	  } // variation up down
	} // loop sty


 
      } // end if(skipLimits==false)
      // } // end if(systcat!=0)

    } // end of if((skipLimits && skipTables)==false)

    std::cout << "  === This is the end of all hope ===" << std::endl;

} // end analisi



  //___________________________________________________________________
  double Analysis_mc::pu_weight ( TH1D *histo, double numberInteractions){
    double nI = numberInteractions;   
    double factore=0;
    factore = histo->GetBinContent(histo->FindBin(nI));
    return factore;
  }

 
double Analysis_mc::displMuoVars(double idispl, double ipt) {
  double ieff = 1.0;
  // 2016
  if(year==0) {
    if     (idispl<0.2) {
      if     (ipt< 6.) ieff = 0.995;
      else if(ipt<10.) ieff = 0.995;
      else if(ipt<20.) ieff = 1.000;
      else             ieff = 0.987;
    }
    else if(idispl<0.5) {
      if     (ipt< 6.) ieff = 1.005;
      else if(ipt<10.) ieff = 1.002;
      else if(ipt<20.) ieff = 1.002;
      else             ieff = 0.991;
    }
    else if(idispl<1.0) {
      if     (ipt< 6.) ieff = 1.018;
      else if(ipt<10.) ieff = 1.007;
      else if(ipt<20.) ieff = 0.985;
      else             ieff = 1.013;
    }
    else {
      if     (ipt< 6.) ieff = 1.008;
      else if(ipt<10.) ieff = 1.021;
      else if(ipt<20.) ieff = 0.976;
      else             ieff = 1.012;
    }
  }
  // 2017
  else if(year==1) {
    if     (idispl<0.2) {
      if     (ipt< 6.) ieff = 0.995;
      else if(ipt<10.) ieff = 0.995;
      else if(ipt<20.) ieff = 1.000;
      else             ieff = 0.987;
    }
    else if(idispl<0.5) {
      if     (ipt< 6.) ieff = 1.005;
      else if(ipt<10.) ieff = 1.002;
      else if(ipt<20.) ieff = 1.002;
      else             ieff = 0.991;
    }
    else if(idispl<1.0) {
      if     (ipt< 6.) ieff = 1.018;
      else if(ipt<10.) ieff = 1.007;
      else if(ipt<20.) ieff = 0.985;
      else             ieff = 1.013;
    }
    else {
      if     (ipt< 6.) ieff = 1.008;
      else if(ipt<10.) ieff = 1.021;
      else if(ipt<20.) ieff = 0.976;
      else             ieff = 1.012;
    }
  }
  // 2018
  else {
    if     (idispl<0.2) {
      if     (ipt< 6.) ieff = 0.994;
      else if(ipt<10.) ieff = 0.996;
      else if(ipt<20.) ieff = 0.991;
      else             ieff = 0.986;
    }
    else if(idispl<0.5) {
      if     (ipt< 6.) ieff = 0.993;
      else if(ipt<10.) ieff = 0.996;
      else if(ipt<20.) ieff = 0.997;
      else             ieff = 1.003;
    }
    else if(idispl<1.0) {
      if     (ipt< 6.) ieff = 0.992;
      else if(ipt<10.) ieff = 0.994;
      else if(ipt<20.) ieff = 1.009;
      else             ieff = 0.999;
    }
    else {
      if     (ipt< 6.) ieff = 1.011;
      else if(ipt<10.) ieff = 1.023;
      else if(ipt<20.) ieff = 0.997;
      else             ieff = 0.995;
    }
  }

  return ieff;
}





/*
//
// ========================================================
//   Write tables (NEED TO CHECK AFTER LAST CHANGES!!!!)
// ========================================================
//
if(skipTables==false) {
std::ofstream tabletexS, tabletexL;
tabletexS.open("tabelle/tables_"+sgn+"_"+cpl+"_short.txt");
tabletexL.open("tabelle/tables_"+sgn+"_"+cpl+"_long.txt");
size_t nsrbins = plots_SR[icoup][0][0][1+isign]->GetNbinsX();
std::vector<double> totconts(nsrbins+3, 0.0);
std::vector<double> totstats(nsrbins+3, 0.0);
std::vector<double> binconts(3, 0.0);
std::vector<double> binstats(3, 0.0);
//
// Write table: signal
// Row header
tabletexL << left << std::setw(2*ntab) << "  signal";
tabletexS << left << std::setw(2*ntab) << "  signal";
for(size_t ibin=0; ibin<nsrbins; ++ibin) {
tabletexL << " & $"   << left << std::setw(ntab/2) <<  plots_SR[icoup][0][0][1+isign]->GetBinContent(ibin+1)
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << plots_SR[icoup][0][0][1+isign]->GetBinError(ibin+1)
<< "$";
// Group by final state
/// >>> WARNING: if bin numbering changes, this needs to be updated!
size_t ibintmp = (ibin<6 ? 0 : (ibin<12 ? 1 : 2));
binconts[ibintmp] += plots_SR[icoup][0][0][1+isign]->GetBinContent(ibin+1);
binstats[ibintmp] += plots_SR[icoup][0][0][1+isign]->GetBinError(ibin+1) * plots_SR[icoup][0][0][1+isign]->GetBinError(ibin+1);
}
//
for(size_t ibintmp=0; ibintmp<3; ++ibintmp) {
tabletexS << " & $"   << left << std::setw(ntab/2)  << binconts[ibintmp]
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << std::sqrt(binstats[ibintmp])
<< "$";
binconts[ibintmp] = 0.;
binstats[ibintmp] = 0.;
}
//
tabletexL << " \\\\\n  \\hline\n";
tabletexS << " \\\\\n  \\hline\n";
 
//
// Write table: backgrounds
for(unsigned bkg=0; bkg<nBkg; ++bkg) {
// Row header
tabletexL << left << std::setw(2*ntab) << ("  "+labelPerProc[bkgNames[bkg]]);
tabletexS << left << std::setw(2*ntab) << ("  "+labelPerProc[bkgNames[bkg]]);
for(size_t ibin=0; ibin<nsrbins; ++ibin) {
tabletexL << " & $"   << left << std::setw(ntab/2)  << plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinContent(ibin+1)
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinError(ibin+1)
<< "$";
// Add to total background
totconts[ibin] += plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinContent(ibin+1);
totstats[ibin] += plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinError(ibin+1) * plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinError(ibin+1);
// Group by final state
/// >>> WARNING: if bin numbering changes, this needs to be updated!
size_t ibintmp = (ibin<6 ? 0 : (ibin<12 ? 1 : 2));
binconts[ibintmp] += plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinContent(ibin+1);
binstats[ibintmp] += plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinError(ibin+1) * plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinError(ibin+1);
// Add to total background!
totconts[nsrbins+ibintmp] += plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinContent(ibin+1);
totstats[nsrbins+ibintmp] += plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinError(ibin+1) * plots_SR[icoup][0][0][1+nSamples_signal+bkg]->GetBinError(ibin+1);
}
//
for(size_t ibintmp=0; ibintmp<3; ++ibintmp) {
tabletexS << " & $"   << left << std::setw(ntab/2) << binconts[ibintmp]
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << std::sqrt(binstats[ibintmp])
<< "$";
binconts[ibintmp] = 0.;
binstats[ibintmp] = 0.;
}
//
tabletexL << " \\\\\n  \\hline\n";
tabletexS << " \\\\\n  \\hline\n";
}

//
// Write table: total background
// Row header
tabletexL << left << std::setw(2*ntab) << "  Total background";
tabletexS << left << std::setw(2*ntab) << "  Total background";
for(size_t ibin=0; ibin<nsrbins; ++ibin) {
tabletexL << " & $"   << left << std::setw(ntab/2) << totconts[ibin]
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << std::sqrt(totstats[ibin])
<< "$";
totconts[ibin] = 0.;
totstats[ibin] = 0.;
}
//
for(size_t ibintmp=0; ibintmp<3; ++ibintmp) {
tabletexS << " & $"   << left << std::setw(ntab/2)  << totconts[nsrbins+ibintmp]
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << std::sqrt(totstats[nsrbins+ibintmp])
<< "$";
totconts[nsrbins+ibintmp] = 0.;
totstats[nsrbins+ibintmp] = 0.;
}
//
tabletexL << " \\\\\n  \\hline\n";
tabletexS << " \\\\\n  \\hline\n";

//
// Write table: data
// Row header
tabletexL << left << std::setw(2*ntab) << "  Observed";
tabletexS << left << std::setw(2*ntab) << "  Observed";
for(size_t ibin=0; ibin<nsrbins; ++ibin) {
tabletexL << " & $"   << left << std::setw(ntab/2)  << 0 //dataYields[0][couplidx[icoup]][6]->GetBinContent(ibin+1)
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << 0 //dataYields[0][couplidx[icoup]][6]->GetBinError(ibin+1)
<< "$";
// Group by final state
/// >>> WARNING: if bin numbering changes, this needs to be updated!
size_t ibintmp = (ibin<6 ? 0 : (ibin<12 ? 1 : 2));
binconts[ibintmp] += 0; //dataYields[0][couplidx[icoup]][6]->GetBinContent(ibin+1);
binstats[ibintmp] += 0; //dataYields[0][couplidx[icoup]][6]->GetBinError(ibin+1) * dataYields[0][couplidx[icoup]][6]->GetBinError(ibin+1);
}
//
for(size_t ibintmp=0; ibintmp<3; ++ibintmp) {
tabletexS << " & $"   << left << std::setw(ntab/2)  << binconts[ibintmp]
<< " \\pm " << left << std::setw(ntab/2) << std::setprecision(2) << std::sqrt(binstats[ibintmp])
<< "$";
binconts[ibintmp] = 0.;
binstats[ibintmp] = 0.;
}
//
tabletexL << " \\\\\n  \\hline\n";
tabletexS << " \\\\\n  \\hline\n";

tabletexS.close();
tabletexL.close();
//
// ========================================================
//
} // end if(skipTables==false)
*/

