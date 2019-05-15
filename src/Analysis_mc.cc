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
Analysis_mc::Analysis_mc(unsigned jaar, const std::string& list, const std::string& directory) : TObject()

{
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
void Analysis_mc::initSample(unsigned jaar,const Sample& samp){ 

  //update current sample
  currentSample = samp;
  sampleFile = samp.getFile();
  sampleFile->cd("blackJackAndHookers");
  fChain = (TTree*) sampleFile->Get("blackJackAndHookers/blackJackAndHookersTree");
  initTree(fChain, samp.isData());
  nEntries = fChain->GetEntries();
  if(!samp.isData()){

    //read sum of simulated event weights
    TH1D* hCounter = new TH1D("hCounter", "Events counter", 1, 0, 1);
    hCounter->Read("hCounter"); 
    double sumSimulatedEventWeights = hCounter->GetBinContent(1);
    delete hCounter;

    //event weights set with lumi depending on sample's era 
    double dataLumi;
    if( jaar == 0 ){
      dataLumi = lumi2016;
    }
    else if ( jaar == 1 ){
      dataLumi = lumi2017;
    }
    else dataLumi = lumi2018;

    // N.B.: getXSec() returns the cross section, or the *re-weighted* cross section in case of V^2 (or ctau) re-weighting
    scale = samp.getXSec()*dataLumi*1000/sumSimulatedEventWeights;       //xSec*lumi divided by total sum of simulated event weights
  }
}
//_______________________________________________________ initialize sample ____
void Analysis_mc::initSample(){ //initialize the next sample in the list 
  initSample(jaar,samples[++currentSampleIndex]);
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
void Analysis_mc::initTree(TTree *tree, const bool isData, unsigned jaar)
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
  fChain->SetBranchAddress("_Flag_goodVertices", &_Flag_goodVertices, &b__Flag_goodVertices);
  fChain->SetBranchAddress("_Flag_HBHENoiseFilter", &_Flag_HBHENoiseFilter, &b__Flag_HBHENoiseFilter);
  fChain->SetBranchAddress("_Flag_HBHENoiseIsoFilter", &_Flag_HBHENoiseIsoFilter, &b__Flag_HBHENoiseIsoFilter);
  fChain->SetBranchAddress("_Flag_EcalDeadCellTriggerPrimitiveFilter", &_Flag_EcalDeadCellTriggerPrimitiveFilter, &b__Flag_EcalDeadCellTriggerPrimitiveFilter);
  fChain->SetBranchAddress("_Flag_BadPFMuonFilter", &_Flag_BadPFMuonFilter, &b__Flag_BadPFMuonFilter);
  fChain->SetBranchAddress("_Flag_BadChargedCandidateFilter", &_Flag_BadChargedCandidateFilter, &b__Flag_BadChargedCandidateFilter);
  //fChain->SetBranchAddress("_updated_ecalBadCalibFilter", &_updated_ecalBadCalibFilter, &b__updated_ecalBadCalibFilter);  
  fChain->SetBranchAddress("_passTrigger_1l", &_passTrigger_1l, &b__passTrigger_1l);   
  fChain->SetBranchAddress("_HLT_IsoMu24", &_HLT_IsoMu24, &b__HLT_IsoMu24);
  if (jaar==0) fChain->SetBranchAddress("_HLT_IsoTkMu24", &_HLT_IsoTkMu24, &b__HLT_IsoTkMu24);
  if (jaar==0) fChain->SetBranchAddress("_HLT_Ele27_WPTight_Gsf", &_HLT_Ele27_WPTight_Gsf, &b__HLT_Ele27_WPTight_Gsf);   
  if (jaar!=0)fChain->SetBranchAddress("_HLT_IsoMu27", &_HLT_IsoMu27, &b__HLT_IsoMu27);  
  if (jaar!=0)fChain->SetBranchAddress("_HLT_Ele32_WPTight_Gsf", &_HLT_Ele32_WPTight_Gsf, &b__HLT_Ele32_WPTight_Gsf);
  if (jaar!=0)fChain->SetBranchAddress("_HLT_Ele35_WPTight_Gsf", &_HLT_Ele35_WPTight_Gsf, &b__HLT_Ele35_WPTight_Gsf);
  if (jaar!=0)fChain->SetBranchAddress("_HLT_Ele32_WPTight_Gsf_L1DoubleEG", &_HLT_Ele32_WPTight_Gsf_L1DoubleEG, &b__HLT_Ele32_WPTight_Gsf_L1DoubleEG);   
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
  fChain->SetBranchAddress("_lPtCorr", _lPtCorr, &b__lPtCorr);
  fChain->SetBranchAddress("_lPtScaleUp", _lPtScaleUp, &b__lPtScaleUp);
  fChain->SetBranchAddress("_lPtScaleDown", _lPtScaleDown, &b__lPtScaleDown);
  fChain->SetBranchAddress("_lPtResUp", _lPtResUp, &b__lPtResUp);
  fChain->SetBranchAddress("_lPtResDown", _lPtResDown, &b__lPtResDown);
  fChain->SetBranchAddress("_lECorr", _lECorr, &b__lECorr);
  fChain->SetBranchAddress("_lEScaleUp", _lEScaleUp, &b__lEScaleUp);
  fChain->SetBranchAddress("_lEScaleDown", _lEScaleDown, &b__lEScaleDown);
  fChain->SetBranchAddress("_lEResUp", _lEResUp, &b__lEResUp);
  fChain->SetBranchAddress("_lEResDown", _lEResDown, &b__lEResDown);
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
  fChain->SetBranchAddress("_jetIsLoose", _jetIsLoose, &b__jetIsLoose);
  fChain->SetBranchAddress("_jetIsTight", _jetIsTight, &b__jetIsTight);
  fChain->SetBranchAddress("_jetIsTightLepVeto", _jetIsTightLepVeto, &b__jetIsTightLepVeto);
  fChain->SetBranchAddress("_jetNeutralHadronFraction", _jetNeutralHadronFraction, &b__jetNeutralHadronFraction);
  fChain->SetBranchAddress("_jetChargedHadronFraction", _jetChargedHadronFraction, &b__jetChargedHadronFraction);
  fChain->SetBranchAddress("_jetNeutralEmFraction", _jetNeutralEmFraction, &b__jetNeutralEmFraction);
  fChain->SetBranchAddress("_jetChargedEmFraction", _jetChargedEmFraction, &b__jetChargedEmFraction);
  fChain->SetBranchAddress("_jetHFHadronFraction", _jetHFHadronFraction, &b__jetHFHadronFraction);
  fChain->SetBranchAddress("_jetHFEmFraction", _jetHFEmFraction, &b__jetHFEmFraction);
  fChain->SetBranchAddress("_met", &_met, &b__met);
  fChain->SetBranchAddress("_metRaw", &_metRaw, &b__metRaw);
  fChain->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
  fChain->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
  fChain->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
  fChain->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
  fChain->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
  fChain->SetBranchAddress("_metRawPhi", &_metRawPhi, &b__metRawPhi);
  fChain->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
  fChain->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
  fChain->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
  fChain->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
  fChain->SetBranchAddress("_metSignificance", &_metSignificance, &b__metSignificance);

  if(!isData){
    fChain->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
    fChain->SetBranchAddress("_weight", &_weight, &b__weight);
    fChain->SetBranchAddress("_lheHTIncoming", &_lheHTIncoming, &b__lheHTIncoming);
    fChain->SetBranchAddress("_ctauHN", &_ctauHN, &b__ctauHN);
    fChain->SetBranchAddress("_nLheTau", &_nLheTau, &b__nLheTau);
    fChain->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
    fChain->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
    fChain->SetBranchAddress("_nPsWeights", &_nPsWeights, &b__nPsWeights);
    fChain->SetBranchAddress("_psWeight", _psWeight, &b__psWeight);
    fChain->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
    fChain->SetBranchAddress("_gen_pdgID", _gen_pdgID, &b__gen_pdgID);
    fChain->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
    fChain->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
    fChain->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
    fChain->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
    fChain->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
    fChain->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
    fChain->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
    fChain->SetBranchAddress("_gen_vertex_x", _gen_vertex_x, &b__gen_vertex_x);
    fChain->SetBranchAddress("_gen_vertex_y", _gen_vertex_y, &b__gen_vertex_y);
    fChain->SetBranchAddress("_gen_vertex_z", _gen_vertex_z, &b__gen_vertex_z);
    fChain->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
    fChain->SetBranchAddress("_gen_lMinDeltaR", _gen_lMinDeltaR, &b__gen_lMinDeltaR);
    fChain->SetBranchAddress("_gen_lPassParentage", _gen_lPassParentage, &b__gen_lPassParentage);    
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
  }
}

//          ================= ================= ================= ================= ================= =================          // 
//          ================= ================= ================= ================= ================= =================          //

//_______________________________________________________ analysis function ____
void Analysis_mc::analisi( unsigned jaar, const std::string& list, const std::string& directory,
			   TString outfilename,
			   int systcat, int systdir
			   ) {

  std::ofstream zero("zero.txt"); 
  std::ofstream one("one.txt");  
  std::ofstream two("two.txt"); 
  std::ofstream three("three.txt"); 
  std::ofstream four("four.txt"); 
  std::ofstream five("five.txt"); 
 


  
  cout<<"in analisi"<<endl;
  cout<<"---------------------------"<<endl;   
  setTDRStyle();
  if(systdir<0) {
    std::cout << " >>> Dummy message (to avoid warnings): systdir " << systdir << std::endl;
  }

  TFile *fout = new TFile(outfilename.Data(), "recreate");
  
  // ------------ pile up -----------------------------------------------//
  TH1D *pileUpWeight[1];

  TFile hfile_pu("/user/mvit/CMSSW_9_4_4/src/HNL_analysis/PU/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root");
  //TFile hfile_pu("/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/ANALYSIS/20190318_MartinasCode/samples.noSync/2016/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root");
  pileUpWeight[0] = (TH1D*)hfile_pu.Get("puw_Run2016Inclusive_central");


  TGraphAsymmErrors *fakeRate_mu[3];
  TGraphAsymmErrors *fakeRate_e[3];
  TGraphAsymmErrors *fakeRate_mumu[3];
  TGraphAsymmErrors *fakeRate_ee[3];
  TGraphAsymmErrors *fakeRate_mue[3];


  TFile hfile1("/user/mvit/CMSSW_9_4_4/src/closure_2016/FR/fake_rate_mu.root");
  fakeRate_mu[0] = (TGraphAsymmErrors*)hfile1.Get("fakeRate_mu_eta1");
  fakeRate_mu[1] = (TGraphAsymmErrors*)hfile1.Get("fakeRate_mu_eta2");
  fakeRate_mu[2] = (TGraphAsymmErrors*)hfile1.Get("fakeRate_mu_eta3");
  TFile hfile2("/user/mvit/CMSSW_9_4_4/src/closure_2016/FR/fake_rate_e.root");
  fakeRate_e[0] = (TGraphAsymmErrors*)hfile2.Get("fakeRate_e_eta1");
  fakeRate_e[1] = (TGraphAsymmErrors*)hfile2.Get("fakeRate_e_eta2");
  fakeRate_e[2] = (TGraphAsymmErrors*)hfile2.Get("fakeRate_e_eta3");
  TFile hfile_dfr1("/user/mvit/CMSSW_9_4_4/src/closure_2016/FR/fake_rate_mumu.root");
  fakeRate_mumu[0]= (TGraphAsymmErrors*)hfile_dfr1.Get("fakeRate_mu_eta1");
  fakeRate_mumu[1]= (TGraphAsymmErrors*)hfile_dfr1.Get("fakeRate_mu_eta2");
  fakeRate_mumu[2]= (TGraphAsymmErrors*)hfile_dfr1.Get("fakeRate_mu_eta3");
  TFile hfile_dfr2("/user/mvit/CMSSW_9_4_4/src/closure_2016/FR/fake_rate_ee.root");
  fakeRate_ee[0]= (TGraphAsymmErrors*)hfile_dfr2.Get("fakeRate_e_eta1");
  fakeRate_ee[1]= (TGraphAsymmErrors*)hfile_dfr2.Get("fakeRate_e_eta2");
  fakeRate_ee[2]= (TGraphAsymmErrors*)hfile_dfr2.Get("fakeRate_e_eta3");
  TFile hfile_dfr3("/user/mvit/CMSSW_9_4_4/src/closure_2016/FR/fake_rate_emu.root");
  fakeRate_mue[0]= (TGraphAsymmErrors*)hfile_dfr3.Get("fakeRate_emu_eta1");
  fakeRate_mue[1]= (TGraphAsymmErrors*)hfile_dfr3.Get("fakeRate_emu_eta2");
  fakeRate_mue[2]= (TGraphAsymmErrors*)hfile_dfr3.Get("fakeRate_emu_eta3");




  
  if (jaar == 0 ) {
   
  }
  else if (jaar == 1 ) {

  }
  else {

  }
  // ------------ b tagging -----------------------------------------------//
  // b-tagging working points (DeepCsv_b + DeepCsv_bb)
  double btagCuts[3][3];
  // selected WP (0: loose; 1: medium; 2: tight)
  BTagEntry::OperatingPoint bwp = BTagEntry::OP_LOOSE;    // = 0
  //BTagEntry::OperatingPoint bwp = BTagEntry::OP_MEDIUM; // = 1
  //BTagEntry::OperatingPoint bwp = BTagEntry::OP_TIGHT;  // = 2

  //  - 2016
  btagCuts[0][0] = bjet_loose_2016; // loose
  btagCuts[0][1] = 0.6321; // medium
  btagCuts[0][2] = 0.8953; // tight
  //  - 2017
  btagCuts[1][0] = bjet_loose_2017; // loose
  btagCuts[1][1] = 0.4941; // medium
  btagCuts[1][2] = 0.8001; // tight
  //  - 2018
  btagCuts[2][0] = bjet_loose_2018; // loose
  btagCuts[2][1] = 0.4184; // medium
  btagCuts[2][2] = 0.7527; // tight

  // B-tagging calibration + reader
  BTagCalibration calib("DeepCSV", (year==0 ? "DeepCSV_2016LegacySF_WP_V1.csv" : (year==1 ? "DeepCSV_94XSF_WP_V4_B_F.csv" : "DeepCSV_102XSF_WP_V1.csv")));
  BTagCalibrationReader reader(bwp,             // working point
			       "central",       // central sys type
			       {"up", "down"}); // other sys types

  reader.load(calib,             // calibration instance
	      BTagEntry::FLAV_B, // b-tag flavor
	      "comb");           // measurement type



  // ------------   samples info -----------------------------------------------//
  std::vector <Sample> samples  = readSampleList(list, directory);

  //std::vector <Sample> samples  = readSampleList(list, directory);
  /*
    if (jaar == 0) {
    const int nSamples = samples.size();
    const int nSamples_eff = 2;
    const int nSamples_signal = 2;
    }
    else if (jaar == 1 ) {
    const int nSamples = samples.size();
    const int nSamples_eff = 2;
    const int nSamples_signal = 2;
    }
    else {
    const int nSamples = samples.size();
    const int nSamples_eff = 2;
    const int nSamples_signal = 2;
    } 
  */
  TH1D* Histos[nDist][nChannel][nCat][nSamples_eff +1];
		       
  for(int i = 0; i < nDist; ++i){
    float BinWidth = (HistMax[i] - HistMin[i])/nBins[i];
    std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
    for(int effsam = 0; effsam < nSamples_eff + 1; ++effsam){
      for(int cat = 0; cat < nCat; ++cat){
	for(int cha = 0; cha < nChannel; ++cha){               
	  Histos[i][cha][cat][effsam] = new TH1D(eff_names[effsam] +"_"+ channelNames[cha] +"_"+ catNames[cat] +"_"+ Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + "; events /" + Yaxis + Units[i], nBins[i], HistMin[i], HistMax[i]);
	  Histos[i][cha][cat][effsam]->Sumw2();
	}
      }
    }
  }

  //Calculate the center of the maximum bin of each histogram
  double maxBinC[nDist];
  for(int i = 0; i < nDist; ++i){
    maxBinC[i] = Histos[i][0][0][0]->GetBinCenter(Histos[i][0][0][0]->GetNbinsX());
  }
  
  // ------------   run over samples -----------------------------------------------//

  for(int sam = 0,effsam = 0; sam < samples.size(); ++sam, ++effsam){

    initSample(jaar,samples[sam]);
    if (samples[sam].isData()) continue;
    //check consistency
    std::cout << "sample initialized: --> " << std::endl;
    std::cout << "fileName: " << samples[sam].getFileName() << "  process name: " << samples[sam].getProcessName() << "   xsec: " << samples[sam].getXSec() << std::endl;
    if(samples[sam].isData()) std::cout << " is Data" << std::endl;
    if(samples[sam].isMC()  ) std::cout << " is MC"   << std::endl;
    if(sam != 0){
      if(samples[sam].getProcessName() == samples[sam-1].getProcessName()) --effsam;     
    }
    if (samples[sam].getProcessName() != "DY") continue;
    // For lifetime re-weighting (hip hip hip hurray)
    double ctauOld(0.), ctauNew(0.), ctWeight(1.);
    /* if(samples[sam].isNewPhysicsSignal()) {
       std::cout << " is signal" << std::endl;
       if(samples[sam].getHNLV2New()>0.) {
       ctauOld = samples[sam].getHNLctau();
       ctauNew = samples[sam].getHNLctauNew();
       std::cout << "  ==> HNL lifetime re-weighting: " << std::endl;
       std::cout << "      (" << samples[sam].getHNLV2()    << ", " << ctauOld
       << ") --> (" << samples[sam].getHNLV2New() << ", " << ctauNew
       << ")" << std::endl;

       ctWeight = (ctauOld/ctauNew) * TMath::Exp(((1./ctauOld)-(1./ctauNew))*_ctauHN);
       }
       }*/

    double progress = 0; 	//For printing progress bar 
    // ------------   run over entries -----------------------------------------------//  
    for (Long64_t it = 0; it < nEntries; ++it){
      GetEntry(samples[sam], it);
      //   std::cout<<"after get tree"<<std::endl;

      //print progess
      /* if(it%100 == 0 && it != 0){
	 progress += (double) (100./nEntries);
	 printProgress(progress);
	 } else if(it == nEntries -1){
	 progress = 1.;
	 printProgress(progress);
	 }*/
      // N.B.: ctWeight = 1 unless it is a ctau-reweighted signal sample
      ctWeight = 1;
      double scal = 0;
      scal = scale * _weight * ctWeight * pu_weight(*&pileUpWeight[0],_nTrueInt);
      bwght = 1.;

      // std::cout<<"after pu"<<std::endl;

      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PARAMETERS AND CUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      std::vector<unsigned> ind;
      double*           conePt = new double[_nL];
      double           _ptReal[_nL];
      double           _EReal[_nL];
      Bool_t            _passedMVA90[_nL];     
      
      unsigned         ind_new_leading=0;
      unsigned         ind_new_p=0;
      unsigned         ind_new_pp=0;
      unsigned*         _isLooseCutBasedElectronWithoutIsolatio= new unsigned[_nL];
      unsigned*         _isOurMedium= new unsigned[_nL];
      unsigned*         _passTimingVeto= new unsigned[_nL];
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
      unsigned        l1=0;
      unsigned        l2=0;
      unsigned        l3=0;
      TLorentzVector  v4l1;
      TLorentzVector  v4l2;
      TLorentzVector  v4l3;
      TLorentzVector  v4l2_naked;
      TLorentzVector  v4l3_naked;
      TLorentzVector  v4l2_propagated;
      TLorentzVector  v4l3_propagated;

      double            _vertex_X=-1;
      double            _vertex_Y=-1;
      double            _vertex_Z=-1;
      double            _vertex_R2D=-1;
      double            _vertex_sR2D=-1;
      double            _vertex_R=-1;
      double            _vertex_sR=-1;
      double            _vertex_chi2=-1;
      double            _vertex_normchi2=-1;
      double _vertex_ndf =-1;

      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //------------------------------------------------------------ lepton selection for FO
      for(unsigned i = 0; i < _nL; ++i){
	_ptReal[i]=_lPt[i];
	_EReal[i] =_lE[i];
      }
      int _lIndex[_nL];
       for(unsigned i = 0; i < _nL; ++i){
	_lIndex[i] = i+1;
      }


      
       //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<""<<std::endl;
       //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<""<<std::endl;
       //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<""<<std::endl;

      

      //select leptons
      //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"event 203"<<"  "<<_eventNb<<std::endl;
      const unsigned lCount = selectLepConeCorr(ind);
      if (lCount < 3) continue;
      //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"after 3 event 203"<<"  "<<_eventNb<<std::endl;

      //------------------------------------------------------------ jet pt variation and nJet and bjet
      /* for (unsigned j =0; j < _nJets ; j++){
	 _jetPt[j]=_jetSmearedPt[j];
	 if(systcat==8) {
	 if(systdir==0) _jetPt[j]=_jetSmearedPt_JECDown[j];	   
	 else _jetPt[j]=_jetSmearedPt_JECUp[j];	   
	 }
	 else if(systcat==9) {
	 if(systdir==0)  _jetPt[j]=_jetSmearedPt_JERDown[j];	  
	 else  _jetPt[j]=_jetSmearedPt_JERUp[j];	  
	 }*/
      for (unsigned j =0; j < _nJets ; j++){
	if(jetIsBJet(j)  && _jetPt[j]<1000. && std::abs(_jetEta[j])<2.4) {
	  double bjetSf = 1.;
	  // b-jet systematics
	  if(systcat==10) {
	    if(systdir==0)  bjetSf = reader.eval_auto_bounds("down", BTagEntry::FLAV_B, std::abs(_jetEta[j]), _jetPt[j]);	  
	    else  bjetSf = reader.eval_auto_bounds("up", BTagEntry::FLAV_B, std::abs(_jetEta[j]), _jetPt[j]);	    
	  }
	  // b-jet central SF
	  else bjetSf = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, std::abs(_jetEta[j]), _jetPt[j]);
	  // Scale the b-veto event weight
	  bwght *= bjetSf;
	}	
      }
      //counting bjet and njet
      for (unsigned j =0; j < _nJets ; j++){
	if (jetIsGood(j)) ++goodjet;
	if (jetIsBJet(j)) ++bjet;
      }
      // ------------ ==================== -----------------------------------------------//
      // ------------   event selection   -----------------------------------------------//
      //assign the l1 index
      ind_new_leading = l1Index(ind);
      if (l1Index(ind) == -1) continue; //in case there are not l1 at all
      //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"after leading event 203"<<"  "<<_eventNb<<"   pt leading: "<< _lPt[ind_new_leading]<<std::endl;

      //check how many displaced there are (displaced --> dxy, common vertex, FO, no l1)
      unsigned displacedC = 0;
      std::vector<TLorentzVector> lepV_displaced;
      std::vector<int> charge_displaced;
      std::vector<unsigned> temp_index;
      
      

      
      int index_to_use_for_l2_l3[2]={0,0};
      //find the right OS pair with min invariant mass
      int min_test= 9999;
      int min_mass=999;
      displacedC=0;
      	  if (_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971)std::cout<<"=========================================================="<<std::endl;

      for(unsigned l = 0; l < lCount; ++l){
		  if (_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971)std::cout<<"      "<<std::endl;

	for(unsigned j = l+1; j < lCount; ++j){
	  
	  // if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"==============="<<std::endl;
	  if (!(_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971))continue;
		  
	  //if (!(_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 )) continue;
	  //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"         first: "<< l<<" ind[l] "<< ind[l]<<"   second "<<j<<" ind[j] "<< ind[j]<<std::endl;
	  //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"         firstI : "<<_lIndex[ind[l]]<<"  second I "<<_lIndex[ind[j]]<<std::endl;

	  
	  //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"inside loop leptons displaced"<<"  "<<_eventNb<<"  pt  "<<_lPt[ind[l]]<<"   "<<_lPt[ind[j]]<<std::endl;
	  //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"dxy: "<<fabs(_dxy[ind[l]])<<"  "<<fabs(_dxy[ind[j]])<<"   reliso "<<_relIso[ind[l]]<<"  "<<_relIso[ind[j]]<<std::endl;

	  //std::cout<<"calling the function lepIsDisplaced with those index: fixed one: "<<ind[l]<< "    looping over the rest"<<std::endl;
	  //if(!lepIsDisplaced(ind[l] , ind_new_leading, ind)) continue;
	  //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"l is dispalced"<<"  "<<_eventNb<<_lPt[ind[l]]<<std::endl;

	  // if(!lepIsDisplaced(ind[j] , ind_new_leading, ind)) continue;
	  //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"j is dispalced"<<"  "<<_eventNb<<_lPt[ind[j]]<<std::endl;

	  // if (_lCharge[ind[l]] == _lCharge[ind[j]]) continue;
	  //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"OS"<<"  "<<_eventNb<<std::endl;

	  if (_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971)std::cout<<"before chioce pair with  "<<ind[l]<<"  "<< ind[j]<<"  "<<std::endl;
	  if (!IsDisplacedPair(ind[l] ,ind[j], ind_new_leading, ind)) continue;
	   if (_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971)std::cout<<"after chioce pair with  "<<ind[l]<<"  "<< ind[j]<<"  "<<std::endl;
	   if (_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971){
	     for(unsigned v = 0; v < _nVFit_os; ++v){
	       if ((_vertices_os[v][0] == (_lIndex[ind[l]] * 100 +  _lIndex[ind[j]] )) ||(_vertices_os[v][0] == (_lIndex[ind[l]] +  _lIndex[ind[j]] *100) )) std::cout<<"was the vertex found:  "<<_vertices_os[v][0]<<"   x   "<< _vertices_os[v][1]<<"    where the lindex were:  "<<_lIndex[ind[l]]<<"  "<<_lIndex[ind[j]]<< std::endl;
	     }

	   }


	   
	  ++displacedC;
	  TLorentzVector temp_displaced1;
	  TLorentzVector temp_displaced2;
	  temp_displaced1.SetPtEtaPhiE(_lPt[ind[l]],_lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
	  temp_displaced2.SetPtEtaPhiE(_lPt[ind[j]],_lEta[ind[j]], _lPhi[ind[j]], _lE[ind[j]]);

	  std::cout<<"**********************   "<< (temp_displaced1+temp_displaced2).M()<<std::endl;
	  std::cout<<"min mass  "<< min_mass<<std::endl;
	  if ( (temp_displaced1+temp_displaced2).M()  < min_mass) {
	    min_mass= (temp_displaced1+temp_displaced2).M();
	    if (_lPt[ind[l]]> _lPt[ind[j]]){
	      index_to_use_for_l2_l3[0] = ind[l];
	      index_to_use_for_l2_l3[1] = ind[j];
	    }
	    if (_lPt[ind[l]] < _lPt[ind[j]]){
	      index_to_use_for_l2_l3[0] = ind[j];
	      index_to_use_for_l2_l3[1] = ind[l];
	    }	    
	  }
	  std::cout<<"       min mass  "<< min_mass<<std::endl;

	  //std::cout<<"mass min: "<<min_mass<<std::endl;
	}//end loop2
      }//end loop1

      if (displacedC < 1) continue;
      //if (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"after displaced event 203"<<"  "<<_eventNb<<std::endl;

      
      //trigger NOT trigger matching!!!!!!
      if (!_passTrigger_1l) continue;

      //f (_eventNb==96541 || _eventNb==113885 || _eventNb==134456 || _eventNb==136224 ) std::cout<<"after trigger event 203"<<"  "<<_eventNb<<std::endl;

      
      if (_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971){
	std::cout<<"================"<<std::endl;
	for(unsigned l = 0; l < lCount; ++l){
	  std::cout<<"leading:   "<< ind_new_leading <<"  pt "<<_lPt[ ind_new_leading]<<std::endl;
	  std::cout<<"is it displcade: "<<lepIsDisplaced(ind[l] , ind_new_leading, ind)<<std::endl;
	  std::cout<<l<<": index "<<ind[l]<< "  pt: "<<_lPt[ind[l]]<<"  relIso: "<< _relIso[ind[l]]<<"  dxy "<< fabs(_dxy[ind[l]])<<"flav: "<< _lFlavor[ind[l]]<<"   ourmedium: "<<muOurMedium(ind[l])<<"   mediumPOG: "<<_lPOGMedium[ind[l]]<<"  "<<_lCharge[ind[l]]<<std::endl;
	}

	std::cout<<"----------->   picked: "<< _lPt[index_to_use_for_l2_l3[0]]<<"   "<<_lPt[index_to_use_for_l2_l3[1]]<<std::endl;
	std::cout<< "was the vertex found:  "<< l2l3_vertex_variable (index_to_use_for_l2_l3[0],index_to_use_for_l2_l3[1])<<"   "<<_vertices[l2l3_vertex_variable (index_to_use_for_l2_l3[0],index_to_use_for_l2_l3[1])][0]<< "  "<<"   x  "<<_vertices[l2l3_vertex_variable (index_to_use_for_l2_l3[0],index_to_use_for_l2_l3[1])][1]<<std::endl;
	for(unsigned l = 0; l < lCount; ++l){
	  for(unsigned j = l+1; j < lCount; ++j){
	    TLorentzVector temp_displaced1;
	    TLorentzVector temp_displaced2;
	    temp_displaced1.SetPtEtaPhiE(_lPt[ind[l]],_lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
	    temp_displaced2.SetPtEtaPhiE(_lPt[ind[j]],_lEta[ind[j]], _lPhi[ind[j]], _lE[ind[j]]);
	    std::cout<<"mass made by: "<< ind[l]<<"  and  "<< ind[j]<<"  is "<< (temp_displaced1+temp_displaced2).M()<<std::endl;
	  }
	}
      }

      
      
      // ------------ changing all the lep info and vertex-----------------------------------------------//
      l1=ind_new_leading;
      l2=index_to_use_for_l2_l3[0];
      l3=index_to_use_for_l2_l3[1];
      v4l1.SetPtEtaPhiE(_lPt[l1],_lEta[l1], _lPhi[l1], _lE[l1]);
      v4l2.SetPtEtaPhiE(_lPt[l2],_lEta[l2], _lPhi[l2], _lE[l2]);
      v4l3.SetPtEtaPhiE(_lPt[l3],_lEta[l3], _lPhi[l3], _lE[l3]);
      v4l2_naked.SetPtEtaPhiE(_ptReal[l2],_lEta[l2], _lPhi[l2], _EReal[l2]);
      v4l3_naked.SetPtEtaPhiE(_ptReal[l3],_lEta[l3], _lPhi[l3], _EReal[l3]);
      flavors_3l[0]=_lFlavor[l1];
      flavors_3l[1]=_lFlavor[l2];
      flavors_3l[2]=_lFlavor[l3];
      charge_3l[0]=_lCharge[l1];
      charge_3l[1]=_lCharge[l2];
      charge_3l[2]=_lCharge[l3];




      

      if (samples[sam].getProcessName() == "DY" )   {    
      zero << Form("%1d %7d %9d\t%+2d (%6.1f)\t%+2d (%6.1f | %6.1f) %1d\t%+2d (%6.1f | %6.1f) %1d\t %6.1f" ,
		   _runNb, _lumiBlock, _eventNb,  
		   (-1)*_lCharge[l1]*(11+2*_lFlavor[l1]),v4l1.Pt(),
		 (-1)*_lCharge[l2]*(11+2*_lFlavor[l2]),v4l2_naked.Pt(),v4l2.Pt(),_lProvenanceCompressed[l2],
		 (-1)*_lCharge[l3]*(11+2*_lFlavor[l3]),v4l3_naked.Pt(),v4l3.Pt(),_lProvenanceCompressed[l3],	
		 _met)<< std::endl;
      }


      //vertex l2l3 info
      int index_l2l3= l2l3_vertex_variable (l2,l3);      
      _vertex_X=_vertices[index_l2l3][1];
      _vertex_Y=_vertices[index_l2l3][2];
      _vertex_Z=_vertices[index_l2l3][3];
      _vertex_chi2=_vertices[index_l2l3][11];
      _vertex_normchi2= _vertices[index_l2l3][11]/_vertices[index_l2l3][10];
      _vertex_ndf =_vertices[index_l2l3][10];



      if (_eventNb==316331001 || _eventNb==300090452 ||_eventNb==279298855 ||_eventNb==258551392 ||_eventNb==111733559 ||_eventNb==108458464 ||_eventNb==66138971){

	std::cout<<"                            "<<l1<<"  "<<l2<<"  "<<l3<<"  "<<v4l1.Pt()<<"  "<<v4l2.Pt()<<"  "<<v4l3.Pt()<<"  "<<index_l2l3<<"  "<<_vertex_X<<"  "<<_vertices[index_l2l3][0]<< std::endl;
      }
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

      if (isCRRun && v4l1.Pt() < 30) continue;
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     sFR and  dRF   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      bool single_fake=false;
      bool Double_fake=false;     
      if ( _closest_l1JetE[l2] ==  _closest_l1JetE[l3] ) Double_fake = true;
      if (!Double_fake) single_fake = true;
      if(Double_fake && _closest_l1JetE[l2] ==0) {
	single_fake = true;
	Double_fake = false;
      }
       // -----------------------------------------------------------//
      if (single_fake && flavors_3l[1] == 1 && v4l2.Pt() < 5) continue;
      if (single_fake && flavors_3l[2] == 1 && v4l3.Pt() < 5) continue;
      if (single_fake && flavors_3l[1] == 0 && v4l2.Pt() < 10) continue;
      if (single_fake && flavors_3l[2] == 0 && v4l3.Pt() < 10) continue;
      if (samples[sam].getProcessName() == "DY" )   {    
      one << Form("%1d %7d %9d\t%+2d (%6.1f)\t%+2d (%6.1f | %6.1f) %1d\t%+2d (%6.1f | %6.1f) %1d\t %6.1f" ,
		   _runNb, _lumiBlock, _eventNb,  
		   (-1)*_lCharge[l1]*(11+2*_lFlavor[l1]),v4l1.Pt(),
		 (-1)*_lCharge[l2]*(11+2*_lFlavor[l2]),v4l2_naked.Pt(),v4l2.Pt(),_lProvenanceCompressed[l2],
		 (-1)*_lCharge[l3]*(11+2*_lFlavor[l3]),v4l3_naked.Pt(),v4l3.Pt(),_lProvenanceCompressed[l3],	
		 _met)<< std::endl;
      }

      // ------------ closest jet info --------------------------------------//
      TLorentzVector  l1Jet[1] ;
      float JEC       ;
      TLorentzVector  lepAwareJet[1] ;
      l1Jet[0].SetPxPyPzE(_closest_l1JetPx[l2],_closest_l1JetPy[l2],_closest_l1JetPz[l2],_closest_l1JetE[l2]);
      JEC             = _closestJEC[l2];
      lepAwareJet[0] = (l1Jet[0] - v4l2_naked - v4l3_naked)*JEC + v4l3_naked + v4l2_naked;  
      double momentum_jet=0.;
      momentum_jet = lepAwareJet[0].Pt();
      if (momentum_jet<10) momentum_jet=12;
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
      if (sideBandRegion) continue;
      if (tightC != 2) continue;
      if (samples[sam].getProcessName() == "DY" )   {    
      two << Form("%1d %7d %9d\t%+2d (%6.1f)\t%+2d (%6.1f | %6.1f) %1d\t%+2d (%6.1f | %6.1f) %1d\t %6.1f" ,
		   _runNb, _lumiBlock, _eventNb,  
		   (-1)*_lCharge[l1]*(11+2*_lFlavor[l1]),v4l1.Pt(),
		 (-1)*_lCharge[l2]*(11+2*_lFlavor[l2]),v4l2_naked.Pt(),v4l2.Pt(),_lProvenanceCompressed[l2],
		 (-1)*_lCharge[l3]*(11+2*_lFlavor[l3]),v4l3_naked.Pt(),v4l3.Pt(),_lProvenanceCompressed[l3],	
		 _met)<< std::endl;
      }

      // ------------------ prompt check for MC ------------------------//
      promptC=0;
      if (_lIsPrompt[l1] || _lProvenanceCompressed[l1]==0) promptC++;
      if (_lIsPrompt[l2] || _lProvenanceCompressed[l2]==0) promptC++;
      if (_lIsPrompt[l3] || _lProvenanceCompressed[l3]==0) promptC++;
      //if (!samples[sam].isData() && promptC!=3) continue;
      // -----------------    applying the FRs    --------------------------------//
      if (sideBandRegion){
	if ( samples[sam].isData()  )scal *= -1;
	if (!samples[sam].isData() )scal  = 1 * scal;
	if (single_fake){
	  if (!_isT[l2]) {
	    double fr = FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
				   _lEta[l2], _lFlavor[l2], _lPt[l2], index_eta,flav_dRF, momentum_jet);
	    scal *= -fr/(1-fr);
	  }
	  if (!_isT[l3]) {
	    double fr = FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
				   _lEta[l3], _lFlavor[l3], _lPt[l3], index_eta,flav_dRF, momentum_jet);
	    scal *= -fr/(1-fr);
	  }	  
	}//sFR
	if (loose_lepton_dFR &&  Double_fake) {
	  double fr = FR_weight (*&fakeRate_mu, *&fakeRate_e, *&fakeRate_mumu,*&fakeRate_ee,*&fakeRate_mue,single_fake, Double_fake,
				 _lEta[l2], _lFlavor[l2], _lPt[l2], index_eta,flav_dRF, momentum_jet);
	  scal *= -fr/(1-fr);
	}
      }//FR
      
      if (single_fake && tightFail_sFR && !_isT[l2] && _relIso[l2] < isolation_tight) continue;
      if (single_fake && tightFail_sFR && !_isT[l3] && _relIso[l3] < isolation_tight) continue;
      if (samples[sam].getProcessName() == "DY" )   {    
      three << Form("%1d %7d %9d\t%+2d (%6.1f)\t%+2d (%6.1f | %6.1f) %1d\t%+2d (%6.1f | %6.1f) %1d\t %6.1f" ,
		   _runNb, _lumiBlock, _eventNb,  
		   (-1)*_lCharge[l1]*(11+2*_lFlavor[l1]),v4l1.Pt(),
		 (-1)*_lCharge[l2]*(11+2*_lFlavor[l2]),v4l2_naked.Pt(),v4l2.Pt(),_lProvenanceCompressed[l2],
		 (-1)*_lCharge[l3]*(11+2*_lFlavor[l3]),v4l3_naked.Pt(),v4l3.Pt(),_lProvenanceCompressed[l3],	
		 _met)<< std::endl;
      }
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     analysis   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      bool internal_conv= true;
      if (_lIsPrompt[l1] && _lMatchPdgId[l1] ==22) internal_conv = false;
      if (_lIsPrompt[l2] && _lMatchPdgId[l2] ==22) internal_conv = false;
      if (_lIsPrompt[l3] && _lMatchPdgId[l3] ==22) internal_conv = false;
      bool external_conv= false;
      if (_lIsPrompt[l1] && _lMatchPdgId[l1] ==22) external_conv = true;
      if (_lIsPrompt[l2] && _lMatchPdgId[l2] ==22) external_conv = true;
      if (_lIsPrompt[l3] && _lMatchPdgId[l3] ==22) external_conv = true;    
      if (samples[sam].getProcessName() == "DY" && !internal_conv) continue;
      if (samples[sam].getFileName() == "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root" && !external_conv) continue;


      if (samples[sam].getProcessName() == "DY" )   {    
      four << Form("%1d %7d %9d\t%+2d (%6.1f)\t%+2d (%6.1f | %6.1f) %1d\t%+2d (%6.1f | %6.1f) %1d\t %6.1f" ,
		   _runNb, _lumiBlock, _eventNb,  
		   (-1)*_lCharge[l1]*(11+2*_lFlavor[l1]),v4l1.Pt(),
		 (-1)*_lCharge[l2]*(11+2*_lFlavor[l2]),v4l2_naked.Pt(),v4l2.Pt(),_lProvenanceCompressed[l2],
		 (-1)*_lCharge[l3]*(11+2*_lFlavor[l3]),v4l3_naked.Pt(),v4l3.Pt(),_lProvenanceCompressed[l3],	
		 _met)<< std::endl;
      }
      // if (photonOverlap (samples[sam])) continue;
      
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
      TVector3 svMpv =secondary_vertex[0]- primary_vertex[0];
      double vtxR     = svMpv.Mag();
      double vtxRvtxPcosAlpha = svMpv.Dot(l2plusl3)/vtxR;
      // -----------------   masses
      double M_3L= (v4l2 + v4l3 + v4l1).M();
      double M_ZPair = (pair[0]+pair[1]).M();
      double M_l2l3 = (v4l2 + v4l3).M();
      double M_3L_combined = (v4l2 + v4l3 + v4l1).M();
      if (Double_fake) M_3L_combined = (v4l2_naked + v4l3_naked + v4l1).M();
      double M_l2l3_combined = (v4l2 + v4l3).M();
      if (Double_fake) M_l2l3_combined = (v4l2_naked + v4l3_naked).M();
      METvec.SetPtEtaPhiE(_met, 0, _metPhi,_met);    
      TLorentzVector to_use_mT;
      to_use_mT.SetPtEtaPhiE(other[0].Pt(),0, other[0].Phi(), other[0].Pt());
      double mT=(to_use_mT+METvec).M();
      // -----------------   function useful  2 --> SR also    --------------------------------//
      // 0 = mmm
      // 1 = mme OS
      // 2 = mme SS
      // 3 = eee
      // 4 = eem OS
      // 5 = eem SS   
      int SR_channel=0;
      SR_channel=channel(flavors_3l, charge_3l);
      if (isSRRun && SR_channel == -1 ) continue;
      //avoid +++ or ---
      if (isSRRun && SR_channel == 0 && charge_3l[0] == charge_3l[1] && charge_3l[0] == charge_3l[2]) continue;
      if (isSRRun && SR_channel == 3 && charge_3l[0] == charge_3l[1] && charge_3l[0] == charge_3l[2]) continue;
      bool less2     =false;
      bool more2_10  =false;
      bool more10    =false;  
      bool less5     =false;
      bool more5     =false;
      if (D2_delta_pv_sv < 2)                         less2   = true;
      if (D2_delta_pv_sv >= 2 && D2_delta_pv_sv < 10) more2_10= true;
      if (D2_delta_pv_sv >= 10 )                      more10  = true;
      if (M_l2l3_combined < 5 )   less5= true;
      if (M_l2l3_combined > 5 )   more5= true;
      //bin histogram SR
      int bin_SR_muonCoupling =0;
      int bin_SR_eleCoupling =0;
      bin_SR_muonCoupling = SR_bin_muon( SR_channel, less2,  more2_10,  more10,  less5,  more5 );
      bin_SR_eleCoupling =  SR_bin_ele( SR_channel, less2,  more2_10,  more10,  less5,  more5 );


      M_3L_combined = (v4l2 + v4l3 + v4l1).M();
      M_l2l3_combined =  (v4l2 + v4l3).M();
      
      bool selection_0=false;
      bool selection_1=false;
      bool selection_2=false;
      bool selection_3=false;
      bool selection_4=false;
      bool selection_5=false;
      bool selection_final=false;
      if (charge_3l[2] != charge_3l[1])                                        selection_0 = true;
      if ( selection_0 && v4l2.DeltaR(v4l3) < 1)                               selection_1 = true;
      if ( selection_1 && bjet == 0 )                                          selection_2 = true;
      if ( selection_2 && M_3L_combined > 45 && M_3L_combined < 85)            selection_3 = true;
      if ( selection_3 && min_delta_phi > 1)                                   selection_4 = true;
      if ( selection_4 && vtxRvtxPcosAlpha > 0.9)                              selection_5 = true;
      if ( selection_5 && M_l2l3_combined < 50)                                selection_final = true;
     
      if (!selection_0) continue;

      
      if (selection_final && samples[sam].getProcessName() == "DY" )   {    
      five << Form("%1d %7d %9d\t%+2d (%6.1f)\t%+2d (%6.1f | %6.1f) %1d\t%+2d (%6.1f | %6.1f) %1d\t %6.1f" ,
		   _runNb, _lumiBlock, _eventNb,  
		   (-1)*_lCharge[l1]*(11+2*_lFlavor[l1]),v4l1.Pt(),
		 (-1)*_lCharge[l2]*(11+2*_lFlavor[l2]),v4l2_naked.Pt(),v4l2.Pt(),_lProvenanceCompressed[l2],
		 (-1)*_lCharge[l3]*(11+2*_lFlavor[l3]),v4l3_naked.Pt(),v4l3.Pt(),_lProvenanceCompressed[l3],	
		 _met)<< std::endl;
      }
      
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<     histogramm   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      double values[nDist] ={static_cast<double>(0) ,static_cast<double>(0) ,
			     v4l1.Pt(),
			     v4l2.Pt(),
			     v4l3.Pt(),
			     M_3L,
			     M_l2l3,
			     M_l2l3,
			     M_l2l3_combined,
			     M_l2l3_combined,
			     M_ZPair,
			     mT,
			     _met,		    
			     static_cast<double>( goodjet),
			     static_cast<double>(bjet),
			     _met,
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
			     D3_delta_pv_sv,
			     D2_delta_pv_sv,
			     D2_delta_pv_sv,
			     D2_delta_pv_sv,
			     momentum_jet, momentum_jet};
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  filling   histogramm   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      unsigned fill = effsam;
      /* bool isDataDrivenBgk= false;
	 if (samples[sam].isData() && tightFail_sFR     && single_fake)     isDataDrivenBgk= true;
	 if (samples[sam].isData() && loose_lepton_dFR  && Double_fake)     isDataDrivenBgk= true;
	 bool isDataYield= false;
	 if (samples[sam].isData() && !tightFail_sFR     && single_fake)     isDataYield= true;
	 if (samples[sam].isData() && tight_lepton_dFR  && Double_fake)      isDataYield= true;
	 if (isDataDrivenBgk) fill = nSamples_eff;
	 if (isDataYield)     fill = 0;
	 if (isDataYield)     scal = 1;
	 if (isDataYield)     continue;*/

     
 

      int channel_bin = -1;
      channel_bin = SR_channel+1;
      if (isSRRun && channel_bin == -1 ) continue;

             
      // ------------------- Histo SR
      if (SR_channel <= 2) {
	//std::cout<<"sr channel: "<< SR_channel<<"  channel bin "<<channel_bin<< "  "<< bin_SR_muonCoupling<<"  check: ("<<flavors_3l[0]<<","<<flavors_3l[1]<<","<<flavors_3l[2]<<")"<< std::endl;
	//Histos[0][SR_channel][cut_bin][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	//Histos[0][6][cut_bin][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_0)      Histos[0][SR_channel][0][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_1)      Histos[0][SR_channel][1][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_2)      Histos[0][SR_channel][2][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_3)      Histos[0][SR_channel][3][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_4)      Histos[0][SR_channel][4][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_5)      Histos[0][SR_channel][5][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_final)  Histos[0][SR_channel][6][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_0)      Histos[0][6][0][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_1)      Histos[0][6][1][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_2)      Histos[0][6][2][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_3)      Histos[0][6][3][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_4)      Histos[0][6][4][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_5)      Histos[0][6][5][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);
	if (selection_final)  Histos[0][6][6][fill] -> Fill(static_cast<double>(bin_SR_muonCoupling), scal);

	
      }
      if (SR_channel > 2) {
	//std::cout<<"sr channel: "<< SR_channel<<"  channel bin "<<channel_bin<< "  "<< bin_SR_eleCoupling<<"  check: ("<<flavors_3l[0]<<","<<flavors_3l[1]<<","<<flavors_3l[2]<<")"<< std::endl;
	if (selection_0)      Histos[0][SR_channel][0][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_1)      Histos[0][SR_channel][1][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_2)      Histos[0][SR_channel][2][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_3)      Histos[0][SR_channel][3][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_4)      Histos[0][SR_channel][4][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_5)      Histos[0][SR_channel][5][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_final)  Histos[0][SR_channel][6][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_0)      Histos[0][7][0][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_1)      Histos[0][7][1][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_2)      Histos[0][7][2][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_3)      Histos[0][7][3][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_4)      Histos[0][7][4][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_5)      Histos[0][7][5][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	if (selection_final)  Histos[0][7][6][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	//Histos[0][SR_channel][cut_bin][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
	//Histos[0][7][cut_bin][fill] -> Fill(static_cast<double>(bin_SR_eleCoupling), scal);
      }
      // ------------------- Histo cut flow
      //Histos[1][SR_channel][0][fill]->Fill(static_cast<double>(cut_bin+1), scal);
      //Histos[1][SR_channel][cut_bin][fill]->Fill(static_cast<double>(cut_bin+1), scal);

      if (selection_0)      Histos[1][SR_channel][0][fill] -> Fill(static_cast<double>(1), scal);
      if (selection_1)      Histos[1][SR_channel][0][fill] -> Fill(static_cast<double>(2), scal);
      if (selection_2)      Histos[1][SR_channel][0][fill] -> Fill(static_cast<double>(3), scal);
      if (selection_3)      Histos[1][SR_channel][0][fill] -> Fill(static_cast<double>(4), scal);
      if (selection_4)      Histos[1][SR_channel][0][fill] -> Fill(static_cast<double>(5), scal);
      if (selection_5)      Histos[1][SR_channel][0][fill] -> Fill(static_cast<double>(6), scal);
      if (selection_final)  Histos[1][SR_channel][0][fill] -> Fill(static_cast<double>(7), scal);
      
      if (SR_channel <= 2){
	//Histos[1][6][0][fill]->Fill(static_cast<double>(cut_bin+1), scal);
	if (selection_0)      Histos[1][6][0][fill] -> Fill(static_cast<double>(1), scal);
	if (selection_1)      Histos[1][6][0][fill] -> Fill(static_cast<double>(2), scal);
	if (selection_2)      Histos[1][6][0][fill] -> Fill(static_cast<double>(3), scal);
	if (selection_3)      Histos[1][6][0][fill] -> Fill(static_cast<double>(4), scal);
	if (selection_4)      Histos[1][6][0][fill] -> Fill(static_cast<double>(5), scal);
	if (selection_5)      Histos[1][6][0][fill] -> Fill(static_cast<double>(6), scal);
	if (selection_final)  Histos[1][6][0][fill] -> Fill(static_cast<double>(7), scal);
      }
      if (SR_channel > 2){
	if (selection_0)      Histos[1][7][0][fill] -> Fill(static_cast<double>(1), scal);
	if (selection_1)      Histos[1][7][0][fill] -> Fill(static_cast<double>(2), scal);
	if (selection_2)      Histos[1][7][0][fill] -> Fill(static_cast<double>(3), scal);
	if (selection_3)      Histos[1][7][0][fill] -> Fill(static_cast<double>(4), scal);
	if (selection_4)      Histos[1][7][0][fill] -> Fill(static_cast<double>(5), scal);
	if (selection_5)      Histos[1][7][0][fill] -> Fill(static_cast<double>(6), scal);
	if (selection_final)  Histos[1][7][0][fill] -> Fill(static_cast<double>(7), scal);
      }     
      // ------------------- all the other histograms
      for(int numero_histo = 0; numero_histo < nDist; ++numero_histo){
	//Histos[numero_histo][SR_channel][cut_bin][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);

	if (selection_0) Histos[numero_histo][SR_channel][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_1) Histos[numero_histo][SR_channel][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_2) Histos[numero_histo][SR_channel][2][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_3) Histos[numero_histo][SR_channel][3][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_4) Histos[numero_histo][SR_channel][4][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_5) Histos[numero_histo][SR_channel][5][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	if (selection_final) Histos[numero_histo][SR_channel][6][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);


	
	if (SR_channel <= 2){
	  if (selection_0) Histos[numero_histo][6][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_1) Histos[numero_histo][6][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_2) Histos[numero_histo][6][2][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_3) Histos[numero_histo][6][3][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_4) Histos[numero_histo][6][4][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_5) Histos[numero_histo][6][5][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_final) Histos[numero_histo][6][6][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  // Histos[numero_histo][6][cut_bin][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	}
	if (SR_channel > 2)  {
	  if (selection_0) Histos[numero_histo][7][0][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_1) Histos[numero_histo][7][1][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_2) Histos[numero_histo][7][2][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_3) Histos[numero_histo][7][3][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_4) Histos[numero_histo][7][4][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_5) Histos[numero_histo][7][5][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  if (selection_final) Histos[numero_histo][7][6][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	  //Histos[numero_histo][7][cut_bin][fill]->Fill(TMath::Min(values[numero_histo], maxBinC[numero_histo]), scal);
	}
      }//end histo
      
    }//end loop over the entries
  }//loop over samples

  std::cout<<"multiboson: "<< Histos[2][6][0][24]->GetSumOfWeights()<<std::endl;
  std::cout<<"dy: "<< Histos[2][6][0][21]->GetSumOfWeights()<<std::endl;


  // TH1D* Histos[nDist][nChannel][nCat][nSamples_eff +1];

  // THIS IS THE UNBLIND PLOT ===>IT HAS TO SILENT IN THE PLOTTING!!!!!!!!!!!!!!!!!!
  TH1D* dataYields[nDist][nChannel][nCat];
  for(unsigned dist = 0; dist < nDist; ++dist){
    for(unsigned cat = 0; cat < nCat; ++cat){
      for(int cha = 0; cha < nChannel; ++cha){               
	if (isSRRun) dataYields[dist][cha][cat] = (TH1D*) Histos[dist][cha][cat][nSamples_signal+1]->Clone();
	if (isCRRun) dataYields[dist][cha][cat] = (TH1D*) Histos[dist][cha][cat][0]->Clone();
      }
    }
  }

  TH1D* bkgYields[nDist][nChannel][nCat][nSamples_eff - nSamples_signal]; //change to nSamples_eff if sig is removed
  for(unsigned dist = 0; dist < nDist; ++dist){
    for(unsigned cat = 0; cat < nCat; ++cat){
      for(int cha = 0; cha < nChannel; ++cha){
	
	for(unsigned effsam1 = nSamples_signal+1; effsam1 < nSamples_eff +1 ; ++effsam1){	  
	  // put_at_zero(*&Histos[dist][cha][cat][effsam1]);
	  
	  bkgYields[dist][cha][cat][effsam1 -nSamples_signal-1] = (TH1D*) Histos[dist][cha][cat][effsam1]->Clone();
	  
	  if(effsam1 > nSamples_signal+1 && effsam1 < nSamples_eff){	  
	    if (isSRRun) dataYields[dist][cha][cat]->Add(bkgYields[dist][cha][cat][effsam1 -nSamples_signal-1]);
	  }	  
	}
      }
    }
  }

  


  /*

    const std::vector< std::string > uncNames = {"JEC_2016", "uncl", "scale", "pileup", "bTag_udsg_2016", "bTag_bc_2016", "isr", "fsr", "prefiring", "WZ_extrapolation",
    "lepton_reco", "muon_id_stat_2016", "electron_id_stat_2016", "lepton_id_syst", "pdf", "scaleXsec", "pdfXsec"};



    const std::string systNames[8]= { "lumi", "pdfAcc", "JEC","pdf","pu",   "stasignal","statDY","statttbar",	"statWJets",	"statmultiboson", "statXgamma","statTTTX", "statNonPrompt"};
    const std::string systDist[8]= {"lnN","lnN", "lnN",	"lnN",	"lnN", 	"lnN", "lnN", "lnN","lnN"};
    const std::string bgkNames[6]= {"DY",  	"ttbar",	"WJets",	"multiboson", 	"Xgamma",    	"TTTX"};



    lumi	lnN	1.025	-	1.025	-	-	1.025	-
    pdfAcc	lnN	1.00184	-	-	-	-	-	-
    JEC	shape	1	1	1	1	1	1	-
    metUncl	shape	1	1	1	1	1	1	-
    scale_elebin4	shape	1	1	1	1	1	1	-
    pdf	shape	1	1	1	1	1	1	-
    pu	shape	1	1	1	1	1	1	-
    btagSF	shape	1	1	1	1	1	1	-
    id_eff	shape	1	1	1	1	1	1	-
    trigeff	shape	1	1	1	1	1	1	-
    fakeEWK	shape	-	-	-	-	-	-	1
    ZZmt	shape	-	1	-	-	-	-	-
    scaleAcc	shape	1	-	-	-	-	-	-
    lifetime	lnN	1.18522	-	-	-	-	-	-
    statSig4	lnN	1	-	-	-	-	-	-
    statZZH4	lnN	-	8.79488	-	-	-	-	-
    stattriboson4	lnN	-	-	1	-	-	-	-
    statWZ4	lnN	-	-	-	1.35968	-	-	-
    statXgamma4	lnN	-	-	-	-	1	-	-
    statTTX4	lnN	-	-	-	-	-	1	-
    statnonPrompt4	lnN	-	-	-	-	-	-	2.89821
    extraZZH	lnN	-	1.1	-	-	-	-	-
    extratriboson	lnN	-	-	1.5	-	-	-	-
    extraWZ	lnN	-	-	-	1.094	-	-	-
    extraXgamma	lnN	-	-	-	-	1.15	-	-
    extraTTX	lnN	-	-	-	-	-	1.5	-
    extranonPrompt	lnN	-	-	-	-	-	-	1.3


  */

  

  TH1D* signals[nSamples_signal];

  for(unsigned dist = 0; dist < nDist; ++dist){
    for(unsigned cat = 0; cat < nCat; ++cat){
      for(int cha = 0; cha < nChannel; ++cha){               
	for (unsigned signal_sample = 0; signal_sample< nSamples_signal; signal_sample++){
	  signals[signal_sample] = (TH1D*) Histos[dist][cha][cat][signal_sample+1]->Clone() ;     
	}

      
	plotDataVSMC(cat,cha,dist,
		     dataYields[dist][cha][cat], bkgYields[dist][cha][cat],
		     eff_names,nSamples_eff -  nSamples_signal -1 ,
		     catNames[cat], channelNames[cha], channelNames[cha]+"_"+ Histnames_ossf[dist]+"_"+catNames[cat],
		     true,
		     2, true, signals,  sigNames, nSamples_signal, false);
      }
    }//end cat
  }//end histo
 

  /*
  // da qui e' la roba per le data card

      
  for (int ii = 0; ii < nSR; ii++){// loop over bin
	
  double bgkield[6]= {0.,0.,0.,0.,0.,0.};
  bgkield[0] = Histos[dist][cat][63] -> GetBinContent(ii+1);
  bgkield[1] = Histos[dist][cat][64] -> GetBinContent(ii+1);
  bgkield[2] = Histos[dist][cat][65] -> GetBinContent(ii+1);
  bgkield[3] = Histos[dist][cat][66] -> GetBinContent(ii+1);
  bgkield[4] = Histos[dist][cat][67] -> GetBinContent(ii+1);
  bgkield[5] = Histos[dist][cat][68] -> GetBinContent(ii+1);


  std::vector<std::vector<double> > systUnc (7,vector<double>(7,0.));

  for (int i =0; i < 7; i++){
  for (int j =0; j < 7; j++){
  if (i != j)  systUnc[i][j] = 0;
  else{	    
  if (i == 1 && Histos[dist][cat][63] -> GetBinContent(ii+1) !=0) systUnc[i][j] = (1+ (Histos[dist][cat][63] -> GetBinError(ii+1)) / (Histos[dist][cat][63] -> GetBinContent(ii+1)));
  if (i == 1 && systUnc[i][j] >= 2)                               systUnc[i][j] = 1.99;
  if (i == 1 && Histos[dist][cat][63] -> GetBinContent(ii+1) ==0) systUnc[i][j] = 1;

  if (i == 2 && Histos[dist][cat][64] -> GetBinContent(ii+1) !=0) systUnc[i][j] = (1+ (Histos[dist][cat][64] -> GetBinError(ii+1)) / (Histos[dist][cat][64] -> GetBinContent(ii+1)));
  if (i == 2 && Histos[dist][cat][64] -> GetBinContent(ii+1) ==0) systUnc[i][j] = 1;

  if (i == 2 && Histos[dist][cat][65] -> GetBinContent(ii+1) !=0) systUnc[i][j] = (1+ (Histos[dist][cat][64] -> GetBinError(ii+1)) / (Histos[dist][cat][64] -> GetBinContent(ii+1)));
  if (i == 2 && systUnc[i][j] >= 2)                               systUnc[i][j] = 1.99;
  if (i == 2 && Histos[dist][cat][65] -> GetBinContent(ii+1) ==0) systUnc[i][j] = 1;

  if (i == 3 && Histos[dist][cat][65] -> GetBinContent(ii+1) !=0) systUnc[i][j] = (1+ (Histos[dist][cat][65] -> GetBinError(ii+1)) / (Histos[dist][cat][65] -> GetBinContent(ii+1)));
  if (i == 3 && systUnc[i][j] >= 2)                               systUnc[i][j] = 1.99;
  if (i == 3 && Histos[dist][cat][65] -> GetBinContent(ii+1) ==0) systUnc[i][j] = 1;

  if (i == 4 && Histos[dist][cat][66] -> GetBinContent(ii+1) !=0) systUnc[i][j] = (1+ (Histos[dist][cat][66] -> GetBinError(ii+1)) / (Histos[dist][cat][66] -> GetBinContent(ii+1)));
  if (i == 4 && systUnc[i][j] >= 2)                               systUnc[i][j] = 1.99;
  if (i == 4 && Histos[dist][cat][66] -> GetBinContent(ii+1) ==0) systUnc[i][j] = 1;

  if (i == 5 && Histos[dist][cat][67] -> GetBinContent(ii+1) !=0) systUnc[i][j] = (1+ (Histos[dist][cat][67] -> GetBinError(ii+1)) / (Histos[dist][cat][67] -> GetBinContent(ii+1)));
  if (i == 5 && systUnc[i][j] >= 2)                               systUnc[i][j] = 1.99;
  if (i == 5 && Histos[dist][cat][67] -> GetBinContent(ii+1) ==0) systUnc[i][j] = 1;

  if (i == 6 && Histos[dist][cat][68] -> GetBinContent(ii+1) !=0) systUnc[i][j] = (1+ (Histos[dist][cat][68] -> GetBinError(ii+1)) / (Histos[dist][cat][68] -> GetBinContent(ii+1)));
  if (i == 6 && systUnc[i][j] >= 2)                               systUnc[i][j] = 1.99;
  if (i == 6 && Histos[dist][cat][68] -> GetBinContent(ii+1) ==0) systUnc[i][j] = 1;

  }
  }
  }

  // IMPORTANTE   da qui e' la roba per le data card
  for (unsigned signal_sample = 0; signal_sample< 62; signal_sample++){	  
  if(signals[signal_sample] -> GetBinContent(ii+1) !=0) systUnc[0][0]= (1+ (signals[signal_sample] -> GetBinError(ii+1)) / (signals[signal_sample] -> GetBinContent(ii+1)));
  if(systUnc[0][0] >=2 )                                systUnc[0][0] = 1.99;	  
  if(signals[signal_sample] -> GetBinContent(ii+1) ==0) systUnc[0][0] = 1;


	  
  printDataCard(   dataYields[dist][cat] ->GetBinContent(ii+1), 
  signals[signal_sample]->GetBinContent(ii+1),
  sigNamespp[signal_sample],
  bgkield,
  6,
  bgkNames,
  systUnc, 7, systNames,systDist,
  sigNamespp[signal_sample]+"_bin"+std::to_string(ii+1)+".txt",
  false, sigNamespp[signal_sample]+"_bin"+std::to_string(ii+1), ii+1);


  }// end loop signal


  }//loop bin
  */
 


 
  

}//END ANALIUSI



//___________________________________________________________________
double Analysis_mc::pu_weight ( TH1D *histo, double numberInteractions){
  double nI = numberInteractions;   
  double factore=0;
  factore = histo->GetBinContent(histo->FindBin(nI));
  return factore;
}

 
