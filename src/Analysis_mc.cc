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
  fChain->SetBranchAddress("_updated_ecalBadCalibFilter", &_updated_ecalBadCalibFilter, &b__updated_ecalBadCalibFilter);  
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

  cout<<"in analisi"<<endl;
  cout<<"---------------------------"<<endl;   
  //setTDRStyle();
  if(systdir<0) {
    std::cout << " >>> Dummy message (to avoid warnings): systdir " << systdir << std::endl;
  }

  TFile *fout = new TFile(outfilename.Data(), "recreate");
  
  // ------------ pile up -----------------------------------------------//
  TH1D *pileUpWeight[1];

  TFile hfile_pu("/user/mvit/CMSSW_9_4_4/src/HNL_analysis/PU/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root");
  pileUpWeight[0] = (TH1D*)hfile_pu.Get("puw_Run2016Inclusive_central");
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
  // ------------   run over samples -----------------------------------------------//  
  for(int sam = 0; sam < samples.size(); ++sam){
    initSample(jaar,samples[sam]);
    //check consistency
    cout<<"sample initialized: --> "<<endl;
    cout<<"fileName: "<<samples[sam].getFileName()<<"  process name: "<< samples[sam].getProcessName()<< "   xsec: "<< samples[sam].getXSec()<<endl;  
    if (samples[sam].isData()) cout<<"is Data"<<endl;
    if (samples[sam].isMC()) cout<<"is MC"<<endl;
    if (samples[sam].isNewPhysicsSignal()) cout<<"is signal"<<endl;
    
    double progress = 0; 	//For printing progress bar 
    // ------------   run over entries -----------------------------------------------//  
    for (Long64_t it = 0; it < nEntries/100; ++it){
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
      double scal = 0;
      scal = scale*_weight * pu_weight(*&pileUpWeight[0],_nTrueInt);
      bwght=1.;

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
      kind[0]=-1;   
      skip_event[0]= -1;
      for (int i =0; i < 3; i++){
	lepton_reco[i].SetPtEtaPhiE(0.,0.,0.,0.);
	lepton_transv[i].SetPtEtaPhiE(0.,0.,0.,0.);
	pair[i].SetPtEtaPhiE(0.,0.,0.,0.);
	flavors_3l[i]=0;
	charge_3l[i]=0;	
      }
      unsigned        l1=0;
      unsigned        l2=0;
      unsigned        l3=0;
      TLorentzVector  v4l1;
      TLorentzVector  v4l2;
      TLorentzVector  v4l3;

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
      //select leptons
      //      std::cout<<"after real"<<std::endl;
	      
      const unsigned lCount = selectLepConeCorr(ind);
      //     std::cout<<"after selct lep"<<std::endl;

      if (lCount < 3) continue;

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

      //check how many displaced there are (displaced --> dxy, common vertex, FO, no l1)
      unsigned displacedC = 0;
      std::vector<TLorentzVector> lepV_displaced;
      std::vector<int> charge_displaced;
      std::vector<unsigned> temp_index;
      
    
      int index_to_use_for_l2_l3[2]={0,0};
      //find the right OS pair with min invariant mass
      int min_test= 9999;
      int min_mass=999; 
      for(unsigned l = 0; l < lCount; ++l){
	for(unsigned j = l+1; j < lCount; ++j){
	  if(!lepIsDisplaced(ind[l] , ind_new_leading, ind)) continue;
	  if(!lepIsDisplaced(ind[j] , ind_new_leading, ind)) continue;
	  if (_lCharge[ind[l]] == _lCharge[ind[j]]) continue;
	  ++displacedC;
	  TLorentzVector temp_displaced1;
	  TLorentzVector temp_displaced2;
	  temp_displaced1.SetPtEtaPhiE(_lPt[ind[l]],_lEta[ind[l]], _lPhi[ind[l]], _lE[ind[l]]);
	  temp_displaced2.SetPtEtaPhiE(_lPt[ind[j]],_lEta[ind[j]], _lPhi[ind[j]], _lE[ind[j]]);
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
	}//end loop2
      }//end loop1
      if (displacedC< 2) continue;
     
      // ------------ changing all the lep info -----------------------------------------------//
      l1=ind_new_leading;
      l2=index_to_use_for_l2_l3[0];
      l3=index_to_use_for_l2_l3[1];
      v4l1.SetPtEtaPhiE(_lPt[l1],_lEta[l1], _lPhi[l1], _lE[l1]);
      v4l2.SetPtEtaPhiE(_lPt[l2],_lEta[l2], _lPhi[l2], _lE[l2]);
      v4l3.SetPtEtaPhiE(_lPt[l3],_lEta[l3], _lPhi[l3], _lE[l3]);
      int index_l2l3= l2l3_vertex_variable (l2,l3);

      _vertex_X=_vertices[index_l2l3][1];
      _vertex_Y=_vertices[index_l2l3][2];
      _vertex_Z=_vertices[index_l2l3][3];
      _vertex_chi2=_vertices[index_l2l3][11];
      _vertex_normchi2= _vertices[index_l2l3][11]/_vertices[index_l2l3][10];
      _vertex_ndf =_vertices[index_l2l3][10];
      
      std::cout<<"vetrex: "<<_vertex_X<<" "<<_vertex_Y<<" "<<_vertex_Z<<" "<<_vertex_ndf<<"    ------> "<<_vertices_os[index_l2l3][1]<< std::endl;
      // ------------ ==================== -----------------------------------------------//
     

      if (!_passTrigger_1l) continue;

      
    }//end loop over the entries
    
  }//loop over samples
  
}//END ANALIUSI



//___________________________________________________________________
double Analysis_mc::pu_weight ( TH1D *histo, double numberInteractions){
  double nI = numberInteractions;   
  double factore=0;
  factore = histo->GetBinContent(histo->FindBin(nI));
  return factore;
}
