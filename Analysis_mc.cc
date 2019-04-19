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
#include "Analysis_mc.h"
#include "TApplication.h"
#include "TColor.h"



//include C++ library classes
using std::cout;
using std::endl;
using std::flush;
using std::ofstream;

//include other parts of the code
#include "tdrstyle.h"
#include "plotCode_new.h"
//#include "Selection.h"

// For b-tagging SFs and variations thereof
#include "BTagCalibrationStandalone.h"
//#include "BTagCalibrationStandalone.cpp"


using namespace std;




static Double_t pigreco= TMath::ACos(-1);
ClassImp(Analysis_mc)

//_______________________________________________________default constructor_____
Analysis_mc::Analysis_mc():TObject()

{
}
//_______________________________________________________ constructor_____
Analysis_mc::Analysis_mc(unsigned jaar, int selezione, std::string FileNameTree_in) : TObject()

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

  if(selezione==0) 
    std::cout << " >>> Dummy message (to avoid warnings): selezione " << selezione
	      << ",  FileNameTree_in " << FileNameTree_in.c_str() << std::endl;
}
//________________________________________________________________distruttore_____
Analysis_mc::~Analysis_mc()	 {
  // destructor
}

//_______________________________________________________ constructor_____
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



//==================================================================
void Analysis_mc::analisi(int selezione, int num_histo_kin,
			  TString outfilename,
			  int systcat, int systdir
                          ) {
    
  cout<<"in analisi"<<endl;
  cout<<"---------------------------"<<endl;   
  //setTDRStyle();
  if(selezione==0 || num_histo_kin==0 || systdir<0) {
    std::cout << " >>> Dummy message (to avoid warnings): selezione " << selezione
	      << ", num_histo_kin " << num_histo_kin << ", systdir " << systdir << std::endl;
  }

  TFile *fout = new TFile(outfilename.Data(), "recreate");

  TH1D *pileUpWeight[1];    
  TFile hfile_pu("/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/ANALYSIS/20190318_MartinasCode/samples.noSync/2016/puWeights_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_Summer16.root");
  pileUpWeight[0] = (TH1D*)hfile_pu.Get("puw_Run2016Inclusive_central");

  // b-tagging working points (DeepCsv_b + DeepCsv_bb)
  double btagCuts[3][3];
  // selected WP (0: loose; 1: medium; 2: tight)
  BTagEntry::OperatingPoint bwp = BTagEntry::OP_LOOSE;    // = 0
  //BTagEntry::OperatingPoint bwp = BTagEntry::OP_MEDIUM; // = 1
  //BTagEntry::OperatingPoint bwp = BTagEntry::OP_TIGHT;  // = 2

  //  - 2016
  btagCuts[0][0] = 0.2217; // loose
  btagCuts[0][1] = 0.6321; // medium
  btagCuts[0][2] = 0.8953; // tight
  //  - 2017
  btagCuts[1][0] = 0.1522; // loose
  btagCuts[1][1] = 0.4941; // medium
  btagCuts[1][2] = 0.8001; // tight
  //  - 2018
  btagCuts[2][0] = 0.1241; // loose
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


  const double glugluToZZkFactor = 1;
  const double WZSF = 1;
  const double ZZSF = 1;
  const double XgammaSF = 1;
  const double low_coupling= 0.0001;
  const double low_coupling_2= 0.00001;
  const double high_coupling = 0.001;
  const double coupling_factor_low= low_coupling / 0.0001;
  const double coupling_factor_low_2= low_coupling_2 / 0.0001;
 
 
  // const int nSamples= 109;
  // const int nSamples_eff = 69;
  // const int nSamples_signal=62;
  const int nSamples = 2;
  const int nSamples_eff = 2;
  const int nSamples_signal = 2;

  const double fl = 1;
    

  const TString fileList[nSamples] = {  "HeavyNeutrino_trilepton_M-4_V-0.00290516780927_e_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.00290516780927_e_massiveAndCKM_LO.root"
    //"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root",  // dummy (???)
					//
					//"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_Summer16.root"
					//"HeavyNeutrino_trilepton_M-5_V-0.01_mu_massiveAndCKM_LO.root"
					//"HeavyNeutrino_trilepton_M-2_V-0.0110905365064_mu_massiveAndCKM_LO.root",
					//"HeavyNeutrino_trilepton_M-4_V-0.0013_mu_massiveAndCKM_LO.root"
					//"HeavyNeutrino_trilepton_M-10_V-0.00244948974278_mu_massiveAndCKM_LO.root"
					// "HeavyNeutrino_trilepton_M-1_V-0.0949736805647_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-1_V-0.212367605816_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-2_V-0.0110905365064_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-2_V-0.0248394846967_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-3_V-0.00707813534767_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-4_V-0.00290516780927_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-5_V-0.00145602197786_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-6_V-0.00202484567313_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-8_V-0.00151327459504_e_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-10_V-0.000756967634711_e_massiveAndCKM_LO.root",
					// //
					// "HeavyNeutrino_trilepton_M-1_V-0.0949736805647_mu_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-2_V-0.0110905365064_mu_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-3_V-0.00707813534767_mu_massiveAndCKM_LO.root",
					// "HeavyNeutrino_trilepton_M-4_V-0.00290516780927_mu_massiveAndCKM_LO.root"
  };
  const TString names[nSamples]   =  {  "total",
					//
					"TT"
					//"M-2_V-0.0110905365064_mu",
					//"M-4_V-0.0013_mu",
					//"M-10_V-0.00244948974278_mu"
					// "M-1_V-0.0949736805647_e",
					// "M-1_V-0.212367605816_e",
					// "M-2_V-0.0110905365064_e",
					// "M-2_V-0.0248394846967_e",
					// "M-3_V-0.00707813534767_e",
					// "M-4_V-0.00290516780927_e",
					// "M-5_V-0.00145602197786_e",
					// "M-6_V-0.00202484567313_e",
					// "M-8_V-0.00151327459504_e",
					// "M-10_V-0.000756967634711_e",
					// //
					// "M-1_V-0.0949736805647_e",
					// "M-2_V-0.0110905365064_e",
					// "M-3_V-0.00707813534767_e",
					// "M-4_V-0.00290516780927_e"
  };

  const TString eff_names[nSamples_eff +1 ] = { "total",
						//
						"TT"
						//"M-2_V-0.0110905365064_mu",
						//"M-4_V-0.0013_mu",
						//"M-10_V-0.00244948974278_mu",
						// "M-1_V-0.0949736805647_e",
						// "M-1_V-0.212367605816_e",
						// "M-2_V-0.0110905365064_e",
						// "M-2_V-0.0248394846967_e",
						// "M-3_V-0.00707813534767_e",
						// "M-4_V-0.00290516780927_e",
						// "M-5_V-0.00145602197786_e",
						// "M-6_V-0.00202484567313_e",
						// "M-8_V-0.00151327459504_e",
						// "M-10_V-0.000756967634711_e",
						// //
						// "M-1_V-0.0949736805647_e",
						// "M-2_V-0.0110905365064_e",
						// "M-3_V-0.00707813534767_e",
						// "M-4_V-0.00290516780927_e",
						"no prompt"
  };

  const double xSections[nSamples]= {0,        
				     //
				     87.315,
				     //fl*5.287e-01,
				     //fl*6.750e-03,
				     //2.521e-02 // M=10 GeV
				     // fl*fl*3.94e+01,
				     // fl*fl*fl*1.97e+02,
				     // fl*5.37e-01,
				     // fl*fl*2.69,
				     // fl*2.06e-01,
				     // fl*3.41e-02,
				     // fl*8.60e-03,
				     // fl*1.66e-02,
				     // fl*9.42e-03,
				     // fl*2.37e-03,
				     // //
				     // fl*fl*3.94e+01,
				     // fl*5.37e-01,
				     // fl*2.06e-01,
				     // fl*3.41e-02
  };

  /*
  const TString fileList[nSamples] = {  "ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root",

					"HeavyNeutrino_trilepton_M-2_V-0.00244948974278_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.00282842712475_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.00316227766017_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.004472135955_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.00547722557505_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.00707106781187_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.00836660026534_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.0141421356237_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.0173205080757_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.01_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.022360679775_mu_massiveAndCKM_LO.root",
					
					"HeavyNeutrino_trilepton_M-3_V-0.00244948974278_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.00282842712475_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.00316227766017_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.004472135955_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.00547722557505_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.00707106781187_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.00836660026534_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.0141421356237_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.0173205080757_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.01_mu_massiveAndCKM_LO.root",

					"HeavyNeutrino_trilepton_M-4_V-0.00244948974278_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.00282842712475_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.00316227766017_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.004472135955_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.00547722557505_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.00707106781187_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.00836660026534_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.0141421356237_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.01_mu_massiveAndCKM_LO.root",
					
					"HeavyNeutrino_trilepton_M-5_V-0.00244948974278_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.00282842712475_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.00316227766017_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.004472135955_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.00547722557505_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.00707106781187_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.00836660026534_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.01_mu_massiveAndCKM_LO.root",
					
					"HeavyNeutrino_trilepton_M-8_V-0.00244948974278_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-8_V-0.00282842712475_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-8_V-0.00316227766017_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-8_V-0.004472135955_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-8_V-0.00547722557505_mu_massiveAndCKM_LO.root",

					"HeavyNeutrino_trilepton_M-10_V-0.00244948974278_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-10_V-0.00282842712475_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-10_V-0.00316227766017_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-10_V-0.004472135955_mu_massiveAndCKM_LO.root",

					"HeavyNeutrino_trilepton_M-15_V-0.00244948974278_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-15_V-0.00282842712475_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-15_V-0.00316227766017_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-15_V-0.004472135955_mu_massiveAndCKM_LO.root",


					"HeavyNeutrino_trilepton_M-20_V-1.93390796058e-05_mu_massiveAndCKM_LO.root",


					"HeavyNeutrino_trilepton_M-1_V-0.0949736805647_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-1_V-0.212367605816_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.0110905365064_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-2_V-0.0248394846967_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-3_V-0.00707813534767_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-4_V-0.00290516780927_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-5_V-0.00145602197786_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-6_V-0.00202484567313_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-8_V-0.00151327459504_mu_massiveAndCKM_LO.root",
					"HeavyNeutrino_trilepton_M-10_V-0.000756967634711_mu_massiveAndCKM_LO.root",



					"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",

					"TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", "TTJets_SingleLeptFromT_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root", "TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root ",
					"ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1.root" ,"ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root" ,"ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1.root" ,"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root" ,"ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1.root",
					
					"WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root",


					"ZH_HToZZ_4LFilter_M125_13TeV_powheg2-minlo-HZJ_JHUgenV6_pythia8.root", "VBF_HToZZTo4L_M125_13TeV_powheg2_JHUgenV6_pythia8.root",

					"WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "WWZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root","ZZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root","WZZ_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "WZG_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root","WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "WGGJets_TuneCUETP8M1_13TeV_madgraphMLM_pythia8.root",
					
					"GluGluToContinToZZTo2e2mu_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo2e2tau_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo2mu2tau_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo4e_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo4mu_13TeV_MCFM701_pythia8.root", "GluGluToContinToZZTo4tau_13TeV_MCFM701_pythia8.root","WWTo2L2Nu_13TeV-powheg.root", "WWTo2L2Nu_DoubleScattering_13TeV-pythia8.root","ZZTo4L_13TeV_powheg_pythia8.root ", "WZTo3LNu_mllmin01_13TeV-powheg-pythia8_ext1.root","WplusH_HToZZTo4L_M125_13TeV_powheg2-minlo-HWJ_JHUgenV6_pythia8.root",


					"WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root","ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root", 

					"TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8.root", "TTGG_0Jets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8.root", "TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root "    , "TTWJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8.root","TTWJetsToQQ_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8_0.root", "TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "TTZToLL_M-1to10_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","tZq_ll_4f_13TeV-amcatnlo-pythia8.root","TTTT_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root", "TTWW_TuneCUETP8M2T4_13TeV-madgraph-pythia8.root", "TTWZ_TuneCUETP8M2T4_13TeV-madgraph-pythia8.root", "TTZZ_TuneCUETP8M2T4_13TeV-madgraph-pythia8.root", "ST_tWll_5f_LO_13TeV_MadGraph_pythia8.root"

  };
       
  const TString names[nSamples]   =  {  "total",      
					"M-2_V-0.00244948974278_mu",
					"M-2_V-0.00282842712475_mu",
					"M-2_V-0.00316227766017_mu",
					"M-2_V-0.004472135955_mu",
					"M-2_V-0.00547722557505_mu",
					"M-2_V-0.00707106781187_mu",
					"M-2_V-0.00836660026534_mu",
					"M-2_V-0.0141421356237_mu",
					"M-2_V-0.0173205080757_mu",
					"M-2_V-0.01_mu",
					"M-2_V-0.022360679775_mu",
					"M-3_V-0.00244948974278_mu",
					"M-3_V-0.00282842712475_mu",
					"M-3_V-0.00316227766017_mu",
					"M-3_V-0.004472135955_mu",
					"M-3_V-0.00547722557505_mu",
					"M-3_V-0.00707106781187_mu",
					"M-3_V-0.00836660026534_mu",
					"M-3_V-0.0141421356237_mu",
					"M-3_V-0.0173205080757_mu",
					"M-3_V-0.01_mu",
					"M-4_V-0.00244948974278_mu",
					"M-4_V-0.00282842712475_mu",
					"M-4_V-0.00316227766017_mu",
					"M-4_V-0.004472135955_mu",
					"M-4_V-0.00547722557505_mu",
					"M-4_V-0.00707106781187_mu",
					"M-4_V-0.00836660026534_mu",
					"M-4_V-0.0141421356237_mu",
					"M-4_V-0.01_mu",
					"M-5_V-0.00244948974278_mu",
					"M-5_V-0.00282842712475_mu",
					"M-5_V-0.00316227766017_mu",
					"M-5_V-0.004472135955_mu",
					"M-5_V-0.00547722557505_mu",
					"M-5_V-0.00707106781187_mu",
					"M-5_V-0.00836660026534_mu",
					"M-5_V-0.01_mu",
					"M-8_V-0.00244948974278_mu",
					"M-8_V-0.00282842712475_mu",
					"M-8_V-0.00316227766017_mu",
					"M-8_V-0.004472135955_mu",
					"M-8_V-0.00547722557505_mu",
					"M-10_V-0.00244948974278_mu",
					"M-10_V-0.00282842712475_mu",
					"M-10_V-0.00316227766017_mu",
					"M-10_V-0.004472135955_mu",
					"M-15_V-0.00244948974278_mu",
					"M-15_V-0.00282842712475_mu",
					"M-15_V-0.00316227766017_mu",
					"M-15_V-0.004472135955_mu",
					"M-20_V-1.93390796058e-05_mu",
					
					"M-1_V-0.0949736805647_mu",
					"M-1_V-0.212367605816_mu",
					"M-2_V-0.0110905365064_mu",
					"M-2_V-0.0248394846967_mu",
					"M-3_V-0.00707813534767_mu",
					"M-4_V-0.00290516780927_mu",
					"M-5_V-0.00145602197786_mu",
					"M-6_V-0.00202484567313_mu",
					"M-8_V-0.00151327459504_mu",
					"M-10_V-0.000756967634711_mu",


    
					"DY",    "DY",   
					"TTbar","TTbar","TTbar","TTbar", "TTbar","TTbar", "TTbar","TTbar",
					"WJets",
					"multiboson","multiboson",
					"multiboson", "multiboson", "multiboson","multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson", "multiboson","multiboson",
					"X+#gamma","X+#gamma",      
					"TT/T + X","TT/T + X","TT/T + X", "TT/T + X","TT/T + X","TT/T + X","TT/T + X","TT/T + X","TT/T + X","TT/T + X","TT/T + X","TT/T + X","TT/T + X"
				
					
  };
  const TString eff_names[nSamples_eff +1 ] = { "total",      
						"M-2_V-0.00244948974278_mu",
						"M-2_V-0.00282842712475_mu",
						"M-2_V-0.00316227766017_mu",
						"M-2_V-0.004472135955_mu",
						"M-2_V-0.00547722557505_mu",
						"M-2_V-0.00707106781187_mu",
						"M-2_V-0.00836660026534_mu",
						"M-2_V-0.0141421356237_mu",
						"M-2_V-0.0173205080757_mu",
						"M-2_V-0.01_mu",
						"M-2_V-0.022360679775_mu",

						"M-3_V-0.00244948974278_mu",
						"M-3_V-0.00282842712475_mu",
						"M-3_V-0.00316227766017_mu",
						"M-3_V-0.004472135955_mu",
						"M-3_V-0.00547722557505_mu",
						"M-3_V-0.00707106781187_mu",
						"M-3_V-0.00836660026534_mu",
						"M-3_V-0.0141421356237_mu",
						"M-3_V-0.0173205080757_mu",
						"M-3_V-0.01_mu",

						"M-4_V-0.00244948974278_mu",
						"M-4_V-0.00282842712475_mu",
						"M-4_V-0.00316227766017_mu",
						"M-4_V-0.004472135955_mu",
						"M-4_V-0.00547722557505_mu",
						"M-4_V-0.00707106781187_mu",
						"M-4_V-0.00836660026534_mu",
						"M-4_V-0.0141421356237_mu",
						"M-4_V-0.01_mu",

						"M-5_V-0.00244948974278_mu",
						"M-5_V-0.00282842712475_mu",
						"M-5_V-0.00316227766017_mu",
						"M-5_V-0.004472135955_mu",
						"M-5_V-0.00547722557505_mu",
						"M-5_V-0.00707106781187_mu",
						"M-5_V-0.00836660026534_mu",
						"M-5_V-0.01_mu",

						"M-8_V-0.00244948974278_mu",
						"M-8_V-0.00282842712475_mu",
						"M-8_V-0.00316227766017_mu",
						"M-8_V-0.004472135955_mu",
						"M-8_V-0.00547722557505_mu",

						"M-10_V-0.00244948974278_mu",
						"M-10_V-0.00282842712475_mu",
						"M-10_V-0.00316227766017_mu",
						"M-10_V-0.004472135955_mu",

						"M-15_V-0.00244948974278_mu",
						"M-15_V-0.00282842712475_mu",
						"M-15_V-0.00316227766017_mu",
						"M-15_V-0.004472135955_mu",


						"M-20_V-1.93390796058e-05_mu",

						"M-1_V-0.0949736805647_mu",
						"M-1_V-0.212367605816_mu",
						"M-2_V-0.0110905365064_mu",
						"M-2_V-0.0248394846967_mu",
						"M-3_V-0.00707813534767_mu",
						"M-4_V-0.00290516780927_mu",
						"M-5_V-0.00145602197786_mu",
						"M-6_V-0.00202484567313_mu",
						"M-8_V-0.00151327459504_mu",
						"M-10_V-0.000756967634711_mu",
  
						"DY",  
						"t#bar{t}",
						"WJets",
						"multiboson", 
						"X+#gamma",    
						"TT/T + X",		
						"no prompt"};

  // 9 0.05173
  //61526.7, wjet
  const double xSections[nSamples]= {0,        

				     fl*2.570e-02,
				     fl*3.436e-02 ,
				     fl*4.295e-02 ,
				     fl*8.619e-02 ,
				     fl*1.290e-01 ,
				     fl*2.151e-01,
				     fl*3.022e-01,
				     fl*8.631e-01 ,
				     fl*1.296,
				     fl*4.295e-01,
				     fl*fl*2.158,

				     fl*2.423e-02,
				     fl*3.238e-02,
				     fl*4.024e-02,
				     fl*8.099e-02,
				     fl*1.211e-01,
				     fl*2.028e-01,
				     fl*2.812e-01,
				     fl*fl*8.070e-01,
				     fl*fl*1.174,
				     fl*4.029e-01,


				     fl*2.40e-02,
				     fl*3.16e-02,
				     fl*3.99e-02,
				     fl*7.90e-02,
				     fl*1.19e-01,
				     fl*2.02e-01,
				     fl*2.80e-01,
				     fl*fl*7.896e-01,
				     fl*3.95e-01,

				     

				     fl*2.401e-02,
				     fl*3.177e-02,
				     fl*3.929e-02,
				     fl*7.967e-02,
				     fl*1.186e-01,
				     fl*2.004e-01,
				     fl*2.798e-01,
				     fl*4.007e-01,


				     fl*2.444e-02,
				     fl*3.286e-02,
				     fl*4.117e-02,
				     fl*8.243e-02,
				     fl*1.240e-01,

				     fl*2.47e-02,
				     fl*3.29e-02,
				     fl*4.15e-02,
				     fl*8.29e-02,

				     fl*2.45e-02,
				     fl*3.26e-02,
				     fl*4.09e-02,
				     fl*8.13e-02,

				     fl*1.442e-06,


				     fl*fl*3.94e+01,
				     fl*fl*fl*1.97e+02,
				     fl*5.37e-01,
				     fl*fl*2.69,
				     fl*2.06e-01,
				     fl*3.41e-02,
				     fl*8.60e-03,
				     fl*1.66e-02,
				     fl*9.42e-03,
				     fl*2.37e-03,


				     18610, 1921.8*3,
				     87.315, 182.17540224, 182.17540224,  3.68064 , 80.95 ,  136.02 , 35.85, 35.85 , 
				     61526.7,

				     0.000652,0.001034,

				     0.2086, 0.1651,  0.01398,0.05565,0.04123 , 0.2147 , 1.711, 0.006699,0.006699,0.006699,0.003339,0.003339,0.003339,12.178, 0.1729, 1.256, 58.59,0.0002339,


				     489,117.864,
				     
				     2.967, 0.01731, 3.697, 0.2043,0.40620, 0.2728,0.0283,0.0942, 0.009103, 0.007834,0.002938, 0.00156, 0.0110
				     
				    

				    
  };
  */
    
  double luminosity = 35.867;
    
  TFile *hfile[nSamples];
  TTree *inputTree[nSamples];
 
  ULong64_t       _runNb;
  ULong64_t       _lumiBlock;
  ULong64_t       _eventNb;
  UChar_t         _nVertex;
  Double_t        _weight = 0.;
  Double_t        _lheHTIncoming;
  Double_t        _ctauHN;
  UInt_t         _nLheWeights;
  _nLheWeights = 110;
  Double_t        _lheWeight[110];   //[_nLheWeights]
  Float_t         _nTrueInt = 0;
  UInt_t          _ttgEventType;
  UInt_t          _zgEventType;
  Double_t        _gen_met;
  Double_t        _gen_metPhi;
  UInt_t          _gen_nPh;
  _gen_nPh = 7;
  Double_t        _gen_phPt[7];   //[_gen_nPh]
  Double_t        _gen_phEta[7];   //[_gen_nPh]
  Double_t        _gen_phPhi[7];   //[_gen_nPh]
  Double_t        _gen_phE[7];   //[_gen_nPh]
  Int_t           _gen_phMomPdg[7];   //[_gen_nPh]
  Bool_t          _gen_phIsPrompt[7];   //[_gen_nPh]
  Double_t        _gen_phMinDeltaR[7];   //[_gen_nPh]
  Bool_t          _gen_phPassParentage[7];   //[_gen_nPh]
  UInt_t          _gen_nL;
  _gen_nL = 20;
  Double_t        _gen_lPt[20];   //[_gen_nL]
  Double_t        _gen_lEta[20];   //[_gen_nL]
  Double_t        _gen_lPhi[20];   //[_gen_nL]
  Double_t        _gen_lE[20];   //[_gen_nL]
  UInt_t          _gen_lFlavor[20];   //[_gen_nL]
  Int_t           _gen_lCharge[20];   //[_gen_nL]
  Int_t           _gen_lMomPdg[20];   //[_gen_nL]
  Bool_t          _gen_lIsPrompt[20];   //[_gen_nL]
 
  UInt_t         _nL;
  _nL = 20;
  UInt_t         _nMu;
  _nMu = 20;
  UInt_t         _nEle;
  _nEle = 20;
  UInt_t         _nLight;
  _nLight = 20;
  UInt_t         _nTau;
  UInt_t         _nVFit;
  _nVFit  = 25;
  UInt_t         _nGoodLeading;
  UInt_t        _lIndex[20];   //[_nL]
  double        _vertices[25][12];
  double        _lDisplaced[25][24];   //[_nVFit]

  Double_t        _lPt[20];   //[_nL]
  Double_t        _lEta[20];   //[_nL]
  Double_t        _lEtaSC[20];   //[_nLight]
  Double_t        _lPhi[20];   //[_nL]
  Double_t        _lE[20];   //[_nL]
  UInt_t          _lFlavor[20];   //[_nL]
  Int_t           _lCharge[20];   //[_nL]
  double        _dxy[20];   //[_nL]
  double        _dz[20];   //[_nL]
  double        _3dIP[20];   //[_nL]
  double        _3dIPSig[20];   //[_nL]
  double        _2dIP[20];   //[_nL]
  double        _2dIPSig[20];   //[_nL]
  double   _lMatchPt[20];
  double   _lMatchEta[20];
  double   _lMatchPhi[20];

  double   _lMatchVertexX[20];
  double   _lMatchVertexY[20];
  double   _lMatchVertexZ[20];
  double   _pvX = 0.;
  double   _pvY = 0.;
  double   _pvZ = 0.;
  double   _pvXErr = 0.;
  double   _pvYErr = 0.;
  double   _pvZErr = 0.;
  double _closestJetCsvV2[20];
  double _closestJetDeepCsv_b[20];
  double           _lMuTime[20];
  double          _lMuTimeErr[20];
  double          _lMuRPCTime[20];
  double          _lMuRPCTimeErr[20];
  int             _lMuTimenDof[20];
  int             _lMuRPCTimenDof[20];
  Float_t         _lElectronMva[20];   //[_nLight]
  Bool_t          _lElectronPassEmu[20];   //[_nLight]
  Bool_t          _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto[20];   //[_nL]
  Bool_t          _lPOGVeto[20];   //[_nL]
  Bool_t          _lPOGLoose[20];   //[_nL]
  Bool_t          _lPOGMedium[20];   //[_nL]
  Bool_t          _lPOGTight[20];   //[_nL]
  Bool_t          _lpassConversionVeto[20];   //[_nL]
  UInt_t          _eleNumberInnerHitsMissing[20];   //[_nL]
  Bool_t          _lGlobalMuon[20];   //[_nL]
  Bool_t          _lTrackerMuon[20];   //[_nL]
  double        _lInnerTrackValidFraction[20];   //[_nL]
  double        _lGlobalTrackNormalizeChi2[20];   //[_nL]
  double        _lCQChi2Position[20];   //[_nL]
  double        _lCQTrackKink[20];   //[_nL]
  UInt_t          _lNumberOfMatchedStation[20];   //[_nL]
  UInt_t          _lNumberOfValidPixelHits[20];   //[_nL]
  UInt_t          _muNumberInnerHits[20];   //[_nL]
  UInt_t          _lTrackerLayersWithMeasurement[20];   //[_nL]
  Bool_t          _lEleIsEB[20];   //[_nL]
  Bool_t          _lEleIsEE[20];   //[_nL]
  double        _lEleSuperClusterOverP[20];   //[_nL]
  double        _lEleEcalEnergy[20];   //[_nL]
  double        _lElefull5x5SigmaIetaIeta[20];   //[_nL]
  double        _lEleDEtaInSeed[20];   //[_nL]
  double        _lEleDeltaPhiSuperClusterTrackAtVtx[20];   //[_nL
  double        _lElehadronicOverEm[20];   //[_nL]
  double        _lEleInvMinusPInv[20];   //[_nL]
  double        _relIso[20];   //[_nLight]
  double        _puCorr[20];   //[_nL]
  double        _absIso03[20];   //[_nL]
  double        _absIso04[20];   //[_nL]
  double        _sumNeutralHadronEt04[20];   //[_nL]
  double        _sumChargedHadronPt04[20];   //[_nL]
  double        _sumPhotonEt04[20];   //[_nL]
  double        _sumNeutralHadronEt03[20];   //[_nL]
  double        _sumChargedHadronPt03[20];   //[_nL]
  double        _sumPhotonEt03[20];   //[_nL]
  double        _trackIso[20];   //[_nL]
  double        _ecalIso[20];   //[_nL]
  double        _hcalIso[20];   //[_nL]
  //double        _deltaBIso[20];   //[_nL]
  double        _ecalPFClusterIso[20];   //[_nL]
  double        _hcalPFClusterIso[20];   //[_nL]
  double        _ptRel[20];   //[_nLight]
  double        _ptRatio[20];   //[_nLight]
  double        _closestJetCsv[20];   //[_nLight]
  UInt_t          _selectedTrackMult[20];   //[_nLight]
  double        _muonSegComp[20];   //[_nMu]
  Bool_t          _tauMuonVeto[20];   //[_nL]
  Bool_t          _tauEleVeto[20];   //[_nL]
  Bool_t          _decayModeFindingNew[20];   //[_nL]
  Bool_t          _tauVLooseMvaNew[20];   //[_nL]
  Bool_t          _tauLooseMvaNew[20];   //[_nL]
  Bool_t          _tauMediumMvaNew[20];   //[_nL]
  Bool_t          _tauTightMvaNew[20];   //[_nL]
  Bool_t          _tauVTightMvaNew[20];   //[_nL]
  Bool_t          _tauVTightMvaOld[20];   //[_nL]
  Bool_t          _lIsPrompt[20];   //[_nL]
  Int_t           _lMatchPdgId[20];   //[_nL]
  UInt_t          _lProvenance[20];   //[_nL]
  UInt_t          _lProvenanceCompressed[20];   //[_nL]
  UInt_t         _lHasTrigger[20];
  UInt_t         _nPh;
  _nPh = 20;
  Double_t        _phPt[20];   //[_nPh]
  Double_t        _phEta[20];   //[_nPh]
  Double_t        _phEtaSC[20];   //[_nPh]
  Double_t        _phPhi[20];   //[_nPh]
  Double_t        _phE[20];   //[_nPh]
  Bool_t          _phCutBasedLoose[20];   //[_nPh]
  Bool_t          _phCutBasedMedium[20];   //[_nPh]
  Bool_t          _phCutBasedTight[20];   //[_nPh]
  Double_t        _phMva[20];   //[_nPh]
  Double_t        _phRandomConeChargedIsolation[20];   //[_nPh]
  Double_t        _phChargedIsolation[20];   //[_nPh]
  Double_t        _phNeutralHadronIsolation[20];   //[_nPh]
  Double_t        _phPhotonIsolation[20];   //[_nPh]
  Double_t        _phSigmaIetaIeta[20];   //[_nPh]
  //Double_t        _phSigmaIetaIphi[20];   //[_nPh]
  Double_t        _phHadronicOverEm[20];   //[_nPh]
  Bool_t          _phPassElectronVeto[20];   //[_nPh]
  Bool_t          _phHasPixelSeed[20];   //[_nPh]
  Bool_t          _phIsPrompt[20];   //[_nPh]
  Int_t           _phMatchMCPhotonAN15165[20];   //[_nPh]
  Int_t           _phMatchMCLeptonAN15165[20];   //[_nPh]
  Int_t           _phMatchPdgId[20];   //[_nPh]
  UInt_t         _nJets;
  _nJets = 20;
  Double_t        _jetPt[20];   //[_nJets]
  Double_t        _jetSmearedPt[20];   //[_nJets]
  Double_t        _jetSmearedPt_JECUp[20];   //[_nJets]
  Double_t        _jetSmearedPt_JECDown[20];   //[_nJets]
  Double_t        _jetSmearedPt_JERUp[20];   //[_nJets]
  Double_t        _jetSmearedPt_JERDown[20];   //[_nJets]
  Double_t        _jetEta[20];   //[_nJets]
  Double_t        _jetPhi[20];   //[_nJets]
  Double_t        _jetE[20];   //[_nJets]
  Double_t        _jetCsvV2[20];   //[_nJets]
  Double_t        _jetDeepCsv_udsg[20];   //[_nJets]
  Double_t        _jetDeepCsv_b[20];   //[_nJets]
  Double_t        _jetDeepCsv_c[20];   //[_nJets]
  Double_t        _jetDeepCsv_bb[20];   //[_nJets]
  UInt_t          _jetHadronFlavor[20];   //[_nJets]
  //UInt_t          _jetId[20];   //[_nJets]
  bool            _jetIsTight[20];
  Double_t        _met = 0.;
  Double_t        _metJECDown;
  Double_t        _metJECUp;
  Double_t        _metUnclDown;
  Double_t        _metUnclUp;
  Double_t        _metPhi = 0.;
  Double_t        _metPhiJECDown;
  Double_t        _metPhiJECUp;
  Double_t        _metPhiUnclDown;
  Double_t        _metPhiUnclUp;
    
  // List of branches
  TBranch        *b__runNb;   //!
  TBranch        *b__lumiBlock;   //!
  TBranch        *b__eventNb;   //!
  TBranch        *b__nVertex;   //!
  TBranch        *b__weight;   //!
  TBranch        *b__lheHTIncoming;   //!
  TBranch        *b__ctauHN;   //!
  TBranch        *b__nLheWeights;   //!
  TBranch        *b__lheWeight;   //!
  TBranch        *b__nTrueInt;   //!
  TBranch        *b__ttgEventType;   //!
  TBranch        *b__zgEventType;   //!
  TBranch        *b__gen_met;   //!
  TBranch        *b__gen_metPhi;   //!
  TBranch        *b__gen_nPh;   //!
  TBranch        *b__gen_phPt;   //!
  TBranch        *b__gen_phEta;   //!
  TBranch        *b__gen_phPhi;   //!
  TBranch        *b__gen_phE;   //!
  TBranch        *b__gen_phMomPdg;   //!
  TBranch        *b__gen_phIsPrompt;   //!
  TBranch        *b__gen_phMinDeltaR;   //!
  TBranch        *b__gen_phPassParentage;   //!
  TBranch        *b__gen_nL;   //!
  TBranch        *b__gen_lPt;   //!
  TBranch        *b__gen_lEta;   //!
  TBranch        *b__gen_lPhi;   //!
  TBranch        *b__gen_lE;   //!
  TBranch        *b__gen_lFlavor;   //!
  TBranch        *b__gen_lCharge;   //!
  TBranch        *b__gen_lMomPdg;   //!
  TBranch        *b__gen_lIsPrompt;   //!
 
  TBranch        *b__nL;   //!
  TBranch        *b__nMu;   //!
  TBranch        *b__nEle;   //!
  TBranch        *b__nLight;   //!
  TBranch        *b__nTau;   //!
  TBranch        *b__nVFit;   //!
  TBranch        *b__nGoodLeading;   //!
  //TBranch        *b__lIndex;   //!
  TBranch        *b__vertices;   //!
  TBranch        *b__lDisplaced;   //!
   TBranch *b__closestJetCsvV2;
  TBranch *b__closestJetDeepCsv_b;

  TBranch        *b__lPt;   //!
  TBranch        *b__lEta;   //!
  TBranch        *b__lEtaSC;   //!
  TBranch        *b__lPhi;   //!
  TBranch        *b__lE;   //!
  TBranch        *b__lFlavor;   //!
  TBranch        *b__lCharge;   //!
  TBranch        *b__dxy;   //!
  TBranch        *b__dz;   //!
  TBranch        *b__3dIP;   //!
  TBranch        *b__3dIPSig;   //!
  TBranch        *b__2dIP;   //!
  TBranch        *b__2dIPSig;   //!
  TBranch        *b__lElectronMva;   //!
  TBranch        *b__lElectronPassEmu;   //!
  TBranch        *b__lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto;   //!
  TBranch        *b__lPOGVeto;   //!
  TBranch        *b__lPOGLoose;   //!
  TBranch        *b__lPOGMedium;   //!
  TBranch        *b__lPOGTight;   //!
  TBranch        *b__lpassConversionVeto;   //!
  TBranch        *b__eleNumberInnerHitsMissing;   //!
  TBranch        *b__lGlobalMuon;   //!
  TBranch        *b__lTrackerMuon;   //!
  TBranch        *b__lInnerTrackValidFraction;   //!
  TBranch        *b__lGlobalTrackNormalizeChi2;   //!
  TBranch        *b__lCQChi2Position;   //!
  TBranch        *b__lCQTrackKink;   //!
  TBranch        *b__lNumberOfMatchedStation;   //!
  TBranch        *b__lNumberOfValidPixelHits;   //!
  TBranch        *b__muNumberInnerHits;   //!
  TBranch        *b__lTrackerLayersWithMeasurement;   //!
  TBranch           *b__lMuTime;
  TBranch          *b__lMuTimeErr;
  TBranch          *b__lMuRPCTime;
  TBranch          *b__lMuRPCTimeErr;
  TBranch             *b__lMuTimenDof;
  TBranch             *b__lMuRPCTimenDof;
  TBranch        *b__lEleIsEB;   //!
  TBranch        *b__lEleIsEE;   //!
  TBranch        *b__lEleSuperClusterOverP;   //!
  TBranch        *b__lEleEcalEnergy;   //!
  TBranch        *b__lElefull5x5SigmaIetaIeta;   //!
  TBranch        *b__lEleDEtaInSeed;   //!
  TBranch        *b__lEleDeltaPhiSuperClusterTrackAtVtx;   //!
  TBranch        *b__lElehadronicOverEm;   //!
  TBranch        *b__lEleInvMinusPInv;   //!
  TBranch        *b__relIso;   //!
  TBranch        *b__lHasTrigger;
  TBranch        *b__puCorr;   //!
  TBranch        *b__absIso03;   //!
  TBranch        *b__absIso04;   //!
  TBranch        *b__sumNeutralHadronEt04;   //!
  TBranch        *b__sumChargedHadronPt04;   //!
  TBranch        *b__sumPhotonEt04;   //!
  TBranch        *b__sumNeutralHadronEt03;   //!
  TBranch        *b__sumChargedHadronPt03;   //!
  TBranch        *b__sumPhotonEt03;   //!
  TBranch        *b__trackIso;   //!
  TBranch        *b__ecalIso;   //!
  TBranch        *b__hcalIso;   //!
  //TBranch        *b__deltaBIso;   //!
  TBranch        *b__ecalPFClusterIso;   //!
  TBranch        *b__hcalPFClusterIso;   //!
  TBranch        *b__ptRel;   //!
  TBranch        *b__ptRatio;   //!
  TBranch        *b__closestJetCsv;   //!
  TBranch        *b__selectedTrackMult;   //!
  TBranch        *b__muonSegComp;   //!
  TBranch        *b__tauMuonVeto;   //!
  TBranch        *b__tauEleVeto;   //!
  TBranch        *b__decayModeFindingNew;   //!
  TBranch        *b__tauVLooseMvaNew;   //!
  TBranch        *b__tauLooseMvaNew;   //!
  TBranch        *b__tauMediumMvaNew;   //!
  TBranch        *b__tauTightMvaNew;   //!
  TBranch        *b__tauVTightMvaNew;   //!
  TBranch        *b__tauVTightMvaOld;   //!
  TBranch        *b__lIsPrompt;   //!
  TBranch        *b__lMatchPdgId;   //!
  TBranch        *b__lProvenance;
  TBranch        *b__lProvenanceCompressed;


  TBranch   *b__lMatchPt;
  TBranch   *b__lMatchEta;
  TBranch   *b__lMatchPhi;

  TBranch   *b__lMatchVertexX;
  TBranch   *b__lMatchVertexY;
  TBranch   *b__lMatchVertexZ;
  TBranch   *b__pvX;
  TBranch   *b__pvY;
  TBranch   *b__pvZ;
  TBranch   *b__pvXErr;
  TBranch   *b__pvYErr;
  TBranch   *b__pvZErr;
  TBranch        *b__nPh;   //!
  TBranch        *b__phPt;   //!
  TBranch        *b__phEta;   //!
  TBranch        *b__phEtaSC;   //!
  TBranch        *b__phPhi;   //!
  TBranch        *b__phE;   //!
  TBranch        *b__phCutBasedLoose;   //!
  TBranch        *b__phCutBasedMedium;   //!
  TBranch        *b__phCutBasedTight;   //!
  TBranch        *b__phMva;   //!
  TBranch        *b__phRandomConeChargedIsolation;   //!
  TBranch        *b__phChargedIsolation;   //!
  TBranch        *b__phNeutralHadronIsolation;   //!
  TBranch        *b__phPhotonIsolation;   //!
  TBranch        *b__phSigmaIetaIeta;   //!
  //TBranch        *b__phSigmaIetaIphi;   //!
  TBranch        *b__phHadronicOverEm;   //!
  TBranch        *b__phPassElectronVeto;   //!
  TBranch        *b__phHasPixelSeed;   //!
  TBranch        *b__phIsPrompt;   //!
  TBranch        *b__phMatchMCPhotonAN15165;   //!
  TBranch        *b__phMatchMCLeptonAN15165;   //!
  TBranch        *b__phMatchPdgId;   //!
  TBranch        *b__nJets;   //!
  TBranch        *b__jetPt;   //!
  TBranch        *b__jetSmearedPt;
  TBranch        *b__jetSmearedPt_JECUp;   //!
  TBranch        *b__jetSmearedPt_JECDown;   //!
  TBranch        *b__jetSmearedPt_JERUp;   //!
  TBranch        *b__jetSmearedPt_JERDown;   //!
  TBranch        *b__jetEta;   //!
  TBranch        *b__jetPhi;   //!
  TBranch        *b__jetE;   //!
  TBranch        *b__jetCsvV2;   //!
  TBranch       *b__jetIsTight;   //!
  TBranch        *b__jetDeepCsv_udsg;   //!
  TBranch        *b__jetDeepCsv_b;   //!
  TBranch        *b__jetDeepCsv_c;   //!
  TBranch        *b__jetDeepCsv_bb;   //!
  TBranch        *b__jetHadronFlavor;   //!
  //TBranch        *b__jetId;   //!
  TBranch        *b__met;   //!
  TBranch        *b__metJECDown;   //!
  TBranch        *b__metJECUp;   //!
  TBranch        *b__metUnclDown;   //!
  TBranch        *b__metUnclUp;   //!
  TBranch        *b__metPhi;   //!
  TBranch        *b__metPhiJECDown;   //!
  TBranch        *b__metPhiJECUp;   //!
  TBranch        *b__metPhiUnclDown;   //!
  TBranch        *b__metPhiUnclUp;   //!
    
  double hcounter[nSamples];

  // Dummy
  for(unsigned i=0; i<20; ++i){
    _lIndex[i] = i+1;
  }

  
  for(int sam = 0; sam < nSamples; ++sam){
        
    cout<<"================================= 1"<<endl;
    cout<<"--------> "<< "name " << names[sam] << endl;
    cout<<"--------> "<< "fileList[sam] " << fileList[sam] << endl;
    // Signal
    //hfile[sam] = new TFile("/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/ANALYSIS/20190318_theoryUncertainties/samples.noSync/"+fileList[sam],"read");
    // Background
    hfile[sam] = new TFile("/Users/trocino/Documents/Work/Analysis/HeavyNeutrino/ANALYSIS/20190318_MartinasCode/samples.noSync/2016/"+fileList[sam],"read");
    hfile[sam]->cd("blackJackAndHookers");
    inputTree[sam] = static_cast<TTree*>(hfile[sam]->Get("blackJackAndHookers/blackJackAndHookersTree"));
    inputTree[sam]->SetBranchAddress("_runNb", &_runNb, &b__runNb);
    inputTree[sam]->SetBranchAddress("_lumiBlock", &_lumiBlock, &b__lumiBlock);
    inputTree[sam]->SetBranchAddress("_eventNb", &_eventNb, &b__eventNb);
    inputTree[sam]->SetBranchAddress("_nVertex", &_nVertex, &b__nVertex);
    inputTree[sam]->SetBranchAddress("_met", &_met, &b__met);
    inputTree[sam]->SetBranchAddress("_metJECDown", &_metJECDown, &b__metJECDown);
    inputTree[sam]->SetBranchAddress("_metJECUp", &_metJECUp, &b__metJECUp);
    inputTree[sam]->SetBranchAddress("_metUnclDown", &_metUnclDown, &b__metUnclDown);
    inputTree[sam]->SetBranchAddress("_metUnclUp", &_metUnclUp, &b__metUnclUp);
    inputTree[sam]->SetBranchAddress("_metPhi", &_metPhi, &b__metPhi);
    inputTree[sam]->SetBranchAddress("_metPhiJECDown", &_metPhiJECDown, &b__metPhiJECDown);
    inputTree[sam]->SetBranchAddress("_metPhiJECUp", &_metPhiJECUp, &b__metPhiJECUp);
    inputTree[sam]->SetBranchAddress("_metPhiUnclDown", &_metPhiUnclDown, &b__metPhiUnclDown);
    inputTree[sam]->SetBranchAddress("_metPhiUnclUp", &_metPhiUnclUp, &b__metPhiUnclUp);
    inputTree[sam]->SetBranchAddress("_nTrueInt", &_nTrueInt, &b__nTrueInt);
    inputTree[sam]->SetBranchAddress("_weight", &_weight, &b__weight);
    inputTree[sam]->SetBranchAddress("_lheHTIncoming", &_lheHTIncoming, &b__lheHTIncoming);
    inputTree[sam]->SetBranchAddress("_ctauHN", &_ctauHN, &b__ctauHN);
    inputTree[sam]->SetBranchAddress("_nLheWeights", &_nLheWeights, &b__nLheWeights);
    inputTree[sam]->SetBranchAddress("_lheWeight", _lheWeight, &b__lheWeight);
    inputTree[sam]->SetBranchAddress("_ttgEventType", &_ttgEventType, &b__ttgEventType);
    inputTree[sam]->SetBranchAddress("_zgEventType", &_zgEventType, &b__zgEventType);
    inputTree[sam]->SetBranchAddress("_gen_met", &_gen_met, &b__gen_met);
    inputTree[sam]->SetBranchAddress("_gen_metPhi", &_gen_metPhi, &b__gen_metPhi);
    inputTree[sam]->SetBranchAddress("_gen_nPh", &_gen_nPh, &b__gen_nPh);
    inputTree[sam]->SetBranchAddress("_gen_phPt", _gen_phPt, &b__gen_phPt);
    inputTree[sam]->SetBranchAddress("_gen_phEta", _gen_phEta, &b__gen_phEta);
    inputTree[sam]->SetBranchAddress("_gen_phPhi", _gen_phPhi, &b__gen_phPhi);
    inputTree[sam]->SetBranchAddress("_gen_phE", _gen_phE, &b__gen_phE);
    inputTree[sam]->SetBranchAddress("_gen_phMomPdg", _gen_phMomPdg, &b__gen_phMomPdg);
    inputTree[sam]->SetBranchAddress("_gen_phIsPrompt", _gen_phIsPrompt, &b__gen_phIsPrompt);
    inputTree[sam]->SetBranchAddress("_gen_nL", &_gen_nL, &b__gen_nL);
    inputTree[sam]->SetBranchAddress("_gen_lPt", _gen_lPt, &b__gen_lPt);
    inputTree[sam]->SetBranchAddress("_gen_lEta", _gen_lEta, &b__gen_lEta);
    inputTree[sam]->SetBranchAddress("_gen_lPhi", _gen_lPhi, &b__gen_lPhi);
    inputTree[sam]->SetBranchAddress("_gen_lE", _gen_lE, &b__gen_lE);
    inputTree[sam]->SetBranchAddress("_gen_lFlavor", _gen_lFlavor, &b__gen_lFlavor);
    inputTree[sam]->SetBranchAddress("_gen_lCharge", _gen_lCharge, &b__gen_lCharge);
    inputTree[sam]->SetBranchAddress("_gen_lMomPdg", _gen_lMomPdg, &b__gen_lMomPdg);
    inputTree[sam]->SetBranchAddress("_gen_lIsPrompt", _gen_lIsPrompt, &b__gen_lIsPrompt);
    
    inputTree[sam]->SetBranchAddress("_nL", &_nL, &b__nL);
    inputTree[sam]->SetBranchAddress("_nMu", &_nMu, &b__nMu);
    inputTree[sam]->SetBranchAddress("_nEle", &_nEle, &b__nEle);
    inputTree[sam]->SetBranchAddress("_nLight", &_nLight, &b__nLight);
    inputTree[sam]->SetBranchAddress("_nTau", &_nTau, &b__nTau);
    inputTree[sam]->SetBranchAddress("_nVFit", &_nVFit, &b__nVFit);
    inputTree[sam]->SetBranchAddress("_nGoodLeading", &_nGoodLeading, &b__nGoodLeading);
    //inputTree[sam]->SetBranchAddress("_lIndex", _lIndex, &b__lIndex);
    inputTree[sam]->SetBranchAddress("_vertices", _vertices, &b__vertices);
    inputTree[sam]->SetBranchAddress("_lDisplaced", _lDisplaced, &b__lDisplaced);
    inputTree[sam]->SetBranchAddress("_lPt", _lPt, &b__lPt);
    inputTree[sam]->SetBranchAddress("_lEta", _lEta, &b__lEta);
    inputTree[sam]->SetBranchAddress("_lEtaSC", _lEtaSC, &b__lEtaSC);
    inputTree[sam]->SetBranchAddress("_lPhi", _lPhi, &b__lPhi);
    inputTree[sam]->SetBranchAddress("_lE", _lE, &b__lE);
    inputTree[sam]->SetBranchAddress("_lFlavor", _lFlavor, &b__lFlavor);
    inputTree[sam]->SetBranchAddress("_lCharge", _lCharge, &b__lCharge);
    inputTree[sam]->SetBranchAddress("_dxy", _dxy, &b__dxy);
    inputTree[sam]->SetBranchAddress("_dz", _dz, &b__dz);
    inputTree[sam]->SetBranchAddress("_3dIP", _3dIP, &b__3dIP);
    inputTree[sam]->SetBranchAddress("_3dIPSig", _3dIPSig, &b__3dIPSig);
    inputTree[sam]->SetBranchAddress("_2dIP", _2dIP, &b__2dIP);
    inputTree[sam]->SetBranchAddress("_2dIPSig", _2dIPSig, &b__2dIPSig);
    //inputTree[sam]->SetBranchAddress("_lElectronMva", _lElectronMva, &b__lElectronMva);
    inputTree[sam]->SetBranchAddress("_lElectronMvaFall17Iso", _lElectronMva, &b__lElectronMva);
    inputTree[sam]->SetBranchAddress("_lElectronPassEmu", _lElectronPassEmu, &b__lElectronPassEmu);
    inputTree[sam]->SetBranchAddress("_lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto", _lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto, &b__lLooseCBwoIsolationwoMissingInnerhitswoConversionVeto);
    inputTree[sam]->SetBranchAddress("_lPOGVeto", _lPOGVeto, &b__lPOGVeto);
    inputTree[sam]->SetBranchAddress("_lPOGLoose", _lPOGLoose, &b__lPOGLoose);
    inputTree[sam]->SetBranchAddress("_lPOGMedium", _lPOGMedium, &b__lPOGMedium);
    inputTree[sam]->SetBranchAddress("_lPOGTight", _lPOGTight, &b__lPOGTight);
    inputTree[sam]->SetBranchAddress("_lElectronPassConvVeto", _lpassConversionVeto, &b__lpassConversionVeto);
    inputTree[sam]->SetBranchAddress("_lElectronMissingHits", _eleNumberInnerHitsMissing, &b__eleNumberInnerHitsMissing);
    inputTree[sam]->SetBranchAddress("_lGlobalMuon", _lGlobalMuon, &b__lGlobalMuon);
    inputTree[sam]->SetBranchAddress("_lTrackerMuon", _lTrackerMuon, &b__lTrackerMuon);
    inputTree[sam]->SetBranchAddress("_lInnerTrackValidFraction", _lInnerTrackValidFraction, &b__lInnerTrackValidFraction);
    inputTree[sam]->SetBranchAddress("_lGlobalTrackNormalizeChi2", _lGlobalTrackNormalizeChi2, &b__lGlobalTrackNormalizeChi2);
    inputTree[sam]->SetBranchAddress("_lCQChi2Position", _lCQChi2Position, &b__lCQChi2Position);
    inputTree[sam]->SetBranchAddress("_lCQTrackKink", _lCQTrackKink, &b__lCQTrackKink);
    inputTree[sam]->SetBranchAddress("_lMuonSegComp", _muonSegComp, &b__muonSegComp);
    inputTree[sam]->SetBranchAddress("_lNumberOfMatchedStation", _lNumberOfMatchedStation, &b__lNumberOfMatchedStation);
    inputTree[sam]->SetBranchAddress("_lNumberOfValidPixelHits", _lNumberOfValidPixelHits, &b__lNumberOfValidPixelHits);
    inputTree[sam]->SetBranchAddress("_muNumberInnerHits", _muNumberInnerHits, &b__muNumberInnerHits);
    inputTree[sam]->SetBranchAddress("_lTrackerLayersWithMeasurement", _lTrackerLayersWithMeasurement, &b__lTrackerLayersWithMeasurement);
    inputTree[sam]->SetBranchAddress("_lMuTime", _lMuTime, &b__lMuTime);
    inputTree[sam]->SetBranchAddress("_lMuTimeErr", _lMuTimeErr, &b__lMuTimeErr);
    inputTree[sam]->SetBranchAddress("_lMuRPCTime", _lMuRPCTime, &b__lMuRPCTime);
    inputTree[sam]->SetBranchAddress("_lMuRPCTimeErr", _lMuRPCTimeErr, &b__lMuRPCTimeErr);
    inputTree[sam]->SetBranchAddress("_lMuTimenDof", _lMuTimenDof, &b__lMuTimenDof);
    inputTree[sam]->SetBranchAddress("_lMuRPCTimenDof", _lMuRPCTimenDof, &b__lMuRPCTimenDof);
    inputTree[sam]->SetBranchAddress("_lEleIsEB", _lEleIsEB, &b__lEleIsEB);
    inputTree[sam]->SetBranchAddress("_lEleIsEE", _lEleIsEE, &b__lEleIsEE);
    inputTree[sam]->SetBranchAddress("_lEleSuperClusterOverP", _lEleSuperClusterOverP, &b__lEleSuperClusterOverP);
    inputTree[sam]->SetBranchAddress("_lEleEcalEnergy", _lEleEcalEnergy, &b__lEleEcalEnergy);
    inputTree[sam]->SetBranchAddress("_lElefull5x5SigmaIetaIeta", _lElefull5x5SigmaIetaIeta, &b__lElefull5x5SigmaIetaIeta);
    inputTree[sam]->SetBranchAddress("_lEleDEtaInSeed", _lEleDEtaInSeed, &b__lEleDEtaInSeed);
    inputTree[sam]->SetBranchAddress("_lEleDeltaPhiSuperClusterTrackAtVtx", _lEleDeltaPhiSuperClusterTrackAtVtx, &b__lEleDeltaPhiSuperClusterTrackAtVtx);
    inputTree[sam]->SetBranchAddress("_lElehadronicOverEm", _lElehadronicOverEm, &b__lElehadronicOverEm);
    inputTree[sam]->SetBranchAddress("_lEleInvMinusPInv", _lEleInvMinusPInv, &b__lEleInvMinusPInv);
    //    inputTree[sam]->SetBranchAddress("_eleNumberInnerHitsMissing", _eleNumberInnerHitsMissing, &b__eleNumberInnerHitsMissing);
    inputTree[sam]->SetBranchAddress("_relIso", _relIso, &b__relIso);
    inputTree[sam]->SetBranchAddress("_lHasTrigger", _lHasTrigger, &b__lHasTrigger);

    inputTree[sam]->SetBranchAddress("_puCorr", _puCorr, &b__puCorr);
    inputTree[sam]->SetBranchAddress("_absIso03", _absIso03, &b__absIso03);
    inputTree[sam]->SetBranchAddress("_absIso04", _absIso04, &b__absIso04);
    inputTree[sam]->SetBranchAddress("_sumNeutralHadronEt04", _sumNeutralHadronEt04, &b__sumNeutralHadronEt04);
    inputTree[sam]->SetBranchAddress("_sumChargedHadronPt04", _sumChargedHadronPt04, &b__sumChargedHadronPt04);
    inputTree[sam]->SetBranchAddress("_sumPhotonEt04", _sumPhotonEt04, &b__sumPhotonEt04);
    inputTree[sam]->SetBranchAddress("_sumNeutralHadronEt03", _sumNeutralHadronEt03, &b__sumNeutralHadronEt03);
    inputTree[sam]->SetBranchAddress("_sumChargedHadronPt03", _sumChargedHadronPt03, &b__sumChargedHadronPt03);
    inputTree[sam]->SetBranchAddress("_sumPhotonEt03", _sumPhotonEt03, &b__sumPhotonEt03);
    inputTree[sam]->SetBranchAddress("_trackIso", _trackIso, &b__trackIso);
    inputTree[sam]->SetBranchAddress("_ecalIso", _ecalIso, &b__ecalIso);
    inputTree[sam]->SetBranchAddress("_hcalIso", _hcalIso, &b__hcalIso);
    //inputTree[sam]->SetBranchAddress("_deltaBIso", _deltaBIso, &b__deltaBIso);
    inputTree[sam]->SetBranchAddress("_ecalPFClusterIso", _ecalPFClusterIso, &b__ecalPFClusterIso);
    inputTree[sam]->SetBranchAddress("_hcalPFClusterIso", _hcalPFClusterIso, &b__hcalPFClusterIso);
    inputTree[sam]->SetBranchAddress("_ptRel", _ptRel, &b__ptRel);
    inputTree[sam]->SetBranchAddress("_ptRatio", _ptRatio, &b__ptRatio);
    inputTree[sam]->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
    inputTree[sam]->SetBranchAddress("_tauMuonVeto", _tauMuonVeto, &b__tauMuonVeto);
    inputTree[sam]->SetBranchAddress("_tauEleVeto", _tauEleVeto, &b__tauEleVeto);
    inputTree[sam]->SetBranchAddress("_decayModeFindingNew", _decayModeFindingNew, &b__decayModeFindingNew);
    inputTree[sam]->SetBranchAddress("_tauVLooseMvaNew", _tauVLooseMvaNew, &b__tauVLooseMvaNew);
    inputTree[sam]->SetBranchAddress("_tauLooseMvaNew", _tauLooseMvaNew, &b__tauLooseMvaNew);
    inputTree[sam]->SetBranchAddress("_tauMediumMvaNew", _tauMediumMvaNew, &b__tauMediumMvaNew);
    inputTree[sam]->SetBranchAddress("_tauTightMvaNew", _tauTightMvaNew, &b__tauTightMvaNew);
    inputTree[sam]->SetBranchAddress("_tauVTightMvaNew", _tauVTightMvaNew, &b__tauVTightMvaNew);
    inputTree[sam]->SetBranchAddress("_tauVTightMvaOld", _tauVTightMvaOld, &b__tauVTightMvaOld);
   
    //    inputTree[sam]->SetBranchAddress("_selectedTrackMult", _selectedTrackMult, &b__selectedTrackMult);
    inputTree[sam]->SetBranchAddress("_lIsPrompt", _lIsPrompt, &b__lIsPrompt);
    inputTree[sam]->SetBranchAddress("_lMatchPdgId", _lMatchPdgId, &b__lMatchPdgId);
    inputTree[sam]->SetBranchAddress("_lProvenance", _lProvenance, &b__lProvenance);
    inputTree[sam]->SetBranchAddress("_lProvenanceCompressed", _lProvenanceCompressed, &b__lProvenanceCompressed);
    inputTree[sam]->SetBranchAddress("_lMatchPt", _lMatchPt, &b__lMatchPt);
    inputTree[sam]->SetBranchAddress("_lMatchEta", _lMatchEta, &b__lMatchEta);
    inputTree[sam]->SetBranchAddress("_lMatchPhi", _lMatchPhi, &b__lMatchPhi);

    inputTree[sam]->SetBranchAddress("_lMatchVertexX", _lMatchVertexX, &b__lMatchVertexX);
    inputTree[sam]->SetBranchAddress("_lMatchVertexY", _lMatchVertexY, &b__lMatchVertexY);
    inputTree[sam]->SetBranchAddress("_lMatchVertexZ", _lMatchVertexZ, &b__lMatchVertexZ);

    inputTree[sam]->SetBranchAddress("_pvX", &_pvX, &b__pvX);
    inputTree[sam]->SetBranchAddress("_pvY", &_pvY, &b__pvY);
    inputTree[sam]->SetBranchAddress("_pvZ", &_pvZ, &b__pvZ);

    inputTree[sam]->SetBranchAddress("_pvXErr", &_pvXErr, &b__pvXErr);
    inputTree[sam]->SetBranchAddress("_pvYErr", &_pvYErr, &b__pvYErr);
    inputTree[sam]->SetBranchAddress("_pvZErr", &_pvZErr, &b__pvZErr);
    inputTree[sam]->SetBranchAddress("_closestJetCsvV2", _closestJetCsvV2, &b__closestJetCsvV2);
    inputTree[sam]->SetBranchAddress("_closestJetDeepCsv_b", _closestJetDeepCsv_b, &b__closestJetDeepCsv_b);
    
    inputTree[sam]->SetBranchAddress("_nPh", &_nPh, &b__nPh);
    inputTree[sam]->SetBranchAddress("_phPt", _phPt, &b__phPt);
    inputTree[sam]->SetBranchAddress("_phEta", _phEta, &b__phEta);
    inputTree[sam]->SetBranchAddress("_phEtaSC", _phEtaSC, &b__phEtaSC);
    inputTree[sam]->SetBranchAddress("_phPhi", _phPhi, &b__phPhi);
    inputTree[sam]->SetBranchAddress("_phE", _phE, &b__phE);
    inputTree[sam]->SetBranchAddress("_phCutBasedLoose", _phCutBasedLoose, &b__phCutBasedLoose);
    inputTree[sam]->SetBranchAddress("_phCutBasedMedium", _phCutBasedMedium, &b__phCutBasedMedium);
    inputTree[sam]->SetBranchAddress("_phCutBasedTight", _phCutBasedTight, &b__phCutBasedTight);
    inputTree[sam]->SetBranchAddress("_phMva", _phMva, &b__phMva);
    inputTree[sam]->SetBranchAddress("_phRandomConeChargedIsolation", _phRandomConeChargedIsolation, &b__phRandomConeChargedIsolation);
    inputTree[sam]->SetBranchAddress("_phChargedIsolation", _phChargedIsolation, &b__phChargedIsolation);
    inputTree[sam]->SetBranchAddress("_phNeutralHadronIsolation", _phNeutralHadronIsolation, &b__phNeutralHadronIsolation);
    inputTree[sam]->SetBranchAddress("_phPhotonIsolation", _phPhotonIsolation, &b__phPhotonIsolation);
    inputTree[sam]->SetBranchAddress("_phSigmaIetaIeta", _phSigmaIetaIeta, &b__phSigmaIetaIeta);
    //inputTree[sam]->SetBranchAddress("_phSigmaIetaIphi", _phSigmaIetaIphi, &b__phSigmaIetaIphi);
    inputTree[sam]->SetBranchAddress("_phHadronicOverEm", _phHadronicOverEm, &b__phHadronicOverEm);
    inputTree[sam]->SetBranchAddress("_phPassElectronVeto", _phPassElectronVeto, &b__phPassElectronVeto);
    inputTree[sam]->SetBranchAddress("_phHasPixelSeed", _phHasPixelSeed, &b__phHasPixelSeed);
    inputTree[sam]->SetBranchAddress("_phIsPrompt", _phIsPrompt, &b__phIsPrompt);
   
    inputTree[sam]->SetBranchAddress("_phMatchPdgId", _phMatchPdgId, &b__phMatchPdgId);
    inputTree[sam]->SetBranchAddress("_nJets", &_nJets, &b__nJets);
    inputTree[sam]->SetBranchAddress("_jetPt", _jetPt, &b__jetPt);
    inputTree[sam]->SetBranchAddress("_jetSmearedPt", _jetSmearedPt, &b__jetSmearedPt);
    inputTree[sam]->SetBranchAddress("_jetSmearedPt_JECUp", _jetSmearedPt_JECUp, &b__jetSmearedPt_JECUp);
    inputTree[sam]->SetBranchAddress("_jetSmearedPt_JECDown", _jetSmearedPt_JECDown, &b__jetSmearedPt_JECDown);
    inputTree[sam]->SetBranchAddress("_jetSmearedPt_JERUp", _jetSmearedPt_JERUp, &b__jetSmearedPt_JERUp);
    inputTree[sam]->SetBranchAddress("_jetSmearedPt_JERDown", _jetSmearedPt_JERDown, &b__jetSmearedPt_JERDown);
    inputTree[sam]->SetBranchAddress("_jetEta", _jetEta, &b__jetEta);
    inputTree[sam]->SetBranchAddress("_jetPhi", _jetPhi, &b__jetPhi);
    inputTree[sam]->SetBranchAddress("_jetE", _jetE, &b__jetE);
    inputTree[sam]->SetBranchAddress("_jetCsvV2", _jetCsvV2, &b__jetCsvV2);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_udsg", _jetDeepCsv_udsg, &b__jetDeepCsv_udsg);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_b", _jetDeepCsv_b, &b__jetDeepCsv_b);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_c", _jetDeepCsv_c, &b__jetDeepCsv_c);
    inputTree[sam]->SetBranchAddress("_jetDeepCsv_bb", _jetDeepCsv_bb, &b__jetDeepCsv_bb);
    inputTree[sam]->SetBranchAddress("_jetHadronFlavor", _jetHadronFlavor, &b__jetHadronFlavor);
    //inputTree[sam]->SetBranchAddress("_jetId", _jetId, &b__jetId);
    inputTree[sam]->SetBranchAddress("_jetIsTight", _jetIsTight, &b__jetIsTight);

    TH1D* _hCounter = new TH1D("hCounter", "Events counter", 5,0,5);
    _hCounter->Read("hCounter");
    hcounter[sam] = _hCounter->GetBinContent(1);

  }//end for on tree
    
    
  //******************* HISTO **********************
  const int nCat=1;  // n. entries in catNames
  const int nSR=12;  // n. search regions 

  const int nDist = 1;  // Number of distributions to plot

  TH1D* Histos[nDist][nCat][nSamples_eff +1];

  // Only to be initialized and filled for QCD and PDF uncertainties
  //  - only one distribution (search regions); nDist = 1 anyway here...
  //  - no need for nCat
  //  - nSamples_eff needed
  //  - QCD scales: need for 6 variations 
  //  - PDFs: need for 100 variations
  std::vector<unsigned> theoSystVars;
  bool runtheosyst = (systcat==1 || systcat==2);
  if(systcat==1) {
    theoSystVars.push_back(2);
    theoSystVars.push_back(3);
    theoSystVars.push_back(4);
    theoSystVars.push_back(5);
    theoSystVars.push_back(7);
    theoSystVars.push_back(9);
  }
  else if(systcat==2) {
    for(unsigned l=10; l<110; ++l)
      theoSystVars.push_back(l);
  }
  const unsigned nTheoVars = theoSystVars.size();
  TH1D* systHistos[nTheoVars][nDist][nSamples_eff+1];

  //const TString catNames[nCat]= {"M2_V0024","M2_V0028","M2_V0031","M2_V0044","M2_V0054","M2_V0070","M2_V0083","M2_V0141","M2_V0173","M2_V01","M2_V0223","M3_V0024","M3_V0028","M3_V0031","M3_V0044","M3_V0054","M3_V0070","M3_V0083","M3_V0141","M3_V0173","M3_V01","M4_V0024","M4_V0028","M4_V0031","M4_V0044","M4_V0054","M4_V0070","M4_V0083","M4_V0141","M4_V01","M5_V0024","M5_V0028","M5_V0031","M5_V0044","M5_V0054","M5_V0070","M5_V0083","M5_V01","M8_V0024","M8_V0028","M8_V0031","M8_V0044","M8_V0054","M10_V0024","M10_V0028","M10_V0031","M10_V0044","M15_V0024","M15_V0028","M15_V0031","M15_V0044","M20_V1933"};
  //const TString catNames[nCat]={"mucoupling","single_fake", "double_fakes", "samejet", "deltaR"};
  const TString catNames[nCat]={"mucoupling"};
  const TString Histnames_ossf[nDist] = {"events" };


  const TString Xaxes[nDist] = {""};


  const TString Units[nDist] = {" "};
  const double HistMin[nDist] = {0.5};

  
  
  const double HistMax[nDist] = { 12.5};
  const int nBins[nDist] =      { 12};


  
			       
  cout<<"------ 1"<<endl;
  for(int i = 0; i < nDist; ++i){
    float BinWidth = (HistMax[i] - HistMin[i])/nBins[i];
    std::ostringstream strs; strs << BinWidth; std::string Yaxis = strs.str();
    for(int effsam = 0; effsam < nSamples_eff + 1; ++effsam){
      for(int cat = 0; cat < nCat; ++cat){               
	Histos[i][cat][effsam] = new TH1D(eff_names[effsam] + catNames[cat] + Histnames_ossf[i] , eff_names[effsam] + catNames[cat] + Histnames_ossf[i] + ";" + Xaxes[i] + "; events /" + Yaxis + Units[i], nBins[i], HistMin[i], HistMax[i]);
	Histos[i][cat][effsam]->Sumw2();
      }

      // Only for theory systs
      if(runtheosyst) {
	for(unsigned sidx=0; sidx<nTheoVars; ++sidx) {
	  systHistos[sidx][i][effsam] = new TH1D(eff_names[effsam] + "_syst_" + theoSystVars[sidx] + "_" + Histnames_ossf[i] , eff_names[effsam] + "_syst_" + theoSystVars[sidx] + "_" + Histnames_ossf[i] + ";" + Xaxes[i] + "; events /" + Yaxis + Units[i], nBins[i], HistMin[i], HistMax[i]);
	  systHistos[sidx][i][effsam]->Sumw2();
	}
      }
    }
  }
  // Calculate the center of the maximum bin of each histogram
  double maxBinC[nDist];
  for(int i = 0; i < nDist; ++i){
    maxBinC[i] = Histos[i][0][0]->GetBinCenter(Histos[i][0][0]->GetNbinsX());
  }
  cout<<"--2"<<endl;


  
  
  


  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> PARAMETERS AND CUTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  const double met_cuts =80;
  const double mlll_cuts = 85;
  const int number_veto_leptons=3;
  const double b_jets_wp= 0.5426;
  const double b_jets_pt= 25;

  const double MVA_cuts_pt15[3] = {0.77, 0.56, 0.48};
  const double MVA_cuts_pt25[3] = {0.52, 0.11, -0.01};
  const double newMVALooseFR[3]= {-0.02, -0.52, -0.52}; 

  

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
  //Loop over all samples
  Double_t Counter[nSamples];
  for(int i = 0; i  < nSamples; ++i){
    Counter[i] = 0;
  }
  Double_t scale[nSamples];
    
  double tot_ttbar=0;
  double heavy_ttbar=0;
  double light_ttbar=0;




  
  



  
  //-------------->  LOOP ON ALL SAMPLES
  for(int sam = 0, effsam = 0; sam < nSamples; ++sam, ++effsam){
    Long64_t nEntries = inputTree[sam]->GetEntries();
    cout<<"----------------"<<endl;
    if(sam != 0){
      cout<<"=========  effsam: "<<effsam<<" "<<names[sam]<<"   "<<sam<<endl;
      if(names[sam] == names[sam -1]) --effsam;
      cout<<"+++++++++  effsam: "<<effsam<<" "<<names[sam]<<"   "<<sam<<endl;

    }
    if(sam >= 0){ // 1 data set
      scale[sam] = xSections[sam]*luminosity*1000/(hcounter[sam]);
    }
    if (effsam == 0) continue;
    std::cout<<"Entries in "<< fileList[sam] <<" "<<nEntries<<std::endl;
    cout << effsam << endl;
    //if (effsam != 13) continue;
    //  if ( fileList[sam] != "TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root") continue;

    
    
    //if (effsam !=1) continue;
    //  if (sam > 7) continue;
        
    double counters_cut[20] ;
    double counters_cut_ossf[20] ;
    double counters_cut_no_ossf[20] ;


    // if (sam <= 12 ){

    for (int i =0; i< 20; i++){
      counters_cut[i] =0.;
      counters_cut_ossf[i] =0.;
      counters_cut_no_ossf[i] =0.;
    }
    // }
    double progress = 0; 	//For printing progress bar
    

bool dy_ttt=false;
  bool dy_ttt_df=false;
  bool dy_ttt_cj=false;
  bool dy_ttt_dr=false;

  bool tt_ttt=false;
  bool tt_ttt_df=false;
  bool tt_ttt_cj=false;
  bool tt_ttt_dr=false;
 
  bool left_ttt=false;
  bool left_ttt_df=false;
  bool left_ttt_cj=false;
  bool left_ttt_dr=false;
    
        
    //--------------> LOOOP ON THE TREE"S ENTRIES
    for (Long64_t it = 0; it < nEntries; ++it){
      inputTree[sam]->GetEntry(it);
      
      //if (it%10000 == 0) cout<<'.'<<flush;
            
      if(it%100 == 0 && it != 0){
	progress += (double) (100./nEntries);
	printProgress(progress);
      } else if(it == nEntries -1){
	progress = 1.;
	printProgress(progress);
      }
            
            
      double scal = 0;
      scal = scale[sam]*_weight * pu_weight(*&pileUpWeight[0],_nTrueInt);
      //
      // Weight for b-jet veto
      double bwght = 1.;
     
            
      if (effsam == 0) continue;
      if (effsam == 69) continue;   
            
      // if(_nL < number_veto_leptons) continue;
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> VECTORS AND VARIABLES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      int               nBjets = 0;
      unsigned*         ind = new unsigned[_nL];	//new indices of good leptons
      //double*           conePt = new double[_nL];
      double            faxtore_FR=0;            
      double            prov_index[3]={-1,-1,-1};
      double            prov_number_tight[1]= {-1};
      double            prov_number_prompt[1]= {-1};
      double            prov_fattore[1]= {-1};
      int               skip_event[1]= {-1};
      double            faxtore[1]= {-100};            
      TLorentzVector    lepton_reco[3];
      int               flavors_3l[3];
      int               charge_3l[3];
      TLorentzVector    sum_3l_rec;	//M_3l
      TLorentzVector    pair [3];
      int               kind[1] ={-1}; // 0 no-ossf
      TLorentzVector    sum_2l_rec_pair; 	//M_2l best Z candidate
      int               event_clas[1]={-1}; // 1* = eee   2* = emm   3* = eem  4* = mmm
      int               check_mt= -1;
      int               check_deltaR= -1;
      int               check_mt_os= -1;
      Double_t          delta_R_max=-1;
      Double_t          delta_R_min=-1;
      Double_t          _mll_min=50000;
      Double_t          _mll_min_os=50000;
      TLorentzVector    lepton_transv[3];
      TLorentzVector    METvec;
      double            m_T=0;
      double            MET=0;
      double            MET_PHI=0;           
      int               nominal=-1;
      int               jec_check=-1;
      int               unclustered_met=-1;
      int               up=-1;
      int               down=-1;
      double            new_met[1]= {0};
      double            new_met_phi[1]= {0};
      int               new_number_bjets[1]= {0};            
      int               new_pt_checks= -1;
      bool              trigger_fired = false;
      bool              low_pt_event=false;
      bool              high_pt_event = false;            
      unsigned*         _SR_low= new unsigned[8];
      double             search_region_fill[1]={-1};
      bool              data_control_region=false;



      unsigned          lCount = 0;	//Count number of FO leptons that are not taus
      unsigned*         _isFO= new unsigned[_nL];
      Bool_t            _passedMVA90[_nL];   
      
     
      unsigned           promptC = 0;
      double             low_mass_pt_base[1];
      
      unsigned*          _isWithTrack= new unsigned[_nL];
      unsigned           wTrack=0;

      double            iV_ls=0;
      double            iV_lt=0;
      double            iV_st=0;

      double            _vertex_X[3];
      double            _vertex_Y[3];
      double            _vertex_Z[3];

      double            _vertex_sX[3];
      double            _vertex_sY[3];
      double            _vertex_sZ[3];
      
      double            _vertex_R2D[3];
      double            _vertex_sR2D[3];
      double            _vertex_R[3];
      double            _vertex_sR[3];
      double            _vertex_chi2[3];
      double            _vertex_normchi2[3];

      unsigned         ind_new_leading=0;
      unsigned         ind_new_p=0;
      unsigned         ind_new_pp=0;

      bool             isAll = false;
      bool             isAll_met_mlll = false;

      unsigned*         _isLooseCutBasedElectronWithoutIsolatio= new unsigned[_nL];
      unsigned*         _isTightCutBasedElectronWithoutIsolatio= new unsigned[_nL];
      unsigned*         _isMediumCutBasedElectronWithoutIsolatio= new unsigned[_nL];
      unsigned*         _isOurMedium= new unsigned[_nL];
      unsigned*         _passTimingVeto= new unsigned[_nL];

      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<            
      const double isolation_loose=1.2;
      const double isolation_tight=0.2;
      //------------------------------------------------------------ selection class
     
      double*  conePt = new double[_nL];

      int goodjet=0;
      int bjet=0;

      for (unsigned j =0; j < _nJets ; j++){
	TLorentzVector jet;
	// jet.SetPtEtaPhiE(_jetPt[j],_jetEta[j],_jetPhi[j],_jetE[j]);
	jet.SetPtEtaPhiE(_jetSmearedPt[j],_jetEta[j],_jetPhi[j],_jetE[j]);
	if(systcat==8) {
	  if(systdir==0) {
	    jet.SetPtEtaPhiE(_jetSmearedPt_JECDown[j],_jetEta[j],_jetPhi[j],_jetE[j]);
	  }
	  else {
	    jet.SetPtEtaPhiE(_jetSmearedPt_JECUp[j],_jetEta[j],_jetPhi[j],_jetE[j]);
	  }
	}
	else if(systcat==9) {
	  if(systdir==0) {
	    jet.SetPtEtaPhiE(_jetSmearedPt_JERDown[j],_jetEta[j],_jetPhi[j],_jetE[j]);
	  }
	  else {
	    jet.SetPtEtaPhiE(_jetSmearedPt_JERUp[j],_jetEta[j],_jetPhi[j],_jetE[j]);
	  }
	}
	bool ciccia = true;
	for(unsigned i = 0; i < _nL; ++i){
	  if (!_lPOGLoose[i]) continue;
	  TLorentzVector leptone;
	  leptone.SetPtEtaPhiE(_lPt[i],_lEta[i],_lPhi[i],_lE[i]);
	  if (leptone.DeltaR(jet) < 0.4) ciccia = false;
	}
	// if (ciccia && _jetIsTight[j] && _jetEta[j] < 2.4 && _jetPt[j] > 25) goodjet++;
	// if (ciccia && _jetIsTight[j] && _jetEta[j] < 2.4 && _jetPt[j] > 25 && _jetCsvV2[j] > 0.546) bjet++;

	if(ciccia && _jetIsTight[j] && _jetEta[j]<2.4 && jet.Pt()>25.) goodjet++;
	if(ciccia && _jetIsTight[j] && _jetEta[j]<2.4 && jet.Pt()>25.
	   && (_jetDeepCsv_b[j]+_jetDeepCsv_bb[j])>btagCuts[year][bwp]) {
	  bjet++;
	  if(jet.Pt()<1000. && std::abs(_jetEta[j])<2.4) {
	    double bjetSf = 1.;
	    // b-jet systematics
	    if(systcat==10) {
	      if(systdir==0) {
		bjetSf = reader.eval_auto_bounds("down", BTagEntry::FLAV_B, std::abs(_jetEta[j]), jet.Pt());
	      }
	      else {
		bjetSf = reader.eval_auto_bounds("up", BTagEntry::FLAV_B, std::abs(_jetEta[j]), jet.Pt());
	      }
	    }
	    // b-jet central SF
	    else {
	      bjetSf = reader.eval_auto_bounds("central", BTagEntry::FLAV_B, std::abs(_jetEta[j]), jet.Pt());
	    }
	    //
	    // Scale the b-veto event weight
	    bwght *= bjetSf;
	  }
	}
      } // end loop on jets: for (unsigned j =0; j < _nJets ; j++)


      // ------------ ==================== -----------------------------------------------//
      // ------------ muon and electron ID -----------------------------------------------//
      for(unsigned i = 0; i < _nL; ++i){
	
	conePt[i]=_lPt[i];
	//conePt[i]= _lPt[i]* maximum (0, _relIso[i] - 0.2);
	//_lPt[i] = _lPt[i]* maximum (0, _relIso[i] - 0.2);
	//_lE[i] = _lE[i] *  maximum (0, _relIso[i] - 0.2);
	if (_lFlavor[i] == 0) {
	  _isLooseCutBasedElectronWithoutIsolatio[i] = true;
	  if (!(_lEleIsEB[i] || _lEleIsEE[i])) _isLooseCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lElefull5x5SigmaIetaIeta[i]                  >= (_lEleIsEB[i] ? 0.11       : 0.0314 ))       _isLooseCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleDEtaInSeed  [i]                          >= (_lEleIsEB[i] ? 0.00477    : 0.00868))       _isLooseCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleDeltaPhiSuperClusterTrackAtVtx [i]       >= (_lEleIsEB[i] ? 0.222      : 0.213  ))       _isLooseCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lElehadronicOverEm[i]                        >= (_lEleIsEB[i] ? 0.298      : 0.101  ))       _isLooseCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleInvMinusPInv[i]                          >= (_lEleIsEB[i] ? 0.241      : 0.14   ))       _isLooseCutBasedElectronWithoutIsolatio[i] =  false;

	  _isMediumCutBasedElectronWithoutIsolatio[i] = true;
	  if (!(_lEleIsEB[i] || _lEleIsEE[i])) _isMediumCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lElefull5x5SigmaIetaIeta[i]                  >= (_lEleIsEB[i] ? 0.00998    : 0.0298 ))       _isMediumCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleDEtaInSeed  [i]                          >= (_lEleIsEB[i] ? 0.00311    : 0.00609))       _isMediumCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleDeltaPhiSuperClusterTrackAtVtx [i]       >= (_lEleIsEB[i] ? 0.103      : 0.045  ))       _isMediumCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lElehadronicOverEm[i]                        >= (_lEleIsEB[i] ? 0.253      : 0.0878  ))      _isMediumCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleInvMinusPInv[i]                          >= (_lEleIsEB[i] ? 0.134      : 0.13   ))       _isMediumCutBasedElectronWithoutIsolatio[i] =  false;
	  
	  _isTightCutBasedElectronWithoutIsolatio[i] = true;
	  if (!(_lEleIsEB[i] || _lEleIsEE[i])) _isTightCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lElefull5x5SigmaIetaIeta[i]                  >= (_lEleIsEB[i] ? 0.00998    : 0.0292 ))       _isTightCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleDEtaInSeed  [i]                          >= (_lEleIsEB[i] ? 0.00308    : 0.00605))       _isTightCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleDeltaPhiSuperClusterTrackAtVtx [i]       >= (_lEleIsEB[i] ? 0.0816      : 0.0394  ))     _isTightCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lElehadronicOverEm[i]                        >= (_lEleIsEB[i] ? 0.0414      : 0.0641  ))     _isTightCutBasedElectronWithoutIsolatio[i] =  false;
	  if(_lEleInvMinusPInv[i]                          >= (_lEleIsEB[i] ? 0.0129      : 0.0129   ))    _isTightCutBasedElectronWithoutIsolatio[i] =  false;



	  _passedMVA90[i]=false;
	  int eta = -1;
	  if(TMath::Abs(_lEta[i]) < 0.8 ) eta = 0;
	  else if(TMath::Abs(_lEta[i]) < 1.479 ) eta = 1;
	  else eta = 2;
	  _passedMVA90[i] = _lElectronMva[i] >  std::min( MVA_cuts_pt15[eta], std::max(MVA_cuts_pt25[eta] , MVA_cuts_pt15[eta] + (MVA_cuts_pt25[eta] - MVA_cuts_pt15[eta])*0.1 *( _lPt[i]-15) ) );
	  
	}// end ele


	
	if (_lFlavor[i] == 1){
	  _isOurMedium[i] = false;
	  bool goodGlob = false;
	  goodGlob = _lGlobalMuon[i] && _lCQChi2Position[i] < 12   &&  _lCQTrackKink[i] < 20;
          _isOurMedium[i] = _lPOGLoose[i] &&  _muonSegComp[i] > (goodGlob ? 0.303 : 0.451);

	  
	  //time veto
	  _passTimingVeto[i] = true;
	  bool cmbok =(_lMuTimenDof[i] >7   );
	  bool rpcok =(_lMuRPCTimenDof[i] >1  && _lMuRPCTimeErr[i]==0);
	  if (rpcok){
	    if ((fabs(_lMuRPCTime[i])>10) &&	!(cmbok && fabs(_lMuTime[i])<10))   _passTimingVeto[i]=false;
	  }
	  else{
	    if (cmbok && (_lMuTime[i]>20 || _lMuTime[i]<-45)) _passTimingVeto[i]=false;
	  }
	}// end muon
      }//end new ID
      // ------------ ==================== -----------------------------------------------//


      // ------------ ==================== -----------------------------------------------//
      // ------------   object selection   -----------------------------------------------//
      for(unsigned l = 0; l < _nL; ++l){
	_isFO[l] = false;
		if (fabs(_dz[l]) > 10) continue;

	if (_lFlavor[l] == 0 && _lPt[l] < 10 ) continue;
	if (_lFlavor[l] == 0 && _relIso[l] > isolation_loose) continue;
	if (_lFlavor[l] == 0 && !_isLooseCutBasedElectronWithoutIsolatio[l]) continue;
	

	if (_lFlavor[l] == 1 && _relIso[l] > isolation_loose)   continue;
	if (_lFlavor[l] == 1 && _lPt[l] < 5 ) continue;
	if (_lFlavor[l] == 1 && !_isOurMedium[l]) continue;
	//if (_lFlavor[l] == 1 && !_passTimingVeto[l]) continue;

	_isFO[l] = true;
	if(_isFO[l] && _lFlavor[l] != 2){
	  ind[lCount] = l;
	  ++lCount;   // number of FO letpons (not taus)
	}
      } 
  
      if(lCount < 3) continue;
      //Order FO leptons by Pt
      unsigned* ordind = new unsigned[lCount];
      std::set<unsigned> usedLep;      
      for(unsigned k =0; k < lCount; ++k){
	//unsigned maxI = 999;
	double maxConePt = 0;
	for(unsigned l = 0; l < lCount; ++l){
	  if(usedLep.find(ind[l]) == usedLep.end()){
	    if(conePt[ind[l]] > maxConePt){
	      maxConePt = conePt[ind[l]];
	      ordind[k] = ind[l];
	    }
	  }
	}
	usedLep.insert(ordind[k]);
      }
      for(unsigned i = 0; i < lCount; ++i){
	ind[i] = ordind[i];
      }
      double pippa0=_lPt[ind[0]];
      double pippa1=_lPt[ind[1]];
      double pippa2=_lPt[ind[2]];

      // ------------ ==================== -----------------------------------------------//


      // ------------ ==================== -----------------------------------------------//
      // ------------   event selection   -----------------------------------------------//
      ind_new_leading = -1;
      int counter_leading=0;
      //*****     ---------------->      leading
      for(unsigned l = 0; l < lCount; ++l){
	if (counter_leading == 0){
	  if ((_lFlavor[ind[l]]== 1 && _lPOGLoose[ind[l]] && _lPOGMedium[ind[l]] && _relIso[ind[l]] < 0.1 && TMath::Abs(_dxy[ind[l]]) < 0.05 && TMath::Abs(_dz[ind[l]]) < 0.1 && fabs(_3dIPSig [ind[l]]) < 4 && _lPt[ind[l]]> 25 )                           ||                          (_lFlavor[ind[l]]== 0 &&  _passedMVA90[ind[l]] && _relIso[ind[l]] < 0.1 && TMath::Abs(_dxy[ind[l]]) < 0.05 && TMath::Abs(_dz[ind[l]]) < 0.1 && fabs(_3dIPSig [ind[l]]) < 4 && _lPt[ind[l]]> 30)) {
	    ++counter_leading;
	    ind_new_leading = ind[l];
	  }
	}
      }
      if (counter_leading <1)continue;
      for (unsigned i =0; i < lCount; i ++){
	if (ind[i] == ind_new_leading) {
	  lepton_reco[0].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);
	  flavors_3l[0]=_lFlavor[ind[i]];
	  charge_3l[0]=_lCharge[ind[i]];
	  conePt[0] =  _lPt[ind[i]];
	}
      }
      //*****     ---------------->      displaced pair
      unsigned displacedC = 0;
      unsigned* _isDisplaced= new unsigned[_nL];
      for(unsigned l = 0; l < lCount; ++l){
	_isDisplaced[ind[l]] = false;
	if (ind[l] == ind_new_leading ) continue;
	int number_possible_vertices= 0;
	if (TMath::Abs(_dxy[ind[l]]) > 0.01 ) {
	  for(unsigned secondo = 0; secondo < lCount; ++secondo){
	    if (ind[l] == ind[secondo]) continue;
	    for(unsigned v = 0; v < _nVFit; ++v){
	      if ((_vertices[v][0] == (_lIndex[ind[l]] * 100 +  _lIndex[ind[secondo]] )) ||(_vertices[v][0] == (_lIndex[ind[l]] +  _lIndex[ind[secondo]] *100) )) number_possible_vertices++;
	    }
	  }
	}//dispalced
	if (TMath::Abs(_dxy[ind[l]]) > 0.01  && number_possible_vertices>0)	_isDisplaced[ind[l]] = true;
      }     
      for(unsigned l = 0; l < lCount; ++l){
	if(_isDisplaced[ind[l]]) ++displacedC;	
      }
      if (displacedC < 2) continue;


      if (flavors_3l[0] == 1 && ((_lHasTrigger[ind_new_leading] & 1)==0 && (_lHasTrigger[ind_new_leading] & 2)==0)) continue;
      if (flavors_3l[0] == 0 && ((_lHasTrigger[ind_new_leading] & 1)==0 )) continue;

  
      //  ----------->      easy case == just 2 lepton
      TLorentzVector   lepton_tobeselected[20];
      int              index_displaced[20];
      int index_l[2];
      int vertex_only_2_found= 0;
      TVector3 l2_position[1];
      TVector3 l3_position[1];
    
      TLorentzVector leptons_propagated[3];
      leptons_propagated[0].SetPtEtaPhiE(lepton_reco[0].Pt(),0, lepton_reco[0].Phi(), lepton_reco[0].Pt());
      leptons_propagated[1].SetPtEtaPhiE(lepton_reco[1].Pt(),0, lepton_reco[1].Phi(), lepton_reco[1].Pt());
      leptons_propagated[2].SetPtEtaPhiE(lepton_reco[2].Pt(),0, lepton_reco[2].Phi(), lepton_reco[2].Pt());
       
      TVector3 position_propagated[2];
      // more than 2 leptons!!! 
      //   if (displacedC > 2){
	
      for (int i =0; i < 20; i ++){
	index_displaced[i] = -50;
	lepton_tobeselected[i].SetPtEtaPhiE( 0,  0, 0, 0);
      }
      int displacedC_check=0;
      for (unsigned i =0; i < lCount; i ++){
	if (_isDisplaced[ind[i]]){	    
	  lepton_tobeselected[displacedC_check].SetPtEtaPhiE( _lPt[ind[i]],  _lEta[ind[i]], _lPhi[ind[i]], _lE[ind[i]]);
	  index_displaced[displacedC_check] = ind[i];
	  displacedC_check++;
	}
      }

      int index_1=-5;
      int index_2= -5;
      for (unsigned h = 0; h <displacedC; h ++ ){
	for (unsigned g = 0; g < displacedC; g ++){
	  if (h != g) {
	    int find_vertex=0;
	    double vertice_st= -1;
	    for(unsigned v = 0; v < _nVFit; ++v){
	      if ((_vertices[v][0] == (_lIndex[index_displaced[h]] * 100 +  _lIndex[index_displaced[g]] )) ||(_vertices[v][0] == (_lIndex[index_displaced[h]] +  _lIndex[index_displaced[g]] *100) ))	{
		find_vertex=1;
		vertice_st= v;
	      }
	    }//vertex
	    if ( find_vertex== 0) continue;

	    TLorentzVector leptons_propagated_pair[3];
	    leptons_propagated_pair[0].SetPtEtaPhiE(0,0,0,0);
	    leptons_propagated_pair[1].SetPtEtaPhiE(_lPt[index_displaced[h]],0,_lPhi[index_displaced[h]], _lPt[index_displaced[h]]);
	    leptons_propagated_pair[2].SetPtEtaPhiE(_lPt[index_displaced[g]],0,_lPhi[index_displaced[g]], _lPt[index_displaced[g]]);
	    for(unsigned v = 0; v < _nVFit; ++v){
	      if (v == vertice_st){
		if (_vertices[v][0] == (_lIndex[index_displaced[h]] * 100 +  _lIndex[index_displaced[g]] )){
		  double x1= _lDisplaced[v][3];
		  double y1= _lDisplaced[v][4];
		  double z1= _lDisplaced[v][5];
		  double x2= _lDisplaced[v][15];
		  double y2= _lDisplaced[v][16];
		  double z2= _lDisplaced[v][17];
		  position_propagated[0].SetXYZ(_lDisplaced[v][0],_lDisplaced[v][1],_lDisplaced[v][2]);
		  position_propagated[1].SetXYZ(_lDisplaced[v][12],_lDisplaced[v][13],_lDisplaced[v][14]);

		  if (_lFlavor[index_displaced[h]] ==1)leptons_propagated_pair[1].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.105*0.105) );
		  if (_lFlavor[index_displaced[h]] ==0)leptons_propagated_pair[1].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.0005*0.0005) );
		  if (_lFlavor[index_displaced[g]] ==1)leptons_propagated_pair[2].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.105*0.105) );
		  if (_lFlavor[index_displaced[g]] ==0)leptons_propagated_pair[2].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.0005*0.0005) );
		}// 203
		if (_vertices[v][0] == (_lIndex[index_displaced[h]] +  _lIndex[index_displaced[g]] *100)){
		  double x1= _lDisplaced[v][3];
		  double y1= _lDisplaced[v][4];
		  double z1= _lDisplaced[v][5];
		  double x2= _lDisplaced[v][15];
		  double y2= _lDisplaced[v][16];
		  double z2= _lDisplaced[v][17];
		  if (_lFlavor[index_displaced[g]] ==1)leptons_propagated_pair[2].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.105*0.105) );
		  if (_lFlavor[index_displaced[g]] ==0)leptons_propagated_pair[2].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.0005*0.0005) );
		  if (_lFlavor[index_displaced[h]] ==1)leptons_propagated_pair[1].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.105*0.105) );
		  if (_lFlavor[index_displaced[h]] ==0)leptons_propagated_pair[1].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.0005*0.0005) );
		  position_propagated[1].SetXYZ(_lDisplaced[v][0],_lDisplaced[v][1],_lDisplaced[v][2]);
		  position_propagated[0].SetXYZ(_lDisplaced[v][12],_lDisplaced[v][13],_lDisplaced[v][14]);
		}// 203
	      }// only st
	    }	 
        
	    if (h == 0 && g ==1) {
	      _mll_min = (leptons_propagated_pair[1]+leptons_propagated_pair[2]).M();
	      index_1 = index_displaced[h];
	      index_2 = index_displaced[g];
	      l2_position[0].SetXYZ(position_propagated[0].X(), position_propagated[0].Y(),position_propagated[0].Z());
	      l3_position[0].SetXYZ(position_propagated[1].X(), position_propagated[1].Y(),position_propagated[1].Z());

	      if (leptons_propagated_pair[1].Pt() < leptons_propagated_pair[2].Pt()) {
		index_1 = index_displaced[g];
		index_2 = index_displaced[h];
		l2_position[0].SetXYZ(position_propagated[1].X(), position_propagated[1].Y(),position_propagated[1].Z());
		l3_position[0].SetXYZ(position_propagated[0].X(), position_propagated[0].Y(),position_propagated[0].Z());
	      }
            
	    }
	    if((leptons_propagated_pair[1]+leptons_propagated_pair[2]).M() < _mll_min) {
	      _mll_min = (leptons_propagated_pair[1]+leptons_propagated_pair[2]).M();
	      index_1 = index_displaced[h];
	      index_2 = index_displaced[g];
	      l2_position[0].SetXYZ(position_propagated[0].X(), position_propagated[0].Y(),position_propagated[0].Z());
	      l3_position[0].SetXYZ(position_propagated[1].X(), position_propagated[1].Y(),position_propagated[1].Z());
	      if (leptons_propagated_pair[1].Pt() < leptons_propagated_pair[2].Pt()) {
		index_1 = index_displaced[g];
		index_2 = index_displaced[h];
		l2_position[0].SetXYZ(position_propagated[1].X(), position_propagated[1].Y(),position_propagated[1].Z());
		l3_position[0].SetXYZ(position_propagated[0].X(), position_propagated[0].Y(),position_propagated[0].Z());
	      }

	    }
	  }// h diverso g
	}// loop1
      }//loop 2
      index_l [0] = index_1;
      index_l[1]  = index_2; 


      if (index_1 != -5 && index_2!=-5){
	lepton_reco[1].SetPtEtaPhiE( _lPt[index_l[0]],  _lEta[index_l[0]], _lPhi[index_l[0]], _lE[index_l[0]]);
	flavors_3l[1]=_lFlavor[index_l[0]];
	charge_3l[1]=_lCharge[index_l[0]];
	conePt[1] =  _lPt[index_l[0]];
	lepton_reco[2].SetPtEtaPhiE( _lPt[index_l[1]],  _lEta[index_l[1]], _lPhi[index_l[1]], _lE[index_l[1]]);
	flavors_3l[2]=_lFlavor[index_l[1]];
	charge_3l[2]=_lCharge[index_l[1]];
	conePt[2] =  _lPt[index_l[1]];
      }
      //}// end >2
      
      if (index_l [0] == -5 || index_l [1]== -5) continue;
      double posx_vtx_good=0;
      double posy_vtx_good=0;
      double posz_vtx_good=0;
      for(unsigned v = 0; v < _nVFit; ++v){
	if ((_vertices[v][0] == (_lIndex[index_l[0]] * 100 +  _lIndex[index_l[1]] )) ||(_vertices[v][0] == (_lIndex[index_l[0]] +  _lIndex[index_l[1]] *100) )) {
	  posx_vtx_good        = _vertices[v][1];
	  posy_vtx_good        = _vertices[v][2];
	  posz_vtx_good        = _vertices[v][3]; 	  
	}
      }//end vertici

      
    
      double vertice_st= -1;
      for(unsigned v = 0; v < _nVFit; ++v){
	if ((_vertices[v][0] == (_lIndex[index_l[0]] * 100 +  _lIndex[index_l[1]] )) ||(_vertices[v][0] == (_lIndex[index_l[0]] +  _lIndex[index_l[1]] *100) )) {
	  vertice_st= v;
	}
      }
      for(unsigned v = 0; v < _nVFit; ++v){
	if (v == vertice_st){
	  if (_vertices[v][0] == (_lIndex[index_l[0]] * 100 +  _lIndex[index_l[1]] )){
	 

	    double x1= _lDisplaced[v][3];
	    double y1= _lDisplaced[v][4];
	    double z1= _lDisplaced[v][5];
	    double x2= _lDisplaced[v][15];
	    double y2= _lDisplaced[v][16];
	    double z2= _lDisplaced[v][17];
	    if (_lFlavor[index_l[0]] ==1)leptons_propagated[1].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.105*0.105) );
	    if (_lFlavor[index_l[0]] ==0)leptons_propagated[1].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.0005*0.0005) );
	    if (_lFlavor[index_l[1]] ==1)leptons_propagated[2].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.105*0.105) );
	    if (_lFlavor[index_l[1]] ==0)leptons_propagated[2].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.0005*0.0005) );
	  }// 203
	  if (_vertices[v][0] == (_lIndex[index_l[0]] +  _lIndex[index_l[1]] *100)){

	     
	    double x1= _lDisplaced[v][3];
	    double y1= _lDisplaced[v][4];
	    double z1= _lDisplaced[v][5];
	    double x2= _lDisplaced[v][15];
	    double y2= _lDisplaced[v][16];
	    double z2= _lDisplaced[v][17];
	    if (_lFlavor[index_l[1]] ==1)leptons_propagated[2].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.105*0.105) );
	    if (_lFlavor[index_l[1]] ==0)leptons_propagated[2].SetXYZT(x1,y1,z1,TMath::Sqrt(x1*x1+y1*y1+z1*z1+0.0005*0.0005) );
	    if (_lFlavor[index_l[0]] ==1)leptons_propagated[1].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.105*0.105) );
	    if (_lFlavor[index_l[0]] ==0)leptons_propagated[1].SetXYZT(x2,y2,z2,TMath::Sqrt(x2*x2+y2*y2+z2*z2+0.0005*0.0005) );	   
	  }// 203
	}// only s
      }



      TLorentzVector lepton_prima[2];
      lepton_prima[0].SetPtEtaPhiE(lepton_reco[1].Pt(), lepton_reco[1].Eta(), lepton_reco[1].Phi(), lepton_reco[1].E());
      lepton_prima[1].SetPtEtaPhiE(lepton_reco[2].Pt(), lepton_reco[2].Eta(), lepton_reco[2].Phi(), lepton_reco[2].E());

      //lepton_reco[1].SetPtEtaPhiE( leptons_propagated[1].Pt(),  leptons_propagated[1].Eta(), leptons_propagated[1].Phi(), leptons_propagated[1].E());
      //lepton_reco[2].SetPtEtaPhiE( leptons_propagated[2].Pt(),  leptons_propagated[2].Eta(), leptons_propagated[2].Phi(), leptons_propagated[2].E());
    
      // ------------ ==================== -----------------------------------------------//
      // ------------   tight selection   -----------------------------------------------//
      // ------------ ==================== -----------------------------------------------//
      // ------------   tight selection   -----------------------------------------------//
      unsigned* _isT= new unsigned[_nL];

      for(unsigned l = 0; l < lCount; ++l){
	_isT[ind[l]] = false;
	if (ind[l] == ind_new_leading) continue;
	if (int(ind[l]) != index_l[1] && int(ind[l]) != index_l[0]) continue;
	
	if(_lFlavor[ind[l]]== 1 && _lPt[ind[l]] < 5) continue; 
	if(_lFlavor[ind[l]]== 0 && _lPt[ind[l]] < 10) continue;

	if (_relIso[ind[l]] > isolation_tight) continue;
	//if (fabs(_dxy[ind[l]]) > 0.01) continue;
	if (_lFlavor[ind[l]] == 0 &&  _relIso[ind[l]] <= isolation_tight && _isLooseCutBasedElectronWithoutIsolatio[ind[l]])  _isT[ind[l]] = true;
	if (_lFlavor[ind[l]] == 1 && _relIso[ind[l]] <= isolation_tight && _isOurMedium[ind[l]] && _passTimingVeto[ind[l]])   _isT[ind[l]] = true;
      }


     
      /////////////////////////////// --------------- ////////////////////////
      /////////////////////////////// --------------- ////////////////////////
      promptC=0;
      if (_lIsPrompt[ind_new_leading]) promptC++;
      if (_lIsPrompt[index_l[1]])promptC++;
      if (_lIsPrompt[index_l[0]])promptC++;


      bool dy_to_be_skiped= true;
      if (fileList[sam] == "DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root" || fileList[sam] == "DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root"){
	if (_lMatchPdgId[ind_new_leading] == 22 && _lIsPrompt[ind_new_leading] )  dy_to_be_skiped = false;
	if (_lMatchPdgId[index_l[1]] == 22 && _lIsPrompt[index_l[1]] )  dy_to_be_skiped = false;
	if (_lMatchPdgId[index_l[0]] == 22 && _lIsPrompt[index_l[0]] )  dy_to_be_skiped = false;
      }
      if (!dy_to_be_skiped) continue;

      bool zgamma_to_be_skiped= false;
      if (fileList[sam] == "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root"){
	if (_lMatchPdgId[ind_new_leading] == 22 && _lIsPrompt[ind_new_leading] )  zgamma_to_be_skiped = true;
	if (_lMatchPdgId[index_l[1]] == 22 && _lIsPrompt[index_l[1]] )  zgamma_to_be_skiped = true;
	if (_lMatchPdgId[index_l[0]] == 22 && _lIsPrompt[index_l[0]] )  zgamma_to_be_skiped = true;
      }
      if (fileList[sam] == "ZGTo2LG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8.root" && !zgamma_to_be_skiped) continue;


      double prov1=-1;
      double prov2=-1;
      double prov3=-1;
      if(flavors_3l[0] == 1) prov1 =  _lProvenanceCompressed[ind_new_leading] ;
      if(flavors_3l[1] == 1) prov2 =  _lProvenanceCompressed[index_l[0]];
      if(flavors_3l[2] == 1) prov3 =  _lProvenanceCompressed[index_l[1]];
      if(flavors_3l[0] == 0) prov1 =  5+_lProvenanceCompressed[ind_new_leading] ;
      if(flavors_3l[1] == 0) prov2 =  5+_lProvenanceCompressed[index_l[0]];
      if(flavors_3l[2] == 0) prov3 =  5+_lProvenanceCompressed[index_l[1]];

       

      for(unsigned v = 0; v < _nVFit; ++v){
	if ((_vertices[v][0] == (_lIndex[ind_new_leading] * 100 +  _lIndex[index_l[0]])) || (_vertices[v][0] == (_lIndex[ind_new_leading]  +  _lIndex[index_l[0]]*100))) {
	  _vertex_X[0]        = _vertices[v][1];
	  _vertex_Y[0]        = _vertices[v][2];
	  _vertex_Z[0]        = _vertices[v][3];
	  _vertex_R2D[0]      = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2]);
	  _vertex_sR2D[0]     = TMath::Sqrt(derivate2_with_sigmaR2D(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2])  + derivate2_with_sigmaR2D(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2])   + 2*_vertices[v][7]*_vertices[v][7] *derivateR2D(_vertices[v][1],_vertices[v][1], _vertices[v][2])*derivateR2D(_vertices[v][2],_vertices[v][1], _vertices[v][2]) );
	  _vertex_R[0]        = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2] + _vertices[v][3]*_vertices[v][3]);
	  _vertex_sR[0]= TMath::Sqrt(derivate2_with_sigmaR(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     derivate2_with_sigmaR(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     derivate2_with_sigmaR(_vertices[v][3], _vertices[v][6],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][7]*_vertices[v][7]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][9]*_vertices[v][9]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][8]*_vertices[v][8]*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3])  );
	  _vertex_chi2[0]     = _vertices[v][11];
	  _vertex_normchi2[0] = _vertices[v][11]/_vertices[v][10];
	}
	if ((_vertices[v][0] == (_lIndex[index_l[0]] * 100 +  _lIndex[index_l[1]] )) ||(_vertices[v][0] == (_lIndex[index_l[0]] +  _lIndex[index_l[1]] *100) )) {
	  _vertex_X[1]        = _vertices[v][1];
	  _vertex_Y[1]        = _vertices[v][2];
	  _vertex_Z[1]        = _vertices[v][3];
	  _vertex_sX[1]        = _vertices[v][4];
	  _vertex_sY[1]        = _vertices[v][5];
	  _vertex_sZ[1]        = _vertices[v][6];
	  _vertex_R2D[1]      = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2]);
	  _vertex_sR2D[1]     = TMath::Sqrt(derivate2_with_sigmaR2D(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2])  + derivate2_with_sigmaR2D(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2])   + 2*_vertices[v][7]*_vertices[v][7] *derivateR2D(_vertices[v][1],_vertices[v][1], _vertices[v][2])*derivateR2D(_vertices[v][2],_vertices[v][1], _vertices[v][2]) );
	  _vertex_R[1]        = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2] + _vertices[v][3]*_vertices[v][3]);
	  _vertex_sR[1]= TMath::Sqrt(derivate2_with_sigmaR(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     derivate2_with_sigmaR(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     derivate2_with_sigmaR(_vertices[v][3], _vertices[v][6],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][7]*_vertices[v][7]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][9]*_vertices[v][9]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][8]*_vertices[v][8]*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3])  );
	  _vertex_chi2[1]     = _vertices[v][11];
	  _vertex_normchi2[1] = _vertices[v][11]/_vertices[v][10];
	}
	if ((_vertices[v][0] == (_lIndex[ind_new_leading] * 100 +  _lIndex[index_l[1]])) || (_vertices[v][0] == (_lIndex[ind_new_leading] +  _lIndex[index_l[1]]*100))) {
	  _vertex_X[2]        = _vertices[v][1];
	  _vertex_Y[2]        = _vertices[v][2];
	  _vertex_Z[2]        = _vertices[v][3];
	  _vertex_R2D[2]      = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2]);
	  _vertex_sR2D[2]     = TMath::Sqrt(derivate2_with_sigmaR2D(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2])  + derivate2_with_sigmaR2D(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2])   + 2*_vertices[v][7]*_vertices[v][7] *derivateR2D(_vertices[v][1],_vertices[v][1], _vertices[v][2])*derivateR2D(_vertices[v][2],_vertices[v][1], _vertices[v][2]) );
	  _vertex_R[2]        = TMath::Sqrt(_vertices[v][1]*_vertices[v][1]+ _vertices[v][2]*_vertices[v][2] + _vertices[v][3]*_vertices[v][3]);
	  _vertex_sR[2]= TMath::Sqrt(derivate2_with_sigmaR(_vertices[v][1], _vertices[v][4],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     derivate2_with_sigmaR(_vertices[v][2], _vertices[v][5],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     derivate2_with_sigmaR(_vertices[v][3], _vertices[v][6],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][7]*_vertices[v][7]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][9]*_vertices[v][9]*derivateR(_vertices[v][1],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3]) +
				     2*_vertices[v][8]*_vertices[v][8]*derivateR(_vertices[v][2],_vertices[v][1], _vertices[v][2], _vertices[v][3])*derivateR(_vertices[v][3],_vertices[v][1], _vertices[v][2], _vertices[v][3])  );
	  _vertex_chi2[2]     = _vertices[v][11];
	  _vertex_normchi2[2] = _vertices[v][11]/_vertices[v][10];
	}	           
      }// end loop vertices


      TVector3 primary_vertex[1];
      TVector3 secondary_vertex[1];
      TVector3 secondary_vertex_mc[2];
      primary_vertex[0].SetXYZ(_pvX,_pvY,_pvZ);
      secondary_vertex[0].SetXYZ(_vertex_X[1],_vertex_Y[1],_vertex_Z[1]);
      secondary_vertex_mc[0].SetXYZ(_lMatchVertexX[index_l[0]],_lMatchVertexY[index_l[0]], _lMatchVertexZ[index_l[0]]);
      secondary_vertex_mc[1].SetXYZ(_lMatchVertexX[index_l[1]],_lMatchVertexY[index_l[1]], _lMatchVertexZ[index_l[1]]);

      double delta_pv_sv=  (primary_vertex[0] - secondary_vertex[0]).Mag();
      double delta_pv_sv_mcl1=  (primary_vertex[0] - secondary_vertex_mc[0]).Mag();
      double delta_pv_sv_mcl2=  (primary_vertex[0] - secondary_vertex_mc[1]).Mag();
      double delta_sv_sv_mcl1=  (secondary_vertex[0] - secondary_vertex_mc[0]).Mag();
      double delta_sv_sv_mcl2=  (secondary_vertex[0] - secondary_vertex_mc[1]).Mag();


       double D2_delta_pv_sv= sqrt(  (primary_vertex[0].X()-secondary_vertex[0].X())*(primary_vertex[0].X()-secondary_vertex[0].X())   +    (primary_vertex[0].Y()-secondary_vertex[0].Y())*(primary_vertex[0].Y()-secondary_vertex[0].Y()) );

      double sigma_delta_pv_sv= TMath::Sqrt(error_ingridient(_pvX,_vertex_X[1],_pvXErr,  _vertex_sX[1]) + error_ingridient(_pvY,_vertex_Y[1],_pvYErr, _vertex_sY[1]) + error_ingridient(_pvZ,_vertex_Z[1],_pvZErr, _vertex_sX[1]) )/delta_pv_sv;

      double sigma_delta_pv_sv_mcl1= TMath::Sqrt(error_ingridient(_pvX,_lMatchVertexX[index_l[0]],_pvXErr, 0) + error_ingridient(_pvY,_lMatchVertexY[index_l[0]],_pvYErr, 0) + error_ingridient(_pvZ,_lMatchVertexZ[index_l[0]],_pvZErr, 0) )/delta_pv_sv_mcl1;
      double sigma_delta_pv_sv_mcl2= TMath::Sqrt(error_ingridient(_pvX,_lMatchVertexX[index_l[1]],_pvXErr, 0) + error_ingridient(_pvY,_lMatchVertexY[index_l[1]],_pvYErr, 0) + error_ingridient(_pvZ,_lMatchVertexZ[index_l[1]],_pvZErr, 0) )/delta_pv_sv_mcl2;

      double sigma_delta_sv_sv_mcl1= TMath::Sqrt(error_ingridient(_vertex_X[1],_lMatchVertexX[index_l[0]],_vertex_sX[1], 0) + error_ingridient(_vertex_Y[1],_lMatchVertexY[index_l[0]],_vertex_sY[1], 0) + error_ingridient(_vertex_Z[1],_lMatchVertexZ[index_l[0]],_vertex_sZ[1], 0) )/delta_sv_sv_mcl1;
      double sigma_delta_sv_sv_mcl2= TMath::Sqrt(error_ingridient(_vertex_X[1],_lMatchVertexX[index_l[1]],_vertex_sX[1], 0) + error_ingridient(_vertex_Y[1],_lMatchVertexY[index_l[1]],_vertex_sY[1], 0) + error_ingridient(_vertex_Z[1],_lMatchVertexZ[index_l[1]],_vertex_sZ[1], 0) )/delta_sv_sv_mcl2;


      if (_lProvenanceCompressed[index_l[0]] == 4){
	delta_pv_sv_mcl1 = -5;
	delta_sv_sv_mcl1 = -5;	 
      }
      if (_lProvenanceCompressed[index_l[1]] == 4){
	delta_pv_sv_mcl2 = -5;
	delta_sv_sv_mcl2 = -5;	 
      }

      double delta_mathced=0;
      if(_lProvenanceCompressed[index_l[0]] == 4 || _lProvenanceCompressed[index_l[1]] == 4)delta_mathced= -5;
      else delta_mathced= (secondary_vertex_mc[0]- secondary_vertex_mc[1]).Mag();

      double ratio_pv_sv_pv=delta_pv_sv/  primary_vertex[0].Mag();    
      double ratio_pv_sv=  delta_pv_sv/sigma_delta_pv_sv;
      double ratio_pv_sv_mcl1=  delta_pv_sv_mcl1/sigma_delta_pv_sv_mcl1;
      double ratio_pv_sv_mcl2=  delta_pv_sv_mcl2/sigma_delta_pv_sv_mcl2;
      double ratio_sv_sv_mcl1=  delta_sv_sv_mcl1/sigma_delta_sv_sv_mcl1;
      double ratio_sv_sv_mcl2=  delta_sv_sv_mcl2/sigma_delta_sv_sv_mcl2;
        
        
      double delta_position_l2=(secondary_vertex[0] - l2_position[0]).Mag();
      double delta_position_l3=(secondary_vertex[0] - l3_position[0]).Mag();


      double cos_stuff=0;
      TVector3 l2plusl3=  (lepton_reco[1] + lepton_reco[2]).Vect().Unit();
      TVector3 svMpv =secondary_vertex[0]- primary_vertex[0];
      double vtxR     = svMpv.Mag();
      double vtxRvtxPcosAlpha = svMpv.Dot(l2plusl3)/vtxR;
 
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   
      //M_3l
      sum_3l_rec.SetPtEtaPhiE(0,0,0,0);
      sum_3l_rec= (lepton_reco[0]+ lepton_reco[1]+lepton_reco[2] );
  

            
      //================== event classification ========================
      // ---------------- > OSSF or NO_OSSF
      ossf_no_ossf( kind, pair,lepton_reco[0], lepton_reco[1], lepton_reco[2], flavors_3l, charge_3l);
      if (kind[0]  == -1) continue;

      //---------------- > M_2l best Z candidate
      sum_2l_rec_pair.SetPtEtaPhiE(0,0,0,0);
      sum_2l_rec_pair= (pair[0]+pair[1] );
      bool ossf_event= false;
      if (kind[0] == 1) ossf_event = true;

      
      METvec.SetPtEtaPhiE(_met, 0, _metPhi,_met);    
      lepton_transv[0].SetPtEtaPhiE(lepton_reco[0].Pt(),0, lepton_reco[0].Phi(), lepton_reco[0].Pt());
      lepton_transv[1].SetPtEtaPhiE(lepton_reco[1].Pt(),0, lepton_reco[1].Phi(), lepton_reco[1].Pt());
      lepton_transv[2].SetPtEtaPhiE(lepton_reco[2].Pt(),0, lepton_reco[2].Phi(), lepton_reco[2].Pt());	

      double mT_met = 0;
      mT_met = (lepton_transv[0] +lepton_transv[1] +lepton_transv[2] + METvec).Mag();
      double mlll_met = 0;
      mlll_met=  (lepton_reco[0] +lepton_reco[1] +lepton_reco[2]).Mag() +_met;
  
    


      // ---------------- > CHANNELS
      class_os( event_clas,  flavors_3l, charge_3l);
      if (event_clas[0] == -1) continue;
      // 1* = eee
      // 2* = emm
      // 3* = eem
      // 4* = mmm
      isAll = true;

      // ---------------- > M_min
      check_mt=-1;    
      _mll_min = (lepton_reco[0]+lepton_reco[1]).M();
      check_mt= 2;    
      if ( (lepton_reco[0]+lepton_reco[2]).M() < _mll_min){
	_mll_min = (lepton_reco[0]+lepton_reco[2]).M();
	check_mt= 1;
      }
      if ( (lepton_reco[1]+lepton_reco[2]).M() < _mll_min) {
	_mll_min = (lepton_reco[1]+lepton_reco[2]).M();
	check_mt= 0;
      }   
      // 2 == ls
      // 1 == lt
      // 0 == st
      // ---------------- > M_min OS
      check_mt_os=-1;
      if ((charge_3l[0] != charge_3l[1] )) {
	_mll_min_os = (lepton_reco[0]+lepton_reco[1]).M();
	check_mt_os= 2;
      }
      if ((charge_3l[0] != charge_3l[2] ) && (lepton_reco[0]+lepton_reco[2]).M() < _mll_min_os){
	_mll_min_os = (lepton_reco[0]+lepton_reco[2]).M();
	check_mt_os= 1;
      }
      if ((charge_3l[1] != charge_3l[2] ) && (lepton_reco[1]+lepton_reco[2]).M() < _mll_min_os) {
	_mll_min_os = (lepton_reco[1]+lepton_reco[2]).M();
	check_mt_os= 0;
      }   
      check_deltaR=-1;
      delta_R_min = lepton_reco[0].DeltaR(lepton_reco[1]);
      check_deltaR=2;
      if (lepton_reco[0].DeltaR(lepton_reco[2]) < delta_R_min) {
	delta_R_min = lepton_reco[0].DeltaR(lepton_reco[2]);
	check_deltaR=1;
      }
      if (lepton_reco[1].DeltaR(lepton_reco[2]) < delta_R_min) {
	delta_R_min = lepton_reco[1].DeltaR(lepton_reco[2]);
	check_deltaR=0;
      }



      double min_delta_phi = 0;
      min_delta_phi = fabs(lepton_reco[0].DeltaPhi(lepton_reco[1]));
      if (fabs(lepton_reco[0].DeltaPhi(lepton_reco[2])) < min_delta_phi)  min_delta_phi = fabs(lepton_reco[0].DeltaPhi(lepton_reco[2]));
      double ration1_deltaphi = fabs(lepton_reco[1].DeltaPhi(lepton_reco[2]))/fabs(lepton_reco[0].DeltaPhi(lepton_reco[1]));
      double ration2_deltaphi = fabs(lepton_reco[0].DeltaPhi(lepton_reco[1]))/fabs(lepton_reco[2].DeltaPhi(lepton_reco[1]));

      double ration1_deltaphimin = fabs(lepton_reco[1].DeltaPhi(lepton_reco[2]))/min_delta_phi;
      double ration2_deltaphimin = min_delta_phi/fabs(lepton_reco[2].DeltaPhi(lepton_reco[1]));

      double ration1_deltaphiRmin = fabs(lepton_reco[1].DeltaR(lepton_reco[2]))/min_delta_phi;
      double ration2_deltaphiRmin = min_delta_phi/fabs(lepton_reco[2].DeltaR(lepton_reco[1]));


      //=============================================================
      
      
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      HISTOGRAMS     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      double pt_cone_leading=          lepton_reco[0].Pt() ;
      double pt_cone_sub_leading=      lepton_reco[1].Pt();
      double pt_cone_trailing=         lepton_reco[2].Pt();
  
      double trackIso_subleading =_trackIso[index_l[0]];
      double trackIso_trailing= _trackIso[index_l[1]];

      if (lepton_reco[1].DeltaR(lepton_reco[2]) < 0.3){
	_trackIso[index_l[1]] =  _trackIso[index_l[1]] - pt_cone_sub_leading;
	_trackIso[index_l[0]] =  _trackIso[index_l[0]] - pt_cone_trailing;
      }
            


      bool selection_0=false;
      bool selection_1=false;
      bool selection_2=false;
      bool selection_3=false;
      bool selection_4=false;
      bool selection_5=false;
      bool selection_final=false;
     
      if (!_isT[index_l[1]] ) continue;
      if (!_isT[index_l[0]] ) continue;

      if (charge_3l[2] != charge_3l[1])                                        selection_0 = true;
      if ( selection_0 && lepton_reco[1].DeltaR(lepton_reco[2]) < 1)           selection_1 = true;
      if ( selection_1 &&  bjet == 0 )                                         selection_2 = true;
      else {
	// If there are b-jets, do not skip the event, just give it a small weight, (1-SF)
	//  N.B.:  the "rejected" events Nr are scaled to Nr' = Nr * SF ==> DNr = Nr' - Nr = (SF-1)*Nr
	//         the "accepted" events Na are corrected to Na' = Na - DNr = Na - (SF-1)*Nr = Na + (1-SF)*Nr
	scal *= (1. - bwght);
	selection_2 = true;
      }
      if ( selection_2 &&  sum_3l_rec.M() > 45 && sum_3l_rec.M() < 85)         selection_3 = true;
      if ( selection_3 && min_delta_phi > 1)                                   selection_4 = true;
      if ( selection_4 &&  vtxRvtxPcosAlpha > 0.9)                             selection_5 = true;
      if ( selection_5 )                                                       selection_final = true;




      bool less2=false;
      bool more2_20=false;
      bool more20=false;
      
      bool less5=false;
      bool more5=false;

      if (D2_delta_pv_sv < 2)                         less2= true;
      if (D2_delta_pv_sv >= 2 && D2_delta_pv_sv < 10) more2_20= true;
      if (D2_delta_pv_sv >= 10 )                      more20= true;

      if ((lepton_reco[1]+lepton_reco[2]).M() < 5 )   less5= true;
      if ((lepton_reco[1]+lepton_reco[2]).M() > 5 )   more5= true;

       
      //======================= event selection =================== 
      if (!selection_5) continue;
      if (flavors_3l[0] !=1) continue;
      if (flavors_3l[1] ==0 && flavors_3l[2] ==0 ) continue;
      if ((lepton_reco[1]+lepton_reco[2]).M() > 50 ) continue;
      //if (delta_pv_sv < 0.5) continue;
      // if (ratio_pv_sv < 15 ) continue;

       int index_gen1=-1;
      int index_gen2=-1;
      int index_gen3=-1;
      TVector3 vector2;
      TVector3 vector3;

      for(unsigned l = 0; l < _gen_nL; ++l){
	if (_lMatchPt[ind_new_leading] == _gen_lPt[l]) index_gen1 = l;
	if (_lMatchPt[index_l[0]] == _gen_lPt[l]) index_gen2 = l;
	if (_lMatchPt[index_l[1]] == _gen_lPt[l]) index_gen3 = l;
      }

      if (index_gen2 !=-1) vector2.SetXYZ(_lMatchVertexX[index_l[0]],_lMatchVertexY[index_l[0]] ,_lMatchVertexZ[index_l[0]]);
      if (index_gen3 !=-1) vector3.SetXYZ(_lMatchVertexX[index_l[1]],_lMatchVertexY[index_l[1]] ,_lMatchVertexZ[index_l[1]]);

      tt_ttt_cj=false;
      tt_ttt_dr=false;
      tt_ttt_df=false;
      tt_ttt=false;

      if (_closestJetCsvV2[index_l[0]] ==  _closestJetCsvV2[index_l[1]]){
	tt_ttt_cj= true;
      }
      if (lepton_reco[1].DeltaR(lepton_reco[2]) < 0.4){
	tt_ttt_dr= true;
      }
      if (_closestJetCsvV2[index_l[0]] ==  _closestJetCsvV2[index_l[1]]){
	if (_lProvenanceCompressed[index_l[0]] !=0 && _lProvenanceCompressed[index_l[1]] !=0) {
	  tt_ttt_df= true;
	}	  
      }
      if (_closestJetCsvV2[index_l[0]] !=  _closestJetCsvV2[index_l[1]]){
	if (_lProvenance[index_l[0]] == 9 && _lProvenance[index_l[1]] == 9 && vector2.DeltaR(vector3) == 0 ) {
	  tt_ttt_df= true;
	}
	if (_lProvenanceCompressed[index_l[0]] == 1 && _lProvenanceCompressed[index_l[1]] == 1 ) {
	  if ((_lProvenance[index_l[0]] + _lProvenance[index_l[1]]) == 18 || (_lProvenance[index_l[0]] + _lProvenance[index_l[1]]) == 17){
	    if ( vector2.DeltaR(vector3) < 1.2 && _gen_lMomPdg[index_gen2] !=0 &&  _gen_lMomPdg[index_gen3] !=0){
	      tt_ttt_df= true;
	    }
	    if ( _gen_lMomPdg[index_gen2] ==0 ||  _gen_lMomPdg[index_gen3] ==0){
	      tt_ttt_df= true;
	    }
	  }
	}
      }
      if (!tt_ttt_df) tt_ttt=true;
      
          

      //======================= binning histogram!    ===================
      int bin_histo[nSR]={0,0,0,0,0,0,0,0,0,0,0,0};
      if (less2){
	if (less5){
	  if (ossf_event)         bin_histo[0]=1;
	  if (!ossf_event)        bin_histo[1]=1;
	}//less5
	if (more5){
	  if (ossf_event)         bin_histo[2]=1;
	  if (!ossf_event)        bin_histo[3]=1;
	}//more5
      }//less2
      
      if (more2_20){
	if (less5){
	  if (ossf_event)         bin_histo[4]=1;
	  if (!ossf_event)        bin_histo[5]=1;
	}//less5
	if (more5){
	  if (ossf_event)         bin_histo[6]=1;
	  if (!ossf_event)        bin_histo[7]=1;
	}//more5
      }//less2
      
      if (more20){
	if (less5){
	  if (ossf_event)         bin_histo[8]=1;
	  if (!ossf_event)        bin_histo[9]=1;
	}//less5
	if (more5){
	  if (ossf_event)         bin_histo[10]=1;
	  if (!ossf_event)        bin_histo[11]=1;
	}//more5
      }//less2
      //==========================================



      
      
            
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            
 
      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  FILLING  HISTOGRAMS  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      double values[nDist] ={static_cast<double>(0)};

      //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      unsigned fill = effsam;
  
      // ------------------- Histo kinematics
      for(int numero_histo = 0; numero_histo < nDist; ++numero_histo){
	for (int ii = 0; ii< nSR; ii++){
	  if (bin_histo[ii] == 1)  {
	    Histos[0][0][fill]->Fill(TMath::Min(static_cast<double>(ii+1), maxBinC[0]), scal);
	    // if (tt_ttt) Histos[0][1][fill]->Fill(TMath::Min(static_cast<double>(ii+1), maxBinC[0]), scal);
	    // if (tt_ttt_df) Histos[0][2][fill]->Fill(TMath::Min(static_cast<double>(ii+1), maxBinC[0]), scal);
	    // if (tt_ttt_cj) Histos[0][3][fill]->Fill(TMath::Min(static_cast<double>(ii+1), maxBinC[0]), scal);
	    // if (tt_ttt_dr) Histos[0][4][fill]->Fill(TMath::Min(static_cast<double>(ii+1), maxBinC[0]), scal);

	    if(runtheosyst) {
	      for(unsigned sidx=0; sidx<nTheoVars; ++sidx) {
		double wghtCorr = _lheWeight[theoSystVars[sidx]];
		systHistos[sidx][0][fill]->Fill(TMath::Min(static_cast<double>(ii+1), maxBinC[0]), scal*wghtCorr);
	      }
	    }
	  }
	}
      }//end histo
      
    
    }
  } 
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  /*
  //Split data and MC histograms for plotting and propagating uncertainties
  TH1D* dataYields[nDist][nCat];
  for(unsigned dist = 0; dist < nDist; ++dist){
    for(unsigned cat = 0; cat < nCat; ++cat){
      dataYields[dist][cat] = (TH1D*) Histos[dist][cat][63]->Clone();
    }
  }
    
  // cout<< "ok 1"<<endl;
    
    
  TH1D* bkgYields[nDist][nCat][nSamples_eff -62]; //change to nSamples_eff if sig is removed
  for(unsigned dist = 0; dist < nDist; ++dist){
    for(unsigned cat = 0; cat < nCat; ++cat){
      for(unsigned effsam1 = 63; effsam1 < nSamples_eff +1 ; ++effsam1){
	
	put_at_zero(*&Histos[dist][cat][effsam1]);

	if(runtheosyst && cat==0) {
	  for(unsigned sidx=0; sidx<nTheoVars; ++sidx) {
	    put_at_zero(*&systHistos[sidx][dist][effsam1]);
	  }
	}

	//	cout<< effsam1<<"   "<<nSamples_eff<<endl;
	bkgYields[dist][cat][effsam1 -63] = (TH1D*) Histos[dist][cat][effsam1]->Clone();

	if(effsam1 > 63 && effsam1 < 69){
	  
	  dataYields[dist][cat]->Add(bkgYields[dist][cat][effsam1 -63]);
	}
      }
    }
  }

  const std::string catNamespp[nCat]= {""};

  const std::string systNames[7]= {"stasignal","statDY","statttbar",	"statWJets",	"statmultiboson", "statXgamma","statTTTX"};
  const std::string systDist[7]= {"lnN","lnN", "lnN",	"lnN",	"lnN", 	"lnN", "lnN"};
  const std::string bgkNames[6]= {"DY",  	"ttbar",	"WJets",	"multiboson", 	"Xgamma",    	"TTTX"};

  const TString sigNames[62] = {"M2_v0024mu","M2_v0028mu","M2_v0031mu","M2_v0044mu","M2_v0054mu","M2_v0070mu","M2_v0083mu","M2_v0141mu","M2_v0173mu","M2_v01mu","M2_v0223mu","M3_v0024mu","M3_v0028mu","M3_v0031mu","M3_v0044mu","M3_v0054mu","M3_v0070mu","M3_v0083mu","M3_v0141mu","M3_v0173mu","M3_v01mu","M4_v0024mu","M4_v0028mu","M4_v0031mu","M4_v0044mu","M4_v0054mu","M4_v0070mu","M4_v0083mu","M4_v0141mu","M4_v01mu","M5_v0024mu","M5_v0028mu","M5_v0031mu","M5_v0044mu","M5_v0054mu","M5_v0070mu","M5_v0083mu","M5_v01mu","M8_v0024mu","M8_v0028mu","M8_v0031mu","M8_v0044mu","M8_v0054mu","M10_v0024mu","M10_v0028mu","M10_v0031mu","M10_v0044mu","M15_v0024mu","M15_v0028mu","M15_v0031mu","M15_v0044mu","M20_v1933mu", "M1_v0949mu","M1_v2123mu","M2_v0110mu","M2_v0248mu","M3_v0071mu","M4_v0029mu","M5_v0014mu","M6_v0020mu","M8_v00154mu","M10_v0007mu"};
  
  const std::string sigNamespp[62] = { "M2_v0024mu","M2_v0028mu","M2_v0031mu","M2_v0044mu","M2_v0054mu","M2_v0070mu","M2_v0083mu","M2_v0141mu","M2_v0173mu","M2_v0100mu","M2_v0223mu","M3_v0024mu","M3_v0028mu","M3_v0031mu","M3_v0044mu","M3_v0054mu","M3_v0070mu","M3_v0083mu","M3_v0141mu","M3_v0173mu","M3_v0100mu","M4_v0024mu","M4_v0028mu","M4_v0031mu","M4_v0044mu","M4_v0054mu","M4_v0070mu","M4_v0083mu","M4_v0141mu","M4_v0100mu","M5_v0024mu","M5_v0028mu","M5_v0031mu","M5_v0044mu","M5_v0054mu","M5_v0070mu","M5_v0083mu","M5_v0100mu","M8_v0024mu","M8_v0028mu","M8_v0031mu","M8_v0044mu","M8_v0054mu","M10_v0024mu","M10_v0028mu","M10_v0031mu","M10_v0044mu","M15_v0024mu","M15_v0028mu","M15_v0031mu","M15_v0044mu","M20_v1933mu", "M1_v0949mu","M1_v2123mu","M2_v0110mu","M2_v0248mu","M3_v0071mu","M4_v0029mu","M5_v0014mu","M6_v0020mu","M8_v00154mu","M10_v0007mu"};


  
  
  
  TH1D* signals[nSamples_signal];
  //Plot the yields as a function of the search region


  
  for(unsigned dist = 0; dist < nDist; ++dist){
    for(unsigned cat = 0; cat < nCat; ++cat){
      for (unsigned signal_sample = 0; signal_sample< 62; signal_sample++){
	signals[signal_sample] = (TH1D*) Histos[dist][cat][signal_sample+1]->Clone() ;     
      }
      plotDataVSMC(cat,dist,dataYields[dist][cat], bkgYields[dist][cat], eff_names, nSamples_eff -  nSamples_signal - 1, Histnames_ossf[dist] + "_" +  catNames[cat], catNames[cat], true, 2, true, signals,  sigNames , nSamples_signal, false);

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


	for (unsigned signal_sample = 0; signal_sample< 62; signal_sample++){	  
	  if(signals[signal_sample] -> GetBinContent(ii+1) !=0) systUnc[0][0]= (1+ (signals[signal_sample] -> GetBinError(ii+1)) / (signals[signal_sample] -> GetBinContent(ii+1)));
	  if(systUnc[0][0] >=2 )                                systUnc[0][0] = 1.99;	  
	  if(signals[signal_sample] -> GetBinContent(ii+1) ==0) systUnc[0][0] = 1;


	  
	  // printDataCard(   dataYields[dist][cat] ->GetBinContent(ii+1), 
	  // 		   signals[signal_sample]->GetBinContent(ii+1),
	  // 		   sigNamespp[signal_sample],
	  // 		   bgkield,
	  // 		   6,
	  // 		   bgkNames,
	  // 		   systUnc, 7, systNames,systDist,
	  // 		   sigNamespp[signal_sample]+"_bin"+std::to_string(ii+1)+".txt",
	  // 		   false, sigNamespp[signal_sample]+"_bin"+std::to_string(ii+1), ii+1);


	}// end loop signal


      }//loop bin
      
    }//end cat
  }//end histo

*/

  fout->cd();
  for(int i=0; i<nDist; ++i) {
    for(int effsam=0; effsam<nSamples_eff; ++effsam) {
      if(effsam==0) continue;
      for(int cat=0; cat<nCat; ++cat) {
	// Uncertainty from MC statistics
	if(systcat==11) {
	  for(int k=1; k<=nBins[i]; ++k) {
	    double newcont = Histos[i][cat][effsam]->GetBinContent(k);
	    if(systdir==0) newcont -= Histos[i][cat][effsam]->GetBinError(k);
	    else           newcont += Histos[i][cat][effsam]->GetBinError(k);
	    Histos[i][cat][effsam]->SetBinContent(k, newcont);
	  }
	}
	Histos[i][cat][effsam]->Write();
      }
      if(runtheosyst) {
	std::vector<double> errorByBin(nBins[i], 0.);
	std::vector<double> meanByBin(nBins[i], 0.);
	for(unsigned sidx=0; sidx<nTheoVars; ++sidx) {
	  systHistos[sidx][i][effsam]->Write();
	  //
	  // Theory uncertainties
	  for(int k=0; k<nBins[i]; ++k) {
	    // QCD uncertainties
	    if(systcat==1) {
	      //std::cout << Histos[i][0][effsam]->GetBinContent(k+1) << " - " << systHistos[sidx][i][effsam]->GetBinContent(k+1) << std::endl;
	      double deltabin = Histos[i][0][effsam]->GetBinContent(k+1);
	      deltabin = deltabin>0. ? std::abs(systHistos[sidx][i][effsam]->GetBinContent(k+1) - deltabin)/deltabin : 0.;
	      if(deltabin>errorByBin[k]) errorByBin[k] = deltabin;
	    }
	    // PDF uncertainties
	    else if(systcat==2) {
	      //std::cout << Histos[i][0][effsam]->GetBinContent(k+1) << " - " << systHistos[sidx][i][effsam]->GetBinContent(k+1) << std::endl;
	      double iadd = Histos[i][0][effsam]->GetBinContent(k+1);
	      iadd = iadd>0. ? systHistos[sidx][i][effsam]->GetBinContent(k+1)/iadd : 0.;
	      meanByBin[k] += iadd;
	      errorByBin[k] += iadd*iadd;
	    }
	  } // end loop on bins: for(unsigned k=0; k<nBins[i]; ++k)
	} // end loop on vars: for(unsigned sidx=0; sidx<nTheoVars; ++sidx)
	//
	// Print theory systs results
	std::cout << " --- " << (systcat==1 ? "QCD" : "PDF") << " uncertainties --- " << std::endl;
	for(int k=0; k<nBins[i]; ++k) {
	  std::cout << "   >> bin " << i+1 << ": ";
	  if(systcat==1) {
	    std::cout << errorByBin[k] << std::endl;
	  }
	  else if(systcat==2) {
	    // Var[x] = [1/(N-1)] * [Sum(xi^2) - (Sum(xi))^2/N]
	    errorByBin[k] = (errorByBin[k] - (meanByBin[k]*meanByBin[k]/nTheoVars))/(nTheoVars-1);
	    errorByBin[k] = sqrt(errorByBin[k]);
	    std::cout << errorByBin[k] << std::endl;
	  }
	} // end for(int k=0; k<nBins[i]; ++k)
      } // end if(runtheosyst)
    } // end for(int effsam=0; effsam<nSamples_eff; ++effsam)
  } // end for(int i=0; i<nDist; ++i)
  fout->Close();

}// end analisi

//void plotDataVSMC(int categoria,int istogramma, TH1D* data, TH1D** bkg, const TString* names, const unsigned nHist, const TString& file, const TString& file2,const bool ylog = false, const unsigned widthopt = 0, const bool plotsig = false, TH1D** signal = nullptr, const TString* signames =nullptr, const unsigned nSig = 0, const bool signorm=false);


//___________________________________________________________________
void Analysis_mc::printDataCard(const double obsYield, const double sigYield, const std::string& sigName, const double* bkgYield, const unsigned nBkg, const std::string* bkgNames, const std::vector<std::vector<double> >& systUnc, const unsigned nSyst, const std::string* systNames, const std::string* systDist, const std::string& cardName, const bool shapeCard, const std::string& shapeFileName,int number_bin){

  //stream for writing card
  std::ofstream card;

  //add .txt to name if no file extension is given
  card.open(cardName + ((cardName.find(".") == std::string::npos) ? ".txt" : "") ); //add .txt to name if no file extension is given

  //define number of channels, background sources and systematics
  card << "imax 1 number of channels \n";
  card << "jmax " << nBkg << " number of backgrounds \n";
  card << "kmax " << nSyst << " number of nuisance parameters (sources of systematical uncertainties) \n";
  card << "---------------------------------------------------------------------------------------- \n";

  //define the channels and the number of observed events
  card << "bin bin1 \n";
  card << "observation " << obsYield << "\n";

  //define all backgrounds and their yields
  card << "---------------------------------------------------------------------------------------- \n";
  if(shapeCard){
    card << "shapes * * " << shapeFileName + ".root  $PROCESS $PROCESS_$SYSTEMATIC\n";
    card << "---------------------------------------------------------------------------------------- \n";
  }
  card << "bin   ";
  for(unsigned proc = 0; proc < nBkg + 1; ++proc){
    card << "	" << "bin1";
  }
  card << "\n";
  card << "process";
  card << "	" << sigName;
  for(unsigned bkg = 0; bkg < nBkg; ++bkg){
    card << "	" << bkgNames[bkg];
  }
  card << "\n";
  card << "process";
  for(unsigned bkg = 0; bkg < nBkg + 1; ++bkg){
    card << "	" << bkg;
  }
  card << "\n";
  card <<	"rate";
  card << "	" << sigYield;
  for(unsigned bkg = 0; bkg < nBkg; ++bkg){
    if(bkgYield[bkg] <= 0) card << "	" << "0.00";
    else card << "	" << bkgYield[bkg];
  }
  card << "\n";
  card << "---------------------------------------------------------------------------------------- \n";

  //define sources of systematic uncertainty, what distibution they follow and how large their effect is
  for(unsigned syst = 0; syst < nSyst; ++syst){
    if (number_bin<10 )card << systNames[syst]+"0"+std::to_string(number_bin) << "	" << systDist[syst];
    if (number_bin>=10)card << systNames[syst]+std::to_string(number_bin) << "	" << systDist[syst];

    for(unsigned proc = 0; proc < nBkg + 1; ++proc){
      card << "	";
      if(systUnc[syst][proc] == 0) card << "-";
      else card << systUnc[syst][proc];
    }
    card << "\n";
  }
  card.close();		
}






//___________________________________________________________________
double Analysis_mc::pu_weight ( TH1D *histo, double numberInteractions){


   
  //double momentum = part.Pt() * maximum( 1, iso - 0.1);
  double nI = numberInteractions;
    
  double factor=0;
  double factore=0;

  factore = histo->GetBinContent(histo->FindBin(nI));
      

  factor = factore;
  return factor;
 
}
//___________________________________________________________________

void Analysis_mc::norm ( TH2D *histo){

  
  if (histo->Integral() != 0 ) {
    const  double norma= 1/histo->Integral();
    histo-> Scale(norma*100);
  }
  else {
    histo -> SetXTitle("zero integral!!!!!");
    histo -> SetYTitle("zero integral!!!!!");

  }
}
//_______________________________________________________ constructor_____
void Analysis_mc::put_at_zero(TH1D *histo){
  for (int i =0; i < histo-> GetNbinsX(); i++){
    if (histo->GetBinContent( i+1)  < 0) histo-> SetBinContent(i+1, 0.0); // what about the error???
  }
}

//_______________________________________________________ constructor_____
void Analysis_mc::put_at_zero2d(TH2D *histo){
  for (int i =0; i < histo-> GetNbinsX(); i++){
    for (int j =0; j < histo-> GetNbinsY(); j++){

      if (histo->GetBinContent( i+1, j+1)  < 0) histo-> SetBinContent(i+1,j+1, 0.0);
    }
  }
}


//==================================================================

void Analysis_mc::find_leptons(int selezione,  unsigned displacedC, TLorentzVector lepton_tobeselected[20], int index_displaced[20], int index_s[2]){
  // ordine in pT
  int index_1=-5;
  int index_2= -5;
  if (selezione == 0 ){
    index_1 = index_displaced[0];
    index_2 = index_displaced[1];	
  }

  Double_t          delta_R_min=-1;
  Double_t          _mll_min=50000;
 
  //min mass
  if (selezione == 1 ){
    for (unsigned h = 0; h <displacedC; h ++ ){
      for (unsigned g = 0; g < displacedC; g ++){
	if (h != g) {
	  if (h == 0 && g ==1) {
	    _mll_min = (lepton_tobeselected[h]+lepton_tobeselected[g]).M();
	    index_1 = index_displaced[h];
	    index_2 = index_displaced[g];
	    if (lepton_tobeselected[h].Pt() < lepton_tobeselected[g].Pt()) {
	      index_1 = index_displaced[g];
	      index_2 = index_displaced[h];
	    }

	  }
	  if((lepton_tobeselected[h]+lepton_tobeselected[g]).M() < _mll_min) {
	    _mll_min = (lepton_tobeselected[h]+lepton_tobeselected[g]).M();
	    index_1 = index_displaced[h];
	    index_2 = index_displaced[g];
	    if (lepton_tobeselected[h].Pt() < lepton_tobeselected[g].Pt()) {
	      index_1 = index_displaced[g];
	      index_2 = index_displaced[h];
	    }

	  }
	}
      }// loop1
    }//loop 2
  }//selezione1

  //min delta
  if (selezione == 2 ){
    for (unsigned h = 0; h <displacedC; h ++ ){
      for (unsigned g = 0; g < displacedC; g ++){
	if (h != g) {
	  if (h == 0 && g ==1) {
	    delta_R_min=  lepton_tobeselected[h].DeltaR(lepton_tobeselected[g]);
	    index_1 = index_displaced[h];
	    index_2 = index_displaced[g];

	    if (lepton_tobeselected[h].Pt() < lepton_tobeselected[g].Pt()) {
	      index_1 = index_displaced[g];
	      index_2 = index_displaced[h];
	    }

	  }
	  if(lepton_tobeselected[h].DeltaR(lepton_tobeselected[g]) < delta_R_min) {
	    delta_R_min=  lepton_tobeselected[h].DeltaR(lepton_tobeselected[g]);
	    index_1 = index_displaced[h];
	    index_2 = index_displaced[g];

	    if (lepton_tobeselected[h].Pt() < lepton_tobeselected[g].Pt()) {
	      index_1 = index_displaced[g];
	      index_2 = index_displaced[h];
	    }
	  }
	}	
      }// loop1
    }//loop 2
  }
  index_s [0] = index_1;
  index_s[1]  = index_2; 

}//end funciton




//==================================================================
double Analysis_mc::maximum(double a, double b){
  double massimo=0;
  massimo = a;
  if (massimo < b) massimo = b;
  return massimo+1;
}
//___________________________________________________________________
double Analysis_mc::error_ingridient(double a, double b, double ea, double eb){
  double result =0;
  result = (a-b)*(a-b) *(ea*ea + eb*eb);
  return result;
}
//___________________________________________________________________
double Analysis_mc::derivateR(double a, double x, double y, double z){
  double result =0;
  result = a /((x*x + y*y + z*z)*(x*x + y*y + z*z));
  return result;
}
//___________________________________________________________________
double Analysis_mc::derivate2_with_sigmaR(double a, double sa, double x, double y, double z){
  double result =0;
  result = (a*a)*(sa*sa) /(x*x + y*y + z*z);
  return result;
}
//___________________________________________________________________
double Analysis_mc::derivateR2D(double a, double x, double y){
  double result =0;
  result = a /((x*x + y*y )*(x*x + y*y ));
  return result;
}
//___________________________________________________________________
double Analysis_mc::derivate2_with_sigmaR2D(double a, double sa, double x, double y){
  double result =0;
  result = (a*a)*(sa*sa) /(x*x + y*y );
  return result;
}


//___________________________________________________________________

void Analysis_mc::from_TGraph_to_TH1D (TGraphAsymmErrors *graph, TH1D *histo, int number_point){
    
  const int numero= number_point;
    
  double x_graph[numero];
  double y_graph[numero];
  for (int i =0; i <number_point; i ++){
    x_graph[i]=0;
    y_graph[i]=0;
  }
  for (int i =0; i <number_point; i ++){
    graph -> GetPoint(i, x_graph[i], y_graph[i]);
    histo->SetBinContent (i+1, x_graph[i],  y_graph[i]);
        
    //cout<<i<<") "<<y_graph[i]<<"  "<< histo->GetBinContent (i+1)<<endl;
        
  }
}
/*
//==================================================================
double Analysis_mc::fakeWeight(const unsigned ind, const int flavors, const double conePt, const double eta, const bool tight, TH2D* frMap, const unsigned lCount){
unsigned nFO = 0;
double fr[4];
for(unsigned l = 0; l < lCount; ++l){
if(!tight[ind[l]]){
fr[nFO] = frMap[flavors[ind[l]]]->GetBinContent(frMap[flavors[ind[l]]]->FindBin(TMath::Min(conePt[l], 99.), fabs(eta[ind[l]])));
++nFO;
}
}
double weight = 1;
for(unsigned f = 0; f < nFO; ++f){
weight *= fr[f]/(1-fr[f]);
}
if(nFO == 2) weight*= -1;
return weight;
}
*/






//==================================================================
double Analysis_mc::FR_factor(TGraphAsymmErrors *fakeRate_mu[3],
			      TGraphAsymmErrors *fakeRate_e[3],
			      double eta,
			      double flavors,
			      double lptcone
			      ){
    
    
  eta = fabs(eta);
    
    
  TH1D *fakeRate_mu_histo[3];
  TH1D *fakeRate_e_histo[3];
  Double_t newBins_mu1[7] = {5,10, 15, 25, 35, 50, 70};
  Double_t newBins_e1[6] = {10, 15, 25, 35, 50, 70};
  fakeRate_mu_histo[0]= new TH1D("fake_rate_mu_histo_eta1","",6,newBins_mu1);
  fakeRate_mu_histo[1]= new TH1D("fake_rate_mu_histo_eta2","",6,newBins_mu1);
  fakeRate_mu_histo[2]= new TH1D("fake_rate_mu_histo_eta3","",6,newBins_mu1);
  fakeRate_e_histo[0]= new TH1D("fake_rate_e_histo_eta1","",5,newBins_e1);
  fakeRate_e_histo[1]= new TH1D("fake_rate_e_histo_eta2","",5,newBins_e1);
  fakeRate_e_histo[2]= new TH1D("fake_rate_e_histo_eta3","",5,newBins_e1);
    
  for (int i=0; i< 3; i++){
    from_TGraph_to_TH1D(*&fakeRate_mu[i],*&fakeRate_mu_histo[i],6);
    from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],5);
  }
    
  //double momentum = part.Pt() * maximum( 1, iso - 0.1);
  double momentum = lptcone;
    
  double factor=0;
  double factore=0;
  if (flavors == 0){
    if (momentum < 49){
      if (eta < 0.8){
	factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(momentum));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(momentum));
      }//eta1
      else {
	factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(momentum));
      }//eta1
    }// <70
    else {
      if (eta < 0.8){
	factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(45));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(45));
      }//eta1
      else {
	factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(45));
      }//eta1
    }
  }//e
    
  if (flavors == 1){
    if (momentum < 49){
      if (eta < 0.8){
	factore = fakeRate_mu_histo[0]->GetBinContent(fakeRate_mu_histo[0]->FindBin(momentum));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_mu_histo[1]->GetBinContent(fakeRate_mu_histo[1]->FindBin(momentum));
      }//eta1
      else {
	factore = fakeRate_mu_histo[2]->GetBinContent(fakeRate_mu_histo[2]->FindBin(momentum));
      }//eta1
    }// <70
    else {
      if (eta < 0.8){
	factore = fakeRate_mu_histo[0]->GetBinContent(fakeRate_mu_histo[0]->FindBin(45));
      }//eta1
      else if (eta > 0.8 && eta<1.479){
	factore = fakeRate_mu_histo[1]->GetBinContent(fakeRate_mu_histo[1]->FindBin(45));
      }//eta1
      else {
	factore = fakeRate_mu_histo[2]->GetBinContent(fakeRate_mu_histo[2]->FindBin(45));
      }//eta1
    }
  }//e
    
  delete fakeRate_mu_histo[0];
  delete fakeRate_mu_histo[1];
  delete fakeRate_mu_histo[2];
  delete fakeRate_e_histo[0];
  delete fakeRate_e_histo[1];
  delete fakeRate_e_histo[2];
    
    
    
  return factore;
    
}





//___________________________________________________________________
void Analysis_mc::class_os(int event_clas[1], int  flavors_3l[3], int  charge_3l[3]){
    
  int ch_lepton1=charge_3l[0];
  int ch_lepton2=charge_3l[1];
  int ch_lepton3=charge_3l[2];
  int fl_lepton1=flavors_3l[0];
  int fl_lepton2=flavors_3l[1];
  int fl_lepton3=flavors_3l[2];
    
    
  // 1* = eee
  // 2* = emm
  // 3* = eem
  // 4* = mmm
    
    
    
    
    
  event_clas[0]=-1;
    
    
  if (ch_lepton1 == ch_lepton2 && ch_lepton1 == ch_lepton3 && ch_lepton3 == ch_lepton2)   event_clas[0]=-1;
    
    
  if (fl_lepton2 == 0 ||  fl_lepton3  == 0 ||  fl_lepton1== 0) {
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 0 ) event_clas[0] = 10; //e e e
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 1 ) event_clas[0] = 30; //e e mu
    if ((fl_lepton2 + fl_lepton3 + fl_lepton1) == 2 ) {
      if (fl_lepton2 == 1 ||  fl_lepton3  == 1 ||  fl_lepton1== 1) event_clas[0]=20; // e mu mu
    }
  }// at least an electron
  else if (fl_lepton2 == 1 &&  fl_lepton3  == 1 &&  fl_lepton1 == 1) {
    event_clas[0] = 40;
  }
  else {
    event_clas[0] =-1;
  }
    
  if (event_clas[0]  == 30){
    if ((fl_lepton1 == 0 && fl_lepton2 == 0 && fl_lepton3 == 1) ) { // e e mu
      if ((ch_lepton1 + ch_lepton2) == 0) event_clas[0] = 3;
    }
    if (fl_lepton1 == 0 && fl_lepton2 == 1 && fl_lepton3 == 0) { // e mu e
      if ((ch_lepton1 + ch_lepton3) == 0) event_clas[0] = 3;
    }
    if (fl_lepton1 == 1 && fl_lepton2 == 0 && fl_lepton3 == 0 ) { // mu e e
      if ((ch_lepton2 + ch_lepton3) == 0) event_clas[0] = 3;
    }
  }
    
    
  if (event_clas[0]  == 20){
    if ((fl_lepton1 == 1 && fl_lepton2 == 1 && fl_lepton3 == 0) ) { // e e mu
      if ((ch_lepton1 + ch_lepton2) == 0) event_clas[0] = 2;
    }
    if (fl_lepton1 == 1 && fl_lepton2 == 0 && fl_lepton3 == 1 ) { // e mu e
      if ((ch_lepton1 + ch_lepton3) == 0) event_clas[0] = 2;
    }
    if (fl_lepton1 == 0 && fl_lepton2 == 1 && fl_lepton3 == 1 ) { // mu e e
      if  ((ch_lepton2 + ch_lepton3) == 0)event_clas[0] = 2;
    }
  }
    
    
  if (event_clas[0]  == 10 ) {
    if ((ch_lepton1 + ch_lepton2 + ch_lepton3) == 1 || (ch_lepton1 + ch_lepton2 + ch_lepton3) == -1)  event_clas[0] =1;
  }
    
  if (event_clas[0]  == 40 ) {
    if ((ch_lepton1 + ch_lepton2 + ch_lepton3) == 1 || (ch_lepton1 + ch_lepton2 + ch_lepton3) == -1)  event_clas[0] =4;
  }
    
    
}



//___________________________________________________________________
void Analysis_mc::ossf_no_ossf(int kind[1],TLorentzVector pair[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int  flavors_3l[3], int  charge_3l[3]){
    
  int ch_lepton1=charge_3l[0];
  int ch_lepton2=charge_3l[1];
  int ch_lepton3=charge_3l[2];
  int fl_lepton1=flavors_3l[0];
  int fl_lepton2=flavors_3l[1];
  int fl_lepton3=flavors_3l[2];
    
    
  kind[0] = -1;
    
  if (ch_lepton1 == ch_lepton2 && ch_lepton1 == ch_lepton3 && ch_lepton3 == ch_lepton2)   kind[0] = -1;
    
    
  // OSSF
  if (     ((ch_lepton1 != ch_lepton2)    && (fl_lepton1 == fl_lepton2))  || ((ch_lepton1 != ch_lepton3)   && (fl_lepton1 == fl_lepton3)) || ((ch_lepton2 != ch_lepton3)  && (fl_lepton3 == fl_lepton2)) ){ // ossf
    //cout<<"in function where kind is 1: "<<kind[0]<<endl;
        
        
    kind[0] = 1;
    double i_m[3]={33333,33333,33333};
    double mass_inv=0;
    int index_inv=100;
    double min_mass=999;
    if ((ch_lepton1 != ch_lepton2)  && (fl_lepton1 == fl_lepton2)) i_m[0]= TMath:: Abs((leep1 + leep2).Mag() - 91.1876);
    if ((ch_lepton1 != ch_lepton3)  && (fl_lepton1 == fl_lepton3)) i_m[1]= TMath:: Abs((leep1 + leep3).Mag() - 91.1876);
    if ((ch_lepton2 != ch_lepton3)  && (fl_lepton3 == fl_lepton2)) i_m[2]= TMath:: Abs((leep2 + leep3).Mag() - 91.1876);
    for (int i =0; i < 3; i++){
      if (i_m[i] == 33333) continue;
      mass_inv = i_m[i];
      if (min_mass > mass_inv ){
	min_mass = mass_inv;
	index_inv = i;
      }
    }
    if (index_inv == 0) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[2].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
            
    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
            
    }
  }// end ossf
  // No_OSSF
  else if (   ((ch_lepton1 + ch_lepton2) == 0  )  || ((ch_lepton1 + ch_lepton3) == 0   ) || ((ch_lepton3 + ch_lepton2) == 0   )   ){
    //cout<<"in function where kind is 0: "<<kind[0]<<endl;
    kind[0] = 0;
    double i_m[3]={33333,33333,33333};
    double mass_inv=0;
    int index_inv=100;
    double min_mass=999;
    if ((ch_lepton1 != ch_lepton2) ) i_m[0]= TMath:: Abs((leep1 + leep2).Mag() - 91.1876);
    if ((ch_lepton1 != ch_lepton3)  ) i_m[1]= TMath:: Abs((leep1 + leep3).Mag() - 91.1876);
    if ((ch_lepton2 != ch_lepton3) ) i_m[2]= TMath:: Abs((leep2 + leep3).Mag() - 91.1876);
    for (int i =0; i < 3; i++){
      if (i_m[i] == 33333) continue;
      mass_inv = i_m[i];
      if (min_mass > mass_inv ){
	min_mass = mass_inv;
	index_inv = i;
      }
    }
    if (index_inv == 0) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[2].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
            
    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      pair[2].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
            
    }
        
        
  }//end no-ossf
    
  /*
    cout<< ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>   "<<kind<<endl;
    cout<<"1: "<<ch_lepton1<<"  "<<fl_lepton1<< "  "<< leep1.Pt()<< endl;
    cout<<"2: "<<ch_lepton2<<"  "<<fl_lepton2<< "  "<< leep2.Pt()<< endl;
    cout<<"3: "<<ch_lepton3<<"  "<<fl_lepton3<< "  "<< leep3.Pt()<< endl;
    cout<<"pair 1: "<<pair1.Pt()<< endl;
    cout<<"pair 2: "<<pair2.Pt()<< endl;
  */
  //cout<<"in function "<<kind[0]<<endl;
}





//___________________________________________________________________
void Analysis_mc::fr_selection(int number, TLorentzVector lepton_fake_order[3],TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int index_leptons[3],  int flavor_leptons[3], int origin_leptons[3],int index_3l[3],  int flavor_3l[3], int origin_3l[3]){
    
  lepton_fake_order[0].SetPtEtaPhiE(0,0,0,0);
  lepton_fake_order[1].SetPtEtaPhiE(0,0,0,0);
  lepton_fake_order[2].SetPtEtaPhiE(0,0,0,0);
  for(int i =0; i< 3; i++){
    index_leptons[i]=  -5;
    flavor_leptons[i]= -5;
    origin_leptons[i]= -5;
  }
    
  if (number == 3) {
    lepton_fake_order[0].SetPtEtaPhiE(leep1.Pt(),leep1.Eta(),leep1.Phi(),leep1.E());
    lepton_fake_order[1].SetPtEtaPhiE(leep2.Pt(),leep2.Eta(),leep2.Phi(),leep2.E());
    lepton_fake_order[2].SetPtEtaPhiE(leep3.Pt(),leep3.Eta(),leep3.Phi(),leep3.E());
    index_leptons[0]=index_3l[0];
    index_leptons[1]=index_3l[1];
    index_leptons[2]=index_3l[2];
    flavor_leptons[0]=flavor_3l[0];
    flavor_leptons[1]=flavor_3l[1];
    flavor_leptons[2]=flavor_3l[2];
    origin_leptons[0]=origin_3l[0];
    origin_leptons[1]=origin_3l[1];
    origin_leptons[2]=origin_3l[2];
  }
  if (number == 0) {
    lepton_fake_order[0].SetPtEtaPhiE(leep1.Pt(),leep1.Eta(),leep1.Phi(),leep1.E());
    lepton_fake_order[1].SetPtEtaPhiE(leep2.Pt(),leep2.Eta(),leep2.Phi(),leep2.E());
    lepton_fake_order[2].SetPtEtaPhiE(leep3.Pt(),leep3.Eta(),leep3.Phi(),leep3.E());
    index_leptons[0]=index_3l[0];
    index_leptons[1]=index_3l[1];
    index_leptons[2]=index_3l[2];
    flavor_leptons[0]=flavor_3l[0];
    flavor_leptons[1]=flavor_3l[1];
    flavor_leptons[2]=flavor_3l[2];
    origin_leptons[0]=origin_3l[0];
    origin_leptons[1]=origin_3l[1];
    origin_leptons[2]=origin_3l[2];
  }
    
    
    
}//end fr
