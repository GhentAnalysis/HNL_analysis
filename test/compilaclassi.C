#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#endif

void compilaclassi(TString myopt="fast"){
  std::string opt = "kg";
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  //gSystem->CompileMacro("tdrstyle.cc"                  , opt.c_str());
  gSystem->CompileMacro("plotCode_new.cc"              , opt.c_str());
  gSystem->CompileMacro("BTagCalibrationStandalone.cpp", opt.c_str());
  gSystem->CompileMacro("Analysis_mc.cc"               , opt.c_str());
  gSystem->CompileMacro("analisi.C"                    , opt.c_str());
  
}
