#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#endif

void compilaclassi(TString myopt="fast"){
  std::string opt = "kg";
  if(myopt.Contains("force")){
    opt = "kfg";
  }
  gSystem->CompileMacro("tdrstyle.cc"                  , opt.c_str());
  gSystem->CompileMacro("src/plotCode_new.cc"              , opt.c_str());
  gSystem->CompileMacro("src/Sample.cc"              , opt.c_str());
  gSystem->CompileMacro("src/objectSelection.cc"              , opt.c_str());
  gSystem->CompileMacro("src/kinematicTools.cc"              , opt.c_str());
  gSystem->CompileMacro("src/eventSelection.cc"              , opt.c_str());
  gSystem->CompileMacro("src/stringTools.cc"              , opt.c_str());
  gSystem->CompileMacro("src/Analysis_mc.cc"               , opt.c_str());



  gSystem->CompileMacro("/bTagging/BTagCalibrationStandalone.cpp", opt.c_str());
  gSystem->CompileMacro("analisi.C"                    , opt.c_str());
  
}
