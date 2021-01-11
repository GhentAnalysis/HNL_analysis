
#include "../interface/Analysis_mc.h"
#include "../interface/kinematicTools.h"


//______________________________________________ele base
bool Analysis_mc::eleIsClean(const unsigned electronIndex) const{
  if( is2016() ){
    return eleIsClean2016(electronIndex);
  }
  else if ( is2017() ){
    return eleIsClean2017(electronIndex);
  }
  else{
    return eleIsClean2018(electronIndex);
  }     
}
//______________________________________________muon loose used for ele cleaning base 
bool Analysis_mc::lepIsLoose(const unsigned leptonIndex) const{
  if( is2016() ){
    return lepIsLoose2016(leptonIndex);
  }
  else if ( is2017() ){
    return lepIsLoose2017(leptonIndex);
  }
  else{
    return lepIsLoose2018(leptonIndex);
  }  
}

//______________________________________________ele loose ID
bool Analysis_mc::eleIsLoose(const unsigned leptonIndex) const{
  if( is2016() ){
    return eleIsLoose2016(leptonIndex);
  }
  else if ( is2017() ){
    return eleIsLoose2017(leptonIndex);
  }
  else{
    return eleIsLoose2018(leptonIndex);
  } 

}



//______________________________________________remove electrons overlapping with muons
bool Analysis_mc::eleIsCleanBase(const unsigned electronIndex, bool (Analysis_mc::*looseMuon)(const unsigned) const) const{
  //make sure this lepton is an electron
  if( !isElectron(electronIndex) ){
    std::cerr << "Error: trying to clean non-electron object from muon overlap." << std::endl;
    return 999;
  }
  //check separation with every muon
  for(unsigned m = 0; m < _nMu; ++m){
    if( !isMu(m) ){
      std::cerr << "Error trying to clean electron from non-muon object" << std::endl;
    }
    if( ( this->*looseMuon )(m) ){
      if( kinematics::deltaR(_lPhi[m], _lEta[m], _lPhi[electronIndex], _lEta[electronIndex]) < 0.05 ){
	return false;
      }
    }
  }
  return true;
}
//______________________________________________ele clean base 16
bool Analysis_mc::eleIsClean2016(const unsigned electronIndex) const{
  return eleIsCleanBase(electronIndex, &Analysis_mc::lepIsLoose2016);
}
//______________________________________________ele clean base 17
bool Analysis_mc::eleIsClean2017(const unsigned electronIndex) const{
  return eleIsCleanBase(electronIndex, &Analysis_mc::lepIsLoose2017);
}
//______________________________________________ele clean base 18
bool Analysis_mc::eleIsClean2018(const unsigned electronIndex) const{
  return eleIsCleanBase(electronIndex, &Analysis_mc::lepIsLoose2018);
}

//______________________________________________loose POG 16
bool Analysis_mc::lepIsLoose2016(const unsigned leptonIndex) const{
  return _lPOGLoose[leptonIndex];
}
//______________________________________________loose POG 17
bool Analysis_mc::lepIsLoose2017(const unsigned leptonIndex) const{
  return _lPOGLoose[leptonIndex];
}
//______________________________________________loose POG 18
bool Analysis_mc::lepIsLoose2018(const unsigned leptonIndex) const{
  return _lPOGLoose[leptonIndex];
}
//______________________________________________ele loose ID
bool Analysis_mc::eleIsLoose2016(const unsigned leptonIndex) const{
  if( !isElectron(leptonIndex)) return false;
  if (!(_lEleIsEB[leptonIndex] || _lEleIsEE[leptonIndex])) return false;
	
  double elehadronicOverEm = 0.;
  double rho = 0.;
  if (std::abs (_lEta[leptonIndex]) < 1.00) rho	= _puCorr[leptonIndex] / 0.1440;
  if (std::abs (_lEta[leptonIndex]) >= 1.00 && std::abs (_lEta[leptonIndex]) < 1.479) rho= _puCorr[leptonIndex] / 0.1562;
  if (std::abs (_lEta[leptonIndex]) >= 1.479 && std::abs (_lEta[leptonIndex]) < 2.00) rho= _puCorr[leptonIndex] / 0.1032;
  if (std::abs (_lEta[leptonIndex]) >= 2.00 && std::abs (_lEta[leptonIndex]) < 2.20) rho= _puCorr[leptonIndex] / 0.0859;
  if (std::abs (_lEta[leptonIndex]) >= 2.20 && std::abs (_lEta[leptonIndex]) < 2.30) rho= _puCorr[leptonIndex] / 0.1116;
  if (std::abs (_lEta[leptonIndex]) >= 2.30 && std::abs (_lEta[leptonIndex]) < 2.40) rho= _puCorr[leptonIndex] / 0.1321;
  if (std::abs (_lEta[leptonIndex]) >= 2.40 && std::abs (_lEta[leptonIndex]) < 2.50) rho= _puCorr[leptonIndex] / 0.1654;
  double _x = _lEleIsEB[leptonIndex] ? 1.12       : 0.5 ;
  double _y = _lEleIsEB[leptonIndex] ? 0.0368       : 0.201 ;
  elehadronicOverEm = 	_lElehadronicOverEm[leptonIndex] - (_x + _y*rho)/_lEleEcalEnergy[leptonIndex];
	
	
  if(_lElefull5x5SigmaIetaIeta[leptonIndex]                  >= (_lEleIsEB[leptonIndex] ? 0.11       : 0.0314 ))       return false;
  if(_lEleDEtaInSeed  [leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.00477    : 0.00868))       return false;
  if(_lEleDeltaPhiSuperClusterTrackAtVtx [leptonIndex]       >= (_lEleIsEB[leptonIndex] ? 0.222      : 0.213  ))       return false;
  if(elehadronicOverEm                                       >= (_lEleIsEB[leptonIndex] ? 0.298      : 0.101  ))       return false;
  if(_lEleInvMinusPInv[leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.241      : 0.14   ))       return false;
  return true;

}
//______________________________________________ele loose ID
bool Analysis_mc::eleIsLoose2017(const unsigned leptonIndex) const{
  if( !isElectron(leptonIndex)) return false;
  if (!(_lEleIsEB[leptonIndex] || _lEleIsEE[leptonIndex])) return false;
	
  double elehadronicOverEm = 0.;
  double rho = 0.;
  if (std::abs (_lEta[leptonIndex]) < 1.00) rho	= _puCorr[leptonIndex] / 0.1440;
  if (std::abs (_lEta[leptonIndex]) >= 1.00 && std::abs (_lEta[leptonIndex]) < 1.479) rho= _puCorr[leptonIndex] / 0.1562;
  if (std::abs (_lEta[leptonIndex]) >= 1.479 && std::abs (_lEta[leptonIndex]) < 2.00) rho= _puCorr[leptonIndex] / 0.1032;
  if (std::abs (_lEta[leptonIndex]) >= 2.00 && std::abs (_lEta[leptonIndex]) < 2.20) rho= _puCorr[leptonIndex] / 0.0859;
  if (std::abs (_lEta[leptonIndex]) >= 2.20 && std::abs (_lEta[leptonIndex]) < 2.30) rho= _puCorr[leptonIndex] / 0.1116;
  if (std::abs (_lEta[leptonIndex]) >= 2.30 && std::abs (_lEta[leptonIndex]) < 2.40) rho= _puCorr[leptonIndex] / 0.1321;
  if (std::abs (_lEta[leptonIndex]) >= 2.40 && std::abs (_lEta[leptonIndex]) < 2.50) rho= _puCorr[leptonIndex] / 0.1654;
  double _x = _lEleIsEB[leptonIndex] ? 1.12       : 0.5 ;
  double _y = _lEleIsEB[leptonIndex] ? 0.0368       : 0.201 ;
  elehadronicOverEm = 	_lElehadronicOverEm[leptonIndex] - (_x + _y*rho)/_lEleEcalEnergy[leptonIndex];
	
  if(_lElefull5x5SigmaIetaIeta[leptonIndex]                  >= (_lEleIsEB[leptonIndex] ? 0.11       : 0.0314 ))       return false;
  if(_lEleDEtaInSeed  [leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.00477    : 0.00868))       return false;
  if(_lEleDeltaPhiSuperClusterTrackAtVtx [leptonIndex]       >= (_lEleIsEB[leptonIndex] ? 0.222      : 0.213  ))       return false;
  if(elehadronicOverEm                                       >= (_lEleIsEB[leptonIndex] ? 0.298      : 0.101  ))       return false;
  if(_lEleInvMinusPInv[leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.241      : 0.14   ))       return false;
  return true;

}
//______________________________________________ele loose ID
bool Analysis_mc::eleIsLoose2018(const unsigned leptonIndex) const{
 if( !isElectron(leptonIndex)) return false;
  if (!(_lEleIsEB[leptonIndex] || _lEleIsEE[leptonIndex])) return false;
	
  double elehadronicOverEm = 0.;
  double rho = 0.;
  if (std::abs (_lEta[leptonIndex]) < 1.00) rho	= _puCorr[leptonIndex] / 0.1440;
  if (std::abs (_lEta[leptonIndex]) >= 1.00 && std::abs (_lEta[leptonIndex]) < 1.479) rho= _puCorr[leptonIndex] / 0.1562;
  if (std::abs (_lEta[leptonIndex]) >= 1.479 && std::abs (_lEta[leptonIndex]) < 2.00) rho= _puCorr[leptonIndex] / 0.1032;
  if (std::abs (_lEta[leptonIndex]) >= 2.00 && std::abs (_lEta[leptonIndex]) < 2.20) rho= _puCorr[leptonIndex] / 0.0859;
  if (std::abs (_lEta[leptonIndex]) >= 2.20 && std::abs (_lEta[leptonIndex]) < 2.30) rho= _puCorr[leptonIndex] / 0.1116;
  if (std::abs (_lEta[leptonIndex]) >= 2.30 && std::abs (_lEta[leptonIndex]) < 2.40) rho= _puCorr[leptonIndex] / 0.1321;
  if (std::abs (_lEta[leptonIndex]) >= 2.40 && std::abs (_lEta[leptonIndex]) < 2.50) rho= _puCorr[leptonIndex] / 0.1654;
  double _x = _lEleIsEB[leptonIndex] ? 1.12       : 0.5 ;
  double _y = _lEleIsEB[leptonIndex] ? 0.0368       : 0.201 ;
  elehadronicOverEm = 	_lElehadronicOverEm[leptonIndex] - (_x + _y*rho)/_lEleEcalEnergy[leptonIndex];
	
  if(_lElefull5x5SigmaIetaIeta[leptonIndex]                  >= (_lEleIsEB[leptonIndex] ? 0.11       : 0.0314 ))       return false;
  if(_lEleDEtaInSeed  [leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.00477    : 0.00868))       return false;
  if(_lEleDeltaPhiSuperClusterTrackAtVtx [leptonIndex]       >= (_lEleIsEB[leptonIndex] ? 0.222      : 0.213  ))       return false;
  if(elehadronicOverEm                                       >= (_lEleIsEB[leptonIndex] ? 0.298      : 0.101  ))       return false;
  if(_lEleInvMinusPInv[leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.241      : 0.14   ))       return false;
  return true;

}



//______________________________________________muon medium ID
bool Analysis_mc::muOurMedium(const unsigned leptonIndex) const{
  if( !isMu(leptonIndex)) return false; 
  bool _isOurMedium = false;
  bool goodGlob = false;
  goodGlob = _lGlobalMuon[leptonIndex] && _lCQChi2Position[leptonIndex] < 12   &&  _lCQTrackKink[leptonIndex] < 20;
  _isOurMedium = _lPOGLoose[leptonIndex] &&  _lMuonSegComp[leptonIndex] > (goodGlob ? 0.303 : 0.451);
  return _isOurMedium;
}
//______________________________________________muon medium ID
bool Analysis_mc::muTimeVeto(const unsigned leptonIndex) const{
  if( !isMu(leptonIndex)) return false;
  
  bool _passTimingVeto = false;
  _passTimingVeto = true;
  bool cmbok =(_lMuTimenDof[leptonIndex] >7   );
  bool rpcok =(_lMuRPCTimenDof[leptonIndex] >1  && _lMuRPCTimeErr[leptonIndex]==0);
  if (rpcok){
    if ((fabs(_lMuRPCTime[leptonIndex])>10) &&	!(cmbok && fabs(_lMuTime[leptonIndex])<10))   _passTimingVeto=false;
  }
  else{
    if (cmbok && (_lMuTime[leptonIndex]>20 || _lMuTime[leptonIndex]<-45)) _passTimingVeto=false;
  }
  return _passTimingVeto;
}

//______________________________________________FO definition!
bool Analysis_mc::lepIsFOBase(const unsigned leptonIndex) const{
  if( isTau(leptonIndex) ) return false;
  //ele
  if( isElectron(leptonIndex) && !eleIsClean(leptonIndex)) return false;
  if( isElectron(leptonIndex) && !eleIsLoose(leptonIndex)) return false;
  if (isElectron(leptonIndex) && _lPt[leptonIndex] < ele_pt) return false;
  if ( isElectron(leptonIndex) && _relIso[leptonIndex] > ele_iso_loose) return false;
  //mu
  if( isMu(leptonIndex) && !muOurMedium(leptonIndex)) return false;
  //if( isMu(leptonIndex) && !muTimeVeto(leptonIndex)) return false;
  if (isMu(leptonIndex) && _lPt[leptonIndex] < mu_pt) return false;
  if ( isMu(leptonIndex) && _relIso[leptonIndex] > mu_iso_loose) return false;

  if (fabs(_dz[leptonIndex]) > 10) return false;
  
  return true;
} 

//______________________________________________tight definition!
bool Analysis_mc::lepIsTightDisplaced(const unsigned leptonIndex) const{
  if (!lepIsFOBase(leptonIndex)) return false;
  //ele
  if ( isElectron(leptonIndex) && _relIso[leptonIndex] > ele_iso_tight) return false;
  //mu
  //if( isMu(leptonIndex) && !muTimeVeto(leptonIndex)) return false;
  if (isMu(leptonIndex) && _relIso[leptonIndex] > mu_iso_tight) return false;

  return true;
}

//______________________________________________tight definition for l1!
bool Analysis_mc::lepIsTightPrompt(const unsigned leptonIndex) const{
  if (!lepIsFOBase(leptonIndex)) return false;
  //Variable
  if (_relIso[leptonIndex] > 0.1)  return false;
  //if (fabs(_3dIPSig[leptonIndex]) > 4)  return false;
  if (fabs(_dz[leptonIndex]) > 0.1)  return false;
  if (fabs(_dxy[leptonIndex]) > 0.05)  return false;
  //ID
  if (isMu(leptonIndex) && !_lPOGMedium[leptonIndex]) return false;
  //if( isMu(leptonIndex) && !muTimeVeto(leptonIndex)) return false;
  if (isElectron(leptonIndex) && !elePassMVA(leptonIndex)) return false;
  //pT
  if( is2016() && isMu(leptonIndex) && _lPt[leptonIndex] < mu_2016) return false;
  if( is2017() && isMu(leptonIndex) && _lPt[leptonIndex] < mu_2017) return false;
  if( is2018() && isMu(leptonIndex) && _lPt[leptonIndex] < mu_2018) return false;
  if( is2016() && isElectron(leptonIndex) && _lPt[leptonIndex] < ele_2016) return false;
  if( is2017() && isElectron(leptonIndex) && _lPt[leptonIndex] < ele_2017) return false;
  if( is2018() && isElectron(leptonIndex) && _lPt[leptonIndex] < ele_2018) return false;
  //if (!lepPromptTriggerMatching(leptonIndex)) return false;

  return true;
}
//______________________________________________trigger matching for prompt leptons!
bool Analysis_mc::lepPromptTriggerMatching(const unsigned leptonIndex) const{
  if (!lepIsFOBase(leptonIndex)) return false;
  //if (!lepIsTightPrompt(leptonIndex)) return false;
  
  if( is2016() && isMu(leptonIndex) && !((_lHasTrigger[leptonIndex] & 1)==0 && (_lHasTrigger[leptonIndex] & 2)==0)) return true;
  if( is2016() && isElectron(leptonIndex) && !((_lHasTrigger[leptonIndex] & 1)==0)) return true;

  // std::cout<< isMu(leptonIndex)<<"   "<< "   "<< "  "<<(_lHasTrigger[leptonIndex] & 1) <<" "<< (_lHasTrigger[leptonIndex] & 2)<<std::endl;
  if( is2017() && isMu(leptonIndex) && !((_lHasTrigger[leptonIndex] & 1)==0 && (_lHasTrigger[leptonIndex] & 2)==0)) return true;
  if( is2017() && isElectron(leptonIndex) && !((_lHasTrigger[leptonIndex] & 1)==0)) return true;

  if( is2018() && isMu(leptonIndex) && !((_lHasTrigger[leptonIndex] & 1)==0 && (_lHasTrigger[leptonIndex] & 2)==0)) return true;
  if( is2018() && isElectron(leptonIndex) && !((_lHasTrigger[leptonIndex] & 1)==0)) return true;

  // If none of the above (to avoid compilation warnings)
  return false;
}

//______________________________________________ele MVA category
double Analysis_mc::findEleMVACategory(const unsigned leptonIndex) const{
    double category = -1;
    if( !isElectron(leptonIndex)) category = -1;
    if (_lPt[leptonIndex] < 10 && fabs(_lEtaSC[leptonIndex]) < 0.800) category =0;
    if (_lPt[leptonIndex] < 10 && fabs(_lEtaSC[leptonIndex]) >= 0.800 && fabs(_lEtaSC[leptonIndex]) < 1.479) category =1;
    if (_lPt[leptonIndex] < 10 && fabs(_lEtaSC[leptonIndex]) >= 1.479) category =2;	
    if (_lPt[leptonIndex] >= 10 && fabs(_lEtaSC[leptonIndex]) < 0.800) category =3;
    if (_lPt[leptonIndex] >= 10 && fabs(_lEtaSC[leptonIndex]) >= 0.800 && fabs(_lEtaSC[leptonIndex]) < 1.479) category =4;
    if (_lPt[leptonIndex] >= 10 && fabs(_lEtaSC[leptonIndex]) >= 1.479) category =5;
    return category;	
}

//______________________________________________ele MVA category
double Analysis_mc::convertMVAInRawMva(const unsigned leptonIndex) const{
    return -0.5 * (std::log((1-_lElectronMvaFall17NoIso[leptonIndex])/(1+_lElectronMvaFall17NoIso[leptonIndex])));
}

//______________________________________________ele MVA ID 2016
bool Analysis_mc::elePassMVA2016(const unsigned leptonIndex) const{
    if( !isElectron(leptonIndex)) return false;
    bool _passedMVA90 =false;
    double mvaRaw = convertMVAInRawMva(leptonIndex);
    double category = findEleMVACategory(leptonIndex);
     if (category == 0)_passedMVA90 = mvaRaw > 2.77072387339 - exp(-_lPt[leptonIndex] / 3.81500912145) * 8.16304860178;
     if (category == 1)_passedMVA90 = mvaRaw >1.85602317813 - exp(-_lPt[leptonIndex] / 2.18697654938) * 11.8568936824;
     if (category == 2)_passedMVA90 = mvaRaw >1.73489307814 - exp(-_lPt[leptonIndex] / 2.0163211971) * 17.013880078;
     if (category == 3)_passedMVA90 = mvaRaw >5.9175992258 - exp(-_lPt[leptonIndex] / 13.4807294538) * 9.31966232685;
     if (category == 4)_passedMVA90 = mvaRaw >5.01598837255 - exp(-_lPt[leptonIndex] / 13.1280451502) * 8.79418193765;
     if (category == 5)_passedMVA90 = mvaRaw >4.16921343208 - exp(-_lPt[leptonIndex] / 13.2017224621) * 9.00720913211;	
     return _passedMVA90;
}
//______________________________________________ele MVA ID 2016
bool Analysis_mc::elePassMVA2017(const unsigned leptonIndex) const{
    if( !isElectron(leptonIndex)) return false;
    bool _passedMVA90 =false;
    double mvaRaw = convertMVAInRawMva(leptonIndex);
    double category = findEleMVACategory(leptonIndex);
     if (category == 0)_passedMVA90 = mvaRaw > 2.77072387339 - exp(-_lPt[leptonIndex] / 3.81500912145) * 8.16304860178;
     if (category == 1)_passedMVA90 = mvaRaw >1.85602317813 - exp(-_lPt[leptonIndex] / 2.18697654938) * 11.8568936824;
     if (category == 2)_passedMVA90 = mvaRaw >1.73489307814 - exp(-_lPt[leptonIndex] / 2.0163211971) * 17.013880078;
     if (category == 3)_passedMVA90 = mvaRaw >5.9175992258 - exp(-_lPt[leptonIndex] / 13.4807294538) * 9.31966232685;
     if (category == 4)_passedMVA90 = mvaRaw >5.01598837255 - exp(-_lPt[leptonIndex] / 13.1280451502) * 8.79418193765;
     if (category == 5)_passedMVA90 = mvaRaw >4.16921343208 - exp(-_lPt[leptonIndex] / 13.2017224621) * 9.00720913211;	
     return _passedMVA90;
}
//______________________________________________ele MVA ID 2016
bool Analysis_mc::elePassMVA2018(const unsigned leptonIndex) const{
    if( !isElectron(leptonIndex)) return false;
    bool _passedMVA90 =false;
    double mvaRaw = convertMVAInRawMva(leptonIndex);
    double category = findEleMVACategory(leptonIndex);
     if (category == 0)_passedMVA90 = mvaRaw > 2.77072387339 - exp(-_lPt[leptonIndex] / 3.81500912145) * 8.16304860178;
     if (category == 1)_passedMVA90 = mvaRaw >1.85602317813 - exp(-_lPt[leptonIndex] / 2.18697654938) * 11.8568936824;
     if (category == 2)_passedMVA90 = mvaRaw >1.73489307814 - exp(-_lPt[leptonIndex] / 2.0163211971) * 17.013880078;
     if (category == 3)_passedMVA90 = mvaRaw >5.9175992258 - exp(-_lPt[leptonIndex] / 13.4807294538) * 9.31966232685;
     if (category == 4)_passedMVA90 = mvaRaw >5.01598837255 - exp(-_lPt[leptonIndex] / 13.1280451502) * 8.79418193765;
     if (category == 5)_passedMVA90 = mvaRaw >4.16921343208 - exp(-_lPt[leptonIndex] / 13.2017224621) * 9.00720913211;	
     return _passedMVA90;
}
//______________________________________________ele MVA ID
bool Analysis_mc::elePassMVA(const unsigned leptonIndex) const{
  if( !isElectron(leptonIndex)) return false;
  bool _passedMVA90 =false;
  if( is2016()) _passedMVA90 = elePassMVA2016(leptonIndex);
  if( is2017()) _passedMVA90 = elePassMVA2017(leptonIndex);
  if( is2018()) _passedMVA90 = elePassMVA2018(leptonIndex);
  return _passedMVA90;
  //return _lPOGMedium[leptonIndex];
}
//______________________________________________jet cleaning
bool Analysis_mc::jetIsCleanBase(const unsigned jetIndex, bool (Analysis_mc::*leptonIsFO)(const unsigned) const) const{
  for(unsigned l = 0; l < _nLight; ++l){
    if( (this->*leptonIsFO)(l)){
      double deltaR = kinematics::deltaR(_lPhi[l], _lEta[l], _jetPhi[jetIndex], _jetEta[jetIndex]);
      if(deltaR < 0.4){
	return false;
      }
    }
  }
  return true;
}
//______________________________________________jet cleaning
bool Analysis_mc::jetIsClean(const unsigned jetIndex) const{
  return jetIsCleanBase(jetIndex, &Analysis_mc::lepIsFOBase);
}
//______________________________________________jet ID
bool Analysis_mc::jetIsGood(const unsigned jetIndex, int pt_variation) const{
  if (fabs(_jetEta[jetIndex]) > 2.4) return false;
  if (!jetIsClean(jetIndex)) return false;
  if (!_jetIsTight[jetIndex] ) return false;
  if (pt_variation < jet_pt_cut) return false;
  return true;
}
//______________________________________________jet ID
bool Analysis_mc::jetIsBJet(const unsigned jetIndex, int pt_variation) const{
  if (!jetIsGood(jetIndex, pt_variation)) return false;
  if (is2016() && deepCSV(jetIndex) < bjet_loose_2016) return false;
  if (is2017() && deepCSV(jetIndex) < bjet_loose_2017) return false;
  if (is2018() && deepCSV(jetIndex) < bjet_loose_2018) return false;
  return true;
}

//______________________________________________jet b tagged
double Analysis_mc::deepCSV(const unsigned jetIndex) const{
  return _jetDeepCsv_b[jetIndex] + _jetDeepCsv_bb[jetIndex];
}
