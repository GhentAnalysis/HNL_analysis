
#include "../interface/Analysis_mc.h"
#include "../interface/kinematicTools.h"


//______________________________________________ele base
bool Analyis_mc::eleIsClean(const unsigned electronIndex) const{
    if( is2016() ){
        return eleIsClean2016(electronIndex);
    }
    if else ( is2017() ){
        return eleIsClean2017(electronIndex);
    }
    else{
      return eleIsClean2018(electronIndex);
    }     
}
//______________________________________________muon loose used for ele cleaning base 
bool Analyis_mc::lepIsLoose(const unsigned leptonIndex) const{
  if( is2016() ){
    return lepIsLoose2016(leptonIndex);
  }
  if else ( is2017() ){
      return lepIsLoose2017(leptonIndex);
    }
  else{
    return lepIsLoose2018(leptonIndex);
  }  
}

//______________________________________________ele loose ID
bool Analyis_mc::eleIsLoose(const unsigned leptonIndex) const{
  if( is2016() ){
    return eleIsLoose2016(leptonIndex);
  }
  if else ( is2017() ){
      return eleIsLoose2017(leptonIndex);
    }
  else{
    return eleIsLoose2018(leptonIndex);
  } 

}



//______________________________________________remove electrons overlapping with muons
bool Analyis_mc::eleIsCleanBase(const unsigned electronIndex, bool (Analyis_mc::*looseMuon)(const unsigned) const) const{
    //make sure this lepton is an electron
    if( !isElectron(electronIndex) ){
        std::cerr << "Error: trying to clean non-electron object from muon overlap." << std::endl;
        return 999;
    }
    //check separation with every muon
    for(unsigned m = 0; m < _nMu; ++m){
        if( !isMuon(m) ){
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
bool Analyis_mc::eleIsClean2016(const unsigned electronIndex) const{
    return eleIsCleanBase(electronIndex, &treeReader::lepIsLoose2016);
}
//______________________________________________ele clean base 17
bool Analyis_mc::eleIsClean2017(const unsigned electronIndex) const{
    return eleIsCleanBase(electronIndex, &treeReader::lepIsLoose2017);
}
//______________________________________________ele clean base 18
bool Analyis_mc::eleIsClean2018(const unsigned electronIndex) const{
    return eleIsCleanBase(electronIndex, &treeReader::lepIsLoose2018);
}

//______________________________________________loose POG 16
bool Analyis_mc::lepIsLoose2016(const unsigned leptonIndex) const{
    return _lPOGLoose[leptonIndex];
}
//______________________________________________loose POG 17
bool Analyis_mc::lepIsLoose2017(const unsigned leptonIndex) const{
   return _lPOGLoose[leptonIndex];
}
//______________________________________________loose POG 18
bool Analyis_mc::lepIsLoose2018(const unsigned leptonIndex) const{
   return _lPOGLoose[leptonIndex];
}
//______________________________________________ele loose ID
bool Analyis_mc::eleIsLoose2016(const unsigned leptonIndex) const{
  if( !isElectron(leptonIndex)) return false;
  
  if (!(_lEleIsEB[leptonIndex] || _lEleIsEE[leptonIndex])) return false;
  if(_lElefull5x5SigmaIetaIeta[leptonIndex]                  >= (_lEleIsEB[leptonIndex] ? 0.11       : 0.0314 ))       return false;
  if(_lEleDEtaInSeed  [leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.00477    : 0.00868))       return false;
  if(_lEleDeltaPhiSuperClusterTrackAtVtx [leptonIndex]       >= (_lEleIsEB[leptonIndex] ? 0.222      : 0.213  ))       return false;
  if(_lElehadronicOverEm[leptonIndex]                        >= (_lEleIsEB[leptonIndex] ? 0.298      : 0.101  ))       return false;
  if(_lEleInvMinusPInv[leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.241      : 0.14   ))       return false;
  return true;

}
//______________________________________________ele loose ID
bool Analyis_mc::eleIsLoose2017(const unsigned leptonIndex) const{
  if( !isElectron(leptonIndex)) return false;
  
  if (!(_lEleIsEB[leptonIndex] || _lEleIsEE[leptonIndex])) return false;
  if(_lElefull5x5SigmaIetaIeta[leptonIndex]                  >= (_lEleIsEB[leptonIndex] ? 0.11       : 0.0314 ))       return false;
  if(_lEleDEtaInSeed  [leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.00477    : 0.00868))       return false;
  if(_lEleDeltaPhiSuperClusterTrackAtVtx [leptonIndex]       >= (_lEleIsEB[leptonIndex] ? 0.222      : 0.213  ))       return false;
  if(_lElehadronicOverEm[leptonIndex]                        >= (_lEleIsEB[leptonIndex] ? 0.298      : 0.101  ))       return false;
  if(_lEleInvMinusPInv[leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.241      : 0.14   ))       return false;
  return true;

}
//______________________________________________ele loose ID
bool Analyis_mc::eleIsLoose2018(const unsigned leptonIndex) const{
  if( !isElectron(leptonIndex)) return false;
  
  if (!(_lEleIsEB[leptonIndex] || _lEleIsEE[leptonIndex])) return false;
  if(_lElefull5x5SigmaIetaIeta[leptonIndex]                  >= (_lEleIsEB[leptonIndex] ? 0.11       : 0.0314 ))       return false;
  if(_lEleDEtaInSeed  [leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.00477    : 0.00868))       return false;
  if(_lEleDeltaPhiSuperClusterTrackAtVtx [leptonIndex]       >= (_lEleIsEB[leptonIndex] ? 0.222      : 0.213  ))       return false;
  if(_lElehadronicOverEm[leptonIndex]                        >= (_lEleIsEB[leptonIndex] ? 0.298      : 0.101  ))       return false;
  if(_lEleInvMinusPInv[leptonIndex]                          >= (_lEleIsEB[leptonIndex] ? 0.241      : 0.14   ))       return false;
  return true;

}



//______________________________________________muon medium ID
bool Analyis_mc::muOurMedium(const unsigned leptonIndex) const{
  if( !isMu(leptonIndex)) return false; 
  bool _isOurMedium = false;
  bool goodGlob = false;
  goodGlob = _lGlobalMuon[leptonIndex] && _lCQChi2Position[leptonIndex] < 12   &&  _lCQTrackKink[leptonIndex] < 20;
  _isOurMedium = _lPOGLoose[leptonIndex] &&  _muonSegComp[leptonIndex] > (goodGlob ? 0.303 : 0.451);
  return _isOurMedium;
}
//______________________________________________muon medium ID
bool Analyis_mc::muTimeVeto(const unsigned leptonIndex) const{
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
bool Analyis_mc::lepIsFOBase(const unsigned leptonIndex) const{
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
bool Analyis_mc::lepIsTightDisplaced(const unsigned leptonIndex) const{
  if (!lepIsFOBase(leptonIndex)) return false;
  //ele
  if ( isElectron(leptonIndex) && _relIso[leptonIndex] > ele_iso_tight) return false;
  //mu
  if( isMu(leptonIndex) && !muTimeVeto(leptonIndex)) return false;
  if (isMu(leptonIndex) && _relIso[leptonIndex] > mu_iso_tight) return false;

  return true;
}

//______________________________________________tight definition for l1!
bool Analyis_mc::lepIsTightPrompt(const unsigned leptonIndex) const{
  if (!lepIsFOBase(leptonIndex)) return false;
  //Variable
  if (_relIso[leptonIndex] > 0.1)  return false;
  if (fabs(_3dIPSig[leptonIndex]) > 4)  return false;
  if (fabs(_dz[leptonIndex]) > 0.1)  return false;
  if (fabs(_dxy[leptonIndex]) > 0.05)  return false;
  //ID
  if (isMu(leptonIndex) && !_lPOGMedium(leptonIndex)) return false;
  if( isMu(leptonIndex) && !muTimeVeto(leptonIndex)) return false;
  if (isElectron(leptonIndex) && !elePassMVA[leptonIndex]) return false;
  //pT
  if( is2016() && isMu(leptonIndex) && _lPt[leptonIndex] < mu_2016) return false;
  if( is2017() && isMu(leptonIndex) && _lPt[leptonIndex] < mu_2017) return false;
  if( is2018() && isMu(leptonIndex) && _lPt[leptonIndex] < mu_2018) return false;
  if( is2016() && isElectron(leptonIndex) && _lPt[leptonIndex] < ele_2016) return false;
  if( is2017() && isElectron(leptonIndex) && _lPt[leptonIndex] < ele_2017) return false;
  if( is2018() && isElectron(leptonIndex) && _lPt[leptonIndex] < ele_2018) return false;
  
  return true;
}
//______________________________________________trigger matching for prompt leptons!
bool Analyis_mc::lepPromptTriggerMatching(const unsigned leptonIndex) const{
  if (!lepIsFOBase(leptonIndex)) return false;
  if (!lepIsTightPrompt(leptonIndex)) return false;
  
  if( is2016() && isMu(leptonIndex) && ((_lHasTrigger[leptonIndex] & 1)==0 || (_lHasTrigger[leptonIndex] & 2)==0)) return true;
  if( is2016() && isElectron(leptonIndex) && (_lHasTrigger[leptonIndex] & 1)==0) return true;

  if( is2017() && isMu(leptonIndex) && ((_lHasTrigger[leptonIndex] & 1)==0 || (_lHasTrigger[leptonIndex] & 2)==0)) return true;
  if( is2017() && isElectron(leptonIndex) && (_lHasTrigger[leptonIndex] & 1)==0) return true;

  if( is2018() && isMu(leptonIndex) && ((_lHasTrigger[leptonIndex] & 1)==0 || (_lHasTrigger[leptonIndex] & 2)==0)) return true;
  if( is2018() && isElectron(leptonIndex) && (_lHasTrigger[leptonIndex] & 1)==0) return true;
}


//______________________________________________ele MVA ID
bool Analyis_mc::elePassMVA(const unsigned leptonIndex) const{
  if( !isElectron(leptonIndex)) return false;
  bool _passedMVA90 =false;
  int eta = -1;
  if(TMath::Abs(_lEta[leptonIndex]) < 0.8 ) eta = 0;
  else if(TMath::Abs(_lEta[leptonIndex]) < 1.479 ) eta = 1;
  else eta = 2;
  _passedMVA90 = __lElectronMvaFall17Iso[leptonIndex] >  std::min( MVA_cuts_pt15[eta], std::max(MVA_cuts_pt25[eta] , MVA_cuts_pt15[eta] + (MVA_cuts_pt25[eta] - MVA_cuts_pt15[eta])*0.1 *( _lPt[leptonIndex]-15) ) );
  return _passedMVA90;
}


//______________________________________________jet cleaning
bool Analyis_mc::jetIsCleanBase(const unsigned jetIndex, bool (Analyis_mc::*leptonIsFO)(const unsigned) const) const{
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
bool Analyis_mc::jetIsClean(const unsigned jetIndex) const{
    return jetIsCleanBase(jetIndex, &Analyis_mc::lepIsFOBase);
}
//______________________________________________jet ID
bool Analyis_mc::jetIsGood(const unsigned jetIndex) const{
  if (fabs(_jetEta[jetIndex]) > 2.4) return false;
  if (!jetIsClean[jetIndex]) return false;
  if (!_jetIsTight[jetIndex] ) return false;
  if (_jetPt[jetIndex] < jet_pt_cut) return false;
  return true;
}
//______________________________________________jet ID
bool Analyis_mc::jetIsBJet(const unsigned jetIndex) const{
  if (!jetIsGood[jetIndex]) return false;
  if (is2016() && deepCSV(jetIndex) < bjet_loose_2016) return false;
  if (is2017() && deepCSV(jetIndex) < bjet_loose_2017) return false;
  if (is2018() && deepCSV(jetIndex) < bjet_loose_2018) return false;
  return true;
}

//______________________________________________jet b tagged
double Analyis_mc::deepCSV(const unsigned jetIndex) const{
    return _jetDeepCsv_b[jetIndex] + _jetDeepCsv_bb[jetIndex];
}
