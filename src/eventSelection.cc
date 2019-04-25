#include "../interface/Analysis_mc.h"

/*
* order objects by their pT
*/
//______________________________________________pt order
void Analyis_mc::orderByPt(std::vector<unsigned>& ind, const double* pt, const unsigned count) const{
    std::vector<std::pair<double, unsigned>> ptMap;
    for(unsigned p = 0; p < count; ++p){
        ptMap.push_back({pt[ind[p]], ind[p]});
    }
    std::sort(ptMap.begin(), ptMap.end(), [](std::pair<double, unsigned>& p1, std::pair<double, unsigned>& p2){ return p1.first > p2.first; } );
    for(unsigned p = 0; p < count; ++p){
        ind[p] = ptMap[p].second;
    }
}
//______________________________________________cone correction
double Analyis_mc::coneCorr(const unsigned ind) const{
    double corr = 1.;
    if( lepIsFOBase(ind) && !lepIsTightDisplaced(ind) ){
      if (isMu(leptonIndex)) corr *=( 1. + std::max(_relIso[ind] - mu_iso_tight, 0.) );
      if (isElectron(leptonIndex)) corr *=( 1. + std::max(_relIso[ind] - ele_iso_tight, 0.) );
    }
    return corr;
}
//______________________________________________cone correction
void Analyis_mc::applyConeCorrection(){
    for(unsigned l = 0; l < _nLight; ++l){
        double coneC = coneCorr(l);
        _lPt[l] *= coneC;
        _lE[l] *= coneC;
    }
}
//______________________________________________order only good leptons
unsigned Analyis_mc::selectLep(std::vector<unsigned>& ind) const{
    ind.clear();
    unsigned lCount = 0;
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsFOBase(l)){
            ++lCount;
            ind.push_back(l);
        }
    }
    if(lCount < 3) return 0;
    orderByPt(ind, _lPt, lCount);
    return lCount;	
}
//______________________________________________order only good leptons with cone corrections
unsigned Analyis_mc::selectLepConeCorr(std::vector<unsigned>& ind){
    //apply cone correction
    applyConeCorrection();
    //select and order cone-corrected leptons
    return selectLep(ind);
}

//______________________________________________find index of l1
int Analyis_mc::l1Index(std::vector<unsigned>& ind){
  int index_leading = -1;
  int counter_leading=0;
  for(unsigned l = 0; l < lCount; ++l){
    if (counter_leading == 0){
      if (lepIsTightPrompt(ind[l]) && lepPromptTriggerMatching(ind[l])) {   
	++counter_leading;
	index_leading = ind[l];

      }//only good tight prompt
    }//only 1
  }//loop light
  return index_leading;
}

//______________________________________________check if they can be called displaced
bool Analyis_mc::lepIsDisplaced(const unsigned leptonIndex, int index_taken_by_l1, std::vector<unsigned>& ind) const{
  int number_found_verteces=-1;
  
  if (leptonIndex == index_taken_by_l1) return false;
  if (!lepIsFOBase(leptonIndex)) return false;
  if (fabs(_dxy[leptonIndex]) > dxy_cut) return false;
  //looking for a common vertex with an other lepton
  for(unsigned sd = 0; sd < lCount; ++sd){
    if (leptonIndex == ind[sd]) continue;
    if (ind[sd] == index_taken_by_l1) continue;
    if (fabs(_dxy[ind[sd]]) > dxy_cut) continue;
    if (!lepIsFOBase(ind[sd])) continue;
    for(unsigned v = 0; v < _nVFit; ++v){
      if (vertex_found(leptonIndex,ind[sd],  _vertices[v][0]) ) ++number_found_verteces;   
    }//loop vertecies
  }//loop second lepton
  if (number_found_verteces <= 0) return false;

  return true;
}

//______________________________________________function di check if 2 indeces make a vertex
bool Analyis_mc::vertex_found(const unsigned leptonIndex1, const unsigned leptonIndex2, int vertex_index) const{
  int Index1 = leptonIndex1+1;
  int Index2 = leptonIndex2+1;
  return (vertex_index == (Index1*100 + Index2) ) || (vertex_index == (Index1 + Index2*100) );
}
















unsigned treeReader::dilFlavorComb(const std::vector<unsigned>& ind) const{
    unsigned flavCount[3] = {0,0,0};
    for(unsigned l = 0; l < 2; ++l) ++flavCount[_lFlavor[ind[l]]];
    if(flavCount[2] == 0){
        if(flavCount[1] == 0) return 0; //ee
        else if(flavCount[1] == 1) return 1; //em
        else return 2; //mm
    } else if(flavCount[2] == 1){
        if(flavCount[1] == 0) return 3; //et
        else return 4; //mt
    }
    return 5; //tt
}


/*
* lepton event selection
*/



//amount of consecutive tight leptons
unsigned treeReader::tightLepCount(const std::vector<unsigned>& ind, const unsigned lCount) const{
    unsigned tightC = 0; 
    for(unsigned l = 0; l < lCount; ++l){
        if(lepIsTight(ind[l])) ++tightC;
        else return tightC;
    }
    return tightC;
}

//lepton pT thresholds
bool treeReader::passPtCuts(const std::vector<unsigned>& ind) const{
    if(_lPt[ind[0]] <= 25) return false;
    if(_lPt[ind[1]] <= 15) return false;
    return true;
}

//check if leptons are prompt
bool treeReader::promptLeptons() const{

    //don't apply this check to data
    if( isData() ){
        return true;
    }
    
    //check whether every MC lepton is prompt
    for(unsigned l = 0; l < _nLight; ++l){
        if(lepIsGood(l) && !_lIsPrompt[l]) return false;
    }
    return true;
}


/*
* number of jets 
*/


/*
* remove overlapping phase-space between different MC samples
*/

//check if lepton comes from matrix-element conversion
bool treeReader::lepFromMEExtConversion(const unsigned leptonIndex) const{
    bool fromConversion = (_lMatchPdgId[leptonIndex] == 22);
    bool promptConversion = (_lIsPrompt[leptonIndex] && _lProvenanceConversion[leptonIndex] == 0);
    return (fromConversion && promptConversion);
}

//reject events with overlapping photon-production phase-space
bool treeReader::photonOverlap(const Sample& samp, const bool mcNonprompt) const{

    //in case of data-driven nonprompt estimation the only sample that is to be cleaned is Wgamma (for overlap with WZ)
    bool isInclusiveSample;
    bool isPhotonSample;
    if( mcNonprompt ){
        isInclusiveSample = (samp.getFileName().find("DYJetsToLL") != std::string::npos) ||
            (samp.getFileName().find("TTTo2L") != std::string::npos) ||
            (samp.getFileName().find("TTJets") != std::string::npos );

        isPhotonSample = (samp.getFileName().find("ZGTo2LG") != std::string::npos) ||
            (samp.getFileName().find("TTGJets") != std::string::npos) ||
            (samp.getFileName().find("WGToLNuG") != std::string::npos);
    } else {
        isInclusiveSample = false;
        isPhotonSample = (samp.getFileName().find("WGToLNuG") != std::string::npos);
    }

    //require inclusive sample to contain no external conversions
    if(isInclusiveSample){
       for(unsigned l = 0; l < _nLight; ++l){
            if(lepIsGood(l) && lepFromMEExtConversion(l) ){
                return true;
            }
        } 
    //require photon samples to have atlease one external conversion
    } if(isPhotonSample){
        bool hasConversion = false;
        for(unsigned l = 0; l < _nLight; ++l){
            if(lepIsGood(l) && lepFromMEExtConversion(l) ){
                hasConversion = true;
            }
        }
        return !(hasConversion);
    }
    return false;
}

bool treeReader::photonOverlap(const bool mcNonprompt) const{
    return photonOverlap(currentSample, mcNonprompt);
}

//overlap removal between inclusive and HT-binned samples
bool treeReader::htOverlap(const Sample& samp) const{
    if(samp.getFileName().find("DYJetsToLL_M-50_Tune") != std::string::npos){
        return _gen_HT > 70.;
    } else if(samp.getFileName().find("DYJetsToLL_M-10_50_Tune") != std::string::npos){
        return _gen_HT > 100.;
    }
    return false;
}

bool treeReader::htOverlap() const{
    return htOverlap(currentSample);
}

/*
* trigger selection
*/

bool treeReader::passSingleLeptonTriggers() const{
    return (_passTrigger_e || _passTrigger_m);
}

bool treeReader::passDileptonTriggers() const{
    return (_passTrigger_ee || _passTrigger_em || _passTrigger_mm);
}

bool treeReader::passTrileptonTriggers() const{
    return (_passTrigger_eee || _passTrigger_eem || _passTrigger_emm || _passTrigger_mmm); 
}

bool treeReader::passTriggerCocktail() const{
    bool pass = passSingleLeptonTriggers() ||
        passDileptonTriggers() ||
        passTrileptonTriggers();
    return pass;
}

bool treeReader::passMETFilters() const{
    return _passMETFilters;
}
