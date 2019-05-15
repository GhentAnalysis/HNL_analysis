#include "../interface/Analysis_mc.h"

/*
* order objects by their pT
*/
//______________________________________________pt order
void Analysis_mc::orderByPt(std::vector<unsigned>& ind, const double* pt, const unsigned count) const{
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
double Analysis_mc::coneCorr(const unsigned ind) const{
    double corr = 1.;
    if( lepIsFOBase(ind) && !lepIsTightDisplaced(ind) ){
      if (isMu(ind)) corr *=( 1. + std::max(_relIso[ind] - mu_iso_tight, 0.) );
      if (isElectron(ind)) corr *=( 1. + std::max(_relIso[ind] - ele_iso_tight, 0.) );
    }
    return corr;
}
//______________________________________________cone correction
void Analysis_mc::applyConeCorrection(){
    for(unsigned l = 0; l < _nLight; ++l){
        double coneC = coneCorr(l);
        _lPt[l] *= coneC;
        _lE[l] *= coneC;
    }
}
//______________________________________________order only good leptons
unsigned Analysis_mc::selectLep(std::vector<unsigned>& ind) const{
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
unsigned Analysis_mc::selectLepConeCorr(std::vector<unsigned>& ind){
    //apply cone correction
    applyConeCorrection();
    //select and order cone-corrected leptons
    return selectLep(ind);
}

//______________________________________________find index of l1
int Analysis_mc::l1Index(const std::vector<unsigned>& ind){
  int index_leading = -1;
  int counter_leading=0;
  for(unsigned l = 0; l < ind.size(); ++l){
    //&& lepPromptTriggerMatching(ind[l])
    if (counter_leading == 0){
      if (lepIsTightPrompt(ind[l]) ) {   
	++counter_leading;
	index_leading = ind[l];

      }//only good tight prompt
    }//only 1
  }//loop light
  return index_leading;
}
//______________________________________________check if they can be called displaced pair
bool Analysis_mc::IsDisplacedPair(const unsigned leptonIndex1,const unsigned leptonIndex2, int index_taken_by_l1, std::vector<unsigned>& ind) const{
  int number_found_verteces=0;
  if (leptonIndex1 == index_taken_by_l1) return false;
  if (leptonIndex2 == index_taken_by_l1) return false;
  if (!lepIsFOBase(leptonIndex1)) return false;
  if (!lepIsFOBase(leptonIndex2)) return false;
  if (fabs(_dxy[leptonIndex1]) < dxy_cut) return false;
  if (fabs(_dxy[leptonIndex2]) < dxy_cut) return false;
  if (_lCharge[leptonIndex1] == _lCharge[leptonIndex2]) return false;

  for(unsigned v = 0; v < _nVFit; ++v){
    if (vertex_found(leptonIndex1,leptonIndex2,  _vertices[v][0]) ) ++number_found_verteces;
    if (vertex_found(leptonIndex1,leptonIndex2,  _vertices[v][0]) ) std::cout<<"        vertex found: "<<_vertices[v][0]<<"  "<<_vertices[v][1]<< "   what was passed: "<<leptonIndex1<<"  "<<leptonIndex2<<std::endl;

  }//loop vertecies
  
  if (number_found_verteces <= 0) return false;
  return true;

}
//______________________________________________check if they can be called displaced
bool Analysis_mc::lepIsDisplaced(const unsigned leptonIndex, int index_taken_by_l1, std::vector<unsigned>& ind) const{
  int number_found_verteces=0;
  //std::cout<<"                   we should run using the fixed one that is: "<<leptonIndex<<std::endl;
  //std::cout<<"      in lepIsdisplaced:"<<std::endl;
  if (leptonIndex == index_taken_by_l1) return false;
  //std::cout<<"      no leading"<<std::endl;

  if (!lepIsFOBase(leptonIndex)) return false;
  //std::cout<<"      si FO"<<std::endl;

  if (fabs(_dxy[leptonIndex]) < dxy_cut) return false;
  //std::cout<<"      si dxy cut"<<std::endl;
  //std::cout<<_lPt[leptonIndex]<<std::endl;
  //looking for a common vertex with an other lepton
  for(unsigned sd = 0; sd < ind.size(); ++sd){
    //std::cout<<ind[sd]<<"            --> in lepIsdisplaced:"<<std::endl;

    if (leptonIndex == ind[sd]) continue;
    //std::cout<<"            --> no the same:"<<std::endl;

    if (ind[sd] == index_taken_by_l1) continue;
    //std::cout<<"            --> noleading"<<std::endl;

    if (fabs(_dxy[ind[sd]]) < dxy_cut) continue;
    //std::cout<<"            --> si displaced"<<std::endl;

    if (!lepIsFOBase(ind[sd])) continue;
    // std::cout<<"            --> si FO"<<std::endl;
    //std::cout<<_lPt[ind[sd]]<<"   "<<ind[sd]<<std::endl;

    for(unsigned v = 0; v < _nVFit; ++v){
      //std::cout<<v<<": "<< _vertices[v][0]<<"               what is supposed to be: "<< leptonIndex<<"  "<<ind[sd]<<std::endl;
      if (vertex_found(leptonIndex,ind[sd],  _vertices[v][0]) ) ++number_found_verteces;
    }//loop vertecies
  }//loop second lepton
  //std::cout<<"number of common vertex: "<<number_found_verteces<<std::endl;

  if (number_found_verteces <= 0) return false;

  return true;
}

//______________________________________________function di check if 2 indeces make a vertex
bool Analysis_mc::vertex_found(const unsigned leptonIndex1, const unsigned leptonIndex2, int vertex_index) const{
  int Index1 = leptonIndex1+1;
  int Index2 = leptonIndex2+1;
  //std::cout<<"================> in the fucking vertex finding!!! "<<"vertex: "<<vertex_index<<"   "<<Index1<<"    "<< Index2<<std::endl;
  //std::cout<<"prova1: "<<(Index1*100 + Index2)<<"   prova 2 "<< (Index1 + Index2*100)<<std::endl;
  bool fuck = false;
  fuck = (vertex_index == (Index1*100 + Index2) ) || (vertex_index == (Index1 + Index2*100) );
  //std::cout<<"fuck "<<fuck<<std::endl;
  return ((vertex_index == (Index1*100 + Index2) ) || (vertex_index == (Index1 + Index2*100) ) );
}
//______________________________________________function di check if 2 indeces make a vertex
int Analysis_mc::l2l3_vertex_variable(const unsigned leptonIndex1, const unsigned leptonIndex2) {
  int Index1 = leptonIndex1+1;
  int Index2 = leptonIndex2+1;

  int indice=-1;
  for(unsigned v = 0; v < _nVFit; ++v){
    if (  (_vertices[v][0]   == (Index1*100 + Index2) ) || (_vertices[v][0] == (Index1 + Index2*100) )) {
      indice = v;
    }
  }//loop vertecies
  return indice;
}

//______________________________________________check if lepton comes from matrix-element conversion
bool Analysis_mc::lepFromMEExtConversion(const unsigned leptonIndex) const{
    bool fromConversion = (_lMatchPdgId[leptonIndex] == 22);
    // bool promptConversion = (_lIsPrompt[leptonIndex] && _lProvenanceConversion[leptonIndex] == 0);
    bool promptConversion = (_lIsPrompt[leptonIndex] && _lProvenanceCompressed[leptonIndex] == 0);
    return (fromConversion && promptConversion);
}
//______________________________________________reject events with overlapping photon-production phase-space
bool Analysis_mc::photonOverlap(const Sample& samp, const bool mcNonprompt) const{

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
            if(lepIsFOBase(l) && lepFromMEExtConversion(l) ){
                return true;
            }
        } 
    //require photon samples to have atlease one external conversion
    } if(isPhotonSample){
        bool hasConversion = false;
        for(unsigned l = 0; l < _nLight; ++l){
            if(lepIsFOBase(l) && lepFromMEExtConversion(l) ){
                hasConversion = true;
            }
        }
        return !(hasConversion);
    }
    return false;
}
//______________________________________________
bool Analysis_mc::photonOverlap(const bool mcNonprompt) const{
    return photonOverlap(currentSample, mcNonprompt);
}










