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
int Analysis_mc::l1Index(std::vector<unsigned>& ind){
  int index_leading = -1;
  int counter_leading=0;
  std::cout<<index_leading<<std::endl;
  for(unsigned l = 0; l < lCount; ++l){
    std::cout<<l<<") "<<std::endl;

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
bool Analysis_mc::lepIsDisplaced(const unsigned leptonIndex, int index_taken_by_l1, std::vector<unsigned>& ind) const{
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
bool Analysis_mc::vertex_found(const unsigned leptonIndex1, const unsigned leptonIndex2, int vertex_index) const{
  int Index1 = leptonIndex1+1;
  int Index2 = leptonIndex2+1;
  return (vertex_index == (Index1*100 + Index2) ) || (vertex_index == (Index1 + Index2*100) );
}







