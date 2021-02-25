
#include "../interface/Analysis_mc.h"

//_____________________________________________

int Analysis_mc::SR_bin_muon(int channel,double D2_delta_pv_sv,  double M_l2l3_combined ){
  int bin = -1;
  double value_mass_cut = 4.0 ;
  double value_displaced_first_cut = 0.5 ;
  double value_displaced_second_cut = 1.5 ;
  double value_displaced_third_cut = 4 ;
	
  bool mass_bool_less = false;	
  bool mass_bool_more = false; 	
  bool disp_bool_first = false;
  bool disp_bool_second = false;
  bool disp_bool_third = false;
  bool disp_bool_fourth = false;

  mass_bool_less = M_l2l3_combined < value_mass_cut;
  mass_bool_more = M_l2l3_combined >= value_mass_cut;
  disp_bool_first =  D2_delta_pv_sv < 	value_displaced_first_cut;
  disp_bool_second =  D2_delta_pv_sv > 	value_displaced_first_cut  && D2_delta_pv_sv < 	value_displaced_second_cut;
  disp_bool_third =  D2_delta_pv_sv  > 	value_displaced_second_cut && D2_delta_pv_sv < 	value_displaced_third_cut ;
  disp_bool_fourth =  D2_delta_pv_sv > 	value_displaced_third_cut;	
		
  if (channel == 0 || channel == 1 || channel == 2  ){
    if (channel == 0){
      if (mass_bool_less && disp_bool_first)    bin =1;
      if (mass_bool_less && disp_bool_second)   bin =2;
      if (mass_bool_less && disp_bool_third)    bin =3;
      if (mass_bool_less && disp_bool_fourth)   bin =4;

      if (mass_bool_more && disp_bool_first)    bin =5;
      if (mass_bool_more && D2_delta_pv_sv > 	value_displaced_first_cut)   bin =6;
      
      //if (mass_bool_more && disp_bool_first)    bin =5;
      //if (mass_bool_more && disp_bool_second)   bin =6;
      //if (mass_bool_more && disp_bool_third)    bin =7;
      //if (mass_bool_more && disp_bool_fourth)   bin =8;
	    
    }
    // numeration vould be different with 24 bins
    if (channel == 1){
      if (mass_bool_less && disp_bool_first)    bin =7;
      if (mass_bool_less && disp_bool_second)   bin =8;
      if (mass_bool_less && disp_bool_third)    bin =9;
      if (mass_bool_less && disp_bool_fourth)   bin =10;
      if (mass_bool_more && disp_bool_first)    bin =11;
      if (mass_bool_more && D2_delta_pv_sv > 	value_displaced_first_cut)   bin =12;
      //if (mass_bool_more && disp_bool_third)    bin =15;
      //if (mass_bool_more && disp_bool_fourth)   bin =16;

    }
    if (channel == 2){
      if (mass_bool_less && disp_bool_first)    bin =13;
      if (mass_bool_less && disp_bool_second)   bin =14;
      if (mass_bool_less && disp_bool_third)    bin =15;
      if (mass_bool_less && disp_bool_fourth)   bin =16;
      if (mass_bool_more && disp_bool_first)    bin =17;
      if (mass_bool_more && D2_delta_pv_sv > 	value_displaced_first_cut)   bin =18;
      //if (mass_bool_more && disp_bool_third)    bin =23;
      //if (mass_bool_more && disp_bool_fourth)   bin =24;
    }
  }
  return bin;
}
//_____________________________________________

int Analysis_mc::SR_bin_ele(int channel,double D2_delta_pv_sv,  double M_l2l3_combined ){
  int bin = -1;
  double value_mass_cut = 4.0 ;
  double value_displaced_first_cut = 0.5 ;
  double value_displaced_second_cut = 1.5 ;
  double value_displaced_third_cut = 4 ;
	
  bool mass_bool_less = false;	
  bool mass_bool_more = false; 	
  bool disp_bool_first = false;
  bool disp_bool_second = false;
  bool disp_bool_third = false;
  bool disp_bool_fourth = false;

  mass_bool_less = M_l2l3_combined < value_mass_cut;
  mass_bool_more = M_l2l3_combined >= value_mass_cut;
   disp_bool_first =  D2_delta_pv_sv < 	value_displaced_first_cut;
  disp_bool_second =  D2_delta_pv_sv > 	value_displaced_first_cut  && D2_delta_pv_sv < 	value_displaced_second_cut;
  disp_bool_third =  D2_delta_pv_sv  > 	value_displaced_second_cut && D2_delta_pv_sv < 	value_displaced_third_cut ;
  disp_bool_fourth =  D2_delta_pv_sv > 	value_displaced_third_cut;	
	
  if (channel == 3 || channel == 4 || channel == 5  ){
    if (channel == 3){
      if (mass_bool_less && disp_bool_first)    bin =1;
      if (mass_bool_less && disp_bool_second)   bin =2;
      if (mass_bool_less && disp_bool_third)    bin =3;
      if (mass_bool_less && disp_bool_fourth)   bin =4;
      if (mass_bool_more && disp_bool_first)    bin =5;
      if (mass_bool_more && D2_delta_pv_sv > 	value_displaced_first_cut)   bin =6;
      //if (mass_bool_more && disp_bool_third)    bin =7;
      //if (mass_bool_more && disp_bool_fourth)   bin =8;
    }
    if (channel == 4){
      if (mass_bool_less && disp_bool_first)    bin =7;
      if (mass_bool_less && disp_bool_second)   bin =8;
      if (mass_bool_less && disp_bool_third)    bin =9;
      if (mass_bool_less && disp_bool_fourth)   bin =10;
      if (mass_bool_more && disp_bool_first)    bin =11;
      if (mass_bool_more && D2_delta_pv_sv > 	value_displaced_first_cut)   bin =12;
      //if (mass_bool_more && disp_bool_third)    bin =15;
      //if (mass_bool_more && disp_bool_fourth)   bin =16;

    }
    if (channel == 5){
      if (mass_bool_less && disp_bool_first)    bin =13;
      if (mass_bool_less && disp_bool_second)   bin =14;
      if (mass_bool_less && disp_bool_third)    bin =15;
      if (mass_bool_less && disp_bool_fourth)   bin =16;
      if (mass_bool_more && disp_bool_first)    bin =17;
      if (mass_bool_more && D2_delta_pv_sv > 	value_displaced_first_cut)   bin =18;
      //if (mass_bool_more && disp_bool_third)    bin =23;
      //if (mass_bool_more && disp_bool_fourth)   bin =24;

    }
  }
  return bin;
}
//_____________________________________________
int Analysis_mc::channel(int  flavors_3l[3], int  charge_3l[3]){
  int canale=-1;
  if (flavors_3l[0] == flavors_3l[1] && flavors_3l[0] == flavors_3l[2] ) {
    if (flavors_3l[0] == 0) canale = 3;
    if (flavors_3l[0] == 1) canale = 0;
  }
  if (flavors_3l[0] == 1){
    if (flavors_3l[1]==1 && flavors_3l[2]==0){
      if (charge_3l[0] != charge_3l[1]) canale = 1;	    
      if (charge_3l[0] == charge_3l[1]) canale = 2;
    }
    if (flavors_3l[2]==1 && flavors_3l[1]==0){
      if (charge_3l[0] != charge_3l[2]) canale = 1;	    
      if (charge_3l[0] == charge_3l[2]) canale = 2;
    }   
  }//mme

  if (flavors_3l[0] == 0){
    if (flavors_3l[1]==0 && flavors_3l[2]==1){
      if (charge_3l[0] != charge_3l[1]) canale = 4;	    
      if (charge_3l[0] == charge_3l[1]) canale = 5;
    }
    if (flavors_3l[2]==0 && flavors_3l[1]==1){
      if (charge_3l[0] != charge_3l[2]) canale = 4;	    
      if (charge_3l[0] == charge_3l[2]) canale = 5;
    }   
  }//eem

  // 0 = mmm
  // 1 = mme OS
  // 2 = mme SS
  // 3 = eee
  // 4 = eem OS
  // 5 = eem SS
  return canale;
  
}




//______________________________________________
bool Analysis_mc::resonanceVeto(int channel, double D2_delta_pv_sv, int  flavors_3l[3], int  charge_3l[3], double M_l2l3_combined ,double M_l1l2_combined, double M_l1l3_combined){
  bool veto = false;
  
  //* veto for resonances!
  bool j_psi_veto_l2l3 = true;
  bool psi_2_veto_l2l3 = true;
  bool omega_veto_l2l3 = true;
  bool phi_veto_l2l3 = true;	
  bool j_psi_veto_l1l2 = true;
  bool psi_2_veto_l1l2 = true;
  bool omega_veto_l1l2 = true;
  bool phi_veto_l1l2 = true;
  bool z_veto_l1l2 = true;
  bool upsilon_veto_l1l2 = true;
  bool upsilon2_veto_l1l2 = true;
  bool upsilon3_veto_l1l2 = true;	
  bool j_psi_veto_l1l3 = true;
  bool psi_2_veto_l1l3 = true;
  bool omega_veto_l1l3 = true;
  bool phi_veto_l1l3 = true;
  bool z_veto_l1l3 = true;
  bool upsilon_veto_l1l3 = true;
  bool upsilon2_veto_l1l3 = true;
  bool upsilon3_veto_l1l3 = true;

      
  if (channel == 0 || channel == 3){
    if (D2_delta_pv_sv < 1.5 && fabs (M_l2l3_combined - 3.0969) < 0.08 ) j_psi_veto_l2l3 = false;
    if (D2_delta_pv_sv < 1.5 && fabs (M_l2l3_combined - 3.6861) < 0.08 ) psi_2_veto_l2l3 = false;
    if (D2_delta_pv_sv < 1.5 && fabs (M_l2l3_combined - 0.7827) < 0.08 ) omega_veto_l2l3 = false;
    if (D2_delta_pv_sv < 1.5 && fabs (M_l2l3_combined - 1.0190) < 0.08 ) phi_veto_l2l3 = false;
		
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 0.7827) < 0.08 ) omega_veto_l1l2 = false;
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 1.0190) < 0.08 ) phi_veto_l1l2 = false;
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 3.0969) < 0.08 ) j_psi_veto_l1l2 = false;
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 3.6861) < 0.08 ) psi_2_veto_l1l2 = false;
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 9.4603) < 0.08 ) upsilon_veto_l1l2 = false;
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 10.0233) < 0.08 ) upsilon2_veto_l1l2 = false;
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 10.3552) < 0.08 ) upsilon3_veto_l1l2 = false;
    if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 91.1876) < 10 )  z_veto_l1l2 = false;

    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 0.7827) < 0.08 ) omega_veto_l1l3 = false;
    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 1.0190) < 0.08 ) phi_veto_l1l3 = false;
    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 3.0969) < 0.08 ) j_psi_veto_l1l3 = false;
    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 3.6861) < 0.08 ) psi_2_veto_l1l3 = false;
    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 9.4603) < 0.08 ) upsilon_veto_l1l3 = false;
    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 10.0233) < 0.08 ) upsilon2_veto_l1l3 = false;
    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 10.3552) < 0.08 ) upsilon3_veto_l1l3 = false;
    if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 91.1876) < 10 )  z_veto_l1l3 = false;
  }
  if (channel == 1 || channel == 4){
    if (flavors_3l[0] == flavors_3l[1]){
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 0.7827) < 0.08 ) omega_veto_l1l2 = false;
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 1.0190) < 0.08 ) phi_veto_l1l2 = false;
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 3.0969) < 0.08 ) j_psi_veto_l1l2 = false;
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 3.6861) < 0.08 ) psi_2_veto_l1l2 = false;
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 9.4603) < 0.08 ) upsilon_veto_l1l2 = false;
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 10.0233) < 0.08 ) upsilon2_veto_l1l2 = false;
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 10.3552) < 0.08 ) upsilon3_veto_l1l2 = false;
      if (charge_3l[0]!= charge_3l[1] && fabs (M_l1l2_combined - 91.1876) < 10 )  z_veto_l1l2 = false;
    }
    if (flavors_3l[0] == flavors_3l[2]){
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 0.7827) < 0.08 ) omega_veto_l1l3 = false;
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 1.0190) < 0.08 ) phi_veto_l1l3 = false;
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 3.0969) < 0.08 ) j_psi_veto_l1l3 = false;
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 3.6861) < 0.08 ) psi_2_veto_l1l3 = false;
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 9.4603) < 0.08 ) upsilon_veto_l1l3 = false;
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 10.0233) < 0.08 ) upsilon2_veto_l1l3 = false;
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 10.3552) < 0.08 ) upsilon3_veto_l1l3 = false;
      if (charge_3l[0]!= charge_3l[2] && fabs (M_l1l3_combined - 91.1876) < 10 )  z_veto_l1l3 = false;
    }
  }

  
  veto=   j_psi_veto_l2l3  && psi_2_veto_l2l3  && omega_veto_l2l3  && phi_veto_l2l3 &&
    j_psi_veto_l1l2  && psi_2_veto_l1l2  && omega_veto_l1l2  && phi_veto_l1l2  && z_veto_l1l2  && upsilon_veto_l1l2  && upsilon2_veto_l1l2  && upsilon3_veto_l1l2 
    && j_psi_veto_l1l3  && psi_2_veto_l1l3  && omega_veto_l1l3  && phi_veto_l1l3  && z_veto_l1l3  && upsilon_veto_l1l3  && upsilon2_veto_l1l3  && upsilon3_veto_l1l3;

  return veto;
}


//_____________________________________________ SF prompt ele
double Analysis_mc::SF_btag_eff(TH2F *sf_btag_eff[3], const double eta, const double pt, const int flav){
   double sfValue = 1;
   if (flav == 0) {   
   int binx = sf_btag_eff[0]->GetXaxis()->FindBin(pt);
   int biny = sf_btag_eff[0]->GetYaxis()->FindBin(fabs(eta));
   sfValue = sf_btag_eff[0]->GetBinContent(binx,biny);
   }
   if (flav == 4) {   
   int binx = sf_btag_eff[1]->GetXaxis()->FindBin(pt);
   int biny = sf_btag_eff[1]->GetYaxis()->FindBin(fabs(eta));
   sfValue = sf_btag_eff[1]->GetBinContent(binx,biny);
   }
   if (flav == 5) {   
   int binx = sf_btag_eff[2]->GetXaxis()->FindBin(pt);
   int biny = sf_btag_eff[2]->GetYaxis()->FindBin(fabs(eta));
   sfValue = sf_btag_eff[2]->GetBinContent(binx,biny);

   }
   return sfValue;	
}

//=============================================================================================================================================================================
//===================================================================          SF leptons          ============================================================================
//_____________________________________________ SF prompt muon ID
double Analysis_mc::SF_prompt_muon(TH2D *muon_sf_histogram[1],TH2F *muon_sf_isoIP_histogram[1], const unsigned leptonIndex){
  double sfValue = 1;
  double sfValue_ID = 1;
  double sfValue_IsoIP = 1;
 
  int binx_ID =0;
  int biny_ID =0;
  if (muon_sf_histogram[0]->GetXaxis()->GetXmax() > 4) { // it means that xaxis has pt so y axis has |eta|
    biny_ID = muon_sf_histogram[0]->GetYaxis()->FindBin(std::abs(_lEta[leptonIndex]));
    if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetXaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetXaxis()->GetNbins()) )   binx_ID =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
    else binx_ID = muon_sf_histogram[0]->GetXaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetXaxis()->GetBinLowEdge(1)));  
  }
  else { // it means that xaxis has eta so y axis has pt
    binx_ID = muon_sf_histogram[0]->GetXaxis()->FindBin(_lEta[leptonIndex]);
    if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetYaxis()->GetNbins()) )   biny_ID =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
    else biny_ID = muon_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1)));   
  }
  sfValue_ID = muon_sf_histogram[0]->GetBinContent(binx_ID,biny_ID);

  //made by kirill has x absolute eta, while y pt
  int binx_IsoIP =0;
  int biny_IsoIP =0;
  binx_IsoIP = muon_sf_isoIP_histogram[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndex]));
  if (_lPt[leptonIndex] > muon_sf_isoIP_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_isoIP_histogram[0]->GetYaxis()->GetNbins()))       biny_IsoIP =  muon_sf_isoIP_histogram[0]->GetYaxis()->GetNbins(); 
  else biny_IsoIP = muon_sf_isoIP_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_isoIP_histogram[0]->GetYaxis()->GetBinLowEdge(1)));
  sfValue_IsoIP = muon_sf_isoIP_histogram[0]->GetBinContent(binx_IsoIP,biny_IsoIP);

  sfValue = sfValue_IsoIP*sfValue_ID;
  if (sfValue == 0) std::cout<<" -------------------  WARNING the SF (ID+ISO+IP) for prompt muon is ZERO  ---------------"<<std::endl;
  return sfValue;	
}
//_____________________________________________ SF prompt muon error
double Analysis_mc::SF_prompt_muon_error(TH2D *muon_sf_histogram_syst[1],TH2D *muon_sf_histogram[1],TH2F *muon_sf_isoIP_histogram[1],TH2F *muon_sf_isoIP_histogram_syst[1],  const unsigned leptonIndex){
  double sfValue = 1;
  double sfValue_ID_1 = 1;
  double sfValue_ID_2 = 1;
  double sfValue_ID = 1;
  double sfValue_IsoIP_1 = 1;
  double sfValue_IsoIP_2 = 1;
  double sfValue_IsoIP = 1;

  int binx_ID_1 =0;
  int biny_ID_1 =0;
  if (muon_sf_histogram_syst[0]->GetXaxis()->GetXmax() > 4) { // it means that xaxis has pt so y axis has |eta|
    biny_ID_1 = muon_sf_histogram_syst[0]->GetYaxis()->FindBin(std::abs(_lEta[leptonIndex]));
    if (_lPt[leptonIndex] > muon_sf_histogram_syst[0]->GetXaxis()->GetBinUpEdge(muon_sf_histogram_syst[0]->GetXaxis()->GetNbins()) )binx_ID_1 =  muon_sf_histogram_syst[0]->GetYaxis()->GetNbins(); 
    else binx_ID_1 = muon_sf_histogram_syst[0]->GetXaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram_syst[0]->GetXaxis()->GetBinLowEdge(1)));  
  }
  else { // it means that xaxis has eta so y axis has pt
    binx_ID_1 = muon_sf_histogram_syst[0]->GetXaxis()->FindBin(_lEta[leptonIndex]);
    if (_lPt[leptonIndex] > muon_sf_histogram_syst[0]->GetYaxis()->GetBinUpEdge(muon_sf_histogram_syst[0]->GetYaxis()->GetNbins()) )biny_ID_1 =  muon_sf_histogram_syst[0]->GetYaxis()->GetNbins(); 
    else biny_ID_1 = muon_sf_histogram_syst[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram_syst[0]->GetYaxis()->GetBinLowEdge(1)));   
  }			
  sfValue_ID_1 = muon_sf_histogram_syst[0]->GetBinErrorLow(binx_ID_1,biny_ID_1);

  int binx_ID_2 =0;
  int biny_ID_2 =0;
  if (muon_sf_histogram[0]->GetXaxis()->GetXmax() > 4) { // it means that xaxis has pt so y axis has |eta|
    biny_ID_2 = muon_sf_histogram[0]->GetYaxis()->FindBin(std::abs(_lEta[leptonIndex]));
    if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetXaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetXaxis()->GetNbins()) )binx_ID_2 =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
    else binx_ID_2 = muon_sf_histogram[0]->GetXaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetXaxis()->GetBinLowEdge(1)));  
  }
  else { // it means that xaxis has eta so y axis has pt
    binx_ID_2 = muon_sf_histogram[0]->GetXaxis()->FindBin(_lEta[leptonIndex]);
    if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetYaxis()->GetNbins()) )biny_ID_2 =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
    else biny_ID_2 = muon_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1)));   
  }			
  sfValue_ID_2 = muon_sf_histogram[0]->GetBinErrorLow(binx_ID_2,biny_ID_2);
  sfValue_ID= std::max(sfValue_ID_1,sfValue_ID_2);


  int binx_IsoIP_1 =0;
  int biny_IsoIP_1 =0;
  binx_IsoIP_1 = muon_sf_isoIP_histogram[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndex]));
  if (_lPt[leptonIndex] > muon_sf_isoIP_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_isoIP_histogram[0]->GetYaxis()->GetNbins()))       biny_IsoIP_1 =  muon_sf_isoIP_histogram[0]->GetYaxis()->GetNbins(); 
  else biny_IsoIP_1 = muon_sf_isoIP_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_isoIP_histogram[0]->GetYaxis()->GetBinLowEdge(1)));
  sfValue_IsoIP_1 = muon_sf_isoIP_histogram[0]->GetBinErrorLow(binx_IsoIP_1,biny_IsoIP_1);
  
  int binx_IsoIP_2 =0;
  int biny_IsoIP_2 =0;
  binx_IsoIP_2 = muon_sf_isoIP_histogram_syst[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndex]));
  if (_lPt[leptonIndex] > muon_sf_isoIP_histogram_syst[0]->GetYaxis()->GetBinUpEdge(muon_sf_isoIP_histogram_syst[0]->GetYaxis()->GetNbins()))       biny_IsoIP_2 =  muon_sf_isoIP_histogram_syst[0]->GetYaxis()->GetNbins(); 
  else biny_IsoIP_2 = muon_sf_isoIP_histogram_syst[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_isoIP_histogram_syst[0]->GetYaxis()->GetBinLowEdge(1)));
  sfValue_IsoIP_2 =muon_sf_isoIP_histogram_syst[0]->GetBinContent(binx_IsoIP_2,biny_IsoIP_2);
  sfValue_IsoIP= std::max(sfValue_IsoIP_1,sfValue_IsoIP_2);


  sfValue = TMath::Sqrt(sfValue_IsoIP*sfValue_IsoIP + sfValue_ID*sfValue_ID);

  return sfValue;		
}
//_____________________________________________ SF prompt muon trigger
double Analysis_mc::SF_trigger_muon(TH2F *muon_sf_histogram[1], const unsigned leptonIndex){
   double sfValue = 1;  
   int binx =0;
   int biny =0;	
   binx = muon_sf_histogram[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndex]));
   if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetYaxis()->GetNbins()) )biny =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
   else biny = muon_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1))); 
   sfValue = muon_sf_histogram[0]->GetBinContent(binx,biny);
   return sfValue;	
}
//_____________________________________________ SF prompt muon trigger error
double Analysis_mc::SF_trigger_muon_error(TH2F *muon_sf_histogram[1], const unsigned leptonIndex){
     double sfValue = 1;  
   int binx =0;
   int biny =0;	
   binx = muon_sf_histogram[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndex]));
   if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetYaxis()->GetNbins()) )biny =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
   else biny = muon_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1))); 
   sfValue = muon_sf_histogram[0]->GetBinErrorLow(binx,biny);	
   return sfValue;			
}
//_____________________________________________ SF prompt ele
double Analysis_mc::SF_prompt_ele(TH2F *ele_sf_histogram[1], const unsigned leptonIndex){
   double sfValue = 1;  	
   int binx = ele_sf_histogram[0]->GetXaxis()->FindBin(_lEtaSC[leptonIndex]);
   int biny = 0.;
  if (_lPt[leptonIndex] > ele_sf_histogram[0]->GetYaxis()->GetBinUpEdge(ele_sf_histogram[0]->GetYaxis()->GetNbins()) )biny =  ele_sf_histogram[0]->GetYaxis()->GetNbins(); 
  else biny = ele_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], ele_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1))); 		
   sfValue = ele_sf_histogram[0]->GetBinContent(binx,biny);	
   return sfValue;	
}

//_____________________________________________ SF prompt ele
double Analysis_mc::SF_prompt_ele_error(TH2F *ele_sf_histogram[1], const unsigned leptonIndex){
   double sfValue = 1;  	
   int binx = ele_sf_histogram[0]->GetXaxis()->FindBin(_lEtaSC[leptonIndex]);
   int biny=0.;	
   if (_lPt[leptonIndex] > ele_sf_histogram[0]->GetYaxis()->GetBinUpEdge(ele_sf_histogram[0]->GetYaxis()->GetNbins()) )biny =  ele_sf_histogram[0]->GetYaxis()->GetNbins(); 
   else biny = ele_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], ele_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1))); 	
   sfValue = ele_sf_histogram[0]->GetBinErrorLow(binx,biny);	
   return sfValue;	
}

//_____________________________________________ SF prompt muon trigger
double Analysis_mc::SF_trigger_ele(TH2F *muon_sf_histogram[1], const unsigned leptonIndex){
   double sfValue = 1;  
   int binx =0;
   int biny =0;	
   binx = muon_sf_histogram[0]->GetXaxis()->FindBin(std::abs(_lEtaSC[leptonIndex]));
   if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetYaxis()->GetNbins()) )biny =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
   else biny = muon_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1))); 
   sfValue = muon_sf_histogram[0]->GetBinContent(binx,biny);
   return sfValue;	
}
//_____________________________________________ SF prompt muon trigger error
double Analysis_mc::SF_trigger_ele_error(TH2F *muon_sf_histogram[1], const unsigned leptonIndex){
     double sfValue = 1;  
   int binx =0;
   int biny =0;	
   binx = muon_sf_histogram[0]->GetXaxis()->FindBin(std::abs(_lEtaSC[leptonIndex]));
   if (_lPt[leptonIndex] > muon_sf_histogram[0]->GetYaxis()->GetBinUpEdge(muon_sf_histogram[0]->GetYaxis()->GetNbins()) )biny =  muon_sf_histogram[0]->GetYaxis()->GetNbins(); 
   else biny = muon_sf_histogram[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndex], muon_sf_histogram[0]->GetYaxis()->GetBinLowEdge(1))); 
   sfValue = muon_sf_histogram[0]->GetBinErrorLow(binx,biny);	
   return sfValue;			
}




//_____________________________________________ displaced mess
double Analysis_mc::displaced_weight (int  flavors_3l[3],int channel,unsigned _lElectronMissingHits_l2, unsigned _lElectronMissingHits_l3, double displ2, double displ3, double sum_pt, double D2_delta_pv_sv, double displEleVars[8], TH2F *sf_sv_effcy_num[1], TH2F *sf_sv_effcy_den[1], TH2F *sf_isoID_nPMuon[1],TH2F *sf_isoID_nPMuon_syst[1],const unsigned leptonIndexl2,const unsigned leptonIndexl3){
  
  double weight =1.;

  if (channel == 3){// ee case, conversion studies as SF 
    int indElel2 = -1;
    int indElel3 = -1;
    if (displ2 <= 8 ) indElel2 =0;
    else if (displ2 > 8  && displ2 <= 15)  indElel2 =1;
    else if (displ2 > 15  && displ2 <= 20) indElel2 =2;
    else if (displ2 > 20  && displ2 <= 25) indElel2 =3;
    else if (displ2 > 25  && displ2 <= 30) indElel2 =4;
    else if (displ2 > 30  && displ2 <= 37) indElel2 =5;
    else if (displ2 > 37  && displ2 <= 50) indElel2 =6;
    else  indElel2 =7;
    if (displ3 <= 8 ) indElel3 =0;
    else if (displ3 > 8  && displ3 <= 15)  indElel3 =1;
    else if (displ3 > 15  && displ3 <= 20) indElel3 =2;
    else if (displ3 > 20  && displ3 <= 25) indElel3 =3;
    else if (displ3 > 25  && displ3 <= 30) indElel3 =4;
    else if (displ3 > 30  && displ3 <= 37) indElel3 =5;
    else if (displ3 > 37  && displ3 <= 50) indElel3 =6;
    else   indElel3 =7; 
    weight *= displEleVars[indElel2];	
    weight *= displEleVars[indElel3];
  }//end ee
 
  if (channel == 0 ){ // mm case, luka's studies on the SV AND per-muon SF(ID and ISO) measured in Zevents by kirill
    int binx_n,biny_n,binx_d,biny_d  =0;
    double xvariable, yvariable = 0.;
    xvariable = D2_delta_pv_sv;
    yvariable = sum_pt;
    if (xvariable > 20) xvariable = 7.;
    if (yvariable > 20) yvariable = 10.;
    binx_n = sf_sv_effcy_num[0] ->GetXaxis()->FindBin(xvariable);
    biny_n = sf_sv_effcy_num[0] ->GetYaxis()->FindBin(yvariable);
    binx_d = sf_sv_effcy_den[0] ->GetXaxis()->FindBin(xvariable);
    biny_d = sf_sv_effcy_den[0] ->GetYaxis()->FindBin(yvariable); 
    weight *=  (sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n))    /    (sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d)) ; //only from luka study
    //std::cout<<"SF function ----> weight   "<<weight<< "   which is "<< sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n)<<"  /  "<<sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d)<<std::endl;
    // std::cout<<"   "<<D2_delta_pv_sv<<"   "<< sum_pt<<std::endl;

    double sfValue_IsoID_l2 =1.;
    int binx_IsoID_l2 =0;
    int biny_IsoID_l2 =0;
    binx_IsoID_l2 = sf_isoID_nPMuon[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl2]));
    if (_lPt[leptonIndexl2] > sf_isoID_nPMuon[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon[0]->GetYaxis()->GetNbins()))       biny_IsoID_l2 =  sf_isoID_nPMuon[0]->GetYaxis()->GetNbins(); 
    else biny_IsoID_l2 = sf_isoID_nPMuon[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl2], sf_isoID_nPMuon[0]->GetYaxis()->GetBinLowEdge(1)));
    sfValue_IsoID_l2 = sf_isoID_nPMuon[0]->GetBinContent(binx_IsoID_l2,biny_IsoID_l2); // SF l2

    double sfValue_IsoID_l3 =1.;
    int binx_IsoID_l3 =0;
    int biny_IsoID_l3 =0;
    binx_IsoID_l3 = sf_isoID_nPMuon[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl3]));
    if (_lPt[leptonIndexl3] > sf_isoID_nPMuon[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon[0]->GetYaxis()->GetNbins()))       biny_IsoID_l3 =  sf_isoID_nPMuon[0]->GetYaxis()->GetNbins();
    else biny_IsoID_l3 = sf_isoID_nPMuon[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl3], sf_isoID_nPMuon[0]->GetYaxis()->GetBinLowEdge(1)));
    sfValue_IsoID_l3 = sf_isoID_nPMuon[0]->GetBinContent(binx_IsoID_l3,biny_IsoID_l3); // SF l3

    weight *=sfValue_IsoID_l2 * sfValue_IsoID_l3; // all the SF for the mm case 
  }

  if (channel == 1 || channel==2 || channel == 4 || channel ==5 ){ // em case, electron SF from conversion, sqrt SV from luka, single muon SF from kirill
    if (flavors_3l[1] == 1 && flavors_3l[2] == 0){ // me
      int binx_n,biny_n,binx_d,biny_d  =0;
      double xvariable, yvariable = 0.;
      xvariable = D2_delta_pv_sv;
      yvariable = sum_pt;
      if (xvariable > 20) xvariable = 7.;
      if (yvariable > 20) yvariable = 10.;
      binx_n = sf_sv_effcy_num[0] ->GetXaxis()->FindBin(xvariable);
      biny_n = sf_sv_effcy_num[0] ->GetYaxis()->FindBin(yvariable);
      binx_d = sf_sv_effcy_den[0] ->GetXaxis()->FindBin(xvariable);
      biny_d = sf_sv_effcy_den[0] ->GetYaxis()->FindBin(yvariable);
      weight *=  TMath::Sqrt((sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n))    /    (sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d))) ;
      
      
      int indElel3 = -1;
      if (displ3 <= 8 ) indElel3 =0;
      else if (displ3 > 8  && displ3 <= 15)  indElel3 =1;
      else if (displ3 > 15  && displ3 <= 20) indElel3 =2;
      else if (displ3 > 20  && displ3 <= 25) indElel3 =3;
      else if (displ3 > 25  && displ3 <= 30) indElel3 =4;
      else if (displ3 > 30  && displ3 <= 37) indElel3 =5;
      else if (displ3 > 37  && displ3 <= 50) indElel3 =6;
      else    indElel3 =7; 
      weight *= displEleVars[indElel3];

      double sfValue_IsoID_l2 =1.;
      int binx_IsoID_l2 =0;
      int biny_IsoID_l2 =0;
      binx_IsoID_l2 = sf_isoID_nPMuon[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl2]));
      if (_lPt[leptonIndexl2] > sf_isoID_nPMuon[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon[0]->GetYaxis()->GetNbins()))       biny_IsoID_l2 =  sf_isoID_nPMuon[0]->GetYaxis()->GetNbins(); 
      else biny_IsoID_l2 = sf_isoID_nPMuon[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl2], sf_isoID_nPMuon[0]->GetYaxis()->GetBinLowEdge(1)));
      sfValue_IsoID_l2 = sf_isoID_nPMuon[0]->GetBinContent(binx_IsoID_l2,biny_IsoID_l2); // SF l2
      weight *=  sfValue_IsoID_l2;
    }
    if (flavors_3l[1] == 0 && flavors_3l[2] == 1){ // em
      int binx_n,biny_n,binx_d,biny_d  =0;
      double xvariable, yvariable = 0.;
      xvariable = D2_delta_pv_sv;
      yvariable = sum_pt;
      if (xvariable > 20) xvariable = 7.;
      if (yvariable > 20) yvariable = 10.;    
      binx_n = sf_sv_effcy_num[0] ->GetXaxis()->FindBin(xvariable);
      biny_n = sf_sv_effcy_num[0] ->GetYaxis()->FindBin(yvariable);
      binx_d = sf_sv_effcy_den[0] ->GetXaxis()->FindBin(xvariable);
      biny_d = sf_sv_effcy_den[0] ->GetYaxis()->FindBin(yvariable);
      weight *=  TMath::Sqrt((sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n))    /    (sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d))) ;
          
      int indElel2 = -1;
      if (displ2 <= 8 ) indElel2 =0;
      else if (displ2 > 8  && displ2 <= 15)  indElel2 =1;
      else if (displ2 > 15  && displ2 <= 20) indElel2 =2;
      else if (displ2 > 20  && displ2 <= 25) indElel2 =3;
      else if (displ2 > 25  && displ2 <= 30) indElel2 =4;
      else if (displ2 > 30  && displ2 <= 37) indElel2 =5;
      else if (displ2 > 37  && displ2 <= 50) indElel2 =6;
      else    indElel2 =7;
      weight *= displEleVars[indElel2];

      double sfValue_IsoID_l3 =1.;
      int binx_IsoID_l3 =0;
      int biny_IsoID_l3 =0;
      binx_IsoID_l3 = sf_isoID_nPMuon[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl3]));
      if (_lPt[leptonIndexl3] > sf_isoID_nPMuon[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon[0]->GetYaxis()->GetNbins()))       biny_IsoID_l3 =  sf_isoID_nPMuon[0]->GetYaxis()->GetNbins();
      else biny_IsoID_l3 = sf_isoID_nPMuon[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl3], sf_isoID_nPMuon[0]->GetYaxis()->GetBinLowEdge(1)));
      sfValue_IsoID_l3 = sf_isoID_nPMuon[0]->GetBinContent(binx_IsoID_l3,biny_IsoID_l3); // SF l3
      weight *=  sfValue_IsoID_l3;
    }
  }
  //std::cout<<"weight at the end: "<<weight<<std::endl;
  return weight;
}
//_____________________________________________ displaced mess
double Analysis_mc::displaced_weight_error (int  flavors_3l[3],int channel,unsigned _lElectronMissingHits_l2, unsigned _lElectronMissingHits_l3,double displ2, double displ3,  double sum_pt, double D2_delta_pv_sv, double displEleVars[8], TH2F *sf_sv_effcy_num[1], TH2F *sf_sv_effcy_den[1], TH2F *sf_isoID_nPMuon[1],TH2F *sf_isoID_nPMuon_syst[1],const unsigned leptonIndexl2,const unsigned leptonIndexl3){
  double weight_error =1.;

  if (channel == 3){ // ee case, error
    double eleele=1;
    int indElel2 = -1;
    int indElel3 = -1;
    if (displ2 <= 8 ) indElel2 =0;
    else if (displ2 > 8  && displ2 <= 15)  indElel2 =1;
    else if (displ2 > 15  && displ2 <= 20) indElel2 =2;
    else if (displ2 > 20  && displ2 <= 25) indElel2 =3;
    else if (displ2 > 25  && displ2 <= 30) indElel2 =4;
    else if (displ2 > 30  && displ2 <= 37) indElel2 =5;
    else if (displ2 > 37  && displ2 <= 50) indElel2 =6;
    else   indElel2 =7;
    if (displ3 <= 8 ) indElel3 =0;
    else if (displ3 > 8  && displ3 <= 15)  indElel3 =1;
    else if (displ3 > 15  && displ3 <= 20) indElel3 =2;
    else if (displ3 > 20  && displ3 <= 25) indElel3 =3;
    else if (displ3 > 25  && displ3 <= 30) indElel3 =4;
    else if (displ3 > 30  && displ3 <= 37) indElel3 =5;
    else if (displ3 > 37  && displ3 <= 50) indElel3 =6;
    else  indElel3 =7; 
    eleele *= displEleVars[indElel2];	
    eleele *= displEleVars[indElel3];
    weight_error = 0.5* std::abs(1 - eleele);



    
  }
  
  if (channel == 0 ){ // mm is the mess
    weight_error = 1.;
    int binx_n,biny_n,binx_d,biny_d  =0;
    double xvariable, yvariable = 0.;
    xvariable = D2_delta_pv_sv;
    yvariable = sum_pt;
    if (xvariable > 20) xvariable = 7.;
    if (yvariable > 20) yvariable = 10.;
    binx_n = sf_sv_effcy_num[0] ->GetXaxis()->FindBin(xvariable);
    biny_n = sf_sv_effcy_num[0] ->GetYaxis()->FindBin(yvariable);
    binx_d = sf_sv_effcy_den[0] ->GetXaxis()->FindBin(xvariable);
    biny_d = sf_sv_effcy_den[0] ->GetYaxis()->FindBin(yvariable);

    double error_sv =  0.5* std::abs(1 -   ( (sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n))    /    (sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d))) );
    //std::cout<<"   "<<D2_delta_pv_sv<<"   "<< sum_pt<<std::endl;
    //std::cout<<"xvariable "<< xvariable<< "   bin: "<< binx_n<<"     yvariable: "<<yvariable<<"  "<< biny_n<<std::endl;
    //std::cout<<"error_sv   "<<error_sv<< "   which is "<< sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n)<<"  /  "<<sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d)<<std::endl;
    
    double sfValue_IsoID_l2 =1.;
    int binx_IsoID_l2 =0;
    int biny_IsoID_l2 =0;
    binx_IsoID_l2 = sf_isoID_nPMuon[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl2]));
    if (_lPt[leptonIndexl2] > sf_isoID_nPMuon[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon[0]->GetYaxis()->GetNbins()))       biny_IsoID_l2 =  sf_isoID_nPMuon[0]->GetYaxis()->GetNbins(); 
    else biny_IsoID_l2 = sf_isoID_nPMuon[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl2], sf_isoID_nPMuon[0]->GetYaxis()->GetBinLowEdge(1)));
    sfValue_IsoID_l2 = sf_isoID_nPMuon[0]->GetBinContent(binx_IsoID_l2,biny_IsoID_l2); // SF l2
    double sfValue_IsoID_l3 =1.;
    int binx_IsoID_l3 =0;
    int biny_IsoID_l3 =0;
    binx_IsoID_l3 = sf_isoID_nPMuon[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl3]));
    if (_lPt[leptonIndexl3] > sf_isoID_nPMuon[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon[0]->GetYaxis()->GetNbins()))       biny_IsoID_l3 =  sf_isoID_nPMuon[0]->GetYaxis()->GetNbins();
    else biny_IsoID_l3 = sf_isoID_nPMuon[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl3], sf_isoID_nPMuon[0]->GetYaxis()->GetBinLowEdge(1)));
    sfValue_IsoID_l3 = sf_isoID_nPMuon[0]->GetBinContent(binx_IsoID_l3,biny_IsoID_l3); // SF l3

    double sfValue_IsoID_l2_error =1.;
    int binx_IsoID_l2_error =0;
    int biny_IsoID_l2_error =0;
    binx_IsoID_l2_error = sf_isoID_nPMuon_syst[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl2]));
    if (_lPt[leptonIndexl2] > sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins()))       biny_IsoID_l2_error =  sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins(); 
    else biny_IsoID_l2_error = sf_isoID_nPMuon_syst[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl2], sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinLowEdge(1)));
    sfValue_IsoID_l2_error = sf_isoID_nPMuon_syst[0]->GetBinContent(binx_IsoID_l2_error,biny_IsoID_l2_error); // SF l2
    double sfValue_IsoID_l3_error =1.;
    int binx_IsoID_l3_error =0;
    int biny_IsoID_l3_error =0;
    binx_IsoID_l3_error = sf_isoID_nPMuon_syst[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl3]));
    if (_lPt[leptonIndexl3] > sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins()))       biny_IsoID_l3_error =  sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins();
    else biny_IsoID_l3_error = sf_isoID_nPMuon_syst[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl3], sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinLowEdge(1)));
    sfValue_IsoID_l3_error = sf_isoID_nPMuon_syst[0]->GetBinContent(binx_IsoID_l3_error,biny_IsoID_l3_error); // SF l3
    
    double error_idsf = sfValue_IsoID_l2*sfValue_IsoID_l3_error +  sfValue_IsoID_l3*sfValue_IsoID_l2_error + sfValue_IsoID_l2_error*sfValue_IsoID_l3_error;
    //std::cout<<"error_idsf   "<<error_idsf<<std::endl;


    
    weight_error = TMath::Sqrt (error_sv*error_sv + error_idsf*error_idsf);
    //std::cout<<"total error: "<< weight_error <<std::endl;
    
  }

  
  if (channel == 1 || channel==2 || channel == 4 || channel ==5 ){ // em case, electron SF from conversion, sqrt SV from luka, single muon SF from kirill
    if (flavors_3l[1] == 1 && flavors_3l[2] == 0){ // me
      double weight_firstpart = 1.;
     int binx_n,biny_n,binx_d,biny_d  =0;
      double xvariable, yvariable = 0.;
      xvariable = D2_delta_pv_sv;
      yvariable = sum_pt;
      if (xvariable > 20) xvariable = 7.;
      if (yvariable > 20) yvariable = 10.;
      binx_n = sf_sv_effcy_num[0] ->GetXaxis()->FindBin(xvariable);
      biny_n = sf_sv_effcy_num[0] ->GetYaxis()->FindBin(yvariable);
      binx_d = sf_sv_effcy_den[0] ->GetXaxis()->FindBin(xvariable);
      biny_d = sf_sv_effcy_den[0] ->GetYaxis()->FindBin(yvariable);
      weight_firstpart *=  TMath::Sqrt((sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n))    /    (sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d))) ;

      
      
      int indElel3 = -1;
      if (displ3 <= 8 ) indElel3 =0;
      else if (displ3 > 8  && displ3 <= 15)  indElel3 =1;
      else if (displ3 > 15  && displ3 <= 20) indElel3 =2;
      else if (displ3 > 20  && displ3 <= 25) indElel3 =3;
      else if (displ3 > 25  && displ3 <= 30) indElel3 =4;
      else if (displ3 > 30  && displ3 <= 37) indElel3 =5;
      else if (displ3 > 37  && displ3 <= 50) indElel3 =6;
      else   indElel3 =7; 
      weight_firstpart *= displEleVars[indElel3];
      double error_sv = 0.5* std::abs(1 - weight_firstpart);

      double sfValue_IsoID_l2_error =1.;
      int binx_IsoID_l2_error =0;
      int biny_IsoID_l2_error =0;
      binx_IsoID_l2_error = sf_isoID_nPMuon_syst[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl2]));
      if (_lPt[leptonIndexl2] > sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins()))       biny_IsoID_l2_error =  sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins(); 
      else biny_IsoID_l2_error = sf_isoID_nPMuon_syst[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl2], sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinLowEdge(1)));
      sfValue_IsoID_l2_error = sf_isoID_nPMuon_syst[0]->GetBinContent(binx_IsoID_l2_error,biny_IsoID_l2_error); // SF l2

      weight_error = TMath::Sqrt (error_sv*error_sv + sfValue_IsoID_l2_error*sfValue_IsoID_l2_error);
    }
    if (flavors_3l[1] == 0 && flavors_3l[2] == 1){ // em
      double weight_firstpart = 1.;
      int binx_n,biny_n,binx_d,biny_d  =0;
      double xvariable, yvariable = 0.;
      xvariable = D2_delta_pv_sv;
      yvariable = sum_pt;
      if (xvariable > 20) xvariable = 7.;
      if (yvariable > 20) yvariable = 10.;    
      binx_n = sf_sv_effcy_num[0] ->GetXaxis()->FindBin(xvariable);
      biny_n = sf_sv_effcy_num[0] ->GetYaxis()->FindBin(yvariable);
      binx_d = sf_sv_effcy_den[0] ->GetXaxis()->FindBin(xvariable);
      biny_d = sf_sv_effcy_den[0] ->GetYaxis()->FindBin(yvariable);
      weight_firstpart *=  TMath::Sqrt((sf_sv_effcy_num[0]->GetBinContent(binx_n,biny_n))    /    (sf_sv_effcy_den[0]->GetBinContent(binx_d,biny_d))) ;


      int indElel2 = -1;
      if (displ2 <= 8 ) indElel2 =0;
      else if (displ2 > 8  && displ2 <= 15)  indElel2 =1;
      else if (displ2 > 15  && displ2 <= 20) indElel2 =2;
      else if (displ2 > 20  && displ2 <= 25) indElel2 =3;
      else if (displ2 > 25  && displ2 <= 30) indElel2 =4;
      else if (displ2 > 30  && displ2 <= 37) indElel2 =5;
      else if (displ2 > 37  && displ2 <= 50) indElel2 =6;
      else  indElel2 =7;
      weight_firstpart *= displEleVars[indElel2];
      double error_sv = 0.5* std::abs(1 - weight_firstpart);

      double sfValue_IsoID_l3_error =1.;
      int binx_IsoID_l3_error =0;
      int biny_IsoID_l3_error =0;
      binx_IsoID_l3_error = sf_isoID_nPMuon_syst[0]->GetXaxis()->FindBin(std::abs(_lEta[leptonIndexl3]));
      if (_lPt[leptonIndexl3] > sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinUpEdge(sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins()))       biny_IsoID_l3_error =  sf_isoID_nPMuon_syst[0]->GetYaxis()->GetNbins();
      else biny_IsoID_l3_error = sf_isoID_nPMuon_syst[0]->GetYaxis()->FindBin(std::max(_lPt[leptonIndexl3], sf_isoID_nPMuon_syst[0]->GetYaxis()->GetBinLowEdge(1)));
      sfValue_IsoID_l3_error = sf_isoID_nPMuon_syst[0]->GetBinContent(binx_IsoID_l3_error,biny_IsoID_l3_error); // SF l3

      weight_error = TMath::Sqrt (error_sv*error_sv + sfValue_IsoID_l3_error*sfValue_IsoID_l3_error);

    }
  }
  
  return weight_error;
}
//_____________________________________________
void Analysis_mc::zCandidate(TLorentzVector pair[2],TLorentzVector other[1], TLorentzVector leep1, TLorentzVector leep2,TLorentzVector leep3, int  flavors_3l[3], int  charge_3l[3]){
  int ch_lepton1=charge_3l[0];
  int ch_lepton2=charge_3l[1];
  int ch_lepton3=charge_3l[2];
  int fl_lepton1=flavors_3l[0];
  int fl_lepton2=flavors_3l[1];
  int fl_lepton3=flavors_3l[2];
  
  // OSSF
  if (     ((ch_lepton1 != ch_lepton2)    && (fl_lepton1 == fl_lepton2))  || ((ch_lepton1 != ch_lepton3)   && (fl_lepton1 == fl_lepton3)) || ((ch_lepton2 != ch_lepton3)  && (fl_lepton3 == fl_lepton2)) ){ // ossf
    //cout<<"in function where kind is 1: "<<kind[0]<<endl;
        
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
      other[0].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      other[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
            
    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      other[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
            
    }
  }// end ossf
  // No_OSSF
  else if (   ((ch_lepton1 + ch_lepton2) == 0  )  || ((ch_lepton1 + ch_lepton3) == 0   ) || ((ch_lepton3 + ch_lepton2) == 0   )   ){
    //cout<<"in function where kind is 0: "<<kind[0]<<endl;
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
      other[0].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
    }
    if (index_inv == 1) {
      pair[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      other[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
            
    }
    if (index_inv == 2) {
      pair[0].SetPtEtaPhiE( leep2.Pt(),  leep2.Eta(), leep2.Phi(), leep2.E());
      pair[1].SetPtEtaPhiE( leep3.Pt(),  leep3.Eta(), leep3.Phi(), leep3.E());
      other[0].SetPtEtaPhiE( leep1.Pt(),  leep1.Eta(), leep1.Phi(), leep1.E());
            
    }
        

  }
}
//______________________________________________
double Analysis_mc::FR_weight(TGraphAsymmErrors *fakeRate_mu_sFR[3],
			      TGraphAsymmErrors *fakeRate_e_sFR[3],
			      TGraphAsymmErrors *fakeRate_mumu_dFR[3],
			      TGraphAsymmErrors *fakeRate_ee_dFR[3],
			      TGraphAsymmErrors *fakeRate_emu_dFR[3],
			      bool   isSFR,
			      bool   isDFR,
			      double etaLepton,
			      double flavorsLepton,
			      double ptLepton,
			      double etaJet,
			      double flavorsJet,
			      double ptJet
			      ) {
  double factor=0;

  //std::cout<< "                   inside function for FR: "<<factor<<std::endl;
  if (isSFR && flavorsLepton == 0 )  factor =  sFR_factor_e  (*&fakeRate_e_sFR ,etaLepton, ptLepton);
  if (isSFR && flavorsLepton == 1 )  factor =  sFR_factor_mu (*&fakeRate_mu_sFR ,etaLepton, ptLepton);

  if (isDFR && flavorsJet == 0 )     factor =  dFR_factor_ee   (*&fakeRate_ee_dFR ,etaJet, ptJet);
  if (isDFR && flavorsJet == 1 )     factor =  dFR_factor_mumu (*&fakeRate_mumu_dFR ,etaJet, ptJet);
  if (isDFR && flavorsJet == 2 )     factor =  dFR_factor_emu  (*&fakeRate_emu_dFR ,etaJet, ptJet);
  //std::cout<< "                   return factor: "<<factor<<std::endl;

  return factor;
}


//==================================================================
double Analysis_mc::dFR_factor_ee(TGraphAsymmErrors *fakeRate_e[3],                             
				  int eta,
				  double lptcone
				  ){
    
  const int nBinMu=4;
  const int nBinMu3=4;
    
  Double_t newBins_e1[nBinMu+1] = {10,25,35,55, 80};
  Double_t newBins_e3[nBinMu3+1] ={10,25,35,55, 80};

  //std::cout<<"                             it is ee: eta and pt  "<<eta<< "   "<<lptcone<< std::endl;
   
  TH1D *fakeRate_e_histo[3]; 
  fakeRate_e_histo[0]= new TH1D("fake_rate_e_histo_eta1","",nBinMu,newBins_e1);
  fakeRate_e_histo[1]= new TH1D("fake_rate_e_histo_eta2","",nBinMu,newBins_e1);
  fakeRate_e_histo[2]= new TH1D("fake_rate_e_histo_eta3","",nBinMu3,newBins_e3);  
  for (int i=0; i< 3; i++){
    if (i ==0 || i ==1) from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],nBinMu);
    if (i ==2) from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],nBinMu3);
  }
  //std::cout<<"                             check first bin eta1: "<< fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(20))<<std::endl;

  double momentum = lptcone;
  double factore=0;
  if (momentum < 80){
    if (eta == 1)  factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(momentum));
    if (eta == 2)  factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(momentum));
    if (eta == 3)  factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(momentum));
  }//eta1
  else {
    if (eta == 1)  factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(70));
    if (eta == 2)  factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(70));
    if (eta == 3)  factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(70));
  }  
  delete fakeRate_e_histo[0];
  delete fakeRate_e_histo[1];
  delete fakeRate_e_histo[2];
  //std::cout<<"                             it is ee: factor "<<factore<<std::endl;

  return factore;  
}

//==================================================================
double Analysis_mc::dFR_factor_emu(TGraphAsymmErrors *fakeRate_e[3],                             
				   int eta,
				   double lptcone
				   ){
    
  const int nBinMu=4;
  const int nBinMu3=4;
    
  Double_t newBins_e1[nBinMu+1] = {10,25,35,55, 80};
  Double_t newBins_e3[nBinMu3+1] = {10,25,35,55, 80};
  // std::cout<<"                             it is emu: eta and pt  "<<eta<< "   "<<lptcone<< std::endl;

  TH1D *fakeRate_e_histo[3]; 
  fakeRate_e_histo[0]= new TH1D("fake_rate_e_histo_eta1","",nBinMu,newBins_e1);
  fakeRate_e_histo[1]= new TH1D("fake_rate_e_histo_eta2","",nBinMu,newBins_e1);
  fakeRate_e_histo[2]= new TH1D("fake_rate_e_histo_eta3","",nBinMu3,newBins_e3);  
  for (int i=0; i< 3; i++){
    if (i ==0 || i ==1) from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],nBinMu);
    if (i ==2) from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],nBinMu3);
  }
  //std::cout<<"                             check first bin eta1: "<< fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(20))<<std::endl;

  double momentum = lptcone;
  double factore=0;
  if (momentum < 80){
    if (eta == 1)  factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(momentum));
    if (eta == 2)  factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(momentum));
    if (eta == 3)  factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(momentum));
  }//eta1
  else {
    if (eta == 1)  factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(70));
    if (eta == 2)  factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(70));
    if (eta == 3)  factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(70));
  }  
  delete fakeRate_e_histo[0];
  delete fakeRate_e_histo[1];
  delete fakeRate_e_histo[2];
  //std::cout<<"                             it is emu: factor "<<factore<<std::endl;

  return factore;  
}

//==================================================================
double Analysis_mc::dFR_factor_mumu(TGraphAsymmErrors *fakeRate_e[3],                             
				    int eta,
				    double lptcone
				    ){
    
  const int nBinMu=4;
  //const int nBinMu3=4;
  Double_t newBins_e1[nBinMu+1] = {10,25,35,55, 80};
  //  std::cout<<"                             it is mumu: eta and pt  "<<eta<< "   "<<lptcone<< std::endl;

  TH1D *fakeRate_e_histo[3];  
  fakeRate_e_histo[0]= new TH1D("fake_rate_e_histo_eta1","",nBinMu,newBins_e1);
  fakeRate_e_histo[1]= new TH1D("fake_rate_e_histo_eta2","",nBinMu,newBins_e1);
  fakeRate_e_histo[2]= new TH1D("fake_rate_e_histo_eta3","",nBinMu,newBins_e1);  
  for (int i=0; i< 3; i++){
    from_TGraph_to_TH1D(*&fakeRate_e[i],*&fakeRate_e_histo[i],nBinMu);
  }

  //std::cout<<"                             check first bin eta1: "<< fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(20))<<std::endl;
  
  double momentum = lptcone;
  double factore=0;
  if (momentum < 80){
    if (eta == 1)  factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(momentum));
    if (eta == 2)  factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(momentum));
    if (eta == 3)  factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(momentum));
  }//eta1
  else {
    if (eta == 1)  factore = fakeRate_e_histo[0]->GetBinContent(fakeRate_e_histo[0]->FindBin(70));
    if (eta == 2)  factore = fakeRate_e_histo[1]->GetBinContent(fakeRate_e_histo[1]->FindBin(70));
    if (eta == 3)  factore = fakeRate_e_histo[2]->GetBinContent(fakeRate_e_histo[2]->FindBin(70));
  }  
  delete fakeRate_e_histo[0];
  delete fakeRate_e_histo[1];
  delete fakeRate_e_histo[2];
  //   std::cout<<"                             it is mumu: factor "<<factore<<std::endl;

  return factore;  
}


//______________________________________________
double Analysis_mc::sFR_factor_e (TGraphAsymmErrors *fakeRate[3],                             
				  double eta,                          
				  double lptcone
				  ){
  eta = fabs(eta);
  TH1D *fakeRate_histo[3];
  Double_t newBins[6] = {10, 15, 25, 35, 50, 70};
  fakeRate_histo[0]= new TH1D("fake_rate_e_histo_eta1","",5,newBins);
  fakeRate_histo[1]= new TH1D("fake_rate_e_histo_eta2","",5,newBins);
  fakeRate_histo[2]= new TH1D("fake_rate_e_histo_eta3","",5,newBins);
  for (int i=0; i< 3; i++){
    from_TGraph_to_TH1D(*&fakeRate[i],*&fakeRate_histo[i],5);
  }
  double momentum = lptcone;
  if (lptcone < 10 ) momentum = 11.;
  double factore=0;
  if (momentum < 70){
    if (eta < 0.8){
      factore = fakeRate_histo[0]->GetBinContent(fakeRate_histo[0]->FindBin(momentum));
    }//eta1
    else if (eta > 0.8 && eta<1.479){
      factore = fakeRate_histo[1]->GetBinContent(fakeRate_histo[1]->FindBin(momentum));
    }//eta1
    else {
      factore = fakeRate_histo[2]->GetBinContent(fakeRate_histo[2]->FindBin(momentum));
    }//eta1
  }// <70
  else {
    if (eta < 0.8){
      factore = fakeRate_histo[0]->GetBinContent(fakeRate_histo[0]->FindBin(60));
    }//eta1
    else if (eta > 0.8 && eta<1.479){
      factore = fakeRate_histo[1]->GetBinContent(fakeRate_histo[1]->FindBin(60));
    }//eta1
    else {
      factore = fakeRate_histo[2]->GetBinContent(fakeRate_histo[2]->FindBin(60));
    }//eta1
  }
 
  delete fakeRate_histo[0];
  delete fakeRate_histo[1];
  delete fakeRate_histo[2];
  return factore;
}


//______________________________________________
double Analysis_mc::sFR_factor_mu (TGraphAsymmErrors *fakeRate[3],                             
				   double eta,                          
				   double lptcone
				   ){
  eta = fabs(eta);
  TH1D *fakeRate_histo[3];
  Double_t newBins[7] = {5,10, 15, 25, 35, 50, 70};
  fakeRate_histo[0]= new TH1D("fake_rate_e_histo_eta1","",6,newBins);
  fakeRate_histo[1]= new TH1D("fake_rate_e_histo_eta2","",6,newBins);
  fakeRate_histo[2]= new TH1D("fake_rate_e_histo_eta3","",6,newBins);
  for (int i=0; i< 3; i++){
    from_TGraph_to_TH1D(*&fakeRate[i],*&fakeRate_histo[i],6);
  }
  double momentum = lptcone;
  double factore=0;
  if (momentum < 70){
    if (eta < 0.8){
      factore = fakeRate_histo[0]->GetBinContent(fakeRate_histo[0]->FindBin(momentum));
    }//eta1
    else if (eta > 0.8 && eta<1.479){
      factore = fakeRate_histo[1]->GetBinContent(fakeRate_histo[1]->FindBin(momentum));
    }//eta1
    else {
      factore = fakeRate_histo[2]->GetBinContent(fakeRate_histo[2]->FindBin(momentum));
    }//eta1
  }// <70
  else {
    if (eta < 0.8){
      factore = fakeRate_histo[0]->GetBinContent(fakeRate_histo[0]->FindBin(60));
    }//eta1
    else if (eta > 0.8 && eta<1.479){
      factore = fakeRate_histo[1]->GetBinContent(fakeRate_histo[1]->FindBin(60));
    }//eta1
    else {
      factore = fakeRate_histo[2]->GetBinContent(fakeRate_histo[2]->FindBin(60));
    }//eta1
  }
 
  delete fakeRate_histo[0];
  delete fakeRate_histo[1];
  delete fakeRate_histo[2];
  return factore;
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
    //std::cout<<"                                                          "<<i<<") "<<x_graph[i]<<"   "<<y_graph[i]<<std::endl;
    histo->SetBinContent (i+1, x_graph[i],  y_graph[i]);        
  }
}




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


//_______________________________________________________ constructor_____
void Analysis_mc::put_at_zero(int anno, int iSystematics,int iVariation,int channel, int option, TH1D *histo){
  //option 1 --> non prompt
  //option 0 --> MC

  double ee_sf =0.;
  double em_sf =0.;
  double mm_sf =0.;
  if (anno == 0){
    ee_sf =0.38;
    em_sf =0.20;
    mm_sf =0.21;
  }
  if (anno == 1){
    ee_sf =0.59;
    em_sf =0.33;
    mm_sf =0.26;
  }
  if (anno == 2){
    ee_sf =0.47;
    em_sf =0.31;
    mm_sf =0.26;
  }
  ee_sf = ee_sf / (1-ee_sf);
  em_sf = em_sf / (1-em_sf);
  mm_sf = mm_sf / (1-mm_sf);

  if (option == 0){
    for (int i =0; i < histo-> GetNbinsX(); i++){
      double error_original =0;
      double error_to_add =0;
      double error_final =0;
      if (histo->GetBinContent( i+1)  < 0.  || std::isnan(histo->GetBinContent( i+1)) || histo->GetBinContent( i+1)  <0.00001) {
	error_original = histo-> GetBinError(i+1);
	error_to_add = histo-> GetBinContent(i+1);
	error_final=TMath::Sqrt(error_original*error_original   +    error_to_add*error_to_add );
	histo-> SetBinContent(i+1, 0.000001);
	histo-> SetBinError(i+1, 0.000001);
      }
    }
  }//option 0

  if (option == 1){
    for (int i =0; i < histo-> GetNbinsX(); i++){
      if (histo->GetBinContent( i+1)  <= 0 ){
	histo-> SetBinContent(i+1, 0.0);
	if (channel == 0 && i < 6) histo->SetBinError(i+1, mm_sf); //mumu
	if (channel == 0 && i >= 6) histo->SetBinError(i+1, em_sf); //mue
	if (channel == 1 && i < 6) histo->SetBinError(i+1, ee_sf); //ee
	if (channel == 1 && i >= 6) histo->SetBinError(i+1, em_sf); //mue
      }     
    }
  }//option 1

  if (option == 2){
    for (int i =0; i < histo-> GetNbinsX(); i++){
      if (histo->GetBinContent( i+1)  == 0.0 ){
	histo-> SetBinContent(i+1, 0.0);
	histo->SetBinError(i+1, 0.40/(1-0.40)); //mumu
      }
    }
  }//option 2
}








