#include <vector>
#include <iostream> 
#include <map>
#include <cmath>
#include "pipIDsubroutines.C"

std::map<int,std::vector<double> > hadronMLEID(int hypothesis, int gpart, int e_index[], Int_t q[], Float_t p[], UChar_t sc_sect[], UChar_t dc_sect[], Float_t sc_t[], Float_t sc_r[], UChar_t sc_pd[], int pip_vvp_strict, int pip_R1fid_strict, int pip_MXcut_strict, int ExpOrSim, Float_t ec_ei[], UChar_t ec_sect[], UChar_t cc_sect[], UShort_t nphe[], Float_t ec_eo[], Float_t cx[], Float_t cy[], Float_t cz[], Float_t b[], Float_t tl1_x[], Float_t tl1_y[],  TLorentzVector V4_H[], int currentrunno, int pim_vvp_strict, int pim_R1fid_strict, int pim_MXcut_strict, std::map<int, std::vector<double> > pr_mean_fit, std::map<int, std::vector<double> > pip_mean_fit, std::map<int, std::vector<double> > kp_mean_fit, std::map<int, std::vector<double> > prot_sigma_fit, std::map<int, std::vector<double> > pip_sigma_fit, std::map<int, std::vector<double> > kp_sigma_fit, double pr_conf, double pip_conf, double kp_conf, double pr_anticonf, double pip_anticonf, double kp_anticonf )
{
  Float_t speed_of_light = 29.9792458; // cm/ns
  Float_t pip_mass = 0.13957018; // GeV
  Float_t pim_mass = 0.13957018; // GeV
  Float_t prot_mass = 0.93827203; // GeV
  Float_t kp_mass = 0.49367; // GeV
  Float_t km_mass = 0.49367; // GeV

  Float_t pi = 3.14159265359;
  Float_t pi180 = pi/180.0;

  std::vector<int> hadrons;
  double temp_p_pr = -1.0;
  double temp_p_pip = -1.0;
  double temp_p_kp = -1.0;

  double temp_index_pr = -1;
  double temp_index_pip = -1;
  double temp_index_kp = -1;

  double temp_conf_pr = -1.0;
  double temp_conf_pip = -1.0;
  double temp_conf_kp = -1.0;

  std::cout << "[hadronMLEID] GOING TO LOOP OVER GPART IN HADRONS " << std::endl;
  std::cout << "[hadronMLEID] GPART " << gpart << std::endl;
  for(int k = 1; k < gpart; k++) // loop over particles
    {
      float cand_time = -1.0;
      float beta_measured_pim = -1.0;
      float beta_measured_km = -1.0;
	
      float beta_measured_pr = -1.0;
      float beta_measured_kp = -1.0;
      float beta_measured_pip = -1.0;
	
      float beta_measured = -1.0;

      int sc_sector = sc_sect[k]; //from 0 to 6
      int dc_sector = dc_sect[k]; //from 0 to 6
      
      std::cout << "[hadronMLEID] SECTOR " << dc_sector-1 << std::endl;
      float corrStartTime = e_sctimeCorr(ExpOrSim, sc_t[e_index[1]], sc_sect[e_index[1]], sc_pd[e_index[1]], currentrunno) - (sc_r[e_index[1]]/speed_of_light);
      float ehtcorrMeasuredTime = h_sctimeCorr(ExpOrSim, sc_t[k], sc_sect[k], sc_pd[k], currentrunno) - corrStartTime;
      beta_measured = (sc_r[k]/ehtcorrMeasuredTime)/speed_of_light;
      
      if( q[k] > 0 && k != e_index[1] && sc_sect[k] != 0 && dc_sect[k] != 0 && goodORbadSCpaddle(dc_sect[k], sc_pd[k]) == 1 ){
	//std::cout << " sector " << sc_sector << std::endl;

	// using Nathans pip fiducial cut routine for positive particles.
	// since it is not specific for strictly pion plus. the pim routine
	// can be used for negative hadron pid
	bool R1fid_pass = pip_R1fid_pass(pip_R1fid_strict, dc_sect[k], tl1_x[k], tl1_y[k]);	

	//double hadron_corr_vz  = getCorrZ(ExpOrSim, vx[k], vy[k], vz[k], p[k]*cx[k], p[k]*cy[k], p[k]*cz[k], dc_sect[k]);
	//double el_corr_vz = getCorrZ(ExpOrSim, vx[e_index[1]], vy[e_index[1]], vz[e_index[1]], p[e_index[1]]*cx[e_index[1]], p[e_index[1]]*cy[e_index[1]], p[e_index[1]]*cz[e_index[1]], dc_sect[e_index[1]]);
	//double delta_vz = el_corr_vz - hadron_corr_vz;
	//bool delta_vertex_cut = std::abs( el_corr_vz - hadron_corr_vz ) < 5.0; 


	if( R1fid_pass ){

	  //1/sqrt(1+pow(mass/p[index],2)); 
	  float CalcBeta_pip = 1./sqrt(1.0 + pow(pip_mass/p[k],2));
	  float CalcBeta_pim = 1./sqrt(1.0 + pow(pim_mass/p[k],2));
	  float CalcBeta_kp = 1./sqrt(1.0 + pow(kp_mass/p[k],2));
	  float CalcBeta_km = 1./sqrt(1.0 + pow(km_mass/p[k],2));
	  double CalcBeta_pr = 1./sqrt(1.0 + pow(prot_mass/p[k],2));

	  std::vector<float> calc_beta_pos;
	  calc_beta_pos.push_back(CalcBeta_pip);
	  calc_beta_pos.push_back(CalcBeta_kp);
	  calc_beta_pos.push_back(CalcBeta_pr);
	
	  std::vector<float> calc_beta_neg;
	  calc_beta_neg.push_back(CalcBeta_pim);
	  calc_beta_neg.push_back(CalcBeta_km);

	  //mle fit values
	  //std::cout << " >> GETTING FIT VALUES " << std::endl;
	  //std:: cout << " >> PROTON MEAN MAP SIZE"  << pr_mean_fit.size()  << std::endl;
	  std::vector<double> prot_mean_a = pr_mean_fit[2];//->second;      
	  std::vector<double> prot_mean_b = pr_mean_fit[1];//)->second;      
	  std::vector<double> prot_mean_c = pr_mean_fit[0];//second;      

	  //std:: cout << " >> PIP MEAN MAP SIZE"  << pip_mean_fit.size()  << std::endl;	
	  std::vector<double> pip_mean_a = pip_mean_fit[2];//->second;      
	  std::vector<double> pip_mean_b = pip_mean_fit[1];//->second;      
	  std::vector<double> pip_mean_c = pip_mean_fit[0];//->second;      

	  std::vector<double> Kp_mean_a = kp_mean_fit[2];//->second;      
	  std::vector<double> Kp_mean_b = kp_mean_fit[1];//->second;      
	  std::vector<double> Kp_mean_c = kp_mean_fit[0];//->second;      

	  //std::cout << " proton sigma size " << prot_sigma_fit.size() << std::endl;
	  std::vector<double> prot_sigma_a = prot_sigma_fit[2];//->second;      
	  std::vector<double> prot_sigma_b = prot_sigma_fit[1];//.find(1)->second;      
	  std::vector<double> prot_sigma_c = prot_sigma_fit[0];//.find(2)->second;      
	
	  std::vector<double> pip_sigma_a = pip_sigma_fit[2];//.find(0)->second;      
	  std::vector<double> pip_sigma_b = pip_sigma_fit[1];//.find(1)->second;      
	  std::vector<double> pip_sigma_c = pip_sigma_fit[0];//.find(2)->second;      

	  std::vector<double> Kp_sigma_a = kp_sigma_fit[2];//.find(0)->second;      
	  std::vector<double> Kp_sigma_b = kp_sigma_fit[1];//.find(1)->second;      
	  std::vector<double> Kp_sigma_c = kp_sigma_fit[0];//.find(2)->second;      

	  double mom = p[k];

	  //std::cout << " PRINTING SIZE " << std::endl;
	  //std::cout << " SIZE " << prot_mean_a.size() << " " << prot_mean_a[dc_sector-1] << " " << prot_mean_b[dc_sector-1] << " " << prot_mean_c[dc_sector-1] << std::endl;

	  //	double mean_prot = prot_mean[0] + prot_mean[1]*mom/sqrt(mom*mom+prot_mean[2]) + prot_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+prot_mean[4]);
	  //	double sigma_prot = prot_sigma[0] + prot_sigma[1]/sqrt(mom) + prot_sigma[2]*exp(-mom);
	  //using delta beta in here
	  double mean_prot = CalcBeta_pr + prot_mean_c[dc_sector-1]*pow(mom,2) + prot_mean_b[dc_sector-1]*mom + prot_mean_a[dc_sector-1];
	  double sigma_prot = prot_sigma_c[dc_sector-1]*pow(mom,2) + prot_sigma_b[dc_sector-1]*mom + prot_sigma_a[dc_sector-1];
	  double deltab_prot = beta_measured - mean_prot;// CalcBeta_pr - mean_prot;

	  std::cout << " [hadronMLE] PROTON " << std::endl;

	  std::cout << " sector " << dc_sector-1 << " parameter 0 " << prot_mean_a[dc_sector-1] << std::endl;
	  std::cout << " sector " << dc_sector-1 << " parameter 1 " << prot_mean_b[dc_sector-1] << std::endl;
	  std::cout << " sector " << dc_sector-1 << " parameter 2 " << prot_mean_c[dc_sector-1] << std::endl;

	  std::cout << " theory " << CalcBeta_pr << std::endl;
	  std::cout << " measured " << beta_measured << std::endl;	 
	  std::cout << " mean " << mean_prot << std::endl;
	  std::cout << " sig " << sigma_prot << std::endl;
	  std::cout << " residual " << deltab_prot << std::endl;
	  double norm = sqrt(1./2.0/3.14159);
	  double prob_prot = norm/(sigma_prot) * exp(-0.5 * pow((deltab_prot)/sigma_prot, 2));
	  double conf_prot = (1.0 - TMath::Erf(fabs(deltab_prot)/sigma_prot/sqrt(2.0)));  // removed beta_measures - mean_prot

	  //double mean_pip = pip_mean[0] + pip_mean[1]*mom/sqrt(mom*mom+pip_mean[2]) + pip_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+pip_mean[4]);
	  //double sigma_pip = pip_sigma[0] + pip_sigma[1]/sqrt(mom) + pip_sigma[2]*exp(-mom);
	  double mean_pip =  CalcBeta_pip + pip_mean_c[dc_sector-1]*pow(mom,2) + pip_mean_b[dc_sector-1]*mom + pip_mean_a[dc_sector-1];
	  double sigma_pip =  pip_sigma_c[dc_sector-1]*pow(mom,2) + pip_sigma_b[dc_sector-1]*mom + pip_sigma_a[dc_sector-1];
	  double deltab_pip = beta_measured - mean_pip;//CalcBeta_pip - mean_pip;

	  std::cout << " [hadronMLE] PION PLUS " << std::endl;
	  std::cout << " theory " << CalcBeta_pip << std::endl;
	  std::cout << " mean " << mean_pip << std::endl;
	  std::cout << " sig " << sigma_pip << std::endl;
	  std::cout << " residual " << deltab_pip << std::endl;
	  double prob_pip = (norm/(sigma_pip) * exp(-0.5 * pow((deltab_pip)/sigma_pip, 2)));
	  double conf_pip = (1.0 - TMath::Erf(fabs(deltab_pip)/sigma_pip/sqrt(2.0))); 
	  
	  //double mean_Kp = Kp_mean[0] + Kp_mean[1]*mom/sqrt(mom*mom+Kp_mean[2]) + Kp_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+Kp_mean[4]);
	  //double sigma_Kp = Kp_sigma[0] + Kp_sigma[1]/sqrt(mom) + Kp_sigma[2]*exp(-mom);
	  double mean_Kp = CalcBeta_kp + Kp_mean_c[dc_sector-1]*pow(mom,2) + Kp_mean_b[dc_sector-1]*mom + Kp_mean_a[dc_sector-1];
	  double sigma_Kp =  Kp_sigma_c[dc_sector-1]*pow(mom,2) + Kp_sigma_b[dc_sector-1]*mom + Kp_sigma_a[dc_sector-1];
	  double deltab_kp = beta_measured - mean_Kp; //CalcBeta_kp - mean_Kp;		  

	  std::cout << " [hadronMLE] Kaon Plus" << std::endl;
	  std::cout << " theory " << CalcBeta_kp << std::endl;
	  std::cout << " mean " << mean_Kp << std::endl;
	  std::cout << " sig " << sigma_Kp << std::endl;
	  std::cout << " residual " << deltab_kp << std::endl;
	  double prob_Kp = norm/(sigma_Kp) * exp(-0.5 * pow((deltab_kp)/sigma_Kp, 2));
	  double conf_Kp = (1.0 - TMath::Erf(fabs(deltab_kp)/sigma_Kp/sqrt(2.0)));        
			    
	  std::cout << "[hadronMLEID] PROB " << prob_prot << " " << prob_pip << " " << prob_Kp << std::endl;
	  std::cout << "[hadronMLEID] CONF " << conf_prot << " " << conf_pip << " " << conf_Kp << std::endl;
	
	  if(prob_prot > prob_pip  && prob_prot > prob_Kp && hypothesis == 2212 ){
	    //if( hypothesis == 2212 ){
	    //std::cout << " proton " << std::endl;
	    // commenting out because I want to keep all proton candidates to cut on later
	     //&& conf_prot > pr_conf  && conf_pip < pr_anticonf  && conf_Kp < pr_anticonf){ 
	  
	    /*cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << CalcBeta_pip << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	    cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	    cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	    */ //cout << endl;
	  

	    //keep only the largest momentum candidates
	    std::cout << " momentum of proton " << p[k] << std::endl;
	    if( p[k] > temp_p_pr ){
	      temp_p_pr = p[k];
	      temp_index_pr = k;
	      temp_conf_pr = conf_prot;	      
	    }
	    std::cout << " max proton momentum " << temp_p_pr << std::endl;
	    std::cout << " good proton index " << temp_index_pr << " with confidence level " << temp_conf_pr << std::endl;
	    //return true;
	  }
	  else if(prob_pip > prob_prot && prob_pip  > prob_Kp  && hypothesis == 211 ){
	    //else if(hypothesis == 211 ){
	    //&& conf_pip > pip_conf && conf_prot < pip_anticonf  && conf_Kp  < pip_anticonf){ 
	    //std::cout << " pion plus " << std::endl;
	    /*cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	    cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	     cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	    */ //cout << endl;

	    std::cout << " momentum of pion plus " << p[k] << std::endl;
	    std::cout << " max pion plus momentum " << temp_p_pip << std::endl;
	    if( p[k] > temp_p_pip ){
	      temp_p_pip = p[k];
	      temp_index_pip = k;
	      temp_conf_pip = conf_pip;
	    }
	    std::cout << " max pion plus momentum " << temp_p_pip << std::endl;
	    std::cout << " good pion plus index " << temp_index_pip << " with confidence level " << temp_conf_pip << std::endl;
	    //return true;
	  }
	  else if(prob_Kp > prob_prot && prob_Kp > prob_pip && hypothesis == 321 ){
	    //else if (hypothesis == 321 ){
	    //std::cout << " KAON plus " << std::endl;
	    // && conf_Kp > kp_conf  && conf_prot < pr_anticonf  && conf_pip < pr_anticonf ){ 
	    /* cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	    cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	    cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	    */ //cout << endl;

	    std::cout << " momentum of kaon plus " << p[k] << std::endl;
	    if( p[k] > temp_p_kp ){
	      temp_p_kp = p[k];
	      temp_index_kp = k;
	      temp_conf_kp = conf_Kp;
	    }
	    std::cout << " max kaon plus " << temp_p_kp << std::endl;
	    std::cout << " good kaon plus index " << temp_index_kp << " with confidence level " << temp_conf_kp << std::endl;

	    //return true;
	  }
	}
      }
    }

  std::cout << "[hadronMLE] final variables for hadron " << std::endl;

  std::cout << "[hadronMLE] proton " << std::endl;
  std::cout << " >> index " << temp_index_pr << " conf " << temp_conf_pr << std::endl;
  
  std::map<int, std::vector<double> > hadron_candidates;
  std::vector<double> proton_classifier_values;
  proton_classifier_values.push_back(temp_index_pr);
  proton_classifier_values.push_back(temp_conf_pr);

  std::cout << "[hadronMLE] pion plus " << std::endl;
  std::cout << " >> index " << temp_index_pip << " conf " << temp_conf_pip << std::endl;

  std::vector<double> pip_classifier_values;
  pip_classifier_values.push_back(temp_index_pip);
  pip_classifier_values.push_back(temp_conf_pip);

  std::cout << "[hadronMLE] kaon plus " << std::endl;
  std::cout << " >> index " << temp_index_kp << " conf " << temp_conf_kp << std::endl;

  std::vector<double> kaonP_classifier_values;
  kaonP_classifier_values.push_back(temp_index_kp);
  kaonP_classifier_values.push_back(temp_conf_kp);

  hadron_candidates[2212] = proton_classifier_values;
  hadron_candidates[211] = pip_classifier_values;
  hadron_candidates[321] = kaonP_classifier_values;

  return hadron_candidates;

}

