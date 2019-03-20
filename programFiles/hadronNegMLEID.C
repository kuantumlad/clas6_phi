#include <vector>
#include <iostream> 
#include <map>
#include <cmath>
//#include "pipIDsubroutines.C"

std::map<int,std::vector<double> > hadronNegMLEID(int hypothesis, int gpart, int e_index[], Int_t q[], Float_t p[], UChar_t sc_sect[], UChar_t dc_sect[], Float_t sc_t[], Float_t sc_r[], UChar_t sc_pd[], int pip_vvp_strict, int pip_R1fid_strict, int pip_MXcut_strict, int ExpOrSim, Float_t ec_ei[], UChar_t ec_sect[], UChar_t cc_sect[], UShort_t nphe[], Float_t ec_eo[], Float_t cx[], Float_t cy[], Float_t cz[], Float_t b[], Float_t tl1_x[], Float_t tl1_y[],  TLorentzVector V4_H[], int currentrunno, int pim_vvp_strict, int pim_R1fid_strict, int pim_MXcut_strict, std::map<int, std::vector<double> > pim_mean_fit, std::map<int, std::vector<double> > km_mean_fit, std::map<int, std::vector<double> > pim_sigma_fit, std::map<int, std::vector<double> > km_sigma_fit, double pim_conf, double km_conf, double pim_anticonf, double km_anticonf )
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
  double temp_p_pim = -1.0;
  double temp_p_km = -1.0;

  double temp_index_pim = -1;
  double temp_index_km = -1;

  double temp_conf_pim = -1.0;
  double temp_conf_km = -1.0;

  std::cout << " [HadronNegMLEDID] GOING TO LOOP OVER GPART IN HADRONS " << std::endl;
  std::cout << " GPART " << gpart << std::endl;
  for(int k = 0; k < gpart; k++) // loop over particles
    {
      float cand_time = -1.0;
      float beta_measured_pim = -1.0;
      float beta_measured_km = -1.0;
		
      float beta_measured = -1.0;

      int sc_sector = sc_sect[k]; //from 0 to 6
      std::cout << " SECTOR " << sc_sector << std::endl;
      float corrStartTime = e_sctimeCorr(ExpOrSim, sc_t[e_index[1]], sc_sect[e_index[1]], sc_pd[e_index[1]], currentrunno) - (sc_r[e_index[1]]/speed_of_light);
      float ehtcorrMeasuredTime = h_sctimeCorr(ExpOrSim, sc_t[k], sc_sect[k], sc_pd[k], currentrunno) - corrStartTime;
      beta_measured = (sc_r[k]/ehtcorrMeasuredTime)/speed_of_light;
      
      if( sc_sector >= 1 && q[k] < 0 && k != e_index[1] && sc_sect[k] != 0 && dc_sect[k] != 0 && goodORbadSCpaddle(dc_sect[k], sc_pd[k]) == 1 ){
	//std::cout << " sector " << sc_sector << std::endl;

	// using Nathans pip fiducial cut routine for positive particles.
	// since it is not specific for strictly pion plus. the pim routine
	// can be used for negative hadron pid
	bool R1fid_pass = pip_R1fid_pass(pip_R1fid_strict, dc_sect[k], tl1_x[k], tl1_y[k]);	
	if( R1fid_pass ){

	  float CalcBeta_pim = 1./sqrt(1.0 + pow(pim_mass/p[k],2));
	  float CalcBeta_km = 1./sqrt(1.0 + pow(km_mass/p[k],2));
	
	  std::vector<float> calc_beta_neg;
	  calc_beta_neg.push_back(CalcBeta_pim);
	  calc_beta_neg.push_back(CalcBeta_km);

	  //mle fit values
	  //std:: cout << " >> PIP MEAN MAP SIZE"  << pip_mean_fit.size()  << std::endl;	
	  std::vector<double> pim_mean_a = pim_mean_fit[0];//->second;      
	  std::vector<double> pim_mean_b = pim_mean_fit[1];//->second;      
	  std::vector<double> pim_mean_c = pim_mean_fit[2];//->second;      

	  std::vector<double> Km_mean_a = km_mean_fit[0];//->second;      
	  std::vector<double> Km_mean_b = km_mean_fit[1];//->second;      
	  std::vector<double> Km_mean_c = km_mean_fit[2];//->second;      

	  std::vector<double> pim_sigma_a = pim_sigma_fit[0];//.find(0)->second;      
	  std::vector<double> pim_sigma_b = pim_sigma_fit[1];//.find(1)->second;      
	  std::vector<double> pim_sigma_c = pim_sigma_fit[2];//.find(2)->second;      

	  std::vector<double> Km_sigma_a = km_sigma_fit[0];//.find(0)->second;      
	  std::vector<double> Km_sigma_b = km_sigma_fit[1];//.find(1)->second;      
	  std::vector<double> Km_sigma_c = km_sigma_fit[2];//.find(2)->second;      

	  double mom = p[k];

	  //std::cout << " PRINTING SIZE " << std::endl;
	  //std::cout << " SIZE " << prot_mean_a.size() << " " << prot_mean_a[sc_sector-1] << " " << prot_mean_b[sc_sector-1] << " " << prot_mean_c[sc_sector-1] << std::endl;

	  //	double mean_prot = prot_mean[0] + prot_mean[1]*mom/sqrt(mom*mom+prot_mean[2]) + prot_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+prot_mean[4]);
	  //	double sigma_prot = prot_sigma[0] + prot_sigma[1]/sqrt(mom) + prot_sigma[2]*exp(-mom);
	  //using delta beta in here
	  //double deltab_prot = beta_measured - CalcBeta_pr;
	  //double mean_prot = prot_mean_c[sc_sector-1]*pow(mom,2) + prot_mean_b[sc_sector-1]*mom + prot_mean_a[sc_sector-1];
	  //double sigma_prot = prot_sigma_c[sc_sector-1]*pow(mom,2) + prot_sigma_b[sc_sector-1]*mom + prot_sigma_a[sc_sector-1];
	  //double prob_prot = (1/(sigma_prot*sqrt(2*3.14159))) * exp(-0.5 * pow((deltab_prot + mean_prot)/sigma_prot, 2));
	  //double conf_prot = (1.0 - TMath::Erf(fabs(deltab_prot + mean_prot)/sigma_prot/sqrt(2.0)));  // removed beta_measures - mean_prot

	  //double mean_pip = pip_mean[0] + pip_mean[1]*mom/sqrt(mom*mom+pip_mean[2]) + pip_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+pip_mean[4]);
	  //double sigma_pip = pip_sigma[0] + pip_sigma[1]/sqrt(mom) + pip_sigma[2]*exp(-mom);
	  double deltab_pim = beta_measured - CalcBeta_pim;
	  double mean_pim =  pim_mean_c[sc_sector-1]*pow(mom,2) + pim_mean_b[sc_sector-1]*mom + pim_mean_a[sc_sector-1];
	  double sigma_pim =  pim_sigma_c[sc_sector-1]*pow(mom,2) + pim_sigma_b[sc_sector-1]*mom + pim_sigma_a[sc_sector-1];
	  double prob_pim = (1/(sigma_pim*sqrt(2*3.14159))) * exp(-0.5 * pow((deltab_pim + mean_pim)/sigma_pim, 2));
	  double conf_pim = (1.0 - TMath::Erf(fabs(deltab_pim + mean_pim)/sigma_pim/sqrt(2.0))); 
	  
	  //double mean_Kp = Kp_mean[0] + Kp_mean[1]*mom/sqrt(mom*mom+Kp_mean[2]) + Kp_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+Kp_mean[4]);
	  //double sigma_Kp = Kp_sigma[0] + Kp_sigma[1]/sqrt(mom) + Kp_sigma[2]*exp(-mom);
	  double deltab_km = beta_measured - CalcBeta_km;		  
	  double mean_Km =  Km_mean_c[sc_sector-1]*pow(mom,2) + Km_mean_b[sc_sector-1]*mom + Km_mean_a[sc_sector-1];
	  double sigma_Km =  Km_sigma_c[sc_sector-1]*pow(mom,2) + Km_sigma_b[sc_sector-1]*mom + Km_sigma_a[sc_sector-1];
	  double prob_Km = (1/(sigma_Km*sqrt(2*3.14159))) * exp(-0.5 * pow((deltab_km + mean_Km)/sigma_Km, 2));
	  double conf_Km = (1.0 - TMath::Erf(fabs(deltab_km + mean_Km)/sigma_Km/sqrt(2.0)));        
			    
	  if( prob_pim  > prob_Km  && hypothesis == -211 ){
	    //&& conf_pip > pip_conf && conf_prot < pip_anticonf  && conf_Kp  < pip_anticonf){ 
	    //std::cout << " pion plus " << std::endl;
	    //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	    //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	    // cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	    //cout << endl;

	    if( p[k] > temp_p_pim ){
	      temp_p_pim = p[k];
	      temp_index_pim = k;
	      temp_conf_pim = conf_pim;
	    }
	    //return true;
	  }
	  if( prob_Km > prob_pim && hypothesis == -321 ){

	    //std::cout << " KAON plus " << std::endl;
	    // && conf_Kp > kp_conf  && conf_prot < pr_anticonf  && conf_pip < pr_anticonf ){ 
	    //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	    //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	    //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	    //cout << endl;
	    if( p[k] > temp_p_km ){
	      temp_p_km = p[k];
	      temp_index_km = k;
	      temp_conf_km = conf_Km;
	    }
	    //return true;
	  }
	}
      }
    }
  std::map<int, std::vector<double> > hadron_candidates;

  std::vector<double> pim_classifier_values;
  pim_classifier_values.push_back(temp_index_pim);
  pim_classifier_values.push_back(temp_conf_pim);

  std::vector<double> kaonM_classifier_values;
  kaonM_classifier_values.push_back(temp_index_km);
  kaonM_classifier_values.push_back(temp_conf_km);

  hadron_candidates[-211] = pim_classifier_values;
  hadron_candidates[-321] = kaonM_classifier_values;

  return hadron_candidates;

}

