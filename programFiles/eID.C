#include "eIDsubroutines.C"
#include "vertexCorr.C"
#include "sctimeCorr.C"

int eID(Int_t gpart, Int_t q[], Float_t p[], UChar_t cc_sect[], UChar_t sc_sect[], UChar_t ec_sect[], UChar_t dc_sect[], Float_t cx[], Float_t cy[], Float_t cz[], Float_t tl1_x[], Float_t tl1_y[], Float_t tl3_x[], Float_t tl3_y[], Float_t tl3_z[], Float_t tl3_cx[], Float_t tl3_cy[], Float_t tl3_cz[], int e_zvertex_strict, Float_t vz[], Float_t vy[], Float_t vx[], int e_ECsampling_strict, std::map<int, std::vector<double> > el_ecsf_mean, std::map<int,std::vector<double> > el_ecsf_sigma, int ExpOrSim, Float_t etot[], int e_ECoutVin_strict, Float_t ec_ei[], Float_t ec_eo[], Float_t ech_x[], Float_t ech_y[], Float_t ech_z[], int e_CCthetaMatching_strict, UShort_t cc_segm[], int e_ECgeometric_strict, int e_R1fid_strict, int e_R3fid_strict, int e_CCphiMatching_strict, UChar_t sc_pd[], int e_CCfiducial_strict, int e_CCnphe_strict, UShort_t nphe[], int el_cut_pass_rate[], std::map<int, std::vector<TH1D*> > h_p, std::map<int, std::vector<TH2D*> > h_ectotp, std::map<int , std::vector<TH2D*> > h_eceo_ecei, std::map<int, std::vector<TH2D*> > h_ec_pos, std::map<int, std::vector<TH1D*> > h_nphe, std::map< int, std::vector<TH2D*> > h_theta_phi, std::map< int, std::vector<TH2D*> > h_dcr1, std::map<int, std::vector<TH2D*> > h_dcr3, std::map<int, std::vector<TH2D*> > h_cc_theta, std::map<int, std::vector<TH1D*> > h_vz, std::map<int, std::vector<TH1D*> > h_corr_vz)
{
  vector<int> Veindex;
  for(int k = 0; k < gpart; k++) // loop over particles
    {
      //std::cout << " charge " << q[k] << std::endl;
      if(q[k] == -1 )//{// && cc_sect[k] != 0 && sc_sect[k] != 0 && ec_sect[k] != 0 && dc_sect[k] != 0 ) /// && goodORbadSCpaddle(dc_sect[k], sc_pd[k]) == 1)
	{
	  //std::cout << " pass charge " << std::endl;
	  //el_individual_cut_pass_rate[0]++;

	  Float_t phi = atan3(cy[k],cx[k])*57.2957795;
	  Float_t relphi = get_rel_phi2(shift180180to30330(atan2(cy[k],cx[k])*57.2957795), dc_sect[k]);
	  Float_t thetaCC = 57.2957795*get_thetaCC(tl3_x[k], tl3_y[k], tl3_z[k], tl3_cx[k], tl3_cy[k], tl3_cz[k]);

	  // check if particle passes eID cuts:
	  Bool_t zvertex_pass =false, ECsampling_pass = false, ECoutVin_pass = false, ECgeometric_pass = false, CCthetaMatching_pass = false, R1fid_pass = false, R3fid_pass = false, CCphiMatching_pass = false, CCfiducial_pass = false, CCnphe_pass = false;
	  //zvertex_pass = ECsampling_pass = ECoutVin_pass = ECgeometric_pass = CCthetaMatching_pass = R1fid_pass = R3fid_pass = CCphiMatching_pass = CCfiducial_pass = 0;
	  //std::cout << " vertex number " << getCorrZ(ExpOrSim, vx[k], vy[k], vz[k], p[k]*cx[k], p[k]*cy[k], p[k]*cz[k], dc_sect[k]) << std::endl;
	  zvertex_pass = e_zvertex_pass(e_zvertex_strict, ExpOrSim, getCorrZ(ExpOrSim, vx[k], vy[k], vz[k], p[k]*cx[k], p[k]*cy[k], p[k]*cz[k], dc_sect[k]));
	  double corr_vz = getCorrZ(ExpOrSim, vx[k], vy[k], vz[k], p[k]*cx[k], p[k]*cy[k], p[k]*cz[k], dc_sect[k]);
	  // fill histograms for all negatively charged particles
	  //std::cout<< " dc sect " <<  (int)dc_sect[k] << std::endl;
	  h_p[0][(int)dc_sect[k]]->Fill(p[k]);
	  h_vz[0][(int)dc_sect[k]]->Fill(vz[k]);
	  h_corr_vz[0][(int)dc_sect[k]]->Fill(corr_vz);
	  h_ectotp[0][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	  h_eceo_ecei[0][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	  h_ec_pos[0][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	  h_nphe[0][(int)cc_sect[k]]->Fill(nphe[k]);
	  h_theta_phi[0][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	  h_dcr1[0][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	  h_dcr3[0][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	  int pmt = cc_segm[k]/1000 -1;
	  int segment = cc_segm[k]%1000/10;	  
	  h_cc_theta[0][(int)sc_sect[k]]->Fill(segment,thetaCC);

	  h_p[0][0]->Fill(p[k]);
	  h_vz[0][0]->Fill(vz[k]);
	  h_corr_vz[0][0]->Fill(corr_vz);
	  h_ectotp[0][0]->Fill(p[k], etot[k]/p[k]);
	  h_eceo_ecei[0][0]->Fill( ec_ei[k], ec_eo[k] );
	  h_ec_pos[0][0]->Fill(ech_x[k], ech_y[k]);
	  h_nphe[0][0]->Fill(nphe[k]);
	  h_theta_phi[0][0]->Fill(relphi,thetaCC);
	  h_dcr1[0][0]->Fill(tl1_x[k], tl1_y[k]);
	  h_dcr3[0][0]->Fill(tl3_x[k], tl3_y[k]);
	  h_cc_theta[0][0]->Fill(segment,thetaCC);

	  
	  if( zvertex_pass ){
	    h_p[1][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[1][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[1][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[1][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[1][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[1][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[1][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[1][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[1][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[1][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[1][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[1][0]->Fill(p[k]);
	    h_vz[1][0]->Fill(vz[k]);                                                                                                                                                                    
 	    h_corr_vz[1][0]->Fill(corr_vz); 
	    h_ectotp[1][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[1][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[1][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[1][0]->Fill(nphe[k]);
	    h_theta_phi[1][0]->Fill(relphi,thetaCC);
	    h_dcr1[1][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[1][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[1][0]->Fill(segment,thetaCC);
	  }	  

	  if( e_CCfiducial_pass(e_CCfiducial_strict, thetaCC, relphi) ){
	    h_p[2][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[2][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[2][(int)dc_sect[k]]->Fill(corr_vz);	 
	    h_ectotp[2][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[2][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[2][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[2][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[2][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[2][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[2][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[2][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[2][0]->Fill(p[k]);
	    h_vz[2][0]->Fill(vz[k]); 
 	    h_corr_vz[2][0]->Fill(corr_vz); 
	    h_ectotp[2][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[2][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[2][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[2][0]->Fill(nphe[k]);
	    h_theta_phi[2][0]->Fill(relphi,thetaCC);
	    h_dcr1[2][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[2][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[2][0]->Fill(segment,thetaCC);
	  }	  

	  if( e_CCphiMatching_pass(e_CCphiMatching_strict, ccphimatching(cc_segm[k], phi)) ){
	    h_p[3][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[3][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[3][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[3][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[3][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[3][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[3][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[3][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[3][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[3][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[3][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[3][0]->Fill(p[k]);
	    h_vz[3][0]->Fill(vz[k]); 
 	    h_corr_vz[3][0]->Fill(corr_vz); 
	    h_ectotp[3][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[3][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[3][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[3][0]->Fill(nphe[k]);
	    h_theta_phi[3][0]->Fill(relphi,thetaCC);
	    h_dcr1[3][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[3][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[3][0]->Fill(segment,thetaCC);
	  }	  

	  if( e_CCthetaMatching_pass(e_CCthetaMatching_strict, ExpOrSim, dc_sect[k], thetaCC, (cc_segm[k]%1000/10)) ){
	    h_p[4][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[4][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[4][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[4][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[4][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[4][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[4][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[4][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[4][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[4][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[4][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[4][0]->Fill(p[k]);
	    h_vz[4][0]->Fill(vz[k]); 
	    h_corr_vz[4][0]->Fill(corr_vz); 
 	    h_ectotp[4][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[4][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[4][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[4][0]->Fill(nphe[k]);
	    h_theta_phi[4][0]->Fill(relphi,thetaCC);
	    h_dcr1[4][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[4][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[4][0]->Fill(segment,thetaCC);

	  }	  

	  if( e_R1fid_pass(e_R1fid_strict, dc_sect[k], tl1_x[k], tl1_y[k]) ){
	    h_p[5][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[5][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[5][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[5][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[5][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[5][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[5][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[5][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[5][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[5][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[5][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[5][0]->Fill(p[k]);
	    h_vz[5][0]->Fill(vz[k]); 
	    h_corr_vz[5][0]->Fill(corr_vz); 
 	    h_ectotp[5][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[5][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[5][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[5][0]->Fill(nphe[k]);
	    h_theta_phi[5][0]->Fill(relphi,thetaCC);
	    h_dcr1[5][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[5][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[5][0]->Fill(segment,thetaCC);
	  }	  

	  if( e_R3fid_pass(e_R3fid_strict, dc_sect[k], tl3_x[k], tl3_y[k]) ){
	    h_p[6][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[6][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[6][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[6][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[6][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[6][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[6][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[6][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[6][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[6][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[6][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[6][0]->Fill(p[k]);
	    h_vz[6][0]->Fill(vz[k]); 
	    h_corr_vz[6][0]->Fill(corr_vz); 
 	    h_ectotp[6][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[6][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[6][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[6][0]->Fill(nphe[k]);
	    h_theta_phi[6][0]->Fill(relphi,thetaCC);
	    h_dcr1[6][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[6][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[6][0]->Fill(segment,thetaCC);
	  }	  

	  if( e_ECgeometric_pass(e_ECgeometric_strict, ech_x[k], ech_y[k], ech_z[k])){
	    h_p[7][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[7][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[7][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[7][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[7][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[7][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[7][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[7][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[7][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[7][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[7][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[7][0]->Fill(p[k]);
	    h_vz[7][0]->Fill(vz[k]); 
	    h_corr_vz[7][0]->Fill(corr_vz); 
 	    h_ectotp[7][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[7][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[7][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[7][0]->Fill(nphe[k]);
	    h_theta_phi[7][0]->Fill(relphi,thetaCC);
	    h_dcr1[7][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[7][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[7][0]->Fill(segment,thetaCC);
	  }	  

 	  if( e_ECoutVin_pass(e_ECoutVin_strict, ec_ei[k])){
	    h_p[8][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[8][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[8][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[8][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[8][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[8][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[8][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[8][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[8][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[8][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[8][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[8][0]->Fill(p[k]);
	    h_vz[8][0]->Fill(vz[k]); 
 	    h_corr_vz[8][0]->Fill(corr_vz); 
	    h_ectotp[8][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[8][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[8][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[8][0]->Fill(nphe[k]);
	    h_theta_phi[8][0]->Fill(relphi,thetaCC);
	    h_dcr1[8][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[8][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[8][0]->Fill(segment,thetaCC);
	  }	  

	  if( e_ECsampling_pass(e_ECsampling_strict, ExpOrSim, dc_sect[k], etot[k], p[k], el_ecsf_mean, el_ecsf_sigma)){
	    h_p[9][(int)dc_sect[k]]->Fill(p[k]);
	    h_vz[9][(int)dc_sect[k]]->Fill(vz[k]);
	    h_corr_vz[9][(int)dc_sect[k]]->Fill(corr_vz);
	    h_ectotp[9][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[9][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[9][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[9][(int)cc_sect[k]]->Fill(nphe[k]);
	    h_theta_phi[9][(int)dc_sect[k]]->Fill(relphi,thetaCC);
	    h_dcr1[9][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[9][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
	    int pmt = cc_segm[k]/1000 -1;
	    int segment = cc_segm[k]%1000/10;	  
	    h_cc_theta[9][(int)sc_sect[k]]->Fill(segment,thetaCC);

	    h_p[9][0]->Fill(p[k]);
	    h_vz[9][0]->Fill(vz[k]); 
 	    h_corr_vz[9][0]->Fill(corr_vz); 
	    h_ectotp[9][0]->Fill(p[k], etot[k]/p[k]);
	    h_eceo_ecei[9][0]->Fill( ec_ei[k], ec_eo[k] );
	    h_ec_pos[9][0]->Fill(ech_x[k], ech_y[k]);
	    h_nphe[9][0]->Fill(nphe[k]);
	    h_theta_phi[9][0]->Fill(relphi,thetaCC);
	    h_dcr1[9][0]->Fill(tl1_x[k], tl1_y[k]);
	    h_dcr3[9][0]->Fill(tl3_x[k], tl3_y[k]);
	    h_cc_theta[9][0]->Fill(segment,thetaCC);
	  }	  
	  
	  if (p[k] > 0.5 && p[k] < 6.0 ){
	    if(zvertex_pass) {
	      ECsampling_pass = e_ECsampling_pass(e_ECsampling_strict, ExpOrSim, dc_sect[k], etot[k], p[k], el_ecsf_mean, el_ecsf_sigma ); el_cut_pass_rate[0]++; 
	      //std::cout << " pass vertex " << std::endl;
	      //}
	      if(zvertex_pass && ECsampling_pass) {
		ECoutVin_pass = e_ECoutVin_pass(e_ECoutVin_strict, ec_ei[k]); el_cut_pass_rate[1]++;
		//std::cout << " pass ecsf " << std::endl;
		//}
		if(zvertex_pass && ECsampling_pass && ECoutVin_pass) {
		  ECgeometric_pass = e_ECgeometric_pass(e_ECgeometric_strict, ech_x[k], ech_y[k], ech_z[k]);  el_cut_pass_rate[2]++;
		  //std::cout << " pass ecoutvin " << std::endl;
		  //}
		  if(zvertex_pass && ECsampling_pass && ECoutVin_pass && ECgeometric_pass) {
		    CCthetaMatching_pass = e_CCthetaMatching_pass(e_CCthetaMatching_strict, ExpOrSim, dc_sect[k], thetaCC, (cc_segm[k]%1000)/10);  el_cut_pass_rate[3]++;
		    //std::cout << " pass ec fiducual " << std::endl;
		  //}
		    if(zvertex_pass && ECsampling_pass && ECoutVin_pass && ECgeometric_pass && CCthetaMatching_pass) {
		      R1fid_pass = e_R1fid_pass(e_R1fid_strict, dc_sect[k], tl1_x[k], tl1_y[k]);  el_cut_pass_rate[4]++; 
		      //std::cout << " pass cctheta " << std::endl;
		      //}
		      if(zvertex_pass && ECsampling_pass && ECoutVin_pass && ECgeometric_pass && CCthetaMatching_pass && R1fid_pass) {
			R3fid_pass = e_R3fid_pass(e_R3fid_strict, dc_sect[k], tl3_x[k], tl3_y[k]);  el_cut_pass_rate[5]++;
			//std::cout << " pass R1fid " << std::endl;
			//}
			if(zvertex_pass && ECsampling_pass && ECoutVin_pass && ECgeometric_pass && CCthetaMatching_pass && R1fid_pass && R3fid_pass) {
			  CCphiMatching_pass = e_CCphiMatching_pass(e_CCphiMatching_strict, ccphimatching(cc_segm[k], phi));  el_cut_pass_rate[6]++;
			  //std::cout << " pass R3fid " << std::endl;
			  //}
			  if(zvertex_pass && ECsampling_pass && ECoutVin_pass && ECgeometric_pass && CCthetaMatching_pass && R1fid_pass && R3fid_pass && CCphiMatching_pass) {
			    CCfiducial_pass = e_CCfiducial_pass(e_CCfiducial_strict, thetaCC, relphi);  el_cut_pass_rate[7]++;
			    //std::cout<< " pass ccphi " << std::endl;
			    //}
			    if(zvertex_pass && ECsampling_pass && ECoutVin_pass && ECgeometric_pass && CCthetaMatching_pass && R1fid_pass && R3fid_pass && CCphiMatching_pass && CCfiducial_pass) { 
			      CCnphe_pass = e_CCNphe_pass(e_CCnphe_strict,nphe[k]);  el_cut_pass_rate[8]++; 
			      //std::cout << " pass ccfid " << std::endl;
			      //}
			      if(zvertex_pass && ECsampling_pass && ECoutVin_pass && ECgeometric_pass && CCthetaMatching_pass && R1fid_pass && R3fid_pass && CCphiMatching_pass && CCfiducial_pass && CCnphe_pass) {
				Veindex.push_back(k);  el_cut_pass_rate[9]++; 
				//std::cout << " pass ccnphe " << std::endl;
				h_p[10][(int)dc_sect[k]]->Fill(p[k]);
				h_vz[10][(int)dc_sect[k]]->Fill(vz[k]);
				h_corr_vz[10][(int)dc_sect[k]]->Fill(corr_vz);
				h_ectotp[10][(int)dc_sect[k]]->Fill(p[k], etot[k]/p[k]);
				h_eceo_ecei[10][(int)ec_sect[k]]->Fill( ec_ei[k], ec_eo[k] );
				h_ec_pos[10][(int)ec_sect[k]]->Fill(ech_x[k], ech_y[k]);
				h_nphe[10][(int)cc_sect[k]]->Fill(nphe[k]);
				h_theta_phi[10][(int)dc_sect[k]]->Fill(relphi,thetaCC);
				h_dcr1[10][(int)dc_sect[k]]->Fill(tl1_x[k], tl1_y[k]);
				h_dcr3[10][(int)dc_sect[k]]->Fill(tl3_x[k], tl3_y[k]);
				int pmt = cc_segm[k]/1000 -1;
				int segment = cc_segm[k]%1000/10;	  
				h_cc_theta[10][(int)sc_sect[k]]->Fill(segment,thetaCC);

				//integrated over all sectors
				h_p[10][0]->Fill(p[k]);
				h_vz[10][0]->Fill(vz[k]); 
				h_corr_vz[10][0]->Fill(corr_vz); 
				h_ectotp[10][0]->Fill(p[k], etot[k]/p[k]);
				h_eceo_ecei[10][0]->Fill( ec_ei[k], ec_eo[k] );
				h_ec_pos[10][0]->Fill(ech_x[k], ech_y[k]);
				h_nphe[10][0]->Fill(nphe[k]);
				h_theta_phi[10][0]->Fill(relphi,thetaCC);
				h_dcr1[10][0]->Fill(tl1_x[k], tl1_y[k]);
				h_dcr3[10][0]->Fill(tl3_x[k], tl3_y[k]);
				h_cc_theta[10][0]->Fill(segment,thetaCC);
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }	
	  }
	}
    } // end of loop over gpart
  
  // %%%%% find highest momentum electron %%%%%
  int electron_index = -123;
  if(Veindex.size() > 0) electron_index = Veindex[0];
  if(Veindex.size() > 1){
    for(unsigned int v = 1; v < Veindex.size(); v++){
      //std::cout << " p electron " << p[electron_index] << std::endl;
      if(p[Veindex[v]] > p[electron_index]) electron_index = Veindex[v];
    }
  }
  
  return electron_index;
}
