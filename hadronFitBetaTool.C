#include <fstream>
#include <iostream>
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "programFiles/functions.C"
#include "programFiles/eID.C"
#include "programFiles/hadronID.C"
#include "programFiles/getGenIndices.C"
#include "MomCorr.C"
#include "CalculatorTools.C"


// change the size of the beta vs p histogram range to momentum of 2.2 (smaller than now and then refit) results will be better


int hadronFitBetaTool( int filestart = 1, int Nfiles = 1){

  TStopwatch *stopwatch = new TStopwatch();
  
  TFile *fIn = new TFile("hadron_beta_out.root");  

  Float_t e_mass = 0.000511; // GeV
  Float_t prot_mass = 0.938272; // GeV
  Float_t pip_mass = 0.13957; // GeV
  Float_t pim_mass = 0.13957; // GeV
  Float_t kp_mass = 0.4937; // GeV
  Float_t km_mass = 0.4937; // GeV
  Float_t speed_of_light = 29.9792458; // cm/ns
  Float_t Beam_Energy = 5.498; // GeV
  Float_t pi = 3.14159265359;
  Float_t pi180 = pi/180.0;
  Float_t pi180_inv = 180.0/pi;
  int ExpOrSim = 1;

  std::vector<TH2F*> h2_betap_pr;
  std::vector<TH2F*> h2_betap_pip;
  std::vector<TH2F*> h2_betap_kp;

  std::vector<TH2F*> h2_betap_pim;
  std::vector<TH2F*> h2_betap_km;

  //create vectors of beta vs p histograms
  for( int s = 0; s < 6; s++ ){
    std::cout << " sector " << s << std::endl;
    //positive hadrons
    h2_betap_pr.push_back( (TH2F*)fIn->Get(Form("h_delta_beta_pr_sector_%d", s+1)) );
    h2_betap_pip.push_back( (TH2F*)fIn->Get(Form("h_delta_beta_pip_sector_%d", s+1)) );
    h2_betap_kp.push_back( (TH2F*)fIn->Get(Form("h_delta_beta_kp_sector_%d", s+1)) );

    //negative hadrons
    h2_betap_pim.push_back( (TH2F*)fIn->Get(Form("h_delta_beta_pim_sector_%d", s+1)) );
    h2_betap_km.push_back( (TH2F*)fIn->Get(Form("h_beta_all_km_sector_%d", s+1)) );

    std::cout << ((TH2F*)fIn->Get(Form("h_delta_beta_pip_sector_%d", s+1)))->GetEntries() << std::endl;
  }

  std::vector<TGraphErrors*> h_beta_mean_pr;
  std::vector<TGraphErrors*> h_beta_mean_pip;
  std::vector<TGraphErrors*> h_beta_mean_kp;

  std::vector<TGraphErrors*> h_beta_mean_km;
  std::vector<TGraphErrors*> h_beta_mean_pim;  
  
  std::vector<TGraphErrors*> h_beta_sig_pr;
  std::vector<TGraphErrors*> h_beta_sig_kp;
  std::vector<TGraphErrors*> h_beta_sig_pip;

  std::vector<TGraphErrors*> h_beta_sig_km;
  std::vector<TGraphErrors*> h_beta_sig_pim;
  
  //create vector of tgraphs for means and sig to fit 
  //fits will be used for the maximum likelihood estimators

  std::vector< std::vector<TCanvas*> > v_can;
  std::vector< int > partID  = { 2212, 211, 321, -211, -321 };

  for( int i = 0; i < partID.size(); i++ ){
    //v_can.push_back( std::vector<TCanvas*> );
    std::vector<TCanvas*> temp_canv;
    int pid = partID[i];
    for( int s = 0; s <= 5; s++ ){
      temp_canv.push_back(new TCanvas(Form("c_p%d_%d",pid,s+1), Form("c_p%d_%d",pid,s+1), 800, 800 ));
      //v_can[i].push_back( new TCanvas(Form("c_p%d_%d",pid,s), Form("c_p%d_%d",pid,s), 800, 800 ));
      //v_can[i][s]->Divide(5,5);
      temp_canv[s]->Divide(5,5);
    }
    v_can.push_back( temp_canv );
  }

  for( int i = 0; i < h2_betap_pr.size(); i++ ){

    int nbins_x_pr = h2_betap_pr[i]->GetNbinsX();
    int nbins_x_pip = h2_betap_pip[i]->GetNbinsX();
    int nbins_x_kp = h2_betap_kp[i]->GetNbinsX();

    int nbins_x_pim = h2_betap_pim[i]->GetNbinsX();
    int nbins_x_km = h2_betap_km[i]->GetNbinsX();

    int ndiv_pr = 5;
    int ndiv_pip = 1;
    int ndiv_kp = 5;

    int ndiv_pim = 1;
    int ndiv_km = 10;//5;

    int nfits_pr = nbins_x_pr / ndiv_pr;
    int nfits_pip = nbins_x_pip / ndiv_pip;
    int nfits_kp = nbins_x_kp / ndiv_kp;

    int nfits_pim = nbins_x_pim / ndiv_pim;
    int nfits_km = nbins_x_km / ndiv_km;

    int start_bin_pr = 0;
    int start_bin_pip = 0;
    int start_bin_kp = 0;
    
    int start_bin_pim = 0;
    int start_bin_km = 0;

    int end_bin_pr = ndiv_pr;
    int end_bin_pip = ndiv_pip;
    int end_bin_kp = ndiv_kp;

    int end_bin_pim = ndiv_pim;
    int end_bin_km = ndiv_km;

    //fill with the fit mean and sigma for the tgraphs

    std::vector<double> temp_center_pr;
    std::vector<double> temp_center_pip;
    std::vector<double> temp_center_kp;

    std::vector<double> temp_center_err_pr;
    std::vector<double> temp_center_err_pip;
    std::vector<double> temp_center_err_kp;

    std::vector<double> temp_fit_means_pr;    
    std::vector<double> temp_fit_means_pip;    
    std::vector<double> temp_fit_means_kp;    

    std::vector<double> temp_fit_sigmas_pr;    
    std::vector<double> temp_fit_sigmas_pip;    
    std::vector<double> temp_fit_sigmas_kp;    

    std::vector<double> temp_fit_means_err_pr;    
    std::vector<double> temp_fit_means_err_pip;    
    std::vector<double> temp_fit_means_err_kp;    

    std::vector<double> temp_fit_sigmas_err_pr;    
    std::vector<double> temp_fit_sigmas_err_pip;    
    std::vector<double> temp_fit_sigmas_err_kp;    

    std::vector<double> temp_center_km;
    std::vector<double> temp_center_pim;
    
    std::vector<double> temp_center_err_km;
    std::vector<double> temp_center_err_pim;
    
    std::vector<double> temp_fit_means_km;
    std::vector<double> temp_fit_means_pim;

    std::vector<double> temp_fit_sigmas_km;
    std::vector<double> temp_fit_sigmas_pim;

    std::vector<double> temp_fit_means_err_km;
    std::vector<double> temp_fit_means_err_pim;

    std::vector<double> temp_fit_sigmas_err_km;
    std::vector<double> temp_fit_sigmas_err_pim;

    std::cout << " PERFORMING FITS N TIMES FOR PROTONS: " << nfits_pr << std::endl;
    
    for( int j = 0; j < nfits_pr; j++ ){
      std::cout << " >> FITTING BETWEEN BINS " << start_bin_pr << " " << end_bin_pr << std::endl;

      TH1D *h_beta =  h2_betap_pr[i]->ProjectionY(Form("pr_betap_projY_s%d_%d",i,j),start_bin_pr, end_bin_pr );
      //TF1 *fit_beta = fitHistogram(h_beta);

      double beta_min = 0.0; double beta_max = 0.0;
      double current_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(start_bin_pr);
      if( current_p < 2.0 ){ beta_min = -0.08; beta_max = 0.08; }
      else if( current_p >= 2.0 ){ beta_min = -0.04; beta_max = 0.04; }

      TF1 *fit_beta = new TF1("fit_temp_pr","gaus",beta_min,beta_max);
      fit_beta->SetParameter(0,h_beta->GetMaximum());
      fit_beta->SetParameter(1,0.0);
      fit_beta->SetParameter(2,0.0);
      h_beta->Fit(fit_beta,"","",beta_min,beta_max);
      
      v_can[0][i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(start_bin_pr) - h2_betap_pr[i]->GetXaxis()->GetBinWidth(start_bin_pr)/2.0;
      double max_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(end_bin_pr) + h2_betap_pr[i]->GetXaxis()->GetBinWidth(end_bin_pr)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);


      start_bin_pr+=ndiv_pr;
      end_bin_pr+=ndiv_pr;           

      //if ( j < 2 || j > 15 ) continue;      
      
      temp_center_pr.push_back(slice_center);
      temp_center_err_pr.push_back(0.0);

      temp_fit_means_err_pr.push_back(mean_err);
      temp_fit_sigmas_err_pr.push_back(sig_err);
    
      temp_fit_means_pr.push_back(mean);
      temp_fit_sigmas_pr.push_back(sig);
     
    
    }

    

    // %%%%%%%%%%% FITTING PIP BETA VS P NOW %%%%%%%%%%%%%%%%%
    std::cout << " PERFORMING FITS N TIMES PIP: " << nfits_pip << std::endl;
    for( int j = 0; j < nfits_pip; j++ ){
      std::cout << " >> FITTING BETWEEN BINS " << start_bin_pip << " " << end_bin_pip << std::endl;


      TH1D *h_beta =  h2_betap_pip[i]->ProjectionY(Form("pip_betap_projY_s%d_%d",i,j),start_bin_pip, end_bin_pip );
      //TF1 *fit_beta = fitHistogram(h_beta);
      TF1 *fit_beta = new TF1("fit_temp_pip","gaus",-0.08,0.08);
      fit_beta->SetParameter(0,h_beta->GetMaximum());
      fit_beta->SetParameter(1,0.0);
      fit_beta->SetParameter(2,0.0);

      double beta_min = 0.0; double beta_max = 0.0;
      double current_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(start_bin_pr);
      if( current_p < 2.5 ){ beta_min = -0.03; beta_max = 0.03; }
      else if( current_p >= 2.5 ){ beta_min = -0.02; beta_max = 0.02; }

      h_beta->Fit(fit_beta,"","",beta_min,beta_max);

      v_can[1][i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_pip[i]->GetXaxis()->GetBinCenter(start_bin_pip) - h2_betap_pip[i]->GetXaxis()->GetBinWidth(start_bin_pip)/2.0;
      double max_p = h2_betap_pip[i]->GetXaxis()->GetBinCenter(end_bin_pip) + h2_betap_pip[i]->GetXaxis()->GetBinWidth(end_bin_pip)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);

      start_bin_pip+=ndiv_pip;
      end_bin_pip+=ndiv_pip;           

      //if ( j < 8 || j > (nfits_pip - 6) ) continue;            
      temp_center_pip.push_back(slice_center);
      temp_center_err_pip.push_back(0.0);

      temp_fit_means_err_pip.push_back(mean_err);
      temp_fit_sigmas_err_pip.push_back(sig_err);
    
      temp_fit_means_pip.push_back(mean);
      temp_fit_sigmas_pip.push_back(sig);
         
    }
    
    
    // %%%%%%%%%%%%%% FITTING KAON PLUS BETA VS P NOW %%%%%%%%%%%%%%%
    std::cout << " PERFORMING FITS N TIMES: " << nfits_kp << std::endl;
    for( int j = 0; j < nfits_kp; j++ ){
      std::cout << " >> FITTING BETWEEN BINS " << start_bin_kp << " " << end_bin_kp << std::endl;

      TH1D *h_beta =  h2_betap_kp[i]->ProjectionY(Form("kp_betap_projY_s%d_%d",i,j),start_bin_kp, end_bin_kp );
      //TF1 *fit_beta = fitHistogram(h_beta);
      TF1 *fit_beta = new TF1("fit_temp","gaus",-0.08,0.08);
      fit_beta->SetParameter(0,h_beta->GetMaximum());
      fit_beta->SetParameter(1,0.0);
      fit_beta->SetParameter(2,0.0);

      double beta_min = 0.0; double beta_max = 0.0;
      double current_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(start_bin_pr);
      if( current_p < 2.0 ){ beta_min = -0.021; beta_max = 0.021; }
      else if( current_p >= 2.0 ){ beta_min = -0.018; beta_max = 0.018; }

      h_beta->Fit(fit_beta,"","", beta_min, beta_max);
     
      v_can[2][i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_kp[i]->GetXaxis()->GetBinCenter(start_bin_kp) - h2_betap_kp[i]->GetXaxis()->GetBinWidth(start_bin_kp)/2.0;
      double max_p = h2_betap_kp[i]->GetXaxis()->GetBinCenter(end_bin_kp) + h2_betap_kp[i]->GetXaxis()->GetBinWidth(end_bin_kp)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);

      start_bin_kp+=ndiv_kp;
      end_bin_kp+=ndiv_kp;           

      //if ( j < 2 || j > 6 ) continue;      
      temp_center_kp.push_back(slice_center);
      temp_center_err_kp.push_back(0.0);

      temp_fit_means_err_kp.push_back(mean_err);
      temp_fit_sigmas_err_kp.push_back(sig_err);
    
      temp_fit_means_kp.push_back(mean);
      temp_fit_sigmas_kp.push_back(sig);
         
    }

    // %%%%%%%%%%%%%% FITTING PION MINUS BETA VS P NOW %%%%%%%%%%%%%%%
    std::cout << " PERFORMING FITS N TIMES: " << nfits_pim << std::endl;
    for( int j = 0; j < nfits_pim; j++ ){
      std::cout << " >> FITTING BETWEEN BINS " << start_bin_pim << " " << end_bin_pim << std::endl;

      TH1D *h_beta =  h2_betap_pim[i]->ProjectionY(Form("pim_betap_projY_s%d_%d",i,j),start_bin_pim, end_bin_pim );
      //TF1 *fit_beta = fitHistogram(h_beta);

      TF1 *fit_beta = new TF1("fit_temp","gaus",-0.08,0.08);
      fit_beta->SetParameter(0,h_beta->GetMaximum());
      fit_beta->SetParameter(1,0.0);
      fit_beta->SetParameter(2,0.0);
      
      double beta_min = 0.0; double beta_max = 0.0;
      double current_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(start_bin_pr);
      if( current_p < 2.0 ){ beta_min = -0.08; beta_max = 0.08; }
      else if( current_p >= 2.0 ){ beta_min = -0.04; beta_max = 0.04; }

      h_beta->Fit(fit_beta,"","",beta_min,beta_max);

      v_can[3][i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_pim[i]->GetXaxis()->GetBinCenter(start_bin_pim) - h2_betap_pim[i]->GetXaxis()->GetBinWidth(start_bin_pim)/2.0;
      double max_p = h2_betap_pim[i]->GetXaxis()->GetBinCenter(end_bin_pim) + h2_betap_pim[i]->GetXaxis()->GetBinWidth(end_bin_pim)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);

      start_bin_pim+=ndiv_pim;
      end_bin_pim+=ndiv_pim;           
      //if ( j < 2 || j > 15 ) continue;      
      
      temp_center_pim.push_back(slice_center);
      temp_center_err_pim.push_back(0.0);

      temp_fit_means_err_pim.push_back(mean_err);
      temp_fit_sigmas_err_pim.push_back(sig_err);
    
      temp_fit_means_pim.push_back(mean);
      temp_fit_sigmas_pim.push_back(sig);
         
    }

    // %%%%%%%%%%%%%% FITTING KAON MINUS BETA VS P NOW %%%%%%%%%%%%%%%
    std::cout << " PERFORMING KAON FITS N TIMES: " << nfits_km << std::endl;
    for( int j = 0; j < nfits_km; j++ ){
      std::cout << " >> FITTING KAON MINUS BETWEEN BINS " << start_bin_km << " " << end_bin_km << std::endl;

      //if ( j <= 7 ) continue;

      TH1D *h_beta =  h2_betap_km[i]->ProjectionY(Form("km_betap_projY_s%d_%d",i,j),start_bin_km, end_bin_km );
      TF1 *fit_beta = fitHistogram(h_beta);

      //TF1 *fit_beta = new TF1("fit_temp","gaus",-0.08,0.08);
      //fit_beta->SetParameter(0,h_beta->GetMaximum());
      //fit_beta->SetParameter(1,0.0);
      //fit_beta->SetParameter(2,0.0);
      //below arguments will fit the range in the histogram...
      //double beta_min = 0.0; double beta_max = 0.0;
      double current_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(start_bin_pr);
      //if( current_p < 2.0 ){ beta_min = -0.08; beta_max = 0.08; }
      //else if( current_p >= 2.0 ){ beta_min = -0.04; beta_max = 0.04; }

      //h_beta->Fit(fit_beta,"","",beta_min,beta_max);

      v_can[4][i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_km[i]->GetXaxis()->GetBinCenter(start_bin_km) - h2_betap_km[i]->GetXaxis()->GetBinWidth(start_bin_km)/2.0;
      double max_p = h2_betap_km[i]->GetXaxis()->GetBinCenter(end_bin_km) + h2_betap_km[i]->GetXaxis()->GetBinWidth(end_bin_km)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);

      start_bin_km+=ndiv_km;
      end_bin_km+=ndiv_km;                
      //if ( j < 2 || j > 10 ) continue;      

      temp_center_km.push_back(slice_center);
      temp_center_err_km.push_back(0.0);

      temp_fit_means_err_km.push_back(mean_err);
      temp_fit_sigmas_err_km.push_back(sig_err);
    
      temp_fit_means_km.push_back(mean);
      temp_fit_sigmas_km.push_back(sig);
     
    
    }
    

    h_beta_mean_pr.push_back( new TGraphErrors(temp_center_pr.size(), &(temp_center_pr[0]), &(temp_fit_means_pr[0]), &(temp_center_err_pr[0]), &(temp_fit_means_err_pr[0])) );
    h_beta_mean_pip.push_back( new TGraphErrors(temp_center_pip.size(), &(temp_center_pip[0]), &(temp_fit_means_pip[0]), &(temp_center_err_pip[0]), &(temp_fit_means_err_pip[0])) );
    h_beta_mean_kp.push_back( new TGraphErrors(temp_center_kp.size(), &(temp_center_kp[0]), &(temp_fit_means_kp[0]), &(temp_center_err_kp[0]), &(temp_fit_means_err_kp[0])) );
    h_beta_mean_pim.push_back( new TGraphErrors(temp_center_pim.size(), &(temp_center_pim[0]), &(temp_fit_means_pim[0]), &(temp_center_err_pim[0]), &(temp_fit_means_err_pim[0])) );
    h_beta_mean_km.push_back( new TGraphErrors(temp_center_km.size(), &(temp_center_km[0]), &(temp_fit_means_km[0]), &(temp_center_err_km[0]), &(temp_fit_means_err_km[0])) );

    h_beta_sig_pr.push_back( new TGraphErrors( temp_center_pr.size(), &(temp_center_pr[0]), &(temp_fit_sigmas_pr[0]), &(temp_center_err_pr[0]),&(temp_fit_sigmas_err_pr[0])) );
    h_beta_sig_pip.push_back( new TGraphErrors( temp_center_pip.size(), &(temp_center_pip[0]), &(temp_fit_sigmas_pip[0]), &(temp_center_err_pip[0]),&(temp_fit_sigmas_err_pip[0])) );
    h_beta_sig_kp.push_back( new TGraphErrors( temp_center_kp.size(), &(temp_center_kp[0]), &(temp_fit_sigmas_kp[0]), &(temp_center_err_kp[0]),&(temp_fit_sigmas_err_kp[0])) );
    h_beta_sig_pim.push_back( new TGraphErrors( temp_center_pim.size(), &(temp_center_pim[0]), &(temp_fit_sigmas_pim[0]), &(temp_center_err_pim[0]),&(temp_fit_sigmas_err_pim[0])) );
    h_beta_sig_km.push_back( new TGraphErrors( temp_center_km.size(), &(temp_center_km[0]), &(temp_fit_sigmas_km[0]), &(temp_center_err_km[0]),&(temp_fit_sigmas_err_km[0])) );
    
  }
    

  TCanvas *c_temp_pr = new TCanvas("c_final_pr", "c_final_pr", 800, 800);
  c_temp_pr->Divide(4,2);

  TCanvas *c_temp_pr_sig = new TCanvas("c_final_sig_pr", "c_final_sig_pr", 800, 800);
  c_temp_pr_sig->Divide(4,2);

  std::cout << " size " << h_beta_mean_pr.size() << std::endl;

  //std::map<std::string, std::vector<double> > prot_beta_fit_mean;
  std::vector<std::string> particle_names = {"PROT","KAONP","KAONM","PIONP","PIONM"};
  std::vector<std::string> par_names = {"_DBETA_MU_A","_DBETA_MU_B","_DBETA_MU_C"};//,"_DBETA_MU_D"};
  std::vector<std::string> par_names_sig = {"_DBETA_SIGMA_A","_DBETA_SIGMA_B","_DBETA_SIGMA_C"};

  //std::vector<double> prot_par_mean;
  //std::vector<double> prot_par_sig;

  std::vector<std::string> protonBetaMeanValues;
  std::vector<std::string> protonBetaSigValues;

  std::vector<std::string> pipBetaMeanValues;
  std::vector<std::string> pipBetaSigValues;

  std::vector<std::string> kpBetaMeanValues;
  std::vector<std::string> kpBetaSigValues;

  std::vector<std::string> pimBetaMeanValues;
  std::vector<std::string> pimBetaSigValues;

  std::vector<std::string> kmBetaMeanValues;
  std::vector<std::string> kmBetaSigValues;

  for( int i = 0; i < par_names.size(); i++){
    protonBetaMeanValues.push_back(particle_names[0]+par_names[i]+" \t" + "6" + "\t" + "v" + "\t");    
    kpBetaMeanValues.push_back(particle_names[1]+par_names[i]+ " \t" + "6" + "\t" + "v" + "\t");    
    kmBetaMeanValues.push_back(particle_names[2]+par_names[i]+" \t" + "6" + "\t" + "v" + "\t");
    pipBetaMeanValues.push_back(particle_names[3]+par_names[i]+ " \t" + "6" + "\t" + "v" + "\t");    
    pimBetaMeanValues.push_back(particle_names[4]+par_names[i]+ " \t" + "6" + "\t" + "v" + "\t");
    
  }

  for( int i = 0; i < par_names_sig.size(); i++ ){
    protonBetaSigValues.push_back(particle_names[0]+par_names_sig[i]+" \t" + "6" + "\t" + "v" + "\t");
    kpBetaSigValues.push_back(particle_names[1]+par_names_sig[i]+ " \t" + "6" + "\t" + "v" + "\t");
    kmBetaSigValues.push_back(particle_names[2]+par_names_sig[i]+" \t" + "6" + "\t" + "v" + "\t");
    pipBetaSigValues.push_back(particle_names[3]+par_names_sig[i]+ " \t" + "6" + "\t" + "v" + "\t");
    pimBetaSigValues.push_back(particle_names[4]+par_names_sig[i]+ " \t" + "6" + "\t" + "v" + "\t");
  }
  std::cout << " preparing string outs " << std::endl;

  for(int s = 0; s < 6; s++ ){
    c_temp_pr->cd(s+1);
  

    /// %%%%% MEAN FIT %%%%%%%
    TF1 *fit_mean_pr = new TF1(Form("fit_mean_pr_s%d",s),"[0] + [1]*x + [2]*x*x",0.6, 2.5) ;//[0] + [1]*x/sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.1, 4.0);
    fit_mean_pr->SetParameter(0,0.0);
    fit_mean_pr->SetParameter(1,0.0);
    fit_mean_pr->SetParameter(2,0.0);
    //fit_mean_pr->SetParameter(3,0.1);
    //fit_mean_pr->SetParameter(4,0.1);

    h_beta_mean_pr[s]->Fit(Form("fit_mean_pr_s%d",s),"R");//,"",-0.08,0.08);
    

    //fill proton beta mean fit results here
    //prot_beta_fit_mean[particle_names[0]+par_names[0]].push_back( fit_mean_pr->GetParameter(0) );

    //std::cout << " >> prot values " << protonBetaMeanValues[0] << std::endl;
    //std::cout << " >> prot values " << protonBetaMeanValues[1] << std::endl;
    //std::cout << " >> prot values " << protonBetaMeanValues[2] << std::endl;
    //std::cout << " >> prot values " << protonBetaMeanValues[3] << std::endl;
    //std::cout << std::to_string(fit_mean_pr->GetParameter(0)) << std::endl;

    h_beta_mean_pr[s]->SetTitle(Form("PR BETA MEAN FOR SECTOR %d", s));
    h_beta_mean_pr[s]->SetMarkerColor(kBlack);
    h_beta_mean_pr[s]->SetMarkerStyle(8);
    h_beta_mean_pr[s]->Draw("AP");       

    fit_mean_pr->SetLineColor(kRed);
    fit_mean_pr->Draw("same");

    c_temp_pr_sig->cd(s+1);

    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_pr = new TF1(Form("fit_sigma_pr_s%d",s),"[0] + [1]*x + [2]*x*x", 0.6, 2.5); // /x + [2]*exp(-x)", 0.1, 2.5);//4.0);
    fit_sigma_pr->SetParameter(0,0.0);
    fit_sigma_pr->SetParameter(1,0.0);
    fit_sigma_pr->SetParameter(2,0.0);//prot_mass);
    //fit_sigma_pr->SetParameter(3,0.1);
    //fit_sigma_pr->SetParameter(4,0.1);

    h_beta_sig_pr[s]->Fit(Form("fit_sigma_pr_s%d",s),"R");

    h_beta_sig_pr[s]->SetTitle(Form("PR BETA SIG FOR SECTOR %d", s));
    h_beta_sig_pr[s]->SetMarkerColor(kBlack);
    h_beta_sig_pr[s]->SetMarkerStyle(8);
    h_beta_sig_pr[s]->Draw("AP");       

    fit_sigma_pr->SetLineColor(kRed);
    fit_sigma_pr->Draw("same");

    //store fit parameters in string to write to text file
    protonBetaMeanValues[0] = protonBetaMeanValues[0] + std::to_string(fit_mean_pr->GetParameter(0)) + "\t";
    protonBetaMeanValues[1] = protonBetaMeanValues[1] + std::to_string(fit_mean_pr->GetParameter(1)) + "\t";
    protonBetaMeanValues[2] = protonBetaMeanValues[2] + std::to_string(fit_mean_pr->GetParameter(2)) + "\t";

    protonBetaSigValues[0] = protonBetaSigValues[0] + std::to_string(fit_sigma_pr->GetParameter(0)) + "\t";
    protonBetaSigValues[1] = protonBetaSigValues[1] + std::to_string(fit_sigma_pr->GetParameter(1)) + "\t";
    protonBetaSigValues[2] = protonBetaSigValues[2] + std::to_string(fit_sigma_pr->GetParameter(2)) + "\t";



  }

  
  TCanvas *c_temp_pip = new TCanvas("c_final_pip", "c_final_pip", 800, 800);
  c_temp_pip->Divide(3,2);

  TCanvas *c_temp_pip_sig = new TCanvas("c_final_sig_pip", "c_final_sig_pip", 800, 800);
  c_temp_pip_sig->Divide(3,2);

  std::cout << " size " << h_beta_mean_pip.size() << std::endl;
  for(int s = 0; s < 6; s++ ){
    c_temp_pip->cd(s+1);

    /// %%%%% MEAN FIT %%%%%%%
    TF1 *fit_mean_pip = new TF1(Form("fit_mean_pip_s%d",s),"[0] + [1]*x + [2]*x*x", 0.2, 3.5); ///sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.1, 4.0);
    fit_mean_pip->SetParameter(0,0.0);
    fit_mean_pip->SetParameter(1,0.010);
    fit_mean_pip->SetParameter(2,0.010);//pip_mass);
    //fit_mean_pip->SetParameter(3,0.0);
    //fit_mean_pip->SetParameter(4,0.10);   

    h_beta_mean_pip[s]->Fit(Form("fit_mean_pip_s%d",s),"R");

    h_beta_mean_pip[s]->SetTitle(Form("PIP BETA MEAN FOR SECTOR %d", s+1));
    h_beta_mean_pip[s]->SetMarkerColor(kBlack);
    h_beta_mean_pip[s]->SetMarkerStyle(8);
    h_beta_mean_pip[s]->Draw("AP");       

    fit_mean_pip->SetLineColor(kRed);
    fit_mean_pip->Draw();

    c_temp_pip_sig->cd(s+1);
    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_pip = new TF1(Form("fit_sigma_pip_s%d",s),"[0] + [1]*x + [2]*x*x", 0.2, 3.5); ////x + [2]*exp(-x)", 0.1, 2.5);// 4.0);
    fit_sigma_pip->SetParameter(0,0.0);
    fit_sigma_pip->SetParameter(1,0.0);//1.0);
    fit_sigma_pip->SetParameter(2,0.0);//pip_mass);
    //fit_sigma_pip->SetParameter(3,0.1);
    //fit_sigma_pip->SetParameter(4,0.1);

    h_beta_sig_pip[s]->Fit(Form("fit_sigma_pip_s%d",s),"R");

    h_beta_sig_pip[s]->SetTitle(Form("PIP BETA SIGMA FOR SECTOR %d", s+1));
    h_beta_sig_pip[s]->SetMarkerColor(kBlack);
    h_beta_sig_pip[s]->SetMarkerStyle(8);
    h_beta_sig_pip[s]->Draw("AP");       

    fit_sigma_pip->SetLineColor(kRed);
    fit_sigma_pip->SetLineStyle(3);
    fit_sigma_pip->Draw();

    //store fit parameters in string to write to text file
    for( int i = 0; i < fit_mean_pip->GetNpar(); i++ ){
      pipBetaMeanValues[i] = pipBetaMeanValues[i] + std::to_string(fit_mean_pip->GetParameter(i)) + "\t";
    }
    for( int i = 0; i < fit_sigma_pip->GetNpar(); i++ ){    
      pipBetaSigValues[i] = pipBetaSigValues[i] + std::to_string(fit_sigma_pip->GetParameter(i)) + "\t";
    }
    

  }
  

  TCanvas *c_temp_kp = new TCanvas("c_final_kp", "c_final_kp", 800, 800);
  c_temp_kp->Divide(4,2);

  TCanvas *c_temp_kp_sig = new TCanvas("c_final_sig_kp", "c_final_sig_kp", 800, 800);
  c_temp_kp_sig->Divide(4,2);

  std::cout << " size " << h_beta_mean_kp.size() << std::endl;
  for(int s = 0; s < 6; s++ ){
    c_temp_kp->cd(s+1);

    /// %%%%% MEAN FIT %%%%%%%%%
    TF1 *fit_mean_kp = new TF1(Form("fit_mean_kp_s%d",s),"[0] + [1]*x + [2]*x*x", 0.38, 1.4);///sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.80, 4.0);
    fit_mean_kp->SetParameter(0,0.0);
    fit_mean_kp->SetParameter(1,0.0);
    fit_mean_kp->SetParameter(2,0.0);//kp_mass);
    //fit_mean_kp->SetParameter(3,0.10);
    //fit_mean_kp->SetParameter(4,0.10);   

    h_beta_mean_kp[s]->Fit(Form("fit_mean_kp_s%d",s),"R");

    h_beta_mean_kp[s]->SetTitle(Form("KP BETA MEAN FOR SECTOR %d", s));
    h_beta_mean_kp[s]->SetMarkerColor(kBlack);
    h_beta_mean_kp[s]->SetMarkerStyle(8);
    h_beta_mean_kp[s]->Draw("AP");      

    fit_mean_kp->SetLineColor(kRed);
    fit_mean_kp->Draw("same"); 

    c_temp_kp_sig->cd(s+1);    
    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_kp = new TF1(Form("fit_sigma_kp_s%d",s),"[0] + [1]*x + [2]*x*x", 0.6, 1.38);///x + [2]*exp(-x)", 0.1, 2.5);
    fit_sigma_kp->SetParameter(0,0.0);
    fit_sigma_kp->SetParameter(1,0.0);
    fit_sigma_kp->SetParameter(2,0.0);//kp_mass);
    //fit_sigma_kp->SetParameter(3,0.0);
    //fit_sigma_kp->SetParameter(4,0.10);
    
    h_beta_sig_kp[s]->Fit(Form("fit_sigma_kp_s%d",s),"R");  

    h_beta_sig_kp[s]->SetTitle(Form("KP BETA SIGMA FOR SECTOR %d", s));
    h_beta_sig_kp[s]->SetMarkerColor(kBlack);
    h_beta_sig_kp[s]->SetMarkerStyle(8);
    h_beta_sig_kp[s]->Draw("AP");      
   
    fit_sigma_kp->SetLineColor(kRed);
    fit_sigma_kp->SetLineStyle(3);
    fit_sigma_kp->Draw("same"); 

    //store fit parameters in string to write to text file
    for( int i = 0; i < par_names.size(); i++ ){
      kpBetaMeanValues[i] = kpBetaMeanValues[i] + std::to_string(fit_mean_kp->GetParameter(i)) + "\t";
    }
    for( int i = 0; i < par_names_sig.size(); i++ ){    
      kpBetaSigValues[i] = kpBetaSigValues[i] + std::to_string(fit_sigma_kp->GetParameter(i)) + "\t";
    }
  }


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Plot Beta vs P fit results for Pion Minus

  TCanvas *c_temp_pim = new TCanvas("c_final_pim", "c_final_pim", 800, 800);
  c_temp_pim->Divide(4,2);

  TCanvas *c_temp_pim_sig = new TCanvas("c_final_sig_pim", "c_final_sig_pim", 800, 800);
  c_temp_pim_sig->Divide(4,2);

  std::cout << " size " << h_beta_mean_pim.size() << std::endl;
  for(int s = 0; s < 6; s++ ){
    c_temp_pim->cd(s+1);

    /// %%%%% MEAN FIT %%%%%%%%%
    TF1 *fit_mean_pim = new TF1(Form("fit_mean_pim_s%d",s),"[0] + [1]*x + [2]*x*x", 0.2, 2.5);///sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.0, 4.0);
    fit_mean_pim->SetParameter(0,0.0);
    fit_mean_pim->SetParameter(1,0.0);
    fit_mean_pim->SetParameter(2,0.0);//pim_mass);
    //fit_mean_pim->SetParameter(3,0.50);
    //fit_mean_pim->SetParameter(4,0.50);

    h_beta_mean_pim[s]->Fit(Form("fit_mean_pim_s%d",s),"R");

    h_beta_mean_pim[s]->SetTitle(Form("PIM BETA MEAN FOR SECTOR %d", s));
    h_beta_mean_pim[s]->SetMarkerColor(kBlack);
    h_beta_mean_pim[s]->SetMarkerStyle(8);
    h_beta_mean_pim[s]->Draw("AP");      

    fit_mean_pim->SetLineColor(kRed);
    fit_mean_pim->Draw("same"); 

    c_temp_pim_sig->cd(s+1);
    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_pim = new TF1(Form("fit_sigma_pim_s%d",s),"[0] + [1]*x + [2]*x*x", 0.2, 2.5); ///x + [2]*exp(-x)", 0.1, 2.5);
    fit_sigma_pim->SetParameter(0,0.0);
    fit_sigma_pim->SetParameter(1,0.0);
    fit_sigma_pim->SetParameter(2,0.0);//pim_mass);
    fit_sigma_pim->SetParameter(3,0.0);
    //fit_sigma_pim->SetParameter(4,0.10);
    
    h_beta_sig_pim[s]->Fit(Form("fit_sigma_pim_s%d",s),"R");  

    h_beta_sig_pim[s]->SetTitle(Form("PIM BETA SIGMA FOR SECTOR %d", s));
    h_beta_sig_pim[s]->SetMarkerColor(kBlack);
    h_beta_sig_pim[s]->SetMarkerStyle(8);
    h_beta_sig_pim[s]->Draw("AP");      
   
    fit_sigma_pim->SetLineColor(kRed);
    fit_sigma_pim->SetLineStyle(3);
    fit_sigma_pim->Draw("same"); 

    //store fit parameters in string to write to text file
    for( int i = 0; i < par_names.size(); i++ ){
      pimBetaMeanValues[i] = pimBetaMeanValues[i] + std::to_string(fit_mean_pim->GetParameter(i)) + "\t";
    }
    for( int i = 0; i < par_names_sig.size(); i++ ){    
      pimBetaSigValues[i] = pimBetaSigValues[i] + std::to_string(fit_sigma_pim->GetParameter(i)) + "\t";
    }


  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Plot KAON MINUS Beta vs P
  TCanvas *c_temp_km = new TCanvas("c_final_km", "c_final_km", 800, 800);
  c_temp_km->Divide(4,2);

  TCanvas *c_temp_km_sig = new TCanvas("c_final_sig_km", "c_final_sig_km", 800, 800);
  c_temp_km_sig->Divide(4,2);

  std::cout << " size " << h_beta_mean_km.size() << std::endl;
  for(int s = 0; s < 6; s++ ){
    c_temp_km->cd(s+1);

    /// %%%%% MEAN FIT %%%%%%%%%
    TF1 *fit_mean_km = new TF1(Form("fit_mean_km_s%d",s),"1.0/sqrt(1.0 + (0.493*0.493)/(x*x)) + [0] + [1]*x + [2]*x*x", 0.5, 3.5);// 0.31, 1.05);///sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.0, 4.0);
    fit_mean_km->SetParameter(0,0.0);
    fit_mean_km->SetParameter(1,0.0);
    fit_mean_km->SetParameter(2,0.0);//km_mass);
    //fit_mean_km->SetParameter(3,0.0);
    //fit_mean_km->SetParameter(4,0.10);   

    h_beta_mean_km[s]->Fit(Form("fit_mean_km_s%d",s),"R");

    h_beta_mean_km[s]->SetTitle(Form("KM BETA MEAN FOR SECTOR %d", s));
    h_beta_mean_km[s]->SetMarkerColor(kBlack);
    h_beta_mean_km[s]->SetMarkerStyle(8);
    h_beta_mean_km[s]->Draw("AP");      

    fit_mean_km->SetLineColor(kRed);
    fit_mean_km->Draw("same"); 

    c_temp_km_sig->cd(s+1);
    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_km = new TF1(Form("fit_sigma_km_s%d",s),"[0] + [1]*x + [2]*x*x", 0.44, 1.1); ///x + [2]*exp(-x)", 0.1, 2.5);
    fit_sigma_km->SetParameter(0,0.0);
    fit_sigma_km->SetParameter(1,0.0);
    fit_sigma_km->SetParameter(2,0.0);//km_mass);
    //fit_sigma_km->SetParameter(3,0.0);
    //fit_sigma_km->SetParameter(4,0.10);
    
    h_beta_sig_km[s]->Fit(Form("fit_sigma_km_s%d",s),"R");  

    h_beta_sig_km[s]->SetTitle(Form("KM BETA SIGMA FOR SECTOR %d", s));
    h_beta_sig_km[s]->SetMarkerColor(kBlack);
    h_beta_sig_km[s]->SetMarkerStyle(8);
    h_beta_sig_km[s]->Draw("AP");      
   
    fit_sigma_km->SetLineColor(kRed);
    fit_sigma_km->SetLineStyle(3);
    fit_sigma_km->Draw("same");

    //store fit parameters in string to write to text file
    for( int i = 0; i < par_names.size(); i++ ){
      kmBetaMeanValues[i] = kmBetaMeanValues[i] + std::to_string(fit_mean_km->GetParameter(i)) + "\t";
    }
    for( int i = 0; i < par_names_sig.size(); i++ ){    
      kmBetaSigValues[i] = kmBetaSigValues[i] + std::to_string(fit_sigma_km->GetParameter(i)) + "\t";
    }
 
  }

  

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Write out fit parameters to text files
  std::string parentDirectory = "./";
  ofstream protMeanOut;
  std::string protonMeanName = parentDirectory+"protonBetaFitParametersCLARY.txt";
  protMeanOut.open(protonMeanName);
  for( int l = 0; l < protonBetaMeanValues.size(); l++ ){
    protMeanOut << protonBetaMeanValues[l] << std::endl;
  }

  ofstream kpMeanOut;
  std::string kpMeanName = parentDirectory+"kpMeanBetaFitParametersCLARY.txt";
  kpMeanOut.open(kpMeanName);
  for( int l = 0; l < kpBetaMeanValues.size(); l++ ){
    kpMeanOut << kpBetaMeanValues[l] << std::endl;
  }

  ofstream kmMeanOut;
  std::string kmMeanName = parentDirectory+"kmMeanBetaFitParametersCLARY.txt";
  kmMeanOut.open(kmMeanName);
  for( int l = 0; l < kmBetaMeanValues.size(); l++ ){
    kmMeanOut << kmBetaMeanValues[l] << std::endl;
  }

  
  ofstream pipMeanOut;
  std::string pipMeanName = parentDirectory+"pipMeanBetaFitParametersCLARY.txt";
  pipMeanOut.open(pipMeanName);
  for( int l = 0; l < pipBetaMeanValues.size(); l++ ){
    pipMeanOut << pipBetaMeanValues[l] << std::endl;
  }

  ofstream pimMeanOut;
  std::string pimMeanName = parentDirectory+"pimMeanBetaFitParametersCLARY.t";
  pimMeanOut.open(pimMeanName);
  for( int l = 0; l < pimBetaMeanValues.size(); l++ ){
    pimMeanOut << pimBetaMeanValues[l] << std::endl;
  }


  ///////////////////////////////

  ofstream protSigOut;
  std::string protonSigName = parentDirectory+"protonSigmaBetaFitPrametersCLARY.txt";
  protSigOut.open(protonSigName);
  for( int l = 0; l < protonBetaSigValues.size(); l++ ){
    protSigOut << protonBetaSigValues[l] << std::endl;
  }

  ofstream kpSigOut;
  std::string kpSigName = parentDirectory+"kpSigmaBetaFitParametersCLARY.txt";
  kpSigOut.open(kpSigName);
  for( int l = 0; l < kpBetaSigValues.size(); l++ ){
    kpSigOut << kpBetaSigValues[l] << std::endl;
  }

  ofstream kmSigOut;
  std::string kmSigName = parentDirectory+"kmSigmaBetaFitParametersCLARY.txt";
  kmSigOut.open(kmSigName);
  for( int l = 0; l < kmBetaSigValues.size(); l++ ){
    kmSigOut << kmBetaSigValues[l] << std::endl;
  }

  
  ofstream pipSigOut;
  std::string pipSigName = parentDirectory+"pipSigmaBetaFitParametersCLARY.txt";
  pipSigOut.open(pipSigName);
  for( int l = 0; l < pipBetaSigValues.size(); l++ ){
    pipSigOut << pipBetaSigValues[l] << std::endl;
  }

  ofstream pimSigOut;
  std::string pimSigName = parentDirectory+"pimSigmaBetaFitParametersCLARY.txt";
  pimSigOut.open(pimSigName);
  for( int l = 0; l < pimBetaSigValues.size(); l++ ){
    pimSigOut << pimBetaSigValues[l] << std::endl;
  }


  return 0;
}
