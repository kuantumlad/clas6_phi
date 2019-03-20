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
#include "programFiles/getGenIndices.C"
#include "MomCorr.C"
#include "CalculatorTools.C"


// change the size of the beta vs p histogram range to momentum of 2.2 (smaller than now and then refit) results will be better


int ecSamplingFractionFit( int filestart = 1, int Nfiles = 1){

  TStopwatch *stopwatch = new TStopwatch();
  
  TFile *fIn = new TFile("ec_sf_out.root");  

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


  std::vector<TH2D*> h_el_ectotp;
  std::vector<TGraphErrors*> g_el_ectotp;
  std::vector<TGraphErrors*> g_el_ectotp_sigma;

  std::vector<TCanvas*> v_can;
  std::vector<TCanvas*> v_can_sigma;
  
  for( int s = 0; s < 6; s++ ){
    v_can.push_back( new TCanvas(Form("can_%d",s),Form("can_%d",s),900,900) );
    v_can_sigma.push_back( new TCanvas(Form("can_sigma_%d",s),Form("can_sigma_%d",s),900,900) );
    h_el_ectotp.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_etotp_pc",s)));
  }

  for( int i = 0; i < h_el_ectotp.size(); i++ ){

    int nbins_x = h_el_ectotp[i]->GetNbinsX();
    int ndiv = 1;

    int nfits = nbins_x/ndiv;

    int start_bin_x = 20;
    int end_bin_x = start_bin_x + ndiv;

    std::vector<double> temp_center;
    std::vector<double> temp_fit_means;
    std::vector<double> temp_fit_sigmas;
    std::vector<double> temp_center_err;
    std::vector<double> temp_fit_means_err;
    std::vector<double> temp_fit_sigmas_err;

    for( int j = 0; j < nfits; j++ ){

      double sf_min = 0;
      double sf_max = 0;
      
      TH1D *h_ecsf = h_el_ectotp[i]->ProjectionY(Form("el_ecsf_projY_s%d_%d",i,j), start_bin_x, end_bin_x );
      TF1 *fit_ecsf = fitHistogram(h_ecsf);

      double current_p = h_el_ectotp[i]->GetXaxis()->GetBinCenter(start_bin_x);

      double min_p = h_el_ectotp[i]->GetXaxis()->GetBinCenter(start_bin_x) - h_el_ectotp[i]->GetXaxis()->GetBinWidth(start_bin_x)/2.0;
      double max_p = h_el_ectotp[i]->GetXaxis()->GetBinCenter(end_bin_x) + h_el_ectotp[i]->GetXaxis()->GetBinWidth(end_bin_x)/2.0;

      double slice_center = (max_p + min_p) / 2.0;
      start_bin_x+=ndiv;
      end_bin_x+=ndiv;

      double mean = fit_ecsf->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_ecsf->GetParameter(2);
      double sig_error = fit_ecsf->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);

      temp_center.push_back(slice_center);
      temp_center_err.push_back(0.0);

      temp_fit_means_err.push_back(mean_err);
      temp_fit_sigmas_err.push_back(sig_err);
    
      temp_fit_means.push_back(mean);
      temp_fit_sigmas.push_back(sig);

    }

    g_el_ectotp.push_back( new TGraphErrors( temp_center.size(), &(temp_center[0]), &(temp_fit_means[0]), &(temp_center_err[0]), &(temp_fit_means_err[0]) ) );
    g_el_ectotp_sigma.push_back( new TGraphErrors( temp_center.size(), &(temp_center[0]), &(temp_fit_sigmas[0]), &(temp_center_err[0]), &(temp_fit_sigmas_err[0]) ) );

  }

  std::vector<std::string> par_names = {"EC_SF_MU_A","EC_SF_MU_B","EC_SF_MU_C","EC_SF_MU_D"};
  std::vector<std::string> par_names_sig = {"EC_SF_SIGMA_A","EC_SF_SIGMA_B","EC_SF_SIGMA_C","EC_SF_MSIGMA_D"};
  

  std::vector<std::string> ec_sf_meanValues;
  std::vector<std::string> ec_sf_sigmaValues;

  for( int i = 0; i < par_names.size(); i++ ){
    ec_sf_meanValues.push_back(par_names[i]+" \t" + "6" + "\t" + "v" + "\t");
    ec_sf_sigmaValues.push_back(par_names_sig[i]+" \t" + "6" + "\t" + "v" + "\t");
  }

  for( int s = 0; s < 6; s++ ){
    v_can[s]->cd();

    TF1 *fit_means_sf = new TF1(Form("fit_means_sf_s%d",s),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0.51, 4.5);

    fit_means_sf->SetParameter(0,0.0);
    fit_means_sf->SetParameter(1,0.0);
    fit_means_sf->SetParameter(2,0.0);
    fit_means_sf->SetParameter(3,0.0);

    g_el_ectotp[s]->Fit(Form("fit_means_sf_s%d",s),"R");

    g_el_ectotp[s]->SetTitle(Form("EC  MEAN FIT SECTOR %d",s));
    g_el_ectotp[s]->SetMarkerColor(kBlack);
    g_el_ectotp[s]->SetMarkerStyle(8);
    g_el_ectotp[s]->Draw("AP");
    
    fit_means_sf->SetLineColor(kRed);
    fit_means_sf->Draw("same");

    for( int i = 0; i < fit_means_sf->GetNpar(); i++ ){
      ec_sf_meanValues[i] = ec_sf_meanValues[i] + " " + std::to_string(fit_means_sf->GetParameter(i)); + "\t";
    }
  }
    

  for( int s = 0; s < 6; s++ ){
    v_can_sigma[s]->cd();

    TF1 *fit_sigmas_sf = new TF1(Form("fit_sigmas_sf_s%d",s),"[0] + [1]*x + [2]*x*x + [3]*x*x*x", 0.51, 4.5);

    fit_sigmas_sf->SetParameter(0,0.0);
    fit_sigmas_sf->SetParameter(1,0.0);
    fit_sigmas_sf->SetParameter(2,0.0);
    fit_sigmas_sf->SetParameter(3,0.0);

    g_el_ectotp_sigma[s]->Fit(Form("fit_sigmas_sf_s%d",s),"R");

    g_el_ectotp_sigma[s]->SetTitle(Form("EC SIGMA FIT SECTOR %d",s));
    g_el_ectotp_sigma[s]->SetMarkerColor(kBlack);
    g_el_ectotp_sigma[s]->SetMarkerStyle(8);
    g_el_ectotp_sigma[s]->Draw("AP");
    
    fit_sigmas_sf->SetLineColor(kRed);
    fit_sigmas_sf->Draw("same");

    for( int i = 0; i < fit_sigmas_sf->GetNpar(); i++ ){
      ec_sf_sigmaValues[i] = ec_sf_sigmaValues[i] + " " + std::to_string(fit_sigmas_sf->GetParameter(i)); + "\t";
    }
  }
    
  ////////////////////////////////////////////////////////////////////////////////
  // print out items
  std::string parentDirectory = "./";
  ofstream ecsfMeanOut;
  std::string ecsfMeanName = parentDirectory+"elECSFMeanFitParametersCLARY.txt";
  ecsfMeanOut.open(ecsfMeanName);
  for( int l = 0; l < ec_sf_meanValues.size(); l++ ){
    ecsfMeanOut << ec_sf_meanValues[l] << std::endl;
  }

  ofstream ecsfSigmaOut;
  std::string ecsfSigmaName = parentDirectory+"elECSFSigmaFitParametersCLARY.txt";
  ecsfSigmaOut.open(ecsfSigmaName);
  for( int l = 0; l < ec_sf_sigmaValues.size(); l++ ){
    ecsfSigmaOut << ec_sf_sigmaValues[l] << std::endl;
  }


  return 0;
}
