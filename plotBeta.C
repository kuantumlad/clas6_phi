#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TH2D.h>
#include <iostream>
#include <vector>
#include <ostream>
#include <fstream>
#include "loadBetaParameters.C"

int plotBeta(){


  std::map<int, std::vector<double> > pr_mean_fit = loadBetaParameters("protonBetaFitParametersCLARY.txt");
  std::map< int , std::vector<double> >::iterator it;
  for( it = pr_mean_fit.begin(); it != pr_mean_fit.end(); it++ ){
    std::cout << " FIT PARAMETER " << it->first << std::endl;
    std::cout << "SIZE OF VECTOR " << (it->second).size()  << std::endl;
    for( std::vector<double>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2){
      std::cout << " SECTOR VALUES " << *it2 << std::endl;
    }

  }

  std::map<int, std::vector<double> > pip_mean_fit = loadBetaParameters("pipMeanBetaFitParametersCLARY.txt");
  std::map<int, std::vector<double> > kp_mean_fit = loadBetaParameters("kpMeanBetaFitParametersCLARY.txt");//beta_fits/kpMeanBetaFitParameters.txt");//kpMeanBetaFitParametersCLARY.txt");
  std::map<int, std::vector<double> > pim_mean_fit = loadBetaParameters("pimMeanBetaFitParametersCLARY.txt");
  std::map<int, std::vector<double> > km_mean_fit = loadBetaParameters("kmMeanBetaFitParametersCLARY.txt");

  std::map<int, std::vector<double> > pr_sigma_fit = loadBetaParameters("protonSigmaBetaFitPrametersCLARY.txt");//protonSigmaBetaFitParametersCLARY.txt");
  std::map<int, std::vector<double> > pip_sigma_fit = loadBetaParameters("pipSigmaBetaFitParametersCLARY.txt");
  std::map<int, std::vector<double> > kp_sigma_fit = loadBetaParameters("kpSigmaBetaFitParametersCLARY.txt");
  std::map<int, std::vector<double> > km_sigma_fit = loadBetaParameters("kmSigmaBetaFitParametersCLARY.txt");

  std::vector<TH2D*> h2_pr_betap;
  std::vector<TH2D*> h2_pip_betap;
  std::vector<TH2D*> h2_kp_betap;
  std::vector<TH2D*> h2_pim_betap;
  std::vector<TH2D*> h2_km_betap;

  for( int s = 0; s < 6; s++ ){

    h2_pr_betap.push_back( new TH2D(Form("h2_pr_betap_s%d",s),Form("h2_pr_betap_s%d",s), 100, 0.0, 3.0, 100, 0.0, 1.2));
    h2_pip_betap.push_back( new TH2D(Form("h2_pip_betap_s%d",s),Form("h2_pip_betap_s%d",s), 100, 0.0, 3.0, 100, 0.0, 1.2));
    h2_kp_betap.push_back( new TH2D(Form("h2_kp_betap_s%d",s),Form("h2_kp_betap_s%d",s), 100, 0.0, 3.0, 100, 0.0, 1.2));
    h2_pim_betap.push_back( new TH2D(Form("h2_pim_betap_s%d",s),Form("h2_pim_betap_s%d",s), 100, 0.0, 3.0, 100, 0.0, 1.2));
    h2_km_betap.push_back( new TH2D(Form("h2_km_betap_s%d",s),Form("h2_km_betap_s%d",s), 100, 0.0, 3.0, 100, 0.0, 1.2));
    
  }
  
  std::vector<double> pr_meanA = pr_mean_fit[0];
  std::vector<double> pr_meanB = pr_mean_fit[1];
  std::vector<double> pr_meanC = pr_mean_fit[2];

  std::vector<double> pip_meanA = pip_mean_fit[0];
  std::vector<double> pip_meanB = pip_mean_fit[1];
  std::vector<double> pip_meanC = pip_mean_fit[2];

  std::vector<double> kp_meanA = kp_mean_fit[0];
  std::vector<double> kp_meanB = kp_mean_fit[1];
  std::vector<double> kp_meanC = kp_mean_fit[2];

  std::vector<double> pim_meanA = pim_mean_fit[0];
  std::vector<double> pim_meanB = pim_mean_fit[1];
  std::vector<double> pim_meanC = pim_mean_fit[2];

  std::vector<double> km_meanA = km_mean_fit[0];
  std::vector<double> km_meanB = km_mean_fit[1];
  std::vector<double> km_meanC = km_mean_fit[2];

  std::vector<double> pr_sigA = pr_sigma_fit[0];

  TCanvas *c_beta_pr = new TCanvas("c_beta_pr","c_beta_pr",900,900);
  c_beta_pr->Divide(3,2);
  TCanvas *c_beta_pip = new TCanvas("c_beta_pip","c_beta_pip",900,900);
  c_beta_pip->Divide(3,2);
  TCanvas *c_beta_kp = new TCanvas("c_beta_kp","c_beta_kp",900,900);
  c_beta_kp->Divide(3,2);
  TCanvas *c_beta_pim = new TCanvas("c_beta_pim","c_beta_pim",900,900);
  c_beta_pim->Divide(3,2);
  TCanvas *c_beta_km = new TCanvas("c_beta_km","c_beta_km",900,900);
  c_beta_km->Divide(3,2);



  for( int s = 0; s < 6; s++ ){
    
    //positive hadrons
    TF1 *fit_beta_pr = new TF1(Form("fit_pr_s%d",s),"1.0/(1 + [0]/x^2)",0.2,5.0);
    fit_beta_pr->SetParameter(0,2*0.938);
    TF1 *fit_pr = new TF1(Form("fit_pr_s%d",s),"1.0/(1 + [0]/x^2) + [1]*x*x + [2]*x + [3]",0.2,5.0);
    fit_pr->SetParameter(0,2*0.938);
    fit_pr->SetParameter(1,pr_meanC[s]);
    fit_pr->SetParameter(2,pr_meanB[s]);
    fit_pr->SetParameter(3,pr_meanA[s]);
    c_beta_pr->cd(s+1);
    fit_beta_pr->SetLineColor(kBlue);
    fit_beta_pr->Draw();
    fit_pr->Draw("same");

    TF1 *fit_beta_pip = new TF1(Form("fit_pip_s%d",s),"1.0/(1 + [0]/x^2)",0.2,5.0);
    fit_beta_pip->SetParameter(0,2*0.1395);
    TF1 *fit_pip = new TF1(Form("fit_pip_s%d",s),"1.0/(1 + [0]/x^2) + [1]*x*x + [2]*x + [3]",0.2,5.0);
    fit_pip->SetParameter(0,2*0.1395);
    fit_pip->SetParameter(1,pip_meanC[s]);
    fit_pip->SetParameter(2,pip_meanB[s]);
    fit_pip->SetParameter(3,pip_meanA[s]);
    c_beta_pip->cd(s+1);
    fit_beta_pip->SetLineColor(kBlue);
    fit_beta_pip->Draw();
    fit_pip->Draw("same");
    
    TF1 *fit_beta_kp = new TF1(Form("fit_kp_s%d",s),"1.0/(1 + [0]/x^2)",0.2,5.0);
    fit_beta_kp->SetParameter(0,2*0.493);
    TF1 *fit_kp = new TF1(Form("fit_kp_s%d",s),"1.0/(1 + [0]/x^2) + [1]*x*x + [2]*x + [3]",0.2,5.0);
    fit_kp->SetParameter(0,2*0.493);
    fit_kp->SetParameter(1,kp_meanC[s]);
    fit_kp->SetParameter(2,kp_meanB[s]);
    fit_kp->SetParameter(3,kp_meanA[s]);
    c_beta_kp->cd(s+1);
    fit_beta_kp->SetLineColor(kBlue);
    fit_beta_kp->Draw();
    fit_kp->Draw("same");

    // negative hadrons
    TF1 *fit_beta_pim = new TF1(Form("fit_pim_s%d",s),"1.0/(1 + [0]/x^2)",0.2,5.0);
    fit_beta_pim->SetParameter(0,2*0.1395);
    TF1 *fit_pim = new TF1(Form("fit_pim_s%d",s),"1.0/(1 + [0]/x^2) + [1]*x*x + [2]*x + [3]",0.2,5.0);
    fit_pim->SetParameter(0,2*0.1395);
    fit_pim->SetParameter(1,pim_meanC[s]);
    fit_pim->SetParameter(2,pim_meanB[s]);
    fit_pim->SetParameter(3,pim_meanA[s]);
    c_beta_pim->cd(s+1);
    fit_beta_pim->SetLineColor(kBlue);
    fit_beta_pim->Draw();
    fit_pim->Draw("same");
    
    TF1 *fit_beta_km = new TF1(Form("fit_km_s%d",s),"1.0/(1 + [0]/x^2)",0.2,5.0);
    fit_beta_km->SetParameter(0,2*0.493);
    TF1 *fit_km = new TF1(Form("fit_km_s%d",s),"1.0/(1 + [0]/x^2) + [1]*x*x + [2]*x + [3]",0.2,5.0);
    fit_km->SetParameter(0,2*0.493);
    fit_km->SetParameter(1,km_meanC[s]);
    fit_km->SetParameter(2,km_meanB[s]);
    fit_km->SetParameter(3,km_meanA[s]);
    c_beta_km->cd(s+1);
    fit_beta_km->SetLineColor(kBlue);
    fit_beta_km->Draw();
    fit_km->Draw("same");



  }



  return 0;
}
