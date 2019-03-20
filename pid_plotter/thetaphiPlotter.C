#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TF1.h>
#include <TAxis.h>

#include <vector>
#include <string>
#include <map>
#include <iostream>


int thetaphiPlotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector< TH2D* > > m_h_thetaphi;
  
  for( int c = 0; c < 11; c++ ){
    std::vector< TH2D* > h_thetaphi;
    for( int s = 0; s <= 6; s++ ){    
      h_thetaphi.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_cut%d_theta_phi",s,c)) );    
    }
    m_h_thetaphi[c]=h_thetaphi;
  }

  
  for( int c = 0; c < 1; c++ ){
    for( int s = 1; s < 2; s++ ){ //start at 1 for sector 1
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,900, 900);
      //gPad->SetLogz();
      //x axis is relative phi (shifted between -30 and 30 for each sector and theta on y-axis from CC
      m_h_thetaphi[c][s]->SetTitle(Form("Angular Distribution: Sector %d",s));
      m_h_thetaphi[c][s]->GetXaxis()->SetTitle("CC#theta(deg)");
      m_h_thetaphi[c][s]->GetYaxis()->SetTitle("#phi(deg)");
      m_h_thetaphi[c][s]->GetYaxis()->SetTitleOffset(1.25);
      m_h_thetaphi[c][s]->Draw("colz");
    }
  }

 

  return 0;
}
