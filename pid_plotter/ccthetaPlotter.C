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


int ccthetaPlotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector< TH2D* > > m_h_cctheta_hit;
  
  for( int c = 0; c < 11; c++ ){
    std::vector< TH2D* > h_cctheta_hit;
    for( int s = 0; s < 6; s++ ){    
      h_cctheta_hit.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_cut%d_cc_theta",s,c)) );    
    }
    m_h_cctheta_hit[c]=h_cctheta_hit;
  }

  
  for( int c = 0; c < 1; c++ ){
    for( int s = 0; s < 1; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,600, 400);
      m_h_cctheta_hit[c][s]->SetTitle(Form("PMT vs CC#theta: Sector %d",s));
      m_h_cctheta_hit[c][s]->GetXaxis()->SetTitle("PMT");
      m_h_cctheta_hit[c][s]->GetYaxis()->SetTitle("CC#theta(deg)");
      m_h_cctheta_hit[c][s]->GetYaxis()->SetTitleOffset(1.25);
      m_h_cctheta_hit[c][s]->Draw("colz");
    }
  }

 

  return 0;
}
