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


int vzPlotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector<TH1D* > > m_h_vz;
  std::map< int, std::vector<TH1D* > > m_h_corr_vz;
  
  for( int c = 0; c < 11; c++ ){
    std::vector< TH1D* > h_vz;
    std::vector< TH1D* > h_corr_vz;
    for( int s = 0; s <= 6; s++ ){    
      h_vz.push_back( (TH1D*)fIn->Get(Form("h_el_s%d_cut%d_vz",s,c)) );    
      h_corr_vz.push_back( (TH1D*)fIn->Get(Form("h_el_s%d_cut%d_corr_vz",s,c)) );    
    }
    m_h_vz[c]=h_vz;
    m_h_corr_vz[c]=h_corr_vz;
  }

  
  for( int c = 0; c < 1; c++ ){
    for( int s = 0; s < 1; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,600, 400);
      m_h_vz[c][s]->SetTitle(Form("Electron Vertex: Sector %d",s));
      m_h_vz[c][s]->GetXaxis()->SetTitle("Zvertex(cm)");
      m_h_vz[c][s]->GetYaxis()->SetTitle("Number of Events");
      m_h_vz[c][s]->GetYaxis()->SetTitleOffset(1.45);
      m_h_vz[c][s]->Draw();
      m_h_corr_vz[c][s]->SetLineColor(kRed);
      m_h_corr_vz[c][s]->Draw("same");
    }
  }

  return 0;
}
