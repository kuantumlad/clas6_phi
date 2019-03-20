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


int pPlotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector<TH1D* > > m_h_p;
  
  for( int c = 0; c < 11; c++ ){
    std::vector< TH1D* > h_p;
    for( int s = 0; s <= 6; s++ ){    
      h_p.push_back( (TH1D*)fIn->Get(Form("h_el_s%d_cut%d_p",s,c)) );    
    }
    m_h_p[c]=h_p;
  }

  
  for( int c = 0; c < 1; c++ ){
    for( int s = 0; s < 1; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,600, 400);
      m_h_p[c][s]->SetTitle(Form("Electron Momentum: Sector %d",s));
      m_h_p[c][s]->GetXaxis()->SetTitle("p(GeV)");
      m_h_p[c][s]->GetYaxis()->SetTitle("Number of Events");
      m_h_p[c][s]->GetYaxis()->SetTitleOffset(1.45);
      m_h_p[c][s]->Draw();
    }
  }

  return 0;
}
