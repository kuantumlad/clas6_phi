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


int eieoPlotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector< TH2D* > > m_h_eieo;
  
  for( int c = 0; c < 11; c++ ){
    std::vector< TH2D* > h_eieo;
    for( int s = 0; s <= 6; s++ ){    
      h_eieo.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_cut%d_ecei_eceo",s,c)) );    
    }
    m_h_eieo[c]=h_eieo;
  }

  
  for( int c = 0; c < 1; c++ ){
    for( int s = 1; s < 2; s++ ){ //start at 1 for sector 1
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,600, 400);
      gPad->SetLogz();
      m_h_eieo[c][s]->SetTitle(Form("Energy Inner vs Outer: Sector %d",s));
      m_h_eieo[c][s]->GetXaxis()->SetTitle("Ein(GeV)");
      m_h_eieo[c][s]->GetYaxis()->SetTitle("Eout(GeV)");
      m_h_eieo[c][s]->GetYaxis()->SetTitleOffset(1.25);
      m_h_eieo[c][s]->Draw("colz");
    }
  }

 

  return 0;
}
