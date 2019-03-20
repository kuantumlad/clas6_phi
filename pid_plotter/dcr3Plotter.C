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


int dcr3Plotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector< TH2D* > > m_h_dcr3_hit;
  
  for( int c = 0; c < 11; c++ ){
    std::vector< TH2D* > h_dcr3_hit;
    for( int s = 0; s < 6; s++ ){    
      h_dcr3_hit.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_cut%d_dcr3",s,c)) );    
    }
    m_h_dcr3_hit[c]=h_dcr3_hit;
  }

  
  for( int c = 0; c < 1; c++ ){
    for( int s = 0; s < 1; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,600, 400);
      m_h_dcr3_hit[0][0]->SetTitle(Form("DCR3: Sector %d",s));
      m_h_dcr3_hit[0][0]->GetXaxis()->SetTitle("Xdc(cm)");
      m_h_dcr3_hit[0][0]->GetYaxis()->SetTitle("Ydc(cm)");
      m_h_dcr3_hit[0][0]->GetYaxis()->SetTitleOffset(1.25);
      m_h_dcr3_hit[c][s]->Draw("colz");
    }
  }

 

  return 0;
}
