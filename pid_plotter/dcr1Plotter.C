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


int dcr1Plotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector< TH2D* > > m_h_dcr1_hit;
  
  for( int c = 0; c < 11; c++ ){
    std::vector< TH2D* > h_dcr1_hit;
    for( int s = 0; s < 6; s++ ){    
      h_dcr1_hit.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_cut%d_dcr1",s,c)) );    
    }
    m_h_dcr1_hit[c]=h_dcr1_hit;
  }
  
  m_h_dcr1_hit[0][0]->SetTitle("DCR1 Hit Position");
  m_h_dcr1_hit[0][0]->GetXaxis()->SetTitle("Xdc(cm)");
  m_h_dcr1_hit[0][0]->GetYaxis()->SetTitle("Ydc(cm)");
  m_h_dcr1_hit[0][0]->GetYaxis()->SetTitleOffset(1.25);
  m_h_dcr1_hit[0][0]->SetMarkerColor(kRed);

  for( int c = 1; c < 2; c++ ){
    TCanvas *can = new TCanvas(Form("can%d",c),Form("can%d",c) ,900, 900);
    m_h_dcr1_hit[0][0]->Draw();       
    m_h_dcr1_hit[c][0]->Draw("same");   
  }

  
  return 0;
}
