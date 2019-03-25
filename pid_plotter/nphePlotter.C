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


int nphePlotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector<TH1D* > > m_h_nphe;
  std::vector<std::string> cut_names = {"allNegatives",
					"Z_VERTEX",
					"CC_FID",
					"CC_PHI",
					"CC_THETA",
					"DC_R1_FID",
					"DC_R3_FID",
					"EC_FID",
					"EC_IN_OUT",
					"EC_SAMPLING",
					"NPHE",
					"cuts"};
  
    
  for( int c = 0; c < 12; c++ ){
    std::vector< TH1D* > h_nphe;
    for( int s = 0; s <= 6; s++ ){    
      h_nphe.push_back( (TH1D*)fIn->Get(Form("h_el_s%d_cut%d_ccnphe",s,c)) );    
    }
    m_h_nphe[c]=h_nphe;
  }

  
  for( int c = 0; c < 12; c++ ){
    for( int s = 0; s <= 6; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,900, 900);
      m_h_nphe[c][s]->SetTitle(Form("NPhe: Sector %d",s));
      m_h_nphe[c][s]->GetXaxis()->SetTitle("Nphe");
      m_h_nphe[c][s]->GetYaxis()->SetTitle("Number of Events");
      m_h_nphe[c][s]->GetYaxis()->SetTitleOffset(1.35);
      m_h_nphe[c][s]->Draw();


    if ( s == 0 ){
	can->SaveAs(Form("img/h1_nphe_%s_all.png",cut_names[c].c_str())); 
      }
      else{
	can->SaveAs(Form("img/h1_nphe_%s_s%d.png",cut_names[c].c_str(), s)); 
      }

    }
  }

  return 0;
}
