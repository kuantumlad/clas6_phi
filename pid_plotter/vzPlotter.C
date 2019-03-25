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
    std::vector< TH1D* > h_vz;
    std::vector< TH1D* > h_corr_vz;
    for( int s = 0; s <= 6; s++ ){    
      h_vz.push_back( (TH1D*)fIn->Get(Form("h_el_s%d_cut%d_vz",s,c)) );    
      h_corr_vz.push_back( (TH1D*)fIn->Get(Form("h_el_s%d_cut%d_corr_vz",s,c)) );    
    }
    m_h_vz[c]=h_vz;
    m_h_corr_vz[c]=h_corr_vz;
  }

  
  for( int c = 0; c < 12; c++ ){
    for( int s = 0; s <= 6; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,600, 400);
      m_h_vz[c][s]->SetTitle(Form("Electron Vertex: Sector %d",s));
      m_h_vz[c][s]->GetXaxis()->SetTitle("Zvertex(cm)");
      m_h_vz[c][s]->GetYaxis()->SetTitle("Number of Events");
      m_h_vz[c][s]->GetYaxis()->SetTitleOffset(1.45);
      m_h_vz[c][s]->Draw();
      m_h_corr_vz[c][s]->SetLineColor(kRed);
      m_h_corr_vz[c][s]->Draw("same");

      if ( s == 0 ){
	can->SaveAs(Form("img/h_vz_%s_all.png",cut_names[c].c_str())); 
      }
      else{
	can->SaveAs(Form("img/h_vz_%s_s%d.png",cut_names[c].c_str(), s)); 
      }

    }
  }

  return 0;
}
