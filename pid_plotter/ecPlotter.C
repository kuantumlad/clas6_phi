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


int ecPlotter(const char* input){

  TFile *fIn = new TFile(input,"");

  std::map< int, std::vector< TH2D* > > m_h_ec_hit;
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
    std::vector< TH2D* > h_ec_hit;
    for( int s = 0; s <= 6; s++ ){    
      h_ec_hit.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_cut%d_ec_pos",s,c)) );    
    }
    m_h_ec_hit[c]=h_ec_hit;
  }

  

  

  for( int c = 0; c < 12; c++ ){
    for( int s = 0; s <=6 ; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,900, 900);

      
      m_h_ec_hit[0][s]->SetTitle("EC Hit Position");
      m_h_ec_hit[0][s]->GetXaxis()->SetTitle("Xec(cm)");
      m_h_ec_hit[0][s]->GetYaxis()->SetTitle("Yec(cm)");
      m_h_ec_hit[0][s]->GetYaxis()->SetTitleOffset(1.25);
      m_h_ec_hit[0][s]->SetMarkerColor(kRed);
      
      m_h_ec_hit[0][s]->Draw();       
      m_h_ec_hit[c][s]->Draw("same");   
      if( s == 0 ){
	can->SaveAs(Form("img/h_ec_fid_%s_all.png",cut_names[c].c_str()));
      }
      else{
	can->SaveAs(Form("img/h_ec_fid_%s_s%d.png",cut_names[c].c_str(), s));
      }
    }
  }

  


  return 0;
}
