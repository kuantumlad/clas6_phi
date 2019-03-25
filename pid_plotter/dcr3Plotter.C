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
    std::vector< TH2D* > h_dcr3_hit;
    for( int s = 0; s <= 6; s++ ){    
      h_dcr3_hit.push_back( (TH2D*)fIn->Get(Form("h_el_s%d_cut%d_dcr3",s,c)) );    
    }
    m_h_dcr3_hit[c]=h_dcr3_hit;
  }

  
  for( int c = 0; c < 12; c++ ){
    for( int s = 0; s <= 6; s++ ){
      TCanvas *can = new TCanvas(Form("can%d_s%d",c,s),Form("can%d_s%d",c,s) ,900, 900);
      m_h_dcr3_hit[c][s]->SetTitle(Form("DCR3: Sector %d",s));
      m_h_dcr3_hit[c][s]->GetXaxis()->SetTitle("Xdc(cm)");
      m_h_dcr3_hit[c][s]->GetYaxis()->SetTitle("Ydc(cm)");
      m_h_dcr3_hit[c][s]->GetYaxis()->SetTitleOffset(1.25);
      m_h_dcr3_hit[c][s]->Draw("colz");

       if ( s == 0 ){
	can->SaveAs(Form("img/h_dcr3_fid_%s_all.png",cut_names[c].c_str())); 
      }
      else{
	can->SaveAs(Form("img/h_dcr3_fid_%s_s%d.png",cut_names[c].c_str(), s)); 
      }

    }
  }

 

  return 0;
}
