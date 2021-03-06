#include <fstream>
#include <iostream>
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "programFiles/functions.C"
#include "programFiles/eID.C"
#include "programFiles/getGenIndices.C"
#include "MomCorr.C"

int ecSamplingFractionFits( int filestart = 1, int Nfiles = 1){

  TStopwatch *stopwatch = new TStopwatch();
  
  TChain *h22chain = new TChain("h22");  
  string firstfilename = "";
  string lastfilename = "";


  // %%%%%%%% read the input files into the TChain %%%%%%%%
  int NtotalFiles = 11625;  
  ifstream filelist;
  filelist.open("programFiles/dataFiles.txt");
  int kStop = Nfiles + filestart;
  if(kStop > NtotalFiles+1) kStop = NtotalFiles+1;
  for(int k = 1; k < kStop; k++){
    string filename;
    filelist>>filename;
    std::cout << " >> ADDING FILE " << filename.c_str() << std::endl;
    if(k >= filestart) h22chain->Add(filename.c_str());
    if(k == filestart) firstfilename = filename;
    if(k == kStop - 1) lastfilename = filename;
  }
  
  firstfilename = firstfilename.substr(39,9); // trim the string down so it's just the run#.subrun# (e.g. 38458.a00)
  lastfilename = lastfilename.substr(39,9);
  cout<<"first file: "<<firstfilename<<" ... last file: "<<lastfilename<<endl;
  // %%%%% end read the input files into a TChain %%%%%

  Float_t e_mass = 0.000511; // GeV
  Float_t prot_mass = 0.938272; // GeV
  Float_t pip_mass = 0.13957; // GeV
  Float_t pim_mass = 0.13957; // GeV
  Float_t kp_mass = 0.4937; // GeV
  Float_t km_mass = 0.4937; // GeV
  Float_t speed_of_light = 29.9792458; // cm/ns
  Float_t Beam_Energy = 5.498; // GeV
  Float_t pi = 3.14159265359;
  Float_t pi180 = pi/180.0;
  Float_t pi180_inv = 180.0/pi;
  int ExpOrSim = 1;

  int e_zvertex_strict = 0;
  int e_ECsampling_strict = 0;
  int e_ECoutVin_strict = 0;
  int e_CCthetaMatching_strict = 0;
  int e_ECgeometric_strict = 0;
  int e_R1fid_strict = 0;
  int e_R3fid_strict = 0;
  int e_CCphiMatching_strict = 0;
  int e_CCfiducial_strict = 0;
  int do_momCorr_e = 0;
  int e_CCnphe_strict = 0;

  MomCorr_e1f *MomCorr = new MomCorr_e1f();

  Int_t gpart;
  Int_t mcentr;
  Int_t mcid;
  Float_t mctheta;
  Float_t mcphi;
  Float_t mcp;
    
  Float_t p[35];
  Int_t q[35];
  Float_t cx[35];
  Float_t cy[35];
  Float_t cz[35];
  Float_t vz[35];
  Float_t vy[35];
  Float_t vx[35];
  Float_t b[35];
    
  UChar_t dc_sect[35];
  Float_t tl1_x[35];
  Float_t tl1_y[35];
  Float_t tl3_x[35]; // used to be dc_xsc in h10
  Float_t tl3_y[35];
  Float_t tl3_z[35];
  Float_t tl3_cx[35]; // used to be dc_cxsc in h10
  Float_t tl3_cy[35];
  Float_t tl3_cz[35];
    
  UChar_t ec_sect[35];
  Float_t etot[35];
  Float_t ec_ei[35];
  Float_t ec_eo[35];
  Float_t ec_t[35];
  Float_t ech_x[35];
  Float_t ech_y[35];
  Float_t ech_z[35];
    
  UChar_t sc_sect[35];
  Float_t sc_t[35];
  Float_t sc_r[35];
  UChar_t sc_pd[35];
    
  UChar_t cc_sect[35];
  UShort_t cc_segm[35];
  UShort_t nphe[35];
  // %%%%% end define variables %%%%%
    
  // %%%%% set branch addresses %%%%%    
  h22chain->SetBranchAddress("gpart", &gpart);
  h22chain->SetBranchAddress("p", p);
  h22chain->SetBranchAddress("q", q);
  h22chain->SetBranchAddress("cx", cx);
  h22chain->SetBranchAddress("cy", cy);
  h22chain->SetBranchAddress("cz", cz);
  h22chain->SetBranchAddress("vz", vz);
  h22chain->SetBranchAddress("vy", vy);
  h22chain->SetBranchAddress("vx", vx);
  h22chain->SetBranchAddress("b", b);
    
  h22chain->SetBranchAddress("dc_sect", dc_sect);
  h22chain->SetBranchAddress("tl1_x", tl1_x);
  h22chain->SetBranchAddress("tl1_y", tl1_y);
  h22chain->SetBranchAddress("tl3_x", tl3_x);
  h22chain->SetBranchAddress("tl3_y", tl3_y);
  h22chain->SetBranchAddress("tl3_z", tl3_z);
  h22chain->SetBranchAddress("tl3_cx", tl3_cx);
  h22chain->SetBranchAddress("tl3_cy", tl3_cy);
  h22chain->SetBranchAddress("tl3_cz", tl3_cz);
    
  h22chain->SetBranchAddress("ec_sect", ec_sect);
  h22chain->SetBranchAddress("etot", etot);
  h22chain->SetBranchAddress("ec_ei", ec_ei);
  h22chain->SetBranchAddress("ec_eo", ec_eo);
  h22chain->SetBranchAddress("ec_t", ec_t);
  h22chain->SetBranchAddress("ech_x", ech_x);
  h22chain->SetBranchAddress("ech_y", ech_y);
  h22chain->SetBranchAddress("ech_z", ech_z);
    
  h22chain->SetBranchAddress("sc_sect", sc_sect);
  h22chain->SetBranchAddress("sc_t", sc_t);
  h22chain->SetBranchAddress("sc_r", sc_r);
  h22chain->SetBranchAddress("sc_pd", sc_pd);
    
  h22chain->SetBranchAddress("cc_sect", cc_sect);
  h22chain->SetBranchAddress("cc_segm", cc_segm);
  h22chain->SetBranchAddress("nphe", nphe);
  // %%%%% end set branch addresses %%%%%
    
  // %%%%% stuff for output root files %%%%%
  string outfilename;
  TFile *outputfile;
  outputfile = new TFile("ec_sf_out.root","RECREATE");

  cout<<"entries: "<<h22chain->GetEntries()/1000<<" thousand"<<endl;

  TLorentzVector V4k(0.0, 0.0, Beam_Energy, Beam_Energy); // x, y, z, t
  TLorentzVector V4ISproton(0.0, 0.0, 0.0, prot_mass); // IS = Initial State

  int el_pass_rate[10] = {0,0,0,0,0,0,0,0,0};

  int limit = h22chain->GetEntries() / 100;
  std::cout << " LIMIT " << limit << std::endl;


  // histograms
  TH2D *h_el_etotp_pc = new TH2D("h_el_etot_p","h_el_etotp",200, 0.0, Beam_Energy, 200, 0.0, 0.50);
  std::vector<TH2D*> h_el_sect_etotp_pc;  
  for( int i = 0; i <=6; i++){
    h_el_sect_etotp_pc.push_back( new TH2D(Form("h_el_s%d_etotp_pc",i),Form("h_el_s%d_etotp_pc",i), 200, 0.0, Beam_Energy, 200, 0.0, 0.50) );
  }
  //end histograms


  for(int i = 0; i < h22chain->GetEntries(); i++){    
    //std::cout << i << std::endl;
    
    h22chain->GetEntry(i);

    if ( i % limit == 0 ){
      double completed = (double)i / (double)h22chain->GetEntries();
      std::cout << " >> COMPLETED " << ceil(completed * 100) << " % --- ON EVENT " << i << std::endl;
    }


    string currentrunno_string = "-123";
    if(ExpOrSim == 1)
      {
	currentrunno_string = h22chain->GetCurrentFile()->GetName();
	currentrunno_string = currentrunno_string.substr(39,5);
    }
    int currentrunno = atoi(currentrunno_string.c_str()); // converts the string to an int
    
    //std::cout << " CURRENT RUN " << currentrunno << std::endl;
    //std::cout << " run string " << currentrunno_string.c_str() << std::endl;

    TVector3 V3_e[2]; // 2 for gen(0) and rec(1)
    TLorentzVector V4_e[2];
    
    // %%%%% electron ID %%%%%
    int e_index[2] = {-123,-123};

    e_index[1] = eID(gpart, q, p, cc_sect, sc_sect, ec_sect, dc_sect, cx, cy, cz, tl1_x, tl1_y, tl3_x, tl3_y, tl3_z, tl3_cx, tl3_cy, tl3_cz, e_zvertex_strict, vz, vy, vx, e_ECsampling_strict, ExpOrSim, etot, e_ECoutVin_strict, ec_ei, ech_x, ech_y, ech_z, e_CCthetaMatching_strict, cc_segm, e_ECgeometric_strict, e_R1fid_strict, e_R3fid_strict, e_CCphiMatching_strict, sc_pd, e_CCfiducial_strict, e_CCnphe_strict, nphe, el_pass_rate);
    //std::cout << " finished electron index is:  " << index[1] << std::endl;
    //%%%%% end electron ID %%%%%    
    
    if(e_index[1] >= 0 ){
      //std::cout << " final electron index " << e_index[1] << std::endl;
      V3_e[1].SetXYZ(p[e_index[1]]*cx[e_index[1]], p[e_index[1]]*cy[e_index[1]], p[e_index[1]]*cz[e_index[1]]);
      V4_e[1].SetXYZT(V3_e[1].X(), V3_e[1].Y(), V3_e[1].Z(), sqrt(V3_e[1].Mag2() + pow(e_mass,2)));
      //V4_e_uncorr[1].SetXYZT(V3_e[1].X(), V3_e[1].Y(), V3_e[1].Z(), sqrt(V3_e[1].Mag2() + pow(e_mass,2)));
      //std::cout << " electron momentum before MOMCORR " << V4_e[1].Px() << "  " << V4_e[1].Py() << " " << V4_e[1].Pz() << " " << V4_e[1].E() << std::endl;
      if(do_momCorr_e && ExpOrSim == 1) V4_e[1] = MomCorr->PcorN(V4_e[1], -1, 11);
      
      double ectot = ec_ei[e_index[1]] + ec_eo[e_index[1]];
      int ec_sector = ec_sect[e_index[1]]-1;

      h_el_etotp_pc->Fill( V4_e[1].P(), ectot/V4_e[1].P() );     
      if( ec_sector >= 0 ){
	h_el_sect_etotp_pc[ec_sector]->Fill( V4_e[1].P(), ectot/V4_e[1].P() );      
      }

    }

  }

  outputfile->Write();
  return 0;
}
