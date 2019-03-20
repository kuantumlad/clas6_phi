{
  gROOT->Reset();
  TStyle *st1 = new TStyle("st1","my style");
  st1->SetOptStat(0);
  st1->SetPalette(1, 0);
  st1->SetOptTitle(0);
  st1->SetDrawBorder(0);
  st1->SetAxisColor(1, "x");
  st1->SetAxisColor(1, "y");
  st1->SetFrameBorderMode(0);

  gStyle->SetCanvasColor(1);
  gStyle->SetFrameFillColor(0);
  gStyle->SetPadTopMargin(0.0);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetPadRightMargin(0.0);
  gStyle->SetOptStat(0);

  st1->SetErrorX(0);

  st1->SetLabelFont(12, "x");
  st1->SetLabelFont(12, "y");

  st1->SetLabelColor(1, "y");
  st1->SetLabelColor(1, "x");

  st1->SetLabelSize(0.06, "y");
  st1->SetLabelSize(0.06, "x");
  st1->SetNdivisions(103, "x");
  st1->SetNdivisions(505, "y");

  st1->cd();
  gROOT->ForceStyle();
}
