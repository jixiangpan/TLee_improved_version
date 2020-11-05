#include "draw.icc"

void plot_FC_new()
{  
  //////////////////////////////////////////////////////////////////////////////////////// Draw style
  
  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kBird);

  double snWidth = 2;

  // use medium bold lines and thick markers
  gStyle->SetLineWidth(snWidth);
  gStyle->SetFrameLineWidth(snWidth);
  gStyle->SetHistLineWidth(snWidth);
  gStyle->SetFuncWidth(snWidth);
  gStyle->SetGridWidth(snWidth);
  gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.0);
  gStyle->SetEndErrorSize(4);
  gStyle->SetEndErrorSize(0);

  ////////////////////////////////////////////////////////////////////////////////////////
  
  TString roostr = "";
  
  roostr = "file_FC_toy_Lee0to5.root";
  TFile *file_fc = new TFile(roostr, "read");
  TTree *tree = (TTree*)file_fc->Get("tree");
  
  // Declaration of leaf types
  Int_t           Lee_strength_scaled100;
  vector<double>  *chi2_null_toy;
  vector<double>  *chi2_gmin_toy;
  vector<double>  *LeeF_gmin_toy;

  // List of branches
  TBranch        *b_Lee_strength_scaled100;   //!
  TBranch        *b_chi2_null_toy;   //!
  TBranch        *b_chi2_gmin_toy;   //!
  TBranch        *b_LeeF_gmin_toy;   //!

  // Set object pointer
  chi2_null_toy = 0;
  chi2_gmin_toy = 0;
  LeeF_gmin_toy = 0;
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree->SetBranchAddress("Lee_strength_scaled100", &Lee_strength_scaled100, &b_Lee_strength_scaled100);
  tree->SetBranchAddress("chi2_null_toy", &chi2_null_toy, &b_chi2_null_toy);
  tree->SetBranchAddress("chi2_gmin_toy", &chi2_gmin_toy, &b_chi2_gmin_toy);
  tree->SetBranchAddress("LeeF_gmin_toy", &LeeF_gmin_toy, &b_LeeF_gmin_toy);

  int entries_tree = tree->GetEntries();
  cout<<endl<<" ---> entries_tree "<<entries_tree<<endl<<endl;
  
}

