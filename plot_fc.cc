#include "draw.icc"

void plot_fc()
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

  roostr = "file_sum_FC_1.root";
  TFile *inputfile = new TFile(roostr, "read");
  TTree *tree = (TTree*)inputfile->Get("tree");
  
  // Declaration of leaf types
  Int_t           Lee_strength_scaled100;
  Double_t        chi2_null_data;
  Double_t        chi2_gmin_data;
  vector<double>  *chi2_null_toy;
  vector<double>  *chi2_gmin_toy;
  vector<double>  *LeeF_gmin_toy;

  // List of branches
  TBranch        *b_Lee_strength_scaled100;   //!
  TBranch        *b_chi2_null_data;   //!
  TBranch        *b_chi2_gmin_data;   //!
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
  tree->SetBranchAddress("chi2_null_data", &chi2_null_data, &b_chi2_null_data);
  tree->SetBranchAddress("chi2_gmin_data", &chi2_gmin_data, &b_chi2_gmin_data);
  tree->SetBranchAddress("chi2_null_toy", &chi2_null_toy, &b_chi2_null_toy);
  tree->SetBranchAddress("chi2_gmin_toy", &chi2_gmin_toy, &b_chi2_gmin_toy);
  tree->SetBranchAddress("LeeF_gmin_toy", &LeeF_gmin_toy, &b_LeeF_gmin_toy);

  /////////////////////////////////

  int entries = tree->GetEntries();

  map<int, vector<double> >map_Lee_vector_dchi2;
  map<int, double>map_data_dchi2;

  double scan_low = 1e6;
  double scan_hgh = 0;
  
  for(int ientry=0; ientry<entries; ientry++) {
    if( ientry%max(1, entries/10)==0 ) cout<<TString::Format(" ---> processing %6.4f, %6d", ientry*1./entries, ientry)<<endl;
    tree->GetEntry( ientry );

    map_data_dchi2[ Lee_strength_scaled100 ] = chi2_null_data - chi2_gmin_data;

    for(int idx=0; idx<(int)chi2_null_toy->size(); idx++) {
      double val_null_chi2 = chi2_null_toy->at(idx);
      double val_gmin_chi2 = chi2_gmin_toy->at(idx);
      map_Lee_vector_dchi2[Lee_strength_scaled100].push_back( val_null_chi2 - val_gmin_chi2 );      
    }// idx

    if( scan_low>Lee_strength_scaled100 ) scan_low = Lee_strength_scaled100;
    if( scan_hgh<Lee_strength_scaled100 ) scan_hgh = Lee_strength_scaled100;    
    
  }// ientry

  cout<<endl<<" ---> size of map "<<map_Lee_vector_dchi2.size()<<endl<<endl;

  map<int, double>map_val_CL;
  TGraph *gh_val_CL = new TGraph();
  TGraph *gh_val_CL_10 = new TGraph();
  
  for(auto it=map_Lee_vector_dchi2.begin(); it!=map_Lee_vector_dchi2.end(); it++) {
    
    int Lee_true_100 = it->first;
    double val_dchi2_data = map_data_dchi2[Lee_true_100];
    int line_eff = 0;    
    int size_vc = map_Lee_vector_dchi2[Lee_true_100].size();
    
    for(int idx=0; idx<size_vc; idx++) {
      double val_dchi2_toy = map_Lee_vector_dchi2[Lee_true_100].at(idx);
      if( val_dchi2_toy<val_dchi2_data ) line_eff++;	
    }// idx

    double val_CL = line_eff*100./size_vc;
    map_val_CL[Lee_true_100] = val_CL;

    gh_val_CL->SetPoint( gh_val_CL->GetN(), Lee_true_100*1./100, val_CL );

    if( Lee_true_100%10==0 && Lee_true_100>0 ) {
      gh_val_CL_10->SetPoint( gh_val_CL_10->GetN(), Lee_true_100*1./100, val_CL );
    }
    
  }// it

  //////////////////////////////////

  scan_low = scan_low/100;
  scan_hgh = scan_hgh/100;
  double step = 0.01;  
  int scan_num = (scan_hgh-scan_low)/step + 1;
  
  cout<<endl<<Form(" ---> scan_low/hgh %6.2f %6.2f", scan_low, scan_hgh)<<endl;
  
  double CL_68 = 68.;
  double CL_95 = 95.;
  double array_CL_68[2] = {0};
  double array_CL_95[2] = {0};

  for(int idx=1; idx<scan_num; idx++) {
    double scan_val_AA = scan_low + (idx-1)*step;
    double scan_val_BB = scan_low + (idx)*step;
    double CL_val_AA = gh_val_CL_10->Eval( scan_val_AA );
    double CL_val_BB = gh_val_CL_10->Eval( scan_val_BB );    
    if( (CL_val_AA>CL_68 && CL_val_BB<CL_68) ) {
      array_CL_68[0] = scan_val_AA;
      break;
    }    
  }// idx  
  for(int idx=scan_num; idx>1; idx--) {
    double scan_val_AA = scan_low + (idx-1)*step;
    double scan_val_BB = scan_low + (idx)*step;    
    double CL_val_AA = gh_val_CL_10->Eval( scan_val_AA );
    double CL_val_BB = gh_val_CL_10->Eval( scan_val_BB );    
    if( (CL_val_BB>CL_68 && CL_val_AA<CL_68) ) {
      array_CL_68[1] = scan_val_BB;
      break;
    }    
  }// idx
  
  for(int idx=1; idx<scan_num; idx++) {
    double scan_val_AA = scan_low + (idx-1)*step;
    double scan_val_BB = scan_low + (idx)*step;
    double CL_val_AA = gh_val_CL_10->Eval( scan_val_AA );
    double CL_val_BB = gh_val_CL_10->Eval( scan_val_BB );    
    if( (CL_val_AA>CL_95 && CL_val_BB<CL_95) ) {
      array_CL_95[0] = scan_val_AA;
      break;
    }    
  }// idx  
  for(int idx=scan_num; idx>1; idx--) {
    double scan_val_AA = scan_low + (idx-1)*step;
    double scan_val_BB = scan_low + (idx)*step;    
    double CL_val_AA = gh_val_CL_10->Eval( scan_val_AA );
    double CL_val_BB = gh_val_CL_10->Eval( scan_val_BB );    
    if( (CL_val_BB>CL_95 && CL_val_AA<CL_95) ) {
      array_CL_95[1] = scan_val_BB;
      break;
    }    
  }// idx
  
  cout<<TString::Format(" ---> 68CL %6.2f %6.2f", array_CL_68[0], array_CL_68[1])<<endl;
  cout<<TString::Format(" ---> 95CL %6.2f %6.2f", array_CL_95[0], array_CL_95[1])<<endl;
  
  //////////////////////////////////
  
  TCanvas *canv_val_CL = new TCanvas("canv_val_CL", "canv_val_CL", 900, 650);
  func_canv_margin(canv_val_CL, 0.15, 0.1, 0.1, 0.15);
  gh_val_CL->Draw("al");
  gh_val_CL->SetLineColor(kBlack);
  func_title_size(gh_val_CL, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(gh_val_CL, "LEE strength", "Confidence level (%)");
  gh_val_CL->GetXaxis()->CenterTitle();
  gh_val_CL->GetYaxis()->CenterTitle();
  
  gh_val_CL_10->Draw("same lp");
  gh_val_CL_10->SetLineColor(kBlue);
  gh_val_CL_10->SetMarkerColor(kBlue);


}

