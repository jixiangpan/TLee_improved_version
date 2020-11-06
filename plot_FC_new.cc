#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TLine.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "THStack.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TBranch.h"

#include "TRandom3.h"
#include "TGaxis.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TText.h"
#include "TLatex.h"

#include "TCanvas.h"
#include "TVirtualPad.h"
#include "TPad.h"
#include "TLegend.h"
#include "TString.h"
#include "TColor.h"

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

  //////////////////////////
  
  bool flag_file_fc = 0;
  int user_val_Lee100 = 200;

  bool flag_file_Asimov = 0;
  
  //////////////////////////////////////////////////////////////////////////////////////// file_fc

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

  map<int, vector<double> >map_fc_Lee_dchi2;
  map<int, vector<double> >map_fc_Lee_gmin;

  double max_Lee100 = 0;
  double min_Lee100 = 1e6;
  
  for(int ientry=0; ientry<entries_tree; ientry++) {    
    tree->GetEntry( ientry );

    int size_vc = chi2_null_toy->size();
    
    for(int idx=0; idx<size_vc; idx++) {
      double val_chi2_pull = chi2_null_toy->at( idx );
      double val_chi2_gmin = chi2_gmin_toy->at( idx );
      double val_Lee_gmin  = LeeF_gmin_toy->at( idx );
      map_fc_Lee_dchi2[ Lee_strength_scaled100 ].push_back( val_chi2_pull - val_chi2_gmin );
      map_fc_Lee_gmin[ Lee_strength_scaled100 ].push_back( val_Lee_gmin );

      if( max_Lee100<=Lee_strength_scaled100 ) max_Lee100 = Lee_strength_scaled100;
      if( min_Lee100>=Lee_strength_scaled100 ) min_Lee100 = Lee_strength_scaled100;
      
    }// idx
    
  }// ientry

  cout<<endl<<" ---> min/max Lee100: "<<min_Lee100<<", "<<max_Lee100<<endl<<endl;
  
  if( flag_file_fc ) {
    
    if( user_val_Lee100<min_Lee100 || user_val_Lee100>max_Lee100 ) exit(1);
    
    sort( map_fc_Lee_dchi2[user_val_Lee100].begin(), map_fc_Lee_dchi2[user_val_Lee100].end() );
    sort( map_fc_Lee_gmin[user_val_Lee100].begin(), map_fc_Lee_gmin[user_val_Lee100].end() );
    
    int size_vc = map_fc_Lee_dchi2[user_val_Lee100].size();
    
    roostr = "h1_user_dchi2";
    TH1D *h1_user_dchi2 = new TH1D(roostr, "",
				   100,
				   map_fc_Lee_dchi2[user_val_Lee100].at(0),
				   map_fc_Lee_dchi2[user_val_Lee100].at(size_vc-1));

    roostr = "h1_user_gmin";
    TH1D *h1_user_gmin = new TH1D(roostr, "", 100, 0, 16);
    
    for(int idx=0; idx<size_vc; idx++) {
      double val_dchi2 = map_fc_Lee_dchi2[user_val_Lee100].at(idx);
      if( fabs(val_dchi2)<1e-4 ) val_dchi2 = 0;
      h1_user_dchi2->Fill( val_dchi2 );

      double val_gmin = map_fc_Lee_gmin[user_val_Lee100].at(idx);
      h1_user_gmin->Fill( val_gmin );
    }

    TCanvas *canv_h1_user_dchi2 = new TCanvas("canv_h1_user_dchi2", "", 900, 650);
    func_canv_margin(canv_h1_user_dchi2, 0.15, 0.1, 0.1, 0.15);
    h1_user_dchi2->Draw("hist");
    h1_user_dchi2->SetLineColor(kBlue);
    func_title_size(h1_user_dchi2, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(h1_user_dchi2, "#Delta#chi^{2} = #chi^{2}_{null} - #chi^{2}_{min}", "Entries");
    h1_user_dchi2->GetXaxis()->CenterTitle();
    h1_user_dchi2->GetYaxis()->CenterTitle();
    h1_user_dchi2->Draw("same axis");
          
    TCanvas *canv_h1_user_gmin = new TCanvas("canv_h1_user_gmin", "", 900, 650);
    func_canv_margin(canv_h1_user_gmin, 0.15, 0.1, 0.1, 0.15);
    h1_user_gmin->Draw("hist");
    h1_user_gmin->SetLineColor(kBlue);
    func_title_size(h1_user_gmin, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(h1_user_gmin, "Best fit of LEE strength", "Entries");
    h1_user_gmin->GetXaxis()->CenterTitle();
    h1_user_gmin->GetYaxis()->CenterTitle();
    h1_user_gmin->Draw("same axis");            
  }

  //////////////////////////////////////////////////////////////////////////////////////// 

  roostr = "file_Asimov.root";
  TFile *file_Asimov = new TFile(roostr, "read");
  TTree *tree_Asimov = (TTree*)file_Asimov->Get("tree_Asimov");
  
  // Declaration of leaf types
  // Int_t           Lee_strength_scaled100;
  vector<double>  *Lee_scan100;
  //vector<double>  *chi2_null_toy;

  // List of branches
  // TBranch        *b_Lee_strength_scaled100;   //!
  TBranch        *b_Lee_scan100;   //!
  //TBranch        *b_chi2_null_toy;   //!

  // Set object pointer
  Lee_scan100 = 0;
  //chi2_null_toy = 0;
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  tree_Asimov->SetBranchAddress("Lee_strength_scaled100", &Lee_strength_scaled100, &b_Lee_strength_scaled100);
  tree_Asimov->SetBranchAddress("Lee_scan100", &Lee_scan100, &b_Lee_scan100);
  tree_Asimov->SetBranchAddress("chi2_null_toy", &chi2_null_toy, &b_chi2_null_toy);

  int entries_Asimov = tree_Asimov->GetEntries();

  map<int, TGraph*>map_graph_Asimov;

  TGraphAsymmErrors *gh_CI68_Asimov = new TGraphAsymmErrors();
  TGraphAsymmErrors *gh_CI95_Asimov = new TGraphAsymmErrors();
  
  
  for( int ientry=0; ientry<entries_Asimov; ientry++ ) {
    if( ientry%max(1, entries_Asimov/10)==0 ) cout<<Form(" ---> processing Asimov %5.2f, %4d", ientry*1./entries_Asimov, ientry)<<endl;
    tree_Asimov->GetEntry( ientry );
    
    map_graph_Asimov[Lee_strength_scaled100] = new TGraph();
    roostr = TString::Format("map_graph_Asimov_%06d", Lee_strength_scaled100);
    map_graph_Asimov[Lee_strength_scaled100]->SetName(roostr);

    int size_vc = Lee_scan100->size();
    for(int idx=0; idx<size_vc; idx++) {
      double val_Lee_scan100 = Lee_scan100->at(idx);
      double val_dchi2_Asimov = chi2_null_toy->at(idx);
      int line_eff = 0;
      
      int size_toy = map_fc_Lee_dchi2[val_Lee_scan100].size();
      for(int itoy=0; itoy<size_toy; itoy++) {
	double val_dchi2_toy = map_fc_Lee_dchi2[val_Lee_scan100].at(itoy);
	if( val_dchi2_toy<=val_dchi2_Asimov ) line_eff++;
      }// itoy

      double val_CL = line_eff*100./size_toy;
      map_graph_Asimov[Lee_strength_scaled100]->SetPoint( map_graph_Asimov[Lee_strength_scaled100]->GetN(), val_Lee_scan100/100, val_CL );      
    }// idx

    ///////////////////// confidence level

    double CI68_low = 0;
    double CI68_hgh = 1e6;    
    double CI68     = 68;
    for(int iscan=min_Lee100; iscan<Lee_strength_scaled100; iscan++) {
      double val_at_low = map_graph_Asimov[Lee_strength_scaled100]->Eval( iscan*1./100 );
      double val_at_hgh = map_graph_Asimov[Lee_strength_scaled100]->Eval( (iscan+1)*1./100 );
      if( val_at_low>=CI68 && val_at_hgh<=CI68 ) {
	CI68_low = iscan*1./100;
	break;
      }
    }
    for(int iscan=max_Lee100; iscan>Lee_strength_scaled100; iscan--) {
      double val_at_low = map_graph_Asimov[Lee_strength_scaled100]->Eval( (iscan-1)*1./100 );
      double val_at_hgh = map_graph_Asimov[Lee_strength_scaled100]->Eval( iscan*1./100 );
      if( val_at_hgh>=CI68 && val_at_low<=CI68 ) {
	CI68_hgh = iscan*1./100;
      }
    }    
    int npoint_CI68 = gh_CI68_Asimov->GetN();
    gh_CI68_Asimov->SetPoint( npoint_CI68, Lee_strength_scaled100*1./100, Lee_strength_scaled100*1./100 );
    gh_CI68_Asimov->SetPointError( npoint_CI68, 0, 0, Lee_strength_scaled100*1./100-CI68_low, CI68_hgh-Lee_strength_scaled100*1./100);
    
    double CI95_low = 0;
    double CI95_hgh = 1e6;    
    double CI95     = 95;
    for(int iscan=min_Lee100; iscan<Lee_strength_scaled100; iscan++) {
      double val_at_low = map_graph_Asimov[Lee_strength_scaled100]->Eval( iscan*1./100 );
      double val_at_hgh = map_graph_Asimov[Lee_strength_scaled100]->Eval( (iscan+1)*1./100 );
      if( val_at_low>=CI95 && val_at_hgh<=CI95 ) {
	CI95_low = iscan*1./100;
	break;
      }
    }
    for(int iscan=max_Lee100; iscan>Lee_strength_scaled100; iscan--) {
      double val_at_low = map_graph_Asimov[Lee_strength_scaled100]->Eval( (iscan-1)*1./100 );
      double val_at_hgh = map_graph_Asimov[Lee_strength_scaled100]->Eval( iscan*1./100 );
      if( val_at_hgh>=CI95 && val_at_low<=CI95 ) {
	CI95_hgh = iscan*1./100;
      }
    }    
    int npoint_CI95 = gh_CI95_Asimov->GetN();
    gh_CI95_Asimov->SetPoint( npoint_CI95, Lee_strength_scaled100*1./100, Lee_strength_scaled100*1./100 );
    gh_CI95_Asimov->SetPointError( npoint_CI95, 0, 0, Lee_strength_scaled100*1./100-CI95_low, CI95_hgh-Lee_strength_scaled100*1./100);
    
  }// ientry

  if( flag_file_Asimov ) {    
    TCanvas *canv_toy_Asimov = new TCanvas("canv_toy_Asimov", "", 900, 650);
    func_canv_margin(canv_toy_Asimov, 0.15, 0.1, 0.1, 0.15);
    map_graph_Asimov[user_val_Lee100]->Draw("al");
    map_graph_Asimov[user_val_Lee100]->SetLineColor(kBlue);
    func_title_size(map_graph_Asimov[user_val_Lee100], 0.05, 0.05, 0.05, 0.05);
    func_xy_title(map_graph_Asimov[user_val_Lee100], "True LEE strength", "Confidence Level (%)");
    map_graph_Asimov[user_val_Lee100]->GetXaxis()->CenterTitle();
    map_graph_Asimov[user_val_Lee100]->GetYaxis()->CenterTitle();          
  }

  //////////////////////////
  
  TCanvas *canv_gh_CI68_Asimov = new TCanvas("canv_gh_CI68_Asimov", "", 900, 650);
  func_canv_margin(canv_gh_CI68_Asimov, 0.15, 0.2, 0.1, 0.15);
  TH1D *h1_CI_Asimov = new TH1D("h1_CI_Asimov", "", 100, min_Lee100/100, max_Lee100/100);
  h1_CI_Asimov->Draw();
  h1_CI_Asimov->SetMinimum(min_Lee100/100);
  h1_CI_Asimov->SetMaximum(max_Lee100/100);
  func_title_size(h1_CI_Asimov, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_CI_Asimov, "True LEE strength", "LEE strength");  
  h1_CI_Asimov->GetXaxis()->CenterTitle();
  h1_CI_Asimov->GetYaxis()->CenterTitle();
  h1_CI_Asimov->GetXaxis()->SetTitleOffset(1.2);
  
  gh_CI95_Asimov->Draw("same 3");
  gh_CI95_Asimov->SetFillColor(kYellow);

  gh_CI68_Asimov->Draw("same 3");
  gh_CI68_Asimov->SetFillColor(kGreen+2);

  TF1 *f1_y2x = new TF1("f1_y2x", "x", 0, 1e6);
  f1_y2x->Draw("same");
  f1_y2x->SetLineColor(kBlack);
  f1_y2x->SetLineStyle(7);
  
  h1_CI_Asimov->Draw("same axis");
  
}

