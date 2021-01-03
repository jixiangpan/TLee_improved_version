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

#include "./src/draw.icc"

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
  int user_val_Lee100 = 100;

  bool flag_file_Asimov = 0;

  double val_best = 0.2881;
  //TString dir_str = "./results_open5e19/";
  //TString dir_str = "./results_fake5/";
  TString dir_str = "./results_fake5_9ch/";


  
  //////////////////////////////////////////////////////////////////////////////////////// file_fc

  roostr = dir_str+"total_fc.root";
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

      double val_dchi2 = val_chi2_pull - val_chi2_gmin;
      if( val_dchi2<0 && fabs(val_dchi2)<1e-6 ) val_dchi2 = 0;
      
      map_fc_Lee_dchi2[ Lee_strength_scaled100 ].push_back( val_dchi2 );
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

    TCanvas *canv_h1_user_dchi2 = new TCanvas("canv_h1_user_dchi2", "canv_h1_user_dchi2", 900, 650);
    func_canv_margin(canv_h1_user_dchi2, 0.15, 0.1, 0.1, 0.15);
    h1_user_dchi2->Draw("hist");
    h1_user_dchi2->SetLineColor(kBlue);
    func_title_size(h1_user_dchi2, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(h1_user_dchi2, "#Delta#chi^{2} = #chi^{2}_{null} - #chi^{2}_{min}", "Entries");
    h1_user_dchi2->GetXaxis()->CenterTitle();
    h1_user_dchi2->GetYaxis()->CenterTitle();
    h1_user_dchi2->Draw("same axis");
    canv_h1_user_dchi2->SaveAs("canv_h1_user_dchi2.png");
      
    TCanvas *canv_h1_user_gmin = new TCanvas("canv_h1_user_gmin", "canv_h1_user_gmin", 900, 650);
    func_canv_margin(canv_h1_user_gmin, 0.15, 0.1, 0.1, 0.15);
    h1_user_gmin->Draw("hist");
    h1_user_gmin->SetLineColor(kBlue);
    func_title_size(h1_user_gmin, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(h1_user_gmin, "Best fit of LEE strength", "Entries");
    h1_user_gmin->GetXaxis()->CenterTitle();
    h1_user_gmin->GetYaxis()->CenterTitle();
    h1_user_gmin->Draw("same axis");
    canv_h1_user_gmin->SaveAs("canv_h1_user_gmin.png");
  }

  //////////////////////////////////////////////////////////////////////////////////////// file_Asimov
  
  roostr = dir_str+"file_Asimov.root";
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
    for(int iscan=min_Lee100; iscan<=Lee_strength_scaled100; iscan++) {
      double val_at_low = map_graph_Asimov[Lee_strength_scaled100]->Eval( iscan*1./100 );
      double val_at_hgh = map_graph_Asimov[Lee_strength_scaled100]->Eval( (iscan+1)*1./100 );
      if( val_at_low>=CI68 && val_at_hgh<=CI68 ) {
        CI68_low = iscan*1./100;
        break;
      }
    }
    for(int iscan=max_Lee100; iscan>=Lee_strength_scaled100; iscan--) {
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
    for(int iscan=min_Lee100; iscan<=Lee_strength_scaled100; iscan++) {
      double val_at_low = map_graph_Asimov[Lee_strength_scaled100]->Eval( iscan*1./100 );
      double val_at_hgh = map_graph_Asimov[Lee_strength_scaled100]->Eval( (iscan+1)*1./100 );
      if( val_at_low>=CI95 && val_at_hgh<=CI95 ) {
        CI95_low = iscan*1./100;
        break;
      }
    }
    for(int iscan=max_Lee100; iscan>=Lee_strength_scaled100; iscan--) {
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
    TCanvas *canv_toy_Asimov = new TCanvas("canv_toy_Asimov", "canv_toy_Asimov", 900, 650);
    func_canv_margin(canv_toy_Asimov, 0.15, 0.1, 0.1, 0.15);
    map_graph_Asimov[user_val_Lee100]->Draw("al");
    map_graph_Asimov[user_val_Lee100]->SetLineColor(kBlue);
    func_title_size(map_graph_Asimov[user_val_Lee100], 0.05, 0.05, 0.05, 0.05);
    func_xy_title(map_graph_Asimov[user_val_Lee100], "True LEE strength", "Confidence Level (%)");
    map_graph_Asimov[user_val_Lee100]->GetXaxis()->CenterTitle();
    map_graph_Asimov[user_val_Lee100]->GetYaxis()->CenterTitle();
    canv_toy_Asimov->SaveAs("canv_toy_Asimov.png");
  }

  //////////////////////////
  
  TCanvas *canv_gh_CI68_Asimov = new TCanvas("canv_gh_CI68_Asimov", "canv_gh_CI68_Asimov", 690, 650);
  func_canv_margin(canv_gh_CI68_Asimov, 0.15, 0.1, 0.1, 0.15);
  TH1D *h1_CI_Asimov = new TH1D("h1_CI_Asimov", "", 100, min_Lee100/100, max_Lee100/100);
  h1_CI_Asimov->Draw();
  h1_CI_Asimov->SetMinimum(min_Lee100/100);
  h1_CI_Asimov->SetMaximum(max_Lee100/100);
  func_title_size(h1_CI_Asimov, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_CI_Asimov, "True LEE strength", "LEE strength");  
  h1_CI_Asimov->GetXaxis()->CenterTitle();
  h1_CI_Asimov->GetYaxis()->CenterTitle();
  h1_CI_Asimov->GetXaxis()->SetNdivisions(509);
  h1_CI_Asimov->GetYaxis()->SetNdivisions(509);  
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
  
  canv_gh_CI68_Asimov->SaveAs("canv_gh_CI68_Asimov.png");
    
  //////////////////////////////////////////////////////////////////////////////////////// file_data

  roostr = dir_str+"file_data.root";
  TFile *file_data = new TFile(roostr, "read");
  TTree *tree_data = (TTree*)file_data->Get("tree_data");

  // Declaration of leaf types
  Double_t        Lee_bestFit_data;
  Double_t        Lee_bestFit_err;
  Double_t        chi2_gmin_data;
  vector<int>     *Lee_scan100_data;
  vector<double>  *chi2_null_scan_data;

  // List of branches
  TBranch        *b_Lee_bestFit_data;   //!
  TBranch        *b_Lee_bestFit_err;   //!
  TBranch        *b_chi2_gmin_data;   //!
  TBranch        *b_Lee_scan100_data;   //!
  TBranch        *b_chi2_null_scan_data;   //!

  // Set object pointer
  Lee_scan100_data = 0;
  chi2_null_scan_data = 0;
  
  // Set branch addresses and branch pointers
  if (!tree_data) return;
  tree_data->SetBranchAddress("Lee_bestFit_data", &Lee_bestFit_data, &b_Lee_bestFit_data);
  tree_data->SetBranchAddress("Lee_bestFit_err", &Lee_bestFit_err, &b_Lee_bestFit_err);
  tree_data->SetBranchAddress("chi2_gmin_data", &chi2_gmin_data, &b_chi2_gmin_data);
  tree_data->SetBranchAddress("Lee_scan100_data", &Lee_scan100_data, &b_Lee_scan100_data);
  tree_data->SetBranchAddress("chi2_null_scan_data", &chi2_null_scan_data, &b_chi2_null_scan_data);

  TGraph *gh_CL_data = new TGraph();
  TGraph *gh_scan_dchi2_data = new TGraph();
  
  double CI68_low = 0;
  double CI68_hgh = 1e6;    
  double CI68     = 68;
    
  double CI95_low = 0;
  double CI95_hgh = 1e6;    
  double CI95     = 95;
    
  for(int ientry=0; ientry<tree_data->GetEntries(); ientry++ ) {
    tree_data->GetEntry( ientry );

    int size_vc = Lee_scan100_data->size();
    for(int idx=0; idx<size_vc; idx++) {
      double val_chi2_null_scan_data = chi2_null_scan_data->at(idx);
      double val_dchi2_data = val_chi2_null_scan_data - chi2_gmin_data;
      if( fabs(val_dchi2_data)<1e-6 ) val_dchi2_data = 0;      
      int val_Lee100 = Lee_scan100_data->at(idx);      
      
      int line_eff = 0;
      int size_toy = map_fc_Lee_dchi2[val_Lee100].size();
      for(int itoy=0; itoy<size_toy; itoy++) {
        double val_dchi2_toy = map_fc_Lee_dchi2[val_Lee100].at(itoy);
        if( val_dchi2_toy<=val_dchi2_data ) line_eff++; 
      }// itoy

      double val_CL = line_eff*100./size_toy;
      gh_CL_data->SetPoint( gh_CL_data->GetN(), val_Lee100*1./100, val_CL );

      gh_scan_dchi2_data->SetPoint( gh_scan_dchi2_data->GetN(), val_Lee100*1./100, val_dchi2_data );
      
    }// idx

    //////////////////////////////////////////////////////////////////////////////
    
    double user_68_low = 0;
    double user_68_hgh = 0;
    for(int iscan=min_Lee100; iscan<=(int)(Lee_bestFit_data*100+0.5); iscan++) {
      double val_at_low = gh_CL_data->Eval( iscan*1./100 );
      double val_at_hgh = gh_CL_data->Eval( (iscan+1)*1./100 );
      if( val_at_low>=CI68 && val_at_hgh<=CI68 ) {
        CI68_low = iscan*1./100;
        user_68_low = iscan*1./100;
        user_68_hgh = (iscan+1)*1./100;
        break;
      }
    }
    //// case A
    if(CI68_low!=0) {      
      for(int idx=0; idx<=10000; idx++) {
        double val_xx = user_68_low + idx*0.0001;
        double val_yy = gh_CL_data->Eval( val_xx );
        //cout<<val_xx<<"\t"<<val_yy<<endl;
        if( val_yy<=CI68 ) {
          CI68_low = val_xx;
          break;
        }
      }
    }
    //// caseB
    if(CI68_low!=0) {
      while( fabs(user_68_low-user_68_hgh)>1e-4 ) {
	double user_68_mid = (user_68_low+user_68_hgh)/2;
	double val_at_mid = gh_CL_data->Eval( user_68_mid );
	if( val_at_mid>=CI68 ) user_68_low = user_68_mid;
	else user_68_hgh = user_68_mid;
      }
      //cout<<" check "<<user_68_low<<"\t"<<user_68_hgh<<endl;
      CI68_low = (user_68_low+user_68_hgh)/2;      
    }

    ////////////////////////
    ////////////////////////
    ////////////////////////
    
    for(int iscan=max_Lee100; iscan>(int)(Lee_bestFit_data*100+0.5); iscan--) {
      double val_at_low = gh_CL_data->Eval( (iscan-1)*1./100 );
      double val_at_hgh = gh_CL_data->Eval( iscan*1./100 );
      if( val_at_hgh>=CI68 && val_at_low<=CI68 ) {
        CI68_hgh = iscan*1./100;
        user_68_low = (iscan-1)*1./100;
        user_68_hgh = (iscan)*1./100;
      }
    }
    //// case A
    for(int idx=0; idx<=10000; idx++) {
      double val_xx = user_68_low + idx*0.0001;
      double val_yy = gh_CL_data->Eval( val_xx );
      //cout<<val_xx<<"\t"<<val_yy<<endl;
      if( val_yy>=CI68 ) {
        CI68_hgh = val_xx;
        break;
      }
    }
    //// caseB
    {
      while( fabs(user_68_low-user_68_hgh)>1e-4 ) {
	double user_68_mid = (user_68_low+user_68_hgh)/2;
	double val_at_mid = gh_CL_data->Eval( user_68_mid );
	if( val_at_mid>=CI68 ) user_68_hgh = user_68_mid;
	else user_68_low = user_68_mid;
      }
      //cout<<" check "<<user_68_low<<"\t"<<user_68_hgh<<endl;
      CI68_hgh = (user_68_low+user_68_hgh)/2;      
    }

    cout<<endl;
    cout<<Form(" ---> data CI68, %4.3f %4.3f", CI68_low, CI68_hgh)<<endl;
    
    //////////////////////////////////////////////////////////////////////////////
    
    double user_95_low = 0;
    double user_95_hgh = 0;
    for(int iscan=min_Lee100; iscan<=(int)(Lee_bestFit_data*100+0.5); iscan++) {
      double val_at_low = gh_CL_data->Eval( iscan*1./100 );      
      double val_at_hgh = gh_CL_data->Eval( (iscan+1)*1./100 );
        
      if( val_at_low>=CI95 && val_at_hgh<=CI95 ) {
        CI95_low = iscan*1./100;
        user_95_low = iscan*1./100;
        user_95_hgh = (iscan+1)*1./100;
        break;
      }
    }
    //// case A
    if(CI95_low!=0) { 
      for(int idx=0; idx<=10000; idx++) {
	double val_xx = user_95_low + idx*0.0001;
	double val_yy = gh_CL_data->Eval( val_xx );
	//cout<<val_xx<<"\t"<<val_yy<<endl;
	if( val_yy<=CI95 ) {
	  CI95_low = val_xx;
	  break;
	}
      }
    }
    //// caseB
    if(CI95_low!=0) {
      while( fabs(user_95_low-user_95_hgh)>1e-4 ) {
	double user_95_mid = (user_95_low+user_95_hgh)/2;
	double val_at_mid = gh_CL_data->Eval( user_95_mid );
	if( val_at_mid>=CI95 ) user_95_low = user_95_mid;
	else user_95_hgh = user_95_mid;
      }
      //cout<<" check "<<user_95_low<<"\t"<<user_95_hgh<<endl;
      CI95_low = (user_95_low+user_95_hgh)/2;      
    }
    
    ////////////////////////
    ////////////////////////
    ////////////////////////
    
    for(int iscan=max_Lee100; iscan>(int)(Lee_bestFit_data*100+0.5); iscan--) {
      double val_at_low = gh_CL_data->Eval( (iscan-1)*1./100 );
      double val_at_hgh = gh_CL_data->Eval( iscan*1./100 );
      if( val_at_hgh>=CI95 && val_at_low<=CI95 ) {
	CI95_hgh = iscan*1./100;
	user_95_low = (iscan-1)*1./100;
	user_95_hgh = (iscan)*1./100;
      }
    }
    //// case A
    for(int idx=0; idx<=10000; idx++) {
      double val_xx = user_95_low + idx*0.0001;
      double val_yy = gh_CL_data->Eval( val_xx );
      //cout<<val_xx<<"\t"<<val_yy<<endl;
      if( val_yy>=CI95 ) {
	CI95_hgh = val_xx;
	break;
      }
    }
    //// caseB
    {
      while( fabs(user_95_low-user_95_hgh)>1e-4 ) {
	double user_95_mid = (user_95_low+user_95_hgh)/2;
	double val_at_mid = gh_CL_data->Eval( user_95_mid );
	if( val_at_mid>=CI95 ) user_95_hgh = user_95_mid;
	else user_95_low = user_95_mid;
      }
      //cout<<" check "<<user_95_low<<"\t"<<user_95_hgh<<endl;
      CI95_hgh = (user_95_low+user_95_hgh)/2;      
    }
    
    cout<<Form(" ---> data CI95, %4.3f %4.3f", CI95_low, CI95_hgh)<<endl;
    cout<<endl;
    
  }// ientry
  
  //////////////////////
  
  TCanvas *canv_gh_CL_data = new TCanvas("canv_gh_CL_data", "canv_gh_CL_data", 900, 650);
  func_canv_margin(canv_gh_CL_data, 0.15, 0.1, 0.1, 0.15);
  gh_CL_data->Draw("al");
  func_title_size(gh_CL_data, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(gh_CL_data, "LEE strength", "Confidence Level (%)");  
  gh_CL_data->GetXaxis()->CenterTitle(); gh_CL_data->GetYaxis()->CenterTitle();
  gh_CL_data->GetXaxis()->SetNdivisions(509); gh_CL_data->GetYaxis()->SetNdivisions(509);  
  gh_CL_data->GetXaxis()->SetTitleOffset(1.2);
  gh_CL_data->GetYaxis()->SetRangeUser(0,110);
  TF1 *f1_CL68 = new TF1("f1_CL68", "68", 0, 1e6);
  f1_CL68->Draw("same"); f1_CL68->SetLineColor(kRed); f1_CL68->SetLineStyle(7);  
  TF1 *f1_CL95 = new TF1("f1_CL95", "95", 0, 1e6);
  f1_CL95->Draw("same"); f1_CL95->SetLineColor(kRed); f1_CL95->SetLineStyle(7);

  canv_gh_CL_data->SaveAs("canv_gh_CL_data.png");
  
  //////////////////////
  
  TCanvas *canv_gh_scan_dchi2_data = new TCanvas("canv_gh_scan_dchi2_data", "canv_gh_scan_dchi2_data", 900, 650);
  func_canv_margin(canv_gh_scan_dchi2_data, 0.15, 0.1, 0.1, 0.15);
  gh_scan_dchi2_data->Draw("al");
  func_title_size(gh_scan_dchi2_data, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(gh_scan_dchi2_data, "LEE strength", "#Delta#chi^{2} = #chi^{2}_{null} - #chi^{2}_{min}");  
  gh_scan_dchi2_data->GetXaxis()->CenterTitle(); gh_scan_dchi2_data->GetYaxis()->CenterTitle();
  gh_scan_dchi2_data->GetXaxis()->SetNdivisions(509); gh_scan_dchi2_data->GetYaxis()->SetNdivisions(509);  
  gh_scan_dchi2_data->GetXaxis()->SetTitleOffset(1.2);  

  canv_gh_scan_dchi2_data->cd(); canv_gh_scan_dchi2_data->Update();
  double x1 = gPad->GetUxmin(); double y1 = gPad->GetUymin();
  double x2 = gPad->GetUxmax(); double y2 = gPad->GetUymax();
  gh_scan_dchi2_data->GetYaxis()->SetRangeUser(0, y2);


  TGraphAsymmErrors *gh_CI95_data = new TGraphAsymmErrors();
  for(int idx=(int)(CI95_low*1000+0.5); idx<=(int)(CI95_hgh*1000+0.5); idx++) {
    double val_Lee = idx*1./1000;
    double val_dchi2 = gh_scan_dchi2_data->Eval( val_Lee );

    int ipoint = gh_CI95_data->GetN();
    gh_CI95_data->SetPoint( ipoint, val_Lee, val_dchi2 );
    gh_CI95_data->SetPointError( ipoint, 0, 0, 1e6, 0 );
  }

  TGraphAsymmErrors *gh_CI68_data = new TGraphAsymmErrors();
  for(int idx=(int)(CI68_low*1000+0.5); idx<=(int)(CI68_hgh*1000+0.5); idx++) {
    double val_Lee = idx*1./1000;
    double val_dchi2 = gh_scan_dchi2_data->Eval( val_Lee );

    int ipoint = gh_CI68_data->GetN();
    gh_CI68_data->SetPoint( ipoint, val_Lee, val_dchi2 );
    gh_CI68_data->SetPointError( ipoint, 0, 0, 1e6, 0 );
  }

  gh_CI95_data->Draw("same 3");
  gh_CI95_data->SetFillColor(kYellow);

  gh_CI68_data->Draw("same 3");
  gh_CI68_data->SetFillColor(kGreen+2);

  gh_scan_dchi2_data->Draw("same l");
  
  canv_gh_scan_dchi2_data->cd();
  canv_gh_scan_dchi2_data->Update();
  double xmin = canv_gh_scan_dchi2_data->GetUxmin();
  double ymin = canv_gh_scan_dchi2_data->GetUymin();
  double xmax = canv_gh_scan_dchi2_data->GetUxmax();
  double ymax = canv_gh_scan_dchi2_data->GetUymax();
  
  roostr = "h2_basic_scan_dchi2";
  TH2D *h2_basic_scan_dchi2 = new TH2D(roostr, "", 100, xmin, xmax, 100, ymin, ymax);
  h2_basic_scan_dchi2->Draw("same axis");
  h2_basic_scan_dchi2->GetXaxis()->SetNdivisions(509);// same as gh_scan_dchi2_data
  h2_basic_scan_dchi2->GetYaxis()->SetNdivisions(509);

  //TLegend *lg_scan_dchi2_data = new TLegend(0.65, 0.46-0.1, 0.94, 0.67-0.1);
  TLegend *lg_scan_dchi2_data = new TLegend(0.15+0.05, 0.46+0.2, 0.5+0.05, 0.67+0.2);
  lg_scan_dchi2_data->AddEntry( "", Form("Best-fit value %4.3f", val_best), "");
  lg_scan_dchi2_data->AddEntry( gh_CI68_data, "68% CI "+TString::Format("(%4.3f, %4.3f)", CI68_low, CI68_hgh), "f");
  lg_scan_dchi2_data->AddEntry( gh_CI95_data, "95% CI "+TString::Format("(%4.3f, %4.3f)", CI95_low, CI95_hgh), "f");
  lg_scan_dchi2_data->Draw();
  lg_scan_dchi2_data->SetBorderSize(0);
  lg_scan_dchi2_data->SetFillStyle(0);
  lg_scan_dchi2_data->SetTextSize(0.05);
  
  canv_gh_scan_dchi2_data->SaveAs("canv_gh_scan_dchi2_data.png");

  //map_graph_Asimov[0]->Draw("al");
}

