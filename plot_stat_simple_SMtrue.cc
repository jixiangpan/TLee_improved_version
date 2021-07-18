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

void plot_stat_simple_SMtrue(double exLEE)
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

  double data_dchi2_LEE2GMIN = exLEE;
  double data_dchi2_LEE2SM   = data_dchi2_LEE2GMIN;
  
  TString roostr = "";
  roostr = "./result_boxopen/file_exLEE_simple.root";
 
  ////////////////////////////////////////////////////////////////////////////////////////

  TFile *file_out = new TFile(roostr, "read");
  TTree *tree = (TTree*)file_out->Get("tree");

  // Declaration of leaf types
  Double_t        chi2_null_null8sm_true8sm;
  Double_t        chi2_gmin_null8sm_true8sm;
  Double_t        chi2_null_null8Lee_true8Lee;
  Double_t        chi2_gmin_null8Lee_true8Lee;

  // List of branches
  TBranch        *b_chi2_null_null8sm_true8sm;   //!
  TBranch        *b_chi2_gmin_null8sm_true8sm;   //!
  TBranch        *b_chi2_null_null8Lee_true8Lee;   //!
  TBranch        *b_chi2_gmin_null8Lee_true8Lee;   //!
  
  // Set branch addresses and branch pointers  
  tree->SetBranchAddress("chi2_null_null8sm_true8sm", &chi2_null_null8sm_true8sm, &b_chi2_null_null8sm_true8sm);
  tree->SetBranchAddress("chi2_gmin_null8sm_true8sm", &chi2_gmin_null8sm_true8sm, &b_chi2_gmin_null8sm_true8sm);
  tree->SetBranchAddress("chi2_null_null8Lee_true8Lee", &chi2_null_null8Lee_true8Lee, &b_chi2_null_null8Lee_true8Lee);
  tree->SetBranchAddress("chi2_gmin_null8Lee_true8Lee", &chi2_gmin_null8Lee_true8Lee, &b_chi2_gmin_null8Lee_true8Lee);

  int entries = tree->GetEntries();  
  //cout<<endl<<" ---> entries "<<entries<<endl<<endl;

  vector<double>vc_dchi2_LEE2GMIN;
  vector<double>vc_dchi2_LEE2SM;
  
  cout<<endl;
  for(int ientry=0; ientry<entries; ientry++) {
    if( ientry%max(entries/10,1)==0 ) cout<<TString::Format(" ---> processing : %4.2f, %8d", ientry*1./entries, ientry)<<endl;
    tree->GetEntry( ientry );

    double dchi2_LEE2SM = chi2_gmin_null8sm_true8sm - chi2_null_null8sm_true8sm;
    //if( dchi2_LEE2SM<0 && fabs(dchi2_LEE2SM)<1e-6 ) dchi2_LEE2SM = 0;
    vc_dchi2_LEE2SM.push_back( dchi2_LEE2SM );
  }
  cout<<endl;

  ///////////////////////////////////////////////////////////////////// LEE2SM
  
  sort( vc_dchi2_LEE2SM.begin(), vc_dchi2_LEE2SM.end() );
  
  int size_vc_dchi2_LEE2SM = vc_dchi2_LEE2SM.size();
  int line_LEE2SM = 0;

  TH1D *h1_pdf_dchi2_LEE2SM = new TH1D("dchi2_LEE2SM", "", 100, -40, 80);

  
  for(int idx=0; idx<size_vc_dchi2_LEE2SM; idx++) {
    double dchi2_LEE2SM = vc_dchi2_LEE2SM.at(idx);
    h1_pdf_dchi2_LEE2SM->Fill( dchi2_LEE2SM );
  }
  h1_pdf_dchi2_LEE2SM->Scale( 1./size_vc_dchi2_LEE2SM );

  double xval_peak = h1_pdf_dchi2_LEE2SM->GetBinCenter( h1_pdf_dchi2_LEE2SM->GetMaximumBin() );
  
  if( data_dchi2_LEE2SM > xval_peak ) {
    for(int idx=0; idx<size_vc_dchi2_LEE2SM; idx++) {
      double dchi2_LEE2SM = vc_dchi2_LEE2SM.at(idx);
      if( dchi2_LEE2SM>=data_dchi2_LEE2SM ) line_LEE2SM++;
    }
  }
  else {
    for(int idx=0; idx<size_vc_dchi2_LEE2SM; idx++) {
      double dchi2_LEE2SM = vc_dchi2_LEE2SM.at(idx);
      if( dchi2_LEE2SM<=data_dchi2_LEE2SM ) line_LEE2SM++;
    }
  }

  
  double pvalue_LEE2SM = line_LEE2SM*1./size_vc_dchi2_LEE2SM;
  double significance_LEE2SM = sqrt( TMath::ChisquareQuantile( 1-pvalue_LEE2SM*2, 1 ) );// one side
  
  cout<<endl<<TString::Format(" ---> simple vs. simple, p-value %8.6f, significance %4.2f sigma",
			      pvalue_LEE2SM, significance_LEE2SM)<<endl<<endl;

  roostr = "canv_h1_pdf_dchi2_LEE2SM";
  TCanvas *canv_h1_pdf_dchi2_LEE2SM = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_pdf_dchi2_LEE2SM, 0.15, 0.1, 0.1, 0.15);
  h1_pdf_dchi2_LEE2SM->Draw("hist f");  
  h1_pdf_dchi2_LEE2SM->SetLineColor(kRed);  
  h1_pdf_dchi2_LEE2SM->SetFillColor(kRed-10); h1_pdf_dchi2_LEE2SM->SetFillStyle(1001);
  func_title_size(h1_pdf_dchi2_LEE2SM, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_pdf_dchi2_LEE2SM, "#Delta#chi^{2} = #chi^{2}_{LEEx=1} - #chi^{2}_{SM}", "PDF");
  h1_pdf_dchi2_LEE2SM->GetXaxis()->CenterTitle(); h1_pdf_dchi2_LEE2SM->GetYaxis()->CenterTitle();

  canv_h1_pdf_dchi2_LEE2SM->cd(); canv_h1_pdf_dchi2_LEE2SM->Update();
  double x1_LEE2SM = gPad->GetUxmin(); double y1_LEE2SM = gPad->GetUymin();
  double x2_LEE2SM = gPad->GetUxmax(); double y2_LEE2SM = gPad->GetUymax();

  TLine *rooline_LEE2SM = new TLine( data_dchi2_LEE2SM, 0, data_dchi2_LEE2SM, y2_LEE2SM );
  rooline_LEE2SM->Draw("same"); rooline_LEE2SM->SetLineWidth(2);
  rooline_LEE2SM->SetLineColor(kBlack); rooline_LEE2SM->SetLineStyle(7);
  
  TPaveText *pt_LEE2SM = new TPaveText( x1_LEE2SM + (x2_LEE2SM-x1_LEE2SM)*0.6, y2_LEE2SM*0.56,
					x1_LEE2SM + (x2_LEE2SM-x1_LEE2SM)*0.8, y2_LEE2SM*0.88, "l");  
  pt_LEE2SM->SetTextSize(0.05); pt_LEE2SM->SetTextFont(42); pt_LEE2SM->SetTextAlign(11);
  pt_LEE2SM->SetBorderSize(0); pt_LEE2SM->SetFillStyle(0);
  pt_LEE2SM->AddText( Form("#Delta#chi^{2} = %5.3f", data_dchi2_LEE2SM) );
  //pt_LEE2SM->AddText( Form("#rightarrow p-value = %5.3f", pvalue_LEE2SM) );
  if( pvalue_LEE2SM>1e-3 ) pt_LEE2SM->AddText( Form("#rightarrow p-value = %5.3f", pvalue_LEE2SM) );
  else {
    double val_log10 = TMath::Log10(pvalue_LEE2SM);
    pt_LEE2SM->AddText( Form("#rightarrow p-value = 10^{%3.1f}", val_log10) );
  }
  pt_LEE2SM->AddText( Form("#rightarrow %3.2f#sigma", significance_LEE2SM) );
  pt_LEE2SM->AddText("Simple vs. simple");
  pt_LEE2SM->Draw();

  roostr = "canv_simple2simple_pdf.png";
  canv_h1_pdf_dchi2_LEE2SM->SaveAs(roostr);

  //////////////

  TH1D *h1_cdf_dchi2_LEE2SM = (TH1D*)h1_pdf_dchi2_LEE2SM->Clone("h1_cdf_dchi2_LEE2SM");
  for(int ibin=1; ibin<=h1_cdf_dchi2_LEE2SM->GetNbinsX(); ibin++) {
    double cdf = 0;
    for(int jbin=1; jbin<=ibin; jbin++) cdf += h1_pdf_dchi2_LEE2SM->GetBinContent(jbin);
    h1_cdf_dchi2_LEE2SM->SetBinContent( ibin, cdf );
  }
  
  roostr = "canv_h1_cdf_dchi2_LEE2SM";
  TCanvas *canv_h1_cdf_dchi2_LEE2SM = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_cdf_dchi2_LEE2SM, 0.15, 0.1, 0.1, 0.15); 
  h1_cdf_dchi2_LEE2SM->Draw("hist"); h1_cdf_dchi2_LEE2SM->SetFillColor(0);
  h1_cdf_dchi2_LEE2SM->SetMinimum(0); h1_cdf_dchi2_LEE2SM->SetMaximum(1.1);
  h1_cdf_dchi2_LEE2SM->SetYTitle("CDF");

  canv_h1_cdf_dchi2_LEE2SM->cd(); canv_h1_cdf_dchi2_LEE2SM->Update();
  double x1_cdf_LEE2SM = gPad->GetUxmin(); double y1_cdf_LEE2SM = gPad->GetUymin();
  double x2_cdf_LEE2SM = gPad->GetUxmax(); double y2_cdf_LEE2SM = gPad->GetUymax();

  TLine *rooline_cdf_LEE2SM = new TLine( data_dchi2_LEE2SM, 0, data_dchi2_LEE2SM, y2_cdf_LEE2SM );
  rooline_cdf_LEE2SM->Draw("same"); rooline_cdf_LEE2SM->SetLineWidth(2);
  rooline_cdf_LEE2SM->SetLineColor(kBlack); rooline_cdf_LEE2SM->SetLineStyle(7);
  
  TPaveText *pt_cdf_LEE2SM = new TPaveText( x1_cdf_LEE2SM + (x2_cdf_LEE2SM-x1_cdf_LEE2SM)*0.6, y2_cdf_LEE2SM*0.48,
					    x1_cdf_LEE2SM + (x2_cdf_LEE2SM-x1_cdf_LEE2SM)*0.8, y2_cdf_LEE2SM*0.8, "l");  
  pt_cdf_LEE2SM->SetTextSize(0.05); pt_cdf_LEE2SM->SetTextFont(42); pt_cdf_LEE2SM->SetTextAlign(11);
  pt_cdf_LEE2SM->SetBorderSize(0); pt_cdf_LEE2SM->SetFillStyle(0);
  pt_cdf_LEE2SM->AddText( Form("#Delta#chi^{2} = %5.3f", data_dchi2_LEE2SM) );
  //pt_cdf_LEE2SM->AddText( Form("#rightarrow p-value = %5.3f", pvalue_LEE2SM) );
  if( pvalue_LEE2SM>1e-3 ) pt_cdf_LEE2SM->AddText( Form("#rightarrow p-value = %5.3f", pvalue_LEE2SM) );
  else {
    double val_log10 = TMath::Log10(pvalue_LEE2SM);
    pt_cdf_LEE2SM->AddText( Form("#rightarrow p-value = 10^{%3.1f}", val_log10) );
  }
  pt_cdf_LEE2SM->AddText( Form("#rightarrow %3.2f#sigma", significance_LEE2SM) );
  pt_cdf_LEE2SM->AddText("Simple vs. simple");
  pt_cdf_LEE2SM->Draw();

  roostr = "canv_simple2simple_cdf.png";
  canv_h1_cdf_dchi2_LEE2SM->SaveAs(roostr);

  
}

