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

void plot_stat_exLEE_exSM(double exLEE, double exSM=0)
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
  double data_dchi2_SM2GMIN  = exSM;
  
  TString roostr = "";
  roostr = "./result_boxopen/file_exLEE_exSM.root";

  bool FLAG_SM2GMIN_is0 = 1;
  if( exSM!=0 ) FLAG_SM2GMIN_is0 = 0;

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
  vector<double>vc_dchi2_SM2GMIN;
  
  cout<<endl;
  for(int ientry=0; ientry<entries; ientry++) {
    if( ientry%max(entries/10,1)==0 ) cout<<TString::Format(" ---> processing : %4.2f, %8d", ientry*1./entries, ientry)<<endl;
    tree->GetEntry( ientry );

    double dchi2_LEE2GMIN = chi2_null_null8Lee_true8Lee - chi2_gmin_null8Lee_true8Lee;
    if( dchi2_LEE2GMIN<0 && fabs(dchi2_LEE2GMIN)<1e-6 ) dchi2_LEE2GMIN = 1e-6;
    if( dchi2_LEE2GMIN>=0 )
      vc_dchi2_LEE2GMIN.push_back( dchi2_LEE2GMIN );

    if( !FLAG_SM2GMIN_is0 ) {
      double dchi2_SM2GMIN = chi2_null_null8sm_true8sm - chi2_gmin_null8sm_true8sm;
      if( dchi2_SM2GMIN<0 && fabs(dchi2_SM2GMIN)<1e-6 ) dchi2_SM2GMIN = 1e-6;
      if( dchi2_SM2GMIN>=0 )
	vc_dchi2_SM2GMIN.push_back( dchi2_SM2GMIN );
    }
  }
  cout<<endl;

  ///////////////////////////////////////////////////////////////////// LEE2GMIN
  
  sort( vc_dchi2_LEE2GMIN.begin(), vc_dchi2_LEE2GMIN.end() );
  
  int size_vc_dchi2_LEE2GMIN = vc_dchi2_LEE2GMIN.size();
  int line_LEE2GMIN = 0;

  TH1D *h1_pdf_dchi2_LEE2GMIN = new TH1D("dchi2_LEE2GMIN", "", 100, 0, 25);
  
  for(int idx=0; idx<size_vc_dchi2_LEE2GMIN; idx++) {
    double dchi2_LEE2GMIN = vc_dchi2_LEE2GMIN.at(idx);
    if( dchi2_LEE2GMIN>=data_dchi2_LEE2GMIN ) line_LEE2GMIN++;

    h1_pdf_dchi2_LEE2GMIN->Fill( dchi2_LEE2GMIN );
  }
  h1_pdf_dchi2_LEE2GMIN->Scale( 1./size_vc_dchi2_LEE2GMIN );
    
  double pvalue_LEE2GMIN = line_LEE2GMIN*1./size_vc_dchi2_LEE2GMIN;
  double significance_LEE2GMIN = sqrt( TMath::ChisquareQuantile( 1-pvalue_LEE2GMIN, 1 ) );// two sides
  
  cout<<endl<<TString::Format(" ---> rejecting LEEx=1, p-value %8.6f, significance %4.2f sigma",
			      pvalue_LEE2GMIN, significance_LEE2GMIN)<<endl<<endl;

  roostr = "canv_h1_pdf_dchi2_LEE2GMIN";
  TCanvas *canv_h1_pdf_dchi2_LEE2GMIN = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_pdf_dchi2_LEE2GMIN, 0.15, 0.1, 0.1, 0.15);
  h1_pdf_dchi2_LEE2GMIN->Draw("hist f");  
  h1_pdf_dchi2_LEE2GMIN->SetLineColor(kBlue);  
  h1_pdf_dchi2_LEE2GMIN->SetFillColor(kBlue-10); h1_pdf_dchi2_LEE2GMIN->SetFillStyle(1001);
  func_title_size(h1_pdf_dchi2_LEE2GMIN, 0.05, 0.05, 0.05, 0.05);
  func_xy_title(h1_pdf_dchi2_LEE2GMIN, "#Delta#chi^{2} = #chi^{2}_{LEEx=1} - #chi^{2}_{min}", "PDF");
  h1_pdf_dchi2_LEE2GMIN->GetXaxis()->CenterTitle(); h1_pdf_dchi2_LEE2GMIN->GetYaxis()->CenterTitle();

  canv_h1_pdf_dchi2_LEE2GMIN->cd(); canv_h1_pdf_dchi2_LEE2GMIN->Update();
  double x1_LEE2GMIN = gPad->GetUxmin(); double y1_LEE2GMIN = gPad->GetUymin();
  double x2_LEE2GMIN = gPad->GetUxmax(); double y2_LEE2GMIN = gPad->GetUymax();

  TLine *rooline_LEE2GMIN = new TLine( data_dchi2_LEE2GMIN, 0, data_dchi2_LEE2GMIN, y2_LEE2GMIN );
  rooline_LEE2GMIN->Draw("same"); rooline_LEE2GMIN->SetLineWidth(2);
  rooline_LEE2GMIN->SetLineColor(kBlack); rooline_LEE2GMIN->SetLineStyle(7);
  
  TPaveText *pt_LEE2GMIN = new TPaveText( x1_LEE2GMIN + (x2_LEE2GMIN-x1_LEE2GMIN)*0.6, y2_LEE2GMIN*0.56,
					  x1_LEE2GMIN + (x2_LEE2GMIN-x1_LEE2GMIN)*0.8, y2_LEE2GMIN*0.88, "l");  
  pt_LEE2GMIN->SetTextSize(0.05); pt_LEE2GMIN->SetTextFont(42); pt_LEE2GMIN->SetTextAlign(11);
  pt_LEE2GMIN->SetBorderSize(0); pt_LEE2GMIN->SetFillStyle(0);
  pt_LEE2GMIN->AddText( Form("#Delta#chi^{2} = %5.3f", data_dchi2_LEE2GMIN) );
  if( pvalue_LEE2GMIN>1e-3 ) pt_LEE2GMIN->AddText( Form("#rightarrow p-value = %5.3f", pvalue_LEE2GMIN) );
  else {
    double val_log10 = TMath::Log10(pvalue_LEE2GMIN);
    pt_LEE2GMIN->AddText( Form("#rightarrow p-value = 10^{%3.1f}", val_log10) );
  }
  pt_LEE2GMIN->AddText( Form("#rightarrow %3.2f#sigma", significance_LEE2GMIN) );
  pt_LEE2GMIN->AddText("Reject LEEx=1");
  pt_LEE2GMIN->Draw();

  roostr = "canv_exLEE_pdf.png";
  canv_h1_pdf_dchi2_LEE2GMIN->SaveAs(roostr);

  canv_h1_pdf_dchi2_LEE2GMIN->SetLogy();
  roostr = "canv_exLEE_pdf_logy.png";
  canv_h1_pdf_dchi2_LEE2GMIN->SaveAs(roostr);

  //////////////

  TH1D *h1_cdf_dchi2_LEE2GMIN = (TH1D*)h1_pdf_dchi2_LEE2GMIN->Clone("h1_cdf_dchi2_LEE2GMIN");
  for(int ibin=1; ibin<=h1_cdf_dchi2_LEE2GMIN->GetNbinsX(); ibin++) {
    double cdf = 0;
    for(int jbin=1; jbin<=ibin; jbin++) cdf += h1_pdf_dchi2_LEE2GMIN->GetBinContent(jbin);
    h1_cdf_dchi2_LEE2GMIN->SetBinContent( ibin, cdf );
  }
  
  roostr = "canv_h1_cdf_dchi2_LEE2GMIN";
  TCanvas *canv_h1_cdf_dchi2_LEE2GMIN = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_cdf_dchi2_LEE2GMIN, 0.15, 0.1, 0.1, 0.15); 
  h1_cdf_dchi2_LEE2GMIN->Draw("hist"); h1_cdf_dchi2_LEE2GMIN->SetFillColor(0);
  h1_cdf_dchi2_LEE2GMIN->SetMinimum(0); h1_cdf_dchi2_LEE2GMIN->SetMaximum(1.1);
  h1_cdf_dchi2_LEE2GMIN->SetYTitle("CDF");

  canv_h1_cdf_dchi2_LEE2GMIN->cd(); canv_h1_cdf_dchi2_LEE2GMIN->Update();
  double x1_cdf_LEE2GMIN = gPad->GetUxmin(); double y1_cdf_LEE2GMIN = gPad->GetUymin();
  double x2_cdf_LEE2GMIN = gPad->GetUxmax(); double y2_cdf_LEE2GMIN = gPad->GetUymax();

  TLine *rooline_cdf_LEE2GMIN = new TLine( data_dchi2_LEE2GMIN, 0, data_dchi2_LEE2GMIN, y2_cdf_LEE2GMIN );
  rooline_cdf_LEE2GMIN->Draw("same"); rooline_cdf_LEE2GMIN->SetLineWidth(2);
  rooline_cdf_LEE2GMIN->SetLineColor(kBlack); rooline_cdf_LEE2GMIN->SetLineStyle(7);
  
  TPaveText *pt_cdf_LEE2GMIN = new TPaveText( x1_cdf_LEE2GMIN + (x2_cdf_LEE2GMIN-x1_cdf_LEE2GMIN)*0.6, y2_cdf_LEE2GMIN*0.48,
					      x1_cdf_LEE2GMIN + (x2_cdf_LEE2GMIN-x1_cdf_LEE2GMIN)*0.8, y2_cdf_LEE2GMIN*0.8, "l");  
  pt_cdf_LEE2GMIN->SetTextSize(0.05); pt_cdf_LEE2GMIN->SetTextFont(42); pt_cdf_LEE2GMIN->SetTextAlign(11);
  pt_cdf_LEE2GMIN->SetBorderSize(0); pt_cdf_LEE2GMIN->SetFillStyle(0);
  pt_cdf_LEE2GMIN->AddText( Form("#Delta#chi^{2} = %5.3f", data_dchi2_LEE2GMIN) );
  if( pvalue_LEE2GMIN>1e-3 ) pt_cdf_LEE2GMIN->AddText( Form("#rightarrow p-value = %5.3f", pvalue_LEE2GMIN) );
  else {
    double val_log10 = TMath::Log10(pvalue_LEE2GMIN);
    pt_cdf_LEE2GMIN->AddText( Form("#rightarrow p-value = 10^{%3.1f}", val_log10) );
  }
  pt_cdf_LEE2GMIN->AddText( Form("#rightarrow %3.2f#sigma", significance_LEE2GMIN) );
  pt_cdf_LEE2GMIN->AddText("Reject LEEx=1");
  pt_cdf_LEE2GMIN->Draw();

  roostr = "canv_exLEE_cdf.png";
  canv_h1_cdf_dchi2_LEE2GMIN->SaveAs(roostr);

  ///////////////////////////////////////////////////////////////////// SM2GMIN

  if( !FLAG_SM2GMIN_is0 ) {
    sort( vc_dchi2_SM2GMIN.begin(), vc_dchi2_SM2GMIN.end() );
  
    int size_vc_dchi2_SM2GMIN = vc_dchi2_SM2GMIN.size();
    int line_SM2GMIN = 0;

    TH1D *h1_pdf_dchi2_SM2GMIN = new TH1D("dchi2_SM2GMIN", "", 100, 0, 25);
  
    for(int idx=0; idx<size_vc_dchi2_SM2GMIN; idx++) {
      double dchi2_SM2GMIN = vc_dchi2_SM2GMIN.at(idx);
      if( dchi2_SM2GMIN>=data_dchi2_SM2GMIN ) line_SM2GMIN++;

      h1_pdf_dchi2_SM2GMIN->Fill( dchi2_SM2GMIN );
    }
    h1_pdf_dchi2_SM2GMIN->Scale( 1./size_vc_dchi2_SM2GMIN );
    
    double pvalue_SM2GMIN = line_SM2GMIN*1./size_vc_dchi2_SM2GMIN;
    double significance_SM2GMIN = sqrt( TMath::ChisquareQuantile( 1-pvalue_SM2GMIN, 1 ) );// two sides
  
    cout<<endl<<TString::Format(" ---> simple vs. simple, p-value %8.6f, significance %4.2f sigma",
				pvalue_SM2GMIN, significance_SM2GMIN)<<endl<<endl;

    roostr = "canv_h1_pdf_dchi2_SM2GMIN";
    TCanvas *canv_h1_pdf_dchi2_SM2GMIN = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_pdf_dchi2_SM2GMIN, 0.15, 0.1, 0.1, 0.15);
    h1_pdf_dchi2_SM2GMIN->Draw("hist f");  
    h1_pdf_dchi2_SM2GMIN->SetLineColor(kRed);  
    h1_pdf_dchi2_SM2GMIN->SetFillColor(kRed-10); h1_pdf_dchi2_SM2GMIN->SetFillStyle(1001);
    func_title_size(h1_pdf_dchi2_SM2GMIN, 0.05, 0.05, 0.05, 0.05);
    func_xy_title(h1_pdf_dchi2_SM2GMIN, "#Delta#chi^{2} = #chi^{2}_{SM} - #chi^{2}_{min}", "PDF");
    h1_pdf_dchi2_SM2GMIN->GetXaxis()->CenterTitle(); h1_pdf_dchi2_SM2GMIN->GetYaxis()->CenterTitle();

    canv_h1_pdf_dchi2_SM2GMIN->cd(); canv_h1_pdf_dchi2_SM2GMIN->Update();
    double x1_SM2GMIN = gPad->GetUxmin(); double y1_SM2GMIN = gPad->GetUymin();
    double x2_SM2GMIN = gPad->GetUxmax(); double y2_SM2GMIN = gPad->GetUymax();

    TLine *rooline_SM2GMIN = new TLine( data_dchi2_SM2GMIN, 0, data_dchi2_SM2GMIN, y2_SM2GMIN );
    rooline_SM2GMIN->Draw("same"); rooline_SM2GMIN->SetLineWidth(2);
    rooline_SM2GMIN->SetLineColor(kBlack); rooline_SM2GMIN->SetLineStyle(7);
  
    TPaveText *pt_SM2GMIN = new TPaveText( x1_SM2GMIN + (x2_SM2GMIN-x1_SM2GMIN)*0.6, y2_SM2GMIN*0.56,
					   x1_SM2GMIN + (x2_SM2GMIN-x1_SM2GMIN)*0.8, y2_SM2GMIN*0.88, "l");  
    pt_SM2GMIN->SetTextSize(0.05); pt_SM2GMIN->SetTextFont(42); pt_SM2GMIN->SetTextAlign(11);
    pt_SM2GMIN->SetBorderSize(0); pt_SM2GMIN->SetFillStyle(0);
    pt_SM2GMIN->AddText( Form("#Delta#chi^{2} = %5.3f", data_dchi2_SM2GMIN) );
    //pt_SM2GMIN->AddText( Form("#rightarrow p-value = %5.3f", pvalue_SM2GMIN) );
    if( pvalue_SM2GMIN>1e-3 ) pt_SM2GMIN->AddText( Form("#rightarrow p-value = %5.3f", pvalue_SM2GMIN) );
    else {
      double val_log10 = TMath::Log10(pvalue_SM2GMIN);
      pt_SM2GMIN->AddText( Form("#rightarrow p-value = 10^{%3.1f}", val_log10) );
    }
    pt_SM2GMIN->AddText( Form("#rightarrow %3.2f#sigma", significance_SM2GMIN) );
    pt_SM2GMIN->AddText("Disfavor SM");
    pt_SM2GMIN->Draw();

    roostr = "canv_exSM_pdf.png";
    canv_h1_pdf_dchi2_SM2GMIN->SaveAs(roostr);

    canv_h1_pdf_dchi2_SM2GMIN->SetLogy();
    roostr = "canv_exSM_pdf_logy.png";
    canv_h1_pdf_dchi2_SM2GMIN->SaveAs(roostr);

    //////////////

    TH1D *h1_cdf_dchi2_SM2GMIN = (TH1D*)h1_pdf_dchi2_SM2GMIN->Clone("h1_cdf_dchi2_SM2GMIN");
    for(int ibin=1; ibin<=h1_cdf_dchi2_SM2GMIN->GetNbinsX(); ibin++) {
      double cdf = 0;
      for(int jbin=1; jbin<=ibin; jbin++) cdf += h1_pdf_dchi2_SM2GMIN->GetBinContent(jbin);
      h1_cdf_dchi2_SM2GMIN->SetBinContent( ibin, cdf );
    }
  
    roostr = "canv_h1_cdf_dchi2_SM2GMIN";
    TCanvas *canv_h1_cdf_dchi2_SM2GMIN = new TCanvas(roostr, roostr, 900, 650);
    func_canv_margin(canv_h1_cdf_dchi2_SM2GMIN, 0.15, 0.1, 0.1, 0.15); 
    h1_cdf_dchi2_SM2GMIN->Draw("hist"); h1_cdf_dchi2_SM2GMIN->SetFillColor(0);
    h1_cdf_dchi2_SM2GMIN->SetMinimum(0); h1_cdf_dchi2_SM2GMIN->SetMaximum(1.1);
    h1_cdf_dchi2_SM2GMIN->SetYTitle("CDF");

    canv_h1_cdf_dchi2_SM2GMIN->cd(); canv_h1_cdf_dchi2_SM2GMIN->Update();
    double x1_cdf_SM2GMIN = gPad->GetUxmin(); double y1_cdf_SM2GMIN = gPad->GetUymin();
    double x2_cdf_SM2GMIN = gPad->GetUxmax(); double y2_cdf_SM2GMIN = gPad->GetUymax();

    TLine *rooline_cdf_SM2GMIN = new TLine( data_dchi2_SM2GMIN, 0, data_dchi2_SM2GMIN, y2_cdf_SM2GMIN );
    rooline_cdf_SM2GMIN->Draw("same"); rooline_cdf_SM2GMIN->SetLineWidth(2);
    rooline_cdf_SM2GMIN->SetLineColor(kBlack); rooline_cdf_SM2GMIN->SetLineStyle(7);
  
    TPaveText *pt_cdf_SM2GMIN = new TPaveText( x1_cdf_SM2GMIN + (x2_cdf_SM2GMIN-x1_cdf_SM2GMIN)*0.6, y2_cdf_SM2GMIN*0.48,
					       x1_cdf_SM2GMIN + (x2_cdf_SM2GMIN-x1_cdf_SM2GMIN)*0.8, y2_cdf_SM2GMIN*0.8, "l");  
    pt_cdf_SM2GMIN->SetTextSize(0.05); pt_cdf_SM2GMIN->SetTextFont(42); pt_cdf_SM2GMIN->SetTextAlign(11);
    pt_cdf_SM2GMIN->SetBorderSize(0); pt_cdf_SM2GMIN->SetFillStyle(0);
    pt_cdf_SM2GMIN->AddText( Form("#Delta#chi^{2} = %5.3f", data_dchi2_SM2GMIN) );
    //pt_cdf_SM2GMIN->AddText( Form("#rightarrow p-value = %5.3f", pvalue_SM2GMIN) );
    if( pvalue_SM2GMIN>1e-3 ) pt_cdf_SM2GMIN->AddText( Form("#rightarrow p-value = %5.3f", pvalue_SM2GMIN) );
    else {
      double val_log10 = TMath::Log10(pvalue_SM2GMIN);
      pt_cdf_SM2GMIN->AddText( Form("#rightarrow p-value = 10^{%3.1f}", val_log10) );
    }
    pt_cdf_SM2GMIN->AddText( Form("#rightarrow %3.2f#sigma", significance_SM2GMIN) );
    pt_cdf_SM2GMIN->AddText("Disfavor SM");
    pt_cdf_SM2GMIN->Draw();

    roostr = "canv_exSM_cdf.png";
    canv_h1_cdf_dchi2_SM2GMIN->SaveAs(roostr);
  }
  
}

