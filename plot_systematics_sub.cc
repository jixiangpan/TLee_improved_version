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

#include "TPrincipal.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

#include "./src/draw.icc"

void plot_systematics_sub()
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////////

  TString roostr = "";

  roostr = "file_collapsed_covariance_matrix.root";
  //roostr = "file_collapsed_covariance_matrix_DetNoRandom.root";
  TFile *roofile_syst = new TFile(roostr, "read");

  TMatrixD* matrix_absolute_flux_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_flux_cov_newworld");
  TMatrixD* matrix_absolute_Xs_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_Xs_cov_newworld");
  TMatrixD* matrix_absolute_detector_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_detector_cov_newworld");
  TMatrixD* matrix_absolute_mc_stat_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_mc_stat_cov_newworld");
  TMatrixD* matrix_absolute_additional_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_additional_cov_newworld");

  map<int, TMatrixD*>matrix_fgx_sub;
  for(int idx=1; idx<=17; idx++) {
    roostr = TString::Format("matrix_sub_flux_geant4_Xs_newworld_%d", idx);
    matrix_fgx_sub[idx] = (TMatrixD*)roofile_syst->Get(roostr);
  }
  map<int, TString>map_name_fgx_sub;
  map_name_fgx_sub[1] = "Skin Effect";
  map_name_fgx_sub[2] = "Horn Current";
  map_name_fgx_sub[3] = "K^{-} Production";
  map_name_fgx_sub[4] = "K^{+} Production";
  map_name_fgx_sub[5] = "K^{0} Production";
  map_name_fgx_sub[6] = "Nucleon Inelastic Xs";
  map_name_fgx_sub[7] = "Nucleon QE Xs";
  map_name_fgx_sub[8] = "Nucleon Total Xs";
  map_name_fgx_sub[9] = "#pi^{-} Production";
  map_name_fgx_sub[10] = "Pion Inelastic Xs";
  map_name_fgx_sub[11] = "Pion QE Xs";
  map_name_fgx_sub[12] = "Pion Total Xs";
  map_name_fgx_sub[13] = "#pi^{+} Production";
  map_name_fgx_sub[14] = "#pi^{+} (Geant4)";
  map_name_fgx_sub[15] = "#pi^{-} (Geant4)";
  map_name_fgx_sub[16] = "proton (Geant4)";
  map_name_fgx_sub[17] = "#nu-Ar interaction (GENIE)";

  map<int, TMatrixD*>matrix_detector_sub;
  for(int idx=1; idx<=10; idx++) {
    roostr = TString::Format("matrix_absolute_detector_sub_cov_newworld_%02d", idx);
    matrix_detector_sub[idx] = (TMatrixD*)roofile_syst->Get(roostr);
  }

  TMatrixD* matrix_pred_newworld = (TMatrixD*)roofile_syst->Get("matrix_pred_newworld");// 1xN
  TMatrixD* matrix_data_newworld = (TMatrixD*)roofile_syst->Get("matrix_data_newworld");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int color_flux     = kRed;
  int color_Xs       = kBlue;
  int color_geant    = kAzure+7;
  int color_detector = kMagenta;
  int color_dirt     = kOrange-3;
  int color_mc_stat  = kGreen+1;
  int color_total    = kBlack;
  
  const int num_ch = 7;
  int bins_ch[num_ch] = {26, 26, 26, 26, 11, 11, 11};  
  int rows = matrix_absolute_additional_cov_newworld->GetNrows();

  // for(int idx=0; idx<rows; idx++) {
  //   double val_flux = (*matrix_absolute_flux_cov_newworld)(idx, idx);
  //   double val_Xs = (*matrix_absolute_Xs_cov_newworld)(idx, idx);
  //   double sum_sub = 0;
  //   for(int jdx=1; jdx<=17; jdx++) sum_sub += (*matrix_fgx_sub[jdx])(idx, idx);
  //   if( fabs(sum_sub-val_flux-val_Xs)>1e-3 ) cout<<" -----------> diff "<<idx<<endl;
  // }

  TLine *line_percentage[num_ch-1];
  for(int idx=0; idx<num_ch-1; idx++) {
    double xx = 0;
    for(int jdx=0; jdx<=idx; jdx++) xx += bins_ch[jdx];
    line_percentage[idx] = new TLine( xx, 0, xx, 110);
    line_percentage[idx]->SetLineStyle(7);
    line_percentage[idx]->SetLineWidth(1);
  }
  
  map<int, TPaveText*>pt_text_ch;
  TString pt_str_ch[num_ch] = {"FC #nu_{e}CC", "PC #nu_{e}CC", "FC #nu_{#mu}CC", "PC #nu_{#mu}CC", "FC CC#pi^{0}", "PC CC#pi^{0}", "NC#pi^{0}"};
  double pt_str_angle[num_ch] = {0, 0, 0, 0,25,25, 25};
  for(int idx=0; idx<num_ch; idx++) {
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += bins_ch[jdx];
    pt_text_ch[idx] = new TPaveText( line_eff+4, 110+2.5, line_eff+4+1, 110+2.5, "l");
    pt_text_ch[idx]->SetTextSize(0.04);
    pt_text_ch[idx]->SetTextFont(42); pt_text_ch[idx]->SetTextAlign(11);
    pt_text_ch[idx]->SetBorderSize(0); pt_text_ch[idx]->SetFillStyle(0);
    pt_text_ch[idx]->AddText( pt_str_ch[idx] );
    ((TText*)pt_text_ch[idx]->GetListOfLines()->Last())->SetTextAngle( pt_str_angle[idx] );        
  }

  TString *axis_label_str = new TString[rows];
  int line_str = 0;  
  map<int, TLine*>map_line_label_xx;
  map<int, TLine*>map_line_label_yy;
  int line_map_line = 0;  
  for(int ich=0; ich<num_ch; ich++) {
    for(int ibin=1; ibin<=bins_ch[ich]; ibin++) {
      line_str++;
      if(ibin==5 || ibin==10 || ibin==15 || ibin==20 || ibin==25) {
	axis_label_str[line_str-1] = TString::Format("%d", ibin*100);

	line_map_line++;
	map_line_label_xx[line_map_line] = new TLine(line_str, 0, line_str, 3);
	map_line_label_yy[line_map_line] = new TLine(0, line_str, 3, line_str);	
      }
      else {
	axis_label_str[line_str-1] ="";
      }
    }
  }


  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// flux
  /////////////////////////////////////////////////////////////////////////////////////////////////////// flux

  int flux_bgn = 1;
  int flux_end = 13;
  
  // map_name_fgx_sub[1] = "Skin Effect";
  // map_name_fgx_sub[2] = "Horn Current";
  // map_name_fgx_sub[3] = "K^{-} Production";
  // map_name_fgx_sub[4] = "K^{+} Production";
  // map_name_fgx_sub[5] = "K^{0} Production";
  // map_name_fgx_sub[6] = "Nucleon Inelastic Xs";
  // map_name_fgx_sub[7] = "Nucleon QE Xs";
  // map_name_fgx_sub[8] = "Nucleon Total Xs";
  // map_name_fgx_sub[9] = "#pi^{-} Production";
  // map_name_fgx_sub[10] = "Pion Inelastic Xs";
  // map_name_fgx_sub[11] = "Pion QE Xs";
  // map_name_fgx_sub[12] = "Pion Total Xs";
  // map_name_fgx_sub[13] = "#pi^{+} Production";
  map<int, int>map_color_flux;
  
  map_color_flux[1] = kRed;
  map_color_flux[4] = kGreen;
  map_color_flux[5] = kBlue;
  map_color_flux[7] = kMagenta;
  map_color_flux[13] = kOrange-3;

  
  //////
  roostr = "h2_sum_flux_abs_cov"; TH2D *h2_sum_flux_abs_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_flux_rel_cov"; TH2D *h2_sum_flux_rel_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_flux_correlation"; TH2D *h2_sum_flux_correlation = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  
  map<int, TH2D*>h2_sub_flux_abs_cov;
  map<int, TH2D*>h2_sub_flux_rel_cov;
  map<int, TH2D*>h2_sub_flux_correlation;  
  for(int idx=flux_bgn; idx<=flux_end; idx++) {
    roostr = TString::Format("h2_sub_flux_abs_cov_%02d", idx); h2_sub_flux_abs_cov[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
    roostr = TString::Format("h2_sub_flux_rel_cov_%02d", idx); h2_sub_flux_rel_cov[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
    roostr = TString::Format("h2_sub_flux_correlation_%02d", idx); h2_sub_flux_correlation[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  }// idx

  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {      
      double cv_i = (*matrix_pred_newworld)(0, ibin-1);
      double cv_j = (*matrix_pred_newworld)(0, jbin-1);
      
      double total_cov_ij = 0;
      double total_cov_i = 0;
      double total_cov_j = 0;

      for(int idx=flux_bgn; idx<=flux_end; idx++) {
	double sub_cov = (*matrix_fgx_sub[idx])(ibin-1, jbin-1);	
	double cov_i = (*matrix_fgx_sub[idx])(ibin-1, ibin-1);
	double cov_j = (*matrix_fgx_sub[idx])(jbin-1, jbin-1);

	total_cov_ij += sub_cov;
	total_cov_i += cov_i;
	total_cov_j += cov_j;
	
	double rel_cov = 0;
	double correlation = 0;
	if(cv_i==0 || cv_j==0) {
	  if( ibin==jbin ) correlation = 1;
	}
	else {
	  rel_cov = sub_cov/cv_i/cv_j;
	  correlation = sub_cov/sqrt(cov_i)/sqrt(cov_j);
	}

	h2_sub_flux_abs_cov[idx]->SetBinContent(ibin, jbin, sub_cov);
	h2_sub_flux_rel_cov[idx]->SetBinContent(ibin, jbin, rel_cov);
	h2_sub_flux_correlation[idx]->SetBinContent(ibin, jbin, correlation);
	
      }// idx

      
      double rel_cov = 0;
      double correlation = 0;      
      if(cv_i==0 || cv_j==0) {
	if( ibin==jbin ) correlation = 1;
      }
      else {
	rel_cov = total_cov_ij/cv_i/cv_j;
	correlation = total_cov_ij/sqrt(total_cov_i)/sqrt(total_cov_j);
      }
            
      h2_sum_flux_abs_cov->SetBinContent( ibin, jbin, total_cov_ij );
      h2_sum_flux_rel_cov->SetBinContent( ibin, jbin, rel_cov );
      h2_sum_flux_correlation->SetBinContent(ibin, jbin, correlation);
      
    }// jbin
  }// ibin


  TH1D *h1_rel_flux = new TH1D("h1_rel_flux", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_rel_flux->SetBinContent( ibin, sqrt(h2_sum_flux_rel_cov->GetBinContent(ibin, ibin)) );
  
  THStack *h1_stack_flux = new THStack("h1_stack_flux", "");
  map<int, TH1D*>h1_frac_sub_flux;
  for(int idx=flux_bgn; idx<=flux_end; idx++) {
    roostr = TString::Format("h1_frac_sub_flux_%02d", idx);
    h1_frac_sub_flux[idx] = new TH1D(roostr, roostr, rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      double total_cov = h2_sum_flux_abs_cov->GetBinContent(ibin, ibin);
      double sub_cov = h2_sub_flux_abs_cov[idx]->GetBinContent(ibin, ibin);
      double val_frac = 0;
      if( total_cov!=0 ) val_frac = sub_cov/total_cov * 100;
      h1_frac_sub_flux[idx]->SetBinContent(ibin, val_frac);      
    }// ibin
    
    //h1_stack_flux->Add( h1_frac_sub_flux[idx] );
    h1_frac_sub_flux[idx]->SetLineColor(kBlack);
    h1_frac_sub_flux[idx]->SetFillColor( map_color_flux[idx] );
  }// idx

  
  for(int idx=flux_bgn; idx<=flux_end; idx++) {
    if( map_color_flux[idx]!=0 ) h1_stack_flux->Add( h1_frac_sub_flux[idx] );
  }

  for(int idx=flux_bgn; idx<=flux_end; idx++) {
    if( map_color_flux[idx]==0 ) h1_stack_flux->Add( h1_frac_sub_flux[idx] );
  }

  ///////////////////////////////
  ///////////////////////////////

  TH2D *h2_basic_stack_flux = new TH2D("h2_basic_stack_flux", "" , rows, 0, rows, 110, 0, 110);
  
  TCanvas *canv_h2_basic_stack_flux = new TCanvas("canv_h2_basic_stack_flux", "canv_h2_basic_stack_flux", 1300, 700);
  func_canv_margin(canv_h2_basic_stack_flux, 0.1, 0.24, 0.11, 0.15);
  h2_basic_stack_flux->Draw();
  func_title_size(h2_basic_stack_flux, 0.06, 0.045, 0.045, 0.045);
  //h2_basic_stack_flux->GetXaxis()->LabelsOption("v R");
  //h2_basic_stack_flux->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_basic_stack_flux, "Reco energy [MeV]", "Syst. percentage [%]");
  h2_basic_stack_flux->GetXaxis()->CenterTitle(); h2_basic_stack_flux->GetYaxis()->CenterTitle();
  h2_basic_stack_flux->GetXaxis()->SetTitleOffset(1.8); h2_basic_stack_flux->GetYaxis()->SetTitleOffset(0.9);

  h1_stack_flux->Draw("same");

  h2_basic_stack_flux->Draw("same axis");
  h2_basic_stack_flux->GetYaxis()->SetTickLength(0.02);
  h2_basic_stack_flux->GetXaxis()->SetTickLength(0);
  h2_basic_stack_flux->GetXaxis()->SetLabelSize(0.065);
  
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_basic_stack_flux->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw();
  }
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_percentage[idx]->Draw();
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }


  int index_flux = 0;
  TLegend *lg_fraction_flux = new TLegend(0.76+0.01, 0.15, 0.98, 0.89);
  lg_fraction_flux->Draw();
  lg_fraction_flux->SetTextSize(0.045);

  lg_fraction_flux->AddEntry("", "Flux Syst", "");
  
  index_flux = 13;
  lg_fraction_flux->AddEntry(h1_frac_sub_flux[index_flux], TString::Format("#color[%d]{%s}", map_color_flux[index_flux], map_name_fgx_sub[index_flux].Data()), "f");  
  index_flux = 5;
  lg_fraction_flux->AddEntry(h1_frac_sub_flux[index_flux], TString::Format("#color[%d]{%s}", map_color_flux[index_flux], map_name_fgx_sub[index_flux].Data()), "f");  
  index_flux = 4;
  lg_fraction_flux->AddEntry(h1_frac_sub_flux[index_flux], TString::Format("#color[%d]{%s}", map_color_flux[index_flux], map_name_fgx_sub[index_flux].Data()), "f");
  index_flux = 7;
  lg_fraction_flux->AddEntry(h1_frac_sub_flux[index_flux], TString::Format("#color[%d]{%s}", map_color_flux[index_flux], map_name_fgx_sub[index_flux].Data()), "f");  
  index_flux = 1;
  lg_fraction_flux->AddEntry(h1_frac_sub_flux[index_flux], TString::Format("#color[%d]{%s}", map_color_flux[index_flux], map_name_fgx_sub[index_flux].Data()), "f");
  index_flux = 2;
  lg_fraction_flux->AddEntry(h1_frac_sub_flux[index_flux], "Others", "f");
  
  canv_h2_basic_stack_flux->SaveAs("canv_h2_basic_stack_flux.png");
  //canv_h2_basic_stack_flux->SaveAs("canv_h2_basic_stack_flux.pdf");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// geant
  /////////////////////////////////////////////////////////////////////////////////////////////////////// geant

  int geant_bgn = 14;
  int geant_end = 16;
  
  map<int, int>map_color_geant;
  
  map_color_geant[14] = kRed;
  map_color_geant[15] = kGreen;
  map_color_geant[16] = kOrange-3;

  //////
  roostr = "h2_sum_geant_abs_cov"; TH2D *h2_sum_geant_abs_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_geant_rel_cov"; TH2D *h2_sum_geant_rel_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_geant_correlation"; TH2D *h2_sum_geant_correlation = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  
  map<int, TH2D*>h2_sub_geant_abs_cov;
  map<int, TH2D*>h2_sub_geant_rel_cov;
  map<int, TH2D*>h2_sub_geant_correlation;  
  for(int idx=geant_bgn; idx<=geant_end; idx++) {
    roostr = TString::Format("h2_sub_geant_abs_cov_%02d", idx); h2_sub_geant_abs_cov[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
    roostr = TString::Format("h2_sub_geant_rel_cov_%02d", idx); h2_sub_geant_rel_cov[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
    roostr = TString::Format("h2_sub_geant_correlation_%02d", idx); h2_sub_geant_correlation[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  }// idx

  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {      
      double cv_i = (*matrix_pred_newworld)(0, ibin-1);
      double cv_j = (*matrix_pred_newworld)(0, jbin-1);
      
      double total_cov_ij = 0;
      double total_cov_i = 0;
      double total_cov_j = 0;

      for(int idx=geant_bgn; idx<=geant_end; idx++) {
	double sub_cov = (*matrix_fgx_sub[idx])(ibin-1, jbin-1);	
	double cov_i = (*matrix_fgx_sub[idx])(ibin-1, ibin-1);
	double cov_j = (*matrix_fgx_sub[idx])(jbin-1, jbin-1);

	total_cov_ij += sub_cov;
	total_cov_i += cov_i;
	total_cov_j += cov_j;
	
	double rel_cov = 0;
	double correlation = 0;
	if(cv_i==0 || cv_j==0) {
	  if( ibin==jbin ) correlation = 1;
	}
	else {
	  rel_cov = sub_cov/cv_i/cv_j;
	  correlation = sub_cov/sqrt(cov_i)/sqrt(cov_j);
	}

	h2_sub_geant_abs_cov[idx]->SetBinContent(ibin, jbin, sub_cov);
	h2_sub_geant_rel_cov[idx]->SetBinContent(ibin, jbin, rel_cov);
	h2_sub_geant_correlation[idx]->SetBinContent(ibin, jbin, correlation);
	
      }// idx

      
      double rel_cov = 0;
      double correlation = 0;      
      if(cv_i==0 || cv_j==0) {
	if( ibin==jbin ) correlation = 1;
      }
      else {
	rel_cov = total_cov_ij/cv_i/cv_j;
	correlation = total_cov_ij/sqrt(total_cov_i)/sqrt(total_cov_j);
      }
            
      h2_sum_geant_abs_cov->SetBinContent( ibin, jbin, total_cov_ij );
      h2_sum_geant_rel_cov->SetBinContent( ibin, jbin, rel_cov );
      h2_sum_geant_correlation->SetBinContent(ibin, jbin, correlation);
      
    }// jbin
  }// ibin


  TH1D *h1_rel_geant = new TH1D("h1_rel_geant", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_rel_geant->SetBinContent( ibin, sqrt(h2_sum_geant_rel_cov->GetBinContent(ibin, ibin)) );
  
  THStack *h1_stack_geant = new THStack("h1_stack_geant", "");
  map<int, TH1D*>h1_frac_sub_geant;
  for(int idx=geant_bgn; idx<=geant_end; idx++) {
    roostr = TString::Format("h1_frac_sub_geant_%02d", idx);
    h1_frac_sub_geant[idx] = new TH1D(roostr, roostr, rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      double total_cov = h2_sum_geant_abs_cov->GetBinContent(ibin, ibin);
      double sub_cov = h2_sub_geant_abs_cov[idx]->GetBinContent(ibin, ibin);
      double val_frac = 0;
      if( total_cov!=0 ) val_frac = sub_cov/total_cov * 100;
      h1_frac_sub_geant[idx]->SetBinContent(ibin, val_frac);      
    }// ibin
    
    //h1_stack_geant->Add( h1_frac_sub_geant[idx] );
    h1_frac_sub_geant[idx]->SetLineColor(kBlack);
    h1_frac_sub_geant[idx]->SetFillColor( map_color_geant[idx] );
  }// idx

  
  for(int idx=geant_bgn; idx<=geant_end; idx++) {
    if( map_color_geant[idx]==0 ) h1_stack_geant->Add( h1_frac_sub_geant[idx] );
  }

  for(int idx=geant_bgn; idx<=geant_end; idx++) {
    if( map_color_geant[idx]!=0 ) h1_stack_geant->Add( h1_frac_sub_geant[idx] );
  }

  ///////////////////////////////
  ///////////////////////////////

  TH2D *h2_basic_stack_geant = new TH2D("h2_basic_stack_geant", "" , rows, 0, rows, 110, 0, 110);
  
  TCanvas *canv_h2_basic_stack_geant = new TCanvas("canv_h2_basic_stack_geant", "canv_h2_basic_stack_geant", 1300, 700);
  func_canv_margin(canv_h2_basic_stack_geant, 0.1, 0.24, 0.11, 0.15);
  h2_basic_stack_geant->Draw();
  func_title_size(h2_basic_stack_geant, 0.06, 0.045, 0.045, 0.045);
  //h2_basic_stack_geant->GetXaxis()->LabelsOption("v R");
  //h2_basic_stack_geant->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_basic_stack_geant, "Reco energy [MeV]", "Syst. percentage [%]");
  h2_basic_stack_geant->GetXaxis()->CenterTitle(); h2_basic_stack_geant->GetYaxis()->CenterTitle();
  h2_basic_stack_geant->GetXaxis()->SetTitleOffset(1.8); h2_basic_stack_geant->GetYaxis()->SetTitleOffset(0.9);

  h1_stack_geant->Draw("same");

  h2_basic_stack_geant->Draw("same axis");
  h2_basic_stack_geant->GetYaxis()->SetTickLength(0.02);
  h2_basic_stack_geant->GetXaxis()->SetTickLength(0);
  h2_basic_stack_geant->GetXaxis()->SetLabelSize(0.065);

  for(int ibin=1; ibin<=rows; ibin++) {
    h2_basic_stack_geant->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw();
  }
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_percentage[idx]->Draw();
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }


  int index_geant = 0;
  TLegend *lg_fraction_geant = new TLegend(0.76+0.01, 0.45, 0.98, 0.89);
  lg_fraction_geant->Draw();
  lg_fraction_geant->SetTextSize(0.045);

  lg_fraction_geant->AddEntry("", "Re-interaction", "");
  
  index_geant = 14;
  lg_fraction_geant->AddEntry(h1_frac_sub_geant[index_geant], TString::Format("#color[%d]{%s}", map_color_geant[index_geant], map_name_fgx_sub[index_geant].Data()), "f");  
  index_geant = 15;
  lg_fraction_geant->AddEntry(h1_frac_sub_geant[index_geant], TString::Format("#color[%d]{%s}", map_color_geant[index_geant], map_name_fgx_sub[index_geant].Data()), "f");  
  index_geant = 16;
  lg_fraction_geant->AddEntry(h1_frac_sub_geant[index_geant], TString::Format("#color[%d]{%s}", map_color_geant[index_geant], map_name_fgx_sub[index_geant].Data()), "f");
 
  canv_h2_basic_stack_geant->SaveAs("canv_h2_basic_stack_geant.png");
  //canv_h2_basic_stack_geant->SaveAs("canv_h2_basic_stack_geant.pdf");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// detector
  /////////////////////////////////////////////////////////////////////////////////////////////////////// detector

  int detector_bgn = 1;
  int detector_end = 10;
  
  map<int, TString>map_name_detector_sub;
  map_name_detector_sub[1] = "LY Down";
  map_name_detector_sub[2] = "LY Rayleigh";
  map_name_detector_sub[3] = "Recomb2";
  map_name_detector_sub[4] = "SCE";
  //map_name_detector_sub[5] = ;
  map_name_detector_sub[6] = "WireMod #theta_{xz}";
  map_name_detector_sub[7] = "WireMod #theta_{yz}";
  map_name_detector_sub[8] = "WireMod x";
  map_name_detector_sub[9] = "WireMod y";
  map_name_detector_sub[10]= "LY Att";
  
  map<int, int>map_color_detector;
  map_color_detector[1] = 2;
  map_color_detector[2] = 3;
  map_color_detector[3] = 4;
  map_color_detector[4] = 5;
  //map_color_detector[5] = ;
  map_color_detector[6] = 6;
  map_color_detector[7] = 7;
  map_color_detector[8] = 8;
  map_color_detector[9] = 9;
  map_color_detector[10] = kOrange-3;
  
  //////
  roostr = "h2_sum_detector_abs_cov"; TH2D *h2_sum_detector_abs_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_detector_rel_cov"; TH2D *h2_sum_detector_rel_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_detector_correlation"; TH2D *h2_sum_detector_correlation = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  
  map<int, TH2D*>h2_sub_detector_abs_cov;
  map<int, TH2D*>h2_sub_detector_rel_cov;
  map<int, TH2D*>h2_sub_detector_correlation;  
  for(int idx=detector_bgn; idx<=detector_end; idx++) {
    roostr = TString::Format("h2_sub_detector_abs_cov_%02d", idx); h2_sub_detector_abs_cov[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
    roostr = TString::Format("h2_sub_detector_rel_cov_%02d", idx); h2_sub_detector_rel_cov[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
    roostr = TString::Format("h2_sub_detector_correlation_%02d", idx); h2_sub_detector_correlation[idx] = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  }// idx

  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {      
      double cv_i = (*matrix_pred_newworld)(0, ibin-1);
      double cv_j = (*matrix_pred_newworld)(0, jbin-1);
      
      double total_cov_ij = 0;
      double total_cov_i = 0;
      double total_cov_j = 0;

      for(int idx=detector_bgn; idx<=detector_end; idx++) {
	if( idx==5 ) continue;
	
	double sub_cov = (*matrix_detector_sub[idx])(ibin-1, jbin-1);	
	double cov_i = (*matrix_detector_sub[idx])(ibin-1, ibin-1);
	double cov_j = (*matrix_detector_sub[idx])(jbin-1, jbin-1);

	total_cov_ij += sub_cov;
	total_cov_i += cov_i;
	total_cov_j += cov_j;
	
	double rel_cov = 0;
	double correlation = 0;
	if(cv_i==0 || cv_j==0) {
	  if( ibin==jbin ) correlation = 1;
	}
	else {
	  rel_cov = sub_cov/cv_i/cv_j;
	  correlation = sub_cov/sqrt(cov_i)/sqrt(cov_j);
	}

	h2_sub_detector_abs_cov[idx]->SetBinContent(ibin, jbin, sub_cov);
	h2_sub_detector_rel_cov[idx]->SetBinContent(ibin, jbin, rel_cov);
	h2_sub_detector_correlation[idx]->SetBinContent(ibin, jbin, correlation);
	
      }// idx

      
      double rel_cov = 0;
      double correlation = 0;      
      if(cv_i==0 || cv_j==0) {
	if( ibin==jbin ) correlation = 1;
      }
      else {
	rel_cov = total_cov_ij/cv_i/cv_j;
	correlation = total_cov_ij/sqrt(total_cov_i)/sqrt(total_cov_j);
      }
            
      h2_sum_detector_abs_cov->SetBinContent( ibin, jbin, total_cov_ij );
      h2_sum_detector_rel_cov->SetBinContent( ibin, jbin, rel_cov );
      h2_sum_detector_correlation->SetBinContent(ibin, jbin, correlation);
      
    }// jbin
  }// ibin


  TH1D *h1_rel_detector = new TH1D("h1_rel_detector", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_rel_detector->SetBinContent( ibin, sqrt(h2_sum_detector_rel_cov->GetBinContent(ibin, ibin)) );
  
  THStack *h1_stack_detector = new THStack("h1_stack_detector", "");
  map<int, TH1D*>h1_frac_sub_detector;
  for(int idx=detector_bgn; idx<=detector_end; idx++) {
    if( idx==5 ) continue;
    
    roostr = TString::Format("h1_frac_sub_detector_%02d", idx);
    h1_frac_sub_detector[idx] = new TH1D(roostr, roostr, rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      double total_cov = h2_sum_detector_abs_cov->GetBinContent(ibin, ibin);
      double sub_cov = h2_sub_detector_abs_cov[idx]->GetBinContent(ibin, ibin);
      double val_frac = 0;
      if( total_cov!=0 ) val_frac = sub_cov/total_cov * 100;
      h1_frac_sub_detector[idx]->SetBinContent(ibin, val_frac);      
    }// ibin
    
    //h1_stack_detector->Add( h1_frac_sub_detector[idx] );
    h1_frac_sub_detector[idx]->SetLineColor(kBlack);
    h1_frac_sub_detector[idx]->SetFillColor( map_color_detector[idx] );
  }// idx

  
  for(int idx=detector_bgn; idx<=detector_end; idx++) {
    if( idx==5 ) continue;
    if( map_color_detector[idx]==0 ) h1_stack_detector->Add( h1_frac_sub_detector[idx] );
  }

  for(int idx=detector_bgn; idx<=detector_end; idx++) {
    if( idx==5 ) continue;
    if( map_color_detector[idx]!=0 ) h1_stack_detector->Add( h1_frac_sub_detector[idx] );
  }

  ///////////////////////////////
  ///////////////////////////////

  TH2D *h2_basic_stack_detector = new TH2D("h2_basic_stack_detector", "" , rows, 0, rows, 110, 0, 110);
  
  TCanvas *canv_h2_basic_stack_detector = new TCanvas("canv_h2_basic_stack_detector", "canv_h2_basic_stack_detector", 1300, 700);
  func_canv_margin(canv_h2_basic_stack_detector, 0.1, 0.24, 0.11, 0.15);
  h2_basic_stack_detector->Draw();
  func_title_size(h2_basic_stack_detector, 0.06, 0.045, 0.045, 0.045);
  //h2_basic_stack_detector->GetXaxis()->LabelsOption("v R");
  //h2_basic_stack_detector->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_basic_stack_detector, "Reco energy [MeV]", "Syst. percentage [%]");
  h2_basic_stack_detector->GetXaxis()->CenterTitle(); h2_basic_stack_detector->GetYaxis()->CenterTitle();
  h2_basic_stack_detector->GetXaxis()->SetTitleOffset(1.8); h2_basic_stack_detector->GetYaxis()->SetTitleOffset(0.9);

  h1_stack_detector->Draw("same");

  h2_basic_stack_detector->Draw("same axis");
  h2_basic_stack_detector->GetYaxis()->SetTickLength(0.02);
  h2_basic_stack_detector->GetXaxis()->SetTickLength(0);
  h2_basic_stack_detector->GetXaxis()->SetLabelSize(0.065);

  for(int ibin=1; ibin<=rows; ibin++) {
    h2_basic_stack_detector->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw();
  }
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_percentage[idx]->Draw();
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }


  int index_detector = 0;
  TLegend *lg_fraction_detector = new TLegend(0.76+0.01, 0.15, 0.98, 0.89);
  lg_fraction_detector->Draw();
  lg_fraction_detector->SetTextSize(0.045);

  lg_fraction_detector->AddEntry("", "Detector Syst", "");

  for(int idx=detector_bgn; idx<=detector_end; idx++) {
    if( idx==5 ) continue;
    index_detector = idx;
    lg_fraction_detector->AddEntry(h1_frac_sub_detector[index_detector], TString::Format("#color[%d]{%s}", map_color_detector[index_detector], map_name_detector_sub[index_detector].Data()), "f");   
  }

  canv_h2_basic_stack_detector->SaveAs("canv_h2_basic_stack_detector.png");
  //canv_h2_basic_stack_detector->SaveAs("canv_h2_basic_stack_detector.pdf");
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  // flux, geant, detector
  // h2_sum_flux_abs_cov
  // h2_sum_flux_rel_cov
  // h2_sum_flux_correlation

  // Xs, MCstat, Dirt
  // TMatrixD* matrix_absolute_Xs_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_Xs_cov_newworld");
  // TMatrixD* matrix_absolute_mc_stat_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_mc_stat_cov_newworld");
  // TMatrixD* matrix_absolute_additional_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_additional_cov_newworld");

  ////////////////////////////
  roostr = "h2_sum_Xs_abs_cov"; TH2D *h2_sum_Xs_abs_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_Xs_rel_cov"; TH2D *h2_sum_Xs_rel_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_Xs_correlation"; TH2D *h2_sum_Xs_correlation = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_Xs_cov_newworld)(ibin-1, jbin-1);
      double cov_i = (*matrix_absolute_Xs_cov_newworld)(ibin-1, ibin-1);
      double cov_j = (*matrix_absolute_Xs_cov_newworld)(jbin-1, jbin-1);

      double cv_i = (*matrix_pred_newworld)(0, ibin-1);
      double cv_j = (*matrix_pred_newworld)(0, jbin-1);
           
      double rel_cov = 0;
      double correlation = 0;      
      if(cv_i==0 || cv_j==0) {
	if( ibin==jbin ) correlation = 1;
      }
      else {
	rel_cov = cov_ij/cv_i/cv_j;
	correlation = cov_ij/sqrt(cov_i)/sqrt(cov_j);
      }

      h2_sum_Xs_abs_cov->SetBinContent(ibin, jbin, cov_ij);
      h2_sum_Xs_rel_cov->SetBinContent(ibin, jbin, rel_cov);
      h2_sum_Xs_correlation->SetBinContent(ibin, jbin, correlation);      
    }
  }
  
  TH1D *h1_rel_Xs = new TH1D("h1_rel_Xs", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_rel_Xs->SetBinContent( ibin, sqrt(h2_sum_Xs_rel_cov->GetBinContent(ibin, ibin)) );
  
  ////////////////////////////
  roostr = "h2_sum_mc_stat_abs_cov"; TH2D *h2_sum_mc_stat_abs_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_mc_stat_rel_cov"; TH2D *h2_sum_mc_stat_rel_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_mc_stat_correlation"; TH2D *h2_sum_mc_stat_correlation = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_mc_stat_cov_newworld)(ibin-1, jbin-1);
      double cov_i = (*matrix_absolute_mc_stat_cov_newworld)(ibin-1, ibin-1);
      double cov_j = (*matrix_absolute_mc_stat_cov_newworld)(jbin-1, jbin-1);

      double cv_i = (*matrix_pred_newworld)(0, ibin-1);
      double cv_j = (*matrix_pred_newworld)(0, jbin-1);
           
      double rel_cov = 0;
      double correlation = 0;      
      if(cv_i==0 || cv_j==0) {
	if( ibin==jbin ) correlation = 1;
      }
      else {
	rel_cov = cov_ij/cv_i/cv_j;
	correlation = cov_ij/sqrt(cov_i)/sqrt(cov_j);
      }

      h2_sum_mc_stat_abs_cov->SetBinContent(ibin, jbin, cov_ij);
      h2_sum_mc_stat_rel_cov->SetBinContent(ibin, jbin, rel_cov);
      h2_sum_mc_stat_correlation->SetBinContent(ibin, jbin, correlation);      
    }
  }
  
  TH1D *h1_rel_mc_stat = new TH1D("h1_rel_mc_stat", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_rel_mc_stat->SetBinContent( ibin, sqrt(h2_sum_mc_stat_rel_cov->GetBinContent(ibin, ibin)) );
  
  ////////////////////////////
  roostr = "h2_sum_dirt_abs_cov"; TH2D *h2_sum_dirt_abs_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_dirt_rel_cov"; TH2D *h2_sum_dirt_rel_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_dirt_correlation"; TH2D *h2_sum_dirt_correlation = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_additional_cov_newworld)(ibin-1, jbin-1);
      double cov_i = (*matrix_absolute_additional_cov_newworld)(ibin-1, ibin-1);
      double cov_j = (*matrix_absolute_additional_cov_newworld)(jbin-1, jbin-1);

      double cv_i = (*matrix_pred_newworld)(0, ibin-1);
      double cv_j = (*matrix_pred_newworld)(0, jbin-1);
           
      double rel_cov = 0;
      double correlation = 0;      
      if(cv_i==0 || cv_j==0) {
	if( ibin==jbin ) correlation = 1;
      }
      else {
	rel_cov = cov_ij/cv_i/cv_j;
	correlation = cov_ij/sqrt(cov_i)/sqrt(cov_j);
      }

      h2_sum_dirt_abs_cov->SetBinContent(ibin, jbin, cov_ij);
      h2_sum_dirt_rel_cov->SetBinContent(ibin, jbin, rel_cov);
      h2_sum_dirt_correlation->SetBinContent(ibin, jbin, correlation);      
    }
  }
  
  TH1D *h1_rel_dirt = new TH1D("h1_rel_dirt", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_rel_dirt->SetBinContent( ibin, sqrt(h2_sum_dirt_rel_cov->GetBinContent(ibin, ibin)) );

  /////////////////////////
  /////////////////////////
  /////////////////////////

  
  // flux, geant, Xs, detector, mc_stat, dirt
  roostr = "h2_sum_total_abs_cov"; TH2D *h2_sum_total_abs_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_total_rel_cov"; TH2D *h2_sum_total_rel_cov = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);
  roostr = "h2_sum_total_correlation"; TH2D *h2_sum_total_correlation = new TH2D(roostr, roostr, rows, 0, rows, rows, 0, rows);

  h2_sum_total_abs_cov->Add( h2_sum_flux_abs_cov );
  h2_sum_total_abs_cov->Add( h2_sum_geant_abs_cov );
  h2_sum_total_abs_cov->Add( h2_sum_Xs_abs_cov );
  h2_sum_total_abs_cov->Add( h2_sum_detector_abs_cov );
  h2_sum_total_abs_cov->Add( h2_sum_mc_stat_abs_cov );
  h2_sum_total_abs_cov->Add( h2_sum_dirt_abs_cov );

  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = h2_sum_total_abs_cov->GetBinContent(ibin, jbin);
      double cov_i  = h2_sum_total_abs_cov->GetBinContent(ibin, ibin);
      double cov_j  = h2_sum_total_abs_cov->GetBinContent(jbin, jbin); 

      double cv_i = (*matrix_pred_newworld)(0, ibin-1);
      double cv_j = (*matrix_pred_newworld)(0, jbin-1);
           
      double rel_cov = 0;
      double correlation = 0;      
      if(cv_i==0 || cv_j==0) {
	if( ibin==jbin ) correlation = 1;
      }
      else {
	rel_cov = cov_ij/cv_i/cv_j;
	correlation = cov_ij/sqrt(cov_i)/sqrt(cov_j);
      }

      h2_sum_total_rel_cov->SetBinContent(ibin, jbin, rel_cov);
      h2_sum_total_correlation->SetBinContent(ibin, jbin, correlation);      
    }
  }
  
  TH1D *h1_rel_total = new TH1D("h1_rel_total", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) h1_rel_total->SetBinContent( ibin, sqrt(h2_sum_total_rel_cov->GetBinContent(ibin, ibin)) );

  ////////////

  THStack *h1_stack_total = new THStack("h1_stack_total", "");
  roostr = "h1_percentage_flux"; TH1D *h1_percentage_flux = new TH1D(roostr, roostr, rows, 0, rows);
  roostr = "h1_percentage_Xs"; TH1D *h1_percentage_Xs = new TH1D(roostr, roostr, rows, 0, rows);
  roostr = "h1_percentage_geant"; TH1D *h1_percentage_geant = new TH1D(roostr, roostr, rows, 0, rows);
  roostr = "h1_percentage_detector"; TH1D *h1_percentage_detector = new TH1D(roostr, roostr, rows, 0, rows);
  roostr = "h1_percentage_mc_stat"; TH1D *h1_percentage_mc_stat = new TH1D(roostr, roostr, rows, 0, rows);
  roostr = "h1_percentage_dirt"; TH1D *h1_percentage_dirt = new TH1D(roostr, roostr, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    double cov_total = h2_sum_total_abs_cov->GetBinContent(ibin, ibin);
    if( cov_total!=0 ) {
      h1_percentage_flux->SetBinContent(ibin, h2_sum_flux_abs_cov->GetBinContent(ibin, ibin)*100./cov_total );
      h1_percentage_Xs->SetBinContent(ibin, h2_sum_Xs_abs_cov->GetBinContent(ibin, ibin)*100./cov_total );
      h1_percentage_geant->SetBinContent(ibin, h2_sum_geant_abs_cov->GetBinContent(ibin, ibin)*100./cov_total );
      h1_percentage_detector->SetBinContent(ibin, h2_sum_detector_abs_cov->GetBinContent(ibin, ibin)*100./cov_total );
      h1_percentage_mc_stat->SetBinContent(ibin, h2_sum_mc_stat_abs_cov->GetBinContent(ibin, ibin)*100./cov_total );
      h1_percentage_dirt->SetBinContent(ibin, h2_sum_dirt_abs_cov->GetBinContent(ibin, ibin)*100./cov_total );      
    }
  }
  
  
  ////////////
  
  TCanvas *canv_h1_rel_total = new TCanvas("canv_h1_rel_total", "canv_h1_rel_total", 1300, 700);
  func_canv_margin(canv_h1_rel_total, 0.1, 0.24, 0.11, 0.15);
  h1_rel_total->Draw();
  h1_rel_total->SetMinimum(0); h1_rel_total->SetMaximum(2);
  func_title_size(h1_rel_total, 0.06, 0.045, 0.045, 0.045);
  func_xy_title(h1_rel_total, "Reco energy [MeV]", "Relative uncertainty");
  h1_rel_total->GetXaxis()->CenterTitle(); h1_rel_total->GetYaxis()->CenterTitle();
  h1_rel_total->GetXaxis()->SetTitleOffset(1.8); h1_rel_total->GetYaxis()->SetTitleOffset(0.9);

  h1_rel_total->Draw(); h1_rel_total->SetLineColor(color_total); h1_rel_total->SetLineWidth(4); 
  h1_rel_total->GetYaxis()->SetTickLength(0.02);
  h1_rel_total->GetXaxis()->SetTickLength(0);
  h1_rel_total->GetXaxis()->SetLabelSize(0.065);
  
  h1_rel_dirt->Draw("same");     h1_rel_dirt->SetLineColor( color_dirt );
  h1_rel_mc_stat->Draw("same");  h1_rel_mc_stat->SetLineColor( color_mc_stat );  
  h1_rel_flux->Draw("same");     h1_rel_flux->SetLineColor( color_flux );
  h1_rel_geant->Draw("same");    h1_rel_geant->SetLineColor( color_geant );
  h1_rel_Xs->Draw("same");       h1_rel_Xs->SetLineColor( color_Xs );    
  h1_rel_detector->Draw("same"); h1_rel_detector->SetLineColor( color_detector );
  
  for(int ibin=1; ibin<=rows; ibin++) { h1_rel_total->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]); }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw(); map_line_label_xx[it_map->first]->SetY2(0.05);
  }
  for(int idx=0; idx<num_ch-1; idx++) { line_percentage[idx]->Draw(); line_percentage[idx]->SetY2(2); }
  for(int idx=0; idx<num_ch; idx++) { pt_text_ch[idx]->Draw(); }

  TLegend *lg_fraction_total = new TLegend(0.76+0.01, 0.45, 0.98, 0.89);
  lg_fraction_total->Draw();
  lg_fraction_total->SetTextSize(0.045);
  lg_fraction_total->AddEntry(h1_rel_total,    TString::Format("#color[%d]{Total}", color_total ), "l");
  lg_fraction_total->AddEntry(h1_rel_flux,     TString::Format("#color[%d]{Flux}", color_flux ), "l");
  lg_fraction_total->AddEntry(h1_rel_Xs,       TString::Format("#color[%d]{#nu-Ar Xs}", color_Xs ), "l");
  lg_fraction_total->AddEntry(h1_rel_geant,    TString::Format("#color[%d]{Reinteraction}", color_geant ), "l");
  lg_fraction_total->AddEntry(h1_rel_detector, TString::Format("#color[%d]{detector}", color_detector ), "l");
  lg_fraction_total->AddEntry(h1_rel_mc_stat,  TString::Format("#color[%d]{MC stat}", color_mc_stat ), "l");
  lg_fraction_total->AddEntry(h1_rel_dirt,     TString::Format("#color[%d]{dirt}", color_dirt ), "l");
  
  canv_h1_rel_total->SaveAs("canv_h1_rel_total.png");
  //canv_h1_rel_total->SaveAs("canv_h1_rel_total.pdf");
  
  ////////////
  
  TCanvas *canv_h1_rel_flux = new TCanvas("canv_h1_rel_flux", "canv_h1_rel_flux", 1300, 700);
  func_canv_margin(canv_h1_rel_flux, 0.1, 0.24, 0.11, 0.15);
  h1_rel_flux->Draw();
  h1_rel_flux->SetMinimum(0); h1_rel_flux->SetMaximum(1);
  func_title_size(h1_rel_flux, 0.06, 0.045, 0.045, 0.045);
  func_xy_title(h1_rel_flux, "Reco energy [MeV]", "Relative uncertainty");
  h1_rel_flux->GetXaxis()->CenterTitle(); h1_rel_flux->GetYaxis()->CenterTitle();
  h1_rel_flux->GetXaxis()->SetTitleOffset(1.8); h1_rel_flux->GetYaxis()->SetTitleOffset(0.9);
  h1_rel_flux->GetYaxis()->SetNdivisions(508);
  
  h1_rel_flux->Draw(); h1_rel_flux->SetLineColor(color_flux);
  h1_rel_flux->GetYaxis()->SetTickLength(0.02);
  h1_rel_flux->GetXaxis()->SetTickLength(0);
  h1_rel_flux->GetXaxis()->SetLabelSize(0.065);
  
  h1_rel_geant->Draw("same");
  h1_rel_Xs->Draw("same");

  for(int ibin=1; ibin<=rows; ibin++) { h1_rel_flux->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]); }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw(); map_line_label_xx[it_map->first]->SetY2(0.03);
  }
  for(int idx=0; idx<num_ch-1; idx++) { line_percentage[idx]->Draw(); line_percentage[idx]->SetY2(1); }
  for(int idx=0; idx<num_ch; idx++) { pt_text_ch[idx]->Draw(); }

  TLegend *lg_rel_flux = new TLegend(0.76+0.01, 0.45+0.2, 0.98, 0.89);
  lg_rel_flux->Draw();
  lg_rel_flux->SetTextSize(0.045);
  lg_rel_flux->AddEntry(h1_rel_flux,     TString::Format("#color[%d]{Flux}", color_flux ), "l");
  lg_rel_flux->AddEntry(h1_rel_Xs,       TString::Format("#color[%d]{#nu-Ar Xs}", color_Xs ), "l");
  lg_rel_flux->AddEntry(h1_rel_geant,    TString::Format("#color[%d]{Reinteraction}", color_geant ), "l");
  
  canv_h1_rel_flux->SaveAs("canv_h1_rel_flux.png");
  //canv_h1_rel_flux->SaveAs("canv_h1_rel_flux.pdf");
    
  ////////////
  
  TCanvas *canv_h1_rel_detector = new TCanvas("canv_h1_rel_detector", "canv_h1_rel_detector", 1300, 700);
  func_canv_margin(canv_h1_rel_detector, 0.1, 0.24, 0.11, 0.15);
  h1_rel_detector->Draw();
  h1_rel_detector->SetMinimum(0); h1_rel_detector->SetMaximum(2);
  func_title_size(h1_rel_detector, 0.06, 0.045, 0.045, 0.045);
  func_xy_title(h1_rel_detector, "Reco energy [MeV]", "Relative uncertainty");
  h1_rel_detector->GetXaxis()->CenterTitle(); h1_rel_detector->GetYaxis()->CenterTitle();
  h1_rel_detector->GetXaxis()->SetTitleOffset(1.8); h1_rel_detector->GetYaxis()->SetTitleOffset(0.9);
  h1_rel_detector->GetYaxis()->SetNdivisions(508);
  
  h1_rel_detector->Draw(); h1_rel_detector->SetLineColor(color_detector);
  h1_rel_detector->GetYaxis()->SetTickLength(0.02);
  h1_rel_detector->GetXaxis()->SetTickLength(0);
  h1_rel_detector->GetXaxis()->SetLabelSize(0.065);
  
  for(int ibin=1; ibin<=rows; ibin++) { h1_rel_detector->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]); }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw(); map_line_label_xx[it_map->first]->SetY2(0.05);
  }
  for(int idx=0; idx<num_ch-1; idx++) { line_percentage[idx]->Draw(); line_percentage[idx]->SetY2(2); }
  for(int idx=0; idx<num_ch; idx++) { pt_text_ch[idx]->Draw(); }

  canv_h1_rel_detector->SaveAs("canv_h1_rel_detector.png");
  //canv_h1_rel_detector->SaveAs("canv_h1_rel_detector.pdf");
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////////// correlation
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  TLine *line_xx_corr[num_ch-1];
  TLine *line_yy_corr[num_ch-1];
  
  for(int idx=0; idx<num_ch-1; idx++) {
    double xx = 0;
    for(int jdx=0; jdx<=idx; jdx++) xx += bins_ch[jdx];
    line_xx_corr[idx] = new TLine( xx, 0, xx, rows); line_xx_corr[idx]->SetLineStyle(7); line_xx_corr[idx]->SetLineWidth(1);
    line_yy_corr[idx] = new TLine( 0, xx, rows, xx); line_yy_corr[idx]->SetLineStyle(7); line_yy_corr[idx]->SetLineWidth(1);
  }
  
  
  map<int, TPaveText*>pt_text_ch_corr;
  //TString pt_str_ch[num_ch] = {"FC #nu_{e}CC", "PC #nu_{e}CC", "FC #nu_{#mu}CC", "PC #nu_{#mu}CC", "FC CC#pi^{0}", "PC CC#pi^{0}", "NC#pi^{0}"};
  //double pt_str_angle[num_ch] = {0, 0, 0, 0,25,25, 25};
  for(int idx=0; idx<num_ch; idx++) {
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += bins_ch[jdx];
    pt_text_ch_corr[idx] = new TPaveText( line_eff+4, rows+2.9, line_eff+4+1, rows+2.9, "l");
    pt_text_ch_corr[idx]->SetTextSize(0.033);
    pt_text_ch_corr[idx]->SetTextFont(42); pt_text_ch_corr[idx]->SetTextAlign(11);
    pt_text_ch_corr[idx]->SetBorderSize(0); pt_text_ch_corr[idx]->SetFillStyle(0);
    pt_text_ch_corr[idx]->AddText( pt_str_ch[idx] );
    ((TText*)pt_text_ch_corr[idx]->GetListOfLines()->Last())->SetTextAngle( pt_str_angle[idx] );        
  }

  /////////////////////////////////// flux
  
  TCanvas *canv_h2_sum_flux_correlation = new TCanvas("canv_h2_sum_flux_correlation", "canv_h2_sum_flux_correlation", 800, 700);
  func_canv_margin(canv_h2_sum_flux_correlation, 0.15, 0.15, 0.1, 0.15);
  canv_h2_sum_flux_correlation->cd();
  canv_h2_sum_flux_correlation->Update();
  h2_sum_flux_correlation->Draw("colz");
  h2_sum_flux_correlation->SetTitle("");
  func_title_size(h2_sum_flux_correlation, 0.05, 0.033, 0.05, 0.033);
  h2_sum_flux_correlation->GetZaxis()->SetLabelSize(0.033);
  h2_sum_flux_correlation->GetZaxis()->SetRangeUser(-1, 1);
  h2_sum_flux_correlation->GetXaxis()->SetTickLength(0);
  h2_sum_flux_correlation->GetYaxis()->SetTickLength(0);
  h2_sum_flux_correlation->GetXaxis()->SetLabelSize(0.05);
  h2_sum_flux_correlation->GetYaxis()->SetLabelSize(0.05);
  func_xy_title(h2_sum_flux_correlation, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_sum_flux_correlation->GetXaxis()->CenterTitle(); h2_sum_flux_correlation->GetYaxis()->CenterTitle();
  h2_sum_flux_correlation->GetXaxis()->SetTitleOffset(2.2); h2_sum_flux_correlation->GetYaxis()->SetTitleOffset(1.9);

  for(int idx=0; idx<num_ch; idx++) { pt_text_ch_corr[idx]->Draw(); pt_text_ch_corr[idx]->SetTextSize(0.033); }  
  for(int idx=0; idx<num_ch-1; idx++) { line_xx_corr[idx]->Draw(); line_yy_corr[idx]->Draw(); }
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_sum_flux_correlation->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_sum_flux_correlation->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {    
    map_line_label_xx[it_map->first]->Draw(); map_line_label_xx[it_map->first]->SetY2(3);
    map_line_label_yy[it_map->first]->Draw();
  }
  
  canv_h2_sum_flux_correlation->SaveAs("canv_h2_sum_flux_correlation.png");
  //canv_h2_sum_flux_correlation->SaveAs("canv_h2_sum_flux_correlation.pdf");
    
  /////////////////////////////////// geant
  
  TCanvas *canv_h2_sum_geant_correlation = new TCanvas("canv_h2_sum_geant_correlation", "canv_h2_sum_geant_correlation", 800, 700);
  func_canv_margin(canv_h2_sum_geant_correlation, 0.15, 0.15, 0.1, 0.15);
  canv_h2_sum_geant_correlation->cd();
  canv_h2_sum_geant_correlation->Update();
  h2_sum_geant_correlation->Draw("colz");
  h2_sum_geant_correlation->SetTitle("");
  func_title_size(h2_sum_geant_correlation, 0.05, 0.033, 0.05, 0.033);
  h2_sum_geant_correlation->GetZaxis()->SetLabelSize(0.033);
  h2_sum_geant_correlation->GetZaxis()->SetRangeUser(-1, 1);
  h2_sum_geant_correlation->GetXaxis()->SetTickLength(0);
  h2_sum_geant_correlation->GetYaxis()->SetTickLength(0);
  h2_sum_geant_correlation->GetXaxis()->SetLabelSize(0.05);
  h2_sum_geant_correlation->GetYaxis()->SetLabelSize(0.05);
  func_xy_title(h2_sum_geant_correlation, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_sum_geant_correlation->GetXaxis()->CenterTitle(); h2_sum_geant_correlation->GetYaxis()->CenterTitle();
  h2_sum_geant_correlation->GetXaxis()->SetTitleOffset(2.2); h2_sum_geant_correlation->GetYaxis()->SetTitleOffset(1.9);

  for(int idx=0; idx<num_ch; idx++) { pt_text_ch_corr[idx]->Draw(); pt_text_ch_corr[idx]->SetTextSize(0.033); }  
  for(int idx=0; idx<num_ch-1; idx++) { line_xx_corr[idx]->Draw(); line_yy_corr[idx]->Draw(); }
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_sum_geant_correlation->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_sum_geant_correlation->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw();
    map_line_label_yy[it_map->first]->Draw();
  }
  
  canv_h2_sum_geant_correlation->SaveAs("canv_h2_sum_geant_correlation.png");
  //canv_h2_sum_geant_correlation->SaveAs("canv_h2_sum_geant_correlation.pdf");
    
  /////////////////////////////////// Xs
  
  TCanvas *canv_h2_sum_Xs_correlation = new TCanvas("canv_h2_sum_Xs_correlation", "canv_h2_sum_Xs_correlation", 800, 700);
  func_canv_margin(canv_h2_sum_Xs_correlation, 0.15, 0.15, 0.1, 0.15);
  canv_h2_sum_Xs_correlation->cd();
  canv_h2_sum_Xs_correlation->Update();
  h2_sum_Xs_correlation->Draw("colz");
  h2_sum_Xs_correlation->SetTitle("");
  func_title_size(h2_sum_Xs_correlation, 0.05, 0.033, 0.05, 0.033);
  h2_sum_Xs_correlation->GetZaxis()->SetLabelSize(0.033);
  h2_sum_Xs_correlation->GetZaxis()->SetRangeUser(-1, 1);
  h2_sum_Xs_correlation->GetXaxis()->SetTickLength(0);
  h2_sum_Xs_correlation->GetYaxis()->SetTickLength(0);
  h2_sum_Xs_correlation->GetXaxis()->SetLabelSize(0.05);
  h2_sum_Xs_correlation->GetYaxis()->SetLabelSize(0.05);
  func_xy_title(h2_sum_Xs_correlation, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_sum_Xs_correlation->GetXaxis()->CenterTitle(); h2_sum_Xs_correlation->GetYaxis()->CenterTitle();
  h2_sum_Xs_correlation->GetXaxis()->SetTitleOffset(2.2); h2_sum_Xs_correlation->GetYaxis()->SetTitleOffset(1.9);

  for(int idx=0; idx<num_ch; idx++) { pt_text_ch_corr[idx]->Draw(); pt_text_ch_corr[idx]->SetTextSize(0.033); }  
  for(int idx=0; idx<num_ch-1; idx++) { line_xx_corr[idx]->Draw(); line_yy_corr[idx]->Draw(); }
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_sum_Xs_correlation->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_sum_Xs_correlation->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw();
    map_line_label_yy[it_map->first]->Draw();
  }
  
  canv_h2_sum_Xs_correlation->SaveAs("canv_h2_sum_Xs_correlation.png");
  //canv_h2_sum_Xs_correlation->SaveAs("canv_h2_sum_Xs_correlation.pdf");
             
  /////////////////////////////////// detector
  
  TCanvas *canv_h2_sum_detector_correlation = new TCanvas("canv_h2_sum_detector_correlation", "canv_h2_sum_detector_correlation", 800, 700);
  func_canv_margin(canv_h2_sum_detector_correlation, 0.15, 0.15, 0.1, 0.15);
  canv_h2_sum_detector_correlation->cd();
  canv_h2_sum_detector_correlation->Update();
  h2_sum_detector_correlation->Draw("colz");
  h2_sum_detector_correlation->SetTitle("");
  func_title_size(h2_sum_detector_correlation, 0.05, 0.033, 0.05, 0.033);
  h2_sum_detector_correlation->GetZaxis()->SetLabelSize(0.033);
  h2_sum_detector_correlation->GetZaxis()->SetRangeUser(-1, 1);
  h2_sum_detector_correlation->GetXaxis()->SetTickLength(0);
  h2_sum_detector_correlation->GetYaxis()->SetTickLength(0);
  h2_sum_detector_correlation->GetXaxis()->SetLabelSize(0.05);
  h2_sum_detector_correlation->GetYaxis()->SetLabelSize(0.05);
  func_xy_title(h2_sum_detector_correlation, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_sum_detector_correlation->GetXaxis()->CenterTitle(); h2_sum_detector_correlation->GetYaxis()->CenterTitle();
  h2_sum_detector_correlation->GetXaxis()->SetTitleOffset(2.2); h2_sum_detector_correlation->GetYaxis()->SetTitleOffset(1.9);

  for(int idx=0; idx<num_ch; idx++) { pt_text_ch_corr[idx]->Draw(); pt_text_ch_corr[idx]->SetTextSize(0.033); }  
  for(int idx=0; idx<num_ch-1; idx++) { line_xx_corr[idx]->Draw(); line_yy_corr[idx]->Draw(); }
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_sum_detector_correlation->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_sum_detector_correlation->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw();
    map_line_label_yy[it_map->first]->Draw();
  }
  
  canv_h2_sum_detector_correlation->SaveAs("canv_h2_sum_detector_correlation.png");
  //canv_h2_sum_detector_correlation->SaveAs("canv_h2_sum_detector_correlation.pdf");
         
  /////////////////////////////////// total
  
  TCanvas *canv_h2_sum_total_correlation = new TCanvas("canv_h2_sum_total_correlation", "canv_h2_sum_total_correlation", 800, 700);
  func_canv_margin(canv_h2_sum_total_correlation, 0.15, 0.15, 0.1, 0.15);
  canv_h2_sum_total_correlation->cd();
  canv_h2_sum_total_correlation->Update();
  h2_sum_total_correlation->Draw("colz");
  h2_sum_total_correlation->SetTitle("");
  func_title_size(h2_sum_total_correlation, 0.05, 0.033, 0.05, 0.033);
  h2_sum_total_correlation->GetZaxis()->SetLabelSize(0.033);
  h2_sum_total_correlation->GetZaxis()->SetRangeUser(-1, 1);
  h2_sum_total_correlation->GetXaxis()->SetTickLength(0);
  h2_sum_total_correlation->GetYaxis()->SetTickLength(0);
  h2_sum_total_correlation->GetXaxis()->SetLabelSize(0.05);
  h2_sum_total_correlation->GetYaxis()->SetLabelSize(0.05);
  func_xy_title(h2_sum_total_correlation, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_sum_total_correlation->GetXaxis()->CenterTitle(); h2_sum_total_correlation->GetYaxis()->CenterTitle();
  h2_sum_total_correlation->GetXaxis()->SetTitleOffset(2.2); h2_sum_total_correlation->GetYaxis()->SetTitleOffset(1.9);

  for(int idx=0; idx<num_ch; idx++) { pt_text_ch_corr[idx]->Draw(); pt_text_ch_corr[idx]->SetTextSize(0.033); }  
  for(int idx=0; idx<num_ch-1; idx++) { line_xx_corr[idx]->Draw(); line_yy_corr[idx]->Draw(); }
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_sum_total_correlation->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_sum_total_correlation->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);  
  }
  for(auto it_map=map_line_label_xx.begin(); it_map!=map_line_label_xx.end(); it_map++) {
    map_line_label_xx[it_map->first]->Draw();
    map_line_label_yy[it_map->first]->Draw();
  }
  
  canv_h2_sum_total_correlation->SaveAs("canv_h2_sum_total_correlation.png");
  //canv_h2_sum_total_correlation->SaveAs("canv_h2_sum_total_correlation.pdf");
        
  
  
}

