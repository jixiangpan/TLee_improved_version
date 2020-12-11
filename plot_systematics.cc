#include "./src/draw.icc"

void plot_systematics()
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
  TMatrixD* matrix_absolute_total_cov_newworld = (TMatrixD*)roofile_syst->Get("matrix_absolute_cov_newworld");

  TMatrixD* matrix_pred_newworld = (TMatrixD*)roofile_syst->Get("matrix_pred_newworld");
  TMatrixD* matrix_data_newworld = (TMatrixD*)roofile_syst->Get("matrix_data_newworld");

  map<int, TMatrixD*>matrix_absolute_detector_sub_cov_newworld;
  for(int idx=1; idx<=9; idx++) {
    if( idx==5 ) continue;
    roostr = TString::Format("matrix_absolute_detector_sub_cov_newworld_%02d", idx);
    matrix_absolute_detector_sub_cov_newworld[idx] = (TMatrixD*)roofile_syst->Get(roostr);
  }
  
  int rows = matrix_absolute_total_cov_newworld->GetNrows();

  for(int ibin=0; ibin<rows; ibin++) {
    double val_flux = (*matrix_absolute_flux_cov_newworld)(ibin, ibin);
    double val_Xs = (*matrix_absolute_Xs_cov_newworld)(ibin, ibin);
    double val_detector = (*matrix_absolute_detector_cov_newworld)(ibin, ibin);
    double val_mc_stat = (*matrix_absolute_mc_stat_cov_newworld)(ibin, ibin);
    double val_additional = (*matrix_absolute_additional_cov_newworld)(ibin, ibin);
    double val_total = (*matrix_absolute_total_cov_newworld)(ibin, ibin);
    // cout<<Form(" ---> %3d, total %12.4f, check %12.4f", ibin+1, val_total,
    // 	       val_flux + val_Xs + val_detector + val_mc_stat + val_additional
    // 	       )<<endl;
  }

  ///////////////////////////////////////////////////////////////

  int color_flux = kRed;
  int color_Xs = kBlue;
  int color_detector = kMagenta;
  int color_additional = kOrange-3;
  int color_mc_stat = kGreen+1;
  int color_total = kBlack;
  
  const int num_ch = 7;
  int nbins_ch[num_ch] = {26, 26, 26, 26, 11, 11, 11};
  
  double *axis_label = new double[rows+1];
  int line_axis = 0;
  
  TString *axis_label_str = new TString[rows];
  int line_str = 0;
  
  for(int ich=0; ich<num_ch; ich++) {

    if( ich==0 ) {
      for(int ibin=0; ibin<=nbins_ch[ich]; ibin++) {
	line_axis++;
	axis_label[line_axis-1] = ibin;
      }
    }
    else {
      line_axis--;
      for(int ibin=0; ibin<=nbins_ch[ich]; ibin++) {
	line_axis++;
	axis_label[line_axis-1] = ibin;
	if( ich==num_ch-1 && ibin==nbins_ch[ich] ) axis_label[line_axis-1] = 1./0;
      }
    }

    for(int ibin=1; ibin<=nbins_ch[ich]; ibin++) {
      line_str++;
      if(ibin==5 || ibin==10 || ibin==15 || ibin==20 || ibin==25) {
	axis_label_str[line_str-1] = TString::Format("%d", ibin*100);
      }
      else {
	axis_label_str[line_str-1] ="";
      }
    }
    
  }// ich

  cout<<endl<<" ---> line_axis "<<line_axis<<endl<<endl;
  //for(int ibin=0; ibin<=rows; ibin++) cout<<" ---> "<<axis_label[ibin]<<endl;

  /*
  = "a" sort by alphabetic order
  = ">" sort by decreasing values
  = "<" sort by increasing values
  = "h" draw labels horizonthal
  = "v" draw labels vertical
  = "u" draw labels up (end of label right adjusted)
  = "d" draw labels down (start of label left adjusted)

  chopt = "R": labels are Right adjusted on tick mark.(default is centered)
  chopt = "L": labels are Left adjusted on tick mark.
  chopt = "C": labels are Centered on tick mark.
  chopt = "M": In the Middle of the divisions.
  */

  //////////////////
  
  double line_xy[num_ch] = {0};
  for(int idx=0; idx<num_ch; idx++) {
    for(int jdx=0; jdx<=idx; jdx++) {
      line_xy[idx] += nbins_ch[jdx];
    }
    cout<<Form(" ---> line xy: %4.1f", line_xy[idx])<<endl;
  }

  TLine *line_root_xx[num_ch-1];
  TLine *line_root_yy[num_ch-1];
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx] = new TLine( line_xy[idx], 0, line_xy[idx], rows );
    line_root_xx[idx]->SetLineWidth(1);
    line_root_xx[idx]->SetLineColor(kBlack);
    line_root_xx[idx]->SetLineStyle(1);
    line_root_yy[idx] = new TLine( 0, line_xy[idx], rows, line_xy[idx]);
    line_root_yy[idx]->SetLineWidth(1);
    line_root_yy[idx]->SetLineColor(kBlack);
    line_root_yy[idx]->SetLineStyle(1);
  }

  //////////////////

  map<int, int>map_line_axis_xy;
  
  int line_axis_xy = 0; 
  for(int ich=0; ich<num_ch; ich++) {
    for(int ibin=1; ibin<=nbins_ch[ich]; ibin++) {
      line_axis_xy++;
      if(ibin%5==0) {
	map_line_axis_xy[line_axis_xy] = line_axis_xy;
      }
    }    
  }// ich

  map<int, TLine*>map_root_line_axis_xx;
  map<int, TLine*>map_root_line_axis_yy;
  
  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff] = new TLine( line_eff*1., 0, line_eff*1., 3 );
    map_root_line_axis_xx[line_eff]->SetLineWidth(2);
    map_root_line_axis_xx[line_eff]->SetLineColor(kBlack);
    
    map_root_line_axis_yy[line_eff] = new TLine( 0, line_eff*1., 3, line_eff*1.);
    map_root_line_axis_yy[line_eff]->SetLineWidth(2);
    map_root_line_axis_yy[line_eff]->SetLineColor(kBlack);    
  }
  
  //////////////////
  
  map<int, TPaveText*>pt_text_ch;
  TString pt_str_ch[num_ch] = {"FC #nu_{e}CC", "PC #nu_{e}CC", "FC #nu_{#mu}CC", "PC #nu_{#mu}CC", "FC CC#pi^{0}", "PC CC#pi^{0}", "NC#pi^{0}"};
  double pt_str_angle[num_ch] = {0,0,0,0,30,30,30};
  for(int idx=0; idx<num_ch; idx++) {
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += nbins_ch[jdx];
    pt_text_ch[idx] = new TPaveText( line_eff+4, rows+3, line_eff+4+1, rows+3, "l");
    pt_text_ch[idx]->SetTextSize(0.03);
    pt_text_ch[idx]->SetTextFont(42); pt_text_ch[idx]->SetTextAlign(11);
    pt_text_ch[idx]->SetBorderSize(0); pt_text_ch[idx]->SetFillStyle(0);
    pt_text_ch[idx]->AddText( pt_str_ch[idx] );
    ((TText*)pt_text_ch[idx]->GetListOfLines()->Last())->SetTextAngle( pt_str_angle[idx] );        
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// flux
  
  roostr = "h2_covariance_flux";
  TH2D *h2_covariance_flux = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_covariance_flux->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_covariance_flux->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  roostr = "h2_correlation_flux";
  TH2D *h2_correlation_flux = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_correlation_flux->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_correlation_flux->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TH1D *h1_flux_relerr = new TH1D("h1_flux_relerr", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_flux_relerr->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_flux_cov_newworld)(ibin-1,jbin-1);      
      double cov_i = (*matrix_absolute_flux_cov_newworld)(ibin-1,ibin-1);
      double cov_j = (*matrix_absolute_flux_cov_newworld)(jbin-1,jbin-1);      

      double val_correlation = cov_ij/sqrt(cov_i*cov_j);
      if( cov_i==0 || cov_j==0 ) val_correlation = 0;
   
      h2_covariance_flux->SetBinContent(ibin, jbin, cov_ij);
      h2_correlation_flux->SetBinContent(ibin, jbin, val_correlation);

      if( ibin==jbin ) {
	double val_cv = (*matrix_pred_newworld)(0, ibin-1);
	double val_relerr = sqrt(cov_ij)/val_cv;
	h1_flux_relerr->SetBinContent(ibin, val_relerr);
      }
      
    }// jbin
  }// ibin
  
  TCanvas *canv_h2_correlation_flux = new TCanvas("canv_h2_correlation_flux", "canv_h2_correlation_flux", 800, 700);
  func_canv_margin(canv_h2_correlation_flux, 0.15, 0.15, 0.1, 0.15);
  h2_correlation_flux->Draw("colz");
  func_title_size(h2_correlation_flux, 0.05, 0.033, 0.05, 0.033);
  h2_correlation_flux->GetZaxis()->SetLabelSize(0.033);
  h2_correlation_flux->GetZaxis()->SetRangeUser(-1, 1);
  h2_correlation_flux->GetXaxis()->LabelsOption("v R");
  h2_correlation_flux->GetXaxis()->SetTickLength(0);
  h2_correlation_flux->GetYaxis()->SetTickLength(0);
  func_xy_title(h2_correlation_flux, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_correlation_flux->GetXaxis()->CenterTitle(); h2_correlation_flux->GetYaxis()->CenterTitle();
  h2_correlation_flux->GetXaxis()->SetTitleOffset(2.2); h2_correlation_flux->GetYaxis()->SetTitleOffset(1.9);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw("same");
    line_root_yy[idx]->Draw("same");
  }

  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw("same");
    map_root_line_axis_yy[line_eff]->Draw("same");
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }

  canv_h2_correlation_flux->SaveAs("canv_h2_correlation_flux.png");
  canv_h2_correlation_flux->SaveAs("canv_h2_correlation_flux.root");
  h2_correlation_flux->SaveAs("h2_correlation_flux.root");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// Xs

  roostr = "h2_covariance_Xs";
  TH2D *h2_covariance_Xs = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_covariance_Xs->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_covariance_Xs->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  roostr = "h2_correlation_Xs";
  TH2D *h2_correlation_Xs = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_correlation_Xs->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_correlation_Xs->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TH1D *h1_Xs_relerr = new TH1D("h1_Xs_relerr", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_Xs_relerr->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_Xs_cov_newworld)(ibin-1,jbin-1);      
      double cov_i = (*matrix_absolute_Xs_cov_newworld)(ibin-1,ibin-1);
      double cov_j = (*matrix_absolute_Xs_cov_newworld)(jbin-1,jbin-1);      

      double val_correlation = cov_ij/sqrt(cov_i*cov_j);
      if( cov_i==0 || cov_j==0 ) val_correlation = 0;
  
      h2_covariance_Xs->SetBinContent(ibin, jbin, cov_ij);
      h2_correlation_Xs->SetBinContent(ibin, jbin, val_correlation);

      if( ibin==jbin ) {
	double val_cv = (*matrix_pred_newworld)(0, ibin-1);
	double val_relerr = sqrt(cov_ij)/val_cv;
	h1_Xs_relerr->SetBinContent(ibin, val_relerr);
      }
      
    }// jbin
  }// ibin
  
  TCanvas *canv_h2_correlation_Xs = new TCanvas("canv_h2_correlation_Xs", "canv_h2_correlation_Xs", 800, 700);
  func_canv_margin(canv_h2_correlation_Xs, 0.15, 0.15, 0.1, 0.15);
  h2_correlation_Xs->Draw("colz");
  func_title_size(h2_correlation_Xs, 0.05, 0.033, 0.05, 0.033);
  h2_correlation_Xs->GetZaxis()->SetLabelSize(0.033);
  h2_correlation_Xs->GetZaxis()->SetRangeUser(-1, 1);
  h2_correlation_Xs->GetXaxis()->LabelsOption("v R");
  h2_correlation_Xs->GetXaxis()->SetTickLength(0);
  h2_correlation_Xs->GetYaxis()->SetTickLength(0);
  func_xy_title(h2_correlation_Xs, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_correlation_Xs->GetXaxis()->CenterTitle(); h2_correlation_Xs->GetYaxis()->CenterTitle();
  h2_correlation_Xs->GetXaxis()->SetTitleOffset(2.2); h2_correlation_Xs->GetYaxis()->SetTitleOffset(1.9);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw("same");
    line_root_yy[idx]->Draw("same");
  }

  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw("same");
    map_root_line_axis_yy[line_eff]->Draw("same");
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }

  canv_h2_correlation_Xs->SaveAs("canv_h2_correlation_Xs.png");
  canv_h2_correlation_Xs->SaveAs("canv_h2_correlation_Xs.root");
  h2_correlation_Xs->SaveAs("h2_correlation_Xs.root");
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////// detector

  roostr = "h2_covariance_detector";
  TH2D *h2_covariance_detector = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_covariance_detector->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_covariance_detector->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  roostr = "h2_correlation_detector";
  TH2D *h2_correlation_detector = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_correlation_detector->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_correlation_detector->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TH1D *h1_detector_relerr = new TH1D("h1_detector_relerr", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_detector_relerr->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_detector_cov_newworld)(ibin-1,jbin-1);      
      double cov_i = (*matrix_absolute_detector_cov_newworld)(ibin-1,ibin-1);
      double cov_j = (*matrix_absolute_detector_cov_newworld)(jbin-1,jbin-1);      

      double val_correlation = cov_ij/sqrt(cov_i*cov_j);
      if( cov_i==0 || cov_j==0 ) val_correlation = 0;
      // if( fabs(val_correlation)-1<1e-4 ) {
      // 	if( val_correlation>0 ) val_correlation = 0.999;
      // 	if( val_correlation<0 ) val_correlation = -0.999;
      // }
      
      h2_covariance_detector->SetBinContent(ibin, jbin, cov_ij);
      h2_correlation_detector->SetBinContent(ibin, jbin, val_correlation);

      if( ibin==jbin ) {
	double val_cv = (*matrix_pred_newworld)(0, ibin-1);
	double val_relerr = sqrt(cov_ij)/val_cv;
	h1_detector_relerr->SetBinContent(ibin, val_relerr);
      }
      
    }// jbin
  }// ibin
  
  TCanvas *canv_h2_correlation_detector = new TCanvas("canv_h2_correlation_detector", "canv_h2_correlation_detector", 800, 700);
  func_canv_margin(canv_h2_correlation_detector, 0.15, 0.15, 0.1, 0.15);
  h2_correlation_detector->Draw("colz");
  func_title_size(h2_correlation_detector, 0.05, 0.033, 0.05, 0.033);
  h2_correlation_detector->GetZaxis()->SetLabelSize(0.033);
  h2_correlation_detector->GetZaxis()->SetRangeUser(-1, 1);
  h2_correlation_detector->GetXaxis()->LabelsOption("v R");
  h2_correlation_detector->GetXaxis()->SetTickLength(0);
  h2_correlation_detector->GetYaxis()->SetTickLength(0);
  func_xy_title(h2_correlation_detector, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_correlation_detector->GetXaxis()->CenterTitle(); h2_correlation_detector->GetYaxis()->CenterTitle();
  h2_correlation_detector->GetXaxis()->SetTitleOffset(2.2); h2_correlation_detector->GetYaxis()->SetTitleOffset(1.9);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw("same");
    line_root_yy[idx]->Draw("same");
  }

  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw("same");
    map_root_line_axis_yy[line_eff]->Draw("same");
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }

  canv_h2_correlation_detector->SaveAs("canv_h2_correlation_detector.png");


  ///////////////////////////////////////////////////

  map<int, TH2D*>h2_covariance_detector_sub;
  map<int, TH2D*>h2_correlation_detector_sub;
  map<int, TH1D*>h1_detector_relerr_sub;
  map<int, TCanvas*>canv_h2_correlation_detector_sub;
    
  for(auto it=matrix_absolute_detector_sub_cov_newworld.begin(); it!=matrix_absolute_detector_sub_cov_newworld.end(); it++) {
    int idx = it->first;

    
    roostr = "h2_covariance_detector_sub";
    roostr = TString::Format("h2_covariance_detector_sub_%02d", idx);
    h2_covariance_detector_sub[idx] = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      h2_covariance_detector_sub[idx]->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
      h2_covariance_detector_sub[idx]->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    }

    roostr = "h2_correlation_detector";
    roostr = TString::Format("h2_correlation_detector_sub_%02d", idx);
    h2_correlation_detector_sub[idx] = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      h2_correlation_detector_sub[idx]->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
      h2_correlation_detector_sub[idx]->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    }

    roostr = TString::Format("h1_detector_relerr_sub_%02d", idx);
    h1_detector_relerr_sub[idx] = new TH1D(roostr, "", rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      h1_detector_relerr_sub[idx]->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    }
  
    for(int ibin=1; ibin<=rows; ibin++) {
      for(int jbin=1; jbin<=rows; jbin++) {
	double cov_ij = (*matrix_absolute_detector_sub_cov_newworld[idx])(ibin-1,jbin-1);      
	double cov_i = (*matrix_absolute_detector_sub_cov_newworld[idx])(ibin-1,ibin-1);
	double cov_j = (*matrix_absolute_detector_sub_cov_newworld[idx])(jbin-1,jbin-1);      

	double val_correlation = cov_ij/sqrt(cov_i*cov_j);
	if( cov_i==0 || cov_j==0 ) val_correlation = 0;

	if( val_correlation+1>-1e-4 && val_correlation+1<1e-4 ) val_correlation = -0.999;
	if( val_correlation-1>-1e-4 && val_correlation-1<1e-4 ) val_correlation = 0.999;
	
	if( cov_i==0 || cov_j==0 ) val_correlation = 0;

	h2_covariance_detector_sub[idx]->SetBinContent(ibin, jbin, cov_ij);
	h2_correlation_detector_sub[idx]->SetBinContent(ibin, jbin, val_correlation);

	if( ibin==jbin ) {
	  double val_cv = (*matrix_pred_newworld)(0, ibin-1);
	  double val_relerr = sqrt(cov_ij)/val_cv;
	  h1_detector_relerr_sub[idx]->SetBinContent(ibin, val_relerr);
	}
      
      }// jbin
    }// ibin

    roostr = TString::Format("canv_h2_correlation_detector_sub_%02d", idx);
    canv_h2_correlation_detector_sub[idx] = new TCanvas(roostr, roostr, 800, 700);
    func_canv_margin(canv_h2_correlation_detector_sub[idx], 0.15, 0.15, 0.1, 0.15);
    h2_correlation_detector_sub[idx]->Draw("colz");
    func_title_size(h2_correlation_detector_sub[idx], 0.05, 0.033, 0.05, 0.033);
    h2_correlation_detector_sub[idx]->GetZaxis()->SetLabelSize(0.033);
    h2_correlation_detector_sub[idx]->GetZaxis()->SetRangeUser(-1, 1);
    h2_correlation_detector_sub[idx]->GetXaxis()->LabelsOption("v R");
    h2_correlation_detector_sub[idx]->GetXaxis()->SetTickLength(0);
    h2_correlation_detector_sub[idx]->GetYaxis()->SetTickLength(0);
    func_xy_title(h2_correlation_detector_sub[idx], "Reco energy [MeV]", "Reco energy [MeV]");
    h2_correlation_detector_sub[idx]->GetXaxis()->CenterTitle(); h2_correlation_detector_sub[idx]->GetYaxis()->CenterTitle();
    h2_correlation_detector_sub[idx]->GetXaxis()->SetTitleOffset(2.2); h2_correlation_detector_sub[idx]->GetYaxis()->SetTitleOffset(1.9);
  
    for(int idx=0; idx<num_ch-1; idx++) {
      line_root_xx[idx]->Draw("same");
      line_root_yy[idx]->Draw("same");
      line_root_xx[idx]->SetLineColor(kRed);
      line_root_yy[idx]->SetLineColor(kRed);
    }

    for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
      int line_eff = it->first;    
      map_root_line_axis_xx[line_eff]->Draw("same");
      map_root_line_axis_yy[line_eff]->Draw("same");
    }

    for(int idx=0; idx<num_ch; idx++) {
      pt_text_ch[idx]->Draw();
    }

    roostr = TString::Format("canv_h2_correlation_detector_sub_%02d.png", idx);
    canv_h2_correlation_detector_sub[idx]->SaveAs(roostr);

  }
  
            
  /////////////////////////////////////////////////////////////////////////////////////////////////////// mc_stat

  roostr = "h2_covariance_mc_stat";
  TH2D *h2_covariance_mc_stat = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_covariance_mc_stat->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_covariance_mc_stat->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  roostr = "h2_correlation_mc_stat";
  TH2D *h2_correlation_mc_stat = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_correlation_mc_stat->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_correlation_mc_stat->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TH1D *h1_mc_stat_relerr = new TH1D("h1_mc_stat_relerr", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_mc_stat_relerr->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_mc_stat_cov_newworld)(ibin-1,jbin-1);      
      double cov_i = (*matrix_absolute_mc_stat_cov_newworld)(ibin-1,ibin-1);
      double cov_j = (*matrix_absolute_mc_stat_cov_newworld)(jbin-1,jbin-1);      

      double val_correlation = cov_ij/sqrt(cov_i*cov_j);
      if( cov_i==0 || cov_j==0 ) val_correlation = 0;
   
      h2_covariance_mc_stat->SetBinContent(ibin, jbin, cov_ij);
      h2_correlation_mc_stat->SetBinContent(ibin, jbin, val_correlation);

      if( ibin==jbin ) {
	double val_cv = (*matrix_pred_newworld)(0, ibin-1);
	double val_relerr = sqrt(cov_ij)/val_cv;
	h1_mc_stat_relerr->SetBinContent(ibin, val_relerr);
      }
      
    }// jbin
  }// ibin
  
  TCanvas *canv_h2_correlation_mc_stat = new TCanvas("canv_h2_correlation_mc_stat", "canv_h2_correlation_mc_stat", 800, 700);
  func_canv_margin(canv_h2_correlation_mc_stat, 0.15, 0.15, 0.1, 0.15);
  h2_correlation_mc_stat->Draw("colz");
  func_title_size(h2_correlation_mc_stat, 0.05, 0.033, 0.05, 0.033);
  h2_correlation_mc_stat->GetZaxis()->SetLabelSize(0.033);
  h2_correlation_mc_stat->GetZaxis()->SetRangeUser(-1, 1);
  h2_correlation_mc_stat->GetXaxis()->LabelsOption("v R");
  h2_correlation_mc_stat->GetXaxis()->SetTickLength(0);
  h2_correlation_mc_stat->GetYaxis()->SetTickLength(0);
  func_xy_title(h2_correlation_mc_stat, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_correlation_mc_stat->GetXaxis()->CenterTitle(); h2_correlation_mc_stat->GetYaxis()->CenterTitle();
  h2_correlation_mc_stat->GetXaxis()->SetTitleOffset(2.2); h2_correlation_mc_stat->GetYaxis()->SetTitleOffset(1.9);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw("same");
    line_root_yy[idx]->Draw("same");
    line_root_xx[idx]->SetLineColor(kBlack);
    line_root_yy[idx]->SetLineColor(kBlack);
  }

  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw("same");
    map_root_line_axis_yy[line_eff]->Draw("same");
    map_root_line_axis_xx[line_eff]->SetLineColor(kBlack);
    map_root_line_axis_yy[line_eff]->SetLineColor(kBlack);
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }

  canv_h2_correlation_mc_stat->SaveAs("canv_h2_correlation_mc_stat.png");
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////// additional

  roostr = "h2_covariance_additional";
  TH2D *h2_covariance_additional = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_covariance_additional->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_covariance_additional->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  roostr = "h2_correlation_additional";
  TH2D *h2_correlation_additional = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_correlation_additional->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_correlation_additional->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TH1D *h1_additional_relerr = new TH1D("h1_additional_relerr", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_additional_relerr->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_additional_cov_newworld)(ibin-1,jbin-1);      
      double cov_i = (*matrix_absolute_additional_cov_newworld)(ibin-1,ibin-1);
      double cov_j = (*matrix_absolute_additional_cov_newworld)(jbin-1,jbin-1);      

      double val_correlation = cov_ij/sqrt(cov_i*cov_j);
      if( cov_i==0 || cov_j==0 ) val_correlation = 0;
 
      h2_covariance_additional->SetBinContent(ibin, jbin, cov_ij);
      h2_correlation_additional->SetBinContent(ibin, jbin, val_correlation);

      if( ibin==jbin ) {
	double val_cv = (*matrix_pred_newworld)(0, ibin-1);
	double val_relerr = sqrt(cov_ij)/val_cv;
	h1_additional_relerr->SetBinContent(ibin, val_relerr);
      }
      
    }// jbin
  }// ibin
  
  TCanvas *canv_h2_correlation_additional = new TCanvas("canv_h2_correlation_additional", "canv_h2_correlation_additional", 800, 700);
  func_canv_margin(canv_h2_correlation_additional, 0.15, 0.15, 0.1, 0.15);
  h2_correlation_additional->Draw("colz");
  func_title_size(h2_correlation_additional, 0.05, 0.033, 0.05, 0.033);
  h2_correlation_additional->GetZaxis()->SetLabelSize(0.033);
  h2_correlation_additional->GetZaxis()->SetRangeUser(-1, 1);
  h2_correlation_additional->GetXaxis()->LabelsOption("v R");
  h2_correlation_additional->GetXaxis()->SetTickLength(0);
  h2_correlation_additional->GetYaxis()->SetTickLength(0);
  func_xy_title(h2_correlation_additional, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_correlation_additional->GetXaxis()->CenterTitle(); h2_correlation_additional->GetYaxis()->CenterTitle();
  h2_correlation_additional->GetXaxis()->SetTitleOffset(2.2); h2_correlation_additional->GetYaxis()->SetTitleOffset(1.9);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw("same");
    line_root_yy[idx]->Draw("same");
    line_root_xx[idx]->SetLineColor(kBlack);
    line_root_yy[idx]->SetLineColor(kBlack);
  }

  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw("same");
    map_root_line_axis_yy[line_eff]->Draw("same");
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }

  canv_h2_correlation_additional->SaveAs("canv_h2_correlation_additional.png");
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////// total

  roostr = "h2_covariance_total";
  TH2D *h2_covariance_total = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_covariance_total->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_covariance_total->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  roostr = "h2_correlation_total";
  TH2D *h2_correlation_total = new TH2D(roostr, "", rows, 0, rows, rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_correlation_total->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    h2_correlation_total->GetYaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TH1D *h1_total_relerr = new TH1D("h1_total_relerr", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_total_relerr->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double cov_ij = (*matrix_absolute_total_cov_newworld)(ibin-1,jbin-1);      
      double cov_i = (*matrix_absolute_total_cov_newworld)(ibin-1,ibin-1);
      double cov_j = (*matrix_absolute_total_cov_newworld)(jbin-1,jbin-1);      

      double val_correlation = cov_ij/sqrt(cov_i*cov_j);
      if( cov_i==0 || cov_j==0 ) val_correlation = 0;
   
      h2_covariance_total->SetBinContent(ibin, jbin, cov_ij);
      h2_correlation_total->SetBinContent(ibin, jbin, val_correlation);

      if( ibin==jbin ) {
	double val_cv = (*matrix_pred_newworld)(0, ibin-1);
	double val_relerr = sqrt(cov_ij)/val_cv;
	h1_total_relerr->SetBinContent(ibin, val_relerr);
      }
      
    }// jbin
  }// ibin
  
  TCanvas *canv_h2_correlation_total = new TCanvas("canv_h2_correlation_total", "canv_h2_correlation_total", 800, 700);
  func_canv_margin(canv_h2_correlation_total, 0.15, 0.15, 0.1, 0.15);
  h2_correlation_total->Draw("colz");
  func_title_size(h2_correlation_total, 0.05, 0.033, 0.05, 0.033);
  h2_correlation_total->GetZaxis()->SetLabelSize(0.033);
  h2_correlation_total->GetZaxis()->SetRangeUser(-1, 1);
  h2_correlation_total->GetXaxis()->LabelsOption("v R");
  h2_correlation_total->GetXaxis()->SetTickLength(0);
  h2_correlation_total->GetYaxis()->SetTickLength(0);
  func_xy_title(h2_correlation_total, "Reco energy [MeV]", "Reco energy [MeV]");
  h2_correlation_total->GetXaxis()->CenterTitle(); h2_correlation_total->GetYaxis()->CenterTitle();
  h2_correlation_total->GetXaxis()->SetTitleOffset(2.2); h2_correlation_total->GetYaxis()->SetTitleOffset(1.9);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw("same");
    line_root_yy[idx]->Draw("same");
  }

  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw("same");
    map_root_line_axis_yy[line_eff]->Draw("same");
  }

  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
  }

  canv_h2_correlation_total->SaveAs("canv_h2_correlation_total.png");
    
  /////////////////////////////////////////////////////////////////////////////////////////////////////// flux_Xs
  
  roostr = "h2_relerr_flux_Xs";
  TH2D *h2_relerr_flux_Xs = new TH2D(roostr, "", rows, 0, rows, 100, 0, 1.5);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_relerr_flux_Xs->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TCanvas *canv_h2_relerr_flux_Xs = new TCanvas("canv_h2_relerr_flux_Xs", "canv_h2_relerr_flux_Xs", 1300, 700);
  func_canv_margin(canv_h2_relerr_flux_Xs, 0.15, 0.2, 0.11, 0.15);
  h2_relerr_flux_Xs->Draw();
  func_title_size(h2_relerr_flux_Xs, 0.06, 0.042, 0.04, 0.04);
  h2_relerr_flux_Xs->GetXaxis()->LabelsOption("v R");
  h2_relerr_flux_Xs->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_relerr_flux_Xs, "Reco energy [MeV]", "Relative error");
  h2_relerr_flux_Xs->GetXaxis()->CenterTitle(); h2_relerr_flux_Xs->GetYaxis()->CenterTitle();
  h2_relerr_flux_Xs->GetXaxis()->SetTitleOffset(1.9); h2_relerr_flux_Xs->GetYaxis()->SetTitleOffset(1.);
   
  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw();
    map_root_line_axis_xx[line_eff]->SetY2(0.03);
  }
  
  h1_flux_relerr->Draw("same hist");
  h1_flux_relerr->SetLineColor(color_flux);
  
  h1_Xs_relerr->Draw("same hist");
  h1_Xs_relerr->SetLineColor(color_Xs);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw();
    line_root_xx[idx]->SetLineStyle(7);
    line_root_xx[idx]->SetY2(1.5);
  }

  canv_h2_relerr_flux_Xs->cd();
  canv_h2_relerr_flux_Xs->Update();
  for(int idx=0; idx<num_ch; idx++) {    
    pt_text_ch[idx]->Draw();
    pt_text_ch[idx]->SetTextSize(0.037);
    pt_text_ch[idx]->SetY1(1.53);
    pt_text_ch[idx]->SetY2(1.53);
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += nbins_ch[jdx];
    pt_text_ch[idx]->SetX1(line_eff+4);
    pt_text_ch[idx]->SetX2(line_eff+4+1);
  }
  
  TLegend *lg_relerr_flux_Xs = new TLegend(0.82, 0.75, 0.93, 0.89);
  lg_relerr_flux_Xs->AddEntry(h1_flux_relerr, "Flux", "l");
  lg_relerr_flux_Xs->AddEntry(h1_Xs_relerr, "Xs", "l");
  lg_relerr_flux_Xs->Draw();
  lg_relerr_flux_Xs->SetTextSize(0.04);
    
  h2_relerr_flux_Xs->Draw("same axis");
  canv_h2_relerr_flux_Xs->SaveAs("canv_h2_relerr_flux_Xs.png");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// detector
  
  roostr = "h2_relerr_detector";
  TH2D *h2_relerr_detector = new TH2D(roostr, "", rows, 0, rows, 100, 0, 2.5);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_relerr_detector->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TCanvas *canv_h2_relerr_detector = new TCanvas("canv_h2_relerr_detector", "canv_h2_relerr_detector", 1300, 700);
  func_canv_margin(canv_h2_relerr_detector, 0.15, 0.2, 0.11, 0.15);
  h2_relerr_detector->Draw();
  func_title_size(h2_relerr_detector, 0.06, 0.042, 0.04, 0.04);
  h2_relerr_detector->GetXaxis()->LabelsOption("v R");
  h2_relerr_detector->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_relerr_detector, "Reco energy [MeV]", "Relative error");
  h2_relerr_detector->GetXaxis()->CenterTitle(); h2_relerr_detector->GetYaxis()->CenterTitle();
  h2_relerr_detector->GetXaxis()->SetTitleOffset(1.9); h2_relerr_detector->GetYaxis()->SetTitleOffset(1.);
   
  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw();
    map_root_line_axis_xx[line_eff]->SetY2(0.04);
  }
  
  h1_detector_relerr->Draw("same hist");
  h1_detector_relerr->SetLineColor(color_detector);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw();
    line_root_xx[idx]->SetLineStyle(7);
    line_root_xx[idx]->SetY2(2.5);
  }

  canv_h2_relerr_detector->cd();
  canv_h2_relerr_detector->Update();
  for(int idx=0; idx<num_ch; idx++) {    
    pt_text_ch[idx]->Draw();
    pt_text_ch[idx]->SetTextSize(0.037);
    pt_text_ch[idx]->SetY1(2.55);
    pt_text_ch[idx]->SetY2(2.55);
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += nbins_ch[jdx];
    pt_text_ch[idx]->SetX1(line_eff+4);
    pt_text_ch[idx]->SetX2(line_eff+4+1);
  }
  
  // TLegend *lg_relerr_detector = new TLegend(0.82, 0.75, 0.93, 0.89);
  // lg_relerr_detector->AddEntry(h1_flux_relerr, "Flux", "l");
  // lg_relerr_detector->AddEntry(h1_Xs_relerr, "Xs", "l");
  // lg_relerr_detector->Draw();
  // lg_relerr_detector->SetTextSize(0.04);
    
  h2_relerr_detector->Draw("same axis");
  canv_h2_relerr_detector->SaveAs("canv_h2_relerr_detector.png");
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// total
  
  roostr = "h2_relerr_total";
  TH2D *h2_relerr_total = new TH2D(roostr, "", rows, 0, rows, 100, 0, 2.5);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_relerr_total->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }

  TCanvas *canv_h2_relerr_total = new TCanvas("canv_h2_relerr_total", "canv_h2_relerr_total", 1300, 700);
  func_canv_margin(canv_h2_relerr_total, 0.15, 0.2, 0.11, 0.15);
  h2_relerr_total->Draw();
  func_title_size(h2_relerr_total, 0.06, 0.042, 0.04, 0.04);
  h2_relerr_total->GetXaxis()->LabelsOption("v R");
  h2_relerr_total->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_relerr_total, "Reco energy [MeV]", "Relative error");
  h2_relerr_total->GetXaxis()->CenterTitle(); h2_relerr_total->GetYaxis()->CenterTitle();
  h2_relerr_total->GetXaxis()->SetTitleOffset(1.9); h2_relerr_total->GetYaxis()->SetTitleOffset(1.);
   
  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw();
    map_root_line_axis_xx[line_eff]->SetY2(0.05);
  }

  h1_total_relerr->Draw("same hist");
  h1_total_relerr->SetLineColor(color_total);
  h1_total_relerr->SetLineWidth(4);
  
  h1_additional_relerr->Draw("same hist");
  h1_additional_relerr->SetLineColor(color_additional);
  
  h1_mc_stat_relerr->Draw("same hist");
  h1_mc_stat_relerr->SetLineColor(color_mc_stat);
  
  h1_flux_relerr->Draw("same hist");
  h1_flux_relerr->SetLineColor(color_flux);
  
  h1_Xs_relerr->Draw("same hist");
  h1_Xs_relerr->SetLineColor(color_Xs);
  
  h1_detector_relerr->Draw("same hist");
  h1_detector_relerr->SetLineColor(color_detector);
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw();
    line_root_xx[idx]->SetLineStyle(7);
    line_root_xx[idx]->SetY2(2.5);
  }

  canv_h2_relerr_total->cd();
  canv_h2_relerr_total->Update();
  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
    pt_text_ch[idx]->SetY1(2.55);
    pt_text_ch[idx]->SetY2(2.55);
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += nbins_ch[jdx];
    pt_text_ch[idx]->SetX1(line_eff+4);
    pt_text_ch[idx]->SetX2(line_eff+4+1);
  }

  TLegend *lg_relerr_total = new TLegend(0.82, 0.5, 0.93, 0.89);
  lg_relerr_total->AddEntry(h1_total_relerr, "Total", "l");
  lg_relerr_total->AddEntry(h1_flux_relerr, "Flux", "l");
  lg_relerr_total->AddEntry(h1_Xs_relerr, "Xs", "l");
  lg_relerr_total->AddEntry(h1_detector_relerr, "Detector", "l");
  lg_relerr_total->AddEntry(h1_mc_stat_relerr, "MC stat", "l");
  lg_relerr_total->AddEntry(h1_additional_relerr, "Dirt", "l");
  lg_relerr_total->Draw();
  lg_relerr_total->SetTextSize(0.04);
    
  h2_relerr_total->Draw("same axis");
  canv_h2_relerr_total->SaveAs("canv_h2_relerr_total.png");

  /////////////////////////////////////////////////////////////////////////////////////////////////////// fraction
  
  THStack *h1_stack_fraction = new THStack("h1_stack_fraction", "");
  
  TH1D *h1_fraction_flux = new TH1D("h1_fraction_flux", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_fraction_flux->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    double val_flux = h1_flux_relerr->GetBinContent(ibin);
    double val_total = h1_total_relerr->GetBinContent(ibin);
    double val_frac = val_flux*val_flux*100./val_total/val_total;
    h1_fraction_flux->SetBinContent( ibin, val_frac );
  }
  h1_fraction_flux->SetFillColor(color_flux);
  h1_fraction_flux->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_fraction_flux);
  
  TH1D *h1_fraction_Xs = new TH1D("h1_fraction_Xs", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_fraction_Xs->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    double val_Xs = h1_Xs_relerr->GetBinContent(ibin);
    double val_total = h1_total_relerr->GetBinContent(ibin);
    double val_frac = val_Xs*val_Xs*100./val_total/val_total;
    h1_fraction_Xs->SetBinContent( ibin, val_frac );
  }
  h1_fraction_Xs->SetFillColor(color_Xs);
  h1_fraction_Xs->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_fraction_Xs);

  TH1D *h1_fraction_detector = new TH1D("h1_fraction_detector", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_fraction_detector->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    double val_detector = h1_detector_relerr->GetBinContent(ibin);
    double val_total = h1_total_relerr->GetBinContent(ibin);
    double val_frac = val_detector*val_detector*100./val_total/val_total;
    h1_fraction_detector->SetBinContent( ibin, val_frac );
  }
  h1_fraction_detector->SetFillColor(color_detector);
  h1_fraction_detector->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_fraction_detector);

  TH1D *h1_fraction_mc_stat = new TH1D("h1_fraction_mc_stat", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_fraction_mc_stat->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    double val_mc_stat = h1_mc_stat_relerr->GetBinContent(ibin);
    double val_total = h1_total_relerr->GetBinContent(ibin);
    double val_frac = val_mc_stat*val_mc_stat*100./val_total/val_total;
    h1_fraction_mc_stat->SetBinContent( ibin, val_frac );
  }
  h1_fraction_mc_stat->SetFillColor(color_mc_stat);
  h1_fraction_mc_stat->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_fraction_mc_stat);

  TH1D *h1_fraction_additional = new TH1D("h1_fraction_additional", "", rows, 0, rows);
  for(int ibin=1; ibin<=rows; ibin++) {
    h1_fraction_additional->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
    double val_additional = h1_additional_relerr->GetBinContent(ibin);
    double val_total = h1_total_relerr->GetBinContent(ibin);
    double val_frac = val_additional*val_additional*100./val_total/val_total;
    h1_fraction_additional->SetBinContent( ibin, val_frac );
  }
  h1_fraction_additional->SetFillColor(color_additional);
  h1_fraction_additional->SetLineColor(kBlack);
  h1_stack_fraction->Add(h1_fraction_additional);

  ///////////////////////////////
   
  roostr = "h2_basic_fraction";
  TH2D *h2_basic_fraction = new TH2D(roostr, "", rows, 0, rows, 110, 0, 110);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_basic_fraction->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  TCanvas *canv_h2_basic_fraction = new TCanvas("canv_h2_basic_fraction", "canv_h2_basic_fraction", 1300, 700);
  func_canv_margin(canv_h2_basic_fraction, 0.15, 0.2, 0.11, 0.15);
  h2_basic_fraction->Draw();
  func_title_size(h2_basic_fraction, 0.06, 0.042, 0.04, 0.04);
  h2_basic_fraction->GetXaxis()->LabelsOption("v R");
  h2_basic_fraction->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_basic_fraction, "Reco energy [MeV]", "Syst. percentage");
  h2_basic_fraction->GetXaxis()->CenterTitle(); h2_basic_fraction->GetYaxis()->CenterTitle();
  h2_basic_fraction->GetXaxis()->SetTitleOffset(1.9); h2_basic_fraction->GetYaxis()->SetTitleOffset(1.);
   
  h1_stack_fraction->Draw("same");

  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw();
    line_root_xx[idx]->SetLineStyle(7);
    line_root_xx[idx]->SetY2(110);
  }

  canv_h2_basic_fraction->cd();
  canv_h2_basic_fraction->Update();
  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
    pt_text_ch[idx]->SetY1(112);
    pt_text_ch[idx]->SetY2(112);
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += nbins_ch[jdx];
    pt_text_ch[idx]->SetX1(line_eff+4);
    pt_text_ch[idx]->SetX2(line_eff+4+1);
  }

  TLegend *lg_fraction_total = new TLegend(0.82, 0.5+0.11, 0.93, 0.89);
  lg_fraction_total->AddEntry(h1_fraction_flux, "Flux", "f");
  lg_fraction_total->AddEntry(h1_fraction_Xs, "Xs", "f");  
  lg_fraction_total->AddEntry(h1_fraction_detector, "Detector", "f");
  lg_fraction_total->AddEntry(h1_fraction_mc_stat, "MC stat", "f");
  lg_fraction_total->AddEntry(h1_fraction_additional, "Dirt", "f");
  lg_fraction_total->Draw();
  lg_fraction_total->SetTextSize(0.04);
      
  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw();
    map_root_line_axis_xx[line_eff]->SetY2(3);
  }

  h2_basic_fraction->Draw("same axis");
  canv_h2_basic_fraction->SaveAs("canv_h2_basic_fraction.png");

  
  /////////////////////////////////////////////////////////////////////////////////////////////////////// fraction_detector
  
  THStack *h1_stack_fraction_detector = new THStack("h1_stack_fraction_detector", "");

  map<int, int>color_sub;
  color_sub[1] = 2;
  color_sub[2] = 3;
  color_sub[3] = 4;
  color_sub[4] = 5;
  //color_sub[5] = ;
  color_sub[6] = 6;
  color_sub[7] = 7;
  color_sub[8] = 8;
  color_sub[9] = 9;
  color_sub[10] = kOrange-3;
  
  map<int, TString>str_sub;
  str_sub[1] = "LY Down";
  str_sub[2] = "LY Rayleigh";
  str_sub[3] = "Recomb2";
  str_sub[4] = "SCE";
  //str_sub[5] = ;
  str_sub[6] = "WireMod #theta_{xz}";
  str_sub[7] = "WireMod #theta_{yz}";
  str_sub[8] = "WireMod x";
  str_sub[9] = "WireMod y";

  
  map<int, TH1D*>h1_fraction_detecor_sub;
  for(auto it=h1_detector_relerr_sub.begin(); it!=h1_detector_relerr_sub.end(); it++) {
    int idx = it->first;

    roostr = TString::Format("h1_fraction_detecor_sub_%02d", idx);
    h1_fraction_detecor_sub[idx] = new TH1D(roostr, "", rows, 0, rows);
    for(int ibin=1; ibin<=rows; ibin++) {
      h1_fraction_detecor_sub[idx]->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
      double val_cc = h1_detector_relerr_sub[idx]->GetBinContent(ibin);
      double val_total = h1_detector_relerr->GetBinContent(ibin);
      double val_frac = val_cc*val_cc*100./val_total/val_total;      
    
      if( val_cc==0 || val_total==0 ) val_frac = 0;

      h1_fraction_detecor_sub[idx]->SetBinContent( ibin, val_frac );
    }
    h1_fraction_detecor_sub[idx]->SetFillColor(color_sub[idx]);
    h1_fraction_detecor_sub[idx]->SetLineColor(kBlack);
    h1_stack_fraction_detector->Add(h1_fraction_detecor_sub[idx]);
  }
  
  ///////////////////////////////
   
  roostr = "h2_basic_fraction_detector";
  TH2D *h2_basic_fraction_detector = new TH2D(roostr, "", rows, 0, rows, 110, 0, 110);
  for(int ibin=1; ibin<=rows; ibin++) {
    h2_basic_fraction_detector->GetXaxis()->SetBinLabel(ibin, axis_label_str[ibin-1]);
  }
  TCanvas *canv_h2_basic_fraction_detector = new TCanvas("canv_h2_basic_fraction_detector", "canv_h2_basic_fraction_detector", 1300, 700);
  func_canv_margin(canv_h2_basic_fraction_detector, 0.15, 0.2, 0.11, 0.15);
  h2_basic_fraction_detector->Draw();
  func_title_size(h2_basic_fraction_detector, 0.06, 0.042, 0.04, 0.04);
  h2_basic_fraction_detector->GetXaxis()->LabelsOption("v R");
  h2_basic_fraction_detector->GetXaxis()->SetTickLength(0);
  func_xy_title(h2_basic_fraction_detector, "Reco energy [MeV]", "Syst. percentage");
  h2_basic_fraction_detector->GetXaxis()->CenterTitle(); h2_basic_fraction_detector->GetYaxis()->CenterTitle();
  h2_basic_fraction_detector->GetXaxis()->SetTitleOffset(1.9); h2_basic_fraction_detector->GetYaxis()->SetTitleOffset(1.);
   
  h1_stack_fraction_detector->Draw("same");
  
  for(int idx=0; idx<num_ch-1; idx++) {
    line_root_xx[idx]->Draw();
    line_root_xx[idx]->SetLineStyle(7);
    line_root_xx[idx]->SetY2(110);
  }

  canv_h2_basic_fraction_detector->cd();
  canv_h2_basic_fraction_detector->Update();
  for(int idx=0; idx<num_ch; idx++) {
    pt_text_ch[idx]->Draw();
    pt_text_ch[idx]->SetY1(112);
    pt_text_ch[idx]->SetY2(112);
    int line_eff = 0;
    for(int jdx=0; jdx<idx; jdx++) line_eff += nbins_ch[jdx];
    pt_text_ch[idx]->SetX1(line_eff+4);
    pt_text_ch[idx]->SetX2(line_eff+4+1);
  }
  
  TLegend *lg_fraction_detector_total = new TLegend(0.82, 0.15, 0.97, 0.89);
  lg_fraction_detector_total->Draw();
  lg_fraction_detector_total->SetTextSize(0.04);
  for(auto it=h1_detector_relerr_sub.begin(); it!=h1_detector_relerr_sub.end(); it++) {
    int idx = it->first;
    lg_fraction_detector_total->AddEntry(h1_fraction_detecor_sub[idx], str_sub[idx], "f");
  }  
  
  for(auto it=map_line_axis_xy.begin(); it!=map_line_axis_xy.end(); it++) {
    int line_eff = it->first;    
    map_root_line_axis_xx[line_eff]->Draw();
    map_root_line_axis_xx[line_eff]->SetY2(3);
  }  
  
  h2_basic_fraction_detector->Draw("same axis");
  canv_h2_basic_fraction_detector->SaveAs("canv_h2_basic_fraction_detector.png");
}

