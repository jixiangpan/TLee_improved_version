#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "WCPLEEANA/TLee.h"

#include "WCPLEEANA/Configure_Lee.h"

#include "TApplication.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
  
  usage:
  
  make clean
  make  
  ./read_TLee_v20 -f 1 -p 1
  
  ---> README:
  ---> Makefile: comment the line "ROOTSYS=/home/xji/data0/software/root_build", if you have your own "ROOTSYS"
  ---> minuit2 is in the ROOT
*/

int main(int argc, char** argv)
{
  TString roostr = "";
  
  double scaleF_POT = 1;
  int ifile = 1;
  
  for(int i=1; i<argc; i++) {    
    if( strcmp(argv[i],"-p")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>scaleF_POT ) ) { cerr<<" ---> Error scaleF_POT !"<<endl; exit(1); }
    }

    if( strcmp(argv[i],"-f")==0 ) {
      stringstream convert( argv[i+1] );
      if(  !( convert>>ifile ) ) { cerr<<" ---> Error ifile !"<<endl; exit(1); }
    }
  }
  
  cout<<endl<<" ---> check, scaleF_POT "<<scaleF_POT<<", ifile "<<ifile<<endl<<endl;

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
  
  TApplication theApp("theApp",&argc,argv);
  
  TLee *Lee_test = new TLee();
    
  ////////// just do it one time in the whole procedure
  
  Lee_test->scaleF_POT = scaleF_POT;
  Lee_test->Set_config_file_directory(config_Lee::spectra_file, config_Lee::flux_Xs_directory,
				      config_Lee::detector_directory, config_Lee::mc_directory);
  Lee_test->Set_Spectra_MatrixCov();
  Lee_test->Set_POT_implement();
  Lee_test->Set_TransformMatrix();

  ////////// can do any times
  
  Lee_test->flag_syst_flux_Xs    = config_Lee::flag_syst_flux_Xs;
  Lee_test->flag_syst_detector   = config_Lee::flag_syst_detector;
  Lee_test->flag_syst_additional = config_Lee::flag_syst_additional;
  Lee_test->flag_syst_mc_stat    = config_Lee::flag_syst_mc_stat;
  
  Lee_test->scaleF_Lee = config_Lee::Lee_strength_for_output_covariance_matrix;
  Lee_test->Set_Collapse();

  TFile *file_collapsed_covariance_matrix = new TFile("file_collapsed_covariance_matrix.root", "recreate");
  
  TTree *tree_config = new TTree("tree", "configure information");
  int flag_syst_flux_Xs = config_Lee::flag_syst_flux_Xs;
  int flag_syst_detector = config_Lee::flag_syst_detector;
  int flag_syst_additional = config_Lee::flag_syst_additional;
  int flag_syst_mc_stat = config_Lee::flag_syst_mc_stat;
  double user_Lee_strength_for_output_covariance_matrix = config_Lee::Lee_strength_for_output_covariance_matrix;
  double user_scaleF_POT = scaleF_POT;
  tree_config->Branch("flag_syst_flux_Xs", &flag_syst_flux_Xs, "flag_syst_flux_Xs/I" );
  tree_config->Branch("flag_syst_detector", &flag_syst_detector, "flag_syst_detector/I" );
  tree_config->Branch("flag_syst_additional", &flag_syst_additional, "flag_syst_additional/I" );
  tree_config->Branch("flag_syst_mc_stat", &flag_syst_mc_stat, "flag_syst_mc_stat/I" );
  tree_config->Branch("user_Lee_strength_for_output_covariance_matrix", &user_Lee_strength_for_output_covariance_matrix,
		      "user_Lee_strength_for_output_covariance_matrix/D" );
  tree_config->Branch("user_scaleF_POT", &user_scaleF_POT, "user_scaleF_POT/D" );
  tree_config->Fill();
  file_collapsed_covariance_matrix->cd();
  tree_config->Write();

  Lee_test->matrix_absolute_cov_newworld.Write("matrix_absolute_cov_newworld");// (bins, bins)
  Lee_test->matrix_absolute_flux_cov_newworld.Write("matrix_absolute_flux_cov_newworld");
  Lee_test->matrix_absolute_Xs_cov_newworld.Write("matrix_absolute_Xs_cov_newworld");
  Lee_test->matrix_absolute_detector_cov_newworld.Write("matrix_absolute_detector_cov_newworld");
  Lee_test->matrix_absolute_mc_stat_cov_newworld.Write("matrix_absolute_mc_stat_cov_newworld");
  Lee_test->matrix_absolute_additional_cov_newworld.Write("matrix_absolute_additional_cov_newworld");
 
  Lee_test->matrix_pred_newworld.Write("matrix_pred_newworld");// (1, bins)
  Lee_test->matrix_data_newworld.Write("matrix_data_newworld");// (1, bins)  
  file_collapsed_covariance_matrix->Close();
  
  //////////////////////////////////////////////////////////////////////////////////////// Goodness of fit
  
  Lee_test->scaleF_Lee = 0;
  Lee_test->Set_Collapse();
 
  bool flag_both_numuCC            = config_Lee::flag_both_numuCC;// 1
  bool flag_CCpi0_FC_by_numuCC     = config_Lee::flag_CCpi0_FC_by_numuCC;// 2
  bool flag_CCpi0_PC_by_numuCC     = config_Lee::flag_CCpi0_PC_by_numuCC;// 3
  bool flag_NCpi0_by_numuCC        = config_Lee::flag_NCpi0_by_numuCC;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = config_Lee::flag_nueCC_PC_by_numuCC_pi0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = config_Lee::flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = config_Lee::flag_nueCC_LowE_FC_by_all;// 7
  bool flag_nueCC_FC_by_all        = config_Lee::flag_nueCC_FC_by_all;// 8
  
  ///////////////////////// gof
  
  if( flag_nueCC_FC_by_all ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 26, matrix_gof_trans.GetNcols()-26, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 8);
  }
  
  ///////////////////////// gof
  
  if( flag_nueCC_LowE_FC_by_all ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*4 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=26*4 + 11*3; ibin++) matrix_gof_trans(ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 8, matrix_gof_trans.GetNcols()-8, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 7);
  }
  
  ///////////////////////// gof
  
  if( flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, (26-8) + 26*3 + 11*3 );// oldworld, newworld
    for( int ibin=1; ibin<=(26-8) + 26*3 + 11*3; ibin++) matrix_gof_trans(8+ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( (26-8), matrix_gof_trans.GetNcols()-(26-8), matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 6);
  }
  
  ///////////////////////// gof
  
  if( flag_nueCC_PC_by_numuCC_pi0) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*3+11*3 );// oldworld, newworld

    for(int ibin=1; ibin<=26*3+11*3; ibin++) matrix_gof_trans(26+ibin-1, ibin-1) = 1;
      
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );
    
    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;
    
    Lee_test->Exe_Goodness_of_fit( 26, matrix_gof_trans.GetNcols()-26, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 5);
  }
  
  ///////////////////////// gof
  
  if( flag_both_numuCC ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 26*2 );// oldworld, newworld
    for( int ibin=1; ibin<=26*2; ibin++ ) matrix_gof_trans(26*2+ibin-1, ibin-1) = 1;
    
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;

    Lee_test->Exe_Goodness_of_fit( 26*2, 0, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 1);
  }

  ///////////////////////// gof
  
  if( flag_CCpi0_FC_by_numuCC ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 11 + 26*2 );// oldworld, newworld

    for(int ibin=1; ibin<=11; ibin++) matrix_gof_trans(26*4+ibin-1, ibin-1) = 1;
    for(int ibin=1; ibin<=26*2; ibin++) matrix_gof_trans(26*2+ibin-1, 11+ibin-1) = 1;
      
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;
    
    Lee_test->Exe_Goodness_of_fit( 11, 26*2, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 2);
  }
  
  ///////////////////////// gof
  
  if( flag_CCpi0_PC_by_numuCC ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 11 + 26*2 );// oldworld, newworld

    for(int ibin=1; ibin<=11; ibin++) matrix_gof_trans(26*4+11+ibin-1, ibin-1) = 1;
    for(int ibin=1; ibin<=26*2; ibin++) matrix_gof_trans(26*2+ibin-1, 11+ibin-1) = 1;
      
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;
    
    Lee_test->Exe_Goodness_of_fit( 11, 26*2, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 3);
  }
  
  ///////////////////////// gof
  
  if( flag_NCpi0_by_numuCC ) {
    TMatrixD matrix_gof_trans( Lee_test->bins_newworld, 11 + 26*2 );// oldworld, newworld

    for(int ibin=1; ibin<=11; ibin++) matrix_gof_trans(26*4+11*2+ibin-1, ibin-1) = 1;
    for(int ibin=1; ibin<=26*2; ibin++) matrix_gof_trans(26*2+ibin-1, 11+ibin-1) = 1;
      
    TMatrixD matrix_gof_trans_T( matrix_gof_trans.GetNcols(), matrix_gof_trans.GetNrows() );
    matrix_gof_trans_T.Transpose( matrix_gof_trans );

    TMatrixD matrix_gof_pred = Lee_test->matrix_pred_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_data = Lee_test->matrix_data_newworld * matrix_gof_trans;
    TMatrixD matrix_gof_syst = matrix_gof_trans_T * (Lee_test->matrix_absolute_cov_newworld) * matrix_gof_trans;
    
    Lee_test->Exe_Goodness_of_fit( 11, 26*2, matrix_gof_pred, matrix_gof_data, matrix_gof_syst, 4);
  }

  //////////////////////////////////////////////////////////////////////////////////////// LEE strength fitting

  if( config_Lee::flag_Lee_strength_data ) {

    Lee_test->Set_measured_data();

    Lee_test->Minimization_Lee_strength_FullCov(2, 0);

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;

    /////////////////////////////////////////
    
    double gmin = Lee_test->minimization_chi2;
    TGraph *gh_scan = new TGraph();
    double slow = 0;
    double shgh = 3;
    int nscan = 50;
    double val_max_dchi2 = 0;
    double step = (shgh-slow)/nscan;
    for(int idx=1; idx<=nscan; idx++) {
      //cout<<" ---> scan "<<idx<<endl;
      if( idx%(max(1, nscan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./nscan, idx)<<endl;
      double val_s = slow + (idx-1)*step;
      Lee_test->Minimization_Lee_strength_FullCov(val_s, 1);
      double val_chi2 = Lee_test->minimization_chi2;
      gh_scan->SetPoint( gh_scan->GetN(), val_s, val_chi2 - gmin);
      if( val_max_dchi2<val_chi2 - gmin ) val_max_dchi2 = val_chi2 - gmin;
    }
    double val_dchi2at1 = gh_scan->Eval(1);
    
    TCanvas *canv_gh_scan = new TCanvas("canv_gh_scan", "canv_gh_scan", 900, 650);
    canv_gh_scan->SetLeftMargin(0.15); canv_gh_scan->SetRightMargin(0.1);
    canv_gh_scan->SetTopMargin(0.1); canv_gh_scan->SetBottomMargin(0.15);    
    gh_scan->Draw("al");
    gh_scan->GetXaxis()->SetTitle("LEE strength"); gh_scan->GetYaxis()->SetTitle("#Delta#chi^{2}");
    gh_scan->GetXaxis()->SetLabelSize(0.05); gh_scan->GetXaxis()->SetTitleSize(0.05);
    gh_scan->GetYaxis()->SetLabelSize(0.05); gh_scan->GetYaxis()->SetTitleSize(0.05);
    gh_scan->GetXaxis()->CenterTitle(); gh_scan->GetYaxis()->CenterTitle();
    gh_scan->GetXaxis()->SetTitleOffset(1.2);
    gh_scan->GetYaxis()->SetRangeUser(0, val_max_dchi2*1.1);
    // canv_gh_scan->cd(); canv_gh_scan->Update();
    // double x1 = gPad->GetUxmin(); double y1 = gPad->GetUymin();
    // double x2 = gPad->GetUxmax(); double y2 = gPad->GetUymax();    
    TLine *lineA_dchi2at1 = new TLine(1, 0, 1, val_dchi2at1);    
    lineA_dchi2at1->Draw("same");
    lineA_dchi2at1->SetLineWidth(2);
    lineA_dchi2at1->SetLineColor(kBlue);
    lineA_dchi2at1->SetLineStyle(7);
    TLine *lineB_dchi2at1 = new TLine(0, val_dchi2at1, 1, val_dchi2at1);    
    lineB_dchi2at1->Draw("same");
    lineB_dchi2at1->SetLineWidth(2);
    lineB_dchi2at1->SetLineColor(kBlue);
    lineB_dchi2at1->SetLineStyle(7);
    auto *tt_text_data = new TLatex( 0.2, val_dchi2at1*1.1, Form("#Delta#chi^{2} = %3.2f", val_dchi2at1) );
    tt_text_data->SetTextAlign(11); tt_text_data->SetTextSize(0.05); tt_text_data->SetTextAngle(0);
    tt_text_data->SetTextFont(42);  tt_text_data->Draw(); tt_text_data->SetTextColor(kBlue);

    canv_gh_scan->SaveAs("canv_gh_scan.png");
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////

  if( 0 ) {
    double chi2_null_null8sm_true8sm  = 0;
    double chi2_gmin_null8sm_true8sm  = 0;
    double chi2_null_null8Lee_true8Lee = 0;
    double chi2_gmin_null8Lee_true8Lee = 0;

    TFile *file_out = new TFile(TString::Format("file_out_%03d.root", ifile), "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("chi2_null_null8sm_true8sm", &chi2_null_null8sm_true8sm, "chi2_null_null8sm_true8sm/D" );
    tree->Branch("chi2_gmin_null8sm_true8sm", &chi2_gmin_null8sm_true8sm, "chi2_gmin_null8sm_true8sm/D" );
    tree->Branch("chi2_null_null8Lee_true8Lee", &chi2_null_null8Lee_true8Lee, "chi2_null_null8Lee_true8Lee/D" );
    tree->Branch("chi2_gmin_null8Lee_true8Lee", &chi2_gmin_null8Lee_true8Lee, "chi2_gmin_null8Lee_true8Lee/D" );

    int N_toy = 100;
        
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();
    
    Lee_test->Set_Variations(N_toy);
    
    for(int itoy=1; itoy<=N_toy; itoy++) {
      if( itoy%max(N_toy/10,1)==0 ) {
	cout<<TString::Format(" ---> processing toy ( total cov ): %4.2f, %6d", itoy*1./N_toy, itoy)<<endl;
      }

      int status_fit = 0;
      
      Lee_test->Set_toy_Variation( itoy );

      ///////////////////////////////////
      
      Lee_test->Minimization_Lee_strength_FullCov(1, 1);
      chi2_null_null8Lee_true8Lee = Lee_test->minimization_chi2;

      Lee_test->Minimization_Lee_strength_FullCov(1, 0);
      chi2_gmin_null8Lee_true8Lee = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;

      if( status_fit!=0 ) continue;
      tree->Fill();
    }

    file_out->cd();
    tree->Write();
    file_out->Close();
    
  }

  ////////////////////////////////////////////////////////////////////////////////////////

  if( 0 ) {    
    Lee_test->Set_measured_data();

    Lee_test->Minimization_Lee_strength_FullCov(1, 1);
    double val_chi2_Lee = Lee_test->minimization_chi2;

    Lee_test->Minimization_Lee_strength_FullCov(0, 1);
    double val_chi2_sm = Lee_test->minimization_chi2;

    double val_dchi2 = val_chi2_Lee - val_chi2_sm;
    
    cout<<endl<<TString::Format(" ---> dchi2 = Lee - sm: %4.2f", val_dchi2)<<endl<<endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////// test, Asimov

  if( 0 ) {
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();
  
    Lee_test->Set_toy_Asimov();

    Lee_test->Minimization_Lee_strength_FullCov(2, 0);

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;
  }

  //////////////////////////////////////////////////////////////////////////////////////// test, variation

  if( 0 ) {
    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();

    Lee_test->Set_Variations( 10 );
    Lee_test->Set_toy_Variation( 4 );
    
    Lee_test->Minimization_Lee_strength_FullCov(2, 0);

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;
  }
  
  //////////////////////////////////////////////////////////////////////////////////////// Feldman-Cousins approach

  if( 0 ) {
    
    Lee_test->Set_measured_data();    
    TMatrixD matrix_data_input_fc = Lee_test->matrix_data_newworld;
    
    // Lee_test->Minimization_Lee_strength_FullCov(1, 0);
    
    // cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
    // 				Lee_test->minimization_chi2,
    // 				Lee_test->minimization_Lee_strength_val,
    // 				Lee_test->minimization_Lee_strength_err
    // 				)<<endl<<endl;
    
    /////////////// Feldman-Cousins
    
    double Lee_true_low = 0;
    double Lee_true_hgh = 5;
    double Lee_step     = 0.02;
    
    //int num_toy = 1;    
    //Lee_test->Exe_Feldman_Cousins(Lee_true_low, Lee_true_hgh, Lee_step, num_toy, ifile);

    /////////////// Asimov
    
    //Lee_test->Exe_Fledman_Cousins_Asimov(Lee_true_low, Lee_true_hgh, Lee_step);

    /////////////// measured data
    
    //Lee_test->Exe_Fiedman_Cousins_Data( matrix_data_input_fc, Lee_true_low, Lee_true_hgh, Lee_step );
    
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  cout<<endl<<" Entrer Ctrl+c to end the program"<<endl<<endl;
  
  //theApp.Run();

  return 0;
}
