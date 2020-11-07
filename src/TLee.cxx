#include "WCPLEEANA/TLee.h"

#include "draw.icc"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace DataBase {
  double x1[101]={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
		  11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		  21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		  31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
		  41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
		  51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
		  61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
		  71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
		  81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
		  91, 92, 93, 94, 95, 96, 97, 98, 99, 100};
  
  double yl[101]={0, 0, 0, 0.856632, 1.70317, 2.51005, 3.32075, 4.14046, 4.9693, 5.80646, 6.65117,
		  7.5025, 8.35978, 9.22237, 10.0898, 10.9615, 11.8372, 12.7165, 13.5992, 14.4849, 15.3734,
		  16.2646, 17.1583, 18.0543, 18.9524, 19.8526, 20.7547, 21.6586, 22.5642, 23.4715, 24.3803,
		  25.2906, 26.2023, 27.1153, 28.0297, 28.9452, 29.8619, 30.7797, 31.6987, 32.6187, 33.5396,
		  34.4616, 35.3845, 36.3083, 37.2329, 38.1584, 39.0847, 40.0118, 40.9396, 41.8682, 42.7975,
		  43.7275, 44.6581, 45.5895, 46.5215, 47.454, 48.3873, 49.321, 50.2554, 51.1903, 52.1257,
		  53.0617, 53.9982, 54.9352, 55.8727, 56.8107, 57.7491, 58.6881, 59.6274, 60.5673, 61.5075,
		  62.4482, 63.3892, 64.3307, 65.2725, 66.2148, 67.1575, 68.1005, 69.0438, 69.9876, 70.9317,
		  71.8761, 72.8209, 73.766, 74.7114, 75.6572, 76.6033, 77.5497, 78.4964, 79.4434, 80.3907,
		  81.3383, 82.2862, 83.2342, 84.1827, 85.1314, 86.0804, 87.0296, 87.9791, 88.9288, 89.8788};

  double yh[101]={1.1478, 2.35971, 3.51917, 4.72422, 5.98186, 7.21064, 8.41858, 9.61053, 10.7896, 11.9582, 13.1179,
		  14.27, 15.4155, 16.5552, 17.6898, 18.8197, 19.9454, 21.0673, 22.1858, 23.3011, 24.4133,
		  25.5229, 26.6299, 27.7346, 28.837, 29.9374, 31.0358, 32.1322, 33.2271, 34.3201, 35.4117,
		  36.5017, 37.5904, 38.6776, 39.7635, 40.8483, 41.9318, 43.0141, 44.0955, 45.1757, 46.2549,
		  47.3331, 48.4104, 49.4868, 50.5623, 51.637, 52.7108, 53.7839, 54.8561, 55.9277, 56.9985,
		  58.0686, 59.1381, 60.2068, 61.275, 62.3425, 63.4094, 64.4757, 65.5415, 66.6066, 67.6713,
		  68.7354, 69.7989, 70.862, 71.9246, 72.9866, 74.0483, 75.1094, 76.1701, 77.2304, 78.2902,
		  79.3496, 80.4085, 81.4672, 82.5253, 83.5831, 84.6406, 85.6976, 86.7542, 87.8105, 88.8665,
		  89.9221, 90.9774, 92.0323, 93.0869, 94.1411, 95.1951, 96.2488, 97.3021, 98.3552, 99.4079,
		  100.46, 101.513, 102.564, 103.616, 104.667, 105.718, 106.769, 107.82, 108.87, 109.92};
}


///////////////////////////////////////////////////////// ccc
void TLee::Exe_Fiedman_Cousins_Data(TMatrixD matrix_fakedata, double Lee_true_low, double Lee_true_hgh, double step)
{
  cout<<endl<<" ---> Exe_Feldman_Cousins_Data"<<endl;
  
  Set_fakedata( matrix_fakedata );
  
  //////////////////
  
  Minimization_Lee_strength_FullCov(1, 0);
  
  cout<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
			minimization_chi2,
			minimization_Lee_strength_val,
			minimization_Lee_strength_err
			)<<endl;
  
  //////////////////

  double Lee_bestFit_data = minimization_Lee_strength_val;
  double Lee_bestFit_err = minimization_Lee_strength_err;
  double chi2_gmin_data = minimization_chi2;
  vector<int>Lee_scan100_data;
  vector<double>chi2_null_scan_data;
  
  TTree *tree_data = new TTree("tree_data", "Feldman-Cousins");

  tree_data->Branch( "Lee_bestFit_data", &Lee_bestFit_data, "Lee_bestFit_data/D" );
  tree_data->Branch( "Lee_bestFit_err", &Lee_bestFit_err, "Lee_bestFit_err/D" );
  tree_data->Branch( "chi2_gmin_data", &chi2_gmin_data, "chi2_gmin_data/D");
  tree_data->Branch( "Lee_scan100_data", &Lee_scan100_data );
  tree_data->Branch( "chi2_null_scan_data", &chi2_null_scan_data );

  int num_scan = int(((Lee_true_hgh-Lee_true_low)/step)+0.5) + 1;
  
  for(int idx=1; idx<=num_scan; idx++ ) {
    if( idx%(max(1, num_scan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./num_scan, idx)<<endl;

    double Lee_strength = Lee_true_low + (idx-1)*step;
    int Lee_strength_scaled100 = (int)(Lee_strength*100 + 0.5);
    
    Minimization_Lee_strength_FullCov(Lee_strength, 1);
    chi2_null_scan_data.push_back( minimization_chi2 );
    Lee_scan100_data.push_back( Lee_strength_scaled100 );
  }
  cout<<endl;
  
  tree_data->Fill();
  
  TFile *file_data = new TFile("file_data.root", "recreate");
  file_data->cd();
  tree_data->Write();
  file_data->Close();      
}
  
void TLee::Exe_Fledman_Cousins_Asimov(double Lee_true_low, double Lee_true_hgh, double step)
{
  cout<<endl<<" ---> Exe_Feldman_Cousins_Asimov"<<endl<<endl;
  
  ///////////////////////

  int Lee_strength_scaled100  = 0;  
  vector<double>Lee_scan100;
  vector<double>chi2_null_toy;
  
  TTree *tree_Asimov = new TTree("tree_Asimov", "Feldman-Cousins");

  tree_Asimov->Branch( "Lee_strength_scaled100", &Lee_strength_scaled100, "Lee_strength_scaled100/I" );
  tree_Asimov->Branch( "Lee_scan100", &Lee_scan100 );
  tree_Asimov->Branch( "chi2_null_toy", &chi2_null_toy );

  int num_scan = int(((Lee_true_hgh-Lee_true_low)/step)+0.5) + 1;
  
  for(int idx=1; idx<=num_scan; idx++ ) {
    if( idx%(max(1, num_scan/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./num_scan, idx)<<endl;

    double Lee_strength = Lee_true_low + (idx-1)*step;
    Lee_strength_scaled100 = (int)(Lee_strength*100 + 0.5);
    
    Lee_scan100.clear();
    chi2_null_toy.clear();

    scaleF_Lee = Lee_strength;
    Set_Collapse();
    Set_toy_Asimov();

    for(int jdx=1; jdx<=num_scan; jdx++) {
      double val_Lee_scan = Lee_true_low + (jdx-1)*step;
      double val_Lee_scan100 = (int)(val_Lee_scan*100 + 0.5);
      
      Minimization_Lee_strength_FullCov(val_Lee_scan, 1);
      double val_chi2_null = minimization_chi2;

      Lee_scan100.push_back( val_Lee_scan100 );
      chi2_null_toy.push_back( val_chi2_null );
    }

    tree_Asimov->Fill();
    
  }// idx

  TFile *file_Asimov = new TFile("file_Asimov.root", "recreate");
  file_Asimov->cd();
  tree_Asimov->Write();
  file_Asimov->Close();
  
}

void TLee::Exe_Feldman_Cousins(double Lee_true_low, double Lee_true_hgh, double step, int num_toy, int ifile)
{
  cout<<endl<<" ---> Exe_Feldman_Cousins"<<endl<<endl;

  ///////////////////////

  int Lee_strength_scaled100  = 0;
  vector<double>chi2_null_toy;
  vector<double>chi2_gmin_toy;
  vector<double>LeeF_gmin_toy;
    
  TTree *tree = new TTree("tree", "Feldman-Cousins");

  tree->Branch( "Lee_strength_scaled100", &Lee_strength_scaled100, "Lee_strength_scaled100/I" );
  tree->Branch( "chi2_null_toy", &chi2_null_toy );
  tree->Branch( "chi2_gmin_toy", &chi2_gmin_toy );
  tree->Branch( "LeeF_gmin_toy", &LeeF_gmin_toy );

  int num_idx = int(((Lee_true_hgh-Lee_true_low)/step)+0.5) + 1;
  
  for(int idx=1; idx<=num_idx; idx++ ) {
    if( idx%(max(1, num_idx/10))==0 ) cout<<Form(" ---> scan %4.2f, %3d", idx*1./num_idx, idx)<<endl;

    double Lee_strength = Lee_true_low + (idx-1)*step;

    /////
    scaleF_Lee = Lee_strength;
    Set_Collapse();
    
    Set_Variations( num_toy );

    Lee_strength_scaled100 = (int)(Lee_strength*100 + 0.5);
    
    chi2_null_toy.clear();
    chi2_gmin_toy.clear();    
    LeeF_gmin_toy.clear();
      
    for(int itoy=1; itoy<=num_toy; itoy++) {
      Set_toy_Variation( itoy );

      Minimization_Lee_strength_FullCov(Lee_strength, 1);
      double val_chi2_null = minimization_chi2;
      
      Minimization_Lee_strength_FullCov(Lee_strength, 0);
      double val_chi2_gmin = minimization_chi2;
      if( minimization_status!=0 ) continue;

      chi2_null_toy.push_back( val_chi2_null );
      chi2_gmin_toy.push_back( val_chi2_gmin );
      LeeF_gmin_toy.push_back( minimization_Lee_strength_val );
      
    }// itoy

    tree->Fill();
    
  }// idx
  
  ///////////////////////
  
  TFile *file_FC = new TFile(Form("file_FC_%06d.root", ifile), "recreate");
  file_FC->cd();
  tree->Write();
  file_FC->Close();
  
}

///////////////////////////////////////////////////////// ccc

void TLee::Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed)
{
  TString roostr = "";
  
  ROOT::Minuit2::Minuit2Minimizer min_Lee( ROOT::Minuit2::kMigrad );
  min_Lee.SetPrintLevel(0);
  min_Lee.SetStrategy(1); //0- cursory, 1- default, 2- thorough yet no more successful
  min_Lee.SetMaxFunctionCalls(500000);
  min_Lee.SetMaxIterations(500000);
  min_Lee.SetTolerance(1e-6); // tolerance*2e-3 = edm precision
  min_Lee.SetPrecision(1e-18); //precision in the target function

  /// set fitting parameters
  ROOT::Math::Functor Chi2Functor_Lee( [&](const double *par) {// FCN
      TString roostr = "";
      double chi2 = 0;
      double Lee_strength = par[0];

      /////////      
      TMatrixD matrix_meas(1, bins_newworld);
      for(int ibin=0; ibin<matrix_meas.GetNcols(); ibin++) {
	matrix_meas(0, ibin) = map_fake_data[ibin];	
      }

      /////////
  
      scaleF_Lee = Lee_strength;
      Set_Collapse();
      
      TMatrixD matrix_pred = matrix_pred_newworld;

      /////////
      TMatrixD matrix_cov_syst = matrix_absolute_cov_newworld;
      
      for(int ibin=0; ibin<matrix_cov_syst.GetNrows(); ibin++) {
	double val_stat_cov = 0;	
	double val_meas = matrix_meas(0, ibin);
	double val_pred = matrix_pred(0, ibin);
	
	if( val_meas==0 ) val_stat_cov = val_pred/2;
	else val_stat_cov = 3./( 1./val_meas + 2./val_pred );	
	if( val_meas==0 && val_pred==0 ) val_stat_cov = 1e-6;
	matrix_cov_syst(ibin, ibin) += val_stat_cov;
      }

      TMatrixD matrix_cov_total = matrix_cov_syst;
      TMatrixD matrix_cov_total_inv = matrix_cov_total;
      matrix_cov_total_inv.Invert();
      
      ////////
      TMatrixD matrix_delta = matrix_pred - matrix_meas;
      TMatrixD matrix_delta_T( matrix_delta.GetNcols(), matrix_delta.GetNrows() );
      matrix_delta_T.Transpose( matrix_delta );
      
      TMatrixD matrix_chi2 = matrix_delta * matrix_cov_total_inv *matrix_delta_T;
      chi2 = matrix_chi2(0,0);
      
      /////////
                  
      return chi2;
      
    },// end of FCN
    1 // number of fitting parameters
    );
  
  min_Lee.SetFunction(Chi2Functor_Lee);
  
  min_Lee.SetVariable( 0, "Lee_strength", Lee_initial_value, 1e-2);
  //min_Lee.SetVariableLowerLimit(0, 0);
  min_Lee.SetLowerLimitedVariable(0, "Lee_strength", Lee_initial_value, 1e-2, 0);
  if( flag_fixed ) {
    min_Lee.SetFixedVariable( 0, "Lee_strength", Lee_initial_value );
  }
  
  /// do the minimization
  min_Lee.Minimize();
  int status_Lee = min_Lee.Status();
  const double *par_Lee = min_Lee.X();
  const double *par_Lee_err = min_Lee.Errors();

  if( status_Lee!=0 ) {
    cerr<<endl<<" -----------> Lee strength fitting failed "<<endl<<endl;
    minimization_status = status_Lee;
  }

  minimization_status = status_Lee;
  minimization_chi2 = min_Lee.MinValue();
  minimization_Lee_strength_val = par_Lee[0];
  minimization_Lee_strength_err = par_Lee_err[0];

  /// MinosError
  // {
  //   min_Lee.SetErrorDef(1);
  //   double minosError_low = 0;
  //   double minosError_hgh = 0;
  //   min_Lee.GetMinosError(0, minosError_low, minosError_hgh);
  //   cout<<TString::Format(" ---> Best fit of Lee strength: %5.2f, (1 sigma) from MinosError: %5.2f %5.2f",
  //  			  minimization_Lee_strength_val, minosError_low, minosError_hgh)<<endl;
  // }
  // {
  //   min_Lee.SetErrorDef(4);
  //   double minosError_low = 0;
  //   double minosError_hgh = 0;
  //   min_Lee.GetMinosError(0, minosError_low, minosError_hgh);
  //   cout<<TString::Format(" ---> Best fit of Lee strength: %5.2f, (2 sigma) from MinosError: %5.2f %5.2f",
  // 			  minimization_Lee_strength_val, minosError_low, minosError_hgh)<<endl;
  // }
  // {
  //   min_Lee.SetErrorDef(9);
  //   double minosError_low = 0;
  //   double minosError_hgh = 0;
  //   min_Lee.GetMinosError(0, minosError_low, minosError_hgh);
  //   cout<<TString::Format(" ---> Best fit of Lee strength: %5.2f, (3 sigma) from MinosError: %5.2f %5.2f",
  // 			  minimization_Lee_strength_val, minosError_low, minosError_hgh)<<endl;
  // }
  
}  

///////////////////////////////////////////////////////// ccc

void TLee::Set_toy_Asimov()
{
  for(int ibin=0; ibin<bins_newworld; ibin++) map_fake_data[ibin] = matrix_pred_newworld(0, ibin);
}

void TLee::Set_toy_Variation(int itoy)
{
  for(int ibin=0; ibin<bins_newworld; ibin++) map_fake_data[ibin] = map_toy_variation[itoy][ibin];
}

void TLee::Set_measured_data()
{
  for(int ibin=0; ibin<bins_newworld; ibin++) map_fake_data[ibin] = matrix_data_newworld(0, ibin);
}
  
void TLee::Set_fakedata(TMatrixD matrix_fakedata)
{
  int cols = matrix_fakedata.GetNcols();
  for(int ibin=0; ibin<cols; ibin++) map_fake_data[ibin] = matrix_fakedata(0, ibin);    
}

///////////////////////////////////////////////////////// ccc

void TLee::Set_Variations(int num_toy)
{
  map_toy_variation.clear();

  //////////////////////////

  TMatrixDSym DSmatrix_cov(bins_newworld);
  for(int ibin=0; ibin<bins_newworld; ibin++) {
    for(int jbin=0; jbin<bins_newworld; jbin++) {
      DSmatrix_cov(ibin, jbin) = matrix_absolute_cov_newworld(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen( DSmatrix_cov );
  TMatrixD matrix_eigenvector = DSmatrix_eigen.GetEigenVectors();
  TVectorD matrix_eigenvalue = DSmatrix_eigen.GetEigenValues();

  for(int itoy=1; itoy<=num_toy; itoy++) {    
    TMatrixD matrix_element(bins_newworld, 1);    
    for(int ibin=0; ibin<bins_newworld; ibin++) {
      if( matrix_eigenvalue(ibin)>=0 ) {
	matrix_element(ibin,0) = rand->Gaus( 0, sqrt( matrix_eigenvalue(ibin) ) );
      }
      else {
	matrix_element(ibin,0) = 0;
      }      
    }
    TMatrixD matrix_variation = matrix_eigenvector * matrix_element;
    for(int ibin=0; ibin<bins_newworld; ibin++) {
      double val_with_syst = matrix_variation(ibin,0) + map_pred_spectrum_newworld_bin[ibin];// key point
      if( val_with_syst<0 ) val_with_syst = 0;
      map_toy_variation[itoy][ibin] = rand->PoissonD( val_with_syst );
    }
  }
    
}

///////////////////////////////////////////////////////// ccc

void TLee::Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst, int index)
{
  TString roostr = "";

  cout<<Form(" ---> Goodness of fit, %2d", index)<<endl;

  int color_no = kRed;
  int color_wi = kBlue;

  /////////////////////////////////////////////////////////////////////////////////////////////

  {
    double val_cov = 0;
    for(int ibin=0; ibin<matrix_syst.GetNrows(); ibin++) {
      val_cov += matrix_syst(ibin, ibin);
    }
    if(val_cov==0) {
      for(int ibin=0; ibin<matrix_syst.GetNrows(); ibin++) {
	matrix_syst(ibin, ibin) = 1e-6;
      } 
    }    
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////
  if( num_X==0 ) {
    TMatrixD matrix_delta = matrix_pred - matrix_data;
    TMatrixD matrix_delta_T( matrix_delta.GetNcols(), matrix_delta.GetNrows() );
    matrix_delta_T.Transpose( matrix_delta );

    int rows = matrix_syst.GetNrows();
    TMatrixD matrix_stat(rows, rows);
    for(int i=0; i<rows; i++) {
      double val_pred = matrix_pred(0, i);
      double val_data = matrix_data(0, i);      
      matrix_stat(i,i) = val_pred;      
      if( val_data==1 ) {
	if( val_pred<0.461 ) {// DocDB-32520, when the prediction is sufficiently low
	  double numerator = pow(val_pred-val_data, 2);
	  double denominator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
	  matrix_stat(i,i) = numerator/denominator;
	}
      }      
      if( (val_pred==val_data) && (val_pred==0) ) matrix_stat(i,i) = 1e-6;
    }
    
    TMatrixD matrix_total_cov(rows, rows);
    matrix_total_cov = matrix_stat + matrix_syst;

    TMatrixD matrix_total_cov_inv = matrix_total_cov;
    matrix_total_cov_inv.Invert();

    ///
    TMatrixD matrix_chi2 = matrix_delta * matrix_total_cov_inv *matrix_delta_T;
    double val_chi2 = matrix_chi2(0,0);
    double p_value = TMath::Prob( val_chi2, rows );
    
    cout<<endl<<TString::Format(" ---> GOF: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.8f",
				val_chi2, rows, val_chi2/rows, p_value
				)<<endl<<endl;

    ///////////////////////////

    roostr = TString::Format("h1_pred_%02d", index);
    TH1D *h1_pred = new TH1D(roostr, "", rows, 0, rows);
    for(int ibin=0; ibin<rows; ibin++) {
      h1_pred->SetBinContent(ibin+1, matrix_pred(0, ibin));
      h1_pred->SetBinError( ibin+1, sqrt(matrix_syst(ibin,ibin)) );
    }

    ////
    roostr = TString::Format("h1_ratio_basic_%02d", index);
    TH1D *h1_ratio_basic = new TH1D(roostr, "", rows, 0, rows);

    ////
    TGraphAsymmErrors *gh_ratio = new TGraphAsymmErrors();
    TGraphAsymmErrors *gh_data = new TGraphAsymmErrors();
    double val_max_data = 0;
    
    for(int ibin=0; ibin<rows; ibin++) {
      double val_data = matrix_data(0, ibin);
      double val_pred = matrix_pred(0, ibin);
      double val_ratio= val_data/val_pred;
      if( val_ratio!=val_ratio || val_ratio==1./0 ) val_ratio = 0;

      double val_data_low = 0;
      double val_data_hgh = 0;
      int idx_data = (int)(val_data+0.5);    
      if( idx_data>100 ) {
	val_data_low = val_data - sqrt(val_data);
	val_data_hgh = val_data + sqrt(val_data);
      }
      else {
	val_data_low = DataBase::yl[ idx_data ];
	val_data_hgh = DataBase::yh[ idx_data ];
      }
      double val_ratio_low = val_ratio - val_data_low/val_pred;
      double val_ratio_hgh = val_data_hgh/val_pred - val_ratio;
       
      int n_point = gh_ratio->GetN();
      double val_x = h1_ratio_basic->GetBinCenter(ibin+1);
      double val_halfw = h1_ratio_basic->GetBinWidth(ibin+1)/2;
      gh_ratio->SetPoint( n_point, val_x, val_ratio );
      gh_ratio->SetPointError( n_point, val_halfw, val_halfw, val_ratio_low, val_ratio_hgh );
      
      gh_data->SetPoint( n_point, val_x, val_data );
      gh_data->SetPointError( n_point, val_halfw, val_halfw, val_data-val_data_low, val_data_hgh-val_data );

      if( val_max_data<val_data_hgh ) val_max_data = val_data_hgh;
    }
    
    double val_max_pred = h1_pred->GetBinContent( h1_pred->GetMaximumBin() ) + h1_pred->GetBinError( h1_pred->GetMaximumBin() );
    
    /////////////////////////// plotting

    roostr = TString::Format("canv_spectra_cmp_%02d", index);
    TCanvas *canv_spectra_cmp = new TCanvas(roostr, roostr, 1000, 950);  
    
    /////////////
    TPad *pad_top = new TPad("pad_top", "pad_top", 0, 0.45, 1, 1);
    func_canv_margin(pad_top, 0.15, 0.1, 0.1, 0.05);
    pad_top->Draw();
    pad_top->cd();
    
    h1_pred->Draw("e2");
    h1_pred->SetMinimum(0);
    if( val_max_data > val_max_pred ) h1_pred->SetMaximum( val_max_data*1.1 );

    TH1D *h1_pred_clone = (TH1D*)h1_pred->Clone("h1_pred_clone");
    
    h1_pred->Draw("same e2");
    h1_pred->SetMarkerStyle(1); h1_pred->SetMarkerColor(kRed);
    h1_pred->SetLineColor(kRed); h1_pred->SetLineWidth(2);
    h1_pred->SetFillColor(kRed); h1_pred->SetFillStyle(3005);
    h1_pred->SetTitle("");
    func_title_size(h1_pred, 0.065, 0.065, 0.065, 0.065);
    func_xy_title(h1_pred, "Bin index", "Entries");
    h1_pred->GetXaxis()->SetNdivisions(506); h1_pred->GetYaxis()->SetNdivisions(506);
    h1_pred->GetXaxis()->CenterTitle(); h1_pred->GetYaxis()->CenterTitle();
    h1_pred->GetXaxis()->SetLabelColor(10); h1_pred->GetXaxis()->SetTitleColor(10);
    h1_pred->GetYaxis()->SetTitleOffset(1.1);
    
    h1_pred_clone->Draw("same hist");
    h1_pred_clone->SetLineColor(kRed);
    h1_pred_clone->SetLineWidth(3);

    gh_data->Draw("same pe");
    gh_data->SetMarkerStyle(20); gh_data->SetMarkerSize(1.2); gh_data->SetMarkerColor(kBlue);
    gh_data->SetLineColor(kBlue); gh_data->SetLineWidth(3);
 
    h1_pred->Draw("same axis");
   
    TLegend *lg_top = new TLegend(0.2, 0.7, 0.45, 0.85);
    lg_top->AddEntry(gh_data, "Data", "lep");
    lg_top->AddEntry(h1_pred, "Prediction", "lf");
    lg_top->Draw();
    lg_top->SetBorderSize(0);
    lg_top->SetFillStyle(0);
    lg_top->SetTextSize(0.065);
    
    pad_top->cd();
    pad_top->Update();
    double x1 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.07;
    double y1 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.6;  
    double x2 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.35;
    double y2 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.7;
    TPaveText *pt = new TPaveText( x1,y1,x2,y2,"l");
    pt->SetTextSize(0.065);
    pt->SetTextFont(42);
    pt->SetTextAlign(11);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d", val_chi2,rows));
    //pt->AddText(TString::Format("P-value: %5.3f", p_value));
    pt->Draw();
    
    /////////////
    canv_spectra_cmp->cd();
  
    TPad *pad_bot = new TPad("pad_bot", "pad_bot", 0, 0, 1, 0.45);
    func_canv_margin(pad_bot, 0.15, 0.1, 0.05, 0.3);
    pad_bot->Draw();
    pad_bot->cd();
  
    h1_ratio_basic->Draw();
    h1_ratio_basic->SetMinimum(0); h1_ratio_basic->SetMaximum(2);
    //h1_ratio_basic->SetLineColor(kRed);
    //h1_ratio_basic->SetLineStyle(7);
    func_title_size(h1_ratio_basic, 0.078, 0.078, 0.078, 0.078);
    func_xy_title(h1_ratio_basic, "Bin index", "Data / Pred");
    h1_ratio_basic->GetXaxis()->SetTickLength(0.05);
    h1_ratio_basic->GetXaxis()->CenterTitle(); h1_ratio_basic->GetYaxis()->CenterTitle(); 
    h1_ratio_basic->GetYaxis()->SetTitleOffset(0.92);
    h1_ratio_basic->GetYaxis()->SetNdivisions(509);
    h1_ratio_basic->Draw("same axis");

    TH1D *h1_pred_rel_error = (TH1D*)h1_pred->Clone("h1_pred_rel_error");
    h1_pred_rel_error->Reset();
    for(int ibin=0; ibin<rows; ibin++) {
      double val_err = h1_pred->GetBinError(ibin+1);
      double val_cv = h1_pred->GetBinContent(ibin+1);
      double rel_err = val_err/val_cv;
      if( val_cv==0 ) rel_err = 0;
      h1_pred_rel_error->SetBinContent(ibin+1, 1);
      h1_pred_rel_error->SetBinError(ibin+1, rel_err);
    }
    h1_pred_rel_error->Draw("same e2");

    TF1 *f1_c1 = new TF1("f1_c1", "1", 0, 1e5);
    f1_c1->Draw("same");
    f1_c1->SetLineStyle(7);
    f1_c1->SetLineWidth(2);
    f1_c1->SetLineColor(kBlack);

    gh_ratio->Draw("same pe");
    gh_ratio->SetMarkerStyle(20);
    gh_ratio->SetMarkerSize(1);
    gh_ratio->SetMarkerColor(kBlue);
    gh_ratio->SetLineColor(kBlue);
    
    h1_ratio_basic->Draw("same axis");
 
    roostr = TString::Format("canv_spectra_cmp_%02d.png", index);
    canv_spectra_cmp->SaveAs(roostr);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  else {    
    TMatrixD matrix_pred_Y(num_Y, 1);
    TMatrixD matrix_data_Y(num_Y, 1);
  
    TMatrixD matrix_pred_X(num_X, 1);
    TMatrixD matrix_data_X(num_X, 1);
  
    for(int ibin=1; ibin<=num_Y; ibin++) {
      matrix_pred_Y(ibin-1, 0) = matrix_pred(0, ibin-1);
      matrix_data_Y(ibin-1, 0) = matrix_data(0, ibin-1);
    }
    
    for(int ibin=1; ibin<=num_X; ibin++) {
      matrix_pred_X(ibin-1, 0) = matrix_pred(0, num_Y+ibin-1);
      matrix_data_X(ibin-1, 0) = matrix_data(0, num_Y+ibin-1);
    }

    TMatrixD matrix_cov_stat(num_Y+num_X, num_Y+num_X);

    TMatrixD matrix_cov_total(num_Y+num_X, num_Y+num_X);
    matrix_cov_total = matrix_cov_stat + matrix_syst;
    for(int idx=1; idx<=num_Y+num_X; idx++) {
      if( matrix_cov_total(idx-1, idx-1)==0 ) matrix_cov_total(idx-1, idx-1) = 1e-6;// case inverse
    }   

    TMatrixD matrix_YY(num_Y, num_Y);
    for(int ibin=1; ibin<=num_Y; ibin++) {
      for(int jbin=1; jbin<=num_Y; jbin++) {
	matrix_YY(ibin-1, jbin-1) = matrix_cov_total(ibin-1, jbin-1);
      }
    }
  
    TMatrixD matrix_XX(num_X, num_X);
    for(int ibin=1; ibin<=num_X; ibin++) {
      for(int jbin=1; jbin<=num_X; jbin++) {
	matrix_XX(ibin-1, jbin-1) = matrix_cov_total(num_Y+ibin-1, num_Y+jbin-1);
	if(ibin==jbin) matrix_XX(ibin-1, jbin-1) += matrix_pred_X(ibin-1, 0);// Pearson's term for statistics
      }
    }
    TMatrixD matrix_XX_inv = matrix_XX;
    matrix_XX_inv.Invert();

    TMatrixD matrix_YX(num_Y, num_X);
    TMatrixD matrix_XY(num_X, num_Y);
    for(int ibin=1; ibin<=num_Y; ibin++) {
      for(int jbin=1; jbin<=num_X; jbin++) {
	matrix_YX(ibin-1, jbin-1) = matrix_cov_total( ibin-1, num_Y+jbin-1 );
      }
    }
    matrix_XY.Transpose( matrix_YX );
      
    /////////////////////////////

    TMatrixD matrix_Y_under_X = matrix_pred_Y + matrix_YX * matrix_XX_inv * (matrix_data_X - matrix_pred_X);
    
    TMatrixD matrix_YY_under_XX = matrix_YY - matrix_YX * matrix_XX_inv * matrix_XY;
  
    ///////////////////////////// goodness of fit, Pearson's format

    TMatrixD matrix_goodness_cov_total_noConstraint(num_Y, num_Y);
    for( int i=0; i<num_Y; i++ ) {
      double val_pred = matrix_pred_Y(i, 0);
      double val_data = matrix_data_Y(i, 0);    
      matrix_goodness_cov_total_noConstraint(i,i) = val_pred;    
      if( val_data==1 ) {
	if( val_pred<0.461 ) {// DocDB-32520, when the prediction is sufficiently low
	  double numerator = pow(val_pred-val_data, 2);
	  double denominator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
	  matrix_goodness_cov_total_noConstraint(i,i) = numerator/denominator;
	}
      }

      if( (val_pred==val_data) && (val_pred==0) ) matrix_goodness_cov_total_noConstraint(i,i) = 1e-6;
    }
  
    matrix_goodness_cov_total_noConstraint = matrix_goodness_cov_total_noConstraint + matrix_YY;
  
    TMatrixD matrix_delta_noConstraint(1, num_Y);
    for(int ibin=1; ibin<=num_Y; ibin++) {
      double val_pred = matrix_pred_Y(ibin-1, 0);
      double val_data = matrix_data_Y(ibin-1, 0);
      matrix_delta_noConstraint(0, ibin-1) = val_pred - val_data;
    }
    TMatrixD matrix_delta_noConstraint_T(num_Y, 1);
    matrix_delta_noConstraint_T.Transpose(matrix_delta_noConstraint);
    TMatrixD matrix_cov_noConstraint_inv = matrix_goodness_cov_total_noConstraint;
    matrix_cov_noConstraint_inv.Invert();
    TMatrixD matrix_chi2_noConstraint = matrix_delta_noConstraint * matrix_cov_noConstraint_inv * matrix_delta_noConstraint_T;
    double val_chi2_noConstraint = matrix_chi2_noConstraint(0,0);
    double p_value_noConstraint = TMath::Prob( val_chi2_noConstraint, num_Y );  
    cout<<endl<<TString::Format(" ---> GOF noConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.8f",
				val_chi2_noConstraint, num_Y, val_chi2_noConstraint/num_Y, p_value_noConstraint
				)<<endl;

    ////////
    ////////
  
    TMatrixD matrix_goodness_cov_total_wiConstraint(num_Y, num_Y);
    for( int i=0; i<num_Y; i++ ) {
      double val_pred = matrix_Y_under_X(i, 0);
      double val_data = matrix_data_Y(i, 0);

      matrix_goodness_cov_total_wiConstraint(i,i) = val_pred;
      if( val_data==1 ) {
	if( val_pred<0.461 ) {// DocDB-32520, when the prediction is sufficiently low
	  double numerator = pow(val_pred-val_data, 2);
	  double dewiminator = 2*( val_pred - val_data + val_data*log(val_data/val_pred) );
	  matrix_goodness_cov_total_wiConstraint(i,i) = numerator/dewiminator;
	}
      }
  
      if( (val_pred==val_data) && (val_pred==0) ) matrix_goodness_cov_total_wiConstraint(i,i) = 1e-6;
    }

    matrix_goodness_cov_total_wiConstraint = matrix_goodness_cov_total_wiConstraint + matrix_YY_under_XX;

    TMatrixD matrix_delta_wiConstraint(1, num_Y);
    for(int ibin=1; ibin<=num_Y; ibin++) {
      double val_pred = matrix_Y_under_X(ibin-1, 0);
      double val_data = matrix_data_Y(ibin-1, 0);
      matrix_delta_wiConstraint(0, ibin-1) = val_pred - val_data;
    }
    TMatrixD matrix_delta_wiConstraint_T(num_Y, 1);
    matrix_delta_wiConstraint_T.Transpose(matrix_delta_wiConstraint);
    TMatrixD matrix_cov_wiConstraint_inv = matrix_goodness_cov_total_wiConstraint;
    matrix_cov_wiConstraint_inv.Invert();
    TMatrixD matrix_chi2_wiConstraint = matrix_delta_wiConstraint * matrix_cov_wiConstraint_inv * matrix_delta_wiConstraint_T;
  
    double val_chi2_wiConstraint = matrix_chi2_wiConstraint(0,0);
    double p_value_wiConstraint = TMath::Prob( val_chi2_wiConstraint, num_Y );  
    cout<<TString::Format(" ---> GOF wiConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.8f",
			  val_chi2_wiConstraint, num_Y, val_chi2_wiConstraint/num_Y, p_value_wiConstraint
			  )<<endl<<endl;
    
    ////////////////////////////////////////
    
    roostr = TString::Format("h1_pred_Y_noConstraint_%02d", index);
    TH1D *h1_pred_Y_noConstraint = new TH1D(roostr, "", num_Y, 0, num_Y);
    for(int ibin=1; ibin<=num_Y; ibin++) {
      h1_pred_Y_noConstraint->SetBinContent( ibin, matrix_pred_Y(ibin-1, 0) );
      double val_err = sqrt( matrix_YY(ibin-1, ibin-1) );
      h1_pred_Y_noConstraint->SetBinError( ibin, val_err );
    }
  
    roostr = TString::Format("h1_pred_Y_wiConstraint_%02d", index);
    TH1D *h1_pred_Y_wiConstraint = new TH1D(roostr, "", num_Y, 0, num_Y);
    for(int ibin=1; ibin<=num_Y; ibin++) {
      h1_pred_Y_wiConstraint->SetBinContent( ibin, matrix_Y_under_X(ibin-1, 0) );
      double val_err = sqrt( matrix_YY_under_XX(ibin-1, ibin-1) );
      h1_pred_Y_wiConstraint->SetBinError( ibin, val_err );
    }

    TGraphAsymmErrors *gh_data = new TGraphAsymmErrors();
    TGraphAsymmErrors *gh_ratio_noConstraint = new TGraphAsymmErrors();    
    TGraphAsymmErrors *gh_ratio_wiConstraint = new TGraphAsymmErrors();    

    double val_max_data = 0;
    
    roostr = TString::Format("h1_ratio_basic_%02d", index);
    TH1D *h1_ratio_basic = new TH1D(roostr, "", num_Y, 0, num_Y);

    roostr = TString::Format("h1_ratio_wi2no_%02d", index);
    TH1D *h1_ratio_wi2no = new TH1D(roostr, "", num_Y, 0, num_Y);
    
    for(int ibin=1; ibin<=num_Y; ibin++) {
      double val_data = matrix_data_Y(ibin-1, 0);
      double val_pred_noConstraint = matrix_pred_Y(ibin-1, 0);
      double val_pred_wiConstraint = matrix_Y_under_X(ibin-1, 0);

      double val_ratio_wi2no = val_pred_wiConstraint/val_pred_noConstraint;
      if( val_ratio_wi2no!=val_ratio_wi2no || val_ratio_wi2no==1./0 ) val_ratio_wi2no = 0;
      h1_ratio_wi2no->SetBinContent(ibin, val_ratio_wi2no);
      
      double val_data_low = 0;
      double val_data_hgh = 0;
      int idx_data = (int)(val_data+0.5);    
      if( idx_data>100 ) {
	val_data_low = val_data - sqrt(val_data);
	val_data_hgh = val_data + sqrt(val_data);
      }
      else {
	val_data_low = DataBase::yl[ idx_data ];
	val_data_hgh = DataBase::yh[ idx_data ];
      }
      
      int n_point = ibin-1;
      double val_x = h1_ratio_basic->GetBinCenter(ibin);
      double val_halfw = h1_ratio_basic->GetBinWidth(ibin)/2;
      
      gh_data->SetPoint( n_point, val_x, val_data );
      gh_data->SetPointError( n_point, val_halfw, val_halfw, val_data-val_data_low, val_data_hgh-val_data );

      if( val_max_data<val_data_hgh ) val_max_data = val_data_hgh;
      
      ///////
      double val_ratio_no = val_data/val_pred_noConstraint;
      double val_ratio_no_low = val_ratio_no - val_data_low/val_pred_noConstraint;
      double val_ratio_no_hgh = val_data_hgh/val_pred_noConstraint - val_ratio_no;
      if( val_ratio_no!=val_ratio_no || val_ratio_no==1./0 ) val_ratio_no = 0;
      gh_ratio_noConstraint->SetPoint( n_point, val_x, val_ratio_no );
      gh_ratio_noConstraint->SetPointError( n_point, val_halfw, val_halfw, val_ratio_no_low, val_ratio_no_hgh );

      ///////
      double val_ratio_wi = val_data/val_pred_wiConstraint;
      double val_ratio_wi_low = val_ratio_wi - val_data_low/val_pred_wiConstraint;
      double val_ratio_wi_hgh = val_data_hgh/val_pred_wiConstraint - val_ratio_wi;
      if( val_ratio_wi!=val_ratio_wi || val_ratio_wi==1./0 ) val_ratio_wi = 0;
      gh_ratio_wiConstraint->SetPoint( n_point, val_x, val_ratio_wi );
      gh_ratio_wiConstraint->SetPointError( n_point, val_halfw, val_halfw, val_ratio_wi_low, val_ratio_wi_hgh );
    }
    
    double val_max_pred =
      h1_pred_Y_noConstraint->GetBinContent( h1_pred_Y_noConstraint->GetMaximumBin() ) +
      h1_pred_Y_noConstraint->GetBinError( h1_pred_Y_noConstraint->GetMaximumBin() );
    
    //////////////////////////////////////// Plotting
    
    roostr = TString::Format("canv_spectra_cmp_%02d", index);
    TCanvas *canv_spectra_cmp = new TCanvas(roostr, roostr, 1000, 950);  
    
    /////////////
    TPad *pad_top = new TPad("pad_top", "pad_top", 0, 0.45, 1, 1);
    func_canv_margin(pad_top, 0.15, 0.1, 0.1, 0.05);
    pad_top->Draw();
    pad_top->cd();
    
    h1_pred_Y_noConstraint->Draw("e2");
    h1_pred_Y_noConstraint->SetMinimum(0);
    if( val_max_data > val_max_pred ) h1_pred_Y_noConstraint->SetMaximum( val_max_data*1.1 );	
    
    TH1D *h1_pred_Y_noConstraint_clone = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_clone");
    
    h1_pred_Y_noConstraint->Draw("same e2");
    h1_pred_Y_noConstraint->SetMarkerStyle(1); h1_pred_Y_noConstraint->SetMarkerColor(color_no);
    h1_pred_Y_noConstraint->SetLineColor(color_no); h1_pred_Y_noConstraint->SetLineWidth(2);
    h1_pred_Y_noConstraint->SetFillColor(color_no); h1_pred_Y_noConstraint->SetFillStyle(3005);
    h1_pred_Y_noConstraint->SetTitle("");
    func_title_size(h1_pred_Y_noConstraint, 0.065, 0.065, 0.065, 0.065);
    func_xy_title(h1_pred_Y_noConstraint, "Bin index", "Entries");
    h1_pred_Y_noConstraint->GetXaxis()->SetNdivisions(506); h1_pred_Y_noConstraint->GetYaxis()->SetNdivisions(506);
    h1_pred_Y_noConstraint->GetYaxis()->CenterTitle(); h1_pred_Y_noConstraint->GetXaxis()->CenterTitle();
    h1_pred_Y_noConstraint->GetXaxis()->SetLabelColor(10);
    h1_pred_Y_noConstraint->GetXaxis()->SetTitleColor(10);
    h1_pred_Y_noConstraint->GetYaxis()->SetTitleOffset(1.1);
    
    h1_pred_Y_noConstraint_clone->Draw("same hist");
    h1_pred_Y_noConstraint_clone->SetLineColor(color_no);
    h1_pred_Y_noConstraint_clone->SetLineWidth(3);
    
    h1_pred_Y_wiConstraint->Draw("same e2");    
    TH1D *h1_pred_Y_wiConstraint_clone = (TH1D*)h1_pred_Y_wiConstraint->Clone("h1_pred_Y_wiConstraint_clone");
    
    h1_pred_Y_wiConstraint->Draw("same e2");
    h1_pred_Y_wiConstraint->SetMarkerStyle(1); h1_pred_Y_wiConstraint->SetMarkerColor(color_wi);
    h1_pred_Y_wiConstraint->SetLineColor(color_wi); h1_pred_Y_wiConstraint->SetLineWidth(2);
    h1_pred_Y_wiConstraint->SetFillColor(color_wi); h1_pred_Y_wiConstraint->SetFillStyle(3004);

    h1_pred_Y_noConstraint_clone->Draw("same hist");
    h1_pred_Y_wiConstraint_clone->Draw("same hist");
    h1_pred_Y_wiConstraint_clone->SetLineColor(color_wi);
    h1_pred_Y_wiConstraint_clone->SetLineWidth(3);

    gh_data->Draw("same pe");
    gh_data->SetMarkerStyle(20); gh_data->SetMarkerSize(1.2); gh_data->SetMarkerColor(kBlack);
    gh_data->SetLineColor(kBlack); gh_data->SetLineWidth(3);

    h1_pred_Y_noConstraint->Draw("same axis");

    TLegend *lg_top = new TLegend(0.5, 0.65, 0.85, 0.85);
    lg_top->AddEntry(gh_data, "Data", "lep");
    lg_top->AddEntry(h1_pred_Y_noConstraint, "Pred no constraint", "lf");
    lg_top->AddEntry(h1_pred_Y_wiConstraint, "Pred wi constraint", "lf");  
    lg_top->Draw();
    lg_top->SetBorderSize(0);
    lg_top->SetFillStyle(0);
    lg_top->SetTextSize(0.065);
    //lg_top->SetTextFont(132);

    pad_top->cd();
    pad_top->Update();
    double x1 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.47;
    double y1 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.45;  
    double x2 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.85;
    double y2 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.65;
    TPaveText *pt = new TPaveText( x1,y1,x2,y2,"l");
    pt->SetTextSize(0.065);
    pt->SetTextFont(42);
    pt->SetTextAlign(11);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d", val_chi2_noConstraint,num_Y));
    ((TText*)pt->GetListOfLines()->Last())->SetTextColor(kRed);
    pt->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d", val_chi2_wiConstraint,num_Y));
    ((TText*)pt->GetListOfLines()->Last())->SetTextColor(kBlue);
    pt->Draw();
    
    /////////////
    canv_spectra_cmp->cd();
  
    TPad *pad_bot = new TPad("pad_bot", "pad_bot", 0, 0, 1, 0.45);
    func_canv_margin(pad_bot, 0.15, 0.1, 0.05, 0.3);
    pad_bot->Draw();
    pad_bot->cd();

    h1_ratio_basic->Draw();
    h1_ratio_basic->SetMinimum(0); h1_ratio_basic->SetMaximum(2);  
    func_title_size(h1_ratio_basic, 0.078, 0.078, 0.078, 0.078);
    func_xy_title(h1_ratio_basic, "Bin index", "Data / Pred");
    h1_ratio_basic->GetXaxis()->SetTickLength(0.05);
    h1_ratio_basic->GetXaxis()->CenterTitle(); h1_ratio_basic->GetYaxis()->CenterTitle(); 
    h1_ratio_basic->GetYaxis()->SetTitleOffset(0.92);
    h1_ratio_basic->GetYaxis()->SetNdivisions(509);
    
    
    TH1D *h1_pred_Y_noConstraint_rel_error = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_rel_error");
    h1_pred_Y_noConstraint_rel_error->Reset();
    for(int ibin=0; ibin<num_Y; ibin++) {
      double val_err = h1_pred_Y_noConstraint->GetBinError(ibin+1);
      double val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin+1);
      double rel_err = val_err/val_cv;
      if( val_cv==0 ) rel_err = 0;
      h1_pred_Y_noConstraint_rel_error->SetBinContent(ibin+1, 1);
      h1_pred_Y_noConstraint_rel_error->SetBinError(ibin+1, rel_err);
    }
    h1_pred_Y_noConstraint_rel_error->Draw("same e2");

    TH1D *h1_pred_Y_wiConstraint_rel_error = (TH1D*)h1_pred_Y_wiConstraint->Clone("h1_pred_Y_wiConstraint_rel_error");
    h1_pred_Y_wiConstraint_rel_error->Reset();
    for(int ibin=0; ibin<num_Y; ibin++) {
      double val_err = h1_pred_Y_wiConstraint->GetBinError(ibin+1);
      //double val_cv = h1_pred_Y_wiConstraint->GetBinContent(ibin+1);
      double val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin+1);
      double rel_err = val_err/val_cv;
      if( val_cv==0 ) rel_err = 0;
      h1_pred_Y_wiConstraint_rel_error->SetBinContent(ibin+1, 1);
      h1_pred_Y_wiConstraint_rel_error->SetBinError(ibin+1, rel_err);
    }
    h1_pred_Y_wiConstraint_rel_error->Draw("same e2");

    TF1 *f1_c1 = new TF1("f1_c1", "1", 0, 1e5);
    f1_c1->Draw("same");
    f1_c1->SetLineStyle(7);
    f1_c1->SetLineWidth(2);
    f1_c1->SetLineColor(kBlack);

    h1_ratio_wi2no->Draw("same hist");
    h1_ratio_wi2no->SetLineColor(kGreen+1);
    h1_ratio_wi2no->SetLineWidth(3);
    pad_bot->Update();
    pad_bot->cd();
    TGaxis *axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),gPad->GetUxmax(),
			      gPad->GetUymax(),0,2,509,"+L");    
    axis->Draw();
    axis->SetLineColor(kGreen+1);
    axis->SetLabelColor(10);
    axis->SetTitleColor(kGreen+1);
    axis->SetTitle("Pred wi/no constraint");
    // axis->SetLabelSize(0.078);
    axis->SetTitleSize(0.078);
    axis->SetTitleOffset(0.15);
    axis->CenterTitle();
      
    gh_ratio_noConstraint->Draw("pe");
    gh_ratio_noConstraint->SetMarkerStyle(20); gh_ratio_noConstraint->SetMarkerSize(1.2);
    gh_ratio_noConstraint->SetMarkerColor(color_no);
    gh_ratio_noConstraint->SetLineColor(color_no); gh_ratio_noConstraint->SetLineWidth(3);
    
    gh_ratio_wiConstraint->Draw("pe");
    gh_ratio_wiConstraint->SetMarkerStyle(20); gh_ratio_wiConstraint->SetMarkerSize(1.2);
    gh_ratio_wiConstraint->SetMarkerColor(color_wi);
    gh_ratio_wiConstraint->SetLineColor(color_wi); gh_ratio_wiConstraint->SetLineWidth(3);
    
    h1_ratio_basic->Draw("same axis");
     
    roostr = TString::Format("canv_spectra_cmp_%02d.png", index);
    canv_spectra_cmp->SaveAs(roostr);
  }// else, num_X!=0
   
}
  
///////////////////////////////////////////////////////// ccc

void TLee::Set_Collapse()
{
  //////////////////////////////////////// pred

  TMatrixD matrix_transform_Lee = matrix_transform;
  for(int ibin=0; ibin<matrix_transform_Lee.GetNrows(); ibin++) {
    for(int jbin=0; jbin<matrix_transform_Lee.GetNcols(); jbin++) {
      if( map_Lee_oldworld.find(ibin)!=map_Lee_oldworld.end() ) matrix_transform_Lee(ibin, jbin) *= scaleF_Lee;
    }
  }
  
  map_pred_spectrum_newworld_bin.clear();
  TMatrixD matrix_pred_oldworld(1, bins_oldworld);
  for(int ibin=0; ibin<bins_oldworld; ibin++) matrix_pred_oldworld(0, ibin) = map_input_spectrum_oldworld_bin[ibin];

  matrix_pred_newworld.Clear();
  matrix_pred_newworld.ResizeTo(1, bins_newworld);
  matrix_pred_newworld = matrix_pred_oldworld * matrix_transform_Lee;
  if( bins_newworld!=matrix_pred_newworld.GetNcols() ) { cerr<<"bins_newworld!=matrix_pred_newworld.GetNcols()"<<endl; exit(1); }
  for(int ibin=0; ibin<bins_newworld; ibin++) map_pred_spectrum_newworld_bin[ibin] = matrix_pred_newworld(0, ibin);
  
  ////////////////////////////////////////
  
  matrix_absolute_cov_oldworld.Clear();
  matrix_absolute_cov_oldworld.ResizeTo( bins_oldworld, bins_oldworld );

  if( flag_syst_flux_Xs ) matrix_absolute_cov_oldworld += matrix_input_cov_flux_Xs;
  if( flag_syst_detector ) matrix_absolute_cov_oldworld += matrix_input_cov_detector;
  if( flag_syst_additional ) matrix_absolute_cov_oldworld += matrix_input_cov_additional;

  TMatrixD matrix_transform_Lee_T( bins_newworld, bins_oldworld );
  matrix_transform_Lee_T.Transpose( matrix_transform_Lee );

  matrix_absolute_cov_newworld.Clear();
  matrix_absolute_cov_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_absolute_cov_newworld = matrix_transform_Lee_T * matrix_absolute_cov_oldworld * matrix_transform_Lee;

  if( flag_syst_mc_stat ) {
    for(int ibin=0; ibin<bins_newworld; ibin++) {
      double val_mc_stat_cov = gh_mc_stat_bin[ibin]->Eval( scaleF_Lee );
      matrix_absolute_cov_newworld(ibin, ibin) += val_mc_stat_cov;
    }
  }
  
}

///////////////////////////////////////////////////////// ccc

void TLee::Set_TransformMatrix()
{
   cout<<endl<<" ---> Set_TransformMatrix"<<endl<<endl;

   ////////////////////////////// correponding to "Set_Spectra_MatrixCov"

}

///////////////////////////////////////////////////////// ccc

void TLee::Set_POT_implement()
{
  cout<<endl<<" ---> Set_POT_implement"<<endl<<endl;
  
  ////////////////////////////// pred

  int line_pred = -1;
  for( auto it_ch=map_input_spectrum_ch_bin.begin(); it_ch!=map_input_spectrum_ch_bin.end(); it_ch++ ) {
    int ich = it_ch->first;
    for(int ibin=0; ibin<(int)map_input_spectrum_ch_bin[ich].size(); ibin++) {
      line_pred++;
      map_input_spectrum_ch_bin[ich][ibin] *= scaleF_POT;
      map_input_spectrum_oldworld_bin[line_pred] *= scaleF_POT;
    }// ibin
  }// ich
  
  ////////////////////////////// data

  int line_data = -1;
  for( auto it_ch=map_data_spectrum_ch_bin.begin(); it_ch!=map_data_spectrum_ch_bin.end(); it_ch++ ) {
    int ich = it_ch->first;
    for( int ibin=0; ibin<(int)map_data_spectrum_ch_bin[ich].size(); ibin++ ) {
      line_data++;
      map_data_spectrum_ch_bin[ich][ibin] *= scaleF_POT;
      map_data_spectrum_newworld_bin[line_data] *= scaleF_POT;
    }// ibin
  }// ich

  matrix_data_newworld.Clear();
  matrix_data_newworld.ResizeTo(1, bins_newworld);
  for(int ibin=0; ibin<bins_newworld; ibin++) {
    matrix_data_newworld(0, ibin) = map_data_spectrum_newworld_bin[ibin];
  }
    
  ////////////////////////////// flux_Xs, detector, additional, mc_stat

  double scaleF_POT2 = scaleF_POT * scaleF_POT;
  
  for(int ibin=0; ibin<bins_oldworld; ibin++) {
    for(int jbin=0; jbin<bins_oldworld; jbin++) {      
      matrix_input_cov_flux_Xs(ibin, jbin) *= scaleF_POT2;
      matrix_input_cov_detector(ibin, jbin) *= scaleF_POT2;
      matrix_input_cov_additional(ibin, jbin) *= scaleF_POT2;      
    }// jbin
  }// ibin

  for(auto it=gh_mc_stat_bin.begin(); it!=gh_mc_stat_bin.end(); it++) {
    int ibin = it->first; //cout<<Form(" ---> check %3d, %3d", ibin, gh_mc_stat_bin[ibin]->GetN())<<endl;
    for(int idx=0; idx<gh_mc_stat_bin[ibin]->GetN(); idx++) {
      double x(0), y(0);
      gh_mc_stat_bin[ibin]->GetPoint(idx, x, y);
      gh_mc_stat_bin[ibin]->SetPoint(idx, x, y*scaleF_POT2);
    }// ipoint
  }// ibin
  
}

///////////////////////////////////////////////////////// ccc

void TLee::Set_config_file_directory(TString spectra_file_, TString flux_Xs_directory_, TString detector_directory_, TString mc_directory_)
{
  cout<<endl<<" ---> Set_config_file_directory"<<endl<<endl;

  spectra_file       = spectra_file_;
  flux_Xs_directory  = flux_Xs_directory_;
  detector_directory = detector_directory_;
  mc_directory       = mc_directory_;

  cout<<Form(" spectra_file       %-10s", spectra_file.Data() )<<endl;
  cout<<Form(" flux_Xs_directory  %-10s", flux_Xs_directory.Data() )<<endl;
  cout<<Form(" detector_directory %-10s", detector_directory.Data() )<<endl;
  cout<<Form(" mc_directory       %-10s", mc_directory.Data() )<<endl;  
}


void TLee::Set_Spectra_MatrixCov()
{
  /// spectra should be consist with matrix-cov order
  
  cout<<endl<<" ---> Set_Spectra_MatrixCov"<<endl<<endl;
  TString roostr = "";

  ////////////////////////////////////// pred
  
  // https://www.phy.bnl.gov/xqian/talks/wire-cell/Leeana/configurations/cov_input.txt  
  map_input_spectrum_ch_str[1] = "nueCC_FC_norm";
  map_input_spectrum_ch_str[2] = "nueCC_PC_norm";
  map_input_spectrum_ch_str[3] = "numuCC_FC_norm";
  map_input_spectrum_ch_str[4] = "numuCC_PC_norm";
  map_input_spectrum_ch_str[5] = "CCpi0_FC_norm";
  map_input_spectrum_ch_str[6] = "CCpi0_PC_norm";
  map_input_spectrum_ch_str[7] = "NCpi0_norm";
  map_input_spectrum_ch_str[8] = "Lee_FC";
  map_input_spectrum_ch_str[9] = "Lee_PC";
  map_input_spectrum_ch_str[10]= "nueCC_FC_ext";
  map_input_spectrum_ch_str[11]= "nueCC_PC_ext";
  map_input_spectrum_ch_str[12]= "numuCC_FC_ext";
  map_input_spectrum_ch_str[13]= "numuCC_PC_ext";
  map_input_spectrum_ch_str[14]= "CCpi0_FC_ext";
  map_input_spectrum_ch_str[15]= "CCpi0_PC_ext";
  map_input_spectrum_ch_str[16]= "NCpi0_ext";

  map_Lee_ch[8] = 1;
  map_Lee_ch[9] = 1;

  //////////////////
  //////////////////
  
  roostr = spectra_file;
  TFile *file_spectra = new TFile(roostr, "read");

  ///
  TMatrixD *mat_collapse = (TMatrixD*)file_spectra->Get("mat_collapse");
  matrix_transform.Clear();
  matrix_transform.ResizeTo( mat_collapse->GetNrows(), mat_collapse->GetNcols() );
  matrix_transform = (*mat_collapse);

  ///
  for(int ich=1; ich<=(int)map_input_spectrum_ch_str.size(); ich++) {
    roostr = TString::Format("histo_%d", ich);
    TH1F *h1_spectrum = (TH1F*)file_spectra->Get(roostr);
    int bins = h1_spectrum->GetNbinsX() + 1;    
    cout<<Form(" %2d  %-20s   bin-num %2d", ich, map_input_spectrum_ch_str[ich].Data(), bins)<<endl;
    
    for(int ibin=1; ibin<=bins; ibin++) map_input_spectrum_ch_bin[ich][ibin-1] = h1_spectrum->GetBinContent(ibin);
  }
  cout<<endl;

  bins_oldworld = 0;
  for(auto it_ch=map_input_spectrum_ch_bin.begin(); it_ch!=map_input_spectrum_ch_bin.end(); it_ch++) {
    int ich = it_ch->first;
      for(int ibin=0; ibin<(int)map_input_spectrum_ch_bin[ich].size(); ibin++) {
	bins_oldworld++;
	int index_oldworld = bins_oldworld - 1;	
	map_input_spectrum_oldworld_bin[ index_oldworld ] = map_input_spectrum_ch_bin[ich][ibin];
	if( map_Lee_ch.find(ich)!=map_Lee_ch.end() ) map_Lee_oldworld[index_oldworld] = 1;
    }// ibin
  }// ich

  ////////////////////////////////////// data

  int line_data = -1;
  bins_newworld = 0;
  for(int ich=1; ich<=7; ich++) {
    roostr = TString::Format("hdata_obsch_%d", ich);
    TH1F *h1_spectrum = (TH1F*)file_spectra->Get(roostr);
    for(int ibin=1; ibin<=h1_spectrum->GetNbinsX()+1; ibin++) {
      map_data_spectrum_ch_bin[ich][ibin-1] = h1_spectrum->GetBinContent(ibin);

      line_data++;
      bins_newworld++;
      map_data_spectrum_newworld_bin[line_data] = map_data_spectrum_ch_bin[ich][ibin-1]; 
    }// ibin
  }// ich

  ////////////////////////////////////////// flux_Xs

  map<int, TFile*>map_file_flux_Xs_frac;  
  map<int, TMatrixD*>map_matrix_flux_Xs_frac;
  
  TMatrixD matrix_flux_Xs_frac(bins_oldworld, bins_oldworld);

  for(int idx=1; idx<=17; idx++) {
    roostr = TString::Format(flux_Xs_directory+"cov_%d.root", idx);
    map_file_flux_Xs_frac[idx] = new TFile(roostr, "read");
    map_matrix_flux_Xs_frac[idx] = (TMatrixD*)map_file_flux_Xs_frac[idx]->Get(TString::Format("frac_cov_xf_mat_%d", idx));
    // cout<<TString::Format(" ---> check: flux and Xs, %2d  ", idx)<<roostr<<endl;
    matrix_flux_Xs_frac += (*map_matrix_flux_Xs_frac[idx]);
  }
  // cout<<endl;
  
  ////////////////////////////////////////// detector
  
  map<int, TString>map_detectorfile_str;
  map_detectorfile_str[1] = detector_directory+"cov_LYDown.root";
  map_detectorfile_str[2] = detector_directory+"cov_LYRayleigh.root";
  map_detectorfile_str[3] = detector_directory+"cov_Recomb2.root";
  map_detectorfile_str[4] = detector_directory+"cov_SCE.root";
  map_detectorfile_str[5] = detector_directory+"cov_WMdEdx.root";
  map_detectorfile_str[6] = detector_directory+"cov_WMThetaXZ.root";
  map_detectorfile_str[7] = detector_directory+"cov_WMThetaYZ.root";
  map_detectorfile_str[8] = detector_directory+"cov_WMX.root";
  map_detectorfile_str[9] = detector_directory+"cov_WMYZ.root";
  map_detectorfile_str[10]= detector_directory+"cov_LYatt.root";
  
  map<int, TFile*>map_file_detector_frac;
  map<int, TMatrixD*>map_matrix_detector_frac;
  TMatrixD matrix_detector_frac(bins_oldworld, bins_oldworld);

  int size_map_detectorfile_str = map_detectorfile_str.size();
  for(int idx=1; idx<=size_map_detectorfile_str; idx++) {
    if(idx==5) continue;
    roostr = map_detectorfile_str[idx];
    map_file_detector_frac[idx] = new TFile(roostr, "read");
    map_matrix_detector_frac[idx] = (TMatrixD*)map_file_detector_frac[idx]->Get(TString::Format("frac_cov_det_mat_%d", idx));
    // cout<<TString::Format(" ---> check: detector, %2d  ", idx)<<roostr<<endl;

    matrix_detector_frac += (*map_matrix_detector_frac[idx]);
  }
  // cout<<endl;

  ////////////////////////////////////////// additional

  TMatrixD *matrix_additional_abs_point = (TMatrixD*)file_spectra->Get("cov_mat_add");
  TMatrixD matrix_additional_abs = (*matrix_additional_abs_point);
    
  //////////////////////////////////////////

  matrix_input_cov_flux_Xs.Clear();
  matrix_input_cov_detector.Clear();
  matrix_input_cov_additional.Clear();
  
  matrix_input_cov_flux_Xs.ResizeTo( bins_oldworld, bins_oldworld );
  matrix_input_cov_detector.ResizeTo( bins_oldworld, bins_oldworld );
  matrix_input_cov_additional.ResizeTo( bins_oldworld, bins_oldworld );

  for(int ibin=0; ibin<bins_oldworld; ibin++) {
    for(int jbin=0; jbin<bins_oldworld; jbin++) {
      double val_i = map_input_spectrum_oldworld_bin[ibin];
      double val_j = map_input_spectrum_oldworld_bin[jbin];
      double val_cov = 0;
      
      val_cov = matrix_flux_Xs_frac(ibin, jbin);
      matrix_input_cov_flux_Xs(ibin, jbin) = val_cov * val_i * val_j;
      
      val_cov = matrix_detector_frac(ibin, jbin);
      matrix_input_cov_detector(ibin, jbin) = val_cov * val_i * val_j;
    }
  }

  matrix_input_cov_additional = matrix_additional_abs;
  
  ////////////////////////////////////////// MC statistics

  map<int, map<int, double> >map_mc_stat_file_bin_Lee;
  map<int, map<int, double> >map_mc_stat_file_bin_mcStat;
  int gbins_mc_stat = 137;
  
  for(int ifile=0; ifile<=99; ifile++) {
    roostr = TString::Format(mc_directory+"%d.log", ifile);
    ifstream InputFile_aa(roostr, ios::in);
    if(!InputFile_aa) { cerr<<" No input-list"<<endl; exit(1); }

    int line = 0;    
    double Lee = 1; double run = 1;
    
    for(int idx=1; idx<=gbins_mc_stat+1; idx++) {            
      int gbin = 0; int lbin = 0; double val_pred = 0; double mc_stat = 0; double nn_stat = 0;
      if(idx==1) { InputFile_aa>>Lee>>run; }
      else {
	InputFile_aa>>gbin>>lbin>>val_pred>>mc_stat>>nn_stat;
	line++;
	map_mc_stat_file_bin_Lee[ifile][line-1] = Lee;
	map_mc_stat_file_bin_mcStat[ifile][line-1] = mc_stat;
      }
    }
  }
  
  /// gh_mc_stat_bin
  for(int ibin=0; ibin<gbins_mc_stat; ibin++) {
    gh_mc_stat_bin[ibin] = new TGraph(); gh_mc_stat_bin[ibin]->SetName(TString::Format("gh_mc_stat_bin_%03d", ibin));
    
    for(auto it=map_mc_stat_file_bin_Lee.begin(); it!=map_mc_stat_file_bin_Lee.end(); it++) {
      int ifile = it->first;
      double Lee = map_mc_stat_file_bin_Lee[ifile][ibin];
      double mc_stat = map_mc_stat_file_bin_mcStat[ifile][ibin];
      gh_mc_stat_bin[ibin]->SetPoint( gh_mc_stat_bin[ibin]->GetN(), Lee, mc_stat );
    }
    
    double x,y;
    gh_mc_stat_bin[ibin]->GetPoint( gh_mc_stat_bin[ibin]->GetN()-1, x, y);
    gh_mc_stat_bin[ibin]->SetPoint( gh_mc_stat_bin[ibin]->GetN(), x+1, y);
  }  
  
}
