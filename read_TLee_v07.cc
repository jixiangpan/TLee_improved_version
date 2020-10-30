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

/// minuit2
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"


/////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "draw.icc"

///////
#include <chrono>
using namespace std::chrono;

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

///////////////////////////////////////////////////////////////////////////////////////////////////////// TLee

class TLee
{
public:
  TLee() {
    rand = new TRandom3(0);
    
    bins_total_old_world = 0;
    
    flag_syst_flux_Xs = false;
    flag_syst_detector = false;
    flag_syst_additional = false;
    flag_syst_mc_stat = false;

  }

  /////////////////////////////////////////////////////// data memeber

  TRandom3 *rand;
  
  double scaleF_POT;
  double scaleF_Lee;

  map<int, map<int,double> >map_input_spectrum_ch_bin;  
  map<int, TString>map_input_spectrum_ch_str;
  map<int, int>map_input_spectrum_ch_binnum;
  map<int, double>map_input_spectrum_global_bin;/// 
  map<int, TH1F*>map_input_h1_spectrum;  
  int bins_total_old_world;

  map<int, map<int, int> >map_old_world_local2global;
  map<int, int>map_old_world_global2localch;
  map<int, int>map_old_world_global2localbin;
  
  map<int, int>map_Lee_ch;
  map<int, int>map_Lee_old_world;

  TMatrixD matrix_trans_wiLee;
  TMatrixD matrix_trans_noLee;
  
  map<int, map<int, double> >map_data_spectrum_ch_bin;
  map<int, TH1F*>map_data_h1_spectrum;
  map<int, double>map_data_spectrum_global_bin;
  
  TMatrixD matrix_input_cov_flux_Xs;
  TMatrixD matrix_input_cov_detector;
  TMatrixD matrix_input_cov_additional;

  TMatrixD matrix_input_cov_flux;
  TMatrixD matrix_input_cov_Xs;
  
  map<int, TGraph*>gh_mc_stat_bin;

  ///////// Lee and POT scaled
  
  map<int, map<int,double> >map_input_spectrum_ch_bin_scaled;
  map<int, double>map_input_spectrum_global_bin_scaled;
  
  map<int, map<int,double> >map_data_spectrum_ch_bin_scaled;
  map<int, double>map_data_spectrum_global_bin_scaled;// -------------------->
  
  TMatrixD matrix_input_cov_flux_Xs_scaled;
  TMatrixD matrix_input_cov_detector_scaled;
  TMatrixD matrix_input_cov_additional_scaled;
  TMatrixD matrix_input_cov_total_scaled;

  TMatrixD matrix_input_cov_flux_scaled;
  TMatrixD matrix_input_cov_Xs_scaled;

  TMatrixD matrix_collapse_absolute_covariance_wiLee_flux;
  TMatrixD matrix_collapse_absolute_covariance_wiLee_Xs;
  TMatrixD matrix_collapse_absolute_covariance_wiLee_detector;
  TMatrixD matrix_collapse_absolute_covariance_wiLee_additional;
  TMatrixD matrix_collapse_absolute_covariance_wiLee_mc_stat;
  
  bool flag_syst_flux_Xs;
  bool flag_syst_detector;
  bool flag_syst_additional;
  bool flag_syst_mc_stat;
  
  int num_elements_collapse;
  int num_global_elements;
  
  map<int, map<int, int> >map_collapse_channel_local2global;
  
  TMatrixD matrix_collapse_absolute_covariance_wiLee;// -------------------->
  map<int, double>map_collapsed_prediction_wiLee;// --------------------> index from 1
  
  TMatrixD matrix_collapse_absolute_covariance_noLee;// -------------------->
  map<int, double>map_collapsed_prediction_noLee;// --------------------> index from 1

  //////////
  map<int, double>map_fake_data;// start from index 0
  map<int, double>map_Asimov_data;// start from index 0
  map<int, double>map_measured_data;// start from index 0

  int N_toy;
  
  map<int, map<int, double > >map_meas_variation_TotalSyst_noLee;// itoy, ibin index from 0
  map<int, map<int, double > >map_meas_variation_TotalSyst_wiLee;// itoy, ibin index from 0
  
  ///
  int minimization_status;
  double minimization_chi2;
  double minimization_Lee_strength_val;
  double minimization_Lee_strength_err;

  /////////////////////////////////////////////////////// function member

  void Set_Spectra_MatrixCov();
  void Set_POT_implement();
  void Set_TransformMatrix();
  
  void Set_Collapse();
  
  void Exe_Goodness_of_Fit(map<int, double>map_h1_pred, map<int, double>map_h1_data, TMatrixD matrix_syst, int index);

  //// Y constrained by X
  void Exe_Constraint_GOF(int num_Y, int num_X, map<int, double>map_h1_pred, map<int, double>map_h1_data, TMatrixD matrix_syst, int index);

  void Set_fake_data_Asimov();
  void Set_Variations(int num_toy);
  void Set_fake_data_Variation(int itoy);
    
  void Set_measured_data();
  
  void Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed);
  
};

///////////////////////////// ccc

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
      TMatrixD matrix_meas(1, num_elements_collapse);
      for(int idx=0; idx<num_elements_collapse; idx++) {
	matrix_meas(0, idx) = map_fake_data[idx];	
      }

      /////////
      TMatrixD matrix_pred(1, num_elements_collapse);

      scaleF_Lee = Lee_strength;
      Set_Collapse();
     
      for(int idx=0; idx<num_elements_collapse; idx++) {
	matrix_pred(0, idx) = map_collapsed_prediction_wiLee[idx+1];// ------------>	
      }

      /////////
      TMatrixD matrix_cov_syst = matrix_collapse_absolute_covariance_wiLee;
      
      for(int idx=0; idx<num_elements_collapse; idx++) {

	double val_stat_cov = 0;	
	double val_meas = matrix_meas(0, idx);
	double val_pred = matrix_pred(0, idx);
	
	if( val_meas==0 ) val_stat_cov = val_pred/2;
	else val_stat_cov = 3./( 1./val_meas + 2./val_pred );	
	if( val_meas==0 && val_pred==0 ) val_stat_cov = 1;
	matrix_cov_syst(idx, idx) += val_stat_cov;
      }

      TMatrixD matrix_cov_total = matrix_cov_syst;
      TMatrixD matrix_cov_total_inv = matrix_cov_total;
      matrix_cov_total_inv.Invert();
      
      ////////
      TMatrixD matrix_delta = matrix_pred - matrix_meas;
      TMatrixD matrix_delta_T( num_elements_collapse, 1 );
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
  // 			  minimization_Lee_strength_val, minosError_low, minosError_hgh)<<endl;
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

///////////////////////////// ccc
void TLee::Set_Variations(int num_toy)
{
  N_toy = num_toy;
  
  map_meas_variation_TotalSyst_noLee.clear();
  map_meas_variation_TotalSyst_wiLee.clear();
    
  ////////////////////////////// wiLee
  ////////////////////////////// wiLee

  TMatrixD matrix_collapse_absolute_covariance_wiLee_temp_rel = matrix_collapse_absolute_covariance_wiLee;
  
  TMatrixDSym DSmatrix_cov_wiLee(num_elements_collapse);
  for(int ibin=0; ibin<num_elements_collapse; ibin++) {
    for(int jbin=0; jbin<num_elements_collapse; jbin++) {
      DSmatrix_cov_wiLee(ibin, jbin) = matrix_collapse_absolute_covariance_wiLee_temp_rel(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen_wiLee( DSmatrix_cov_wiLee );
  TMatrixD matrix_eigenvector_wiLee = DSmatrix_eigen_wiLee.GetEigenVectors();
  TVectorD matrix_eigenvalue_wiLee = DSmatrix_eigen_wiLee.GetEigenValues();
  // TMatrixD matrix_eigenvector_wiLee_T(num_elements_collapse, num_elements_collapse);
  // matrix_eigenvector_wiLee_T.Transpose( matrix_eigenvector_wiLee );
  // TMatrixD matrix_cov_wiLee_diag = matrix_eigenvector_wiLee_T * matrix_collapse_absolute_covariance_wiLee_temp_rel * matrix_eigenvector_wiLee;
  // for(int i=0; i<num_elements_collapse; i++) {
  //   for(int j=0; j<num_elements_collapse; j++) {
  //     if( fabs(matrix_cov_wiLee_diag(i,j))<1e-9 ) matrix_cov_wiLee_diag(i,j) = 0;
  //   }
  // }
   
  double *array_variation_wiLee = new double[num_elements_collapse];
  for(int itoy=1; itoy<=N_toy; itoy++) {    
    TMatrixD matrix_element(num_elements_collapse, 1);    
    for(int j=0; j<num_elements_collapse; j++) {
      if( matrix_eigenvalue_wiLee(j)>=0 ) {
	matrix_element(j,0) = rand->Gaus( 0, sqrt( matrix_eigenvalue_wiLee(j) ) );
      }
      else {
	matrix_element(j,0) = 0;
      }      
    }
    TMatrixD matrix_variation = matrix_eigenvector_wiLee * matrix_element;
    for(int j=0; j<num_elements_collapse; j++) {
      array_variation_wiLee[j] = matrix_variation(j,0);
      double val_with_syst = array_variation_wiLee[j] + map_collapsed_prediction_wiLee[j+1];// key point
      if( val_with_syst<0 ) val_with_syst = 0;
      map_meas_variation_TotalSyst_wiLee[itoy][j] = rand->PoissonD( val_with_syst );
    }
  }
  
  delete[]  array_variation_wiLee;
  
  ////////////////////////////// noLee
  ////////////////////////////// noLee
  /*
  TMatrixD matrix_collapse_absolute_covariance_noLee_temp_rel = matrix_collapse_absolute_covariance_noLee;
  
  TMatrixDSym DSmatrix_cov_noLee(num_elements_collapse);
  for(int ibin=0; ibin<num_elements_collapse; ibin++) {
    for(int jbin=0; jbin<num_elements_collapse; jbin++) {
      DSmatrix_cov_noLee(ibin, jbin) = matrix_collapse_absolute_covariance_noLee_temp_rel(ibin, jbin);
    }
  }
  TMatrixDSymEigen DSmatrix_eigen_noLee( DSmatrix_cov_noLee );
  TMatrixD matrix_eigenvector_noLee = DSmatrix_eigen_noLee.GetEigenVectors();
  TVectorD matrix_eigenvalue_noLee = DSmatrix_eigen_noLee.GetEigenValues();
  // TMatrixD matrix_eigenvector_noLee_T(num_elements_collapse, num_elements_collapse);
  // matrix_eigenvector_noLee_T.Transpose( matrix_eigenvector_noLee );
  // TMatrixD matrix_cov_noLee_diag = matrix_eigenvector_noLee_T * matrix_collapse_absolute_covariance_noLee_temp_rel * matrix_eigenvector_noLee;
  // for(int i=0; i<num_elements_collapse; i++) {
  //   for(int j=0; j<num_elements_collapse; j++) {
  //     if( fabs(matrix_cov_noLee_diag(i,j))<1e-9 ) matrix_cov_noLee_diag(i,j) = 0;
  //   }
  // }
   
  double *array_variation_noLee = new double[num_elements_collapse];
  for(int itoy=1; itoy<=N_toy; itoy++) {    
    TMatrixD matrix_element(num_elements_collapse, 1);    
    for(int j=0; j<num_elements_collapse; j++) {
      if( matrix_eigenvalue_noLee(j)>=0 ) {
	matrix_element(j,0) = rand->Gaus( 0, sqrt( matrix_eigenvalue_noLee(j) ) );
      }
      else {
	matrix_element(j,0) = 0;
      }      
    }
    TMatrixD matrix_variation = matrix_eigenvector_noLee * matrix_element;
    for(int j=0; j<num_elements_collapse; j++) {
      array_variation_noLee[j] = matrix_variation(j,0);
      double val_noth_syst = array_variation_noLee[j] + map_collapsed_prediction_noLee[j+1];// key point
      if( val_noth_syst<0 ) val_noth_syst = 0;
      map_meas_variation_TotalSyst_noLee[itoy][j] = rand->PoissonD( val_noth_syst );
    }
  }
  
  delete[]  array_variation_noLee;
  */
}

///////////////////////////// ccc

void TLee::Set_measured_data()
{
  map_fake_data.clear();
  map_measured_data.clear();

  int rows = matrix_collapse_absolute_covariance_wiLee.GetNrows();
  for(int ibin=1; ibin<=rows; ibin++) {
    map_measured_data[ibin-1] = map_data_spectrum_global_bin_scaled[ibin];
    map_fake_data[ibin-1] = map_measured_data[ibin-1];
  }  
}

void TLee::Set_fake_data_Variation(int itoy)
{
  for(int idx=0; idx<num_elements_collapse; idx++) {
    map_fake_data[idx] = map_meas_variation_TotalSyst_wiLee[itoy][idx];
  } 
}

void TLee::Set_fake_data_Asimov()
{
  map_fake_data.clear();
  map_Asimov_data.clear();
  
  int rows = matrix_collapse_absolute_covariance_wiLee.GetNrows();
  for(int ibin=1; ibin<=rows; ibin++) {
    map_Asimov_data[ibin-1] = map_collapsed_prediction_wiLee[ibin];
    map_fake_data[ibin-1] = map_Asimov_data[ibin-1];
  }  
}

///////////////////////////// ccc

void TLee::Exe_Constraint_GOF(int num_Y, int num_X, map<int, double>map_h1_pred, map<int, double>map_h1_data, TMatrixD matrix_syst, int index)
{  
  TString roostr = "";

  int size_map_h1_pred = map_h1_pred.size();
  int size_map_h1_data = map_h1_data.size();
  if( size_map_h1_pred!=size_map_h1_data ) {
    cerr<<" Error: Exe_Constraint_GOF: size_map_h1_data!=size_map_h1_pred"<<endl; exit(1);
  }
    
  if( num_Y+num_X!=size_map_h1_pred ) {
    cerr<<" Error: Exe_Constraint_GOF: num_Y+num_X!=size_map_h1_pred"<<endl; exit(1);
  }

  //////////////////////////////////

  TMatrixD matrix_pred_Y(num_Y, 1);
  TMatrixD matrix_data_Y(num_Y, 1);
  
  TMatrixD matrix_pred_X(num_X, 1);
  TMatrixD matrix_data_X(num_X, 1);
  
  for(int ibin=1; ibin<=num_Y; ibin++) {
    matrix_pred_Y(ibin-1, 0) = map_h1_pred[ibin];
    matrix_data_Y(ibin-1, 0) = map_h1_data[ibin];
  }

  for(int ibin=1; ibin<=num_X; ibin++) {
    matrix_pred_X(ibin-1, 0) = map_h1_pred[num_Y+ibin];
    matrix_data_X(ibin-1, 0) = map_h1_data[num_Y+ibin];
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
      if(ibin==jbin) matrix_XX(ibin-1, jbin-1) += matrix_pred_X(ibin-1, 0);// Pearson's term
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
    //cout<<Form(" ---> %3d, pred/data %10.6f %10.6f", ibin, val_pred, val_data)<<endl;
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
                              )<<endl<<endl;

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
  cout<<endl<<TString::Format(" ---> GOF wiConstraint: chi2 %6.2f, ndf %3d, chi2/ndf %6.2f, p-value %10.8f",
                              val_chi2_wiConstraint, num_Y, val_chi2_wiConstraint/num_Y, p_value_wiConstraint
                              )<<endl<<endl;
  
  
  /////////////////////////////

  roostr = TString::Format("h1_data_%02d", index);
  TH1D *h1_data = new TH1D(roostr, roostr, num_Y, 0, num_Y);
  for(int ibin=1; ibin<=num_Y; ibin++) {
    h1_data->SetBinContent( ibin, matrix_data_Y(ibin-1, 0) );
  }
  
  roostr = TString::Format("h1_pred_Y_noConstraint_%02d", index);
  TH1D *h1_pred_Y_noConstraint = new TH1D(roostr, roostr, num_Y, 0, num_Y);
  for(int ibin=1; ibin<=num_Y; ibin++) {
    h1_pred_Y_noConstraint->SetBinContent( ibin, matrix_pred_Y(ibin-1, 0) );
    double val_err = sqrt( matrix_cov_total(ibin-1, ibin-1) );
    h1_pred_Y_noConstraint->SetBinError( ibin, val_err );
  }
  
  roostr = TString::Format("h1_pred_Y_wiConstraint_%02d", index);
  TH1D *h1_pred_Y_wiConstraint = new TH1D(roostr, roostr, num_Y, 0, num_Y);
  for(int ibin=1; ibin<=num_Y; ibin++) {
    h1_pred_Y_wiConstraint->SetBinContent( ibin, matrix_Y_under_X(ibin-1, 0) );
    double val_err = sqrt( matrix_YY_under_XX(ibin-1, ibin-1) );
    h1_pred_Y_wiConstraint->SetBinError( ibin, val_err );
  }
  
  ///////////////

  TH1D *h1_ratio_basic = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_ratio_basic");
  h1_ratio_basic->Reset();
  h1_ratio_basic->SetTitle("");

  TGraphAsymmErrors *gh_ratio_noConstraint = new TGraphAsymmErrors();
  TGraphAsymmErrors *gh_data = new TGraphAsymmErrors();
  
  for(int ibin=1; ibin<=h1_ratio_basic->GetNbinsX(); ibin++) {    
    double val_data = h1_data->GetBinContent(ibin);
    double val_pred = h1_pred_Y_noConstraint->GetBinContent(ibin);
    double val_ratio = val_data/val_pred;

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
    
    int n_point = gh_ratio_noConstraint->GetN();
    double val_x = h1_ratio_basic->GetBinCenter(ibin);
    double val_halfw = h1_ratio_basic->GetBinWidth(ibin)/2;
    gh_ratio_noConstraint->SetPoint( n_point, val_x, val_ratio );
    gh_ratio_noConstraint->SetPointError( n_point, val_halfw, val_halfw, val_ratio_low, val_ratio_hgh );

    gh_data->SetPoint( n_point, val_x, val_data );
    gh_data->SetPointError( n_point, val_halfw, val_halfw, val_data-val_data_low, val_data_hgh-val_data );
  }

  TH1D *h1_pred_Y_noConstraint_rel_error = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_rel_error");
  h1_pred_Y_noConstraint_rel_error->Reset();
  h1_pred_Y_noConstraint_rel_error->SetTitle("");
  for(int ibin=1; ibin<=h1_pred_Y_noConstraint_rel_error->GetNbinsX(); ibin++) {
    double val_err = h1_pred_Y_noConstraint->GetBinError(ibin);
    double val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin);
    double rel_err = val_err/val_cv;
    if( val_cv==0 ) rel_err = 0;
    h1_pred_Y_noConstraint_rel_error->SetBinContent(ibin, 1);
    h1_pred_Y_noConstraint_rel_error->SetBinError(ibin, rel_err);
  }

  ///////
  
  TGraphAsymmErrors *gh_ratio_wiConstraint = new TGraphAsymmErrors();
  
  for(int ibin=1; ibin<=h1_ratio_basic->GetNbinsX(); ibin++) {    
    double val_data = h1_data->GetBinContent(ibin);
    double val_pred = h1_pred_Y_wiConstraint->GetBinContent(ibin);
    double val_ratio = val_data/val_pred;

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
    
    int n_point = gh_ratio_wiConstraint->GetN();
    double val_x = h1_ratio_basic->GetBinCenter(ibin);
    double val_halfw = h1_ratio_basic->GetBinWidth(ibin)/2;
    gh_ratio_wiConstraint->SetPoint( n_point, val_x, val_ratio );
    gh_ratio_wiConstraint->SetPointError( n_point, val_halfw, val_halfw, val_ratio_low, val_ratio_hgh );  
  }

  TH1D *h1_pred_Y_wiConstraint_rel_error = (TH1D*)h1_pred_Y_wiConstraint->Clone("h1_pred_Y_wiConstraint_rel_error");
  h1_pred_Y_wiConstraint_rel_error->Reset();
  h1_pred_Y_wiConstraint_rel_error->SetTitle("");
  for(int ibin=1; ibin<=h1_pred_Y_wiConstraint_rel_error->GetNbinsX(); ibin++) {
    double val_err = h1_pred_Y_wiConstraint->GetBinError(ibin);
    double val_cv = h1_pred_Y_wiConstraint->GetBinContent(ibin);
    val_cv = h1_pred_Y_noConstraint->GetBinContent(ibin);
    double rel_err = val_err/val_cv;
    if( val_cv==0 ) rel_err = 0;
    h1_pred_Y_wiConstraint_rel_error->SetBinContent(ibin, 1);
    h1_pred_Y_wiConstraint_rel_error->SetBinError(ibin, rel_err);
  }

  ///////////////////////////////////////////////////////////////////////// Plotting  
  
  double ymax_pred = h1_pred_Y_noConstraint->GetMaximum();
  double ymax_data = h1_data->GetMaximum();
  double ymax_user = (ymax_pred>ymax_data) ? ymax_pred : ymax_data;
  
  roostr = TString::Format("canv_spectra_cmp_%02d", index);
  TCanvas *canv_spectra_cmp = new TCanvas(roostr, roostr, 1000, 950);  
  //func_canv_margin(canv_spectra_cmp, 0.15, 0.1,0.1,0.15);

  /////////////
  TPad *pad_top = new TPad("pad_top", "pad_top", 0, 0.45, 1, 1);
  func_canv_margin(pad_top, 0.15, 0.1, 0.1, 0.05);
  pad_top->Draw();
  pad_top->cd();

  ////////
  h1_pred_Y_noConstraint->Draw("e2");
  if( ymax_pred<ymax_data ) h1_pred_Y_noConstraint->SetMaximum(ymax_user*1.3);
  h1_pred_Y_noConstraint->SetMinimum(0);
  h1_pred_Y_noConstraint->SetTitle("");
  func_title_size(h1_pred_Y_noConstraint, 0.065, 0.065, 0.065, 0.065);
  func_xy_title(h1_pred_Y_noConstraint, "Index", "Entries");
  h1_pred_Y_noConstraint->GetXaxis()->SetNdivisions(506);
  h1_pred_Y_noConstraint->GetYaxis()->SetNdivisions(508);
  h1_pred_Y_noConstraint->GetYaxis()->CenterTitle();
  h1_pred_Y_noConstraint->GetXaxis()->CenterTitle();
  h1_pred_Y_noConstraint->GetXaxis()->SetLabelColor(10);
  h1_pred_Y_noConstraint->GetXaxis()->SetTitleColor(10);  
  TH1D *h1_pred_Y_noConstraint_clone = (TH1D*)h1_pred_Y_noConstraint->Clone("h1_pred_Y_noConstraint_clone");
  h1_pred_Y_noConstraint->Draw("same e2");
  h1_pred_Y_noConstraint->SetMarkerStyle(1);
  h1_pred_Y_noConstraint->SetMarkerColor(kRed);
  h1_pred_Y_noConstraint->SetLineColor(kRed);
  h1_pred_Y_noConstraint->SetLineWidth(2);
  h1_pred_Y_noConstraint->SetFillColor(kRed);
  h1_pred_Y_noConstraint->SetFillStyle(3005);
  h1_pred_Y_noConstraint_clone->Draw("same hist");
  h1_pred_Y_noConstraint_clone->SetLineColor(kRed);
  h1_pred_Y_noConstraint_clone->SetLineWidth(2);  
  h1_pred_Y_noConstraint_clone->Draw("same axis");

  ///////
  h1_pred_Y_wiConstraint->Draw("same e2"); 
  TH1D *h1_pred_Y_wiConstraint_clone = (TH1D*)h1_pred_Y_wiConstraint->Clone("h1_pred_Y_wiConstraint_clone");
  h1_pred_Y_wiConstraint->Draw("same e2");
  h1_pred_Y_wiConstraint->SetMarkerStyle(1);
  h1_pred_Y_wiConstraint->SetMarkerColor(kBlue);
  h1_pred_Y_wiConstraint->SetLineColor(kBlue);
  h1_pred_Y_wiConstraint->SetLineWidth(2);
  h1_pred_Y_wiConstraint->SetFillColor(kBlue);
  h1_pred_Y_wiConstraint->SetFillStyle(3004);
  h1_pred_Y_wiConstraint_clone->Draw("same hist");
  h1_pred_Y_wiConstraint_clone->SetLineColor(kBlue);
  h1_pred_Y_wiConstraint_clone->SetLineWidth(2);  
  h1_pred_Y_wiConstraint_clone->Draw("same axis");

  
  TF1 *f1_top = new TF1("f1_top", "0", 0, 1e6);
  f1_top->Draw("same");
  f1_top->SetLineColor(kBlack);
  f1_top->SetLineStyle(7);
  
  gh_data->Draw("same pe");
  gh_data->SetMarkerStyle(20);
  gh_data->SetMarkerSize(1);
  gh_data->SetMarkerColor(kBlack);
  gh_data->SetLineColor(kBlack);

  TLegend *lg_user = new TLegend(0.5, 0.65, 0.85, 0.85);
  lg_user->AddEntry(gh_data, "Data", "lep");
  lg_user->AddEntry(h1_pred_Y_noConstraint, "Pred no constraint", "lf");
  lg_user->AddEntry(h1_pred_Y_wiConstraint, "Pred wi constraint", "lf");  
  lg_user->Draw();
  lg_user->SetBorderSize(0);
  lg_user->SetFillStyle(0);
  lg_user->SetTextSize(0.065);
  // //lg_user->SetTextFont(132);

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
  pt->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d, p-value: %5.3f", val_chi2_noConstraint,num_Y, p_value_noConstraint));
  ((TText*)pt->GetListOfLines()->Last())->SetTextColor(kRed);
  pt->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d, p-value: %5.3f", val_chi2_wiConstraint,num_Y, p_value_wiConstraint));
  ((TText*)pt->GetListOfLines()->Last())->SetTextColor(kBlue);
  pt->Draw();

  /////////////
  canv_spectra_cmp->cd();
  
  TPad *pad_bot = new TPad("pad_bot", "pad_bot", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot, 0.15, 0.1, 0.05, 0.3);
  pad_bot->Draw();
  pad_bot->cd();
  
  h1_ratio_basic->Draw();
  h1_ratio_basic->SetMinimum(0);
  h1_ratio_basic->SetMaximum(2);
  //h1_ratio_basic->SetLineColor(kRed);
  //h1_ratio_basic->SetLineStyle(7);
  func_title_size(h1_ratio_basic, 0.078, 0.078, 0.078, 0.078);
  func_xy_title(h1_ratio_basic, "Index", "Data / Pred");
  h1_ratio_basic->GetXaxis()->SetTickLength(0.05);
  h1_ratio_basic->GetXaxis()->CenterTitle();
  h1_ratio_basic->GetYaxis()->CenterTitle(); 
  h1_ratio_basic->GetYaxis()->SetTitleOffset(0.92);
  h1_ratio_basic->GetYaxis()->SetNdivisions(509);
  h1_ratio_basic->Draw("same axis");

  TF1 *f1_const = new TF1("f1_const", "1", 0, 1e6);
  f1_const->Draw("same");
  f1_const->SetLineColor(kBlack);
  f1_const->SetLineStyle(7);
  
  h1_pred_Y_noConstraint_rel_error->Draw("same e2");
  h1_pred_Y_noConstraint_rel_error->SetMarkerStyle(1);
  h1_pred_Y_noConstraint_rel_error->SetMarkerColor(kRed);
  h1_pred_Y_noConstraint_rel_error->SetLineColor(kRed);
  h1_pred_Y_noConstraint_rel_error->SetLineWidth(2);
  h1_pred_Y_noConstraint_rel_error->SetFillColor(kRed);
  h1_pred_Y_noConstraint_rel_error->SetFillStyle(3005);
  
  h1_pred_Y_wiConstraint_rel_error->Draw("same e2");
  h1_pred_Y_wiConstraint_rel_error->SetMarkerStyle(1);
  h1_pred_Y_wiConstraint_rel_error->SetMarkerColor(kBlue);
  h1_pred_Y_wiConstraint_rel_error->SetLineColor(kBlue);
  h1_pred_Y_wiConstraint_rel_error->SetLineWidth(2);
  h1_pred_Y_wiConstraint_rel_error->SetFillColor(kBlue);
  h1_pred_Y_wiConstraint_rel_error->SetFillStyle(3004);
  
  gh_ratio_noConstraint->Draw("same pe");
  gh_ratio_noConstraint->SetMarkerStyle(20);
  gh_ratio_noConstraint->SetMarkerSize(1);
  gh_ratio_noConstraint->SetMarkerColor(kRed);
  gh_ratio_noConstraint->SetLineColor(kRed);

  gh_ratio_wiConstraint->Draw("same pe");
  gh_ratio_wiConstraint->SetMarkerStyle(20);
  gh_ratio_wiConstraint->SetMarkerSize(1);
  gh_ratio_wiConstraint->SetMarkerColor(kBlue);
  gh_ratio_wiConstraint->SetLineColor(kBlue);
  
  roostr = TString::Format("canv_spectra_cmp_%02d.png", index);
  canv_spectra_cmp->SaveAs(roostr);

  ///////////////////////////////////////////////////// sub pictures, noConstraint
  ///////////////////////////////////////////////////// sub pictures, noConstraint

  roostr = TString::Format("canv_spectrum_noConstraint_%02d", index);
  TCanvas *canv_spectrum_noConstraint = new TCanvas(roostr, roostr, 1000, 950);  
  
  TPad *pad_top_noConstaint = new TPad("pad_top_noConstaint", "pad_top_noConstaint", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_noConstaint, 0.15, 0.1, 0.1, 0.05);
  pad_top_noConstaint->Draw();
  pad_top_noConstaint->cd();
  h1_pred_Y_noConstraint->Draw("e2");
  h1_pred_Y_noConstraint->SetMinimum(0);
  h1_pred_Y_noConstraint_clone->Draw("same hist");
  f1_top->Draw("same");  
  gh_data->Draw("same pe");
  h1_pred_Y_noConstraint_clone->Draw("same axis");  
  TLegend *lg_user_noConstaint = new TLegend(0.5, 0.65, 0.85, 0.85);
  lg_user_noConstaint->AddEntry(gh_data, "Data", "lep");
  lg_user_noConstaint->AddEntry(h1_pred_Y_noConstraint, "Pred no constraint", "lf");
  lg_user_noConstaint->Draw();
  lg_user_noConstaint->SetBorderSize(0);
  lg_user_noConstaint->SetFillStyle(0);
  lg_user_noConstaint->SetTextSize(0.065);
  // //lg_user_noConstaint->SetTextFont(132);
  pad_top_noConstaint->cd();
  pad_top_noConstaint->Update();
  x1 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.47;
  y1 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.45;  
  x2 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.85;
  y2 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.65;
  TPaveText *pt_noConstraint = new TPaveText( x1,y1,x2,y2,"l");
  pt_noConstraint->SetTextSize(0.065);
  pt_noConstraint->SetTextFont(42);
  pt_noConstraint->SetTextAlign(11);
  pt_noConstraint->SetBorderSize(0);
  pt_noConstraint->SetFillStyle(0);
  pt_noConstraint->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d, p-value: %5.3f", val_chi2_noConstraint,num_Y, p_value_noConstraint));
  ((TText*)pt_noConstraint->GetListOfLines()->Last())->SetTextColor(kRed);
  pt_noConstraint->Draw();
  
  canv_spectrum_noConstraint->cd();  
  TPad *pad_bot_noConstraint = new TPad("pad_bot_noConstraint", "pad_bot_noConstraint", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_noConstraint, 0.15, 0.1, 0.05, 0.3);
  pad_bot_noConstraint->Draw();
  pad_bot_noConstraint->cd();
  h1_ratio_basic->Draw();
  f1_const->Draw("same");
  h1_pred_Y_noConstraint_rel_error->Draw("same e2");
  gh_ratio_noConstraint->Draw("same pe");
  h1_ratio_basic->Draw("same axis");
  
  roostr = TString::Format("canv_spectrum_noConstraint_%02d.png", index);
  canv_spectrum_noConstraint->SaveAs(roostr);

  ///////////////////////////////////////////////////// sub pictures, wiConstraint
  ///////////////////////////////////////////////////// sub pictures, wiConstraint

  roostr = TString::Format("canv_spectrum_wiConstraint_%02d", index);
  TCanvas *canv_spectrum_wiConstraint = new TCanvas(roostr, roostr, 1000, 950);  
  
  TPad *pad_top_wiConstaint = new TPad("pad_top_wiConstaint", "pad_top_wiConstaint", 0, 0.45, 1, 1);
  func_canv_margin(pad_top_wiConstaint, 0.15, 0.1, 0.1, 0.05);
  pad_top_wiConstaint->Draw();
  pad_top_wiConstaint->cd();
  h1_pred_Y_wiConstraint->Draw("e2");
  h1_pred_Y_wiConstraint->SetMinimum(0);
  h1_pred_Y_wiConstraint->SetTitle("");
  h1_pred_Y_wiConstraint_clone->Draw("same hist");
  f1_top->Draw("same");  
  gh_data->Draw("same pe");
  h1_pred_Y_wiConstraint_clone->Draw("same axis");  
  TLegend *lg_user_wiConstaint = new TLegend(0.5, 0.65, 0.85, 0.85);
  lg_user_wiConstaint->AddEntry(gh_data, "Data", "lep");
  lg_user_wiConstaint->AddEntry(h1_pred_Y_wiConstraint, "Pred wi constraint", "lf");
  lg_user_wiConstaint->Draw();
  lg_user_wiConstaint->SetBorderSize(0);
  lg_user_wiConstaint->SetFillStyle(0);
  lg_user_wiConstaint->SetTextSize(0.065);
  // //lg_user_wiConstaint->SetTextFont(132);
  pad_top_wiConstaint->cd();
  pad_top_wiConstaint->Update();
  x1 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.47;
  y1 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.45;  
  x2 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.85;
  y2 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.65;
  TPaveText *pt_wiConstraint = new TPaveText( x1,y1,x2,y2,"l");
  pt_wiConstraint->SetTextSize(0.065);
  pt_wiConstraint->SetTextFont(42);
  pt_wiConstraint->SetTextAlign(11);
  pt_wiConstraint->SetBorderSize(0);
  pt_wiConstraint->SetFillStyle(0);
  pt_wiConstraint->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d, p-value: %5.3f", val_chi2_wiConstraint,num_Y, p_value_wiConstraint));
  ((TText*)pt_wiConstraint->GetListOfLines()->Last())->SetTextColor(kBlue);
  pt_wiConstraint->Draw();
  
  canv_spectrum_wiConstraint->cd();  
  TPad *pad_bot_wiConstraint = new TPad("pad_bot_wiConstraint", "pad_bot_wiConstraint", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot_wiConstraint, 0.15, 0.1, 0.05, 0.3);
  pad_bot_wiConstraint->Draw();
  pad_bot_wiConstraint->cd();
  h1_ratio_basic->Draw();
  f1_const->Draw("same");
  h1_pred_Y_wiConstraint_rel_error->Draw("same e2");
  gh_ratio_wiConstraint->Draw("same pe");
  h1_ratio_basic->Draw("same axis");
  
  roostr = TString::Format("canv_spectrum_wiConstraint_%02d.png", index);
  canv_spectrum_wiConstraint->SaveAs(roostr);
  
  ///////////////////////////////////////////////////// 
  ///////////////////////////////////////////////////// 

  //h1_pred_Y_wiConstraint->

  int bins_h1_pred_wi2no = h1_pred_Y_wiConstraint->GetNbinsX();
  roostr = TString::Format("h1_pred_wi2no_%02d", index);
  TH1D *h1_pred_wi2no = new TH1D(roostr, "", bins_h1_pred_wi2no, 0, bins_h1_pred_wi2no);
  for(int ibin=1; ibin<=bins_h1_pred_wi2no; ibin++) {
    double val_wi = h1_pred_Y_wiConstraint->GetBinContent(ibin);
    double val_no = h1_pred_Y_noConstraint->GetBinContent(ibin);
    double val_wi2no = val_wi/val_no;
    h1_pred_wi2no->SetBinContent(ibin, val_wi2no);
  }

  roostr = TString::Format("canv_h1_pred_wi2no_%02d", index);
  TCanvas *canv_h1_pred_wi2no = new TCanvas(roostr, roostr, 900, 650);
  func_canv_margin(canv_h1_pred_wi2no, 0.15, 0.1, 0.1, 0.15);
  canv_h1_pred_wi2no->SetGridy();
  h1_pred_wi2no->Draw("hist");
  h1_pred_wi2no->SetMinimum(0.5);
  h1_pred_wi2no->SetMaximum(1.5);
  h1_pred_wi2no->SetStats(0);
  h1_pred_wi2no->SetLineColor(kBlue);
  func_title_size(h1_pred_wi2no, 0.05, 0.05, 0.05, 0.05);
  h1_pred_wi2no->SetXTitle("Index");
  h1_pred_wi2no->GetXaxis()->CenterTitle();
  h1_pred_wi2no->GetXaxis()->SetTitleOffset(1.2);
  h1_pred_wi2no->SetYTitle("Pred wi/wo");
  h1_pred_wi2no->GetYaxis()->CenterTitle();
  h1_pred_wi2no->GetYaxis()->SetNdivisions(508);
  h1_pred_wi2no->Draw("same axis");
  roostr = TString::Format("canv_h1_pred_wi2no_%02d.png", index);
  canv_h1_pred_wi2no->SaveAs(roostr);
}

///////////////////////////// ccc

void TLee::Exe_Goodness_of_Fit(map<int, double>map_h1_pred, map<int, double>map_h1_data, TMatrixD matrix_syst, int index)
{
  cout<<endl<<" ---> processing Goodness of Fit (GOF)"<<endl<<endl;
  TString roostr = "";

  roostr = TString::Format("h1_pred_%02d", index);
  int size_pred = map_h1_pred.size();
  TH1D *h1_pred = new TH1D(roostr, roostr, size_pred, 0, size_pred);
  for(int ibin=1; ibin<=size_pred; ibin++) h1_pred->SetBinContent(ibin, map_h1_pred[ibin]);
   
  roostr = TString::Format("h1_data_%02d", index);
  int size_data = map_h1_data.size();
  TH1D *h1_data = new TH1D(roostr, roostr, size_data, 0, size_data);
  for(int ibin=1; ibin<=size_data; ibin++) h1_data->SetBinContent(ibin, map_h1_data[ibin]);
     
  for(int ibin=1; ibin<=h1_pred->GetNbinsX(); ibin++) {
    h1_pred->SetBinError( ibin, sqrt(matrix_syst(ibin-1, ibin-1)) );
  }

  ////////
  int rows = matrix_syst.GetNrows();
  
  TMatrixD matrix_pred(1, rows);
  TMatrixD matrix_data(1, rows);
  
  for(int ibin=1; ibin<=rows; ibin++) {
    double val_pred = h1_pred->GetBinContent( ibin );
    matrix_pred(0, ibin-1) = val_pred;
    double val_data = h1_data->GetBinContent( ibin );
    matrix_data(0, ibin-1) = val_data;    
  }

  TMatrixD matrix_delta = matrix_pred - matrix_data;
  TMatrixD matrix_delta_T( rows, 1 );
  matrix_delta_T.Transpose( matrix_delta );
  
  TMatrixD matrix_stat(rows, rows);
  for(int i=0; i<rows; i++) {
    double val_pred = matrix_pred(0, i);
    double val_data = matrix_data(0, i);
    
    matrix_stat(i,i) = val_pred;

    if( matrix_data(0, i)==1 ) {
      if( matrix_pred(0, i)<0.461 ) {// DocDB-32520, when the prediction is sufficiently low
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

  ///////////////

  TH1D *h1_ratio_basic = (TH1D*)h1_pred->Clone("h1_ratio_basic");
  h1_ratio_basic->Reset();
  h1_ratio_basic->SetTitle("");

  TGraphAsymmErrors *gh_ratio = new TGraphAsymmErrors();
  TGraphAsymmErrors *gh_data = new TGraphAsymmErrors();
  
  for(int ibin=1; ibin<=h1_ratio_basic->GetNbinsX(); ibin++) {    
    double val_data = h1_data->GetBinContent(ibin);
    double val_pred = h1_pred->GetBinContent(ibin);
    double val_ratio = val_data/val_pred;

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
    double val_x = h1_ratio_basic->GetBinCenter(ibin);
    double val_halfw = h1_ratio_basic->GetBinWidth(ibin)/2;
    gh_ratio->SetPoint( n_point, val_x, val_ratio );
    gh_ratio->SetPointError( n_point, val_halfw, val_halfw, val_ratio_low, val_ratio_hgh );

    gh_data->SetPoint( n_point, val_x, val_data );
    gh_data->SetPointError( n_point, val_halfw, val_halfw, val_data-val_data_low, val_data_hgh-val_data );
  }

  ////////////////
  
  TH1D *h1_pred_rel_error = (TH1D*)h1_pred->Clone("h1_pred_rel_error");
  h1_pred_rel_error->Reset();
  h1_pred_rel_error->SetTitle("");
  for(int ibin=1; ibin<=h1_pred_rel_error->GetNbinsX(); ibin++) {
    double val_err = h1_pred->GetBinError(ibin);
    double val_cv = h1_pred->GetBinContent(ibin);
    double rel_err = val_err/val_cv;
    if( val_cv==0 ) rel_err = 0;
    h1_pred_rel_error->SetBinContent(ibin, 1);
    h1_pred_rel_error->SetBinError(ibin, rel_err);
  }

  ///////////////////////////////////////////////////////////////////////// Plotting
  
  double ymax_pred = h1_pred->GetMaximum();
  double ymax_data = h1_data->GetMaximum();
  double ymax_user = (ymax_pred>ymax_data) ? ymax_pred : ymax_data;
  
  roostr = TString::Format("canv_spectra_cmp_%02d", index);
  TCanvas *canv_spectra_cmp = new TCanvas(roostr, roostr, 1000, 950);  
  //func_canv_margin(canv_spectra_cmp, 0.15, 0.1,0.1,0.15);

  /////////////
  TPad *pad_top = new TPad("pad_top", "pad_top", 0, 0.45, 1, 1);
  func_canv_margin(pad_top, 0.15, 0.1, 0.1, 0.05);
  pad_top->Draw();
  pad_top->cd();
  
  h1_pred->Draw("e2");
  //h1_pred->SetMaximum(ymax_user*1.3);
  //h1_pred->SetMinimum(1e-3);
  h1_pred->SetTitle("");
  func_title_size(h1_pred, 0.065, 0.065, 0.065, 0.065);
  func_xy_title(h1_pred, "Index", "Entries");
  h1_pred->GetXaxis()->SetNdivisions(506);
  h1_pred->GetYaxis()->SetNdivisions(506);
  h1_pred->GetYaxis()->CenterTitle();
  h1_pred->GetXaxis()->CenterTitle();

  h1_pred->GetXaxis()->SetLabelColor(10);
  h1_pred->GetXaxis()->SetTitleColor(10);
  
  TH1D *h1_pred_clone = (TH1D*)h1_pred->Clone("h1_pred_clone");

  TF1 *f1_top = new TF1("f1_top", "0", 0, 1e6);
  f1_top->Draw("same");
  f1_top->SetLineColor(kBlack);
  f1_top->SetLineStyle(7);
  
  h1_pred->Draw("same e2");
  h1_pred->SetMarkerStyle(1);
  h1_pred->SetMarkerColor(kRed);
  h1_pred->SetLineColor(kRed);
  h1_pred->SetLineWidth(2);
  h1_pred->SetFillColor(kRed);
  h1_pred->SetFillStyle(3005);

  h1_pred_clone->Draw("same hist");
  h1_pred_clone->SetLineColor(kRed);
  h1_pred_clone->SetLineWidth(2);  
  h1_pred_clone->Draw("same axis");

  gh_data->Draw("same pe");
  gh_data->SetMarkerStyle(20);
  gh_data->SetMarkerSize(1);
  gh_data->SetMarkerColor(kBlue);
  gh_data->SetLineColor(kBlue);

  TLegend *lg_user = new TLegend(0.6, 0.7, 0.85, 0.85);
  lg_user->AddEntry(gh_data, "Data", "lep");
  lg_user->AddEntry(h1_pred, "Prediction", "lf");
  lg_user->Draw();
  lg_user->SetBorderSize(0);
  lg_user->SetFillStyle(0);
  lg_user->SetTextSize(0.065);
  //lg_user->SetTextFont(132);

  pad_top->cd();
  pad_top->Update();
  double x1 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.6;
  double y1 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.5;  
  double x2 = gPad->GetUxmin() + ( gPad->GetUxmax()-gPad->GetUxmin() ) * 0.85;
  double y2 = gPad->GetUymin() + ( gPad->GetUymax()-gPad->GetUymin() ) * 0.7;
  TPaveText *pt = new TPaveText( x1,y1,x2,y2,"l");
  pt->SetTextSize(0.065);
  pt->SetTextFont(42);
  pt->SetTextAlign(11);
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->AddText(TString::Format("#chi^{2}/ndf: %4.2f/%d", val_chi2,rows));
  pt->AddText(TString::Format("P-value: %5.3f", p_value));
  pt->Draw();
  // roostr = "canv_matrix";
  // TCanvas *canv_matrix = new TCanvas(roostr, roostr, 900, 650);  
  // func_canv_margin(canv_matrix, 0.15, 0.2,0.1,0.15);
  // matrix_syst.Draw("colz");


  /////////////
  canv_spectra_cmp->cd();
  
  TPad *pad_bot = new TPad("pad_bot", "pad_bot", 0, 0, 1, 0.45);
  func_canv_margin(pad_bot, 0.15, 0.1, 0.05, 0.3);
  pad_bot->Draw();
  pad_bot->cd();
  
  h1_ratio_basic->Draw();
  h1_ratio_basic->SetMinimum(0);
  h1_ratio_basic->SetMaximum(2);
  //h1_ratio_basic->SetLineColor(kRed);
  //h1_ratio_basic->SetLineStyle(7);
  func_title_size(h1_ratio_basic, 0.078, 0.078, 0.078, 0.078);
  func_xy_title(h1_ratio_basic, "Index", "Data / Pred");
  h1_ratio_basic->GetXaxis()->SetTickLength(0.05);
  h1_ratio_basic->GetXaxis()->CenterTitle();
  h1_ratio_basic->GetYaxis()->CenterTitle(); 
  h1_ratio_basic->GetYaxis()->SetTitleOffset(0.92);
  h1_ratio_basic->GetYaxis()->SetNdivisions(509);
  h1_ratio_basic->Draw("same axis");

  TF1 *f1_const = new TF1("f1_const", "1", 0, 1e6);
  f1_const->Draw("same");
  f1_const->SetLineColor(kBlack);
  f1_const->SetLineStyle(7);
  
  h1_pred_rel_error->Draw("same e2");
  h1_pred_rel_error->SetMarkerStyle(1);
  h1_pred_rel_error->SetMarkerColor(kRed);
  h1_pred_rel_error->SetLineColor(kRed);
  h1_pred_rel_error->SetLineWidth(2);
  h1_pred_rel_error->SetFillColor(kRed);
  h1_pred_rel_error->SetFillStyle(3005);
  
  gh_ratio->Draw("same pe");
  gh_ratio->SetMarkerStyle(20);
  gh_ratio->SetMarkerSize(1);
  gh_ratio->SetMarkerColor(kBlue);
  gh_ratio->SetLineColor(kBlue);

  roostr = TString::Format("canv_spectra_cmp_%02d.png", index);
  canv_spectra_cmp->SaveAs(roostr);
}

///////////////////////////// ccc

void TLee::Set_Collapse()
{
  ////////
  TMatrixD matrix_pred_old_world(1, num_global_elements);
  for(int ibin=1; ibin<=num_global_elements; ibin++) {
    matrix_pred_old_world(0, ibin-1) = map_input_spectrum_global_bin[ibin];
  }

  ////////
  TMatrixD matrix_trans_wiLee_scaleF_Lee = matrix_trans_wiLee;
  int rows = matrix_trans_wiLee_scaleF_Lee.GetNrows();
  int cols = matrix_trans_wiLee_scaleF_Lee.GetNcols();
  for(int ibin=0; ibin<rows; ibin++) {
    for(int jbin=0; jbin<cols; jbin++) {
      if( map_Lee_old_world[ibin]==1 ) {
	matrix_trans_wiLee_scaleF_Lee(ibin, jbin) *= scaleF_Lee;
      }
    }
  }
  TMatrixD matrix_trans_wiLee_scaleF_Lee_T(cols, rows);
  matrix_trans_wiLee_scaleF_Lee_T.Transpose(matrix_trans_wiLee_scaleF_Lee);

  ///////////////////////////////////////////// wiLee

  TMatrixD matrix_pred_wiLee = matrix_pred_old_world * matrix_trans_wiLee_scaleF_Lee;

  matrix_input_cov_total_scaled.Clear();
  matrix_input_cov_total_scaled.ResizeTo(num_global_elements, num_global_elements);

  /// ttt  
  // TMatrixD matrix_input_cov_flux_Xs_scaled;
  // TMatrixD matrix_input_cov_detector_scaled;
  // TMatrixD matrix_input_cov_additional_scaled;
  // TMatrixD matrix_input_cov_total_scaled;

  // TMatrixD matrix_input_cov_flux_scaled;
  // TMatrixD matrix_input_cov_Xs_scaled;

  matrix_input_cov_flux_Xs_scaled.Clear();
  matrix_input_cov_flux_Xs_scaled.ResizeTo(num_global_elements, num_global_elements);  
  matrix_input_cov_detector_scaled.Clear();
  matrix_input_cov_detector_scaled.ResizeTo(num_global_elements, num_global_elements);  
  matrix_input_cov_additional_scaled.Clear();
  matrix_input_cov_additional_scaled.ResizeTo(num_global_elements, num_global_elements);
  
  matrix_input_cov_flux_scaled.Clear();
  matrix_input_cov_flux_scaled.ResizeTo(num_global_elements, num_global_elements);  
  matrix_input_cov_Xs_scaled.Clear();
  matrix_input_cov_Xs_scaled.ResizeTo(num_global_elements, num_global_elements);
  

  if( flag_syst_flux_Xs ) {
    matrix_input_cov_total_scaled += matrix_input_cov_flux_Xs;

    matrix_input_cov_flux_Xs_scaled = matrix_input_cov_flux_Xs;
    matrix_input_cov_flux_scaled = matrix_input_cov_flux;
    matrix_input_cov_Xs_scaled = matrix_input_cov_Xs;
  }  
  if( flag_syst_detector ) {
    matrix_input_cov_total_scaled += matrix_input_cov_detector;

    matrix_input_cov_detector_scaled = matrix_input_cov_detector;
  }
  if( flag_syst_additional ) {
    matrix_input_cov_total_scaled += matrix_input_cov_additional;

    matrix_input_cov_additional_scaled = matrix_input_cov_additional;
  }
  
  matrix_collapse_absolute_covariance_wiLee.Clear();
  matrix_collapse_absolute_covariance_wiLee.ResizeTo(cols, cols);
  matrix_collapse_absolute_covariance_wiLee = matrix_trans_wiLee_scaleF_Lee_T * matrix_input_cov_total_scaled * matrix_trans_wiLee_scaleF_Lee;

  // TMatrixD matrix_collapse_absolute_covariance_wiLee_flux;
  // TMatrixD matrix_collapse_absolute_covariance_wiLee_Xs;
  // TMatrixD matrix_collapse_absolute_covariance_wiLee_detector;
  // TMatrixD matrix_collapse_absolute_covariance_wiLee_additional;
  // TMatrixD matrix_collapse_absolute_covariance_wiLee_mc_stat;
  
  matrix_collapse_absolute_covariance_wiLee_flux.Clear();
  matrix_collapse_absolute_covariance_wiLee_flux.ResizeTo(cols, cols);
  matrix_collapse_absolute_covariance_wiLee_flux = matrix_trans_wiLee_scaleF_Lee_T * matrix_input_cov_flux_scaled * matrix_trans_wiLee_scaleF_Lee;
  
  matrix_collapse_absolute_covariance_wiLee_Xs.Clear();
  matrix_collapse_absolute_covariance_wiLee_Xs.ResizeTo(cols, cols);
  matrix_collapse_absolute_covariance_wiLee_Xs = matrix_trans_wiLee_scaleF_Lee_T * matrix_input_cov_Xs_scaled * matrix_trans_wiLee_scaleF_Lee;
  
  matrix_collapse_absolute_covariance_wiLee_detector.Clear();
  matrix_collapse_absolute_covariance_wiLee_detector.ResizeTo(cols, cols);
  matrix_collapse_absolute_covariance_wiLee_detector = matrix_trans_wiLee_scaleF_Lee_T * matrix_input_cov_detector_scaled * matrix_trans_wiLee_scaleF_Lee;
  
  matrix_collapse_absolute_covariance_wiLee_additional.Clear();
  matrix_collapse_absolute_covariance_wiLee_additional.ResizeTo(cols, cols);
  matrix_collapse_absolute_covariance_wiLee_additional = matrix_trans_wiLee_scaleF_Lee_T * matrix_input_cov_additional_scaled * matrix_trans_wiLee_scaleF_Lee;
    
  matrix_collapse_absolute_covariance_wiLee_mc_stat.Clear();
  matrix_collapse_absolute_covariance_wiLee_mc_stat.ResizeTo(cols, cols);
  
  if( flag_syst_mc_stat ) {
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      double val_mc_stat = gh_mc_stat_bin[ibin]->Eval( scaleF_Lee );
      matrix_collapse_absolute_covariance_wiLee(ibin-1, ibin-1) += val_mc_stat;

      matrix_collapse_absolute_covariance_wiLee_mc_stat(ibin-1, ibin-1) = val_mc_stat;
    }
  }

  for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
    map_collapsed_prediction_wiLee[ibin] = matrix_pred_wiLee(0, ibin-1);
  }

  ///////////////////////////////////////////// noLee
  /*
  TMatrixD matrix_trans_noLee_scaleF_Lee = matrix_trans_noLee;
  TMatrixD matrix_trans_noLee_scaleF_Lee_T(cols, rows);
  
  matrix_trans_noLee_scaleF_Lee_T.Transpose(matrix_trans_noLee_scaleF_Lee);

  TMatrixD matrix_pred_noLee = matrix_pred_old_world * matrix_trans_noLee_scaleF_Lee;

  matrix_input_cov_total_scaled.Clear();
  matrix_input_cov_total_scaled.ResizeTo(num_global_elements, num_global_elements);
  if( flag_syst_flux_Xs ) matrix_input_cov_total_scaled += matrix_input_cov_flux_Xs;
  if( flag_syst_detector ) matrix_input_cov_total_scaled += matrix_input_cov_detector;
  if( flag_syst_additional ) matrix_input_cov_total_scaled += matrix_input_cov_additional;

  matrix_collapse_absolute_covariance_noLee.Clear();
  matrix_collapse_absolute_covariance_noLee.ResizeTo(cols, cols);
  matrix_collapse_absolute_covariance_noLee = matrix_trans_noLee_scaleF_Lee_T * matrix_input_cov_total_scaled * matrix_trans_noLee_scaleF_Lee;

  if( flag_syst_mc_stat ) {
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      double val_mc_stat = gh_mc_stat_bin[ibin]->Eval( scaleF_Lee );
      matrix_collapse_absolute_covariance_noLee(ibin-1, ibin-1) += val_mc_stat;
    }
  }

  for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
    map_collapsed_prediction_noLee[ibin] = matrix_pred_noLee(0, ibin-1);
  }
  */  
}

///////////////////////////// ccc

void TLee::Set_TransformMatrix()
{
  num_elements_collapse = 0;
    
  map<int, int>map_collapse_channel_bins;// new-world channels
  int line_map = 0;
  line_map++; map_collapse_channel_bins[line_map] = map_input_spectrum_ch_binnum[1];
  line_map++; map_collapse_channel_bins[line_map] = map_input_spectrum_ch_binnum[2];
  line_map++; map_collapse_channel_bins[line_map] = map_input_spectrum_ch_binnum[3];
  line_map++; map_collapse_channel_bins[line_map] = map_input_spectrum_ch_binnum[4];
  line_map++; map_collapse_channel_bins[line_map] = map_input_spectrum_ch_binnum[5];
  line_map++; map_collapse_channel_bins[line_map] = map_input_spectrum_ch_binnum[6];
  line_map++; map_collapse_channel_bins[line_map] = map_input_spectrum_ch_binnum[7];

  map_collapse_channel_local2global.clear();  
  int size_map = map_collapse_channel_bins.size();  
  for(int idx=1; idx<=size_map; idx++) {
    for(int ibin=1; ibin<=map_collapse_channel_bins[idx]; ibin++) {
      num_elements_collapse++;
      map_collapse_channel_local2global[idx][ibin] = num_elements_collapse;
    }    
  }

  ///////////////////////////////////

  matrix_trans_wiLee.Clear();
  matrix_trans_wiLee.ResizeTo( num_global_elements, num_elements_collapse );// old_world, new_world
  matrix_trans_noLee.Clear();
  matrix_trans_noLee.ResizeTo( num_global_elements, num_elements_collapse );// old_world, new_world

  map<int, int>map_wiLee_old2new;
  map_wiLee_old2new[1] = 1;
  map_wiLee_old2new[2] = 2;
  map_wiLee_old2new[3] = 3;
  map_wiLee_old2new[4] = 4;
  map_wiLee_old2new[5] = 5;
  map_wiLee_old2new[6] = 6;
  map_wiLee_old2new[7] = 7;
  map_wiLee_old2new[8] = 1;
  map_wiLee_old2new[9] = 2;
  map_wiLee_old2new[10]= 1;
  map_wiLee_old2new[11]= 2;
  map_wiLee_old2new[12]= 3;
  map_wiLee_old2new[13]= 4;
  map_wiLee_old2new[14]= 5;
  map_wiLee_old2new[15]= 6;
  map_wiLee_old2new[16]= 7;
  
  for(auto it=map_wiLee_old2new.begin(); it!=map_wiLee_old2new.end(); it++) {
    int old_world_channel = it->first;
    int new_world_channel = it->second;
    for(int ibin=1; ibin<=map_input_spectrum_ch_binnum[old_world_channel]; ibin++ ) {

      ///
      int old_world_global_bin = map_old_world_local2global[old_world_channel][ibin];
      int new_world_global_bin = map_collapse_channel_local2global[new_world_channel][ibin];
      
      ///
      matrix_trans_wiLee(old_world_global_bin-1, new_world_global_bin-1) = 1;
      
      ///
      if( map_Lee_ch.find(old_world_channel)==map_Lee_ch.end() ) {
        matrix_trans_noLee(old_world_global_bin-1, new_world_global_bin-1) = 1;
      }
      
    }// ibin
  }// it

  cout<<endl<<" ---> Finish TransformMatrix"<<endl<<endl;
}

///////////////////////////// ccc

void TLee::Set_POT_implement()
{
  double scaleF_POT2 = scaleF_POT * scaleF_POT;
  
  /////////////// map_input_spectrum_ch_bin_scaled
  
  for(auto it=map_input_spectrum_ch_bin.begin(); it!=map_input_spectrum_ch_bin.end(); it++) {
    int ch = it->first;
    int bins = map_input_spectrum_ch_bin[ch].size();
    
    for(int ibin=1; ibin<=bins; ibin++) {
      int global_bin = map_old_world_local2global[ch][ibin];
      double scaleF = scaleF_POT;
      map_input_spectrum_ch_bin[ch][ibin] = map_input_spectrum_ch_bin[ch][ibin] * scaleF;
      map_input_spectrum_global_bin[global_bin] = map_input_spectrum_global_bin[global_bin] * scaleF;      
    }// ibin
  }// ch

  /////////////// map_data_spectrum_ch_bin_scaled
  
  for(auto it=map_data_spectrum_ch_bin.begin(); it!=map_data_spectrum_ch_bin.end(); it++) {
    int ch = it->first;
    int bins = map_data_spectrum_ch_bin[ch].size();
    for(int ibin=1; ibin<=bins; ibin++) {
      double scaleF = scaleF_POT;
      map_data_spectrum_ch_bin_scaled[ch][ibin] = map_data_spectrum_ch_bin[ch][ibin] * scaleF;
    }
  }

  for(auto it=map_data_spectrum_global_bin.begin(); it!=map_data_spectrum_global_bin.end(); it++) {
    int ibin = it->first;
    double scaleF = scaleF_POT;
    map_data_spectrum_global_bin_scaled[ibin] = map_data_spectrum_global_bin[ibin] * scaleF;
  }

  //////////////// flux&Xs, detector, additional
  ////////////////
  
  int rows = matrix_input_cov_flux_Xs.GetNrows();
  num_global_elements = rows;

  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {      
      matrix_input_cov_flux_Xs(ibin-1, jbin-1) = matrix_input_cov_flux_Xs(ibin-1, jbin-1) * scaleF_POT2;
      matrix_input_cov_detector(ibin-1, jbin-1) = matrix_input_cov_detector(ibin-1, jbin-1) * scaleF_POT2;
      matrix_input_cov_additional(ibin-1, jbin-1) = matrix_input_cov_additional(ibin-1, jbin-1) * scaleF_POT2;

      matrix_input_cov_flux(ibin-1, jbin-1) = matrix_input_cov_flux(ibin-1, jbin-1) * scaleF_POT2;
      matrix_input_cov_Xs(ibin-1, jbin-1) = matrix_input_cov_Xs(ibin-1, jbin-1) * scaleF_POT2;
    }
  }

  //////////////// MC stat
  ////////////////

  int size_gh = gh_mc_stat_bin.size();
  for(int ibin=1; ibin<=size_gh; ibin++) {
    for(int idx=1; idx<=gh_mc_stat_bin[ibin]->GetN(); idx++) {
      double x(0), y(0);
      gh_mc_stat_bin[ibin]->GetPoint(idx-1, x, y);
      gh_mc_stat_bin[ibin]->SetPoint(idx-1, x, y*scaleF_POT2);
    }    
  }

  cout<<" ---> Finish POT_implement"<<endl<<endl;
}

///////////////////////////// ccc

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
  
  roostr = "./data_framework/merge_all.root";
  TFile *file_spectra = new TFile(roostr, "read");
  int size_map_input_spectrum_ch_str = map_input_spectrum_ch_str.size();
 
  for(int idx=1; idx<=size_map_input_spectrum_ch_str; idx++) {
    roostr = TString::Format("histo_%d", idx);
    TH1F *h1_spectrum = (TH1F*)file_spectra->Get(roostr);
    int bins = h1_spectrum->GetNbinsX() + 1;    
    cout<<Form(" %2d  %-20s   bin-num %2d   %-50s", idx, map_input_spectrum_ch_str[idx].Data(), bins, h1_spectrum->GetTitle())<<endl;
    
    map_input_spectrum_ch_binnum[idx] = bins;
    bins_total_old_world += bins;
 
    for(int ibin=1; ibin<=bins; ibin++) {
      map_input_spectrum_ch_bin[idx][ibin] = h1_spectrum->GetBinContent(ibin);
/*
      if( idx==1 || idx==8 ) {
	if( ibin==1 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.5;
	if( ibin==2 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.34;
	if( ibin==3 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.28;
	if( ibin==4 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.23;
	if( ibin==5 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.24;
	if( ibin==6 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.19;
	if( ibin==7 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.21;
	if( ibin==8 ) map_input_spectrum_ch_bin[idx][ibin] *= 1.12;
      }
*/      
    }

    map_input_h1_spectrum[idx] = (TH1F*)file_spectra->Get(roostr);
  }
  cout<<endl;

  int global_bin = 0;
  int size_map_input_spectrum_ch_bin = map_input_spectrum_ch_bin.size();
  for(int idx=1; idx<=size_map_input_spectrum_ch_bin; idx++) {
    for(int ibin=1; ibin<=map_input_spectrum_ch_binnum[idx]; ibin++ ) {
      global_bin++;
      map_input_spectrum_global_bin[global_bin] = map_input_spectrum_ch_bin[idx][ibin];
      
      map_old_world_local2global[idx][ibin] = global_bin;
      map_old_world_global2localch[global_bin] = idx;
      map_old_world_global2localbin[global_bin] = ibin;

      if( map_Lee_ch.find(idx)!=map_Lee_ch.end() ) {
	map_Lee_old_world[global_bin-1] = 1;
      }
      
    }
  }
  
  ////////////////////////////////////// data

  int line_data = 0;
  for(int idx=1; idx<=7; idx++) {
    roostr = TString::Format("hdata_obsch_%d", idx);
    TH1F *h1_spectrum = (TH1F*)file_spectra->Get(roostr);
    for(int ibin=1; ibin<=h1_spectrum->GetNbinsX()+1; ibin++) {
      map_data_spectrum_ch_bin[idx][ibin] = h1_spectrum->GetBinContent(ibin);

      line_data++;
      map_data_spectrum_global_bin[line_data] = h1_spectrum->GetBinContent(ibin);
    }

    map_data_h1_spectrum[idx] = (TH1F*)file_spectra->Get(roostr);
  }
  
  ////////////////////////////////////////// flux_Xs

  map<int, TFile*>map_file_flux_Xs_frac;  
  map<int, TMatrixD*>map_matrix_flux_Xs_frac;
  
  TMatrixD matrix_flux_Xs_frac(bins_total_old_world, bins_total_old_world);
  TMatrixD matrix_flux_frac(bins_total_old_world, bins_total_old_world);
  TMatrixD matrix_Xs_frac(bins_total_old_world, bins_total_old_world);

  for(int idx=1; idx<=17; idx++) {
    roostr = TString::Format("./data_framework/flux_Xs/cov_%d.root", idx);
    map_file_flux_Xs_frac[idx] = new TFile(roostr, "read");
    map_matrix_flux_Xs_frac[idx] = (TMatrixD*)map_file_flux_Xs_frac[idx]->Get(TString::Format("frac_cov_xf_mat_%d", idx));
    cout<<TString::Format(" ---> check: flux and Xs, %2d  ", idx)<<roostr<<endl;

    matrix_flux_Xs_frac += (*map_matrix_flux_Xs_frac[idx]);

    if( idx<=13 ) {// flux
      matrix_flux_frac += (*map_matrix_flux_Xs_frac[idx]);
    }
    else {// geant4 & GENIE
      matrix_Xs_frac += (*map_matrix_flux_Xs_frac[idx]);
    }
  }
  cout<<endl;

  ////////////////////////////////////////// detector
  
  map<int, TString>map_detectorfile_str;
  map_detectorfile_str[1] = "./data_framework/det/cov_LYDown.root";
  map_detectorfile_str[2] = "./data_framework/det/cov_LYRayleigh.root";
  map_detectorfile_str[3] = "./data_framework/det/cov_Recomb2.root";
  map_detectorfile_str[4] = "./data_framework/det/cov_SCE.root";
  map_detectorfile_str[5] = "./data_framework/det/cov_WMdEdx.root";
  map_detectorfile_str[6] = "./data_framework/det/cov_WMThetaXZ.root";
  map_detectorfile_str[7] = "./data_framework/det/cov_WMThetaYZ.root";
  map_detectorfile_str[8] = "./data_framework/det/cov_WMX.root";
  map_detectorfile_str[9] = "./data_framework/det/cov_WMYZ.root";
  
  map<int, TFile*>map_file_detector_frac;
  map<int, TMatrixD*>map_matrix_detector_frac;
  TMatrixD matrix_detector_frac(bins_total_old_world, bins_total_old_world);

  int size_map_detectorfile_str = map_detectorfile_str.size();
  for(int idx=1; idx<=size_map_detectorfile_str; idx++) {
    if(idx==5) continue;
    roostr = map_detectorfile_str[idx];
    map_file_detector_frac[idx] = new TFile(roostr, "read");
    map_matrix_detector_frac[idx] = (TMatrixD*)map_file_detector_frac[idx]->Get(TString::Format("frac_cov_det_mat_%d", idx));
    cout<<TString::Format(" ---> check: detector, %2d  ", idx)<<roostr<<endl;

    matrix_detector_frac += (*map_matrix_detector_frac[idx]);
  }
  cout<<endl;

  ////////////////////////////////////////// additional

  TMatrixD *matrix_additional_abs_point = (TMatrixD*)file_spectra->Get("cov_mat_add");
  TMatrixD matrix_additional_abs = (*matrix_additional_abs_point);

  //////////////////////////////////////////
  //////////////////////////////////////////

  int rows = bins_total_old_world;
  TMatrixD matrix_flux_Xs_abs(rows, rows);
  TMatrixD matrix_detector_abs(rows, rows);

  TMatrixD matrix_flux_abs(rows, rows);
  TMatrixD matrix_Xs_abs(rows, rows);
  
  for(int ibin=1; ibin<=rows; ibin++) {
    for(int jbin=1; jbin<=rows; jbin++) {
      double val_i = map_input_spectrum_global_bin[ibin];
      double val_j = map_input_spectrum_global_bin[jbin];
      double val_cov = 0;

      val_cov = matrix_flux_Xs_frac(ibin-1, jbin-1);
      matrix_flux_Xs_abs(ibin-1, jbin-1) = val_cov * val_i * val_j;      

      val_cov = matrix_detector_frac(ibin-1, jbin-1);
      matrix_detector_abs(ibin-1, jbin-1) = val_cov * val_i * val_j;

      
      val_cov = matrix_flux_frac(ibin-1, jbin-1);
      matrix_flux_abs(ibin-1, jbin-1) = val_cov * val_i * val_j;      

      val_cov = matrix_Xs_frac(ibin-1, jbin-1);
      matrix_Xs_abs(ibin-1, jbin-1) = val_cov * val_i * val_j;      

    }
  }

  matrix_input_cov_flux_Xs.Clear();
  matrix_input_cov_flux_Xs.ResizeTo(rows, rows);
  matrix_input_cov_detector.Clear();
  matrix_input_cov_detector.ResizeTo(rows, rows);
  matrix_input_cov_additional.Clear();
  matrix_input_cov_additional.ResizeTo(rows, rows);

  matrix_input_cov_flux.Clear();
  matrix_input_cov_flux.ResizeTo(rows, rows);
  matrix_input_cov_Xs.Clear();
  matrix_input_cov_Xs.ResizeTo(rows, rows);
  
  matrix_input_cov_flux_Xs = matrix_flux_Xs_abs;
  matrix_input_cov_detector = matrix_detector_abs;
  matrix_input_cov_additional = matrix_additional_abs;

  matrix_input_cov_flux = matrix_flux_abs;
  matrix_input_cov_Xs = matrix_Xs_abs;
  
  ////////////////////////////////////////// MC statistics

  map<int, map<int, double> >map_mc_stat_file_bin_Lee;
  map<int, map<int, double> >map_mc_stat_file_bin_mcStat;
  int gbins_mc_stat = 137;
  
  for(int ifile=0; ifile<=99; ifile++) {
    roostr = TString::Format("./data_framework/mc_stat/%d.log", ifile);
    ifstream InputFile_aa(roostr, ios::in);
    if(!InputFile_aa) { cerr<<" No input-list"<<endl; exit(1); }

    int line = 0;    
    double Lee = 1;
    double run = 1;
    
    for(int idx=1; idx<=gbins_mc_stat+1; idx++) {            
      int gbin = 0;
      int lbin = 0;
      double val_pred = 0;
      double mc_stat = 0;
      double nn_stat = 0;
      if(idx==1) {
	InputFile_aa>>Lee>>run;
      }
      else {
	InputFile_aa>>gbin>>lbin>>val_pred>>mc_stat>>nn_stat;
	line++;
	map_mc_stat_file_bin_Lee[ifile][line] = Lee;
	map_mc_stat_file_bin_mcStat[ifile][line] = mc_stat;
      }
    }
  }
  
  /// gh_mc_stat_bin
  for(int ibin=1; ibin<=gbins_mc_stat; ibin++) {
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


/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


//
// usage: root -l read_TLee.cc+(1, 1)
//

//int main(int argc, char** argv)
void read_TLee_v07(double scalePOT, int nfile)
{
  TString roostr = "";
  
  double scaleF_POT = 1;
  int ifile = 0;
  
  scaleF_POT = scalePOT;  
  ifile = nfile;
  /*
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
  */
  cout<<endl<<" ---> check, scaleF_POT "<<scaleF_POT<<", ifile "<<ifile<<endl<<endl;
  
  ////////////////////////////////////// Draw style
  
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
  
  TLee *Lee_test = new TLee();
 
  Lee_test->Set_Spectra_MatrixCov();  // these functions should be done only one time after the "TLee" definition
  Lee_test->scaleF_POT = scaleF_POT;

  
  //Lee_test->scaleF_POT = 69./5.327231;
  //Lee_test->scaleF_POT = 121./5.327231;
  
  
  Lee_test->Set_POT_implement();
  Lee_test->Set_TransformMatrix();
  
  ////////////////////////////////////////////////////////////////////////////////// Goodness of fit, fg
    
  Lee_test->scaleF_Lee           = 0;
  
  Lee_test->flag_syst_flux_Xs    = 1;
  Lee_test->flag_syst_detector   = 1;
  Lee_test->flag_syst_additional = 1;
  Lee_test->flag_syst_mc_stat    = 1;
  Lee_test->Set_Collapse();

  //////////////////////
  if( 0 )
    {

      // TMatrixD matrix_collapse_absolute_covariance_wiLee_flux;
      // TMatrixD matrix_collapse_absolute_covariance_wiLee_Xs;
      // TMatrixD matrix_collapse_absolute_covariance_wiLee_detector;
      // TMatrixD matrix_collapse_absolute_covariance_wiLee_additional;
      // TMatrixD matrix_collapse_absolute_covariance_wiLee_mc_stat;

      // TMatrixD matrix_collapse_absolute_covariance_wiLee;// -------------------->
      // map<int, double>map_collapsed_prediction_wiLee;// --------------------> index from 1
    
      int rows = Lee_test->num_elements_collapse;
    
      TMatrixD matrix_flux_correlation(rows, rows);
      TMatrixD matrix_Xs_correlation(rows, rows);
      TMatrixD matrix_detector_correlation(rows, rows);
      TMatrixD matrix_additional_correlation(rows, rows);
      TMatrixD matrix_mc_stat_correlation(rows, rows);
      TMatrixD matrix_total_correlation(rows, rows);

      TMatrixD matrix_flux_frac(rows, rows);
      TMatrixD matrix_Xs_frac(rows, rows);
      TMatrixD matrix_detector_frac(rows, rows);
      TMatrixD matrix_additional_frac(rows, rows);
      TMatrixD matrix_mc_stat_frac(rows, rows);
      TMatrixD matrix_total_frac(rows, rows);

    
      TH1D *h1_flux_frac = new TH1D("h1_flux_frac", "", rows, 0, rows);
      TH1D *h1_Xs_frac = new TH1D("h1_Xs_frac", "", rows, 0, rows);
      TH1D *h1_detector_frac = new TH1D("h1_detector_frac", "", rows, 0, rows);
      TH1D *h1_additional_frac = new TH1D("h1_additional_frac", "", rows, 0, rows);
      TH1D *h1_mc_stat_frac = new TH1D("h1_mc_stat_frac", "", rows, 0, rows);
      TH1D *h1_total_frac = new TH1D("h1_total_frac", "", rows, 0, rows);

      TH1D *h1_flux_frac2 = new TH1D("h1_flux_frac2", "", rows, 0, rows);
      TH1D *h1_Xs_frac2 = new TH1D("h1_Xs_frac2", "", rows, 0, rows);
      TH1D *h1_detector_frac2 = new TH1D("h1_detector_frac2", "", rows, 0, rows);
      TH1D *h1_additional_frac2 = new TH1D("h1_additional_frac2", "", rows, 0, rows);
      TH1D *h1_mc_stat_frac2 = new TH1D("h1_mc_stat_frac2", "", rows, 0, rows);
      TH1D *h1_total_frac2 = new TH1D("h1_total_frac2", "", rows, 0, rows);

      TH1D *h1_CV = new TH1D("h1_CV", "", rows, 0, rows);
    
      for(int ibin=1; ibin<=Lee_test->num_elements_collapse; ibin++) {    
	double val_flux       = Lee_test->matrix_collapse_absolute_covariance_wiLee_flux(ibin-1, ibin-1);
	double val_Xs         = Lee_test->matrix_collapse_absolute_covariance_wiLee_Xs(ibin-1, ibin-1);
	double val_detector   = Lee_test->matrix_collapse_absolute_covariance_wiLee_detector(ibin-1, ibin-1);
	double val_additional = Lee_test->matrix_collapse_absolute_covariance_wiLee_additional(ibin-1, ibin-1);
	double val_mc_stat    = Lee_test->matrix_collapse_absolute_covariance_wiLee_mc_stat(ibin-1, ibin-1);
	double val_total      = Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin-1, ibin-1);

	double val_check = val_flux + val_Xs + val_detector + val_additional + val_mc_stat;
	double val_CV = Lee_test->map_collapsed_prediction_wiLee[ibin];
	//cout<<TString::Format(" ---> check %3d, CV %9.3f, check/total %18.4f %18.4f, diff %9.6f", ibin, val_CV, val_check, val_total, val_check-val_total)<<endl;

      }


      for(int ibin=1; ibin<=Lee_test->num_elements_collapse; ibin++) {
	for(int jbin=1; jbin<=Lee_test->num_elements_collapse; jbin++) {
	
	  double cv_i = Lee_test->map_collapsed_prediction_wiLee[ibin];
	  double cv_j = Lee_test->map_collapsed_prediction_wiLee[jbin];

	  ////////// total
	  double cov_total     = Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin-1, jbin-1);
	  double total_sigma_i = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin-1, ibin-1) );
	  double total_sigma_j = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee(jbin-1, jbin-1) );
	  double total_correlation = cov_total/total_sigma_i/total_sigma_j;
	  double total_frac = cov_total/cv_i/cv_j;
	  matrix_total_correlation(ibin-1, jbin-1) = total_correlation;
	  matrix_total_frac(ibin-1, jbin-1) = total_frac;
	  if(ibin==jbin) {
	    h1_total_frac->SetBinContent(ibin, sqrt(total_frac));
	    h1_total_frac2->SetBinContent(ibin, total_frac/total_frac);
	  }
	
	  ////////// Xs
	  double cov_Xs     = Lee_test->matrix_collapse_absolute_covariance_wiLee_Xs(ibin-1, jbin-1);
	  double Xs_sigma_i = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_Xs(ibin-1, ibin-1) );
	  double Xs_sigma_j = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_Xs(jbin-1, jbin-1) );
	  double Xs_correlation = cov_Xs/Xs_sigma_i/Xs_sigma_j;
	  double Xs_frac = cov_Xs/cv_i/cv_j;
	  matrix_Xs_correlation(ibin-1, jbin-1) = Xs_correlation;
	  matrix_Xs_frac(ibin-1, jbin-1) = Xs_frac;
	  if(ibin==jbin) {
	    h1_Xs_frac->SetBinContent(ibin, sqrt(Xs_frac));
	    h1_Xs_frac2->SetBinContent(ibin, Xs_frac/total_frac);
	  }

	  ////////// flux
	  double cov_flux     = Lee_test->matrix_collapse_absolute_covariance_wiLee_flux(ibin-1, jbin-1);
	  double flux_sigma_i = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_flux(ibin-1, ibin-1) );
	  double flux_sigma_j = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_flux(jbin-1, jbin-1) );
	  double flux_correlation = cov_flux/flux_sigma_i/flux_sigma_j;
	  double flux_frac = cov_flux/cv_i/cv_j;
	  matrix_flux_correlation(ibin-1, jbin-1) = flux_correlation;
	  matrix_flux_frac(ibin-1, jbin-1) = flux_frac;
	  if(ibin==jbin) {
	    h1_flux_frac->SetBinContent(ibin, sqrt(flux_frac)); // -------------------------> setting
	    h1_flux_frac2->SetBinContent(ibin, flux_frac/total_frac);
	  }

	  ////////// detector
	  double cov_detector     = Lee_test->matrix_collapse_absolute_covariance_wiLee_detector(ibin-1, jbin-1);
	  double detector_sigma_i = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_detector(ibin-1, ibin-1) );
	  double detector_sigma_j = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_detector(jbin-1, jbin-1) );
	  double detector_correlation = cov_detector/detector_sigma_i/detector_sigma_j;
	  double detector_frac = cov_detector/cv_i/cv_j;
	  matrix_detector_correlation(ibin-1, jbin-1) = detector_correlation;
	  matrix_detector_frac(ibin-1, jbin-1) = detector_frac;
	  if(ibin==jbin) {
	    h1_detector_frac->SetBinContent(ibin, sqrt(detector_frac));
	    h1_detector_frac2->SetBinContent(ibin, detector_frac/total_frac);
	  }

	  ////////// additional
	  double cov_additional     = Lee_test->matrix_collapse_absolute_covariance_wiLee_additional(ibin-1, jbin-1);
	  double additional_sigma_i = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_additional(ibin-1, ibin-1) );
	  double additional_sigma_j = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_additional(jbin-1, jbin-1) );
	  double additional_correlation = cov_additional/additional_sigma_i/additional_sigma_j;
	  double additional_frac = cov_additional/cv_i/cv_j;
	  matrix_additional_correlation(ibin-1, jbin-1) = additional_correlation;
	  matrix_additional_frac(ibin-1, jbin-1) = additional_frac;
	  if(ibin==jbin) {
	    h1_additional_frac->SetBinContent(ibin, sqrt(additional_frac));
	    h1_additional_frac2->SetBinContent(ibin, additional_frac/total_frac);
	  }

	  ////////// mc_stat
	  double cov_mc_stat     = Lee_test->matrix_collapse_absolute_covariance_wiLee_mc_stat(ibin-1, jbin-1);
	  double mc_stat_sigma_i = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_mc_stat(ibin-1, ibin-1) );
	  double mc_stat_sigma_j = sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee_mc_stat(jbin-1, jbin-1) );
	  double mc_stat_correlation = cov_mc_stat/mc_stat_sigma_i/mc_stat_sigma_j;
	  double mc_stat_frac = cov_mc_stat/cv_i/cv_j;
	  matrix_mc_stat_correlation(ibin-1, jbin-1) = mc_stat_correlation;
	  matrix_mc_stat_frac(ibin-1, jbin-1) = mc_stat_frac;
	  if(ibin==jbin) {
	    h1_mc_stat_frac->SetBinContent(ibin, sqrt(mc_stat_frac));
	    h1_mc_stat_frac2->SetBinContent(ibin, mc_stat_frac/total_frac);
	  }

	  ////////// CV
	  if(ibin==jbin) {
	    double val_CV = Lee_test->map_collapsed_prediction_wiLee[ibin];
	    h1_CV->SetBinContent( ibin, val_CV );
	  
	    h1_CV->SetBinError( ibin, sqrt( cov_total ) ); // -------------------------> setting
	  }
	
	
	}// jbin
      }// ibin


      ////////////////////////////////////////


    

      int color_flux = kRed;
      int color_Xs = kBlue;
      int color_detector = kMagenta;
      int color_mc_stat = kGreen+1;
      int color_additional = kCyan;
    
    
    
      TCanvas *canv_h1_CV = new TCanvas("canv_h1_CV", "canv_h1_CV", 1200, 650);
      canv_h1_CV->SetLogy();
      func_canv_margin(canv_h1_CV, 0.15, 0.1, 0.05, 0.15);
    
      TH1D *h1_CV_clone = (TH1D*)h1_CV->Clone("h1_CV_clone");
        
      h1_CV->Draw("e2");
      h1_CV->SetMinimum(1e-4);
      h1_CV->SetMaximum(5e3);
      h1_CV->SetMarkerStyle(1);
      h1_CV->SetFillColor(kRed);
      func_title_size(h1_CV, 0.05, 0.05, 0.05, 0.05);
      func_xy_title(h1_CV, "Bin index", "Entries");
      h1_CV->GetYaxis()->CenterTitle();
      h1_CV->GetXaxis()->CenterTitle();

      h1_CV_clone->SetLineColor(kBlack);
      h1_CV_clone->Draw("same hist");

      canv_h1_CV->Update();
      double ymax = gPad->GetUymax();
      double ymin = gPad->GetUymin();
    
      double xx_array[6] = {26, 26*2, 26*3, 26*4, 26*4+11, 26*4+11*2};
      TLine *line_array[6];
      for(int idx=1; idx<=6; idx++) {
	double xx_val = xx_array[idx-1];
	line_array[idx-1] = new TLine( xx_val, 1e-4, xx_val, 5e3 );
	line_array[idx-1]->Draw("same");
	line_array[idx-1]->SetLineColor(kBlack);
	line_array[idx-1]->SetLineStyle(7);
      }

    
      TH1D *h1_data = new TH1D("h1_data", "", rows, 0, rows);
      for(int ibin=1; ibin<=rows; ibin++) {
	double value = Lee_test->map_data_spectrum_global_bin[ibin];
	h1_data->SetBinContent(ibin, value);
	h1_data->SetBinError(ibin, sqrt(value) );
      }

      h1_data->Draw("same pe");
      h1_data->SetLineColor(kBlue);
      h1_data->SetMarkerColor(kBlue);
    
      h1_CV_clone->Draw("same hist");
      h1_CV_clone->SetLineColor(kCyan);
    
      canv_h1_CV->SaveAs("canv_h1_CV_channels.png");


      
	TCanvas *canv_h1_total_frac = new TCanvas("canv_h1_total_frac", "canv_h1_total_frac", 1200, 650);
	//canv_h1_total_frac->SetLogy();
	func_canv_margin(canv_h1_total_frac, 0.15, 0.1, 0.05, 0.15);    
	h1_total_frac->Draw("hist");
	h1_total_frac->SetLineColor(kBlack);    
	func_title_size(h1_total_frac, 0.05, 0.05, 0.05, 0.05);
	func_xy_title(h1_total_frac, "Bin index", "Rel. Error");
	h1_total_frac->GetYaxis()->CenterTitle();
	h1_total_frac->GetXaxis()->CenterTitle();
	h1_total_frac->SetLineWidth(4);
	h1_total_frac->SetMaximum(2.5);
    	h1_total_frac->SetMinimum(0);
    
	h1_additional_frac->Draw("same");
	h1_additional_frac->SetLineColor(color_additional);
      			
	h1_mc_stat_frac->Draw("same");
	h1_mc_stat_frac->SetLineColor(color_mc_stat);
      		
	h1_detector_frac->Draw("same");
	h1_detector_frac->SetLineColor(color_detector);
      		
	h1_flux_frac->Draw("same");
	h1_flux_frac->SetLineColor(color_flux);
      			
	h1_Xs_frac->Draw("same");
	h1_Xs_frac->SetLineColor(color_Xs);
    
	h1_Xs_frac->Draw("same axis");
    
	canv_h1_total_frac->Update();
	ymax = gPad->GetUymax();
	ymin = gPad->GetUymin();
    
	TLine *line_array_frac[6];
	for(int idx=1; idx<=6; idx++) {
	double xx_val = xx_array[idx-1];
	line_array_frac[idx-1] = new TLine( xx_val, ymin, xx_val, ymax );
	line_array_frac[idx-1]->Draw("same");
	line_array_frac[idx-1]->SetLineColor(kBlack);
	line_array_frac[idx-1]->SetLineStyle(7);
	}
	canv_h1_total_frac->SaveAs("canv_h1_total_frac.png");

       
	///////////////////////////////////////////////////////////

	TCanvas *canv_h1_total2 = new TCanvas("canv_h1_total2", "canv_h1_total2", 1200, 650);
	//canv_h1_total2->SetLogy();
	func_canv_margin(canv_h1_total2, 0.15, 0.1, 0.05, 0.15);

	TH1D *h1_total2 = new TH1D("h1_total2", "", rows, 0, rows);
	h1_total2->Draw("hist");
	h1_total2->SetLineColor(kBlack);    
	func_title_size(h1_total2, 0.05, 0.05, 0.05, 0.05);
	func_xy_title(h1_total2, "Bin index", "Percentage (%)");
	h1_total2->GetYaxis()->CenterTitle();
	h1_total2->GetXaxis()->CenterTitle();
	h1_total2->SetLineWidth(4);
	h1_total2->SetMaximum(120);
	h1_total2->SetMinimum(0);

	// auto hs = new THStack("hs","Test of palette colored lego stack");

	h1_flux_frac2->Scale(100);
	h1_Xs_frac2->Scale(100);
	h1_detector_frac2->Scale(100);
	h1_mc_stat_frac2->Scale(100);
	h1_additional_frac2->Scale(100);
    
	// hs->Add( h1_flux_frac2 );
	// hs->Add( h1_Xs_frac2 );
	// hs->Add( h1_detector_frac2 );
	// hs->Add( h1_mc_stat_frac2 );
	// hs->Add( h1_additional_frac2 );

	h1_flux_frac2->SetFillColor(color_flux);
	h1_Xs_frac2->SetFillColor(color_Xs);
	h1_detector_frac2->SetFillColor(color_detector);
	h1_mc_stat_frac2->SetFillColor(color_mc_stat);
	h1_additional_frac2->SetFillColor(color_additional);

	h1_flux_frac2->SetLineColor(color_flux);
	h1_Xs_frac2->SetLineColor(color_Xs);
	h1_detector_frac2->SetLineColor(color_detector);
	h1_mc_stat_frac2->SetLineColor(color_mc_stat);
	h1_additional_frac2->SetLineColor(color_additional);

	// hs->Draw("same pfc nostack");

	h1_additional_frac2->Add( h1_mc_stat_frac2 );
	h1_additional_frac2->Add( h1_flux_frac2 );
	h1_additional_frac2->Add( h1_Xs_frac2 );
	h1_additional_frac2->Add( h1_detector_frac2 );
	h1_additional_frac2->Draw("same f hist");
    
	h1_mc_stat_frac2->Add( h1_flux_frac2 );
	h1_mc_stat_frac2->Add( h1_Xs_frac2 );
	h1_mc_stat_frac2->Add( h1_detector_frac2 );
	h1_mc_stat_frac2->Draw("same f hist");
    
	h1_detector_frac2->Add( h1_flux_frac2 );
	h1_detector_frac2->Add( h1_Xs_frac2 );
	h1_detector_frac2->Draw("same f hist");
    
	h1_Xs_frac2->Add( h1_flux_frac2 );
	h1_Xs_frac2->Draw("same f hist");

	h1_flux_frac2->Draw("same f hist");
    
	h1_flux_frac2->SetLineColor(kBlack);
	h1_Xs_frac2->SetLineColor(kBlack);
	h1_detector_frac2->SetLineColor(kBlack);
	h1_mc_stat_frac2->SetLineColor(kBlack);
	h1_additional_frac2->SetLineColor(kBlack);
    
	h1_flux_frac2->SetLineWidth(1);
	h1_Xs_frac2->SetLineWidth(1);
	h1_detector_frac2->SetLineWidth(1);
	h1_mc_stat_frac2->SetLineWidth(1);
	h1_additional_frac2->SetLineWidth(1);
    
	canv_h1_total2->Update();
	ymax = gPad->GetUymax();
	ymin = gPad->GetUymin();
    
	//TLine *line_array_frac[6];
	for(int idx=1; idx<=6; idx++) {
	double xx_val = xx_array[idx-1];
	line_array_frac[idx-1] = new TLine( xx_val, ymin, xx_val, ymax );
	line_array_frac[idx-1]->Draw("same");
	line_array_frac[idx-1]->SetLineColor(kBlack);
	line_array_frac[idx-1]->SetLineStyle(7);
	}


    
	h1_total2->Draw("same axis");

    
    
	canv_h1_total2->SaveAs("canv_h1_total2_tack.png");
      
      ///////////////////////////
      
	TCanvas *canv_total_correlation = new TCanvas("canv_total_correlation", "canv_total_correlation", 900, 850);
	func_canv_margin(canv_total_correlation, 0.15, 0.15, 0.15, 0.15);

	TH2D *h2_total_correlation = new TH2D("h2_total_correlation", "", rows, 0, rows, rows, 0, rows);
	h2_total_correlation->Draw("colz");
	h2_total_correlation->GetZaxis()->SetRangeUser(-1,1);
	h2_total_correlation->GetZaxis()->SetLabelSize(0.05);
    
	matrix_total_correlation.Draw("same colz");
    
	func_title_size(h2_total_correlation, 0.05, 0.05, 0.05, 0.05);
	func_xy_title(h2_total_correlation, "Bin index", "Bin index");
	h2_total_correlation->GetYaxis()->CenterTitle();
	h2_total_correlation->GetXaxis()->CenterTitle();

	TLine *line_array_row[6];
	for(int idx=1; idx<=6; idx++) {
	double xx = xx_array[idx-1];
	line_array_row[idx-1] = new TLine( xx, 0, xx, rows );
	line_array_row[idx-1]->Draw("same");
	line_array_row[idx-1]->SetLineColor(kRed);
	}

	TLine *line_array_col[6];
	for(int idx=1; idx<=6; idx++) {
	double xx = xx_array[idx-1];
	line_array_col[idx-1] = new TLine( 0, xx, rows, xx );
	line_array_col[idx-1]->Draw("same");
	line_array_col[idx-1]->SetLineColor(kRed);
	}
    
	canv_total_correlation->SaveAs("canv_total_correlation.png");
      

    }


  

  //////////////////////
  //////////////////////
  
  bool flag_numuCC_FC = 0;// 1
  bool flag_numuCC_PC = 0;// 2
  bool flag_CCpi0_FC  = 0;// 3
  bool flag_CCpi0_PC  = 0;// 4
  bool flag_NCpi0     = 0;// 5
  bool flag_nueCC_PC  = 0;// 6
  bool flag_nueCC_HghE_FC = 0;// 7
  bool flag_nueCC_LowE_FC = 0;// 8
  bool flag_nueCC_LowE_FC_wopi0    = 0;// 9
  bool flag_nueCC_LowE_FC_combined = 0;// 10
  bool flag_both_numuCC            = 0;// 11

  bool flag_nueCC_LowE_FC_by_numuCC     = 0;// 16
  bool flag_nueCC_LowE_FC_by_numuCC_pi0 = 0;// 17

  

  if( flag_nueCC_LowE_FC_by_numuCC_pi0 ) {

    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;
           
    //////////
    TMatrixD matrix_trans(Lee_test->num_elements_collapse, 8 + 26*2+33);// old_world, new_world
    for(int ibin=1; ibin<=8; ibin++) {
      matrix_trans(ibin-1, ibin-1) = 1;
    }

    for(int ibin=1; ibin<=26*2+33; ibin++) {
      matrix_trans(26*2+ibin-1, 8+ibin-1) = 1;
    }

    TMatrixD matrix_trans_T( 8 + 26*2+33, Lee_test->num_elements_collapse);
    matrix_trans_T.Transpose( matrix_trans );

    ///
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;    

    ///
    TMatrixD matrix_pred_old_world(Lee_test->num_elements_collapse, 1);
    for(int ibin=1; ibin<=Lee_test->num_elements_collapse; ibin++) {
      matrix_pred_old_world(ibin-1, 0) = Lee_test->map_collapsed_prediction_wiLee[ibin];
    }    
    TMatrixD matrix_pred_new_world = matrix_trans_T * matrix_pred_old_world;
    for(int ibin=1; ibin<=8 + 26*2+33; ibin++) {
      map_h1_pred[ibin] = matrix_pred_new_world(ibin-1, 0);
    }

    ///
    TMatrixD matrix_data_old_world(Lee_test->num_elements_collapse, 1);
    for(int ibin=1; ibin<=Lee_test->num_elements_collapse; ibin++) {
      matrix_data_old_world(ibin-1, 0) = Lee_test->map_data_spectrum_global_bin_scaled[ibin];
    }    
    TMatrixD matrix_data_new_world = matrix_trans_T * matrix_data_old_world;
    for(int ibin=1; ibin<=8 + 26*2+33; ibin++) {
      map_h1_data[ibin] = matrix_data_new_world(ibin-1, 0);
    }
    
    Lee_test->Exe_Constraint_GOF(8, 26*2+33, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 17);
    
  }

  
  if( flag_nueCC_LowE_FC_by_numuCC ) {

    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;
           
    //////////
    TMatrixD matrix_trans(Lee_test->num_elements_collapse, 8+18 + 26*2);// old_world, new_world
    for(int ibin=1; ibin<=8+18; ibin++) {
      matrix_trans(ibin-1, ibin-1) = 1;
    }

    for(int ibin=1; ibin<=26; ibin++) {
      matrix_trans(26*2+ibin-1, 8+18+ibin-1) = 1;
    }

    for(int ibin=1; ibin<=26; ibin++) {
      matrix_trans(26*3+ibin-1, 8+18+26+ibin-1) = 1;
    }

    TMatrixD matrix_trans_T( 8+18 + 26*2, Lee_test->num_elements_collapse);
    matrix_trans_T.Transpose( matrix_trans );

    ///
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;    

    ///
    TMatrixD matrix_pred_old_world(Lee_test->num_elements_collapse, 1);
    for(int ibin=1; ibin<=Lee_test->num_elements_collapse; ibin++) {
      matrix_pred_old_world(ibin-1, 0) = Lee_test->map_collapsed_prediction_wiLee[ibin];
    }    
    TMatrixD matrix_pred_new_world = matrix_trans_T * matrix_pred_old_world;
    for(int ibin=1; ibin<=8+18 + 26*2; ibin++) {
      map_h1_pred[ibin] = matrix_pred_new_world(ibin-1, 0);
    }

    ///
    TMatrixD matrix_data_old_world(Lee_test->num_elements_collapse, 1);
    for(int ibin=1; ibin<=Lee_test->num_elements_collapse; ibin++) {
      matrix_data_old_world(ibin-1, 0) = Lee_test->map_data_spectrum_global_bin_scaled[ibin];
    }    
    TMatrixD matrix_data_new_world = matrix_trans_T * matrix_data_old_world;
    for(int ibin=1; ibin<=8+18 + 26*2; ibin++) {
      map_h1_data[ibin] = matrix_data_new_world(ibin-1, 0);
    }
    
    Lee_test->Exe_Constraint_GOF(8+18, 26*2, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 16);
    
  }

  
  if( flag_both_numuCC ) {

    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;

    int nn = 26*2;
    
    for(int ibin=1; ibin<=nn; ibin++) {
      map_h1_pred[ibin] = Lee_test->map_collapsed_prediction_wiLee[ ibin + nn ];
      map_h1_data[ibin] = Lee_test->map_data_spectrum_global_bin_scaled[ ibin + nn ];
    }
    
    TMatrixD matrix_collapse_absolute_cov(nn,nn);

    for(int ibin=1; ibin<=nn; ibin++) {
      for(int jbin=1; jbin<=nn; jbin++) {
	matrix_collapse_absolute_cov(ibin-1, jbin-1) = Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin+nn-1, jbin+nn-1);
      }
    }
     
    Lee_test->Exe_Goodness_of_Fit(map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 11);
  }
  
  if( flag_numuCC_FC ) {    
    ///////// set collapsed local2global index
    int num_elements_collapse = 0;  
    map<int, int>map_collapse_channel_bins;// new-world channels
    int line_map = 0;

    /// --------------------------> 
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[3];

    ///
    map<int, map<int, int> >map_collapse_channel_local2global;  
    int size_map = map_collapse_channel_bins.size();
    for(int idx=1; idx<=size_map; idx++) {
      for(int ibin=1; ibin<=map_collapse_channel_bins[idx]; ibin++) {
	num_elements_collapse++;
	map_collapse_channel_local2global[idx][ibin] = num_elements_collapse;
      }    
    }
    cout<<endl<<" ---> num_elements_collapse "<<num_elements_collapse<<endl<<endl;
        
    TMatrixD matrix_trans( Lee_test->num_elements_collapse, num_elements_collapse );// old_world, new_world
    TMatrixD matrix_trans_T( num_elements_collapse, Lee_test->num_elements_collapse );
    
    map<int, double>map_collapsed_pred;
    map_collapsed_pred.clear();

    ///////// set transform matrix --------------------------> 
    map<int, int>map_oldworld2neworld;
    map_oldworld2neworld[3] = 1;

    /////////
    for(auto it=map_oldworld2neworld.begin(); it!=map_oldworld2neworld.end(); it++) {
      int old_world_channel = it->first;
      int new_world_channel = it->second;

      int bins = Lee_test->map_collapse_channel_local2global[old_world_channel].size();
      for(int ibin=1; ibin<=bins; ibin++ ) {
	int old_world_global_bin = Lee_test->map_collapse_channel_local2global[old_world_channel][ibin];
	int new_world_global_bin = map_collapse_channel_local2global[new_world_channel][ibin];      
 	
	matrix_trans(old_world_global_bin-1, new_world_global_bin-1) = 1;
	map_collapsed_pred[ new_world_global_bin ] += Lee_test->map_collapsed_prediction_wiLee[old_world_global_bin];      
      }// ibin
    }// it

    matrix_trans_T.Transpose( matrix_trans );
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;
    
    ////////// set data collapsed --------------------------> 
    map<int, double>map_collapsed_data;
    map_collapsed_data.clear();
    int old_ch = 0;
    int new_ch = 0;
    int bins_ch = 0;
    
    old_ch = 3; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = map_collapsed_pred[ibin];
      map_h1_data[ibin] = map_collapsed_data[ibin];
    }

    Lee_test->Exe_Goodness_of_Fit(map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 1);  
  }

  //////////////////////
  //////////////////////
  
  if( flag_numuCC_PC ) {    
    ///////// set collapsed local2global index
    int num_elements_collapse = 0;  
    map<int, int>map_collapse_channel_bins;// new-world channels
    int line_map = 0;

    /// --------------------------> 
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[4];

    ///
    map<int, map<int, int> >map_collapse_channel_local2global;  
    int size_map = map_collapse_channel_bins.size();
    for(int idx=1; idx<=size_map; idx++) {
      for(int ibin=1; ibin<=map_collapse_channel_bins[idx]; ibin++) {
	num_elements_collapse++;
	map_collapse_channel_local2global[idx][ibin] = num_elements_collapse;
      }    
    }
    cout<<endl<<" ---> num_elements_collapse "<<num_elements_collapse<<endl<<endl;
        
    TMatrixD matrix_trans( Lee_test->num_elements_collapse, num_elements_collapse );// old_world, new_world
    TMatrixD matrix_trans_T( num_elements_collapse, Lee_test->num_elements_collapse );
    
    map<int, double>map_collapsed_pred;
    map_collapsed_pred.clear();

    ///////// set transform matrix --------------------------> 
    map<int, int>map_oldworld2neworld;
    map_oldworld2neworld[4] = 1;

    /////////
    for(auto it=map_oldworld2neworld.begin(); it!=map_oldworld2neworld.end(); it++) {
      int old_world_channel = it->first;
      int new_world_channel = it->second;

      int bins = Lee_test->map_collapse_channel_local2global[old_world_channel].size();
      for(int ibin=1; ibin<=bins; ibin++ ) {
	int old_world_global_bin = Lee_test->map_collapse_channel_local2global[old_world_channel][ibin];
	int new_world_global_bin = map_collapse_channel_local2global[new_world_channel][ibin];      
 	
	matrix_trans(old_world_global_bin-1, new_world_global_bin-1) = 1;
	map_collapsed_pred[ new_world_global_bin ] += Lee_test->map_collapsed_prediction_wiLee[old_world_global_bin];      
      }// ibin
    }// it

    matrix_trans_T.Transpose( matrix_trans );
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;
    
    ////////// set data collapsed --------------------------> 
    map<int, double>map_collapsed_data;
    map_collapsed_data.clear();
    int old_ch = 0;
    int new_ch = 0;
    int bins_ch = 0;
    
    old_ch = 4; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = map_collapsed_pred[ibin];
      map_h1_data[ibin] = map_collapsed_data[ibin];
    }

    Lee_test->Exe_Goodness_of_Fit(map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 2);  
  }

  ////////////////////// 
  //////////////////////
  
  if( flag_CCpi0_FC ) {    
    ///////// set collapsed local2global index
    int num_elements_collapse = 0;  
    map<int, int>map_collapse_channel_bins;// new-world channels
    int line_map = 0;

    /// --------------------------> 
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[5];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[3];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[4];

    ///
    map<int, map<int, int> >map_collapse_channel_local2global;  
    int size_map = map_collapse_channel_bins.size();
    for(int idx=1; idx<=size_map; idx++) {
      for(int ibin=1; ibin<=map_collapse_channel_bins[idx]; ibin++) {
	num_elements_collapse++;
	map_collapse_channel_local2global[idx][ibin] = num_elements_collapse;
      }    
    }
    cout<<endl<<" ---> num_elements_collapse "<<num_elements_collapse<<endl<<endl;
        
    TMatrixD matrix_trans( Lee_test->num_elements_collapse, num_elements_collapse );// old_world, new_world
    TMatrixD matrix_trans_T( num_elements_collapse, Lee_test->num_elements_collapse );
    
    map<int, double>map_collapsed_pred;
    map_collapsed_pred.clear();

    ///////// set transform matrix --------------------------> 
    map<int, int>map_oldworld2neworld;
    map_oldworld2neworld[5] = 1;
    map_oldworld2neworld[3] = 2;
    map_oldworld2neworld[4] = 3;

    /////////
    for(auto it=map_oldworld2neworld.begin(); it!=map_oldworld2neworld.end(); it++) {
      int old_world_channel = it->first;
      int new_world_channel = it->second;

      int bins = Lee_test->map_collapse_channel_local2global[old_world_channel].size();
      for(int ibin=1; ibin<=bins; ibin++ ) {
	int old_world_global_bin = Lee_test->map_collapse_channel_local2global[old_world_channel][ibin];
	int new_world_global_bin = map_collapse_channel_local2global[new_world_channel][ibin];      
 	
	matrix_trans(old_world_global_bin-1, new_world_global_bin-1) = 1;
	map_collapsed_pred[ new_world_global_bin ] += Lee_test->map_collapsed_prediction_wiLee[old_world_global_bin];      
      }// ibin
    }// it

    matrix_trans_T.Transpose( matrix_trans );
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;
    
    ////////// set data collapsed --------------------------> 
    map<int, double>map_collapsed_data;
    map_collapsed_data.clear();
    int old_ch = 0;
    int new_ch = 0;
    int bins_ch = 0;
    
    old_ch = 5; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 3; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 4; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = map_collapsed_pred[ibin];
      map_h1_data[ibin] = map_collapsed_data[ibin];
    }

    int num_Y = map_collapse_channel_bins[1];
    int num_X = 0;
    for(int idx=2; idx<=size_map; idx++) {
      num_X += map_collapse_channel_bins[idx];
    }
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 3);
  }

  ////////////////////// 
  //////////////////////
  
  if( flag_CCpi0_PC ) {    
    ///////// set collapsed local2global index
    int num_elements_collapse = 0;  
    map<int, int>map_collapse_channel_bins;// new-world channels
    int line_map = 0;

    /// --------------------------> 
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[6];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[3];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[4];

    ///
    map<int, map<int, int> >map_collapse_channel_local2global;  
    int size_map = map_collapse_channel_bins.size();
    for(int idx=1; idx<=size_map; idx++) {
      for(int ibin=1; ibin<=map_collapse_channel_bins[idx]; ibin++) {
	num_elements_collapse++;
	map_collapse_channel_local2global[idx][ibin] = num_elements_collapse;
      }    
    }
    cout<<endl<<" ---> num_elements_collapse "<<num_elements_collapse<<endl<<endl;
        
    TMatrixD matrix_trans( Lee_test->num_elements_collapse, num_elements_collapse );// old_world, new_world
    TMatrixD matrix_trans_T( num_elements_collapse, Lee_test->num_elements_collapse );
    
    map<int, double>map_collapsed_pred;
    map_collapsed_pred.clear();

    ///////// set transform matrix --------------------------> 
    map<int, int>map_oldworld2neworld;
    map_oldworld2neworld[6] = 1;
    map_oldworld2neworld[3] = 2;
    map_oldworld2neworld[4] = 3;

    /////////
    for(auto it=map_oldworld2neworld.begin(); it!=map_oldworld2neworld.end(); it++) {
      int old_world_channel = it->first;
      int new_world_channel = it->second;

      int bins = Lee_test->map_collapse_channel_local2global[old_world_channel].size();
      for(int ibin=1; ibin<=bins; ibin++ ) {
	int old_world_global_bin = Lee_test->map_collapse_channel_local2global[old_world_channel][ibin];
	int new_world_global_bin = map_collapse_channel_local2global[new_world_channel][ibin];      
 	
	matrix_trans(old_world_global_bin-1, new_world_global_bin-1) = 1;
	map_collapsed_pred[ new_world_global_bin ] += Lee_test->map_collapsed_prediction_wiLee[old_world_global_bin];      
      }// ibin
    }// it

    matrix_trans_T.Transpose( matrix_trans );
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;
    
    ////////// set data collapsed --------------------------> 
    map<int, double>map_collapsed_data;
    map_collapsed_data.clear();
    int old_ch = 0;
    int new_ch = 0;
    int bins_ch = 0;
    
    old_ch = 6; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 3; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 4; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = map_collapsed_pred[ibin];
      map_h1_data[ibin] = map_collapsed_data[ibin];
    }

    int num_Y = map_collapse_channel_bins[1];
    int num_X = 0;
    for(int idx=2; idx<=size_map; idx++) {
      num_X += map_collapse_channel_bins[idx];
    }
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 4);
  }

  ////////////////////// 
  //////////////////////
  
  if( flag_NCpi0 ) {    
    ///////// set collapsed local2global index
    int num_elements_collapse = 0;  
    map<int, int>map_collapse_channel_bins;// new-world channels
    int line_map = 0;

    /// --------------------------> 
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[7];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[3];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[4];

    ///
    map<int, map<int, int> >map_collapse_channel_local2global;  
    int size_map = map_collapse_channel_bins.size();
    for(int idx=1; idx<=size_map; idx++) {
      for(int ibin=1; ibin<=map_collapse_channel_bins[idx]; ibin++) {
	num_elements_collapse++;
	map_collapse_channel_local2global[idx][ibin] = num_elements_collapse;
      }    
    }
    cout<<endl<<" ---> num_elements_collapse "<<num_elements_collapse<<endl<<endl;
        
    TMatrixD matrix_trans( Lee_test->num_elements_collapse, num_elements_collapse );// old_world, new_world
    TMatrixD matrix_trans_T( num_elements_collapse, Lee_test->num_elements_collapse );
    
    map<int, double>map_collapsed_pred;
    map_collapsed_pred.clear();

    ///////// set transform matrix --------------------------> 
    map<int, int>map_oldworld2neworld;
    map_oldworld2neworld[7] = 1;
    map_oldworld2neworld[3] = 2;
    map_oldworld2neworld[4] = 3;

    /////////
    for(auto it=map_oldworld2neworld.begin(); it!=map_oldworld2neworld.end(); it++) {
      int old_world_channel = it->first;
      int new_world_channel = it->second;

      int bins = Lee_test->map_collapse_channel_local2global[old_world_channel].size();
      for(int ibin=1; ibin<=bins; ibin++ ) {
	int old_world_global_bin = Lee_test->map_collapse_channel_local2global[old_world_channel][ibin];
	int new_world_global_bin = map_collapse_channel_local2global[new_world_channel][ibin];      
 	
	matrix_trans(old_world_global_bin-1, new_world_global_bin-1) = 1;
	map_collapsed_pred[ new_world_global_bin ] += Lee_test->map_collapsed_prediction_wiLee[old_world_global_bin];      
      }// ibin
    }// it

    matrix_trans_T.Transpose( matrix_trans );
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;
    
    ////////// set data collapsed --------------------------> 
    map<int, double>map_collapsed_data;
    map_collapsed_data.clear();
    int old_ch = 0;
    int new_ch = 0;
    int bins_ch = 0;
    
    old_ch = 7; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 3; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 4; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = map_collapsed_pred[ibin];
      map_h1_data[ibin] = map_collapsed_data[ibin];
    }

    int num_Y = map_collapse_channel_bins[1];
    int num_X = 0;
    for(int idx=2; idx<=size_map; idx++) {
      num_X += map_collapse_channel_bins[idx];
    }
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 5);
  }

  ////////////////////// 
  //////////////////////
  
  if( flag_nueCC_PC ) {    
    ///////// set collapsed local2global index
    int num_elements_collapse = 0;  
    map<int, int>map_collapse_channel_bins;// new-world channels
    int line_map = 0;

    /// --------------------------> 
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[2];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[3];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[4];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[5];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[6];
    line_map++; map_collapse_channel_bins[line_map] = Lee_test->map_input_spectrum_ch_binnum[7];

    ///
    map<int, map<int, int> >map_collapse_channel_local2global;  
    int size_map = map_collapse_channel_bins.size();
    for(int idx=1; idx<=size_map; idx++) {
      for(int ibin=1; ibin<=map_collapse_channel_bins[idx]; ibin++) {
	num_elements_collapse++;
	map_collapse_channel_local2global[idx][ibin] = num_elements_collapse;
      }    
    }
    cout<<endl<<" ---> num_elements_collapse "<<num_elements_collapse<<endl<<endl;
        
    TMatrixD matrix_trans( Lee_test->num_elements_collapse, num_elements_collapse );// old_world, new_world
    TMatrixD matrix_trans_T( num_elements_collapse, Lee_test->num_elements_collapse );
    
    map<int, double>map_collapsed_pred;
    map_collapsed_pred.clear();

    ///////// set transform matrix --------------------------> 
    map<int, int>map_oldworld2neworld;
    map_oldworld2neworld[2] = 1;
    map_oldworld2neworld[3] = 2;
    map_oldworld2neworld[4] = 3;
    map_oldworld2neworld[5] = 4;
    map_oldworld2neworld[6] = 5;
    map_oldworld2neworld[7] = 6;


    // for(int ibin=1; ibin<=Lee_test->num_elements_collapse; ibin++) cout<<Form(" ---> check %3d, %12.4f, %12.4f",
    // 									      ibin,
    // 									      Lee_test->map_collapsed_prediction_wiLee[ibin],
    // 									      sqrt( Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin-1, ibin-1) )
    // 									      )<<endl;
    
    /////////
    for(auto it=map_oldworld2neworld.begin(); it!=map_oldworld2neworld.end(); it++) {
      int old_world_channel = it->first;
      int new_world_channel = it->second;

      int bins = Lee_test->map_collapse_channel_local2global[old_world_channel].size();
      for(int ibin=1; ibin<=bins; ibin++ ) {
	int old_world_global_bin = Lee_test->map_collapse_channel_local2global[old_world_channel][ibin];
	int new_world_global_bin = map_collapse_channel_local2global[new_world_channel][ibin];      
 	
	matrix_trans(old_world_global_bin-1, new_world_global_bin-1) = 1;
	map_collapsed_pred[ new_world_global_bin ] += Lee_test->map_collapsed_prediction_wiLee[old_world_global_bin];      
      }// ibin
    }// it

    matrix_trans_T.Transpose( matrix_trans );
    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;
    
    ////////// set data collapsed --------------------------> 
    map<int, double>map_collapsed_data;
    map_collapsed_data.clear();
    int old_ch = 0;
    int new_ch = 0;
    int bins_ch = 0;
    
    old_ch = 2; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 3; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 4; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 5; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 6; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    old_ch = 7; new_ch = map_oldworld2neworld[old_ch]; bins_ch = Lee_test->map_collapse_channel_local2global[old_ch].size();    
    for(int ibin=1; ibin<=bins_ch; ibin++) {
      int extra_new = 0; for(int idx=1; idx<new_ch; idx++) extra_new += map_collapse_channel_bins[idx];
      map_collapsed_data[extra_new+ibin] += Lee_test->map_data_spectrum_ch_bin_scaled[old_ch][ibin];
    }
     
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = map_collapsed_pred[ibin];
      map_h1_data[ibin] = map_collapsed_data[ibin];
    }

    int num_Y = map_collapse_channel_bins[1];
    int num_X = 0;
    for(int idx=2; idx<=size_map; idx++) {
      num_X += map_collapse_channel_bins[idx];
    }
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 6);
  }

  ////////////////////// 
  //////////////////////
  
  if( flag_nueCC_HghE_FC ) {
       
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;

    int bins_LowE = 8;    
    int num_elements_collapse = Lee_test->matrix_collapse_absolute_covariance_wiLee.GetNrows() - bins_LowE;
    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = Lee_test->map_collapsed_prediction_wiLee[ibin + bins_LowE];
      map_h1_data[ibin] = Lee_test->map_data_spectrum_global_bin_scaled[ibin + bins_LowE];
    }
    
    int num_Y = Lee_test->map_input_spectrum_ch_binnum[1] - bins_LowE;
    int num_X = num_elements_collapse - num_Y;
    
    TMatrixD matrix_collapse_absolute_cov(num_elements_collapse, num_elements_collapse);
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      for(int jbin=1; jbin<=num_elements_collapse; jbin++) {
	matrix_collapse_absolute_cov(ibin-1, jbin-1) = Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin+bins_LowE-1, jbin+bins_LowE-1);
      }
    }
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 7);
  }

  ////////////////////// 
  //////////////////////
  
  if( flag_nueCC_LowE_FC ) {
    
    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;

    int bins_LowE = 8;    
    int num_elements_collapse = Lee_test->matrix_collapse_absolute_covariance_wiLee.GetNrows();
    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = Lee_test->map_collapsed_prediction_wiLee[ibin];
      map_h1_data[ibin] = Lee_test->map_data_spectrum_global_bin_scaled[ibin];
    }
    
    int num_Y = bins_LowE;
    int num_X = num_elements_collapse - num_Y;
    
    TMatrixD matrix_collapse_absolute_cov = Lee_test->matrix_collapse_absolute_covariance_wiLee;
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 8);  
  }

  ////////////////////// 
  //////////////////////
  
  if( flag_nueCC_LowE_FC_wopi0 ) {

    ////////// 
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;

    int bins_LowE = 8;    
    int num_elements_collapse = Lee_test->matrix_collapse_absolute_covariance_wiLee.GetNrows() - 11*3;
    
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      map_h1_pred[ibin] = Lee_test->map_collapsed_prediction_wiLee[ibin];
      map_h1_data[ibin] = Lee_test->map_data_spectrum_global_bin_scaled[ibin];
    }
    
    int num_Y = bins_LowE;
    int num_X = num_elements_collapse - num_Y;
    
    TMatrixD matrix_collapse_absolute_cov(num_elements_collapse, num_elements_collapse);
    for(int ibin=1; ibin<=num_elements_collapse; ibin++) {
      for(int jbin=1; jbin<=num_elements_collapse; jbin++) {
	matrix_collapse_absolute_cov(ibin-1, jbin-1) = Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin-1, jbin-1);
      }
    }
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 9);  
  }

  ////////////////////// 
  //////////////////////
    
  if( flag_nueCC_LowE_FC_combined ) {

    int bins_LowE = 8;
    
    int rows_old_world = Lee_test->matrix_collapse_absolute_covariance_wiLee.GetNrows();
    int rows_new_world = rows_old_world - (bins_LowE-1);// combined first 8 bins
    TMatrixD matrix_trans(rows_old_world, rows_new_world);// old_world, new_world
    TMatrixD matrix_trans_T(rows_new_world, rows_old_world);

    for(int iold=1; iold<=rows_old_world; iold++) {
      int inew = 0;      
      if(iold<=bins_LowE) {
	inew = 0;
	matrix_trans(iold-1, inew) = 1;	
      }
      else {
	inew = iold-bins_LowE;
	matrix_trans(iold-1, inew) = 1;
      }
    }

    matrix_trans_T.Transpose( matrix_trans );

    TMatrixD matrix_collapse_absolute_cov = matrix_trans_T * (Lee_test->matrix_collapse_absolute_covariance_wiLee) * matrix_trans;
    
    TMatrixD matrix_pred_old(rows_old_world, 1);
    TMatrixD matrix_data_old(rows_old_world, 1);
    for(int ibin=1; ibin<=rows_old_world; ibin++) {
      matrix_pred_old(ibin-1, 0) = Lee_test->map_collapsed_prediction_wiLee[ibin];
      matrix_data_old(ibin-1, 0) = Lee_test->map_data_spectrum_global_bin_scaled[ibin];
    }
    TMatrixD matrix_pred_new = matrix_trans_T * matrix_pred_old;
    TMatrixD matrix_data_new = matrix_trans_T * matrix_data_old;
    map<int, double>map_h1_pred;
    map<int, double>map_h1_data;
    for(int ibin=1; ibin<=rows_new_world; ibin++) {
      map_h1_pred[ibin] = matrix_pred_new(ibin-1, 0);
      map_h1_data[ibin] = matrix_data_new(ibin-1, 0);
    }

    int num_Y = 1;
    int num_X = rows_new_world - num_Y;
    
    Lee_test->Exe_Constraint_GOF(num_Y, num_X, map_h1_pred, map_h1_data, matrix_collapse_absolute_cov, 10);  
  }
  
  ////////////////////////////////////////////////////////////////////////////////// Lee strength fitting, lf

  
  if( 0 ) {
    //Lee_test->scaleF_POT = 69.5/5.0549670;
    Lee_test->scaleF_Lee = 1; 
    Lee_test->Set_Collapse();

    //Lee_test->Set_fake_data_Asimov();
    Lee_test->Set_measured_data();

    Lee_test->Minimization_Lee_strength_FullCov(1, 0);// scaleF_Lee, fixed or not

    double gmin = Lee_test->minimization_chi2;
  
    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
				Lee_test->minimization_chi2,
				Lee_test->minimization_Lee_strength_val,
				Lee_test->minimization_Lee_strength_err
				)<<endl<<endl;

    TGraph *gh_scan = new TGraph();
    double slow = 0;
    double shgh = 3;
    int nscan = 50;
    double step = (shgh-slow)/nscan;
    for(int idx=1; idx<=nscan; idx++) {
      cout<<" ---> scan "<<idx<<endl;
      double val_s = slow + (idx-1)*step;
      Lee_test->Minimization_Lee_strength_FullCov(val_s, 1);
      double val_chi2 = Lee_test->minimization_chi2;
      gh_scan->SetPoint( gh_scan->GetN(), val_s, val_chi2 - gmin);
    }

    cout<<" ---> check "<<gh_scan->Eval(1.)<<endl;
    
    gh_scan->Draw("al");
    gh_scan->GetXaxis()->SetTitle("LEE strength");
    gh_scan->GetYaxis()->SetTitle("#Delta#chi^{2}");
  }
  
  ////////////////////////////////////////////////////////////////////////////////// variation

  if( 1 ){
    
    Lee_test->scaleF_Lee = 1;//////////////////// true8Lee
    Lee_test->Set_Collapse();
    
    Lee_test->Set_fake_data_Asimov();
    // for(int ibin=0; ibin<Lee_test->matrix_collapse_absolute_covariance_wiLee.GetNrows(); ibin++)
    //   cout<<Form(" ---> %3d %12.4f %12.4f",
    // 		 ibin+1,
    // 		 Lee_test->map_fake_data[ibin],
    // 		 sqrt(Lee_test->matrix_collapse_absolute_covariance_wiLee(ibin, ibin))
    // 		 )<<endl;

    
    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// scaleF_Lee, fixed or not
    
    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
  				Lee_test->minimization_chi2,
  				Lee_test->minimization_Lee_strength_val,
  				Lee_test->minimization_Lee_strength_err
  				)<<endl<<endl;
  }



  
  if(0) {
    
    //Lee_test->scaleF_POT = 1;

    //Lee_test->scaleF_POT = 5.327231;
    //Lee_test->scaleF_POT = 69./5.327231;
    //Lee_test->scaleF_POT = 121./5.327231;
    
    ///////

    Lee_test->scaleF_Lee = 1;
    Lee_test->Set_Collapse();
    Lee_test->Set_fake_data_Asimov(); 
    Lee_test->Minimization_Lee_strength_FullCov(0, 1);
    double chi2_null8sm = Lee_test->minimization_chi2;
    
    Lee_test->scaleF_Lee = 0;
    Lee_test->Set_Collapse();
    Lee_test->Set_fake_data_Asimov(); 
    Lee_test->Minimization_Lee_strength_FullCov(1, 1);
    double chi2_null8Lee = Lee_test->minimization_chi2;

    cout<<endl<<" ---> check null@sm / Lee "<<sqrt(chi2_null8sm)<<"\t"<<sqrt(chi2_null8Lee)<<endl<<endl;
    
  }

  
  if(0) {
    Lee_test->scaleF_Lee = 1;    
    Lee_test->Set_Collapse();
    Lee_test->Set_Variations( 2 );
    Lee_test->Set_fake_data_Variation( 1 );
    Lee_test->Minimization_Lee_strength_FullCov(2, 0);// scaleF_Lee, fixed or not

    cout<<endl<<TString::Format(" ---> Best fit of Lee strength: chi2 %6.2f, %5.2f +/- %5.2f",
  				Lee_test->minimization_chi2,
  				Lee_test->minimization_Lee_strength_val,
  				Lee_test->minimization_Lee_strength_err
  				)<<endl<<endl;

  }

  

  
  if( 0 ) {
    // Lee_test->scaleF_Lee = 1; 
    // Lee_test->Set_Collapse();
    // Lee_test->Set_Variations( 2 );
    // Lee_test->Set_fake_data_Variation( 2 );
    // Lee_test->Minimization_Lee_strength_FullCov(2, 0);// scaleF_Lee, fixed or not

 
    double chi2_null_null8sm_true8sm  = 0;
    double chi2_gmin_null8sm_true8sm  = 0;
    double chi2_null_null8sm_true8Lee = 0;
    double chi2_gmin_null8sm_true8Lee = 0;
    // save Lee_strength at best fit
    // save fake spectrum

    double Lee_chi2_gmin_null8sm_true8sm = 0;
    double Lee_chi2_gmin_null8sm_true8Lee = 0;
    vector<double>vc_data;
    
    
    TFile *file_out = new TFile(TString::Format("file_out_%03d.root", ifile), "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("chi2_null_null8sm_true8sm", &chi2_null_null8sm_true8sm, "chi2_null_null8sm_true8sm/D" );
    tree->Branch("chi2_gmin_null8sm_true8sm", &chi2_gmin_null8sm_true8sm, "chi2_gmin_null8sm_true8sm/D" );
    tree->Branch("chi2_null_null8sm_true8Lee", &chi2_null_null8sm_true8Lee, "chi2_null_null8sm_true8Lee/D" );
    tree->Branch("chi2_gmin_null8sm_true8Lee", &chi2_gmin_null8sm_true8Lee, "chi2_gmin_null8sm_true8Lee/D" );
    
    tree->Branch("Lee_chi2_gmin_null8sm_true8sm", &Lee_chi2_gmin_null8sm_true8sm, "Lee_chi2_gmin_null8sm_true8sm/D" );
    tree->Branch("Lee_chi2_gmin_null8sm_true8Lee", &Lee_chi2_gmin_null8sm_true8Lee, "Lee_chi2_gmin_null8sm_true8Lee/D" );
    tree->Branch("vc_data", &vc_data);

    /////////
    Lee_test->scaleF_Lee = 0;//////////////////// true8sm
    Lee_test->Set_Collapse();
    Lee_test->Set_measured_data();

    Lee_test->Minimization_Lee_strength_FullCov(0, 1);// strength, fix_flag
    double chi2_H0SM_data = Lee_test->minimization_chi2;
    
    Lee_test->Minimization_Lee_strength_FullCov(1, 1);// strength, fix_flag
    double chi2_H1SM_data = Lee_test->minimization_chi2;

    double delta_chi2 = chi2_H1SM_data - chi2_H0SM_data;
    cout<<" ---> check delta_chi2: "<<delta_chi2<<endl;
    
    
    
    int N_toy = 1000;
    
    Lee_test->scaleF_Lee = 0;//////////////////// true8sm
    Lee_test->Set_Collapse();
    Lee_test->Set_Variations(N_toy);

    for(int itoy=1; itoy<=N_toy; itoy++) {
      if( itoy%max(N_toy/10,1)==0 )
	cout<<TString::Format(" ---> processing toy ( total cov ): %4.2f, %6d", itoy*1./N_toy, itoy)<<endl;

    
      int status_fit = 0;

      Lee_test->Set_fake_data_Variation(itoy);
      //Lee_test->Set_fake_data_Asimov();
      
      Lee_test->Minimization_Lee_strength_FullCov(1, 1);
      chi2_null_null8sm_true8sm = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;
      
      Lee_test->Minimization_Lee_strength_FullCov(0, 1);
      chi2_gmin_null8sm_true8sm = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;

      // Lee_test->scaleF_Lee = 1;//////////////////// true8Lee
      // Lee_test->Set_Collapse();
      // Lee_test->Set_Variations(1);
      // Lee_test->Set_fake_data_Variation(1);
      
      // Lee_test->Minimization_Lee_strength_FullCov(0, 1);
      // chi2_null_null8sm_true8Lee = Lee_test->minimization_chi2;
      // status_fit += Lee_test->minimization_status;
      
      // Lee_test->Minimization_Lee_strength_FullCov(1,0);
      // chi2_gmin_null8sm_true8Lee = Lee_test->minimization_chi2;
      // status_fit += Lee_test->minimization_status;

      if( status_fit!=0 ) continue;

      tree->Fill();
    
    }
    
    file_out->cd();
    tree->Write();
    file_out->Close();

  }


  

  
  if( 0 ) {
    // Lee_test->scaleF_Lee = 1; 
    // Lee_test->Set_Collapse();
    // Lee_test->Set_Variations( 2 );
    // Lee_test->Set_fake_data_Variation( 2 );
    // Lee_test->Minimization_Lee_strength_FullCov(2, 0);// scaleF_Lee, fixed or not

 
    double chi2_null_null8sm_true8sm  = 0;
    double chi2_gmin_null8sm_true8sm  = 0;
    double chi2_null_null8sm_true8Lee = 0;
    double chi2_gmin_null8sm_true8Lee = 0;
    // save Lee_strength at best fit
    // save fake spectrum

    double Lee_chi2_gmin_null8sm_true8sm = 0;
    double Lee_chi2_gmin_null8sm_true8Lee = 0;
    vector<double>vc_data;
    
    
    TFile *file_out = new TFile(TString::Format("file_out_%03d.root", ifile), "recreate");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("chi2_null_null8sm_true8sm", &chi2_null_null8sm_true8sm, "chi2_null_null8sm_true8sm/D" );
    tree->Branch("chi2_gmin_null8sm_true8sm", &chi2_gmin_null8sm_true8sm, "chi2_gmin_null8sm_true8sm/D" );
    tree->Branch("chi2_null_null8sm_true8Lee", &chi2_null_null8sm_true8Lee, "chi2_null_null8sm_true8Lee/D" );
    tree->Branch("chi2_gmin_null8sm_true8Lee", &chi2_gmin_null8sm_true8Lee, "chi2_gmin_null8sm_true8Lee/D" );
    
    tree->Branch("Lee_chi2_gmin_null8sm_true8sm", &Lee_chi2_gmin_null8sm_true8sm, "Lee_chi2_gmin_null8sm_true8sm/D" );
    tree->Branch("Lee_chi2_gmin_null8sm_true8Lee", &Lee_chi2_gmin_null8sm_true8Lee, "Lee_chi2_gmin_null8sm_true8Lee/D" );
    tree->Branch("vc_data", &vc_data);

    
    
    int N_toy = 100;    
    
    for(int itoy=1; itoy<=N_toy; itoy++) {
      if( itoy%(N_toy/10)==0 )
	cout<<TString::Format(" ---> processing toy ( total cov ): %4.2f, %6d", itoy*1./N_toy, itoy)<<endl;

      int status_fit = 0;

      Lee_test->scaleF_Lee = 0;//////////////////// true8sm
      Lee_test->Set_Collapse();
      Lee_test->Set_Variations(1);
      Lee_test->Set_fake_data_Variation(1);
      
      Lee_test->Minimization_Lee_strength_FullCov(0, 1);
      chi2_null_null8sm_true8sm = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;
      
      Lee_test->Minimization_Lee_strength_FullCov(1,0);
      chi2_gmin_null8sm_true8sm = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;

      Lee_test->scaleF_Lee = 1;//////////////////// true8Lee
      Lee_test->Set_Collapse();
      Lee_test->Set_Variations(1);
      Lee_test->Set_fake_data_Variation(1);
      
      Lee_test->Minimization_Lee_strength_FullCov(0, 1);
      chi2_null_null8sm_true8Lee = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;
      
      Lee_test->Minimization_Lee_strength_FullCov(1,0);
      chi2_gmin_null8sm_true8Lee = Lee_test->minimization_chi2;
      status_fit += Lee_test->minimization_status;

      if( status_fit!=0 ) continue;

      tree->Fill();
    }
    
    file_out->cd();
    tree->Write();
    file_out->Close();
    
  }

  ///////////////////////////////////////////////////////////
  
  //return 0;
  
}
