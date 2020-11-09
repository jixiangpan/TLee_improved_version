#ifndef TLee_ana
#define TLee_ana

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

/// minuit2
#include "Math/Functor.h"
#include "Minuit2/Minuit2Minimizer.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//#include "draw.icc"

///////////////////////////////////////////////////////////////////////////////////////////////////////// TLee

class TLee {
public:
  TLee() {
    rand = new TRandom3(0);
    flag_individual_cov_newworld = true;
  }

  /////////////////////////////////////////////////////// data memeber

  TRandom3 *rand;
  
  double scaleF_POT;
  double scaleF_Lee;

  bool flag_syst_flux_Xs;
  bool flag_syst_detector;
  bool flag_syst_additional;
  bool flag_syst_mc_stat;

  TString spectra_file;
  TString flux_Xs_directory;
  TString detector_directory;
  TString mc_directory;
  
  /// bin index starts from 0, channel index from 1
  map<int, map<int, double> >map_input_spectrum_ch_bin;
  map<int, TString>map_input_spectrum_ch_str;
  map<int, double>map_input_spectrum_oldworld_bin;
  int bins_oldworld;
  int bins_newworld;

  map<int, int>map_Lee_ch;
  map<int, int>map_Lee_oldworld;

  map<int, map<int, double> >map_data_spectrum_ch_bin;
  map<int, double>map_data_spectrum_newworld_bin;
  TMatrixD matrix_data_newworld;
  
  TMatrixD matrix_input_cov_flux_Xs;  
  TMatrixD matrix_input_cov_flux;
  TMatrixD matrix_input_cov_Xs;
  TMatrixD matrix_input_cov_detector;
  map<int, TMatrixD>matrix_input_cov_detector_sub;
  TMatrixD matrix_input_cov_additional;

  map<int, TGraph*>gh_mc_stat_bin;

  TMatrixD matrix_transform;

  map<int, double>map_pred_spectrum_newworld_bin;
  TMatrixD matrix_pred_newworld;
  
  TMatrixD matrix_absolute_cov_oldworld;
  TMatrixD matrix_absolute_cov_newworld;

  bool flag_individual_cov_newworld;
  TMatrixD matrix_absolute_flux_cov_newworld;
  TMatrixD matrix_absolute_Xs_cov_newworld;
  TMatrixD matrix_absolute_detector_cov_newworld;
  map<int, TMatrixD>matrix_absolute_detector_sub_cov_newworld;
  TMatrixD matrix_absolute_mc_stat_cov_newworld;
  TMatrixD matrix_absolute_additional_cov_newworld;

  map<int, map<int, double> >map_toy_variation;
  
  map<int, double>map_fake_data;
  
  int minimization_status;
  double minimization_chi2;
  double minimization_Lee_strength_val;
  double minimization_Lee_strength_err;

  /////////////////////////////////////////////////////// function member

  void Set_config_file_directory(TString spectra_file_, TString flux_Xs_directory_, TString detector_directory_, TString mc_directory_);
  
  void Set_Spectra_MatrixCov();
  void Set_POT_implement();
  void Set_TransformMatrix();

  void Set_Collapse();

  // Y constrained by X
  void Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst, int index);

  void Set_Variations(int num_toy);
  
  void Set_toy_Asimov();
  void Set_toy_Variation(int itoy);
  void Set_measured_data();
  void Set_fakedata(TMatrixD matrix_fakedata);// 1, bins
  
  // minimization
  void Minimization_Lee_strength_FullCov(double Lee_initial_value, bool flag_fixed);
  
  // Feldman-Cousins approach
  void Exe_Feldman_Cousins(double Lee_true_low, double Lee_true_hgh, double step, int num_toy, int ifile);
  void Exe_Fledman_Cousins_Asimov(double Lee_true_low, double Lee_true_hgh, double step);
  void Exe_Fiedman_Cousins_Data(TMatrixD matrix_fakedata, double Lee_true_low, double Lee_true_hgh, double step);
};

#endif
