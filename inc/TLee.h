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
//#include "Math/Functor.h"
//#include "Minuit2/Minuit2Minimizer.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

//#include "draw.icc"

///////////////////////////////////////////////////////////////////////////////////////////////////////// TLee
/*  
  Set_Spectra_MatrixCov():
      map_Lee_ch[8] = 1;
      map_Lee_ch[9] = 1;  
      if( map_Lee_ch.find(ich)!=map_Lee_ch.end() ) map_Lee_oldworld[index_oldworld] = 1;
*/


class TLee {
public:
  TLee() {
    rand = new TRandom3(0);
  }

  /////////////////////////////////////////////////////// data memeber

  TRandom3 *rand;
  
  double scaleF_POT;
  double scaleF_Lee;

  bool flag_syst_flux_Xs;
  bool flag_syst_detector;
  bool flag_syst_additional;
  bool flag_syst_mc_stat;
  
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
  TMatrixD matrix_input_cov_detector;
  TMatrixD matrix_input_cov_additional;

  map<int, TGraph*>gh_mc_stat_bin;

  TMatrixD matrix_transform;

  map<int, double>map_pred_spectrum_newworld_bin;
  TMatrixD matrix_pred_newworld;
  
  TMatrixD matrix_absolute_cov_oldworld;
  TMatrixD matrix_absolute_cov_newworld;
  
  /////////////////////////////////////////////////////// function member

  void Set_Spectra_MatrixCov();
  void Set_POT_implement();
  void Set_TransformMatrix();

  void Set_Collapse();

  // Y constrained by X
  void Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst, int index);
};

#endif
