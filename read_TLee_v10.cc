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

///////////////////////////////////////////////////////// ccc

void TLee::Exe_Goodness_of_fit(int num_Y, int num_X, TMatrixD matrix_pred, TMatrixD matrix_data, TMatrixD matrix_syst, int index)
{
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

  }
  /////////////////////////////////////////////////////////////////////////////////////////////
  else {
    
  }
  
}
  
///////////////////////////////////////////////////////// ccc

void TLee::Set_Collapse()
{
  //////////////////////////////////////// pred

  map_pred_spectrum_newworld_bin.clear();
  TMatrixD matrix_pred_oldworld(1, bins_oldworld);
  for(int ibin=0; ibin<bins_oldworld; ibin++) matrix_pred_oldworld(0, ibin) = map_input_spectrum_oldworld_bin[ibin];

  matrix_pred_newworld.Clear();
  matrix_pred_newworld.ResizeTo(1, bins_newworld);
  matrix_pred_newworld = matrix_pred_oldworld * matrix_transform;
  if( bins_newworld!=matrix_pred_newworld.GetNcols() ) { cerr<<"bins_newworld!=matrix_pred_newworld.GetNcols()"<<endl; exit(1); }
  for(int ibin=0; ibin<bins_newworld; ibin++) map_pred_spectrum_newworld_bin[ibin] = matrix_pred_newworld(0, ibin);
  
  ////////////////////////////////////////
  
  matrix_absolute_cov_oldworld.Clear();
  matrix_absolute_cov_oldworld.ResizeTo( bins_oldworld, bins_oldworld );

  if( flag_syst_flux_Xs ) matrix_absolute_cov_oldworld += matrix_input_cov_flux_Xs;
  if( flag_syst_detector ) matrix_absolute_cov_oldworld += matrix_input_cov_detector;
  if( flag_syst_additional ) matrix_absolute_cov_oldworld += matrix_input_cov_additional;

  TMatrixD matrix_transform_T( bins_newworld, bins_oldworld );
  matrix_transform_T.Transpose( matrix_transform );

  matrix_absolute_cov_newworld.Clear();
  matrix_absolute_cov_newworld.ResizeTo(bins_newworld, bins_newworld);
  matrix_absolute_cov_newworld = matrix_transform_T * matrix_absolute_cov_oldworld * matrix_transform;

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
    for(int idx=0; idx<ibin; idx++) {
      double x(0), y(0);
      gh_mc_stat_bin[ibin]->GetPoint(idx, x, y);
      gh_mc_stat_bin[ibin]->SetPoint(idx, x, y*scaleF_POT2);
    }// ipoint
  }// ibin
  
}

///////////////////////////////////////////////////////// ccc

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
  
  roostr = "./data_framework/merge_all.root";
  TFile *file_spectra = new TFile(roostr, "read");

  ///
  TMatrixD *mat_collapse = (TMatrixD*)file_spectra->Get("mat_collapse");
  matrix_transform.Clear();
  matrix_transform.ResizeTo( mat_collapse->GetNrows(), mat_collapse->GetNcols() );
  matrix_transform = (*mat_collapse);

  ///
  for(int ich=1; ich<=16; ich++) {
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
    roostr = TString::Format("./data_framework/flux_Xs/cov_%d.root", idx);
    map_file_flux_Xs_frac[idx] = new TFile(roostr, "read");
    map_matrix_flux_Xs_frac[idx] = (TMatrixD*)map_file_flux_Xs_frac[idx]->Get(TString::Format("frac_cov_xf_mat_%d", idx));
    // cout<<TString::Format(" ---> check: flux and Xs, %2d  ", idx)<<roostr<<endl;
    matrix_flux_Xs_frac += (*map_matrix_flux_Xs_frac[idx]);
  }
  // cout<<endl;
  
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
    roostr = TString::Format("./data_framework/mc_stat/%d.log", ifile);
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

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


//
// usage: root -l read_TLee.cc+(1, 1)
//

//int main(int argc, char** argv)
void read_TLee_v10(double scalePOT, int nfile)
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
  
  TLee *Lee_test = new TLee();

  ////////// just do it one time in the whole procedure
  Lee_test->scaleF_POT = scaleF_POT;
  Lee_test->Set_Spectra_MatrixCov();
  Lee_test->Set_POT_implement();
  Lee_test->Set_TransformMatrix();

  ////////// can do any times
  Lee_test->flag_syst_flux_Xs    = 1;
  Lee_test->flag_syst_detector   = 1;
  Lee_test->flag_syst_additional = 1;
  Lee_test->flag_syst_mc_stat    = 1;
  Lee_test->scaleF_Lee           = 0;
  Lee_test->Set_Collapse();
  
  //////////////////////////////////////////////////////////////////////////////////////// Goodness of fit
  
  Lee_test->scaleF_Lee = 0;
  Lee_test->Set_Collapse();
 
  bool flag_both_numuCC = 1;
  



  
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
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  // return 0;
}
