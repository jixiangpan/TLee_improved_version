#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include "stdlib.h"
using namespace std;

#include<map>

#include "TLee.h"

#include "TApplication.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////// MAIN //////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  TString roostr = "";
  
  double scaleF_POT = 1;
  int ifile = 0;
  
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
 
  bool flag_both_numuCC        = 0;// 1
  bool flag_CCpi0_FC_by_numuCC = 1;// 2


  ///////////////////////// gg
  
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

  ///////////////////////// gg
  
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
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  theApp.Run();
  
  return 0;
}
