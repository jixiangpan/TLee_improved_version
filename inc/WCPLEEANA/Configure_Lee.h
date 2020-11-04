namespace config_Lee
{
  ////////// 
  
  TString spectra_file = "./data_framework/merge_all_new_weights.root";
  TString flux_Xs_directory = "./data_framework/flux_Xs/";
  TString detector_directory = "./data_framework/det/";
  TString mc_directory = "./data_framework/mc_stat/";

  /*
    void TLee::Set_Spectra_MatrixCov()
    (*) map_input_spectrum_ch_str
    (*) map_Lee_ch
    (*) for(int ich=1; ich<=7; ich++), data
    (*) for(int idx=1; idx<=17; idx++), flux_Xs
    (*) map_detectorfile_str
    (*) int gbins_mc_stat = 137;
    (*) for(int ifile=0; ifile<=99; ifile++), mc_sat
    
    void TLee::Set_Collapse()
    (*) double val_mc_stat_cov = gh_mc_stat_bin[ibin]->Eval( scaleF_Lee ); mc_stat
  */
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs    = 1;
  bool flag_syst_detector   = 1;
  bool flag_syst_additional = 1;
  bool flag_syst_mc_stat    = 1;

  double Lee_strength_for_output_covariance_matrix = 0;
  
  ////////// goodness of fit
 
  bool flag_both_numuCC            = 0;// 1
  bool flag_CCpi0_FC_by_numuCC     = 0;// 2
  bool flag_CCpi0_PC_by_numuCC     = 0;// 3
  bool flag_NCpi0_by_numuCC        = 0;// 4
  bool flag_nueCC_PC_by_numuCC_pi0 = 0;// 5
  bool flag_nueCC_HghE_FC_by_numuCC_pi0_nueFC = 0;// 6, HghE>800 MeV
  bool flag_nueCC_LowE_FC_by_all   = 0;// 7
  bool flag_nueCC_FC_by_all        = 0;// 8

  ////////// Lee strength fitting -- data

  bool flag_Lee_strength_data = 0;

  //////////
}
