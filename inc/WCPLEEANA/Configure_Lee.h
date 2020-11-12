namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes
  
  TString spectra_file = "./data_framework_Doc33131/merge_all.root";
  TString flux_Xs_directory = "./data_framework_Doc33131/flux_Xs/";
  TString detector_directory = "./data_framework_Doc33131/det_both/";
  TString mc_directory = "./data_framework_Doc33131/mc_stat/";
  
  /*
  TString spectra_file = "";
  TString flux_Xs_directory = "";
  TString detector_directory = "";
  TString mc_directory = "";
  */
  
  /*
  /// some places may need to be changed when use different file-formats
  void TLee::Set_Spectra_MatrixCov()
  (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  (*) map_Lee_ch                     -----> tag LEE channels
  (*) for(int ich=1; ich<=7; ich++)  -----> data channels
  (*) for(int idx=1; idx<=17; idx++) -----> flux_Xs files
  (*) map_detectorfile_str           -----> detector files
  (*) gbins_mc_stat                  -----> number of bins in mc_stat
  (*) for(int ifile=0; ifile<=99; ifile++) ----->  mc_stat files
  */
  
  ////////// display graphics flag

  bool flag_display_graphics = 0;
  
  ////////// systematics flag
  
  bool flag_syst_flux_Xs    = 1;
  bool flag_syst_detector   = 1;
  bool flag_syst_additional = 1;
  bool flag_syst_mc_stat    = 1;

  double Lee_strength_for_outputfile_covariance_matrix = 0;
  
  bool flag_plotting_systematics   = 0;
  
  ////////// goodness of fit
  
  double Lee_strength_for_GoF      = 0;
    
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
