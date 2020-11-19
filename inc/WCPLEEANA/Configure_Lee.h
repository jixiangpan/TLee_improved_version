namespace config_Lee
{
  ////////// input files for spectra and covariance matrixes
  
  /* TString spectra_file = "./data_framework_Doc33131/merge_all.root"; */
  /* TString flux_Xs_directory = "./data_framework_Doc33131/flux_Xs/"; */
  /* TString detector_directory = "./data_framework_Doc33131/det_both/"; */
  /* TString mc_directory = "./data_framework_Doc33131/mc_stat/"; */

  TString spectra_file = "./data_framework_fake7/merge_fake7.root";
  TString flux_Xs_directory = "./data_framework_fake7/XsFlux/";
  TString detector_directory = "./data_framework_fake7/DetVar_11stat/";
  TString mc_directory = "./data_framework_fake7/mc_stat/";

  int channels_observation = 7;// data channels (=hdata_obsch_# in spectra_file above)
                               // which is equal to the channels after collapse

  int syst_cov_flux_Xs_begin = 1;// files in flux_Xs_directory above
  int syst_cov_flux_Xs_end   = 17;
 
  int syst_cov_mc_stat_begin = 0;// files in mc_directory above
  int syst_cov_mc_stat_end   = 99;
   
  /*
  /// some places may need to be changed when use different file-formats
  void TLee::Set_Spectra_MatrixCov()
  (*) map_input_spectrum_ch_str      -----> prediction channels before collapse
  (*) map_Lee_ch                     -----> tag LEE channels
  (*) map_detectorfile_str           -----> detector files
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
