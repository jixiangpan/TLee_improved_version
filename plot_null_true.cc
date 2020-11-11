#include "./src/draw.icc"

void plot_null_true()
{  
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
  
  TString roostr = "";

  int color_fill = kRed;
  
  roostr = "file_null8Lee_true8Lee.root";
  roostr = "file_null8Lee_true8sm.root";
  roostr = "file_null8sm_Lee2sm.root";
  
  TFile *file_out = new TFile(roostr, "read");
  TTree *tree = (TTree*)file_out->Get("tree");

  // Declaration of leaf types
  Double_t        chi2_null_null8sm_true8sm;
  Double_t        chi2_gmin_null8sm_true8sm;
  Double_t        chi2_null_null8Lee_true8Lee;
  Double_t        chi2_gmin_null8Lee_true8Lee;

  // List of branches
  TBranch        *b_chi2_null_null8sm_true8sm;   //!
  TBranch        *b_chi2_gmin_null8sm_true8sm;   //!
  TBranch        *b_chi2_null_null8Lee_true8Lee;   //!
  TBranch        *b_chi2_gmin_null8Lee_true8Lee;   //!

  // Set branch addresses and branch pointers  
  tree->SetBranchAddress("chi2_null_null8sm_true8sm", &chi2_null_null8sm_true8sm, &b_chi2_null_null8sm_true8sm);
  tree->SetBranchAddress("chi2_gmin_null8sm_true8sm", &chi2_gmin_null8sm_true8sm, &b_chi2_gmin_null8sm_true8sm);
  tree->SetBranchAddress("chi2_null_null8Lee_true8Lee", &chi2_null_null8Lee_true8Lee, &b_chi2_null_null8Lee_true8Lee);
  tree->SetBranchAddress("chi2_gmin_null8Lee_true8Lee", &chi2_gmin_null8Lee_true8Lee, &b_chi2_gmin_null8Lee_true8Lee);

  int entries = tree->GetEntries();
  
  cout<<endl<<" ---> entries "<<entries<<endl<<endl;

  vector<double>vc_dchi2_null8Lee_true8Lee;
  
  for(int ientry=0; ientry<entries; ientry++) {
    if( ientry%max(entries/10,1)==0 ) cout<<TString::Format(" ---> processing : %4.2f, %6d", ientry*1./entries, ientry)<<endl;
    tree->GetEntry(ientry);

    double dchi2_null8Lee_true8Lee = chi2_null_null8Lee_true8Lee - chi2_gmin_null8Lee_true8Lee;
    vc_dchi2_null8Lee_true8Lee.push_back( dchi2_null8Lee_true8Lee );
  }

  int size_vc = vc_dchi2_null8Lee_true8Lee.size();
  sort( vc_dchi2_null8Lee_true8Lee.begin(), vc_dchi2_null8Lee_true8Lee.end() );

  double dchi2_data = 5.35;
  int line_eff = 0;

  for(int idx=0; idx<size_vc; idx++) {
    double val_dchi2 = vc_dchi2_null8Lee_true8Lee.at(idx);
    if( val_dchi2>dchi2_data ) line_eff++;      
  }

  double pValue = line_eff*1./size_vc;// double side pValue
  double val_sigma = sqrt( TMath::ChisquareQuantile( 1-pValue, 1 ) );
  val_sigma = sqrt( TMath::ChisquareQuantile( 1-pValue*2, 1 ) );

  
  cout<<endl<<TString::Format(" ---> pValue %6.4f, #sigma %4.2f, chi2 %4.2f", pValue, val_sigma, val_sigma*val_sigma)<<endl<<endl;
  
  ////////////////////////////////////////////

  TH1D *h1_dchi2 = new TH1D("h1_dchi2", "", 100, vc_dchi2_null8Lee_true8Lee.at(int(size_vc*0.01)), vc_dchi2_null8Lee_true8Lee.at(int(size_vc*1)-1) + 5 );

  for(int idx=0; idx<size_vc; idx++) {
    double val_dchi2 = vc_dchi2_null8Lee_true8Lee.at(idx);
    h1_dchi2->Fill( val_dchi2 );
  }

  h1_dchi2->Scale( 1./size_vc );

  TCanvas *canv_h1_dchi2 = new TCanvas("canv_h1_dchi2", "canv_h1_dchi2", 900, 650);
  //canv_h1_dchi2->SetLogy();
  canv_h1_dchi2->SetLeftMargin(0.15); canv_h1_dchi2->SetRightMargin(0.1);
  canv_h1_dchi2->SetTopMargin(0.1); canv_h1_dchi2->SetBottomMargin(0.15);    
  h1_dchi2->Draw("hist f");
  h1_dchi2->GetYaxis()->SetTitle("PDF"); h1_dchi2->GetXaxis()->SetTitle("#Delta#chi^{2} = #chi^{2}_{LEE} - #chi^{2}_{SM}");
  h1_dchi2->GetXaxis()->SetLabelSize(0.05); h1_dchi2->GetXaxis()->SetTitleSize(0.05);
  h1_dchi2->GetYaxis()->SetLabelSize(0.05); h1_dchi2->GetYaxis()->SetTitleSize(0.05);
  h1_dchi2->GetXaxis()->CenterTitle(); h1_dchi2->GetYaxis()->CenterTitle();
  h1_dchi2->GetXaxis()->SetTitleOffset(1.2);
  h1_dchi2->GetYaxis()->SetNdivisions(508);
  
  h1_dchi2->SetLineColor(color_fill);
  h1_dchi2->SetFillStyle(3004);
  h1_dchi2->SetFillColor(color_fill);
  
  canv_h1_dchi2->cd(); canv_h1_dchi2->Update();
  double x1 = gPad->GetUxmin(); double y1 = gPad->GetUymin();
  double x2 = gPad->GetUxmax(); double y2 = gPad->GetUymax();

  TLine *lineA_dchi2at1 = new TLine(dchi2_data, 0, dchi2_data, y2);    
  lineA_dchi2at1->Draw("same");
  lineA_dchi2at1->SetLineWidth(2);
  lineA_dchi2at1->SetLineColor(kBlack);
  lineA_dchi2at1->SetLineStyle(7);
  
  // auto *tt_text_data = new TLatex( dchi2_data+0.5, y2*0.5, Form("p-value = %5.3f #rightarrow %3.2f#sigma", pValue, val_sigma) );
  // tt_text_data->SetTextAlign(11); tt_text_data->SetTextSize(0.05); tt_text_data->SetTextAngle(0);
  // tt_text_data->SetTextFont(42);  tt_text_data->Draw(); tt_text_data->SetTextColor(kBlack);
  
  
  //TPaveText *pt = new TPaveText( dchi2_data+0.5, y2*0.35, dchi2_data+0.5, y2*0.6,"l");
  TPaveText *pt = new TPaveText( x2*0.3, y2*0.35, x2*0.3, y2*0.6,"l");
  
  pt->SetTextSize(0.05);
  pt->SetTextFont(42);
  pt->SetTextAlign(11);
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->AddText( Form("#Delta#chi^{2} = %4.2f", dchi2_data) );
  //((TText*)pt->GetListOfLines()->Last())->SetTextColor(kBlack);
  pt->AddText( Form("#rightarrow p-value = %5.3f", pValue) );
  pt->AddText( Form("#rightarrow %3.2f#sigma", val_sigma) );
  
  pt->Draw();

  h1_dchi2->Draw("same axis");
  //canv_h1_dchi2->SaveAs("canv_h1_null8Lee_true8Lee.png");
  //canv_h1_dchi2->SaveAs("canv_h1_null8Lee_true8sm.png");
  canv_h1_dchi2->SaveAs("canv_h1_true8sm_Lee2sm.png");
  
  ////////////////////////////////////////////
  
  
}

