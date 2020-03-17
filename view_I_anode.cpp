Double_t FittingF_GE (Double_t *x, Double_t *par) {
  //par[0] - gauss ampl
  //par[1] - gauss mean
  //par[2] - gause sigmna
  //par[3] - exp start t
  //par[4] - exp ampl
  //par[5] - exp decay time
  return par[0]*std::exp(-0.5*std::pow((x[0]-par[1])/par[2], 2))
      + (x[0] > par[3] ? par[4]*std::exp(-1.0*x[0]/par[5]) : 0);
}

void view_I_anode(void) {
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  TH1D* hist_V_drift = new TH1D ("Drift velocity [m/s]","Drift velocity [m/s]",300,0, 2e4);
  double Dtime_left = 2.4e-6, Dtime_right = 4.0e-5;
  TH1D* hist_T_drift = new TH1D ("Drift time [s]","Drift time [s]",160, Dtime_left, Dtime_right);
  TH1D* hist_current = new TH1D ("Induced current", "Induced current", Dtime_right/2e-8, 0, Dtime_right);

  int DEF_W = 900, DEF_H = 700;
  std::string fname1("Output/v17.5/eData_7.0Td_L1.8e-2m_VN.root");
  double En_start;
  double En_collision;
  double En_finish;
  double En_avr;
  double pos_start;
  double pos_finish;
  double delta_x;
  double time_start;
  double delta_time;
  double delta_time_full;

  double theta_start;
  double theta_collision;
  double theta_finish;
  double delta_theta;
  double delta_l;
  short process; //enum ProcessType : short {None = 0, Elastic = 1, Resonance = 2, Overflow = 3};

  long int num_of_elastic = 0, num_of_resonance=0, num_of_overflow = 0, num_of_electrons =0;
  long int num_of_ions = 0, num_of_excitations=0, num_of_disso_attach = 0, num_of_NBrS = 0;

  TFile * file = new TFile (fname1.c_str());
  if (0==file) {
    std::cout<<"File not opened"<<std::endl;
    return;
  }
  if (!file->IsOpen()){
    std::cout<<"File not opened"<<std::endl;
    return;
  }
  std::cout<<"Opened file: "<<fname1<<std::endl;
  TTree * tree = (TTree*) file->Get("ElectronHistory");

  tree->SetBranchAddress("energy_initial", &En_start);
  tree->SetBranchAddress("energy_final", &En_finish);
  tree->SetBranchAddress("energy_coll", &En_collision);
  tree->SetBranchAddress("energy_average", &En_avr);

  tree->SetBranchAddress("time_initial", &time_start);
  tree->SetBranchAddress("time_delta", &delta_time);
  tree->SetBranchAddress("time_delta_full", &delta_time_full);
  tree->SetBranchAddress("process_type", &process);

  tree->SetBranchAddress("theta_initial", &theta_start);
  tree->SetBranchAddress("theta_coll", &theta_collision);
  tree->SetBranchAddress("theta_final", &theta_finish);
  tree->SetBranchAddress("theta_delta", &delta_theta);
  tree->SetBranchAddress("path_delta", &delta_l);

  tree->SetBranchAddress("position_initial",&pos_start);
  tree->SetBranchAddress("position_final",&pos_finish);
  tree->SetBranchAddress("position_delta",&delta_x);

  unsigned long int _end_ = tree->GetEntries();
  double prev_time = 0, prev_pos = 0; 	
  for (unsigned long int i=0; i!=_end_; ++i) {
    tree->GetEntry(i);
	if (time_start==0.0) {
		prev_time = time_start+delta_time_full;
		prev_pos = pos_finish;
	} else {
		hist_current->Fill(0.5*(time_start+delta_time_full+prev_time), (pos_finish-prev_pos)/(time_start+delta_time_full-prev_time));
	}
	prev_time = time_start+delta_time_full;
	prev_pos = pos_finish;
	if (((i==(_end_-1))||time_start==0.0)&&i!=0) {
        ++num_of_electrons;
		if (i!=(_end_-1))
        	tree->GetEntry(i-1);
        hist_V_drift->Fill(pos_finish/(time_start+delta_time_full));
        hist_T_drift->Fill(time_start+delta_time_full);
    }
  }
  TCanvas *c0 = new TCanvas ("Drift time", "Drift time", DEF_W, DEF_H);
  //gPad->SetLogy();
  TF1 *ff = new TF1("fit", FittingF_GE, Dtime_left, Dtime_right, 6);
  ff->SetParNames("Gaus amplitude", "Gaus mean", "Gaus sigma", "Slow t0[s]", "Slow amplitude", "Slow tau[s]");
  ff->SetParLimits(0, 100, 3000);
  ff->SetParLimits(1, 1.55e-6, 1.65e-6);
  ff->SetParLimits(2, 5e-8, 6e-7);
  ff->SetParLimits(3, 1.7e-6, 1.775e-6);
  ff->SetParLimits(4, 100, 2e4);
  ff->SetParLimits(5, 3e-7, 8e-7);
  hist_T_drift->Draw();
  c0->Update();
  //hist_T_drift->Fit(ff);
  TLine* line = new TLine();
  line->SetX1(ff->GetParameter(3));
  line->SetX2(ff->GetParameter(3));
  line->SetY1(c0->GetUymin());
  line->SetY2(c0->GetUymax());
  line->SetLineColor(kRed);
  //hist_T_drift->Draw();
  //line->Draw("same");
  //ff->Draw("same");
  TCanvas *c1 = new TCanvas ("Drift velocity", "Drift velocity", DEF_W, DEF_H);
  hist_V_drift->Draw();
  TCanvas *c2 = new TCanvas ("Induced_current", "Induced_current", DEF_W, DEF_H);
  hist_current->Draw("hist");
}

 


