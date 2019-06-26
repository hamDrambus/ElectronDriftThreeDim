void comp_I_anode(void) {
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  double Dtime_left = 2.4e-6, Dtime_right = 12e-6;
  int NN = Dtime_right/2e-8;
  TH1D* hist_current0 = new TH1D ("Induced current 0", "Induced current 0", NN, 0, Dtime_right);
  TH1D* hist_current1 = new TH1D ("Induced current 1", "Induced current 1", NN, 0, Dtime_right);
  TH1D* hist_current2 = new TH1D ("Induced current 2", "Induced current 2", NN, 0, Dtime_right);
  TH1D* hist_current3 = new TH1D ("Induced current 3", "Induced current 3", NN, 0, Dtime_right);

  int DEF_W = 900, DEF_H = 700;
  std::string fname0("Output/v17.1/eData_7.0Td_L1.8e-2m_VN.root");
  std::string fname1("Output/v17.2/eData_7.0Td_L1.8e-2m_VN.root");
  std::string fname2("Output/v17.3/eData_7.0Td_L1.8e-2m_VN.root");
  std::string fname3("Output/v17.4/eData_7.0Td_L1.8e-2m_VN.root");
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

  double max_val = 0;
  for (int nhist = 0; nhist<4; ++nhist) {
    TH1D* hist = NULL;
    std::string fname;
    switch (nhist) 
    {
		case 0: {
		    hist = hist_current0;
		    fname = fname0;
		    break;
		}
		case 1: {
		    hist = hist_current1;
		    fname = fname1;
		    break;
		}
		case 2: {
		    hist = hist_current2;
		    fname = fname2;
		    break;
		}
    case 3: {
		    hist = hist_current3;
		    fname = fname3;
		    break;
		}
    }
    TFile * file = new TFile (fname.c_str());
    if (0==file) {
      std::cout<<"File \""<<fname<<"\" not opened"<<std::endl;
      continue;
    }
    if (!file->IsOpen()){
      std::cout<<"File \""<<fname<<"\" not opened"<<std::endl;
      continue;
    }
    std::cout<<"Opened file: \""<<fname<<"\""<<std::endl;
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
	    if (time_start!=0.0) {
		    hist->Fill(0.5*(time_start+delta_time_full+prev_time), (pos_finish-prev_pos)/(time_start+delta_time_full-prev_time));
	    }
	    prev_time = time_start+delta_time_full;
	    prev_pos = pos_finish;
      for (int bin = 1, bin_end = hist->GetNbinsX()+1; bin!=bin_end; ++bin) {
		    max_val = std::max(max_val, (double) hist->GetBinContent(bin));
	    }
    }
  }
  max_val*=1.1;
	gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("");
  TCanvas *c0 = new TCanvas ("Induced_current", "Induced_current", DEF_W, DEF_H);
  //gPad->SetLogy();
  c0->SetGrid();
	c0->SetTicks();
  TLegend *legend = new TLegend( 0.55, 0.65, 0.9, 0.9);
	//legend->SetHeader("");
	legend->SetMargin(0.25);
	TH2F* frame = new TH2F("frame", "Induced current", 500, 0, Dtime_right, 500, 0, max_val);
	frame->GetXaxis()->SetTitle("t [s]");
	frame->GetYaxis()->SetTitle("arb.");
	frame->Draw();
	
	hist_current0->SetLineWidth(2);
	hist_current0->SetLineColor(kRed);//(kYellow-3);
	hist_current0->Draw("hist cpsame");
	hist_current1->SetLineWidth(2);
	hist_current1->SetLineColor(kBlack);//(kRed);
	hist_current1->Draw("hist cpsame");
	hist_current2->SetLineWidth(2);
	hist_current2->SetLineColor(kBlue);
	hist_current2->Draw("hist cpsame");
  hist_current3->SetLineWidth(2);
	hist_current3->SetLineColor(kGreen);
	hist_current3->Draw("hist cpsame");
	
	legend->AddEntry(hist_current0, (std::string("7.0 Td, 18mm, XS_DA = 1e-18 cm^{2}")).c_str(), "l");
	legend->AddEntry(hist_current1, (std::string("7.0 Td, 18mm, XS_DA = 3e-18 cm^{2}")).c_str(), "l");
	legend->AddEntry(hist_current2, (std::string("7.0 Td, 18mm, XS_DA = 6e-18 cm^{2}")).c_str(), "l");
	legend->AddEntry(hist_current3, (std::string("7.0 Td, 18mm, XS_DA = 1e-17 cm^{2}")).c_str(), "l");
	//legend->AddEntry(histE_4, (std::string("7.0 Td, W=3e-2 eV, Loss=150x")).c_str(), "l");
	
	frame->Draw("sameaxis");
	legend->Draw("same");
	c0->Update();
}

 


