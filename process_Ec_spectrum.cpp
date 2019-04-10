 {
	int DEF_W = 900, DEF_H = 700;
	double EN_MIN_=0;
	double EN_MAX_=14;
	int NN = 600;
	TH1D* histE_0 = new TH1D ("EnergyC 7.0 Td [eV]","EnergyC 7.0 Td [eV]",NN,EN_MIN_, EN_MAX_);
	TH1D* histE_1 = new TH1D ("EnergyC 7.0 Td [eV] modified","EnergyC 7.0 Td [eV] modified",NN,EN_MIN_, EN_MAX_);
	histE_0->SetStats(false);
	histE_1->SetStats(false);
	//gauss peak
	double En_artifical_peak = 11.103;
	double sigma_artifical_peak = 0.005;
	double rel_area_artifical_peak = 0.05;	
	
	std::string fname_0("Output/v12.1/eData_7.0Td");
	std::string out_fname0 ("Output/v12.1/Ec_forms/7.0Td.dat"); //Make sure this folder exists before launching this script!
	std::string out_fname0_1 ("Output/v12.1/Ec_forms/7.0Td_mod.dat");
	
	double En_collision;
	double delta_time;

	double max_val = 0;
	for (int n_hist =0; n_hist<1; ++n_hist) {
	    TH1D* histE = NULL;
		TFile * file = 0;
		std::string filename = fname_0;
		filename= filename+".root";
		file = new TFile (filename.c_str());
		if (0==file ){
		    std::cout<<"File not opened"<<std::endl;
		    break;
		}
		if (!file->IsOpen()){
		    std::cout<<"File not opened"<<std::endl;
		    break;
		}
		std::cout<<"Opened "<<filename<<std::endl;
		TTree * tree = (TTree*) file->Get("ElectronHistory");

		tree->SetBranchAddress("energy_coll", &En_collision);
		tree->SetBranchAddress("time_delta", &delta_time);
		
		unsigned long int _end_ = tree->GetEntries();
		for (unsigned long int i=0;i!=_end_;++i){
		    tree->GetEntry(i);
		    //histE->Fill(std::fabs(En_collision), delta_time);
		    histE_0->Fill(std::fabs(En_collision));
		}
	    
	    //normalization of spectrum;
	    double Norm =0;
	    for (int bin = 1, bin_end = histE_0->GetNbinsX()+1; bin!=bin_end; ++bin) {
			Norm+=histE_0->GetBinContent(bin)*histE_0->GetBinWidth(bin);
	    }
		double Norm2 = 1 + rel_area_artifical_peak;
		double amplitude = rel_area_artifical_peak/sigma_artifical_peak/sqrt(2*M_PI);
	    for (int bin = 1, bin_end = histE_0->GetNbinsX()+1; bin!=bin_end; ++bin) {
			double En = histE_0->GetBinCenter(bin);
			double val = histE_0->GetBinContent(bin)/Norm;
			histE_0->SetBinContent(bin, val);
			double extra = amplitude*exp(-0.5*pow( (En - En_artifical_peak)/sigma_artifical_peak , 2));
			extra = std::max(extra, 0.0);			
			histE_1->SetBinContent(bin, (val + extra)/Norm2);
			max_val = std::max(max_val, (double) histE_0->GetBinContent(bin));
			max_val = std::max(max_val, (double) histE_1->GetBinContent(bin));
	    }
	}
	max_val*=1.1;
	gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("");
	TCanvas *c_ = new TCanvas ("Ec spectra 3D", "Ec spectra 3D", DEF_W, DEF_H);
	c_->SetGrid();
	c_->SetTicks();
	//gPad->SetLogy();//Log Y
	TH2F* frame = new TH2F( "frame", "Collision E spectra 3D", 500, EN_MIN_, EN_MAX_, 500, 0, max_val);
	frame->GetXaxis()->SetTitle("Ee [eV]");
	frame->GetYaxis()->SetTitle("");
	frame->Draw();
	
	histE_0->SetLineWidth(2);
	histE_0->SetLineColor(kRed);//(kYellow-3);
	histE_0->Draw("hist cpsame");
	
	histE_1->SetLineWidth(2);
	histE_1->SetLineColor(kBlack);//(kYellow-3);
	histE_1->Draw("hist cpsame");

	std::ofstream str0, str1;
	str0.open(out_fname0, std::ios_base::trunc);
	str1.open(out_fname0_1, std::ios_base::trunc);
	str0<<"//\""<<fname_0+".root\" data"<<std::endl;
	str0<<"//E[eV]\tSpectrum value"<<std::endl;
	str1<<"//\""<<fname_0+".root\" data"<<std::endl;
	str1<<"//E[eV]\tModified spectrum value"<<std::endl;
	for (int bin = 1, bin_end = histE_0->GetNbinsX()+1; bin!=bin_end; ++bin) {
		double En = histE_0->GetBinCenter(bin);
		double val0 = histE_0->GetBinContent(bin);
		double val1 = histE_1->GetBinContent(bin);
		str0<<En<<"\t"<<val0<<std::endl;
		str1<<En<<"\t"<<val1<<std::endl;
	}
	str0.close();
	str1.close();
	frame->Draw("sameaxis");
	c_->Update();
}



