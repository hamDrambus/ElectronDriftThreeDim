 {
	double EN_LOSS_MAX_=0.0007;
	double EN_MAX_=0.25;
	int NN = 300;
	TH1D* hist_dE_0 = new TH1D ("Energy increase in resonance 1dim [eV]","Energy increase in resonance 1dim [eV]",NN,-EN_MAX_, EN_MAX_);
	TH1D* hist_dE_1 = new TH1D ("Energy increase in resonance 3dim [eV]","Energy increase in resonance 3dim [eV]",NN,-EN_MAX_, EN_MAX_);
	TH1D* hist_dE_2 = new TH1D ("Energy increase in resonance 3dim uniform [eV]","Energy increase in resonance 3dim uniform [eV]",NN,-EN_MAX_, EN_MAX_);
	TH1D* hist_dE_abs_0 = new TH1D ("Energy abs increase in resonance 1dim [eV]","Energy abs increase in resonance 1dim [eV]",NN,0.0, EN_MAX_);
	TH1D* hist_dE_abs_1 = new TH1D ("Energy abs increase in resonance 3dim [eV]","Energy abs increase in resonance 3dim [eV]",NN,0.0, EN_MAX_);
	TH1D* hist_dE_abs_2 = new TH1D ("Energy abs increase in resonance 3dim uniform [eV]","Energy abs increase in resonance 3dim uniform [eV]",NN,0.0, EN_MAX_);
	TH1D* hist_Eloss_0 = new TH1D ("Energy loss in resonance 1dim [eV]","Energy loss in resonance 1dim [eV]",NN,0, EN_LOSS_MAX_);
	TH1D* hist_Eloss_1 = new TH1D ("Energy loss in resonance 3dim [eV]","Energy loss in resonance 3dim [eV]",NN,0, EN_LOSS_MAX_);
	TH1D* hist_Eloss_2 = new TH1D ("Energy loss in resonance 3dim uniform [eV]","Energy loss in resonance 3dim unifrom [eV]",NN,0, EN_LOSS_MAX_);
	hist_dE_0->SetStats(false);
	hist_dE_1->SetStats(false);
	hist_dE_2->SetStats(false);
	hist_dE_abs_0->SetStats(false);
	hist_dE_abs_1->SetStats(false);
	hist_dE_abs_2->SetStats(false);
	hist_Eloss_0->SetStats(false);
	hist_Eloss_1->SetStats(false);
	hist_Eloss_2->SetStats(false);

	int DEF_W = 900, DEF_H = 700;
	std::vector<double> DRIFT_DISTANCE = {3e-3, 3e-3, 3e-3, 3e-3};
	std::string fname_0("../OneDimSimulation/Output/v13.1/eData_7.0Td");
	std::string fname_1("Output/v08.1/eData_7.0Td");
	std::string fname_2("Output/v11.1/eData_7.0Td");
	
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
	short process; //enum ProcessType : short {None = 0, Elastic = 1, Resonance = 2, Overflow = 3};

	double max_val_dE = 0;
	double max_val_dE_abs = 0;
	double max_val_Eloss = 0;
	for (int nhist = 0; nhist<3;++nhist) {
	    TH1D* hist_dE = NULL;
	    TH1D* hist_dE_abs = NULL;
	    TH1D* hist_Eloss = NULL;
	    std::string fname;
	    switch (nhist) 
	    {
		case 0: {
		    hist_dE = hist_dE_0;
		    hist_dE_abs = hist_dE_abs_0;
		    hist_Eloss = hist_Eloss_0;
		    fname = fname_0;
		    break;
		}
		case 1: {
		    hist_dE = hist_dE_1;
		    hist_dE_abs = hist_dE_abs_1;
		    hist_Eloss = hist_Eloss_1;
		    fname = fname_1;
		    break;
		}
		case 2: {
		    hist_dE = hist_dE_2;
		    hist_dE_abs = hist_dE_abs_2;
		    hist_Eloss = hist_Eloss_2;
		    fname = fname_2;
		    break;
		}
	    }
	    if (0==hist_dE){
		std::cout<<"No histogram"<<std::endl;
		continue;
	    }

		TFile * file = 0;
		std::string filename = fname;
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
		std::cout<<"Opened file: "<<filename<<std::endl;
		TTree * tree = (TTree*) file->Get("ElectronHistory");
		
		tree->SetBranchAddress("energy_initial", &En_start);
		tree->SetBranchAddress("energy_final", &En_finish);
		tree->SetBranchAddress("energy_coll", &En_collision);
		tree->SetBranchAddress("energy_average", &En_avr);
		
		tree->SetBranchAddress("time_initial", &time_start);
		tree->SetBranchAddress("time_delta", &delta_time);
		tree->SetBranchAddress("time_delta_full", &delta_time_full);
		tree->SetBranchAddress("process_type", &process);
		
		tree->SetBranchAddress("position_initial",&pos_start);
		tree->SetBranchAddress("position_final",&pos_finish);
		tree->SetBranchAddress("position_delta",&delta_x);
		
		
		unsigned long int _end_ = tree->GetEntries();
		for (unsigned long int i=0;i!=_end_;++i){
		    tree->GetEntry(i);
		    //if ((1==process)||(2==process)) {
		    if ((std::fabs(En_collision)<11.5)&&((std::fabs(En_collision)>10.8))) {
			if (0==nhist) {//1dim
				hist_dE->Fill(std::fabs(En_collision) - std::fabs(En_start));
			} else {
				hist_dE->Fill(En_collision - En_start);
			}
			hist_dE_abs->Fill(std::fabs(En_collision - En_start));
			hist_Eloss->Fill(std::fabs(En_collision)- std::fabs(En_finish));
		    }
		}
	    for (int bin = 1, bin_end = hist_dE->GetNbinsX()+1; bin!=bin_end; ++bin) {
		max_val_dE = std::max(max_val_dE, (double) hist_dE->GetBinContent(bin));
	    }
	    for (int bin = 1, bin_end = hist_dE->GetNbinsX()+1; bin!=bin_end; ++bin) {
		max_val_dE_abs = std::max(max_val_dE_abs, (double) hist_dE_abs->GetBinContent(bin));
	    }
	    for (int bin = 1, bin_end = hist_Eloss->GetNbinsX()+1; bin!=bin_end; ++bin) {
		max_val_Eloss = std::max(max_val_Eloss, (double) hist_Eloss->GetBinContent(bin));
	    }
	}
	max_val_dE*=1.1;
	max_val_dE_abs*=1.1;
	max_val_Eloss*=1.1;
	gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("");
	TCanvas *c_dE = new TCanvas ("E increase", "E increase", DEF_W, DEF_H);
	c_dE->SetGrid();
	c_dE->SetTicks();
	TLegend *legend_dE = new TLegend(0.55, 0.7, 0.9, 0.9);
	//legend->SetHeader("");
	legend_dE->SetMargin(0.2);
	TH2F* frame_dE = new TH2F( "frame_dE", "E increase", 500, -EN_MAX_, EN_MAX_, 500, 0, max_val_dE);
	frame_dE->GetXaxis()->SetTitle("dE [eV]");
	frame_dE->GetYaxis()->SetTitle("");
	frame_dE->Draw();
	
	hist_dE_0->SetLineWidth(2);
	hist_dE_0->SetLineColor(kRed);
	hist_dE_0->Draw("csame");
	hist_dE_1->SetLineWidth(2);
	hist_dE_1->SetLineColor(kBlack);
	hist_dE_1->Draw("csame");
	hist_dE_2->SetLineWidth(2);
	hist_dE_2->SetLineColor(kBlue);
	hist_dE_2->Draw("csame");
	
	legend_dE->AddEntry(hist_dE_0, (std::string("7.0 Td 1D Mean = ")+std::to_string(hist_dE_0->GetMean())).c_str(), "l");
	legend_dE->AddEntry(hist_dE_1, (std::string("7.0 Td 3D Mean = ")+std::to_string(hist_dE_1->GetMean())).c_str(), "l");
	legend_dE->AddEntry(hist_dE_2, (std::string("7.0 Td 3D fixed Mean = ")+std::to_string(hist_dE_2->GetMean())).c_str(), "l");
	
	frame_dE->Draw("sameaxis");
	legend_dE->Draw("same");
	c_dE->Update();

	TCanvas *c_dE_abs = new TCanvas ("E absolute increase", "E absolute increase", DEF_W, DEF_H);
	c_dE_abs->SetGrid();
	c_dE_abs->SetTicks();
	TLegend *legend_dE_abs = new TLegend(0.55, 0.7, 0.9, 0.9);
	//legend->SetHeader("");
	legend_dE_abs->SetMargin(0.2);
	TH2F* frame_dE_abs = new TH2F( "frame_dE_abs", "E absolute increase", 500, 0, EN_MAX_, 500, 0, max_val_dE_abs);
	frame_dE_abs->GetXaxis()->SetTitle("|dE| [eV]");
	frame_dE_abs->GetYaxis()->SetTitle("");
	frame_dE_abs->Draw();
	
	hist_dE_abs_0->SetLineWidth(2);
	hist_dE_abs_0->SetLineColor(kRed);
	hist_dE_abs_0->Draw("csame");
	hist_dE_abs_1->SetLineWidth(2);
	hist_dE_abs_1->SetLineColor(kBlack);
	hist_dE_abs_1->Draw("csame");
	hist_dE_abs_2->SetLineWidth(2);
	hist_dE_abs_2->SetLineColor(kBlue);
	hist_dE_abs_2->Draw("csame");
	
	legend_dE_abs->AddEntry(hist_dE_abs_0, (std::string("7.0 Td 1D Mean =")+std::to_string(hist_dE_abs_0->GetMean())).c_str(), "l");
	legend_dE_abs->AddEntry(hist_dE_abs_1, (std::string("7.0 Td 3D Mean = ")+std::to_string(hist_dE_abs_1->GetMean())).c_str(), "l");
	legend_dE_abs->AddEntry(hist_dE_abs_2, (std::string("7.0 Td 3D fixed Mean = ")+std::to_string(hist_dE_abs_2->GetMean())).c_str(), "l");
	
	frame_dE_abs->Draw("sameaxis");
	legend_dE_abs->Draw("same");
	c_dE_abs->Update();

	TCanvas *c_Eloss = new TCanvas ("E loss", "E loss", DEF_W, DEF_H);
	c_Eloss->SetGrid();
	c_Eloss->SetTicks();
	TLegend *legend_Eloss = new TLegend(0.55, 0.7, 0.9, 0.9);
	//legend->SetHeader("");
	legend_Eloss->SetMargin(0.2);
	TH2F* frame_Eloss = new TH2F( "frame_Eloss", "E loss", 500, 0, EN_LOSS_MAX_, 500, 0, max_val_Eloss);
	frame_Eloss->GetXaxis()->SetTitle("E loss [eV]");
	frame_Eloss->GetYaxis()->SetTitle("");
	frame_Eloss->Draw();
	
	hist_Eloss_0->SetLineWidth(2);
	hist_Eloss_0->SetLineColor(kRed);
	hist_Eloss_0->Draw("csame");
	hist_Eloss_1->SetLineWidth(2);
	hist_Eloss_1->SetLineColor(kBlack);
	hist_Eloss_1->Draw("csame");
	hist_Eloss_2->SetLineWidth(2);
	hist_Eloss_2->SetLineColor(kBlue);
	hist_Eloss_2->Draw("csame");
	
	legend_Eloss->AddEntry(hist_Eloss_0, (std::string("7.0 Td 1D Mean =")+std::to_string(hist_Eloss_0->GetMean())).c_str(), "l");
	legend_Eloss->AddEntry(hist_Eloss_1, (std::string("7.0 Td 3D Mean =")+std::to_string(hist_Eloss_1->GetMean())).c_str(), "l");
	legend_Eloss->AddEntry(hist_Eloss_2, (std::string("7.0 Td 3D fixed Mean =")+std::to_string(hist_Eloss_2->GetMean())).c_str(), "l");
	
	frame_Eloss->Draw("sameaxis");
	legend_Eloss->Draw("same");
	c_Eloss->Update();
}



