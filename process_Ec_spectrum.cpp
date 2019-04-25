 {
	int DEF_W = 900, DEF_H = 700;
	double EN_MIN_=0;
	double EN_MAX_=14;
	int NN = 800;
	TH1D* histE_0 = new TH1D ("EnergyC 7.0 Td [eV]","EnergyC 7.0 Td [eV]",NN, EN_MIN_, EN_MAX_);
	TGraph* graph = new TGraph();
	histE_0->SetStats(false);
	//gauss peak
	double En_artifical_peak = 11.103;
	double sigma_artifical_peak = 0.05;
	double rel_area_artifical_peak = 0.10;	
	std::vector<double> Es, Emods, Vs, Vmods;
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
			Es.push_back(En);
			Vs.push_back(val);
	    }
		double En = 0;
		double Val = 0;
		int ind = 0;
		double dE_last = Es[Es.size()-1] - Es[Es.size()-2];
		if (Es.back()<EN_MAX_) { //ensuring that [ind+1] in the following code is always valid.
			Es.push_back(EN_MAX_);
			Vs.push_back(0);
		}
		if (Es.back()==EN_MAX_) {
			Es.push_back(EN_MAX_ + dE_last);
			Vs.push_back(0);
		}
		Val = Vs.front();
		Es.insert(Es.begin(), En);
		Vs.insert(Vs.begin(), Val);
		En = Es.front();
		while (En< EN_MAX_) {
			double extra = amplitude*exp(-0.5*pow( (En - En_artifical_peak)/sigma_artifical_peak , 2));
			extra = std::max(extra, 0.0);		
			Emods.push_back(En);	
			Vmods.push_back((Val + extra)/Norm2);
			//stepping is not trivial, because it is necessary to decrease dE near peak.			
			bool hist_step = true;			
			if ((Es[ind+1] > (En_artifical_peak - 3*sigma_artifical_peak)) 
				&& (En < (En_artifical_peak + 3*sigma_artifical_peak)))	{
				double dE_p = sigma_artifical_peak/15;
				double dE_h = Es[ind+1]-En;
				if (dE_h>dE_p)
					hist_step = false;
				En+= (hist_step ? dE_h : dE_p);
				Val = Vs[ind] + (En - Es[ind])*(Vs[ind+1] - Vs[ind])/(Es[ind+1]-Es[ind]);
				if (hist_step)
					++ind; 			
			} else {
				En = Es[ind+1];
				Val = Vs[ind+1];
				++ind;
			}
		}
		double extra = amplitude*exp(-0.5*pow( (En - En_artifical_peak)/sigma_artifical_peak , 2));
		extra = std::max(extra, 0.0);
		Emods.push_back(En);	
		Vmods.push_back((Val + extra)/Norm2);
		for (int i = 0, i_end_ = Vmods.size(); i!=i_end_; ++i)
			max_val = std::max(max_val, Vmods[i]);
		for (int i = 0, i_end_ = Vs.size(); i!=i_end_; ++i)
			max_val = std::max(max_val, Vs[i]);
		
	}

	max_val*=1.05;
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
	
	graph->Set(Vmods.size());
	for (int i = 0, i_end_ = Vmods.size(); i!=i_end_; ++i)
		graph->SetPoint(i, Emods[i], Vmods[i]);
	graph->SetLineWidth(2);
	graph->SetLineColor(kBlack);//(kYellow-3);
	graph->Draw("C");

	std::ofstream str0;
	str0.open(out_fname0, std::ios_base::trunc);
	str0<<"//\""<<fname_0+".root\" data"<<std::endl;
	str0<<"//E[eV]\tSpectrum value"<<std::endl;
	for (int i = 0, i_end_ = Vs.size(); i!=i_end_; ++i)
		if (Es[i]<=EN_MAX_)
			str0<<Es[i]<<"\t"<<Vs[i]<<std::endl;
	str0.close();
	str0.open(out_fname0_1, std::ios_base::trunc);
	str0<<"//\""<<fname_0+".root\" data"<<std::endl;
	str0<<"//E[eV]\tModified spectrum value"<<std::endl;
	for (int i = 0, i_end_ = Vmods.size(); i!=i_end_; ++i)
		str0<<Emods[i]<<"\t"<<Vmods[i]<<std::endl;
	str0.close();
	frame->Draw("sameaxis");
	c_->Update();

}



