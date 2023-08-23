void read_spectrum(TH1D* hist, std::string fname, bool average, bool renormalized_f) {
	double En_collision;
	double delta_time;
	TFile * file = 0;
	std::string filename = fname;
	filename= filename+".root";
	file = new TFile (filename.c_str());
	if (0==file) {
	    std::cout<<"File \""<<filename<<"\" not opened"<<std::endl;
	    return;
	}
	if (!file->IsOpen()){
	    std::cout<<"File \""<<filename<<"\" not opened"<<std::endl;
	    return;
	}
	std::cout<<"Opened "<<filename<<std::endl;
	TTree * tree = (TTree*) file->Get("ElectronHistory");

	tree->SetBranchAddress("energy_coll", &En_collision);
	tree->SetBranchAddress("time_delta", &delta_time);

	unsigned long int _end_ = tree->GetEntries();
	for (unsigned long int i=0;i!=_end_;++i){
	    tree->GetEntry(i);
		if (average)
		    hist->Fill(std::fabs(En_collision), delta_time*(renormalized_f ? 1.0/std::sqrt(En_collision): 1.0));
	    else
			hist->Fill(std::fabs(En_collision), renormalized_f ? 1.0/std::sqrt(En_collision): 1.0);
	}
}

double integrate(TH1D *hist, double from, double to, double renormalized_f) {
	double Int = 0;
    for (int bin = 1, bin_end = hist->GetNbinsX()+1; bin!=bin_end; ++bin) {
	    if (hist->GetBinCenter(bin)>from && hist->GetBinCenter(bin)<to) {
           Int += hist->GetBinContent(bin) * hist->GetBinWidth(bin) * (renormalized_f ? std::sqrt(hist->GetBinCenter(bin)): 1.0);
	    }
    }
	return Int;
}

void to_ascii(std::vector<TH1D*> hists, std::vector<std::string> Tds, std::string filename) {
	std::ofstream str;
	str.open(filename, std::ios_base::trunc);
	if (!str.is_open()) {
		std::cout<<"File \""<<filename<<"\" not opened"<<std::endl;
	    return;
	}
	if (hists.empty())
		return;
	str<<"//E[eV]";
	for (int hh = 0, hh_end_ = Tds.size(); hh!=hh_end_; ++hh)
		str<<"\t"<<Tds[hh]<<" Td";
	str<<std::endl;
	for (int bin = 1, bin_end = hists[0]->GetNbinsX()+1; bin!=bin_end; ++bin) {
		str<<hists[0]->GetBinCenter(bin);
		for (int hh = 0, hh_end_ = Tds.size(); hh!=hh_end_; ++hh) {
			str<<"\t"<<hists[hh]->GetBinContent(bin);
		}
		str<<std::endl;
	}
}

int export_e_distributions (void) {
	int DEF_W = 900, DEF_H = 700;
	double EN_MIN_=0;
	double EN_MAX_=15;
	int NN = 150; //dE= 0.01
	bool linear = true;
	bool renormalized_f = false;
	std::vector<Color_t> palette_major = {kBlack, kRed, kBlue, kGreen, kYellow + 2, kMagenta, kOrange + 7};
	std::string framename_Ec = std::string("e- E before collision distributions in gaseous Ar");
	std::string framename_Eavg = std::string("e- E average distributions in gaseous Ar");
	std::vector<std::string> Tds = {"5.1"};
	std::string folder = "Output/v18.1_e_distr_normal/";
	std::vector<std::string> files = {"eData_5.1Td_SN"};
	std::string output_Ec = folder + "Table_Ec_5.1Td" + (renormalized_f ? "_maxwell" : "") + ".dat";
	std::string output_Eavg = folder + "Table_Eavg_5.1Td" + (renormalized_f ? "_maxwell" : "") + ".dat";
	//std::string output_Ec = folder + "Table_Ec" + (renormalized_f ? "_maxwell" : "") + ".dat";
	//std::string output_Eavg = folder + "Table_Eavg" + (renormalized_f ? "_maxwell" : "") + ".dat";
	std::vector<TH1D*> hists_Ec, hists_Eavg;
	double max_val_Ec = 0, max_val_Eavg = 0;

	for (int hh = 0, hh_end_ = files.size(); hh!=hh_end_; ++hh) {
		std::string filename = folder + files[hh];
		std::string histname1 = "Collision energy for " + Tds[hh] + " [eV]";
		std::string histname2 = "Average energy for " + Tds[hh] + " [eV]";
		hists_Ec.push_back(new TH1D (histname1.c_str(), histname1.c_str(), NN, EN_MIN_, EN_MAX_));
		hists_Eavg.push_back(new TH1D (histname2.c_str(), histname2.c_str(), NN, EN_MIN_, EN_MAX_));
		read_spectrum(hists_Ec[hh], filename, false, renormalized_f);
		read_spectrum(hists_Eavg[hh], filename, true, renormalized_f);
		//normalize:
		double integral = integrate(hists_Ec[hh], EN_MIN_, EN_MAX_, renormalized_f);
		hists_Ec[hh]->Scale(1.0/integral);
		integral = integrate(hists_Eavg[hh], EN_MIN_, EN_MAX_, renormalized_f);
		hists_Eavg[hh]->Scale(1.0/integral);
		max_val_Ec = std::max(max_val_Ec, hists_Ec[hh]->GetBinContent(hists_Ec[hh]->GetMaximumBin()));
		max_val_Eavg = std::max(max_val_Eavg, hists_Eavg[hh]->GetBinContent(hists_Eavg[hh]->GetMaximumBin()));
	}
	max_val_Ec*= linear ? 1.2 : 2;
	//max_val_Ec = std::min(max_val_Ec, linear ? 1.0 : 10);
	max_val_Eavg*= linear ? 1.2 : 2;
	//max_val_Eavg = std::min(max_val_Eavg, linear ? 1.0 : 10);
	gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("");
	TCanvas *c_1 = new TCanvas ((std::string(" ") + framename_Ec).c_str(), (std::string(" ") + framename_Ec).c_str(), DEF_W, DEF_H);
	c_1->SetGrid(); c_1->SetTicks(); c_1->ToggleEventStatus(); c_1->ToggleToolBar();
	if (!linear)
		c_1->SetLogy();
	TLegend *legend1 = new TLegend(0.55, 0.65, 0.9, 0.9);
	//legend1->SetHeader("");
	legend1->SetMargin(0.25);
	TH2F* frame1 = new TH2F("frame1", framename_Ec.c_str(), 500, EN_MIN_, EN_MAX_, 500, linear ? 0 : 1e-5, max_val_Ec);
	frame1->GetXaxis()->SetTitle("Energy [eV]");
	if (renormalized_f)
		frame1->GetYaxis()->SetTitle("f [eV^{-3/2}]");
	else
		frame1->GetYaxis()->SetTitle("f [eV^{-1}]");
	frame1->Draw();

	for (int hh = 0, hh_end_ = hists_Ec.size(); hh!=hh_end_; ++hh) {
		hists_Ec[hh]->SetLineWidth(2);
		hists_Ec[hh]->SetLineColor(palette_major[hh]);
		hists_Ec[hh]->Draw("hist Lsame");
    }

	for (int hh = 0, hh_end_ = hists_Ec.size(); hh!=hh_end_; ++hh)
		legend1->AddEntry(hists_Ec[hh], (std::string("E/N = ") + Tds[hh]).c_str(), "l");
	frame1->Draw("sameaxis");
	legend1->Draw("same");
	c_1->Update();

	TCanvas *c_2 = new TCanvas ((std::string(" ") + framename_Eavg).c_str(), (std::string(" ") + framename_Eavg).c_str(), DEF_W, DEF_H);
	c_2->SetGrid(); c_2->SetTicks(); c_2->ToggleEventStatus(); c_2->ToggleToolBar();
	if (!linear)
		c_2->SetLogy();
	TLegend *legend2 = new TLegend(0.55, 0.65, 0.9, 0.9);
	//legend2->SetHeader("");
	legend2->SetMargin(0.25);
	TH2F* frame2 = new TH2F("frame1", framename_Eavg.c_str(), 500, EN_MIN_, EN_MAX_, 500, linear ? 0 : 1e-5, max_val_Eavg);
	frame2->GetXaxis()->SetTitle("Energy [eV]");
	if (renormalized_f)
		frame2->GetYaxis()->SetTitle("f [eV^{-3/2}]");
	else
		frame2->GetYaxis()->SetTitle("f [eV^{-1}]");
	frame2->Draw();

	for (int hh = 0, hh_end_ = hists_Eavg.size(); hh!=hh_end_; ++hh) {
		hists_Eavg[hh]->SetLineWidth(2);
		hists_Eavg[hh]->SetLineColor(palette_major[hh]);
		hists_Eavg[hh]->Draw("hist Lsame");
    }

	for (int hh = 0, hh_end_ = hists_Eavg.size(); hh!=hh_end_; ++hh)
		legend2->AddEntry(hists_Eavg[hh], (std::string("E/N = ") + Tds[hh] + " Td").c_str(), "l");
	frame2->Draw("sameaxis");
	legend2->Draw("same");
	c_2->Update();

	to_ascii(hists_Ec, Tds, output_Ec);
	to_ascii(hists_Eavg, Tds, output_Eavg);
	return 0;
}
