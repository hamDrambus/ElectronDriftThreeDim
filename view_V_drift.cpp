{
	std::vector<std::string> Tds = {"0.3", "1.0", "2.0", "3.0", "4.0", "4.5", "5.0", "5.5", "6.0", "6.5", "7.0", "7.5", "8.3"};
	std::string prefix = "Output/v10.2/";
	std::string Vd_data_fname = prefix + "V_drift.txt";
	std::ofstream str;
	str.open(Vd_data_fname, std::ios_base::trunc);
	if (!str.is_open())
		std::cout<<"File \""<<Vd_data_fname<<"\" not opened"<<std::endl;
	str<<"Td\tVd [m/s]\tErr"<<std::endl;
	int DEF_W = 900, DEF_H = 700;
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.9);
	TCanvas *c_0 = new TCanvas ("Drift velocity", "Drift velocity", DEF_W, DEF_H);
    
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
    
	for (std::size_t i=0, i_end_=Tds.size(); i!=i_end_; ++i) {
		std::string fname = prefix +"eData_" + Tds[i] + "Td.root";
		TFile * file = new TFile (fname.c_str());
		if (0==file) {
		    std::cout<<"File \""<<fname<<"\" not opened"<<std::endl;
		} else {
		    if (!file->IsOpen()){
				std::cout<<"File \""<<fname<<"\" not opened"<<std::endl;
		    } else {
		        std::cout<<"Opened file: "<<fname<<std::endl;
				double Td_val = std::stod(Tds[i]);
				std::string title = std::string("Drift velocity ") + Tds[i] +"Td [m/s]";
				TH1D* hist_V_drift = new TH1D (title.c_str(), title.c_str(), 300, 0, 4.5e3+ (1e4-4.5e3)*Td_val/7.0);
				long int num_of_elastic = 0, num_of_overflow = 0, num_of_electrons = 0;
				long int num_of_ions = 0, num_of_excitations=0;		        
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
		        for (unsigned long int i=0; i!=_end_; ++i) {
		            tree->GetEntry(i);
		            switch (process) {
		              case 0: {
		              ++num_of_elastic;
		              break;
		              }
		              case 1: {
		              ++num_of_ions;
		              break;
		              }
		              case -1: {
		              break;
		              }
		              case -2: {
		              ++num_of_overflow;
		              break;
		              }
		              default: {
		              ++num_of_excitations;
		              }
		            }
		            if (((i==(_end_-1))||time_start==0.0)&&i!=0) {
		                ++num_of_electrons;
		                tree->GetEntry(i-1);
		                hist_V_drift->Fill(pos_finish/(time_start+delta_time_full));
		            }
		        }
		        tree = (TTree*) file->Get("ElectronProcessCounters");
		        if (NULL!=tree) {
		            Long64_t *processes_counters_ = 0;
		            Short_t *processes_IDs_ = 0;
		            char **processes_legends_ = 0;
		            UInt_t processes_size_;
		            _end_ = tree->GetEntries();
		            if (_end_>1) {
		                std::cout<<"Warning: more that 1 processs counters in TTree \"ElectronProcessCounters\""<<std::endl;    
		            }
		            for (unsigned long int i=0; i!=_end_; ++i) {
		                TBranch * branch_size = tree->GetBranch("proc_size");
		                branch_size->SetAddress(&processes_size_);
		                branch_size->GetEntry(i);
		                if (0 != processes_counters_) {
		                    delete [] processes_counters_;
		                    delete [] processes_IDs_;
		                    delete [] processes_legends_;
		                }
		                processes_counters_ = new Long64_t [processes_size_];
		                processes_IDs_ = new Short_t [processes_size_];
		                processes_legends_ = new char * [processes_size_];
		                TBranch * branch_ID = tree->GetBranch("proc_IDs");
		                TBranch * branch_N = tree->GetBranch("proc_Ns");
		                TBranch * branch_name = tree->GetBranch("proc_names");
		                branch_ID->SetAddress(processes_IDs_);
		                branch_N->SetAddress(processes_counters_);
		                branch_name->SetAddress(processes_legends_);
		                branch_ID->GetEntry(i);
		                branch_N->GetEntry(i);
		                branch_name->GetEntry(i);
		                num_of_elastic = 0, num_of_overflow = 0, num_of_ions = 0, num_of_excitations=0;
		                std::cout<<"ProcSize: "<<processes_size_<<std::endl;
		                for (int proc = 0; proc!=processes_size_; ++proc) {
				            switch (processes_IDs_[proc]) {
				              case 0:{
				              num_of_elastic=processes_counters_[proc];
				              break;
				              }
				              case 1: {
				              num_of_ions=processes_counters_[proc];
				              break;
				              }
				              case -1:{
				              break;
				              }
				              case -2:{
				              num_of_overflow=processes_counters_[proc];
				              break;
				              }
				              default: {
				              num_of_excitations+=processes_counters_[proc]; //!!!+=
				              }
				            }
		                }
		            }
		        } else {
		        std::cout<<"Warning: no \"ElectronProcessCounters\" TTree, using process counting from \"ElectronHistory\""<<std::endl;
		        }
				std::cout<<"Num of overflow: "<<num_of_overflow<<std::endl;
				std::cout<<"Num of excitations: "<<num_of_excitations<<std::endl;
				std::cout<<"Num of ionizations: "<<num_of_ions<<std::endl;
				std::cout<<"Num of elastic: "<<num_of_elastic<<std::endl;
				std::cout<<"Num of electrons: "<<num_of_electrons<<std::endl;
				c_0->Clear();
				hist_V_drift->Draw();
				title = prefix + Tds[i]+"Td_Vdrift.png";
				c_0->Update();
				c_0->SaveAs(title.c_str(), "png");
				str<<Td_val<<"\t"<<hist_V_drift->GetMean()<<"\t"<<hist_V_drift->GetStdDev()<<std::endl;
				hist_V_drift->Delete();
			}
		}
		file->Delete();
	}
	str.close();
}


