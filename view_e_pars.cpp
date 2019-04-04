{
    gStyle->SetStatY(0.9);
    gStyle->SetStatX(0.9);
    TH1D* histE = new TH1D ("EnergyC [eV]","EnergyC [eV]",300,0, 15);
    TH1D* hist_dE = new TH1D ("Energy increase in resonance [eV]","Energy increase in resonance [eV]",300,-0.2, 0.2);
    TH1D* hist_dE_abs = new TH1D ("Energy abs increase in resonance [eV]","Energy abs increase in resonance [eV]",300,0.0, 0.2);
    TH1D* hist_Eloss = new TH1D ("Energy loss in resonance [eV]","Energy loss in resonance [eV]",500,0, 0.003);
    TH1D* hist_dT = new TH1D ("dTime [s]","dTime [s]",1000, 0, 2e-11);
	TH1D* hist_Tdelay = new TH1D ("Time delay [s]","Time delay [s]",600, 0, 1e-12);
    TH1D* hist_dl = new TH1D ("dL [m]","dL [m]",300, 0, 2e-6);
    TH1D* hist_V_drift = new TH1D ("Drift velocity [m/s]","Drift velocity [m/s]",300,0, 1e4);
    TH1D* hist_T_drift = new TH1D ("Drift time [s]","Drift time [s]",300,4e-7, 7e-7);
    TH1D* histEAvr = new TH1D ("Energy average", "Energy average", 300, 0, 15);
    TH1D* hist_theta = new TH1D ("Scatter angle near 10 eV","Scatter angle near 10 eV",300, 0, 3.1416);
    //TH1D* hist_Ey = new TH1D ("Ey [eV]", "Ey [eV]",300, 0, 15);
    //TH2D* hist_E_Ey = new TH2D ("E Ey","E Ey",300, 0, 5, 300, 0, 5);
    
    int DEF_W = 900, DEF_H = 700;
    double E_at_time = 1e-11;
    double dt = 6e-12;
    double DRIFT_DISTANCE = 3e-3;
    std::string fname1("Output/v12.6/eData_7.0Td.root");
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
    long int num_of_ions = 0, num_of_excitations=0;
    
    TFile * file = new TFile (fname1.c_str());
    if (0==file) {
        std::cout<<"File not opened"<<std::endl;
    } else {
        if (!file->IsOpen()){
        std::cout<<"File not opened"<<std::endl;
        } else {
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
            for (unsigned long int i=0;i!=_end_;++i){
                tree->GetEntry(i);
                histE->Fill(En_collision);
                hist_dT->Fill(delta_time);
				double delay = delta_time_full-delta_time;
				if (delay>0)
					hist_Tdelay->Fill(delay);
                hist_dl->Fill(delta_l);
                histEAvr->Fill(En_collision, delta_time);
                //histEy->Fill(En_start*std::sin(theta_start)*std::sin(theta_start));
                //hist_E_Ey->Fill(En_start, En_start*std::sin(theta_start)*std::sin(theta_start));
                if ((std::fabs(En_collision)<11.5)&&((std::fabs(En_collision)>10.8))) {
                //if ((1==process)||(2==process)) {
                    hist_dE->Fill(std::fabs(En_collision) - std::fabs(En_start));
                    hist_dE_abs->Fill(std::fabs(En_collision - En_start));
                    hist_Eloss->Fill(std::fabs(En_collision)- std::fabs(En_finish));
                }
                if ((En_collision<10.05)&&(En_collision>9.95)) {
                    hist_theta->Fill(delta_theta);
                }
                switch (process) {
                  case 0:{
                  ++num_of_elastic;
                  break;
                  }
                  case 2: {
                  ++num_of_ions;
                  break;
                  }
                  case -1:{
                  break;
                  }
                  case -2:{
                  ++num_of_overflow;
                  break;
                  }
                  default: {
                  ++num_of_excitations;
                  }
                }
                //if (!((pos_finish<DRIFT_DISTANCE)&&(pos_finish>-10*DRIFT_DISTANCE))) {
                if (((i==(_end_-1))||time_start==0.0)&&i!=0) {
                    ++num_of_electrons;
                    tree->GetEntry(i-1);
                    hist_V_drift->Fill(pos_finish/(time_start+delta_time_full));
                    hist_T_drift->Fill(time_start+delta_time_full);
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
                for (unsigned long int i=0;i!=_end_;++i) {
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
                    num_of_elastic = 0, num_of_resonance=0, num_of_overflow = 0, num_of_ions = 0, num_of_excitations=0;
                    std::cout<<"ProcSize: "<<processes_size_<<std::endl;
                    for (int proc = 0; proc!=processes_size_; ++proc) {
                    switch (processes_IDs_[proc]) {
                      case 0:{
                      num_of_elastic=processes_counters_[proc];
                      break;
                      }
                      case 1:
                      case 2: {
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

            TCanvas *c_ = new TCanvas ("e energy before collision", "e energy before collision", DEF_W, DEF_H);
            histE->Draw();
            TCanvas *c_2 = new TCanvas ("Delta time", "Delta time", DEF_W, DEF_H);
			c_2->SetLogy();
            hist_dT->Draw();
			TCanvas *c_2_5 = new TCanvas ("Time delay", "Time delay", DEF_W, DEF_H);
			c_2_5->SetLogy();
            hist_Tdelay->Draw();
            TCanvas *c_3 = new TCanvas ("Drift velocity", "Drift velocity", DEF_W, DEF_H);
            hist_V_drift->Draw();
            TCanvas *c_3_5 = new TCanvas ("Drift time", "Drift time", DEF_W, DEF_H);
            hist_T_drift->Draw();
            TCanvas *c_4 = new TCanvas ("Delta L", "Delta L", DEF_W, DEF_H);
            hist_dl->Draw();
            TCanvas *c_5 = new TCanvas ("e energy average", "e energy average", DEF_W, DEF_H);
            histEAvr->Draw();
            TCanvas *c_6 = new TCanvas ("Scatter angle at 10eV", "Scatter angle at 10eV", DEF_W, DEF_H);
            hist_theta->Draw();
            TCanvas *c_7 = new TCanvas ("e delta E in resonance [eV]", "e delta E in resonance [eV]", DEF_W, DEF_H);
            hist_dE->Draw();
            TCanvas *c_8 = new TCanvas ("e energy loss in resonance", "e energy loss in resonance", DEF_W, DEF_H);
            hist_Eloss->Draw();
            TCanvas *c_9 = new TCanvas ("e delta E abs in resonance [eV]", "e delta E abs in resonance [eV]", DEF_W, DEF_H);
            hist_dE_abs->Draw();
            //TCanvas *c_14 = new TCanvas ("E Ey [eV]", "E Ey[eV]", DEF_W, DEF_H);
            //hist_E_Ey->Draw();
            //TCanvas *c_15 = new TCanvas ("Ey [eV]", "Ey[eV]", DEF_W, DEF_H);
            //hist_Ey->Draw();         
            std::cout<<"Drift distance: "<<DRIFT_DISTANCE<<std::endl;
            std::cout<<"Num of overflow: "<<num_of_overflow<<std::endl;
            std::cout<<"Num of excitations: "<<num_of_excitations<<std::endl;
            std::cout<<"Num of resonance: "<<num_of_resonance<<std::endl;
            std::cout<<"Num of ionizations: "<<num_of_ions<<std::endl;
            std::cout<<"Num of elastic: "<<num_of_elastic<<std::endl;
            std::cout<<"Num of electrons: "<<num_of_electrons<<std::endl;
        }
    }
}

 


