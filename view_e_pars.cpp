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

void view_e_pars(void) {
  gStyle->SetStatY(0.9);
  gStyle->SetStatX(0.9);
  TH1D* histEc = new TH1D ("EnergyC [eV]","EnergyC [eV]",300,0, 14);
  TH1D* histEi = new TH1D ("EnergyInit [eV]","EnergyInit [eV]",300,0, 14);
  TH1D* histEf = new TH1D ("EnergyFin [eV]","EnergyFin [eV]",300,0, 14);
  TH1D* hist_dE = new TH1D ("Energy increase in resonance [eV]","Energy increase in resonance [eV]",300,-0.2, 0.2);
  TH1D* hist_dE_abs = new TH1D ("Energy abs increase in resonance [eV]","Energy abs increase in resonance [eV]",300,0.0, 0.2);
  TH1D* hist_Eloss = new TH1D ("Energy loss in resonance [eV]","Energy loss in resonance [eV]",500, 0, 0.002);
  TH1D* hist_dT = new TH1D ("dTime [s]","dTime [s]",1000, 0, 2e-11);
  TH1D* hist_Tdelay = new TH1D ("Time delay [s]","Time delay [s]",600, 0, 1e-5);
  TH1D* hist_dl = new TH1D ("dL [m]","dL [m]",300, 0, 2e-6);
  TH1D* hist_V_drift = new TH1D ("Drift velocity [m/s]","Drift velocity [m/s]",300,0, 1e4);
  double Dtime_left = 1.2e-6, Dtime_right = 6.6e-6;
  TH1D* hist_T_drift = new TH1D ("Drift time [s]","Drift time [s]",160, Dtime_left, Dtime_right);
  TH1D* histEAvr = new TH1D ("Energy average", "Energy average", 300, 0, 15);
  TH1D* hist_theta = new TH1D ("Scatter angle near 10 eV","Scatter angle near 10 eV",300, 0, 3.1416);
  TH1D* hist_theta_i = new TH1D ("Initial angle distribution","Initial angle distribution",300, 0, 3.1416);
  TH1D* hist_theta_c = new TH1D ("Collision angle distribution","Collision angle distribution",300, 0, 3.1416);
  TH1D* hist_theta_f = new TH1D ("Final angle distribution","Final angle distribution",300, 0, 3.1416);
  //TH1D* hist_Ey = new TH1D ("Ey [eV]", "Ey [eV]",300, 0, 15);
  //TH2D* hist_E_Ey = new TH2D ("E Ey","E Ey",300, 0, 5, 300, 0, 5);

  int DEF_W = 900, DEF_H = 700;
  double E_at_time = 1e-11;
  double dt = 6e-12;
  double DRIFT_DISTANCE = 3e-3;
  //std::string fname1("Optimization/eData_7.0Td.root");
  std::string fname1("Output/v16.3/eData_7.0Td_VN.root");
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
  for (unsigned long int i=0;i!=_end_;++i) {
    tree->GetEntry(i);
    histEc->Fill(En_collision);
    histEi->Fill(En_start);
    histEf->Fill(En_finish);
    hist_dT->Fill(delta_time);
    double delay = delta_time_full-delta_time;
    if (delay!=0)
      hist_Tdelay->Fill(delay);
    hist_dl->Fill(delta_l);
    histEAvr->Fill(En_collision, delta_time);
    hist_theta_i->Fill(theta_start);
    hist_theta_c->Fill(theta_collision);
    hist_theta_f->Fill(theta_finish);
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
  if (NULL==tree) {
    std::cout<<"Warning: no \"ElectronProcessCounters\" TTree, using process counting from \"ElectronHistory\""<<std::endl;
  } else {
    num_of_elastic = 0, num_of_resonance=0, num_of_overflow = 0, num_of_ions = 0, num_of_excitations=0;
    Long64_t *processes_counters_ = 0;
    Short_t *processes_IDs_ = 0;
    char **processes_legends_ = 0;
    char * particle_name = new char [100];
    UInt_t particle_ID = -1;
    UInt_t processes_size_ = 0;
    _end_ = tree->GetEntries();
    for (unsigned long int i=0; i!=_end_; ++i) {
      TBranch * branch_size = tree->GetBranch("proc_size");
      TBranch * branch_particle_name = tree->GetBranch("particle_name");
      TBranch * branch_particle_ID = tree->GetBranch("particle_ID");
      branch_size->SetAddress(&processes_size_);
      branch_particle_name->SetAddress(particle_name);
      branch_particle_ID->SetAddress(&particle_ID);
      branch_size->GetEntry(i);
      branch_particle_name->GetEntry(i);
      branch_particle_ID->GetEntry(i);
      std::string name(particle_name);
      std::cout<<"Particle #"<<branch_particle_ID<<" \""<<name<<"\":"<<std::endl;
      std::cout<<"Processes number: "<<processes_size_<<std::endl;
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
      for (int proc = 0; proc!=processes_size_; ++proc) {
        std::cout<<"Process "<<processes_IDs_[proc]<<": "<<processes_counters_[proc]<<" counts"<<std::endl;
        if (name == "Argon Van der Waals Molecule") {
          switch (processes_IDs_[proc]) {
            case 0:{
              num_of_elastic += processes_counters_[proc];
              break;
            }
            case 47: {
              num_of_disso_attach = processes_counters_[proc];                         
              break;                          
            }
            case -1:{
              break;
            }
            case -2:{
              num_of_overflow += processes_counters_[proc];
              break;
            }
          }
        }
        if (name == "Argon") {
          switch (processes_IDs_[proc]) {
            case 0:{
              num_of_elastic += processes_counters_[proc];
              break;
            }
            case 1: {
              num_of_NBrS = processes_counters_[proc];
            }
            case 2: {
              num_of_ions = processes_counters_[proc];
              break;
            }
            case -1:{
            break;
            }
            case -2:{
              num_of_overflow += processes_counters_[proc];
              break;
            }
            default: {
              num_of_excitations += processes_counters_[proc]; //!!!+=
            }
          }
        }
      }
    }
  }
  TCanvas *c_i = new TCanvas ("e energy initial", "e energy initial", DEF_W, DEF_H);
  histEi->Draw();
  TCanvas *c_c = new TCanvas ("e energy before collision", "e energy before collision", DEF_W, DEF_H);
  histEc->Draw();
  TCanvas *c_f = new TCanvas ("e energy after collision", "e energy after collision", DEF_W, DEF_H);
  histEf->Draw();
  TCanvas *c_2 = new TCanvas ("Delta time", "Delta time", DEF_W, DEF_H);
  c_2->SetLogy();
  hist_dT->Draw();
  TCanvas *c_2_5 = new TCanvas ("Time delay", "Time delay", DEF_W, DEF_H);
  c_2_5->SetLogy();
  hist_Tdelay->Draw();
  gStyle->SetOptFit();
  TCanvas *c_3 = new TCanvas ("Drift velocity", "Drift velocity", DEF_W, DEF_H);
  hist_V_drift->Draw();
  TCanvas *c_3_5 = new TCanvas ("Drift time", "Drift time", DEF_W, DEF_H);
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
  c_3_5->Update();
  hist_T_drift->Fit(ff);
  TLine* line = new TLine();
  line->SetX1(ff->GetParameter(3));
  line->SetX2(ff->GetParameter(3));
  line->SetY1(c_3_5->GetUymin());
  line->SetY2(c_3_5->GetUymax());
  line->SetLineColor(kRed);
  hist_T_drift->Draw();
  line->Draw("same");
  ff->Draw("same");
  TCanvas *c_4 = new TCanvas ("Delta L", "Delta L", DEF_W, DEF_H);
  hist_dl->Draw();
  TCanvas *c_5 = new TCanvas ("e energy average", "e energy average", DEF_W, DEF_H);
  histEAvr->Draw();
  TCanvas *c_6 = new TCanvas ("Scatter angle at 10eV", "Scatter angle at 10eV", DEF_W, DEF_H);
  hist_theta->Draw();
  TCanvas *c_6i = new TCanvas ("Angle initial", "Angle initial", DEF_W, DEF_H);
  hist_theta_i->Draw();
  TCanvas *c_6c = new TCanvas ("Angle before collision", "Angle before collision", DEF_W, DEF_H);
  hist_theta_c->Draw();
  TCanvas *c_6f = new TCanvas ("Angle final", "Angle final", DEF_W, DEF_H);
  hist_theta_f->Draw();
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
  std::cout<<"Num of e dissos. attachment: "<<num_of_disso_attach<<std::endl;
  std::cout<<"Num of NBrS radiation: "<<num_of_NBrS<<std::endl;
  std::cout<<"Num of excitations: "<<num_of_excitations<<std::endl;
  std::cout<<"Num of resonance: "<<num_of_resonance<<std::endl;
  std::cout<<"Num of ionizations: "<<num_of_ions<<std::endl;
  std::cout<<"Num of elastic: "<<num_of_elastic<<std::endl;
  std::cout<<"Num of electrons: "<<num_of_electrons<<std::endl;
}

 


