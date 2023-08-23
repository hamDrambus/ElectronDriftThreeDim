std::string dbl_to_str (double val, int precision)
{
	std::stringstream ss;
	ss<<std::fixed<<std::setprecision(precision)<<val;
	return ss.str();
}

struct table_row {
  long int num_of_ions;
  long int num_of_exc_3p5_4s; //VUV
  long int num_of_exc_3p5_4p; //VUV and NIR (1 excitation is 2 photons)

  double vuv_yield, vuv_yield_error; //Y/N (Number of photons per 1 drift electron per 1 cm drift, devided by concentration)
  double nir_yield, nir_yield_error;
};

void view_e_procs(void) {
  int Num_of_electrons = 120;
  double DRIFT_DISTANCE = 3e-1; //in cm
  double Temp = 163; //Kelvins
  double k_boltzmann = 1.38e-023; //SI
  double pressure = 0.6*101.5e3; //Pa
  double N = 1e-6*pressure/(k_boltzmann*Temp); //in cm-3;
  std::vector<double> Tds = {4.0, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 8.0, 8.5, 9.0, 10.0};
  std::vector<std::string> fnames;
  std::vector<table_row> results;
  std::string output = "Output/v18.3_VUV_and_NIR_yields/yields.txt";
  for (std::size_t i = 0, i_end_ = Tds.size(); i!=i_end_; ++i) {
    fnames.push_back(std::string("Output/v18.3_VUV_and_NIR_yields/eData_")+dbl_to_str(Tds[i], 1)+"Td_PN.root");
    table_row temp;
    temp.num_of_ions = 0;
    temp.num_of_exc_3p5_4s = 0;
    temp.num_of_exc_3p5_4p = 0;
    temp.vuv_yield = temp.vuv_yield_error = 0.0;
    temp.nir_yield = temp.nir_yield_error = 0.0;
    results.push_back(temp);
  }

  for (std::size_t f = 0, f_end_ = fnames.size(); f!=f_end_; ++f) {
    TFile * file = new TFile (fnames[f].c_str());
    if (0==file) {
      std::cout<<"File not opened"<<std::endl;
      continue;
    }
    if (!file->IsOpen()) {
      std::cout<<"File not opened"<<std::endl;
      file->Delete();
      continue;
    }
    std::cout<<"Opened file: "<<fnames[f]<<std::endl;

    TTree * tree = (TTree*) file->Get("ElectronProcessCounters");
    if (NULL==tree) {
      std::cout<<"Error: no \"ElectronProcessCounters\" TTree"<<std::endl;
      file->Delete();
      continue;
    }
    Long64_t *processes_counters_ = 0;
    Short_t *processes_IDs_ = 0;
    //char **processes_legends_ = 0; //Does not work
    char * particle_name = new char [100];
    UInt_t particle_ID = -1;
    UInt_t processes_size_ = 0;
    unsigned long int _end_ = tree->GetEntries(); //Entry is info for each particle
    for (unsigned long int i=0; i!=_end_; ++i) { //per particle
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
      std::cout<<"Particle #"<<particle_ID<<" \""<<name<<"\":"<<std::endl;
      std::cout<<"Processes number: "<<processes_size_<<std::endl;
      if (name != "Argon") //Not an argon particle
          continue;
      if (0 != processes_counters_) {
          delete [] processes_counters_;
          delete [] processes_IDs_;
      }
      processes_counters_ = new Long64_t [processes_size_];
      processes_IDs_ = new Short_t [processes_size_];
      TBranch * branch_ID = tree->GetBranch("proc_IDs");
      TBranch * branch_N = tree->GetBranch("proc_Ns");
      branch_ID->SetAddress(processes_IDs_);
      branch_N->SetAddress(processes_counters_);
      branch_ID->GetEntry(i);
      branch_N->GetEntry(i);
      for (int proc = 0; proc!=processes_size_; ++proc) {
        if (processes_IDs_[proc] == 2) //Ionization
          results[f].num_of_ions += processes_counters_[proc];
        if (processes_IDs_[proc] >= 3 && processes_IDs_[proc] <= 6) //3p5 4s excitations
          results[f].num_of_exc_3p5_4s += processes_counters_[proc];
        if (processes_IDs_[proc] >= 7 && processes_IDs_[proc] <= 16) //3p5 4p excitations
          results[f].num_of_exc_3p5_4p += processes_counters_[proc];
      }
    }
    results[f].nir_yield = results[f].num_of_exc_3p5_4p;
    results[f].nir_yield_error = std::sqrt(results[f].num_of_exc_3p5_4p); //poisson statistics
    results[f].nir_yield *= 1.0/(Num_of_electrons * DRIFT_DISTANCE * N);
    results[f].nir_yield_error *= 1.0/(Num_of_electrons * DRIFT_DISTANCE * N);

    results[f].vuv_yield = results[f].num_of_exc_3p5_4p + results[f].num_of_exc_3p5_4s;
    results[f].vuv_yield_error = std::sqrt(results[f].vuv_yield); //poisson statistics
    results[f].vuv_yield *= 1.0/(Num_of_electrons * DRIFT_DISTANCE * N);
    results[f].vuv_yield_error *= 1.0/(Num_of_electrons * DRIFT_DISTANCE * N);
  }
  std::ofstream str;
  str.open(output, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cout<<"Could not open output file \""<<output<<"\""<<std::endl;
    return;
  }
  str<<"//Td\tNion\tN4s\tN4p\tVUV_Y/N\tVUV_Y/N_err\tNIR_Y/N\tNIR_Y/N_err"<<std::endl;
  for (std::size_t f = 0, f_end_ = fnames.size(); f!=f_end_; ++f) {
    str<<dbl_to_str(Tds[f], 1)<<"\t"<<results[f].num_of_ions<<"\t"<<results[f].num_of_exc_3p5_4s<<"\t"
        <<results[f].num_of_exc_3p5_4p<<"\t"<<results[f].vuv_yield<<"\t"<<results[f].vuv_yield_error<<"\t"
        <<results[f].nir_yield<<"\t"<<results[f].nir_yield_error<<std::endl;
  }
}
