#include "Manager.h"

Manager::Manager(ArDataTables *Ar_tables, UInt_t RandomSeed) : is_ready_(false), eField_(0), Concentration_(0), Coefficient_(0),
	skip_counter_(0), ArTables_(Ar_tables), skipping_early_events(true), num_of_events(0), num_of_10millions(0)
{
	random_generator_ = new TRandom1(RandomSeed); //TRandom3 is a default. TRandom1 is slower but better
	start_seed_ = RandomSeed;
	processes_size_ = ArTables_->ArAllData_.ArExper_.max_process_ID + 5; //5 = elastic + 2resonances + None + Overthrow
	processes_counters_ = new Long64_t [processes_size_];
	processes_IDs_ = new Short_t [processes_size_];
	processes_legends_ = new char * [processes_size_];
	for (int i = 0, i_end_ = processes_size_; i!=i_end_; ++i) {
		processes_IDs_[i] = i-2;
		processes_counters_[i] = 0;
		switch (processes_IDs_[i]) {
		case (Event::Overflow): {
			processes_legends_[i] = c_str_cp(std::string("\"Overflow (energy exceeds ")+std::to_string(EN_MAXIMUM_) + " eV)\"");
			break;
		}
		case (Event::None): {
			processes_legends_[i] = c_str_cp(std::string("\"None\""));
			break;
		}
		case (Event::Elastic): {
			processes_legends_[i] = c_str_cp(std::string("\"Elastic scattering\""));
			break;
		}
		case (Event::Resonance_3o2): {
			processes_legends_[i] = c_str_cp(std::string("\"Feshbach resonance 3/2\""));
			break;
		}
		case (Event::Resonance_1o2): {
			processes_legends_[i] = c_str_cp(std::string("\"Feshbach resonance 1/2\""));
			break;
		}
		default: {
			InelasticProcess *p = ArTables_->ArAllData_.ArExper_.FindInelastic(processes_IDs_[i]-Event::Ionization);
			if (NULL!=p) {
				processes_legends_[i] = c_str_cp(p->get_name());
			} else {
				processes_legends_[i] = c_str_cp(std::string("\"?Unknown?\""));
			}
		}
		}
	}
	InitTree();
}

Manager::~Manager()
{
	delete [] processes_counters_;
	delete [] processes_IDs_;
	for (int i = 0, i_end_ = processes_size_; i!=i_end_; ++i) {
		delete [] processes_legends_[i];
	}
	delete [] processes_legends_;
	random_generator_->Delete();
	if (sim_data_)
		sim_data_->Delete();
	if (processes_data_)
		processes_data_->Delete();
}

void Manager::InitTree (void)
{
	sim_data_ = new TTree("ElectronHistory", "ElectronHistory");
	sim_data_->Branch("process_type", &event_.process);

	sim_data_->Branch("time_initial", &event_.time_start);
	sim_data_->Branch("time_delta", &event_.delta_time);
	//sim_data_->Branch("time_delta_full", &event_.delta_time_full);

	if (0==SKIP_HISTORY_) { //for spectra
		sim_data_->Branch("energy_initial", &event_.En_start);
		sim_data_->Branch("energy_coll", &event_.En_collision);
		sim_data_->Branch("energy_final", &event_.En_finish);
		//sim_data_->Branch("energy_average", &event_.En_avr);
	} else { //for V drift
		sim_data_->Branch("position_final",&event_.pos_finish);
		//sim_data_->Branch("position_initial",&event_.pos_start);
		//sim_data_->Branch("position_delta",&event_.delta_x);
	}
	sim_data_->Branch("theta_initial", &event_.theta_start);
	sim_data_->Branch("theta_coll", &event_.theta_collision);
	sim_data_->Branch("theta_final", &event_.theta_finish);
	sim_data_->Branch("theta_delta", &event_.delta_theta);
	sim_data_->Branch("path_delta", &event_.delta_l);

	/*sim_data_->Branch("deb_log_rand",&event_.deb_log_rand);
	sim_data_->Branch("deb_solver_y_left",&event_.deb_solver_y_left);
	sim_data_->Branch("deb_solver_y_right",&event_.deb_solver_y_right);
	sim_data_->Branch("deb_solver_E_left",&event_.deb_solver_E_left);
	sim_data_->Branch("deb_solver_E_right",&event_.deb_solver_E_right);*/
	processes_data_ = new TTree("ElectronProcessCounters", "ElectronProcessCounters");
	processes_data_->Branch("proc_size", &processes_size_, "proc_size/i");
	processes_data_->Branch("proc_IDs", processes_IDs_, "proc_IDs[proc_size]/S");
	processes_data_->Branch("proc_Ns", processes_counters_, "proc_Ns[proc_size]/L");
	processes_data_->Branch("proc_names", processes_legends_, "proc_names[proc_size]/C");
}

void Manager::Clear(void)
{
	is_ready_ = false;
	sim_data_->Delete();
	skip_counter_ = 0;
	skipping_early_events = true;
	num_of_events = 0;
	num_of_10millions = 0;
	for (int i = 0, i_end_ = processes_size_; i!=i_end_; ++i) {
		processes_IDs_[i] = i-2;
		processes_counters_[i] = 0;
	}
	InitTree();
}

void Manager::SetParameters(double Concentr /*in SI*/, double E /*in SI*/)
{
	Concentration_ = std::fabs(Concentr);
	eField_ = std::fabs(E);
	is_ready_ = true;
	if ((0 == Concentration_) || (0 == eField_)) {
		std::cout << "Manager::SetParameters(): Can't have 0 values" << std::endl;
		is_ready_ = false;
	}
	Coefficient_ = 1e20/* *e_charge_SIconst*/ * eField_ / Concentration_; // me/mu is set to 1.0
	//1e20 is because XS is expressed in 1e-20 m^2, not in m^2
}

void Manager::SetParameters(double T /*in K*/, double Pressure /*in SI*/, double E /*in SI*/)
{
	if (0 == T) {
		std::cout << "Manager::SetParameters(): Can't have 0 temperature" << std::endl;
		is_ready_ = false;
		return;
	}
	SetParameters(Pressure / (T*boltzmann_SIconst), E);
}

//Int XS(e)*sqrt(e/(e-Eny))*de
//'to' is always > 'from'
long double Manager::XS_integral(long double from, long double to, long double Eny, Event &event)
{
	double E = from, E_prev = from;
	long double dx;
	long double Int = 0;
	if ((from-Eny)/from<1e-6) {//irregularity case
		E = 1e-6*from + Eny;
		Int+=0.5*(ArTables_->TotalCrossSection(E)+ArTables_->TotalCrossSection(from))*sqrt(E_prev*(E_prev-Eny));
		from = E;
		E_prev = from;
	}
	ColoredRange energy_range_ = ColoredInterval (0, 0.1, 2e-4) + ColoredInterval (0.1, 1, 8e-4) +
		ColoredInterval (1, 10, 1e-2) + ColoredInterval (10, EN_MAXIMUM_, 0.02) +
		ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/5) +	//coarse area
		ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/5) +	//coarse area
		ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) + 	//fine area
		ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80) +	//fine area
		ColoredInterval (11.5, EN_MAXIMUM_, 0.003);
	energy_range_.Trim(from, to);
	int i_end_ = energy_range_.NumOfIndices();
	event.deb_N_integral = i_end_;
	for (int i=0; i!=i_end_;++i) {
		E = energy_range_.Value(i);
		Int+=ArTables_->TotalCrossSection(E)*sqrt(E/(E-Eny))*(E-E_prev);
		E_prev = E;
	}
	return Int;
}

long double Manager::XS_integral_table(long double from, long double to, long double Eny, Event &event)
{
	return (*ArTables_->integral_table_)(to, Eny) - (*ArTables_->integral_table_)(from, Eny);
}

long double Manager::XS_integral_for_test(long double from, long double to, long double Eny, long double dE)
{
	double E = from, E_prev = from;
	long double dx;
	long double Int = 0;
	if ((from - Eny) / from<1e-6) {//irregularity case
		E = std::min(1e-6*from + Eny, to);
		Int += 0.5*(ArTables_->TotalCrossSection(E) + ArTables_->TotalCrossSection(from))*sqrt(E*(E - Eny));
		from = E;
		E_prev = from;
	}
	if (from==to)
		return Int;
	while (E<to) {
		Int += ArTables_->TotalCrossSection(E)*sqrt(E / (E - Eny))*dE;
		E_prev = E;
		E += dE;
	}
	E = to;
	Int += ArTables_->TotalCrossSection(E)*sqrt(E / (E - Eny))*(E - E_prev);
	return Int;
}

void Manager::Initialize(Event &event)
{
	if (!is_ready_)
		return;
	skip_counter_=0;
	skipping_early_events = true;
	num_of_events = 0;
	num_of_10millions = 0;
	//processes_counters; //preserve
	start_seed_ = random_generator_->GetSeed();
	event.CrossSections.resize(ArTables_->ArAllData_.ArExper_.max_process_ID + 3, 0); //3==Elastic + Resonances
	event.CrossSectionsSum.resize(ArTables_->ArAllData_.ArExper_.max_process_ID + 3, 0);
	event.En_start = 1 + 3*random_generator_->Uniform();
	//event.En_start = 8;
	event.En_collision = 0;
	event.En_finish = 0;
	event.pos_start = 0;
	event.pos_finish = 0;
	event.delta_x = 0;
	event.time_start = 0;
	event.delta_time = 0;
	event.delta_time_full = 0;
	event.process = Event::None;

	event.theta_start = 0;
	event.theta_collision = 0;
	event.theta_finish = 0;
	event.delta_theta = 0;
	event.delta_l = 0;

	event_ = event;
}

void Manager::Solve (long double LnR, Event &event) //TODO: tabulate
{
	event.deb_log_rand = LnR;
	long double Eny = event.En_start*sin(event.theta_start)*sin(event.theta_start);
	long double e_start = event.En_start;
	short case_ = 0; //0 - normal Vx_start>0, 1 - Vx_start<0, no sign change, 2 - Vx_start<0, sign is changed
	if (event.theta_start>M_PI/2) {
		double TURN_INT = XS_integral(Eny, event.En_start, Eny, event);
		if (LnR>TURN_INT) { //Vx changes its sign.
			e_start = Eny;
			LnR = LnR - TURN_INT;
			case_ = 2;
		} else { //Vx is always < 0.
			LnR = TURN_INT - LnR; //energy decrease from En_start to E_coll is changed to energy increase from Eny to E_coll
			e_start = Eny;
			case_ = 1;
		}
	}
	//approximate right point. By default it is EN_MAXIMUM_, but there are a lot of extra calculations in this case
	long double left = e_start;
	long double right = e_start + 10*LnR/(ArTables_->TotalCrossSection(left)*sqrt(left/(left-Eny)));
	right = std::min((double)right, EN_MAXIMUM_);
	long double I_max = XS_integral(e_start, right, Eny, event);
	if (I_max < LnR) {
		I_max = XS_integral(e_start, EN_MAXIMUM_, Eny, event);
		right = EN_MAXIMUM_;
	}
	long double e_finish = right;
	event.En_collision = right;
	long double prev_solution = e_start;
	long double f_left = -LnR, f_right = I_max-LnR, f_new;
	double convergence_criteria = 2e-4*std::min(e_start, e_finish - prev_solution);
	if (I_max < LnR) {
		event.process = Event::Overflow;
		event.deb_solver_y_left = 0;
		event.deb_solver_y_right = LnR-I_max;
		event.deb_solver_E_left = left;
		event.deb_solver_E_right = right;
		event.deb_solver_E_delta = 0;
	} else {
		while (convergence_criteria < std::fabs(e_finish - prev_solution)) {
			prev_solution = e_finish;
			e_finish = (left*f_right - right*f_left) / (f_right - f_left);
			f_new = XS_integral(e_start, e_finish, Eny, event) - LnR;
			if (f_new < 0) {
				left = e_finish;
				f_left = f_new;
			} else {
				right = e_finish;
				f_right = f_new;
			}
			convergence_criteria = 2e-4*(std::min(std::max(e_start, (long double)0.1*XS_EL_EN_MINIMUM_),
					std::max(e_finish, (long double)0.1*XS_EL_EN_MINIMUM_)));//std::max(2e-6, std::fabs(5e-4*event.En_finish));
		}
		event.En_collision = e_finish;
		event.deb_solver_y_left = f_left;
		event.deb_solver_y_right = f_right;
		event.deb_solver_E_left = left;
		event.deb_solver_E_right = right;
		event.deb_solver_E_delta = e_finish - prev_solution;
	}
	event.En_collision = e_finish;
	if (1==case_) {
		event.theta_collision = acos(-sqrt((event.En_collision-Eny)/event.En_collision));
	} else {
		event.theta_collision = acos(sqrt((event.En_collision-Eny)/event.En_collision));
	}
}

void Manager::Solve_table (long double LnR, Event &event)
{
	event.deb_log_rand = LnR;
	long double Eny = event.En_start*sin(event.theta_start)*sin(event.theta_start);
	long double e_start = event.En_start;
	short case_ = 0; //0 - normal Vx_start>0, 1 - Vx_start<0, no sign change, 2 - Vx_start<0, sign is changed
	double INT = XS_integral_table(Eny, event.En_start, Eny, event);
	if (event.theta_start>M_PI/2) {
		if (LnR>INT) { //Vx changes its sign.
			INT = LnR - INT; //INT(Eny, Eny)===0;
			case_ = 2;
		} else { //Vx is always < 0.
			INT = INT - LnR; //energy decrease from En_start to E_coll is changed to energy increase from Eny to E_coll
			case_ = 1;
		}
	} else {
		INT = INT + LnR;
	}
	event.En_collision = ArTables_->integral_table_->find_E(Eny, INT); //TODO: add debug info - uncertainty and overflow - status output variable
	if (isnan(event.En_collision))
		event.En_collision = -1;
	if (event.En_collision<0) {
		std::cout<<"Manager::Solve_table::Error: found E<0 in table. Eny= "<<Eny<<" Ei= "<<event.En_start<<" Int= "<<INT<<std::endl;
		event.En_collision = std::min(0.01*event.En_start, 0.001) + event.En_start;
	}
	if (1==case_) {
		event.theta_collision = acos(-sqrt((event.En_collision-Eny)/event.En_collision));
	} else {
		event.theta_collision = acos(sqrt((event.En_collision-Eny)/event.En_collision));
	}
	INT+=1;
}

//Solve with high accuracy integral
void Manager::Solve_test (long double LnR, Event &event)
{
	event.deb_log_rand = LnR;
	long double Eny = event.En_start*sin(event.theta_start)*sin(event.theta_start);
	long double e_start = event.En_start;
	short case_ = 0; //0 - normal Vx_start>0, 1 - Vx_start<0, no sign change, 2 - Vx_start<0, sign is changed
	if (event.theta_start>M_PI/2) {
		double TURN_INT = XS_integral(Eny, event.En_start, Eny, event);
		if (LnR>TURN_INT) { //Vx changes its sign.
			e_start = Eny;
			LnR = LnR - TURN_INT;
			case_ = 2;
		} else { //Vx is always < 0.
			LnR = TURN_INT - LnR; //energy decrease from En_start to E_coll is changed to energy increase from Eny to E_coll
			e_start = Eny;
			case_ = 1;
		}
	}
	//approximate right point. By default it is EN_MAXIMUM_, but there are a lot of extra calculations in this case
	long double left = e_start;
	long double right = e_start + 10*LnR/(ArTables_->TotalCrossSection(left)*sqrt(left/(left-Eny)));
	right = std::min((double)right, EN_MAXIMUM_);
	long double I_max = XS_integral_for_test(e_start, right, Eny, 1e-5);
	if (I_max < LnR) {
		I_max = XS_integral_for_test(e_start, EN_MAXIMUM_, Eny, 1e-5);
		right = EN_MAXIMUM_;
	}
	long double e_finish = right;
	event.En_collision = right;
	long double prev_solution = e_start;
	long double f_left = -LnR, f_right = I_max-LnR, f_new;
	double convergence_criteria = 2e-4*std::min(e_start, e_finish - prev_solution);
	if (I_max < LnR) {
		event.process = Event::Overflow;
		event.deb_solver_y_left = 0;
		event.deb_solver_y_right = LnR-I_max;
		event.deb_solver_E_left = left;
		event.deb_solver_E_right = right;
		event.deb_solver_E_delta = 0;
	} else {
		while (convergence_criteria < std::fabs(e_finish - prev_solution)) {
			prev_solution = e_finish;
			e_finish = (left*f_right - right*f_left) / (f_right - f_left);
			f_new = XS_integral_for_test(e_start, e_finish, Eny, 1e-5) - LnR;
			if (f_new < 0) {
				left = e_finish;
				f_left = f_new;
			} else {
				right = e_finish;
				f_right = f_new;
			}
			convergence_criteria = 2e-4*(std::min(std::max(e_start, (long double)0.1*XS_EL_EN_MINIMUM_),
					std::max(e_finish, (long double)0.1*XS_EL_EN_MINIMUM_)));//std::max(2e-6, std::fabs(5e-4*event.En_finish));
		}
		event.En_collision = e_finish;
		event.deb_solver_y_left = f_left;
		event.deb_solver_y_right = f_right;
		event.deb_solver_E_left = left;
		event.deb_solver_E_right = right;
		event.deb_solver_E_delta = e_finish - prev_solution;
	}
	event.En_collision = e_finish;
	if (1==case_) {
		event.theta_collision = acos(-sqrt((event.En_collision-Eny)/event.En_collision));
	} else {
		event.theta_collision = acos(sqrt((event.En_collision-Eny)/event.En_collision));
	}
}

//Indefinite integral of sqrt(1+x*2*cos_th+x*x)dx
//see wolfram mathematica or gradshteyn
long double Manager::Path_integral (long double x, long double cos_th)
{
	long double R = sqrt(1+x*2*cos_th+x*x);
	long double out = 0.5*(cos_th+x)*R;
	out -=0.5*(cos_th*cos_th-1)*log(2*cos_th+2*(R + x));
	return out;
}

void Manager::DoStepLength(Event &event)
{
	if (!is_ready_)
		return;

	long double L = - log(random_generator_->Uniform());
	L *= Coefficient_; //Calculated once for fixed parameters;
	//solving L = XS_integral(Ei, Ec) for Ec===E collision.
	Solve_table(L, event);

	//Energy is in eV
	long double vel_0 = sqrt(2.0*e_charge_SIconst*event.En_start / e_mass_SIconst)*cos(event.theta_start);
	long double vel_1 = sqrt(2.0*e_charge_SIconst*event.En_collision / e_mass_SIconst)*cos(event.theta_collision);
	event.delta_time = (vel_1 - vel_0)*e_mass_SIconst / (e_charge_SIconst * eField_); //in s
	event.delta_x = (event.En_collision - event.En_start) / eField_;
	event.pos_finish = event.pos_start + event.delta_x;
	event.En_avr = event.En_start + eField_*vel_0*event.delta_time / 4.0 + e_charge_SIconst*std::pow(eField_*event.delta_time, 2) / (6 * e_mass_SIconst);
	//Total path in meters:
	//calculated as Int{sqrt((dy/dt)^2 + (dx/dt)^2 ) dt} from 0 to delta_time
	long double up_limit = event.delta_time*eField_*sqrt(e_charge_SIconst/(2*e_mass_SIconst*event.En_start));
	event.delta_l = Path_integral(up_limit, event.theta_start) - Path_integral(0, event.theta_start);
	event.delta_l*=2*event.En_start/eField_;
	//!!!ENERGY CUT!!! TODO:remove
	event.En_collision = (event.En_collision<EN_CUT_) ? EN_CUT_: event.En_collision;
	//!!!ENERGY CUT!!! TODO:remove
	event_ = event;
}

void Manager::DoScattering(Event &event)
{
	if (!is_ready_)
		return;
	for (int i=0, end_ = event.CrossSections.size(); i!=end_; ++i) {
		event.CrossSections[i] = ArTables_->CrossSection(std::fabs(event.En_collision), i);//process type == index in cross section array.
		event.CrossSectionsSum[i] = std::max(event.CrossSections[i], 0.0) + ((i==0) ? 0.0 : event.CrossSectionsSum[i-1]);
	}
	for (int i=0, end_ = event.CrossSections.size(); i!=end_; ++i)
		event.CrossSectionsSum[i] /= event.CrossSectionsSum[end_-1];

	bool is_overflow = (event.process==Event::Overflow);
	double R1 = random_generator_->Uniform();
	for (int i=0, end_ = event.CrossSections.size(); i!=end_; ++i)
		if (R1<event.CrossSectionsSum[i]) {
			event.process = i;
			break;
		}

	double R2 = random_generator_->Uniform();
	event.delta_theta = ArTables_->generate_Theta (event.En_collision, event.process, R2);
	double phi = random_generator_->Uniform()*2.0*M_PI;
	double cos_th_f = std::cos(event.delta_theta)*std::cos(event.theta_collision) + std::sin(event.delta_theta)*std::sin(event.theta_collision)*std::cos(phi);
	if (cos_th_f<-1.0) //just in case of precision problems
		cos_th_f = -1.0;
	if (cos_th_f>1.0)
		cos_th_f = 1.0;
	event.theta_finish = std::acos(cos_th_f);// event.theta_collision + (R3<0.5 ? event.delta_theta: -event.delta_theta); - v7.x - 2D implementation
	long double gamma_f = e_mass_eVconst/Ar_mass_eVconst;
	double EnergyLoss = 2*(1-cos(event.delta_theta))*event.En_collision*gamma_f /pow(1 + gamma_f, 2);
	switch (event.process) {
		case (Event::Resonance_3o2): {
			EnergyLoss *=RESONANCE_EN_LOSS_FACTOR_;
#ifdef RESONANCE_EN_LOSS_
			EnergyLoss = RESONANCE_EN_LOSS_;
#endif
			break;
		}
		case (Event::Resonance_1o2): {
			EnergyLoss *=RESONANCE_EN_LOSS_FACTOR_;
#ifdef RESONANCE_EN_LOSS_
			EnergyLoss = RESONANCE_EN_LOSS_;
#endif
			break;
		}
		case (Event::Elastic): {
			break;
		}
		case (Event::Ionization): {
			InelasticProcess *p = ArTables_->ArAllData_.ArExper_.FindInelastic(event.process-Event::Ionization);
			if (NULL!=p) {
				EnergyLoss = (event.En_collision > 0) ? p->get_En_thresh() : -p->get_En_thresh();
				EnergyLoss += (event.En_collision - EnergyLoss)/2.0; //Consider that residual energy is equally divided between 2 electrons
			}
			break;
		}
		default: {
			InelasticProcess *p = ArTables_->ArAllData_.ArExper_.FindInelastic(event.process-Event::Ionization);
			if (NULL!=p)
				EnergyLoss = (event.En_collision > 0) ? p->get_En_thresh() : -p->get_En_thresh();
		}
	}
	event.En_finish = event.En_collision - EnergyLoss;
	event.delta_time_full = event.delta_time;
	if (event.process == Event::Resonance_3o2) {
		event.delta_time_full += random_generator_->Exp(h_bar_eVconst/Width_3o2_);
	}
	if (event.process == Event::Resonance_1o2) {
		event.delta_time_full += random_generator_->Exp(h_bar_eVconst/Width_1o2_);
	}
	if (is_overflow)
		event.process = Event::Overflow;
	event_ = event;
}

void Manager::PostStepAction(Event &event)
{
	if (!is_ready_)
		return;
	for (int i = 0, i_end_ = processes_size_; i!=i_end_; ++i) {
		if (processes_IDs_[i] == event.process) {
			++processes_counters_[i];
			break;
		}
	}
	if (event.pos_start>=DRIFT_DISTANCE_HISTORY) {
		skipping_early_events = false;
		skip_counter_=0;
	}
	if (skip_counter_>SKIP_HISTORY_)
		skip_counter_=0;
	if ((0==num_of_events)||IsFinished(event)) {
		sim_data_->Fill();
		++skip_counter_;
		++num_of_events;
		return;
	}
	if ((num_of_10millions+1)<=(num_of_events/10000000)) {
		std::cout<<"Processed "<<(num_of_10millions+1)<<"e7 collisions. Position: "<<event.pos_finish<<std::endl;
		++num_of_10millions;
	}
	if (!skipping_early_events) {
		if ((0==skip_counter_)||(event.process>Event::Elastic)) {
			sim_data_->Fill();
			skip_counter_=0;
		}
		++skip_counter_;
	}
	++num_of_events;
}

void Manager::DoGotoNext(Event &event)
{
	if (Event::None!=event.process) { //== for the very first event
		event.time_start += event.delta_time_full;
		event.En_start = event.En_finish;
		event.pos_start = event.pos_finish;
		event.process = Event::None;
		event.delta_time = 0;
		event.delta_x = 0;
		event.delta_time_full = 0;
		event.theta_start = event.theta_finish;
	}
	event_ = event;
}

void Manager::DoStep(Event &event)
{
	if (!is_ready_)
		return;
	DoGotoNext(event);
	DoStepLength(event);
	DoScattering(event);
	PostStepAction(event);
}

bool Manager::IsFinished(Event &event)
{
	if (!is_ready_)
		return true;
	//return sim_data_->GetEntries()>=100;
	return !(event.pos_finish < DRIFT_DISTANCE_);
	//return !((event.pos_finish < DRIFT_DISTANCE_)&&(event.pos_finish> -10*DRIFT_DISTANCE_));
}

void Manager::LoopSimulation(void)
{
	if (!is_ready_)
		return;
	Initialize(event_);
	while (!IsFinished(event_)) {
		DoStep(event_);
	}
}

void Manager::WriteHistory(std::string root_fname)
{
	TFile *file = new TFile(root_fname.c_str(), "RECREATE");
	file->cd();
	std::cout<<"Event number: "<<sim_data_->GetEntries()<<std::endl;
	sim_data_->Write("", TObject::kOverwrite);
	processes_data_->Fill();
	processes_data_->Write("", TObject::kOverwrite);
	file->Close();
}

void Manager::Test(void)
{
	std::cout << "Testing Manager->Solve():" << std::endl;
	std::cout << "Testing integral calculations: " << std::endl;
	Event dummy;
	std::cout << "\tEy=0, Ei=1eV:" << std::endl;
	std::cout << "\tEf\tXS_integral()\tdE=1e-3\tdE=1e-5\tTable" << std::endl;
	std::cout << "\t"<<2.0<<"\t" << XS_integral(1, 2, 0, dummy) << "\t" << XS_integral_for_test(1, 2, 0, 1e-3) 
		<< "\t" << XS_integral_for_test(1, 2, 0, 1e-5) << "\t" <<XS_integral_table(1,2,0,dummy)<< std::endl;
	std::cout << "\t" << 1.2 << "\t" << XS_integral(1, 1.2, 0, dummy) << "\t" << XS_integral_for_test(1, 1.2, 0, 1e-3)
		<< "\t" << XS_integral_for_test(1, 1.2, 0, 1e-5) << "\t" << XS_integral_table(1, 1.2, 0, dummy) << std::endl;
	std::cout << "\t" << 1.01 << "\t" << XS_integral(1, 1.01, 0, dummy) << "\t" << XS_integral_for_test(1, 1.01, 0, 1e-3)
		<< "\t" << XS_integral_for_test(1, 1.01, 0, 1e-5) << "\t" << XS_integral_table(1, 1.01, 0, dummy) << std::endl;

	std::cout <<std::endl << "\tEy=0, Ei=11.05eV:" << std::endl;
	std::cout << "\tEf=\tXS_integral()\tdE=1e-3\tdE=1e-5\tTable" << std::endl;
	std::cout << "\t" << 11.3 << "\t" << XS_integral(11.05, 11.3, 0, dummy) << "\t" << XS_integral_for_test(11.05, 11.3, 0, 1e-3)
		<< "\t" << XS_integral_for_test(11.05, 11.3, 0, 1e-5) << "\t" << XS_integral_table(11.05, 11.3, 0, dummy) << std::endl;
	std::cout << "\t" << 11.12 << "\t" << XS_integral(11.05, 11.12, 0, dummy) << "\t" << XS_integral_for_test(11.05, 11.12, 0, 1e-3)
		<< "\t" << XS_integral_for_test(11.05, 11.12, 0, 1e-5) << "\t" << XS_integral_table(11.05, 11.12, 0, dummy) << std::endl;
	std::cout << "\t" << 11.07 << "\t" << XS_integral(11.05, 11.07, 0, dummy) << "\t" << XS_integral_for_test(11.05, 11.07, 0, 1e-3)
		<< "\t" << XS_integral_for_test(11.05, 11.07, 0, 1e-5) << "\t" << XS_integral_table(11.05, 11.07, 0, dummy) << std::endl;

	std::cout << std::endl << "\tEy=2.5, Ei=5.0eV:" << std::endl;
	std::cout << "\tEf=\tXS_integral()\tdE=1e-3\tdE=1e-5\tTable" << std::endl;
	std::cout << "\t" << 6.0 << "\t" << XS_integral(5.0, 6.0, 2.5, dummy) << "\t" << XS_integral_for_test(5.0, 6.0, 2.5, 1e-3)
		<< "\t" << XS_integral_for_test(5.0, 6.0, 2.5, 1e-5) << "\t" << XS_integral_table(5.0, 6.0, 2.5, dummy) << std::endl;
	std::cout << "\t" << 5.1 << "\t" << XS_integral(5.0, 5.1, 2.5, dummy) << "\t" << XS_integral_for_test(5.0, 5.1, 2.5, 1e-3)
		<< "\t" << XS_integral_for_test(5.0, 5.1, 2.5, 1e-5) << "\t" << XS_integral_table(5.0, 5.1, 2.5, dummy) << std::endl;
	
	std::cout << std::endl << "\tEy=5.0, Ei=5.0eV:"<<std::endl;
	std::cout << "\tEf=\tXS_integral()\tdE=1e-3\tdE=1e-5\tTable" << std::endl;
	std::cout << "\t" << 6.0 << "\t" << XS_integral(5.0, 6.0, 5.0, dummy) << "\t" << XS_integral_for_test(5.0, 6.0, 5.0, 1e-3)
		<< "\t" << XS_integral_for_test(5.0, 6.0, 5.0, 1e-5) << "\t" << XS_integral_table(5.0, 6.0, 5.0, dummy) << std::endl;
	std::cout << "\t" << 5.1 << "\t" << XS_integral(5.0, 5.1, 5.0, dummy) << "\t" << XS_integral_for_test(5.0, 5.1, 5.0, 1e-3)
		<< "\t" << XS_integral_for_test(5.0, 5.1, 5.0, 1e-5) << "\t" << XS_integral_table(5.0, 5.1, 5.0, dummy) << std::endl;
	std::cout << "\t" << 5.01 << "\t" << XS_integral(5.0, 5.01, 5.0, dummy) << "\t" << XS_integral_for_test(5.0, 5.01, 5.0, 1e-3)
		<< "\t" << XS_integral_for_test(5.0, 5.01, 5.0, 1e-5) << "\t" << XS_integral_table(5.0, 5.01, 5.0, dummy) << std::endl;
	std::cout << std::endl;

	std::cout << "Testing Solvers:" << std::endl;
	std::cout << "\tEi = 5eV, theta = Pi/2, LnR = 0.3" << std::endl;
	std::cout << "Ec=\tOld\tdE=1e-5\tTable" << std::endl;
	double Ec0, Ec1, Ec2;
	dummy.En_start = 5;
	dummy.theta_start = M_PI / 2;
	Solve(0.3, dummy);
	Ec0 = dummy.En_collision;
	Solve_test(0.3, dummy);
	Ec1 = dummy.En_collision;
	Solve_table(0.3, dummy);
	Ec2 = dummy.En_collision;
	std::cout << "\t"<<Ec0<<"\t"<<Ec1<<"\t"<<Ec2<< std::endl<<std::endl;

	std::cout << "\tEi = 1eV, theta = Pi/4, LnR = 0.01" << std::endl;
	std::cout << "Ec=\tOld\tdE=1e-5\tTable" << std::endl;
	dummy.En_start = 1;
	dummy.theta_start = M_PI / 4;
	Solve(0.01, dummy);
	Ec0 = dummy.En_collision;
	Solve_test(0.01, dummy);
	Ec1 = dummy.En_collision;
	Solve_table(0.01, dummy);
	Ec2 = dummy.En_collision;
	std::cout << "\t" << Ec0 << "\t" << Ec1 << "\t" << Ec2 << std::endl << std::endl;

	std::cout << "\tEi = 5eV, theta = Pi, LnR = 0.1" << std::endl;
	std::cout << "Ec=\tOld\tdE=1e-5\tTable" << std::endl;
	dummy.En_start = 5;
	dummy.theta_start = M_PI;
	Solve(0.1, dummy);
	Ec0 = dummy.En_collision;
	Solve_test(0.1, dummy);
	Ec1 = dummy.En_collision;
	Solve_table(0.1, dummy);
	Ec2 = dummy.En_collision;
	std::cout << "\t" << Ec0 << "\t" << Ec1 << "\t" << Ec2 << std::endl << std::endl;

	std::cout << "\tEi = 5eV, theta = 0, LnR = 0.1" << std::endl;
	std::cout << "Ec=\tOld\tdE=1e-5\tTable" << std::endl;
	dummy.En_start = 5;
	dummy.theta_start = 0;
	Solve(0.1, dummy);
	Ec0 = dummy.En_collision;
	Solve_test(0.1, dummy);
	Ec1 = dummy.En_collision;
	Solve_table(0.1, dummy);
	Ec2 = dummy.En_collision;
	std::cout << "\t" << Ec0 << "\t" << Ec1 << "\t" << Ec2 << std::endl << std::endl;

	std::cout << "===============================================" << std::endl << std::endl << std::endl;
}
