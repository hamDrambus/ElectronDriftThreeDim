#include "Manager.h"

Manager::Manager(const Mixture *Ar_tables) :
	skip_counter_(0), material_(Ar_tables), skipping_early_events(true), num_of_events(0)
#ifndef _NO_CERN_ROOT
	, sim_data_(NULL), processes_data_(NULL)
#endif //_NO_CERN_ROOT
{
	switch (gSettings.ProgConsts()->random_generator) {
	case (ProgramConstants::GeneratorClass::TRand1): {
#ifndef _NO_CERN_ROOT
		random_generator_ = new TRandom1();
		boost_random_generator = NULL;
		break;
#endif //_NO_CERN_ROOT
	}
	case (ProgramConstants::GeneratorClass::TRand2): {
#ifndef _NO_CERN_ROOT
		random_generator_ = new TRandom2();
		boost_random_generator = NULL;
		break;
#endif //_NO_CERN_ROOT
	}
	case (ProgramConstants::GeneratorClass::TRand3): {
#ifndef _NO_CERN_ROOT
		random_generator_ = new TRandom3();
		boost_random_generator = NULL;
		break;
#endif //_NO_CERN_ROOT
	}
	case (ProgramConstants::GeneratorClass::BOOST_hellekalek1995): {
		boost_random_generator = new boost::random::hellekalek1995();
#ifndef _NO_CERN_ROOT
		random_generator_ = NULL;
#endif //_NO_CERN_ROOT
		break;
	}
	default: {
#ifndef _NO_CERN_ROOT
		random_generator_ = NULL;
		boost_random_generator = NULL;
#endif //_NO_CERN_ROOT
	}
	}
	std::size_t N_particles = gParticleTable.GetNParticle();
	processes_counters_.resize(N_particles);
	processes_IDs_.resize(N_particles);
	processes_legends_.resize(N_particles);
	const Particle *electron = gParticleTable.GetParticle("electron");
	for (std::size_t ptcl = 0; ptcl!= N_particles; ++ptcl) {
		processes_counters_[ptcl].push_back(0);
		processes_IDs_[ptcl].push_back(Event::Overflow);
		processes_legends_[ptcl].push_back(std::string("\"Overflow (energy exceeds ")+std::to_string(gSettings.ProgConsts()->maximal_energy) + " eV)\"");
		processes_counters_[ptcl].push_back(0);
		processes_IDs_[ptcl].push_back(Event::None);
		processes_legends_[ptcl].push_back("\"None\"");
		std::size_t proc_size = gParticleTable.GetParticle(ptcl)->GetProcessSize(electron);
		for (std::size_t proc = 0; proc!= proc_size; ++proc) {
			processes_counters_[ptcl].push_back(0);
			processes_IDs_[ptcl].push_back(proc);
			processes_legends_[ptcl].push_back(gParticleTable.GetParticle(ptcl)->GetProcessName(electron, proc));
		}
	}
	InitTree();
}

Manager::~Manager()
{
	if (boost_random_generator) {
		delete boost_random_generator;
	}
#ifndef _NO_CERN_ROOT
	if (random_generator_) {
		delete random_generator_;
	}
	if (sim_data_) {
		delete sim_data_;
	}
	if (processes_data_) {
		delete processes_data_;
	}
#endif //_NO_CERN_ROOT
}

bool Manager::isReady(void) const
{
	return (boost::none != Concentration_ && boost::none != Coefficient_ && boost::none != eField_
		&& boost::none != initial_seed_ && boost::none!=Drift_distance_ &&
#ifndef _NO_CERN_ROOT
		(NULL!=random_generator_ || NULL != boost_random_generator)
#else //_NO_CERN_ROOT
		NULL != boost_random_generator
#endif //_NO_CERN_ROOT
		&& boost::none!= run_index_);
}

void Manager::InitTree (void)
{
#ifndef _NO_CERN_ROOT
	if (NULL== sim_data_)
		sim_data_ = new TTree("ElectronHistory", "ElectronHistory");
	const ProgramConstants *sets = gSettings.ProgConsts();

	if (sets->recorded_values.end()!=sets->recorded_values.find("process_type"))//only true values are stored
		sim_data_->Branch("process_type", &event_.process);
	if (sets->recorded_values.end()!=sets->recorded_values.find("time_initial"))
		sim_data_->Branch("time_initial", &event_.time_start);
	if (sets->recorded_values.end()!=sets->recorded_values.find("time_delta"))
		sim_data_->Branch("time_delta", &event_.delta_time);
	if (sets->recorded_values.end()!=sets->recorded_values.find("time_delta_full"))
		sim_data_->Branch("time_delta_full", &event_.delta_time_full);

	if (sets->recorded_values.end()!=sets->recorded_values.find("energy_initial"))
		sim_data_->Branch("energy_initial", &event_.En_start);
	if (sets->recorded_values.end()!=sets->recorded_values.find("energy_collision"))
		sim_data_->Branch("energy_coll", &event_.En_collision);
	if (sets->recorded_values.end()!=sets->recorded_values.find("energy_final"))
		sim_data_->Branch("energy_final", &event_.En_finish);
	if (sets->recorded_values.end()!=sets->recorded_values.find("energy_average"))
		sim_data_->Branch("energy_average", &event_.En_avr);

	if (sets->recorded_values.end()!=sets->recorded_values.find("position_initial"))
		sim_data_->Branch("position_initial",&event_.pos_start);
	if (sets->recorded_values.end()!=sets->recorded_values.find("position_delta"))
		sim_data_->Branch("position_delta",&event_.delta_x);
	if (sets->recorded_values.end()!=sets->recorded_values.find("position_final"))
		sim_data_->Branch("position_final",&event_.pos_finish);
	if (sets->recorded_values.end()!=sets->recorded_values.find("path_delta"))
		sim_data_->Branch("path_delta", &event_.delta_l);

	if (sets->recorded_values.end()!=sets->recorded_values.find("theta_initial"))
		sim_data_->Branch("theta_initial", &event_.theta_start);
	if (sets->recorded_values.end()!=sets->recorded_values.find("theta_collision"))
		sim_data_->Branch("theta_coll", &event_.theta_collision);
	if (sets->recorded_values.end()!=sets->recorded_values.find("theta_delta"))
		sim_data_->Branch("theta_delta", &event_.delta_theta);
	if (sets->recorded_values.end()!=sets->recorded_values.find("theta_final"))
		sim_data_->Branch("theta_final", &event_.theta_finish);

	if (sets->recorded_values.end()!=sets->recorded_values.find("photon_energy"))
		sim_data_->Branch("photon_energy", &event_.photon_En);

	if (sets->recorded_values.end()!=sets->recorded_values.find("deb_log_rand"))
		sim_data_->Branch("deb_log_rand",&event_.deb_log_rand);
	if (sets->recorded_values.end()!=sets->recorded_values.find("deb_solver_y_left"))
		sim_data_->Branch("deb_solver_y_left",&event_.deb_solver_y_left);
	if (sets->recorded_values.end()!=sets->recorded_values.find("deb_solver_y_right"))
		sim_data_->Branch("deb_solver_y_right",&event_.deb_solver_y_right);
	if (sets->recorded_values.end()!=sets->recorded_values.find("deb_solver_E_left"))
		sim_data_->Branch("deb_solver_E_left",&event_.deb_solver_E_left);
	if (sets->recorded_values.end()!=sets->recorded_values.find("deb_solver_E_right"))
		sim_data_->Branch("deb_solver_E_right",&event_.deb_solver_E_right);

	if (processes_data_) {;
		delete processes_data_;
	}
	processes_data_ = NULL;
#else //_NO_CERN_ROOT
	std::cerr << "Manager::InitTree: Warning: ROOT is turned off" << std::endl;
#endif //_NO_CERN_ROOT
}

void Manager::Clear(void)
{
#ifndef _NO_CERN_ROOT
	if (sim_data_) {
		delete sim_data_;
	}
	sim_data_ = NULL;
	if (processes_data_) {
		delete processes_data_;
	}
	processes_data_ = NULL;
#endif //_NO_CERN_ROOT
	skip_counter_ = 0;
	skipping_early_events = (gSettings.ProgConsts()->drift_distance_ignore_history ? true : false);
	num_of_events = 0;
	eField_ = boost::none;
	Concentration_ = boost::none;
	Coefficient_ = boost::none;
	Drift_distance_ = boost::none;
	initial_seed_ = boost::none;
	for (std::size_t i = 0, i_end_ = processes_counters_.size(); i!=i_end_; ++i)
		for (std::size_t j = 0, j_end_ = processes_counters_[i].size(); j!=j_end_; ++j)
			processes_counters_[i][j] = 0;
	InitTree();
}

void Manager::setParameters(double Concentr /*in SI*/, double E /*in SI*/, double drift_distance /*in m*/)
{
	Concentr = std::fabs(Concentr);
	E = std::fabs(E);
	if ((0 == Concentr) || (0 == E)) {
		std::cout << "Manager::SetParameters(): Can't have 0 values" << std::endl;
	}
	Concentration_ = Concentr;
	eField_ = E;
	if (drift_distance < 0) {
		std::cout << "Manager::SetParameters(): Warning: negative drift distance, using its absolute value." << std::endl;
	}
	Drift_distance_ = std::fabs(drift_distance);
	Coefficient_ = 1e20/* *e_charge_SIconst*/ * *eField_ / *Concentration_; // me/mu is set to 1.0
	//1e20 is because XS is expressed in 1e-20 m^2, not in m^2
}

void Manager::setParameters(double T /*in K*/, double Pressure /*in SI*/, double E /*in SI*/, double drift_distance /*in m*/)
{
	if (0 == T) {
		std::cout << "Manager::SetParameters(): Can't have 0 temperature" << std::endl;
		return;
	}
	setParameters(Pressure / (T*gSettings.PhysConsts()->boltzmann_SI), E, drift_distance);
}

bool Manager::setInitialSeed(unsigned long int seed)
{
	initial_seed_ = seed;
	#ifndef _NO_CERN_ROOT
	if (NULL != random_generator_) {
		switch (gSettings.ProgConsts()->random_generator) {
		case (ProgramConstants::GeneratorClass::TRand1): {
			UInt_t seedlist[2] = { *initial_seed_,0 }; //looked into TRandom1.cxx source code.
			//Simple SetSeed ===SetSeed2 set TRandom1 to fixed state, but this state is not defined by GetSeed().
			//in short GetSeed()===F(SetSeed()) which is defined, but its impossible to reproduce TRandom1 using
			((TRandom1 *)random_generator_)->SetSeeds(seedlist);
			//->SetSeed2(*initial_seed_);
			return true;
		}
		case (ProgramConstants::GeneratorClass::TRand2): {
			random_generator_->SetSeed(*initial_seed_);
			return true;
		}
		case (ProgramConstants::GeneratorClass::TRand3): {
			random_generator_->SetSeed(*initial_seed_);
			return true;
		}
		}
	}
	if (NULL != boost_random_generator && ProgramConstants::GeneratorClass::TRand1 == gSettings.ProgConsts()->random_generator) {
		boost_random_generator->seed(*initial_seed_);
		return true;
	}
#else //_NO_CERN_ROOT
	if (NULL != boost_random_generator) {
		boost_random_generator->seed(*initial_seed_);
		return true;
	}
#endif //_NO_CERN_ROOT
	initial_seed_ = boost::none;
	return false;
}

boost::optional<unsigned long int> Manager::getInitialSeed(void) const
{
	return initial_seed_;
}

bool Manager::setRunIndex(std::size_t index)
{
	run_index_ = index;
	return true;
}

boost::optional<std::size_t> Manager::getRunIndex(void) const
{
	return run_index_;
}

//Int XS(e)*sqrt(e/(e-Eny))*de
//'to' is always > 'from'
/*
long double Manager::XS_integral(long double from, long double to, long double Eny, Event &event)
{
	double E = from, E_prev = from;
	long double Int = 0;
	if ((from-Eny)/from<1e-6) {//irregularity case
		E = 1e-6*from + Eny;
		Int+=0.5*(ArTables_->TotalCrossSection(E)+ArTables_->TotalCrossSection(from))*sqrt(E_prev*(E_prev-Eny));
		from = E;
		E_prev = from;
	}
	ColoredRange energy_range_ = ColoredInterval (0, 0.1, 2e-4) + ColoredInterval (0.1, 1, 8e-4) +
		ColoredInterval (1, 10, 1e-2) + ColoredInterval (10, gSettings.PhysConsts()->XS_el_En_maximum, 0.02) +
		ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/5) +	//coarse area
		ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/5) +	//coarse area
		ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) + 	//fine area
		ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80) +		//fine area
		ColoredInterval (11.5, gSettings.PhysConsts()->XS_el_En_maximum, 0.003);
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
*/

double Manager::Uniform()
{
	//not called when not ready, so no checks
#ifndef _NO_CERN_ROOT
	if (NULL != random_generator_)
		return random_generator_->Uniform();
#endif //_NO_CERN_ROOT
	boost::random::uniform_01<boost::random::hellekalek1995> dist(*boost_random_generator);
	return dist();
}

long double Manager::XS_integral_for_test(long double from, long double to, long double Eny, long double dE)
{
	double E = from, E_prev = from;
	long double Int = 0;
	const Particle * particle = gParticleTable.GetParticle(ELECTRON_NAME);
	if ((from - Eny) / from<1e-6) {//irregularity case
		E = std::min(1e-6*from + Eny, to);
		Int += 0.5*(material_->GetCrossSection(particle, E) + material_->GetCrossSection(particle, from))*sqrt(E*(E - Eny));
		from = E;
		E_prev = from;
	}
	if (from==to)
		return Int;
	while (E<to) {
		Int += material_->GetCrossSection(particle, E)*sqrt(E / (E - Eny))*dE;
		E_prev = E;
		E += dE;
	}
	E = to;
	Int += material_->GetCrossSection(particle, E)*sqrt(E / (E - Eny))*(E - E_prev);
	return Int;
}

void Manager::Initialize(void)
{
	if (!isReady()) {
		std::cout << "Manager::Initialize: some of the parameters are not initialized, exiting" << std::endl;
		return;
	}
	skip_counter_ = 0;
	skipping_early_events = (gSettings.ProgConsts()->drift_distance_ignore_history ? true : false);
	num_of_events = 0;
	skip_time_residue_ = 0;
	//processes_counters; //preserve
#ifndef _NO_CERN_ROOT
	if (NULL != random_generator_) {
		e_first_seed_ = random_generator_->GetSeed();
	} else
#endif //_NO_CERN_ROOT
	{
		std::stringstream ss;
		ss << *boost_random_generator;
		std::string str = ss.str();
		e_first_seed_ = std::stoul(str);
	}
}

void Manager::InitializeEvent()
{
	if (!isReady()) {
		std::cout << "Manager::Initialize: some of the parameters are not initialized, exiting" << std::endl;
		return;
	}
	const boost::optional<PDF_routine> *Ec_spec = &(gSettings.ProgConsts()->run_specifics[*run_index_].Ec_spectrum);
	if (*Ec_spec) {
		event_.En_collision = (**Ec_spec)(Uniform());
		event_.En_start= 0;
	} else {
		event_.En_start = 1 + 3*Uniform();
		event_.En_collision = 0;
	}
	event_.En_finish = 0;
	event_.pos_start = 0;
	event_.pos_finish = 0;
	event_.delta_x = 0;
	event_.time_start = 0;
	event_.delta_time = 0;
	event_.delta_time_full = 0;
	event_.process = Event::None;
	event_.particle_ID = gParticleTable.GetParticleID(ELECTRON_NAME);

	event_.theta_start = 0;
	event_.theta_collision = 0;
	event_.theta_finish = 0;
	event_.delta_theta = 0;
	event_.delta_l = 0;

	event_.photon_En = 0;
}

//TODO: maybe move solving to material
void Manager::Solve (long double LnR, Event &event)
{
	event.deb_log_rand = LnR;
	long double Eny = event.En_start*sin(event.theta_start)*sin(event.theta_start);
	long double e_start = event.En_start;
	short case_ = 0; //0 - normal Vx_start>0, 1 - Vx_start<0, no sign change, 2 - Vx_start<0, sign is changed
	const Particle* particle = gParticleTable.GetParticle(event.particle_ID);
	if (event.theta_start>M_PI/2) {
		double TURN_INT = material_->GetUntabXSIntegral(particle, Eny, event.En_start, Eny);
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
	long double right = e_start + 10*LnR/(material_->GetCrossSection(particle, left)*sqrt(left/(left-Eny)));
	right = std::min((double)right, gSettings.ProgConsts()->maximal_energy);
	long double I_max = material_->GetUntabXSIntegral(particle, e_start, right, Eny);
	if (I_max < LnR) {
		I_max = material_->GetUntabXSIntegral(particle, e_start, gSettings.ProgConsts()->maximal_energy, Eny);
		right = gSettings.ProgConsts()->maximal_energy;
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
			f_new = material_->GetUntabXSIntegral(particle, e_start, e_finish, Eny) - LnR;
			if (f_new < 0) {
				left = e_finish;
				f_left = f_new;
			} else {
				right = e_finish;
				f_right = f_new;
			}
			convergence_criteria = 2e-4*(std::min(std::max(e_start, (long double)0.1*gSettings.PhysConsts()->XS_el_En_minimum),
					std::max(e_finish, (long double)0.1*gSettings.PhysConsts()->XS_el_En_minimum)));//std::max(2e-6, std::fabs(5e-4*event.En_finish));
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
	if (isnan(event.theta_collision)) {
		event.theta_collision = 0;
	}
}

void Manager::Solve_table (long double LnR, Event &event)
{
	event.deb_log_rand = LnR;
	long double Eny = event.En_start*sin(event.theta_start)*sin(event.theta_start);
	short case_ = 0; //0 - normal Vx_start>0, 1 - Vx_start<0, no sign change, 2 - Vx_start<0, sign is changed
	const Particle* particle = gParticleTable.GetParticle(event.particle_ID);
	double INT = material_->GetXSIntegral(particle, Eny, event.En_start, Eny);
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
	event.En_collision = material_->GetEnergyFromXSIntegral(particle, Eny, INT); //TODO: add debug info - uncertainty and overflow - status output variable
	if (event.En_collision<0) {
		std::cout<<"Manager::Solve_table::Error: found E<0 in table. Eny= "<<Eny<<" Ei= "<<event.En_start<<" Int= "<<INT<<std::endl;
		//event.En_collision = material_->GetEnergyFromXSIntegral(particle, Eny, INT);
		event.En_collision = std::min(0.01*event.En_start, 0.001) + event.En_start;
		event.En_collision = std::min(event.En_collision, gSettings.ProgConsts()->maximal_energy);
	}
	if (1==case_) {
		event.theta_collision = acos(-sqrt((event.En_collision-Eny)/event.En_collision));
	} else {
		event.theta_collision = acos(sqrt((event.En_collision-Eny)/event.En_collision));
	}
	if (isnan(event.theta_collision)) {
		event.theta_collision = 0;
		//event.En_collision = material_->GetEnergyFromXSIntegral(particle, Eny, INT);
	}
}

//Solve with high accuracy integral
void Manager::Solve_test (long double LnR, Event &event)
{
	event.deb_log_rand = LnR;
	long double Eny = event.En_start*sin(event.theta_start)*sin(event.theta_start);
	long double e_start = event.En_start;
	short case_ = 0; //0 - normal Vx_start>0, 1 - Vx_start<0, no sign change, 2 - Vx_start<0, sign is changed
	const Particle* particle = gParticleTable.GetParticle(event.particle_ID);
	if (event.theta_start > M_PI/2) {
		double TURN_INT = material_->GetUntabXSIntegral(particle, Eny, event.En_start, Eny);
		if (LnR > TURN_INT) { //Vx changes its sign.
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
	long double right = e_start + 10*LnR/(material_->GetUntabCrossSection(particle, left)*sqrt(left/(left-Eny)));
	right = std::min((double)right, gSettings.ProgConsts()->maximal_energy);
	long double I_max = XS_integral_for_test(e_start, right, Eny, 1e-5);
	if (I_max < LnR) {
		I_max = XS_integral_for_test(e_start, gSettings.ProgConsts()->maximal_energy, Eny, 1e-5);
		right = gSettings.ProgConsts()->maximal_energy;
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
			convergence_criteria = 2e-4*(std::min(std::max(e_start, (long double)0.1*gSettings.PhysConsts()->XS_el_En_minimum),
					std::max(e_finish, (long double)0.1*gSettings.PhysConsts()->XS_el_En_minimum)));//std::max(2e-6, std::fabs(5e-4*event.En_finish));
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

void Manager::DoStepLength(void) //In case Ec spectrum is fixed, stepping is done reversed in time
{
	if (!isReady())
		return;
	bool reverse_mode = (event_.En_collision!=0);
	Event step_pars = event_;
	if (reverse_mode) {
		step_pars.theta_start = M_PI - event_.theta_collision;
		step_pars.En_start = event_.En_collision;
	}
	long double L = - log(Uniform());
	L *= *Coefficient_; //Calculated once for fixed parameters;
	//solving L = XS_integral(Ei, Ec) for Ec===E collision.
	Solve_table(L, step_pars);

	if (reverse_mode) {
		event_.En_start = step_pars.En_collision;
		event_.theta_start = M_PI - step_pars.theta_collision;
	} else {
		event_.En_collision = step_pars.En_collision;
		event_.theta_collision = step_pars.theta_collision;
	}
	double eom = gSettings.PhysConsts()->e_charge_SI / gSettings.PhysConsts()->e_mass_SI;
	long double vel_0 = sqrt(2.0*eom*event_.En_start)*cos(event_.theta_start);
	long double vel_1 = sqrt(2.0*eom*event_.En_collision)*cos(event_.theta_collision);
	event_.delta_time = (vel_1 - vel_0) / (eom * *eField_); //in s
	event_.delta_x = (event_.En_collision - event_.En_start) / *eField_;
	if (reverse_mode) {
		event_.delta_x *= -1;
		event_.delta_time *= 1;
	}
	event_.pos_finish = event_.pos_start + event_.delta_x;
	event_.En_avr = event_.En_start + *eField_*vel_0*event_.delta_time / 4.0 + eom*std::pow(*eField_*event_.delta_time, 2) / 6.0;
	//Total path in meters:
	//calculated as Int{sqrt((dy/dt)^2 + (dx/dt)^2 ) dt} from 0 to delta_time
	long double up_limit = event_.delta_time**eField_*sqrt(eom/(2*event_.En_start));
	event_.delta_l = Path_integral(up_limit, event_.theta_start) - Path_integral(0, event_.theta_start);
	event_.delta_l*=2*step_pars.En_start/ *eField_;
}

void Manager::DoScattering(void)
{
	if (!isReady())
		return;
	bool is_overflow = (event_.process==Event::Overflow);
	const Particle* incident = gParticleTable.GetParticle(event_.particle_ID);
	double R2 = Uniform(); //R1 is used in DoStepLength
	const Particle* target = material_->GenerateScatteringParticle(incident, event_.En_collision, R2);
	//TODO: implement exception system
	if (NULL == target) {
		std::cerr<<"Attempting again"<<std::endl;
		R2 = Uniform(); //R1 is used in DoStepLength
		target = material_->GenerateScatteringParticle(incident, event_.En_collision, R2);
		if (NULL == target) {
			std::cerr<<"Fallback to the most abundant element in the mixture"<<std::endl;
			target = material_->GetDominatingParticle(incident, event_.En_collision);
			if (NULL == target) {
				//throw;
				return;
			}
		}
	}
	double R3 = Uniform();
	event_.process = target->GenerateProcess(incident, event_.En_collision, R3);
	if (Event::None == event_.process) {
		//throw
		return;
	}
	double R4 = Uniform();
	event_.delta_theta = target->GenerateScatterAngle(incident, event_.En_collision, event_.process, R4);
	double phi = Uniform()*2.0*M_PI;
	double cos_th_f = std::cos(event_.delta_theta)*std::cos(event_.theta_collision) + std::sin(event_.delta_theta)*std::sin(event_.theta_collision)*std::cos(phi);
	if (cos_th_f<-1.0) //just in case of precision problems
		cos_th_f = -1.0;
	if (cos_th_f>1.0)
		cos_th_f = 1.0;
	event_.theta_finish = std::acos(cos_th_f);
	if (std::isnan(event_.theta_finish)) {
		event_.theta_finish = 0;
	}
	double R5 = Uniform();
	double EnergyLoss = target->GenerateEnergyLoss(incident, event_.En_collision, event_.delta_theta, event_.process, R5);
	event_.photon_En = target->GeneratePhoton(incident, event_.En_collision, event_.delta_theta, event_.process, R5); //same random value must be used!
	event_.En_finish = event_.En_collision - EnergyLoss;
	double R6 = Uniform();
	double time_delay = target->GenerateTimeDelay(incident, event_.En_collision, event_.delta_theta, event_.process, R6);
	event_.delta_time_full = event_.delta_time + time_delay;
	std::vector<const Particle*> outputs = target->GetFinalStates(incident, event_.En_collision, event_.delta_theta, event_.process);
	if (outputs.size()<2) {
		//throw
		return;
	}
	event_.particle_ID_finish = gParticleTable.GetParticleID(outputs[1]);
	event_.particle_ID_target = gParticleTable.GetParticleID(outputs[0]);
	if (is_overflow)
		event_.process = Event::Overflow;
}

void Manager::PostStepAction(void)
{
	if (!isReady())
		return;
	bool do_write = false;
	//Operations which are always performed regardless of writing configuration (every nth event, every x seconds, ignore events before e reached fixed Z):
	std::size_t e_ID = gParticleTable.GetParticleID(ELECTRON_NAME);
	for (int i = 0, i_end_ = processes_counters_.size(); i != i_end_; ++i) {
		if (i == event_.particle_ID_target && e_ID == event_.particle_ID) {
			for (int j = 0, j_end_ = processes_counters_[i].size(); j!=j_end_; ++j) {
				if (processes_IDs_[i][j] == event_.process) {
					++processes_counters_[i][j];
					break;
				}
			}
		}
	}
	if ((0==num_of_events)||IsFinished()) //always record first and last drift event
		do_write = true;
	//End of operations which are always performed
	//Logic of different writing modes:
	if (gSettings.ProgConsts()->skip_history_time) { //ignores drift_distance_ignore_history and skip_history_rate
		if (event_.delta_time_full >=*gSettings.ProgConsts()->skip_history_time) {
			//Event current = event_;
			double delta = event_.delta_time_full;
			while (delta >=*gSettings.ProgConsts()->skip_history_time) {
				event_.delta_time_full = *gSettings.ProgConsts()->skip_history_time;
				if (event_.delta_time > event_.delta_time_full)
					event_.delta_time_full = event_.delta_time;
				sim_data_->Fill();
				delta -= event_.delta_time_full;
				event_.time_start += event_.delta_time_full;
				event_.delta_time = 0;
			}
			event_.delta_time_full = delta;
		}
		skip_time_residue_ += event_.delta_time_full;
		if (skip_time_residue_ >= *gSettings.ProgConsts()->skip_history_time) {
			skip_time_residue_ -= *gSettings.ProgConsts()->skip_history_time;
			do_write = true;
		}
	} else {
		if (skipping_early_events) {  //if true (set at initialization) all events with position less than
			//gSettings.ProgConsts()->drift_distance_ignore_history (if present) are not recored.
			if (event_.pos_start >= gSettings.ProgConsts()->drift_distance_ignore_history) {
				skipping_early_events = false;
				skip_counter_ = 0;
			}
		}
		if (skip_counter_ >= gSettings.ProgConsts()->skip_history_rate || boost::none == gSettings.ProgConsts()->skip_history_rate)
			skip_counter_ = 0;
		if (!skipping_early_events && 0 == skip_counter_)
			do_write = true;
	}
	//End of logic
#ifndef _NO_CERN_ROOT
	//prev_event_ = event_;
	//prev_event_saved_ = false;
	if (do_write) {
		sim_data_->Fill();
		prev_event_saved_ = true;
	}
#endif //_NO_CERN_ROOT
	++num_of_events; //May overflow
	++skip_counter_;
}

void Manager::DoGotoNext(void)
{
	if (Event::None != event_.process) { //== for the very first event
		const boost::optional<PDF_routine> *Ec_spec = &(gSettings.ProgConsts()->run_specifics[*run_index_].Ec_spectrum);
		if (*Ec_spec) {
			event_.En_collision = (**Ec_spec)(Uniform());
			double Rand = Uniform()*2.0 - 1.0;
			event_.theta_collision = std::acos(Rand); //Actually thete_start is uniform but theta collision is not (but close)
			event_.En_start = 0;
		} else {
			event_.En_start = event_.En_finish;
			event_.theta_start = event_.theta_finish;
			event_.En_collision = 0;
		}
		event_.time_start += event_.delta_time_full;
		event_.pos_start = event_.pos_finish;
		event_.process = Event::None;
		event_.delta_time = 0;
		event_.delta_x = 0;
		event_.delta_time_full = 0;
		event_.photon_En = 0;
		event_.particle_ID = event_.particle_ID_finish;
	}
}

void Manager::DoStep(void)
{
	if (!isReady())
		return;
	DoGotoNext();
	DoStepLength(); //In case Ec spectrum is fixed, stepping is done reversed in time
	DoScattering();
	PostStepAction();
}

bool Manager::IsFinished(void)
{
	if (!isReady())
		return true;
	//return sim_data_->GetEntries()>=100;
	return !(event_.pos_finish < Drift_distance_);
}

void Manager::LoopSimulation(void)
{
	Initialize();
	if (!isReady())
		return;
	if (!material_->isValid()) {
		std::cerr << "Manager::LoopSimulation: invalid material" << std::endl;
		return;
	}
	InitializeEvent();
	while (!IsFinished()) {
		DoStep();
	}
}

void Manager::WriteHistory(std::string root_fname)
{
#ifndef _NO_CERN_ROOT
	TFile *file = new TFile(root_fname.c_str(), "RECREATE");
	file->cd();
	std::cout<<"Event number: "<<sim_data_->GetEntries()<<std::endl;
	sim_data_->Write("", TObject::kOverwrite);
	if (NULL != processes_data_) {
		delete processes_data_;
	}
	processes_data_ = new TTree("ElectronProcessCounters", "ElectronProcessCounters");

	//ROOT Trees require pointers, not a standard containers
	Long64_t *processes_counters; //per process_ID
	Short_t *processes_IDs;
	char **processes_legends; //per process_ID
	UInt_t processes_size; //particle_ID
	UInt_t particle_ID;
	char * particle_name;
	//std::vector<std::string> *processes_names;

	processes_data_->Branch("particle_ID", &particle_ID, "particle_ID/i");
	processes_data_->Branch("proc_size", &processes_size, "proc_size/i");
	TBranch *brParticleName = processes_data_->Branch("particle_name", particle_name, "particle_name/C");
	TBranch *brProcessesID = processes_data_->Branch("proc_IDs", processes_IDs, "proc_IDs[proc_size]/S");
	TBranch *brProcessesCounters = processes_data_->Branch("proc_Ns", processes_counters, "proc_Ns[proc_size]/L");
	TBranch *brProcessesName = processes_data_->Branch("proc_names", processes_legends, "proc_names[proc_size]/C");
	//TBranch *brProcessesName1 = processes_data_->Branch("proc_names1", &processes_names);

	for (std::size_t i = 0, i_end_ = processes_counters_.size(); i!=i_end_; ++i) {
		particle_ID = i;
		particle_name = c_str_cp(gParticleTable.GetParticle(i)->GetName());
		processes_size = processes_counters_[i].size();
		processes_counters = new Long64_t [processes_size];
		processes_IDs = new Short_t [processes_size];
		processes_legends = new char* [processes_size];
		//processes_names = new std::vector<std::string> [processes_size];
		for (std::size_t proc = 0; proc != processes_size; ++proc) {
			processes_counters[proc] = processes_counters_[i][proc];
			processes_IDs[proc] = processes_IDs_[i][proc];
			processes_legends[proc] = c_str_cp(processes_legends_[i][proc]);
			//(*processes_names)[proc] = processes_legends_[i][proc];
		}
		brParticleName->SetAddress(particle_name);
		brProcessesID->SetAddress(processes_IDs);
		brProcessesCounters->SetAddress(processes_counters);
		brProcessesName->SetAddress(processes_legends);
		//brProcessesName1->SetAddress(&processes_names);
		processes_data_->Fill();
		for (std::size_t proc = 0; proc != processes_size; ++proc) {
			delete [] processes_legends[proc];
		}
		delete [] processes_counters;
		delete [] processes_IDs;
		delete [] processes_legends;
		//delete [] processes_names;
		delete [] particle_name;
	}

	processes_data_->Write("", TObject::kOverwrite);
	file->Close(); //should delete both trees
	delete file;
	//delete processes_data_; //This ROOT memory management, argh!
	processes_data_ = NULL;
	sim_data_ = NULL;
#else //_NO_CERN_ROOT
	std::cerr << "Manager::WriteHistory: Warning: ROOT is turned off" << std::endl;
#endif //_NO_CERN_ROOT
}

void Manager::WriteLegend(std::string legend_fname)
{
	std::ofstream str;
	open_output_file(legend_fname, str, std::ios_base::trunc);
	if (!str.is_open())
		return;
	for (std::size_t i = 0, i_end_ = processes_counters_.size(); i!=i_end_; ++i) {
		UInt_t processes_size = processes_counters_[i].size();
		str<<"Particle_ID: "<<i<<"\tName: \""<<gParticleTable.GetParticle(i)->GetName()<<"\"\tN of procesess: "<<processes_size<<std::endl;
		str<<"Processes:"<<std::endl;
		for (std::size_t proc = 0; proc != processes_size; ++proc) {
			str<<"ID: "<<processes_IDs_[i][proc]<<"\tName: "<<processes_legends_[i][proc]<<"\tCounter: "<<processes_counters_[i][proc]<<std::endl;
		}
		str<<"=========================================="<<std::endl;

	}
	str.close();
}

/*
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
*/
