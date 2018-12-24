#include "Manager.h"

Manager::Manager(ArDataTables *Ar_tables, UInt_t RandomSeed) : is_ready_(false), eField_(0), Concentration_(0), Coefficient_(0),
	skip_counter_(0), ArTables_(Ar_tables)
{
	random_generator_ = new TRandom1(RandomSeed); //TRandom3 is a default. TRandom1 is slower but better
	InitTree();
}

void Manager::InitTree (void)
{
	sim_data_ = new TTree("ElectronHistory", "ElectronHistory");
	sim_data_->Branch("energy_initial", &event_.En_start);
	sim_data_->Branch("energy_coll", &event_.En_collision);
	sim_data_->Branch("energy_final", &event_.En_finish);
	sim_data_->Branch("energy_average", &event_.En_avr);

	sim_data_->Branch("time_initial", &event_.time_start);
	sim_data_->Branch("time_delta", &event_.delta_time);
	sim_data_->Branch("time_delta_full", &event_.delta_time_full);
	sim_data_->Branch("process_type", &event_.process);

	sim_data_->Branch("position_initial",&event_.pos_start);
	sim_data_->Branch("position_final",&event_.pos_finish);
	sim_data_->Branch("position_delta",&event_.delta_x);

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
}

void Manager::Clear(void)
{
	is_ready_ = false;
	sim_data_->Delete();
	skip_counter_ = 0;
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

void Manager::Initialize(Event &event)
{
	if (!is_ready_)
		return;
	skip_counter_=0;
	event.CrossSections.resize(ArTables_->ArAllData_.ArExper_.max_process_ID + 2, 0); //2==Elastic + Resonance
	event.CrossSectionsSum.resize(ArTables_->ArAllData_.ArExper_.max_process_ID + 2, 0);
	event.En_start = 1 + 3*random_generator_->Uniform();
	//event.En_start = 3;
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
			convergence_criteria = 2e-4*(std::min(e_start, e_finish));//std::max(2e-6, std::fabs(5e-4*event.En_finish));
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
	Solve(L, event);

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
	event.theta_finish = event.delta_theta + event.theta_start;
	if (event.theta_finish>M_PI)
		event.theta_finish = 2*M_PI - event.theta_finish;
	long double gamma_f = e_mass_eVconst/Ar_mass_eVconst;
	double EnergyLoss = 2*(1-cos(event.delta_theta))*event.En_collision*gamma_f /pow(1 + gamma_f, 2);
	switch (event.process) {
		case (Event::Resonance_3o2): {
			EnergyLoss *=RESONANCE_EN_LOSS_FACTOR_;
			break;
		}
		case (Event::Resonance_1o2): {
			EnergyLoss *=RESONANCE_EN_LOSS_FACTOR_;
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
	if ((0==skip_counter_)||(IsFinished(event))||(event.process>Event::Elastic)) {
		sim_data_->Fill();
		skip_counter_=0;
	}
	++skip_counter_;
	if (skip_counter_>SKIP_HISTORY_)
		skip_counter_=0;
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
	file->Close();
	std::ofstream str;
	root_fname.pop_back();
	root_fname.pop_back();
	root_fname.pop_back();
	root_fname.pop_back();
	root_fname.pop_back();
	str.open(root_fname+"process_legend.txt", std::ios_base::trunc);
	str<<"//tree->process value | meaning"<<std::endl;
	str<<Event::Overflow<<"\t\"Overflow (energy exceeds "<<EN_MAXIMUM_<<" eV)\""<<std::endl;
	str<<Event::None<<"\t\"None\""<<std::endl;
	str<<Event::Elastic<<"\t\"Elastic scattering\""<<std::endl;
	str<<Event::Resonance_3o2<<"\t\"Feshbach resonance 3/2\""<<std::endl;
	str<<Event::Resonance_1o2<<"\t\"Feshbach resonance 1/2\""<<std::endl;
	for (int proc = Event::Ionization, end_ = ArTables_->ArAllData_.ArExper_.max_process_ID + Event::Ionization; proc!=end_; ++proc) {
		InelasticProcess *p = ArTables_->ArAllData_.ArExper_.FindInelastic(proc-Event::Ionization);
		if (NULL!=p)
			str<<proc<<"\t"<<p->get_name()<<std::endl;
	}
	str.close();
}

void Manager::Test(void)
{
	int Order = ArTables_->getOrder();
	int Nused = ArTables_->getNused();

	std::string fname_XS_3_4 = "tests/Int_XS_3_4.txt";
	std::string fname_XS_high_3_4 = "tests/Int_XS_high_3_4.txt";
	std::string fname_XS_1_2 = "tests/Int_XS_1_2.txt";
	std::string fname_XS_high_1_2 = "tests/Int_XS_high_1_2.txt";

	ArTables_->setOrder(3);
	ArTables_->setNused(4);
	std::ofstream str;
	str.open(fname_XS_3_4, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	double E;
	for (unsigned int l=0; l<50001;++l) {
		E = -EN_MAXIMUM_ + l*2*EN_MAXIMUM_/50000.0;
		//str<< E<<"\t"<<XS_integral(-EN_MAXIMUM_, E)<<std::endl;
	}
	str.close();

	str.open(fname_XS_high_3_4, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	double Int = 0;
	E = -EN_MAXIMUM_;
	double dE = 1e-5;
	while (E<EN_MAXIMUM_) {
		str<< E<<"\t"<<Int<<std::endl;
		Int+=ArTables_->XS_elastic(std::fabs(E))*dE;
		E+=dE;
	}
	str.close();

	ArTables_->setOrder(1);
	ArTables_->setNused(2);
	str.open(fname_XS_1_2, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	for (unsigned int l=0; l<50001;++l) {
		E = -EN_MAXIMUM_ + l*2*EN_MAXIMUM_/50000.0;
		//str<< E<<"\t"<<XS_integral(-EN_MAXIMUM_, E)<<std::endl;
	}
	str.close();

	str.open(fname_XS_high_1_2, std::ios_base::trunc);
	str<<"//E[eV]\tInt{XS_elastic}"<<std::endl;
	Int = 0;
	E = -EN_MAXIMUM_;
	dE = 1e-5;
	while (E<EN_MAXIMUM_) {
		str<< E<<"\t"<<Int<<std::endl;
		Int+=ArTables_->XS_elastic(std::fabs(E))*dE;
		E+=dE;
	}
	str.close();

	std::string name = std::string("tests/test_XS_integral.sc");
	str.open(name, std::ios_base::trunc);
	str<<"plot '"<<fname_XS_3_4<<"' u 1:2 title 'I(XS) 3,4'"<<std::endl;
	str<<"replot '"<<fname_XS_high_3_4 <<"' u 1:2 w line lc rgb \"#000000\" title 'I(XS) high acc. 3,4'"<<std::endl;
	str<<"replot '"<<fname_XS_1_2<<"' u 1:2 title 'I(XS) 1,2'"<<std::endl;
	str<<"replot '"<<fname_XS_high_1_2 <<"' u 1:2 w line lc rgb \"#FF0000\" title 'I(XS) high acc. 1,2'"<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);

	ArTables_->setOrder(Order);
	ArTables_->setNused(Nused);
}
