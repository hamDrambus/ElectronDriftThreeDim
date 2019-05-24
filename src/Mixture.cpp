#include "Mixture.h"

Mixture::Mixture(std::string mixture_name):
	integral_table_fname(gSettings.ProgConsts()->tabulated_data_folder + "cross_integrals.dat")
{
	ensure_file(integral_table_fname);
	integral_table_ = new FunctionTable; //shared between processes - read only
	is_valid_ = false;
	name_ = mixture_name;
}

Mixture::~Mixture() {};

void Mixture::AddComponent(const Particle* part, double fraction) {
	if (NULL==part)
		return;
	bool duplicate = false;
	for (auto c = components_.begin(), c_end_ = components_.end(); c != c_end_; ++c) {
		if ((c->first == part) || (c->first->GetName() == part->GetName())) {
			c->second = std::max(fraction, 0.0);
			duplicate = true;
			break;
		}
	}
	if (!duplicate)
		components_.push_back(std::pair<const Particle*, double>(part, std::max(fraction, 0.0)));
}

void Mixture::AddComponent(std::string name, double fraction) {
	const Particle* part = gParticleTable.GetParticle(name);
	if (NULL==part) {
		std::cerr<<"Error Mixture(\"" << name_ << "\")::AddComponent: no such particle ("<<name<<")"<<std::endl;
		return;
	}
	AddComponent(part, fraction);
}

//Fallback function. Useless but just in case
const Particle* Mixture::GetDominatingParticle(const Particle *incident, double E) const
{
	if (!isValid())
		return NULL;
	bool use_concentr = false;
	if (NULL == incident) {
		use_concentr = true;
	}
	std::size_t c_end_= components_.size();
	std::size_t c_max = 0;
	double v_max = 0;
	if (!use_concentr) {
		for (std::size_t c = 0; c != c_end_; ++c) {
			double val = 0;
			for (unsigned int pr = 0, pr_end_ = components_[c].first->GetProcessSize(incident); pr != pr_end_; ++pr)
				val += components_[c].first->GetCrossSection(incident, E, pr) * components_[c].second;
			if (val>v_max) {
				v_max = val;
				c_max = c;
			}
		}
		if (0 == v_max) {
			use_concentr = true;
			c_max = 0;
			v_max = 0;
		}
	}
	if (use_concentr) {
		for (std::size_t c = 0; c != c_end_; ++c) {
			double val = components_[c].second;
			if (val>v_max) {
				v_max = val;
				c_max = c;
			}
		}
		if (0 == v_max) {
			return NULL;
		};
	}
	return components_[c_max].first;
}

void Mixture::Prepare(void) {
	while (true) {
		for (auto i = components_.begin(), i_end_ = components_.end(); i != i_end_; ++i)
			if (i->second <= 0) {
				std::cerr << "Mixture(\"" << name_ << "\")::Prepare: Warning! Negative or zero concentration of \"" \
					<< i->first->GetName() << "\" is removed." << std::endl;
				components_.erase(i);
				continue;
			}
		break;
	}
	double Norm = 0;
	for (auto i = components_.begin(), i_end_ = components_.end(); i != i_end_; ++i)
		Norm += i->second;
	for (auto i = components_.begin(), i_end_ = components_.end(); i != i_end_; ++i)
		i->second /= Norm;
	is_valid_ = !components_.empty();
	bool valid_children = true;
	for (auto c = components_.begin(), c_end_ = components_.end(); c != c_end_; ++c) {
		if (!c->first->isValid()) {
			valid_children = false;
			break;
		}
	}
	is_valid_ = is_valid_ && valid_children;
	if (!isValid()) {
		std::cerr<<"Mixture(\"" << name_ << "\")::Prepare: Mixture in invalid"<<std::endl;
		return;
	}
	std::ifstream inp;
	inp.open(integral_table_fname, std::ios_base::binary);
	if (!inp.is_open()) { //TODO: add identification of table corresponding to certain mixture
		std::cout<<"Failed to load \""<<integral_table_fname<<"\""<<std::endl;
		if (generate_integral_table())
			integral_table_->write(integral_table_fname);
		else
			is_valid_ = false;
	} else {
		integral_table_->read(inp);
		inp.close();
		if (integral_table_->is_empty()) {
			std::cout<<"Failed to load \""<<integral_table_fname<<"\""<<std::endl;
			if (generate_integral_table())
				integral_table_->write(integral_table_fname);
			else
				is_valid_ = false;
		}
	}
}

bool Mixture::generate_integral_table(void)
{
	if (!isValid()) {
		return false;
	}
	const Particle *electron = gParticleTable.GetParticle("electron");
	if (NULL == electron) {
		std::cout<<"Mixture(\"" << name_ << "\"): Error: could not find electron particle"<<std::endl;
		return false;
	}
	std::cout<<"Mixture(\"" << name_ << "\"): Generating cross section integrals..."<<std::endl;
	ColoredRange energy_Y_range;
	for (auto c = components_.begin(), c_end_ = components_.end(); c != c_end_; ++c)
		energy_Y_range += c->first->GetEnergies();
	energy_Y_range.Trim(0, gSettings.ProgConsts()->maximal_energy);
	for (long int Ey_i=0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i!=Ey_ind_end_;++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		ColoredRange energy_range = energy_Y_range;
		energy_range = energy_range + ColoredInterval(Ey, 1.1*Ey, Ey/500.0)
				+ ColoredInterval(1.1*Ey, 10*Ey, Ey/100.0) + ColoredInterval(10*Ey, 50*Ey, Ey/10.0);
		energy_range.Trim(Ey, gSettings.ProgConsts()->maximal_energy);
		long double Int = 0;
		double E = Ey, E_prev = Ey;
		for (long int E_i=0, E_i_end_=energy_range.NumOfIndices(); E_i!=E_i_end_;++E_i) {
			E = energy_range.Value(E_i);
			if (E<=Ey) {
				E_prev = E;
				integral_table_->push(Ey, Ey, 0);
				continue;
			}
			if (E_prev == Ey) {//irregularity case
				Int += 0.5*(GetCrossSection(electron, E) + GetCrossSection(electron, E))*sqrt(E*(E - Ey));
			} else {
				Int+=GetCrossSection(electron, E)*sqrt(E/(E-Ey))*(E-E_prev);
			}
			E_prev = E;
			integral_table_->push(E, Ey, Int);
		}
	}
	return true;
}

double Mixture::GetCrossSection(const Particle *incident, double E) const
{
	if (!isValid()) {
		std::cout<<"Mixture(\"" << name_ << "\")::GetCrossSection: Error: invalid mixture"<<std::endl;
		return 0;
	}
	double out = 0;
	for (auto c = components_.begin(), c_end_ = components_.end(); c != c_end_; ++c)
		for (unsigned int pr = 0, pr_end_ = c->first->GetProcessSize(incident); pr != pr_end_; ++pr)
			out += c->first->GetCrossSection(incident, E, pr) * c->second;
	return out;
}

long double Mixture::GetXSIntegral(const Particle *incident, long double En_from, long double En_to, long double En_y) const
{
	if (NULL == incident)
		return 0;
	if (incident->GetName()=="electron")
		return (*integral_table_)(En_to, En_y) - (*integral_table_)(En_from, En_y);
	else
		return 0; //TODO: Implement general cases
}

long double Mixture::GetEnergyFromXSIntegral(const Particle *incident, long double En_y, long double integral_value) const
{
	if (NULL == incident)
		return DBL_MAX;
	if (incident->GetName()=="electron")
		return integral_table_->find_E(En_y, integral_value);
	else
		return DBL_MAX; //TODO: Implement general cases
}

const Particle* Mixture::GenerateScatteringParticle(const Particle *incident, double E, double Rand) const
{
	if (!isValid()) {
		std::cout<<"Mixture(\"" << name_ << "\")::GenerateScatteringParticle: Error: invalid mixture"<<std::endl;
		return NULL;
	}
	std::size_t c_end_= components_.size();
	std::vector<double> CrossSections(c_end_, 0.0), CrossSectionsSum(c_end_, 0.0);
	//TODO: allocating memory each time is quite expensive. Similar issue for Particle. Need to create cache
	//All sizes and particle interations are static (so far), so vectors can be resized for each incident particle
	for (std::size_t c = 0; c != c_end_; ++c) {
		for (unsigned int pr = 0, pr_end_ = components_[c].first->GetProcessSize(incident); pr != pr_end_; ++pr)
			CrossSections[c] += components_[c].first->GetCrossSection(incident, E, pr) * components_[c].second;
		CrossSectionsSum[c] = std::max(CrossSections[c], 0.0) + ((c==0) ? 0.0 : CrossSectionsSum[c - 1]);
	}
	for (std::size_t c = 0; c != c_end_; ++c) {
		CrossSectionsSum[c] /= CrossSectionsSum[c_end_ - 1];
		if (Rand < CrossSectionsSum[c]) {
			return components_[c].first;
		}
	}
	std::cerr <<"Mixture(\"" << name_ << "\")::GenerateScatteringParticle: Error: failed to select particle from "<<c_end_<<" candidates"<<std::endl;
	return NULL;
}

double Mixture::GetUntabCrossSection(const Particle *incident, double E) const
{
	return GetCrossSection(incident, E);
}
long double Mixture::GetUntabXSIntegral(const Particle *incident, long double En_from, long double En_to, long double En_y) const
{
	if (NULL == incident)
		return 0;
	if (incident->GetName()=="electron") {
		double E = En_from, E_prev = En_from;
		long double Int = 0;
		if ((En_from-En_y)/En_from<1e-6) {//irregularity case
			E = 1e-6*En_from + En_y;
			Int+=0.5*(GetCrossSection(incident, E) + GetCrossSection(incident, En_from))*sqrt(E_prev*(E_prev-En_y));
			En_from = E;
			E_prev = En_from;
		}
		ColoredRange energy_range_;
		for (auto c = components_.begin(), c_end_ = components_.end(); c != c_end_; ++c)
			energy_range_ += c->first->GetEnergies();
		energy_range_.Trim(0, gSettings.ProgConsts()->maximal_energy);
		energy_range_.Trim(En_from, En_to);
		int i_end_ = energy_range_.NumOfIndices();
		for (int i=0; i!=i_end_;++i) {
			E = energy_range_.Value(i);
			Int+=GetCrossSection(incident, E)*sqrt(E/(E-En_y))*(E-E_prev);
			E_prev = E;
		}
		return Int;
	} else
		return GetXSIntegral(incident, En_from, En_to, En_y);
}

long double Mixture::GetUntabEnergyFromXSIntegral(const Particle *incident, long double En_y, long double integral_value) const
{
	return GetEnergyFromXSIntegral(incident, En_y, integral_value);
	/* TODO?:
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
	right = std::min((double)right, gSettings.ProgConsts()->maximal_energy);
	long double I_max = XS_integral(e_start, right, Eny, event);
	if (I_max < LnR) {
		I_max = XS_integral(e_start, gSettings.ProgConsts()->maximal_energy, Eny, event);
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
			f_new = XS_integral(e_start, e_finish, Eny, event) - LnR;
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
	*/
}

const Particle* Mixture::GenerateUntabScatteringParticle(const Particle *incident, double E, double Rand) const
{
	return GenerateScatteringParticle(incident, E, Rand);
}

