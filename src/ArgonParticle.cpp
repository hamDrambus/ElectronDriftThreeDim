#include "ArgonParticle.h"

EnergyScanner::EnergyScanner(ScanType type): i(0), type_(type)
{
	/*ElasticXS, Resonance_3o2_XS, Resonance_1o2_XS, ResonancesXS,
		Resonance_3o2_DiffXS, Resonance_1o2_DiffXS, ResonancesDiffXS,
		DiffXS,	InelasticXS, ElasticResXS, XSIntegral, PlotElastic,
		PlotResonance_3o2, PlotResonance_1o2, PlotResonances,
		PlotDiffXS, PlotInelastic, PlotElasticResXS, PlotAllXS*/
	switch (type_) {
	case (ElasticXS): { //used for table construction
		//from 1e-3 eV to 0.1 eV with step 5e-4 eV, etc.
		energy_range_ =ColoredInterval (0, gSettings.PhysConsts()->XS_el_En_minimum, 1e-4) +
				ColoredInterval (gSettings.PhysConsts()->XS_el_En_minimum, 0.1, 5e-4) + ColoredInterval (0.1, 1, 1e-3) +
				ColoredInterval (1, 10, 5e-2) + ColoredInterval (10, gSettings.PhysConsts()->XS_el_En_maximum, 0.1)+
				ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/2) + 	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/2) + 	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/30) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/30); 		//fine area
		energy_range_.Trim(0, gSettings.ProgConsts()->maximal_energy);
		break;
	}
	case (InelasticXS): {
		energy_range_ = ColoredInterval (11.5, gSettings.ProgConsts()->maximal_energy, 0.01);
		break;
	}
	case (ElasticResXS): {
		energy_range_ = ColoredInterval (gSettings.PhysConsts()->XS_el_En_minimum*0.1, 0.1, 5e-4) + ColoredInterval (0.1, 1, 1e-3) +
				ColoredInterval (1, 10, 5e-2) + ColoredInterval (10, gSettings.PhysConsts()->XS_el_En_maximum, 0.1) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/2) + 		//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/2) + 		//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/30) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/30);	//fine area
		energy_range_.Trim(0, gSettings.ProgConsts()->maximal_energy);
		break;
	}
	case (XSIntegral): {
		energy_range_ = ColoredInterval (0, 0.1, 2e-4) + ColoredInterval (0.1, 1, 8e-4) +
				ColoredInterval (1, 10, 1e-2) + ColoredInterval (10, gSettings.PhysConsts()->XS_el_En_maximum, 0.02) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80) +	//fine area
				ColoredInterval (11.5, gSettings.PhysConsts()->XS_el_En_maximum, 0.003);
		energy_range_.Trim(0, gSettings.ProgConsts()->maximal_energy);
		break;
	}
	case (PlotElastic): { //not used
		//from 1e-3 eV to 0.1 eV with step 5e-4 eV, etc.
		energy_range_ = ColoredInterval (gSettings.PhysConsts()->XS_el_En_minimum, 0.1, 3e-4) + ColoredInterval (0.1, 1, 9e-4);
		break;
	}
	case (PlotResonance_3o2): {
		energy_range_ = ColoredInterval (En_3o2_ - 110*Width_3o2_, En_3o2_ + 110*Width_3o2_, Width_3o2_/3) +//coarse area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80); 		//fine area
		energy_range_.Trim(0, gSettings.PhysConsts()->XS_el_En_maximum);
		break;
	}
	case (PlotResonance_1o2): {
		energy_range_ = ColoredInterval (En_1o2_ - 110*Width_1o2_, En_1o2_ + 110*Width_1o2_, Width_1o2_/3) +//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80); 		//fine area
		energy_range_.Trim(0, gSettings.PhysConsts()->XS_el_En_maximum);
		break;
	}
	case (PlotResonances): { //used for plotting total and differential XS in resonances region
		energy_range_ = ColoredInterval (En_1o2_ - 110*Width_1o2_, En_1o2_ + 110*Width_1o2_, Width_1o2_/3) + 	//coarse area
				ColoredInterval (En_3o2_ - 110*Width_3o2_, En_3o2_ + 110*Width_3o2_, Width_3o2_/3) + 			//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) + 		//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80); 			//fine area
		energy_range_.Trim(0, gSettings.PhysConsts()->XS_el_En_maximum);
		break;
	}
	case (PlotDiffXS): {  //used for plotting test of differential to total cross sections
		energy_range_ = ColoredInterval (1e-4, 0.1, 7e-4) +
				ColoredInterval (0.1, 1, 7e-3) +
				ColoredInterval (1, gSettings.PhysConsts()->XS_el_En_maximum, 0.086);
		energy_range_.Trim(0, gSettings.ProgConsts()->maximal_energy);
		break;
	}
	case (PlotInelasticXS): {
		energy_range_ = ColoredInterval (11.5, gSettings.PhysConsts()->XS_el_En_maximum, 0.007) +
				ColoredInterval (gSettings.PhysConsts()->XS_el_En_maximum, 100, 0.1);
		break;
	}
	case (PlotElasticResXS): {
		energy_range_ = ColoredInterval (gSettings.PhysConsts()->XS_el_En_minimum, 0.1, 3e-4) + ColoredInterval (0.1, 1, 9e-4) +
				ColoredInterval (1, 10, 7e-3) + ColoredInterval (10, gSettings.PhysConsts()->XS_el_En_maximum, 0.016) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) +//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80);//fine area
		energy_range_.Trim(0, gSettings.ProgConsts()->maximal_energy);
		break;
	}
	case (PlotAllXS): {
		energy_range_ = ColoredInterval (gSettings.PhysConsts()->XS_el_En_minimum*0.1, 0.1, 2e-4) + ColoredInterval (0.1, 1, 7e-4) +
				ColoredInterval (1, 10, 7e-3) + ColoredInterval (10, gSettings.PhysConsts()->XS_el_En_maximum, 0.016) +
				ColoredInterval (En_1o2_ - 110*Width_1o2_, En_1o2_ + 110*Width_1o2_, Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 110*Width_3o2_, En_3o2_ + 110*Width_3o2_, Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) +//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80) +//fine area
				ColoredInterval (11.5, gSettings.PhysConsts()->XS_el_En_maximum, 0.003);
		energy_range_.Trim(0, gSettings.ProgConsts()->maximal_energy);
		break;
	}
	}
}
long double EnergyScanner::Next(int& err)
{
	if (i<energy_range_.NumOfIndices()) {
		long double val = energy_range_.Value(i);
		++i;
		if (DBL_MAX==val){
			err = 1;
			Reset();
		} else
			err = 0;
		return val;
	}
	err = 1;
	Reset();
	return DBL_MAX;
}
void EnergyScanner::Reset(void)
{
	i = 0;
}


InelasticProcess::InelasticProcess(std::string name, unsigned int ID, double En, double F, std::vector<double> &Ens, std::vector<double> &XSs,
	ArExperimental *ArExper):
		name_(name), ID_(ID), En_threshold_(En), Oscillator_strength_(F), exp_XS_(Ens, XSs, 3, 4), ArExper_(ArExper)
{}

double InelasticProcess::operator ()(double E) const //returns cross section in 1e-16 cm^2
{
	if (E<En_threshold_)
		return 0;
	std::size_t sz = exp_XS_.size();
	if (0==sz)
		goto bb;
	if (E>exp_XS_.getX(sz-1))
		goto bb;
	return Exp_XS(E);
	bb:
	double out = BB_XS(E);
	return out>0 ? out : 0;
}

double InelasticProcess::BB_XS(double E) const
{
	double gamma = (E+gSettings.PhysConsts()->e_mass_eV)/ gSettings.PhysConsts()->Ar_mass_eV;
	double gamma2 = gamma*gamma;
	double beta2 = 1.0-1.0/gamma2;
	double BBconst = 8.0*M_PI* gSettings.PhysConsts()->a_bohr_SI *gSettings.PhysConsts()->a_bohr_SI
		*gSettings.PhysConsts()->Ry_eV*gSettings.PhysConsts()->Ry_eV / gSettings.PhysConsts()->e_mass_eV;
	return (Oscillator_strength_/(En_threshold_*beta2))*log(beta2*gamma2*gSettings.PhysConsts()->e_mass_eV /(2.0*En_threshold_)-beta2)
		*BBconst*E/(E+En_threshold_+ArExper_->E_Ionization);
}

double InelasticProcess::Exp_XS(double E) const
{
	if (E<En_threshold_)
		return 0;
	double out = exp_XS_(E, E);
	return out>0 ? out : 0;
}

double InelasticProcess::get_En_thresh (void) const
{	return En_threshold_;}
std::string InelasticProcess::get_name (void) const
{	return name_;}
unsigned int InelasticProcess::get_ID (void) const
{	return ID_;}

ArExperimental::ArExperimental(void): total_elastic_cross(3, 5) /*fit by 3rd order polynomial*/, max_process_ID(0), E_Ionization(0), is_valid_(false)
{
	std::cout << "Reading Ar exprimental data..." << std::endl;
	std::ifstream inp;
	std::string line, word;
	unsigned int MAX_L = 6; //L == (MAX_L-1) is the last one //TODO: move this "magic" number to either setings or phases' input data file.
	inp.open(gSettings.ProgConsts()->elastic_XS_fname);
	if (!inp.is_open()) {
		std::cerr << "ArExperimental::ArExperimental: Failed to open file \"" << gSettings.ProgConsts()->elastic_XS_fname << "\"" << std::endl;
		goto fail_load;
	}
	while (!inp.eof()) {
		std::getline(inp, line);
		if (line.size()>=2)
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(std::stod(word));
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double XS = std::stod(word); //Cross section is stored in 1e-20 m^2 in both files and this program.
		total_elastic_cross.insert(k, XS);
	}
	inp.close();
	for (int l = 0; l < MAX_L; ++l) {
		if (l < 4) {
			phase_shifts_pos_.push_back(DataVector(3, 4)); /*interpolation by 3rd order polynomial because there's data for 11eV*/
			phase_shifts_pos_.back().use_rightmost(true);
			phase_shifts_neg_.push_back(DataVector(3, 4));
			phase_shifts_neg_.back().use_rightmost(true);
		} else {
			phase_shifts_pos_.push_back(DataVector(3, 5)); /*fit by 3rd order polynomial TODO: tests, tests*/
			phase_shifts_pos_.back().use_rightmost(true);
			phase_shifts_neg_.push_back(DataVector(3, 4));
			phase_shifts_neg_.back().use_rightmost(true);
		}
	}
	inp.open(gSettings.ProgConsts()->elastic_XS_phaseshift_fname);
	if (!inp.is_open()) {
		std::cerr << "ArExperimental::ArExperimental: Failed to open file \"" << gSettings.ProgConsts()->elastic_XS_phaseshift_fname << "\"" << std::endl;
		goto fail_load;
	}
	while (!inp.eof()) {
		std::getline(inp, line);
		bool positive = true;
		if (line.size()>=2) {
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
		}
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double k = std::stod(word);
		//now process string "double\tdouble\tdouble\t..."
		//\t\t means that the double is omitted in the table
		std::vector<double> vals;
		std::vector<bool> is_d;//is double presented in the table
		vals.resize(MAX_L, 0);
		is_d.resize(MAX_L, false);
		for (int l=0;l<MAX_L;++l) {
			if (line.empty())
				break;
			if (line[0]=='\t') { //value is omitted in the table
				line.erase(line.begin());
				continue;
			} //now line does not start with \t
			word = strtoken(line,"\t"); //removes one \t after value
			vals[l] = std::stod(word);
			is_d[l] = true;
		}
		line.clear();
		std::getline(inp, line);
		if (line.size()>=2) {
			if ((line[0]!='\t')||(line[1]!='\t')) //Not a line with negative phase shifts.
				continue;
			line.erase(line.begin(), line.begin()+5);//remove four '\t' for k and L=0 (\t'k'\t'L=0'\t but k and L are substituted for \t for this line
		}
		std::vector<double> vals_1;
		std::vector<bool> is_d_1;
		vals_1.resize(MAX_L-1,0);
		is_d_1.resize(MAX_L-1, false);
		for (int l=0; l<(MAX_L-1); ++l) {
			if (line.empty())
				break;
			if (line[0]=='\t') { //value is omitted in the table
				line.erase(line.begin());
				continue;
			}
			//now line does not start with \t
			word = strtoken(line,"\t"); //removes one \t after value
			vals_1[l] = std::stod(word);
			is_d_1[l] = true;
		}
		for (int l=0; l<MAX_L; ++l) {
			if (is_d[l]) {
				if (0==l)
					phase_shifts_neg_[l].insert(k, vals[l]);
				phase_shifts_pos_[l].insert(k, vals[l]);
			}
		}
		for (int l=0; l<(MAX_L-1); ++l) {
			if (is_d_1[l])
				phase_shifts_neg_[l+1].insert(k, vals_1[l]);
		}
	}
	inp.close();
	//loading ionization and excitation data
	inp.open(gSettings.ProgConsts()->ionization_XS_fname);
	if (!inp.is_open()) {
		std::cerr << "ArExperimental::ArExperimental: Failed to open file \"" << gSettings.ProgConsts()->ionization_XS_fname << "\"" << std::endl;
		goto fail_load;
	}
	read_inelastic(inp, ionizations);
	inp.close();
	if (ionizations.size()) {
		E_Ionization = DBL_MAX;
		for (int i =0, end_ = ionizations.size();i!=end_;++i)
			E_Ionization = std::min(ionizations[i].get_En_thresh(), E_Ionization);
	}
	inp.open(gSettings.ProgConsts()->excitation_XS_fname);
	if (!inp.is_open()) {
		std::cerr << "ArExperimental::ArExperimental: Failed to open file \"" << gSettings.ProgConsts()->excitation_XS_fname << "\"" << std::endl;
		goto fail_load;
	}
	read_inelastic(inp, excitations);
	inp.close();
	if (excitations.size()) {
		First_excitation_En = DBL_MAX;
		for (int i = 0, end_ = excitations.size(); i != end_; ++i)
			First_excitation_En = std::min(excitations[i].get_En_thresh(), First_excitation_En);
	}
	is_valid_ = true;
	std::cout << "Finished reading Ar exprimental data." << std::endl;
	return;
fail_load: //TODO:implement through try-catch (create bad_data or something exception)
	std::cerr << "Failed to load Ar experimental data" << std::endl;
	total_elastic_cross.clear();
	phase_shifts_neg_.clear();
	phase_shifts_pos_.clear();
	ionizations.clear();
	excitations.clear();
	is_valid_ = false;
	return;
}

bool ArExperimental::isValid(void) const
{	return is_valid_; }

void ArExperimental::read_inelastic(std::ifstream &inp, std::vector<InelasticProcess> &to)
{
	short read_section = -1; //-1 - none, 0 - header, 1 - first data line, 2 - second data line
	std::string proc_name;
	double En_thresh, F=0, XS_units;
	std::vector<double> E_vals, XS_vals;
	int line_count = 0;
	std::string line, word;
	while (!inp.eof()) {
		std::getline(inp, line);
		++line_count;
		if (line.size()>=2) {
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
		}
		word = strtoken(line, "\t ,");
		if (word.empty())
			continue;
		if (word[0]=='"') {
			if (read_section>=0) {//write previous data and header
				if (E_vals.size()!=XS_vals.size()) {
					std::size_t new_sz = std::min(E_vals.size(),XS_vals.size());
					std::cout<<"Line "<<line_count<<" warning! "<<proc_name<<" values size mismatch. Trimming data to size: "<<new_sz<<std::endl;
					if (E_vals.size()>new_sz)
						E_vals.erase(E_vals.begin()+new_sz, E_vals.end());
					if (XS_vals.size()>new_sz)
						XS_vals.erase(XS_vals.begin()+new_sz, XS_vals.end());
				}
				to.push_back(InelasticProcess(proc_name, max_process_ID, En_thresh, F, E_vals, XS_vals, this));
				F=0;
				E_vals.clear();
				XS_vals.clear();
				++max_process_ID;
				read_section = -1;
			}
			proc_name = word;
			word = strtoken(line, "\t ,");
			if (word.empty()) {
				std::cout<<"Line "<<line_count<<" warning! No energy specified. Line is ignored"<<std::endl;
				continue;
			}
			En_thresh = std::stod(word);
			word = strtoken(line, "\t ,");
			if (word.empty()) {
				std::cout<<"Line "<<line_count<<" warning! No cross section units specified. Line is ignored"<<std::endl;
				continue;
			}
			XS_units = std::stod(word);
			XS_units/=1e-16; //All cross sections in the program are in 1e-16 cm^2
			word = strtoken(line, "\t ,");
			if (!word.empty())
				F = std::stod(word);
			read_section = 0;
		} else {
			if ((-1==read_section)||(2==read_section)) {
				std::cout<<"Line "<<line_count<<" warning! line ignored"<<std::endl;
				continue;
			}
			while (!word.empty()) {
				double val = std::stod(word);
				if (read_section==0)
					E_vals.push_back(val);
				else
					XS_vals.push_back(val*XS_units);
				word = strtoken(line, "\t ,");
			}
			++read_section;
		}
	}
	if (read_section>=0) {//write previous data and header
		if (E_vals.size()!=XS_vals.size()){
			std::size_t new_sz = std::min(E_vals.size(),XS_vals.size());
			std::cout<<"Line "<<line_count<<" warning! "<<proc_name<<" values size mismatch. Trimming data to size: "<<new_sz<<std::endl;
			if (E_vals.size()>new_sz)
				E_vals.erase(E_vals.begin()+new_sz, E_vals.end());
			if (XS_vals.size()>new_sz)
				XS_vals.erase(XS_vals.begin()+new_sz, XS_vals.end());
		}
		to.push_back(InelasticProcess(proc_name, max_process_ID, En_thresh, F, E_vals, XS_vals, this));
		++max_process_ID;
	}
}

const InelasticProcess * ArExperimental::FindInelastic(short ID) const
{
	for (int i =0, end_ = ionizations.size(); i!=end_; ++i) {
		if (ID == ionizations[i].get_ID()) {
			return &(ionizations[i]);
		}
	}
	for (int i =0, end_ = excitations.size(); i!=end_; ++i) {
		if (ID == excitations[i].get_ID()) {
			return &(excitations[i]);
		}
	}
	return NULL;
}

unsigned int ArExperimental::max_L (long double k) const
{
	return std::max(phase_shifts_pos_.size(), phase_shifts_neg_.size())-1;
}

void ArExperimental::phase_shift (long double k, unsigned int l, long double &ps_pos,  long double &ps_neg) const
{
	if (l>max_L(k)) {
		ps_pos = 0;
		ps_neg = 0;
		return;
	}
	if (l<phase_shifts_pos_.size()) {
		ps_pos = phase_shifts_pos_[l](k, k);
		if (l<phase_shifts_neg_.size())
			ps_neg = phase_shifts_neg_[l](k, k);
		else
			ps_neg = ps_pos;
	} else {
		ps_neg = phase_shifts_neg_[l](k, k);
		ps_pos = ps_neg;
	}
}

ArAllData::ArAllData(void)
{}

//k is in atomic units
void ArAllData::argon_phase_values_exp(long double k, unsigned int l, long double &ps_p, long double &ps_n) const
{
	long double angle_pos, angle_neg;
	ArExper_.phase_shift(k, l, angle_pos, angle_neg);
	ps_p = angle_pos;
	ps_n = angle_neg;
	if (1 == l) {
		long double E = std::pow(k /gSettings.PhysConsts()->a_h_bar_2eM_e_SI, 2.0);//recalculation from atomic units to energy
		long double cot3o2 = -2 * (E - gSettings.PhysConsts()->En_3o2) / gSettings.PhysConsts()->Width_3o2;
		long double cot1o2 = -2 * (E - gSettings.PhysConsts()->En_1o2) / gSettings.PhysConsts()->Width_1o2;
		long double cos3o2 = cot3o2 / std::sqrt(1 + cot3o2*cot3o2);
		long double cos1o2 = cot1o2 / std::sqrt(1 + cot1o2*cot1o2);
		ps_p += std::acos(cos3o2);
		ps_n += std::acos(cos1o2);
	}
}

//k is in atomic units
void ArAllData::argon_phase_values_MERT5(long double k, unsigned int l, long double &ps_p, long double &ps_n) const
{
	//see Kurokawa Phys. Rev. A84 2011, MERT5+ fit http://dx.doi.org/10.1103/PhysRevA.84.062717
	double A = gSettings.PhysConsts()->MERT5_A;
	double D = gSettings.PhysConsts()->MERT5_D;
	double F = gSettings.PhysConsts()->MERT5_F;
	double G = gSettings.PhysConsts()->MERT5_G;
	double A1 = gSettings.PhysConsts()->MERT5_A1;
	double H = gSettings.PhysConsts()->MERT5_H;
	double alpha_d = gSettings.PhysConsts()->MERT5_alpha_d;
	double alpha_q = gSettings.PhysConsts()->MERT5_alpha_q;
	long double tan = 0;
	if (0 == l) {
		tan = -A*k*(1 + 4 * alpha_d*k*k*log(k) / 3) - M_PI*alpha_d*k*k / 3 + D*pow(k, 3) + F* pow(k, 4);
		tan /= (1 + G*pow(k, 3));
	}
	else {
		unsigned int l2 = 2 * l;
		long double al = M_PI / ((l2 + 3)*(l2 + 1)*(l2 - 1));
		long double bl = M_PI*(15 * pow(l2 + 1, 4) - 140 * pow(l2 + 1, 2) + 128) / (pow((l2 + 3)*(l2 + 1)*(l2 - 1), 3)*(l2 + 5)*(l2 - 3));
		long double cl = 3 * al / ((l2 + 5)*(l2 - 3));
		tan = al*alpha_d*k*k + (bl*pow(alpha_d, 2) + cl*alpha_q)*pow(k, 4);
		if (1 == l)
			tan += H*pow(k, 5) - A1*pow(k, 3);
	}
	ps_p = ps_n = std::atan(tan);
}

//In seconds. Considered non zero only for L = 1 and J = 1/2 and 3/2 (J and L are expressed in halves)
long double ArAllData::argon_time_delay_j(long double E, int J, int L) const
{
	if (L!=2) //Only Feshbach resonance is taken into the account
		return 0;
	if (((J != L - 1) && (J != L + 1)) || J < 0) {
		std::cerr << "ArAllData::argon_time_delay_j: Error: No momentum J=" << J << "/2 for e(s=1/2)-Ar(s=0) scattering with L =" << L / 2 << std::endl;
		return 0;
	}
	if (J == L-1) {
		return 2*gSettings.PhysConsts()->h_bar_eVs/
			(1 + std::pow(2 * (E - gSettings.PhysConsts()->En_1o2) / gSettings.PhysConsts()->Width_1o2, 2))/ gSettings.PhysConsts()->Width_1o2;
	} else {
		return 2 * gSettings.PhysConsts()->h_bar_eVs /
			(1 + std::pow(2 * (E - gSettings.PhysConsts()->En_3o2) / gSettings.PhysConsts()->Width_3o2, 2)) / gSettings.PhysConsts()->Width_3o2;
	}
}

//Probability on scatter at certain angle going through total momentum J/2 and orbital momentum L/2. (J&L are expressed in halves)
//Can be negative.
//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
//E in eV
long double ArAllData::argon_scatter_probability_j(long double E, long double theta, int J, int L, int mode) const
{
	if (L % 2 != 0 || L<0) {
		std::cerr << "ArAllData::argon_scatter_probability_j: Error: No orbital momentum L=" << L << "/2 for e(s=1/2)-Ar(s=0) scattering"<< std::endl;
		return 0;
	}
	if (((J != L+1)&&(J != L - 1)) ||J<0) {
		std::cerr << "ArAllData::argon_scatter_probability_j: Error: No momentum J=" << J << "/2 for e(s=1/2)-Ar(s=0) scattering with L ="<<L/2<< std::endl;
		return 0;
	}
	int L0 = L / 2;
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (gSettings.PhysConsts()->phases_En_minimum > E)
			E = gSettings.PhysConsts()->phases_En_minimum;
		k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
		if (E<gSettings.PhysConsts()->phases_En_threshold) {
			L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
			phase_values = &ArAllData::argon_phase_values_MERT5;
		}
		else {
			L_MAX = ArExper_.max_L(k);
			phase_values = &ArAllData::argon_phase_values_exp;
		}
		break;
	}
	}
	long double prob = 0;
	LegendrePolynom P1, P2;
	AssociatedLegendrePolynom AP1, AP2;
	long double cos_th = cos(theta);
	if (L0 > L_MAX)
		return prob;
	long double ph_l0_p = 0; //p - positive, J = L+1/2
	long double ph_l0_n = 0; //n - negative, J = L-1/2
	((*this).*phase_values)(k, L0, ph_l0_p, ph_l0_n);
	ph_l0_p *= 2;
	ph_l0_n *= 2;
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double ph_l_p = ph_l0_p;
		long double ph_l_n = ph_l0_n;
		if (l != L0) {
			((*this).*phase_values)(k, l, ph_l_p, ph_l_n);
			ph_l_p *= 2;
			ph_l_n *= 2;
		}
		long double F_l0_l = 0;
		long double G_l0_l = 0;
		if (J == L + 1) {
			F_l0_l = (L0 + 1)*((l + 1)*cos(ph_l0_p-ph_l_p) - (2*l+1)*cos(ph_l0_p) - (l+1)*cos(ph_l_p) + l*cos(ph_l0_p-ph_l_n) - l*cos(ph_l_n) + 2*l+1);
			G_l0_l = cos(ph_l0_p-ph_l_p) - cos(ph_l0_p-ph_l_n) - cos(ph_l_p) + cos(ph_l_n);
		} else {
			F_l0_l = (L0)*((l + 1)*cos(ph_l0_n-ph_l_p) - (2*l+1)*cos(ph_l0_n) - (l+1)*cos(ph_l_p) + l*cos(ph_l0_n-ph_l_n) - l*cos(ph_l_n) + 2*l+1);
			G_l0_l = -1*(cos(ph_l0_n-ph_l_p) - cos(ph_l0_n-ph_l_n) - cos(ph_l_p) + cos(ph_l_n));
		}
		prob += F_l0_l*P1(cos_th, l)*P2(cos_th, L0) + G_l0_l*AP1(cos_th, l, 1)*AP2(cos_th, L0, 1);
	}
	//cross = prob * M_PI / (2.0*pow(k, 2)) *a_bohr_SIconst*a_bohr_SIconst;
	return prob;
}

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
//mode = 3 - standard formula, but using argon_scatter_probability_j function (slower version, for testing argon_scatter_probability_j)
//E in eV
long double ArAllData::argon_cross_elastic_diff(long double E, long double theta, int mode) const {
	//different formulas are used for E<0.24eV and E>0.24eV!
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units
	if (1!=mode&&2!=mode)
	{
		if (gSettings.PhysConsts()->phases_En_minimum>E)
			k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(gSettings.PhysConsts()->phases_En_minimum);
	}
	long double cross = argon_scatter_spin_flip_amplitude_sq(E, theta, mode);
	cross += argon_scatter_spin_nonflip_amplitude_sq(E, theta, mode);
	cross *= M_PI / (2.0*pow(k, 2));
	cross = std::max(cross, (long double)0);
	return cross*gSettings.PhysConsts()->a_bohr_SI*gSettings.PhysConsts()->a_bohr_SI; //const is multiplied by 1e10
}

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
long double ArAllData::argon_scatter_spin_flip_amplitude_sq(long double E, long double theta, int mode) const
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (gSettings.PhysConsts()->phases_En_minimum>E)
			E = gSettings.PhysConsts()->phases_En_minimum;
		k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
		if (E< gSettings.PhysConsts()->phases_En_threshold) {
			L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
			phase_values = &ArAllData::argon_phase_values_MERT5;
		} else {
			L_MAX = ArExper_.max_L(k);
			phase_values = &ArAllData::argon_phase_values_exp;
		}
		break;
	}
	}
	long double cross = 0;
	AssociatedLegendrePolynom AP1, AP2;
	long double cos_th = cos(theta);
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double l_p = 0; //p - positive
		long double l_n = 0; //n - negative
		((*this).*phase_values)(k, l, l_p, l_n);
		l_p *= 2; //all cosines are taken from double angles
		l_n *= 2;
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double f_p = l_p;
			long double f_n = l_n;
			if (l != f) {
				((*this).*phase_values)(k, f, f_p, f_n);
				f_p *= 2;
				f_n *= 2;
			}
			long double G_l_f = ((0 == l) ? 0 :
				cos(l_p - f_p) - cos(l_n - f_p) - cos(l_p - f_n) + cos(l_n - f_n));
			if (l != f) { //Nondiagonal sum is reduced because f starts from l instead of 0.
				G_l_f *= 2.0;
			}
			cross += G_l_f*AP1(cos_th, l, 1)*AP2(cos_th, f, 1);
		}
	}
	cross = std::max(cross, (long double)0); //can be less than 0 only because of computational uncertainties
	return cross; //|b|^2
}

long double ArAllData::argon_scatter_spin_nonflip_amplitude_sq(long double E, long double theta, int mode) const
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (gSettings.PhysConsts()->phases_En_minimum>E)
			E = gSettings.PhysConsts()->phases_En_minimum;
		k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
		if (E<gSettings.PhysConsts()->phases_En_threshold) {
			L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
			phase_values = &ArAllData::argon_phase_values_MERT5;
		} else {
			L_MAX = ArExper_.max_L(k);
			phase_values = &ArAllData::argon_phase_values_exp;
		}
		break;
	}
	}
	long double cross = 0;
	LegendrePolynom P1, P2;
	long double cos_th = cos(theta);
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double l_p = 0; //p - positive
		long double l_n = 0; //n - negative
		((*this).*phase_values)(k, l, l_p, l_n);
		l_p *= 2; //all cosines are taken from double angles
		l_n *= 2;
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double f_p = l_p;
			long double f_n = l_n;
			if (l != f) {
				((*this).*phase_values)(k, f, f_p, f_n);
				f_p *= 2;
				f_n *= 2;
			}
			long double F_l_f = (l + 1)*(f + 1)*cos(l_p - f_p) + (l + 1)*(f)*cos(l_p - f_n) + (l)*(f + 1)*cos(l_n - f_p) +
				l*f*cos(l_n - f_n) - (l + 1)*(2 * f + 1)*cos(l_p) - (2 * l + 1)*(f + 1)*cos(f_p) +
				(2 * l + 1)*(2 * f + 1) - (l)*(2 * f + 1)*cos(l_n) - (2 * l + 1)*(f)*cos(f_n);
			if (l != f) { //Nondiagonal sum is reduced because f starts from l instead of 0.
				F_l_f *= 2.0;
			}
			cross += F_l_f*P1(cos_th, l)*P2(cos_th, f);
		}
	}
	cross = std::max(cross, (long double)0); //can be less than 0 only because of computational uncertainties
	return cross; //|a|^2
}

//mode = 0 or unexpected - standard formula used in simulation, called by default. Using smoothing between MERT5 and experimental XS
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_XS_)
//mode = 2 - calculate using extrapolation of experimental total XS (normally used only between THRESH_E_XS_ and EN_MAXIMUM_) + resonances
//mode = 3 - calculate using experimental phases (for testing only) + resonances
//mode = 4 - same as 0, but no smoothing
//below XS_EL_EN_MINIMUM_ linear extrapolation to XS(0)=7.491 (1e-20 m^2) (in default mode)
//TODO: add no Ramseur minimum case
long double ArAllData::argon_cross_elastic(long double E, int mode) const//Tabulation of this function must be done carefully. Do not forget near 0 case.
{
	long double cross = 0;
	switch (mode) {
	case 1: {
		long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units is following:
														// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
		unsigned int L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		for (unsigned int l = 0; l <= L_MAX; ++l) {
			long double phase_l_p = 0;
			long double phase_l_n = 0;
			argon_phase_values_MERT5(k, l, phase_l_p, phase_l_n);
			cross += (l + 1)*sin(phase_l_p)*sin(phase_l_p) + l*sin(phase_l_n)*sin(phase_l_n);
		}
		cross *= 4.0*M_PI / pow(k, 2);
		cross *= gSettings.PhysConsts()->a_bohr_SI*gSettings.PhysConsts()->a_bohr_SI;
		break;
	}
	case 2: {
		//this is extrapolation of experimental data
		cross = ArExper_.total_elastic_cross(gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E), gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E));
		//this part is calculated using phase shifts
		cross += argon_cross_resonance_1o2(E);
		cross += argon_cross_resonance_3o2(E);
		break;
	}
	case 3: {
		long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units
		unsigned int L_MAX = ArExper_.max_L(k);;
		for (unsigned int l = 0; l <= L_MAX; ++l) {
			long double phase_l_p = 0;
			long double phase_l_n = 0;
			argon_phase_values_exp(k, l, phase_l_p, phase_l_n); //includes Feshbach resonances
			cross += (l + 1)*sin(phase_l_p)*sin(phase_l_p) + l*sin(phase_l_n)*sin(phase_l_n);
		}
		cross *= 4.0*M_PI / pow(k, 2);
		cross *= gSettings.PhysConsts()->a_bohr_SI*gSettings.PhysConsts()->a_bohr_SI;
		break;
	}
	case 4: {
		if (gSettings.PhysConsts()->no_ramsauer_minimum) {
			double E1 = gSettings.PhysConsts()->ramsauer_minimum_En;
			double XS1 = argon_cross_elastic(E1, 1); //MERT5
			if (E<E1)
				return E*XS1 / E1; //linear increase from (0,0). No Ramsauer minimum.
		}
		if (E<gSettings.PhysConsts()->XS_el_En_minimum) {
			double XS0 = gSettings.PhysConsts()->XS_el_at_0_En;
			cross = XS0 + E*(argon_cross_elastic(gSettings.PhysConsts()->XS_el_En_minimum, 1) - XS0) / gSettings.PhysConsts()->XS_el_En_minimum;
			//linear behavior. MERT5 is used at XS_EL_EN_MINIMUM_
		}
		else {
			if (E < gSettings.PhysConsts()->XS_el_En_thresold) {
				cross = argon_cross_elastic(E, 1); //MERT5
			}
			else {
				cross = argon_cross_elastic(E, 2); //Experiment
			}
		}
		break;
	}
	default: {
		if (gSettings.PhysConsts()->no_ramsauer_minimum) {
			double E1 = gSettings.PhysConsts()->ramsauer_minimum_En;
			double XS1 = argon_cross_elastic(E1, 1); //MERT5
			if (E<E1)
				return E*XS1 / E1; //linear increase from (0,0). No Ramsauer minimum.
		}
		if (E<gSettings.PhysConsts()->XS_el_En_minimum) {
			double XS0 = gSettings.PhysConsts()->XS_el_at_0_En;
			cross = XS0 + E*(argon_cross_elastic(gSettings.PhysConsts()->XS_el_En_minimum, 1) - XS0) / gSettings.PhysConsts()->XS_el_En_minimum;
			//linear behavior. MERT5 is used at XS_el_En_minimum
		} else {
			double smooth_d_En = std::fabs(gSettings.PhysConsts()->XS_el_En_smooth_width);
			double E_l = gSettings.PhysConsts()->XS_el_En_thresold - smooth_d_En;
			double E_r = gSettings.PhysConsts()->XS_el_En_thresold + smooth_d_En;
			if (E <=E_l) {
				cross = argon_cross_elastic(E, 1);//MERT5
				break;
			}
			if (E >E_r) {
				cross = argon_cross_elastic(E, 2);//Experiment
				break;
			}
			cross = ((E_r - E)*argon_cross_elastic(E, 1) + (E - E_l)*argon_cross_elastic(E, 2)) / (E_r - E_l);
			//linear smoothing between MERT5 and experiment
		}
		break;
	}
	}
	return cross;
}

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
long double ArAllData::argon_delay_spin_flip (long double E, long double theta, int mode) const
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (gSettings.PhysConsts()->phases_En_minimum>E)
			E = gSettings.PhysConsts()->phases_En_minimum;
		k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
		if (E<gSettings.PhysConsts()->phases_En_threshold) {
			L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
			phase_values = &ArAllData::argon_phase_values_MERT5;
		} else {
			L_MAX = ArExper_.max_L(k);
			phase_values = &ArAllData::argon_phase_values_exp;
		}
		break;
	}
	}
	long double cos_sum = 0;
	long double sin_sum = 0;
	long double cos_delayed_sum =0;
	long double sin_delayed_sum =0;
	AssociatedLegendrePolynom AP1;
	long double cos_th = cos(theta);
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double l_p = 0; //p - positive
		long double l_n = 0; //n - negative
		((*this).*phase_values)(k, l, l_p, l_n);
		l_p *= 2; //all cosines and sines are taken from double angles
		l_n *= 2;
		cos_sum +=AP1(cos_th, l, 1)*cos(l_p);
		cos_sum +=-1*AP1(cos_th, l, 1)*cos(l_n);
		sin_sum +=AP1(cos_th, l, 1)*sin(l_p);
		sin_sum +=-1*AP1(cos_th, l, 1)*sin(l_n);
		if (l!=0) {
			sin_delayed_sum+= 2*argon_time_delay_j(E, 2*l -1, 2*l) * -1*AP1(cos_th, l, 1) * sin(l_n);
			cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l -1, 2*l) * -1*AP1(cos_th, l, 1) *cos(l_n);
		}
		sin_delayed_sum+= 2*argon_time_delay_j(E, 2*l +1, 2*l) * AP1(cos_th, l, 1) *sin(l_p);
		cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l +1, 2*l) * AP1(cos_th, l, 1) *cos(l_p);
	}
	long double cos_2 = cos_sum*cos_sum;
	long double sin_2 = sin_sum*sin_sum;
	return cos_2/(cos_2 + sin_2) * (sin_delayed_sum * sin_sum/cos_2 + cos_delayed_sum /cos_sum);
}

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
long double ArAllData::argon_delay_spin_nonflip (long double E, long double theta, int mode) const
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (gSettings.PhysConsts()->phases_En_minimum>E)
			E = gSettings.PhysConsts()->phases_En_minimum;
		k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
		if (E< gSettings.PhysConsts()->phases_En_threshold) {
			L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
			phase_values = &ArAllData::argon_phase_values_MERT5;
		} else {
			L_MAX = ArExper_.max_L(k);
			phase_values = &ArAllData::argon_phase_values_exp;
		}
		break;
	}
	}
	long double cos_sum = 0;
	long double sin_sum = 0;
	long double cos_delayed_sum =0;
	long double sin_delayed_sum =0;
	LegendrePolynom P1;
	long double cos_th = cos(theta);
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double l_p = 0; //p - positive
		long double l_n = 0; //n - negative
		((*this).*phase_values)(k, l, l_p, l_n);
		l_p *= 2; //all cosines are taken from double angles
		l_n *= 2;
		cos_sum +=(l+1)*P1(cos_th, l)*cos(l_p);
		cos_sum +=(l)*P1(cos_th, l)*cos(l_n);
		cos_sum -=(2*l+1)*P1(cos_th, l);
		sin_sum +=(l+1)*P1(cos_th, l)*sin(l_p);
		sin_sum +=(l)*P1(cos_th, l)*sin(l_n);
		if (l!=0) {
			sin_delayed_sum+= 2*argon_time_delay_j(E, 2*l -1, 2*l) * P1(cos_th, l) * (l)*sin(l_n);
			cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l -1, 2*l) * P1(cos_th, l) * (l)*cos(l_n);
		}
		sin_delayed_sum+= 2*argon_time_delay_j(E, 2*l +1, 2*l) * P1(cos_th, l) * (l+1)*sin(l_p);
		cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l +1, 2*l) * P1(cos_th, l) * (l+1)*cos(l_p);
	}
	long double cos_2 = cos_sum*cos_sum;
	long double sin_2 = sin_sum*sin_sum;
	return cos_2/(cos_2 + sin_2) * (sin_delayed_sum * sin_sum/cos_2 + cos_delayed_sum /cos_sum);
}

long double ArAllData::argon_delay_spin_nonflip_prob (long double E, long double theta, int mode) const
{
	long double B = argon_scatter_spin_flip_amplitude_sq(E, theta, mode);
	long double A = argon_scatter_spin_nonflip_amplitude_sq(E, theta, mode);
	return A/(A+B);
}

long double ArAllData::argon_ResNBrS_spectrum(long double W, long double E) const //Normalization constant is arbitrary. Not dependant on E
{
	return std::exp(-0.5*std::pow((W - (ArExper_.First_excitation_En - gSettings.PhysConsts()->En_1o2)) / gSettings.PhysConsts()->Width_1o2, 2))
		+ std::exp(-0.5*std::pow((W - (ArExper_.First_excitation_En - gSettings.PhysConsts()->En_3o2)) / gSettings.PhysConsts()->Width_3o2, 2));
}

long double ArAllData::argon_ResNBrS_XS(long double E) const //Normalization constant is taken from "global_definitions.h"
{
	double width_factor = std::pow(gSettings.PhysConsts()->Width_3o2 / 2.0, 2) /
		(std::pow(gSettings.PhysConsts()->Width_3o2 / 2.0, 2) + std::pow((E - gSettings.PhysConsts()->En_3o2) / 2.0, 2));
	width_factor += std::pow(gSettings.PhysConsts()->Width_1o2 / 2.0, 2) /
		(std::pow(gSettings.PhysConsts()->Width_1o2 / 2.0, 2) + std::pow((E - gSettings.PhysConsts()->En_1o2) / 2.0, 2));
	return width_factor*gSettings.PhysConsts()->resonance_NBrS_XS;
}

long double ArAllData::argon_back_scatter_prob(long double E) const
{
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	if (E< gSettings.PhysConsts()->phases_En_threshold) {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
	} else {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
	}
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross = 0;
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double phase_l_p = 0;
		long double phase_l_n = 0;
		((*this).*phase_values)(k, l, phase_l_p, phase_l_n);
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double phase_f_p = phase_l_p;
			long double phase_f_n = phase_l_n;
			if (l != f) {
				((*this).*phase_values)(k, f, phase_f_p, phase_f_n);
			}
			long double cos_l_f = cos(phase_l_p - phase_f_p);
			W += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl(l, f, -1, 0, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
		}
		cross += (2 * l + 1)*sin(phase_l_p)*sin(phase_l_p);
	}
	return W / (2 * cross);
}
//energy loss and input are in eV
long double ArAllData::argon_TM_forward(long double E) const
{
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units is following:
																	  // k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	if (E< gSettings.PhysConsts()->phases_En_threshold) {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
	} else {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
	}
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross = 0;
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double phase_l_p = 0;
		long double phase_l_n = 0;
		((*this).*phase_values)(k, l, phase_l_p, phase_l_n);
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double phase_f_p = phase_l_p;
			long double phase_f_n = phase_l_n;
			if (l != f) {
				((*this).*phase_values)(k, f, phase_f_p, phase_f_n);
			}
			long double cos_l_f = cos(phase_l_p - phase_f_p);
			W += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl_transf(l, f, 0, 1, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl(l, f, 0, 1, 1e-5);
		}
	}
	return W / cross;
}

long double ArAllData::argon_TM_backward(long double E) const
{
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) const = 0;
	unsigned int L_MAX = 0;
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E); //recalculation from energy to atomic units is following:
																	  // k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	if (E< gSettings.PhysConsts()->phases_En_threshold) {
		L_MAX = gSettings.PhysConsts()->MERT5_Lmax;
		phase_values = &ArAllData::argon_phase_values_MERT5;
	} else {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
	}
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross = 0;
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double phase_l_p = 0;
		long double phase_l_n = 0;
		((*this).*phase_values)(k, l, phase_l_p, phase_l_n);
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double phase_f_p = phase_l_p;
			long double phase_f_n = phase_l_n;
			if (l != f) {
				((*this).*phase_values)(k, f, phase_f_p, phase_f_n);
			}
			long double cos_l_f = cos(phase_l_p - phase_f_p);
			W += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl_transf(l, f, -1, 0, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl(l, f, -1, 0, 1e-5);
		}
	}
	return W / cross;
}

long double ArAllData::argon_cross_resonance_3o2(long double E) const
{
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
	long double ps_p_0, ps_n_0;
	ArExper_.phase_shift(k, 1, ps_p_0, ps_n_0);
	long double ps_p = ps_p_0, ps_n = ps_n_0;
	long double cot3o2 = -2 * (E - gSettings.PhysConsts()->En_3o2) / gSettings.PhysConsts()->Width_3o2;
	long double cot1o2 = -2 * (E - gSettings.PhysConsts()->En_1o2) / gSettings.PhysConsts()->Width_1o2;
	long double cos3o2 = cot3o2 / std::sqrt(1 + cot3o2*cot3o2);
	long double cos1o2 = cot1o2 / std::sqrt(1 + cot1o2*cot1o2);
	ps_p += std::acos(cos3o2);
	ps_n += std::acos(cos1o2);

	long double cross = 2 * (pow(sin(ps_p), 2) - pow(sin(ps_p_0), 2));
	cross *= 4.0*M_PI / pow(k, 2);
	return cross * gSettings.PhysConsts()->a_bohr_SI * gSettings.PhysConsts()->a_bohr_SI;
}

long double ArAllData::argon_cross_resonance_1o2(long double E) const
{
	long double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
	long double ps_p_0, ps_n_0;
	ArExper_.phase_shift(k, 1, ps_p_0, ps_n_0);
	long double ps_p = ps_p_0, ps_n = ps_n_0;
	long double cot3o2 = -2 * (E - gSettings.PhysConsts()->En_3o2) / gSettings.PhysConsts()->Width_3o2;
	long double cot1o2 = -2 * (E - gSettings.PhysConsts()->En_1o2) / gSettings.PhysConsts()->Width_1o2;
	long double cos3o2 = cot3o2 / std::sqrt(1 + cot3o2*cot3o2);
	long double cos1o2 = cot1o2 / std::sqrt(1 + cot1o2*cot1o2);
	ps_p += std::acos(cos3o2);
	ps_n += std::acos(cos1o2);

	long double cross = (pow(sin(ps_n), 2) - pow(sin(ps_n_0), 2));
	cross *= 4.0*M_PI / pow(k, 2);
	return cross * gSettings.PhysConsts()->a_bohr_SI * gSettings.PhysConsts()->a_bohr_SI;
}

ArgonParticle::ArgonParticle(void) : Particle(), ArAllData_(),
	total_cross_elastic_fname(gSettings.ProgConsts()->tabulated_data_folder + "total_cross_section_elastic.dat"),
	total_cross_fname(gSettings.ProgConsts()->tabulated_data_folder + "total_cross_section.dat"),
	theta_table_fname(gSettings.ProgConsts()->tabulated_data_folder + "theta_probabilities.dat"),
	time_delay_spin_nonflip_prob_fname(gSettings.ProgConsts()->tabulated_data_folder + "time_delay_spin_nonflip_probabilities.dat"),
	time_delay_spin_flip_fname(gSettings.ProgConsts()->tabulated_data_folder + "time_delay_spin_flip.dat"),
	time_delay_spin_nonflip_fname(gSettings.ProgConsts()->tabulated_data_folder + "time_delay_spin_nonflip.dat"),
	total_resonance_NBrS_spectrum_fname(gSettings.ProgConsts()->tabulated_data_folder + "ResNBrS_spectrum.dat")
{
	name_ = ARGON_NAME;
	mass_ = gSettings.PhysConsts()->Ar_mass_eV;
	width_ = 0;
	std::cout<<"Loading Ar data tables..."<<std::endl;

	if (!ArAllData_.ArExper_.isValid()) {
		std::cout << GetName() << "::ArgonParticle: Error: Invalid experimental data" << std::endl;
		is_valid_ = false;
		std::cout << "Aborting loading data tables" << std::endl;
		return;
	}
	XS_En_sweeper_ = ColoredInterval (0, 0.001, 1e-4) + ColoredInterval (0.001, 0.01, 3e-4) + ColoredInterval (0.01, 0.1, 5e-4) +
			ColoredInterval (0, gSettings.PhysConsts()->XS_el_En_maximum, 1e-3) +
			ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/5) +	//coarse area
			ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/5) +	//coarse area
			ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) + 	//fine area
			ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80);		//fine area
	XS_En_sweeper_.Trim(0, gSettings.ProgConsts()->maximal_energy);

	std::map<unsigned int, std::string> entry;
	entry[0u] = "\"Elastic scattering\"";
	entry[1u] = "\"Resonance NBrS\"";
	for (short ID = 0; ID != ArAllData_.ArExper_.max_process_ID; ++ID) {
		const InelasticProcess *p = ArAllData_.ArExper_.FindInelastic(ID);
		if (NULL!=p) {
			entry[2 + ID] = p->get_name();
		} else {
			entry[2 + ID] = std::string("\"?Unknown?\"");
		}
	}
	processes_[ELECTRON_NAME] = entry;

	is_valid_ = true;
	//shared between threads - read only
	theta_table_.reset(new FunctionTable);
	time_delay_spin_flip_table_.reset(new FunctionTable);
	time_delay_spin_nonflip_table_.reset(new FunctionTable);
	time_delay_spin_nonflip_prob_table_.reset(new FunctionTable);

	ensure_file(total_cross_elastic_fname);
	ensure_file(total_cross_fname);
	ensure_file(theta_table_fname);
	ensure_file(time_delay_spin_nonflip_prob_fname);
	ensure_file(time_delay_spin_flip_fname);
	ensure_file(time_delay_spin_nonflip_fname);
	ensure_file(total_resonance_NBrS_spectrum_fname);

	std::ifstream inp;
	std::ofstream str;
	inp.open(total_cross_elastic_fname); //Cross section is stored in 1e-20 m^2 in both files and this program
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << total_cross_elastic_fname << "\"" << std::endl;
		if (generate_total_cross_elastic_table())
			total_cross_elastic_.write(total_cross_elastic_fname, "E[eV]\tXS elastic [1e-20 m^2]");
		else
			is_valid_ = false;
	} else {
		total_cross_elastic_.read(inp);
		inp.close();
		if (total_cross_elastic_.size() < total_cross_elastic_.getNused()) {
			std::cout << "Failed to load \"" << total_cross_elastic_fname << "\"" << std::endl;
			if (generate_total_cross_elastic_table())
				total_cross_elastic_.write(total_cross_elastic_fname, "E[eV]\tXS elastic [1e-20 m^2]");
			else
				is_valid_ = false;
		}
	}

	inp.open(total_cross_fname); //Cross section is stored in 1e-20 m^2 in both files and this program
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << total_cross_fname << "\"" << std::endl;
		if (generate_total_cross_table())
			total_cross_.write(total_cross_fname, "E[eV]\tXS elastic [1e-20 m^2]");
		else
			is_valid_ = false;
	}
	else {
		total_cross_.read(inp);
		inp.close();
		if (total_cross_.size() < total_cross_.getNused()) {
			std::cout << "Failed to load \"" << total_cross_fname << "\"" << std::endl;
			if (generate_total_cross_table())
				total_cross_.write(total_cross_fname, "E[eV]\tXS elastic [1e-20 m^2]");
			else
				is_valid_ = false;
		}
	}

	inp.open(theta_table_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout<<"Failed to load \""<<theta_table_fname<<"\""<<std::endl;
		if (generate_theta_table())
			theta_table_->write(theta_table_fname);
		else
			is_valid_ = false;
	} else {
		theta_table_->read(inp);
		inp.close();
		if (theta_table_->is_empty()) {
			std::cout<<"Failed to load \""<<theta_table_fname<<"\""<<std::endl;
			if (generate_theta_table())
				theta_table_->write(theta_table_fname);
			else
				is_valid_ = false;
		}
	}

	inp.open(time_delay_spin_nonflip_prob_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << time_delay_spin_nonflip_prob_fname << "\"" << std::endl;
		if (generate_time_delay_spin_nonflip_prob_table())
			time_delay_spin_nonflip_prob_table_->write(time_delay_spin_nonflip_prob_fname);
		else
			is_valid_ = false;
	} else {
		time_delay_spin_nonflip_prob_table_->read(inp);
		inp.close();
		if (time_delay_spin_nonflip_prob_table_->is_empty()) {
			std::cout << "Failed to load \"" << time_delay_spin_nonflip_prob_fname << "\"" << std::endl;
			if (generate_time_delay_spin_nonflip_prob_table())
				time_delay_spin_nonflip_prob_table_->write(time_delay_spin_nonflip_prob_fname);
			else
				is_valid_ = false;
		}
	}

	inp.open(time_delay_spin_nonflip_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << time_delay_spin_nonflip_fname << "\"" << std::endl;
		if (generate_time_delay_spin_nonflip_table())
			time_delay_spin_nonflip_table_->write(time_delay_spin_nonflip_fname);
		else
			is_valid_ = false;
	} else {
		time_delay_spin_nonflip_table_->read(inp);
		inp.close();
		if (time_delay_spin_nonflip_table_->is_empty()) {
			std::cout << "Failed to load \"" << time_delay_spin_nonflip_fname << "\"" << std::endl;
			if (generate_time_delay_spin_nonflip_table())
				time_delay_spin_nonflip_table_->write(time_delay_spin_nonflip_fname);
			else
				is_valid_ = false;
		}
	}

	inp.open(time_delay_spin_flip_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << time_delay_spin_flip_fname << "\"" << std::endl;
		if (generate_time_delay_spin_flip_table())
			time_delay_spin_flip_table_->write(time_delay_spin_flip_fname);
		else
			is_valid_ = false;
	} else {
		time_delay_spin_flip_table_->read(inp);
		inp.close();
		if (time_delay_spin_flip_table_->is_empty()) {
			std::cout << "Failed to load \"" << time_delay_spin_flip_fname << "\"" << std::endl;
			if (generate_time_delay_spin_flip_table())
				time_delay_spin_flip_table_->write(time_delay_spin_flip_fname);
			else
				is_valid_ = false;
		}
	}

	inp.open(total_resonance_NBrS_spectrum_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << total_resonance_NBrS_spectrum_fname << "\"" << std::endl;
		if (generate_ResNBrS_spectrum_table())
			total_resonance_NBrS_spectrum_.write(total_resonance_NBrS_spectrum_fname, "W[eV]\tFprob");
		else
			is_valid_ = false;
	} else {
		total_resonance_NBrS_spectrum_.read(inp);
		inp.close();
		if (total_resonance_NBrS_spectrum_.size()==0) {
			std::cout << "Failed to load \"" << total_resonance_NBrS_spectrum_fname << "\"" << std::endl;
			if (generate_ResNBrS_spectrum_table())
				total_resonance_NBrS_spectrum_.write(total_resonance_NBrS_spectrum_fname, "W[eV]\tFprob");
			else
				is_valid_ = false;
		}
	}
	std::cout<<"Finished loading Ar data tables"<<std::endl;
}

ArgonParticle::~ArgonParticle() {};

bool ArgonParticle::isValid(void) const
{
	return is_valid_ && ArAllData_.ArExper_.isValid();
}

void ArgonParticle::read_data (std::ifstream &inp, DataVector &data, long double y_factor)
{
	std::string line, word;
	while (!inp.eof()&&inp.is_open()) {
		std::getline(inp, line);
		if (line.size()>=2)
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double E = std::stod(word);
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double XS = std::stod(word)*y_factor;
		data.push_back(E, XS);
	}
}

bool ArgonParticle::generate_total_cross_elastic_table(void)
{
	std::cout << GetName() << ": Calculating total elastic cross section..." << std::endl;
	total_cross_elastic_.clear();
	total_cross_elastic_.setOrder(1); //interpolation with 1st order polynomial using 2 points
	total_cross_elastic_.setNused(2);
	if (!ArAllData_.ArExper_.isValid()) {
		std::cout << GetName() << ": Error: Invalid experimental data" << std::endl;
		return false;
	}
	double E, cross;
	for (long i = 0, i_end_ = XS_En_sweeper_.NumOfIndices(); i != i_end_; ++i) {
		E = XS_En_sweeper_.Value(i);
		cross = ArAllData_.argon_cross_elastic(E);
		total_cross_elastic_.push_back(E, cross);
	}
	return true;
}

bool ArgonParticle::generate_total_cross_table(void)
{
	std::cout << GetName() << ": Calculating total cross section..." << std::endl;
	total_cross_.clear();
	total_cross_.setOrder(1); //interpolation with 1st order polynomial using 2 points
	total_cross_.setNused(2);
	if (!ArAllData_.ArExper_.isValid()) {
		std::cerr << GetName() << ": Error: Invalid experimental data" << std::endl;
		return false;
	}
	bool use_table = true;
	bool skip_warning = false;
	if (total_cross_elastic_.size() < (total_cross_elastic_.getOrder() + 1)) {
		std::cerr << GetName() << ": Warning: Invalid elastic cross section table" << std::endl;
		use_table = false;
	}
	auto procs = processes_.find(ELECTRON_NAME);
	if (processes_.end() == procs) {
		std::cerr << GetName() << "::generate_total_cross_table: Error: unsupported particle \"" << ELECTRON_NAME << "\"" << std::endl;
		return false;
	}
	for (long i = 0, i_end_ = XS_En_sweeper_.NumOfIndices(); i != i_end_; ++i) {
		double E = XS_En_sweeper_.Value(i);
		double cross = 0;
		for (auto proc = procs->second.begin(), proc_end_ = procs->second.end(); proc != proc_end_; ++proc) {
			if (0 == proc->first) {
				if (use_table)
					cross += total_cross_elastic_(E, E);
				else 
					cross += ArAllData_.argon_cross_elastic(E);
				continue;
			}
			if (1 == proc->first) {
				cross += ArAllData_.argon_ResNBrS_XS(E);
				continue;
			}
			short ID = proc->first - 2;
			const InelasticProcess *p = ArAllData_.ArExper_.FindInelastic(ID);
			if (NULL != p) {
				cross += (*p)(E);
				continue;
			}
			if (!skip_warning) {
				std::cerr << GetName() << "::generate_total_cross_table: Warning: process #" << proc->first << " cross section not found" << std::endl;
				skip_warning = true;
			}
		}
		total_cross_.push_back(E, cross);
	}
	return true;
}

bool ArgonParticle::generate_theta_table(void)
{
	std::cout<<"Generating theta probability function tables..."<<std::endl;
	if (!ArAllData_.ArExper_.isValid()) {
		std::cout << "Error: Invalid experimental data" << std::endl;
		return false;
	}
	//TODO: use EnergyScanner
	ColoredRange energy_range = ColoredInterval (0, 0.01, 1e-4) + ColoredInterval (0, gSettings.PhysConsts()->XS_el_En_maximum, 1e-3) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/10) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/10) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/200) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/200);		//fine area
	energy_range.Trim(0, gSettings.ProgConsts()->maximal_energy);
	std::vector<double> diff_XS, F, thetas;
	diff_XS.resize(gSettings.ProgConsts()->angle_discretization, 0);
	thetas.resize(gSettings.ProgConsts()->angle_discretization, 0);
	F.resize(gSettings.ProgConsts()->angle_discretization, 0);
	for (long int E_i=0, E_ind_end_ = energy_range.NumOfIndices(); E_i!=E_ind_end_; ++E_i) {
		double E = energy_range.Value(E_i);
		std::size_t i_end_ = diff_XS.size();
		for (std::size_t i=0; i!=i_end_; ++i) {
			thetas[i] = i*M_PI/(double)(i_end_-1);
			F[i]=0;
			diff_XS[i]=ArAllData_.argon_cross_elastic_diff(E, thetas[i], 0);
			if (i!=0) {
				F[i] = F[i-1] + 0.5*(diff_XS[i]+diff_XS[i-1])*(thetas[i]-thetas[i-1]);//Integral
			}
		}
		double norm = 1.0 / F[i_end_ - 1];
		for (std::size_t i=0; i!=i_end_; ++i)
			F[i] *= norm;//normalize probability function
		for (std::size_t th_i=0, th_i_end_=thetas.size(); th_i!=th_i_end_;++th_i)
			theta_table_->push(E, thetas[th_i], F[th_i]);
	}
	return true;
}

bool ArgonParticle::generate_time_delay_spin_flip_table(void)
{
	std::cout << "Generating spin flip time delay tables..." << std::endl;
	if (!ArAllData_.ArExper_.isValid()) {
		std::cout << "Error: Invalid experimental data" << std::endl;
		return false;
	}
	//TODO: use EnergyScanner
	ColoredRange energy_Y_range =
		ColoredInterval(En_3o2_ - 30 * Width_3o2_, En_3o2_ + 30 * Width_3o2_, Width_3o2_ / 10) +	//coarse area
		ColoredInterval(En_3o2_ - 15 * Width_3o2_, En_3o2_ + 15 * Width_3o2_, Width_3o2_ / 200)+		//fine area
		ColoredInterval(En_1o2_ - 30 * Width_1o2_, En_1o2_ + 30 * Width_1o2_, Width_1o2_ / 10) +	//coarse area
		ColoredInterval(En_1o2_ - 15 * Width_1o2_, En_1o2_ + 15 * Width_1o2_, Width_1o2_ / 200);		//fine area
	energy_Y_range.Trim(0, gSettings.ProgConsts()->maximal_energy);
	for (long int Ey_i = 0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i != Ey_ind_end_; ++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		for (std::size_t th_i = 0, th_i_end_ = gSettings.ProgConsts()->angle_discretization; th_i != th_i_end_; ++th_i) {
			double theta = th_i*M_PI / (th_i_end_ - 1);
			time_delay_spin_flip_table_->push(Ey, theta, ArAllData_.argon_delay_spin_flip(Ey, theta));
		}
	}
	return true;
}

bool ArgonParticle::generate_time_delay_spin_nonflip_table(void)
{
	std::cout << "Generating spin nonflip time delay tables..." << std::endl;
	if (!ArAllData_.ArExper_.isValid()) {
		std::cout << "Error: Invalid experimental data" << std::endl;
		return false;
	}
	//TODO: use EnergyScanner
	ColoredRange energy_Y_range =
		ColoredInterval(En_3o2_ - 30 * Width_3o2_, En_3o2_ + 30 * Width_3o2_, Width_3o2_ / 10) +	//coarse area
		ColoredInterval(En_3o2_ - 15 * Width_3o2_, En_3o2_ + 15 * Width_3o2_, Width_3o2_ / 200)+		//fine area
		ColoredInterval(En_1o2_ - 30 * Width_1o2_, En_1o2_ + 30 * Width_1o2_, Width_1o2_ / 10) +	//coarse area
		ColoredInterval(En_1o2_ - 15 * Width_1o2_, En_1o2_ + 15 * Width_1o2_, Width_1o2_ / 200);		//fine area
	energy_Y_range.Trim(0, gSettings.ProgConsts()->maximal_energy);
	for (long int Ey_i = 0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i != Ey_ind_end_; ++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		for (std::size_t th_i = 0, th_i_end_ = gSettings.ProgConsts()->angle_discretization; th_i != th_i_end_; ++th_i) {
			double theta = th_i*M_PI / (th_i_end_ - 1);
			time_delay_spin_nonflip_table_->push(Ey, theta, ArAllData_.argon_delay_spin_nonflip(Ey, theta));
		}
	}
	return true;
}

bool ArgonParticle::generate_time_delay_spin_nonflip_prob_table(void)
{
	std::cout << "Generating spin nonflip time delay probability tables..." << std::endl;
	if (!ArAllData_.ArExper_.isValid()) {
		std::cout << "Error: Invalid experimental data" << std::endl;
		return false;
	}
	//TODO: use EnergyScanner
	ColoredRange energy_Y_range =
		ColoredInterval(En_3o2_ - 30 * Width_3o2_, En_3o2_ + 30 * Width_3o2_, Width_3o2_ / 10) +	//coarse area
		ColoredInterval(En_3o2_ - 15 * Width_3o2_, En_3o2_ + 15 * Width_3o2_, Width_3o2_ / 200)+		//fine area
		ColoredInterval(En_1o2_ - 30 * Width_1o2_, En_1o2_ + 30 * Width_1o2_, Width_1o2_ / 10) +	//coarse area
		ColoredInterval(En_1o2_ - 15 * Width_1o2_, En_1o2_ + 15 * Width_1o2_, Width_1o2_ / 200);		//fine area
	energy_Y_range.Trim(0, gSettings.ProgConsts()->maximal_energy);
	for (long int Ey_i = 0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i != Ey_ind_end_; ++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		for (std::size_t th_i = 0, th_i_end_ = gSettings.ProgConsts()->angle_discretization; th_i != th_i_end_; ++th_i) {
			double theta = th_i*M_PI / (th_i_end_ - 1);
			time_delay_spin_nonflip_prob_table_->push(Ey, theta, ArAllData_.argon_delay_spin_nonflip_prob(Ey, theta));
		}
	}
	return true;
}

bool ArgonParticle::generate_ResNBrS_spectrum_table(void)
{
	std::cout << "Generating resonance NBrS spectrum probability tables..." << std::endl;
	if (!ArAllData_.ArExper_.isValid()) {
		std::cout << "Error: Invalid experimental data" << std::endl;
		return false;
	}
	total_resonance_NBrS_spectrum_.clear();
	total_resonance_NBrS_spectrum_.setOrder(1);
	total_resonance_NBrS_spectrum_.setNused(2);
	total_resonance_NBrS_spectrum_.use_leftmost(true);
	total_resonance_NBrS_spectrum_.use_rightmost(true);
	total_resonance_NBrS_spectrum_.set_leftmost(0);
	total_resonance_NBrS_spectrum_.set_rightmost(0);

	double En_center1 = ArAllData_.ArExper_.First_excitation_En - En_1o2_;
	double En_center2 = ArAllData_.ArExper_.First_excitation_En - En_3o2_;

	ColoredRange energy_range =
		ColoredInterval(En_center1 - 4 * Width_1o2_, En_center1 + 4 * Width_1o2_, Width_1o2_ / 100) +
		ColoredInterval(En_center2 - 4 * Width_3o2_, En_center1 + 4 * Width_3o2_, Width_3o2_ / 100);
	std::vector <double> Fprobs;
	Fprobs.resize(energy_range.NumOfIndices(), 0);
	double Int = 0;
	double W_prev = 0;
	double Spec_prev = 0;
	for (long int E_i = 0, E_i_end_ = energy_range.NumOfIndices(); E_i != E_i_end_; ++E_i) {
		double W = energy_range.Value(E_i);
		double Spec = ArAllData_.argon_ResNBrS_spectrum(W, 0);
		Fprobs[E_i] = 0;
		if (E_i != 0) {
			Fprobs[E_i] = Fprobs[E_i - 1] + 0.5*(Spec + Spec_prev)*(W - W_prev);//Integral
		}
		W_prev = W;
		Spec_prev = Spec;
	}
	for (long int E_i = 0, E_i_end_ = energy_range.NumOfIndices(); E_i != E_i_end_; ++E_i) {
		Fprobs[E_i] /= Fprobs[E_i_end_ - 1];//normalize probability function
		double W = energy_range.Value(E_i);
		total_resonance_NBrS_spectrum_.push_back(W, Fprobs[E_i]);
	}
	return true;
}

double ArgonParticle::GetCrossSection(const Particle *target, double E) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GetCrossSection: Error: NULL target" << std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end() == procs) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported target particle \"" << target->GetName() << "\"" << std::endl;
		return 0;
	}
	if (target->GetName() == ELECTRON_NAME) {
		double output = total_cross_(E, E);
		if (std::isnan(output))
			output = 0;
		return total_cross_(E, E);
	}
	return 0;
}

double ArgonParticle::GetCrossSection(const Particle *target, double E, unsigned int process) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GetCrossSection: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) {
			return total_cross_elastic_(E, E);
		}
		if (1 == process) {
			return ArAllData_.argon_ResNBrS_XS(E);
		}
		short ID = process - 2;
		const InelasticProcess *p = ArAllData_.ArExper_.FindInelastic(ID);
		if (NULL!=p)
			return (*p)(E);
	}
	return 0;
}

double ArgonParticle::GetCrossSection(const Particle *target, double E, double theta, unsigned int process) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GetCrossSection: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) {
			return ArAllData_.argon_cross_elastic_diff(E, theta);
		}
		return GetCrossSection(target, E, process) / (4*M_PI); //considered isotropic for other processes
	}
	return 0;
}

std::vector<const Particle*> ArgonParticle::GetFinalStates(const Particle *target, double E, double theta, unsigned int process) const
{
	std::vector<const Particle*> out;
	if (NULL == target) {
		std::cerr << GetName() << "::GetFinalStates: Error: NULL target"<<std::endl;
		out.push_back(NULL);
		return out;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GetFinalStates: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		out.push_back(NULL);
		return out;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GetFinalStates: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		out.push_back(NULL);
		return out;
	}
	out.push_back(this);
	out.push_back(target);
	return out;
}

//Faster algorithm than in Particle is used.
int ArgonParticle::GenerateProcess(const Particle *target, double E, double Rand) const {
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateProcess: Error: NULL target" << std::endl;
		return -1;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end() == procs) {
		std::cerr << GetName() << "::GenerateProcess: Error: unsupported target particle \"" << target->GetName() << "\"" << std::endl;
		return -1;
	}
	std::size_t n_procs = procs->second.size();
	double CrossSectionTotal = 1.0/GetCrossSection(target, E);
	double ProbabilitySum = 0;
	for (std::size_t ind = 0; ind != n_procs; ++ind) {
		ProbabilitySum += std::max(GetCrossSection(target, E, ind) * CrossSectionTotal, 0.0);
		if (Rand < ProbabilitySum) {
			return ind;
		}
	}
	//Because GetCrossSection(target, E) and GetCrossSection(target, E, ind) are calculated
	//using different algorithms the sum of XS!=XS_total. So there will be some events when ProbabilitySum'll be less than 1.
	//Elastic scaterring (0 proceess) is chosen in such cases.
	return 0;
}

double ArgonParticle::GenerateScatterAngle(const Particle *target, double E, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: NULL target"<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) { //elastic
			return theta_table_->getY(E, Rand);
		}
		//considered isotropic for other processes
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	Rand = Rand*2.0 - 1.0;
	return std::acos(Rand);
}

double ArgonParticle::GenerateEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateEnergyLoss: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateEnergyLoss: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateEnergyLoss: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) {
			long double gamma_f = target->GetMass() / GetMass();
			double EnergyLoss = 2*(1-cos(theta))*E*gamma_f /pow(1 + gamma_f, 2);
			if (gSettings.PhysConsts()->resonance_En_loss) {
				double width_factor = std::pow(Width_3o2_/2.0, 2)/(std::pow(Width_3o2_/2.0, 2) + std::pow((E-En_3o2_)/2.0, 2));
				width_factor += std::pow(Width_1o2_/2.0, 2)/(std::pow(Width_1o2_/2.0, 2) + std::pow((E-En_1o2_)/2.0, 2));
				width_factor *= 0.5**gSettings.PhysConsts()->resonance_En_loss;
				EnergyLoss += width_factor;
			}
			return EnergyLoss;
		}
		if (1 == process) { //same as GeneratePhoton, will return the same energy for the same random numbers
			for (std::size_t i = 0, i_end_ = total_resonance_NBrS_spectrum_.size(); i != i_end_; ++i) {
				if (Rand <= total_resonance_NBrS_spectrum_.getY(i)) {
					if (0 == i)
						return total_resonance_NBrS_spectrum_.getX(i);
					double x0 = total_resonance_NBrS_spectrum_.getX(i - 1);
					double x1 = total_resonance_NBrS_spectrum_.getX(i);
					double y0 = total_resonance_NBrS_spectrum_.getY(i - 1);
					double y1 = total_resonance_NBrS_spectrum_.getY(i);
					return x0 + (x1 - x0)*(Rand - y0) / (y1 - y0);
				}
			}
			return total_resonance_NBrS_spectrum_.getX(total_resonance_NBrS_spectrum_.size()-1);
		}
		short ID = process - 2;
		const InelasticProcess *p = ArAllData_.ArExper_.FindInelastic(ID);
		if (NULL!=p) {
			if ((std::string::npos != p->get_name().find("ION")) || (std::string::npos != p->get_name().find("Ion"))
					|| (std::string::npos != p->get_name().find("ion"))) {
				double EnergyLoss = p->get_En_thresh();
				EnergyLoss += (E - EnergyLoss)/2.0; //Consider that residual energy is equally divided between 2 electrons
				return EnergyLoss;
			}
			return p->get_En_thresh();
		}
	}
	return 0;
}

//For neutral bremsstrahlung or deexcitation. Different from energy loss
double ArgonParticle::GeneratePhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GeneratePhoton: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GeneratePhoton: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GeneratePhoton: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) {
			return 0;
		}
		if (1 == process) {
			for (std::size_t i = 0, i_end_ = total_resonance_NBrS_spectrum_.size(); i != i_end_; ++i) {
				if (Rand <= total_resonance_NBrS_spectrum_.getY(i)) {
					if (0 == i)
						return total_resonance_NBrS_spectrum_.getX(i);
					double x0 = total_resonance_NBrS_spectrum_.getX(i - 1);
					double x1 = total_resonance_NBrS_spectrum_.getX(i);
					double y0 = total_resonance_NBrS_spectrum_.getY(i - 1);
					double y1 = total_resonance_NBrS_spectrum_.getY(i);
					return x0 + (x1 - x0)*(Rand - y0) / (y1 - y0);
				}
			}
			return total_resonance_NBrS_spectrum_.getX(total_resonance_NBrS_spectrum_.size()-1);
		}
		short ID = process - 2;
		const InelasticProcess *p = ArAllData_.ArExper_.FindInelastic(ID);
		if (NULL!=p) {
			if ((std::string::npos != p->get_name().find("ION")) || (std::string::npos != p->get_name().find("Ion")) || (std::string::npos != p->get_name().find("ion")))
				return 0; //Ionization
			return 2 * M_PI * gSettings.PhysConsts()->h_bar_eVs * gSettings.PhysConsts()->light_speed_SI / gSettings.PhysConsts()->Ar_primal_line_nm * 1e-9;
		}
	}
	return 0;
}

double ArgonParticle::GenerateTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) {
			switch (gSettings.PhysConsts()->time_delay_model) {
			case (PhysicalConstants::TimeDelayMode::Precise) : {
				double P = (*time_delay_spin_nonflip_prob_table_)(E, theta);
				P = std::max(P, 0.0);
				if (Rand < P) {
					return (*time_delay_spin_nonflip_table_)(E, theta);
				} else {
					return (*time_delay_spin_flip_table_)(E, theta);
				}
				return 0;
			}
			case (PhysicalConstants::TimeDelayMode::Rough) : {
				double tau = ArAllData_.argon_time_delay_j(E, 3, 2);
				tau += ArAllData_.argon_time_delay_j(E, 1, 2);
				return tau;
			}
			case (PhysicalConstants::TimeDelayMode::None) : {
				return 0;
			}
			}
		}
		return 0;
	}
	return 0;
}

//Untabulated functions:
unsigned int ArgonParticle::GenerateUntabProcess(const Particle *target, double E, double Rand) const
{
	return GenerateProcess(target, E, Rand);
}

double ArgonParticle::GenerateUntabScatterAngle(const Particle *target, double E, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: NULL target"<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) { //elastic
			std::vector<double> diff_XS, thetas;
			diff_XS.resize(gSettings.ProgConsts()->angle_discretization, 0);
			thetas.resize(gSettings.ProgConsts()->angle_discretization, 0);
			for (std::size_t i=0, i_end_=diff_XS.size(); i!=i_end_; ++i) {
				thetas[i] = i*M_PI/(double)(i_end_-1);
				diff_XS[i]=ArAllData_.argon_cross_elastic_diff(E, thetas[i], 0);
			}
			PDF_routine prob_function(thetas, diff_XS);
			return prob_function(Rand);
		}
		//considered isotropic for other processes
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	Rand = Rand*2.0 - 1.0;
	return std::acos(Rand);
}

double ArgonParticle::GenerateUntabEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GenerateEnergyLoss(target, E, theta, process, Rand);
}

double ArgonParticle::GenerateUntabPhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GeneratePhoton(target, E, theta, process, Rand);
}

double ArgonParticle::GenerateUntabTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()== ELECTRON_NAME) {
		if (0 == process) {
			switch (gSettings.PhysConsts()->time_delay_model) {
			case (PhysicalConstants::TimeDelayMode::Precise) : {
				double P = ArAllData_.argon_delay_spin_nonflip_prob(E, theta);
				if (Rand < P) {
					return ArAllData_.argon_delay_spin_nonflip(E, theta);
				} else {
					return ArAllData_.argon_delay_spin_flip(E, theta);
				}
				return 0;
			}
			case (PhysicalConstants::TimeDelayMode::Rough) : {
				double tau = ArAllData_.argon_time_delay_j(E, 3, 2);
				tau += ArAllData_.argon_time_delay_j(E, 1, 2);
				return tau;
			}
			case (PhysicalConstants::TimeDelayMode::None) : {
				return 0;
			}
			}
		}
		return 0;
	}
	return 0;
}
