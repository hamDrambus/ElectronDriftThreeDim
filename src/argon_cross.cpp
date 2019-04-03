#include "argon_cross.h"

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
		energy_range_ =ColoredInterval (0, XS_EL_EN_MINIMUM_, 1e-4) + ColoredInterval (XS_EL_EN_MINIMUM_, 0.1, 5e-4) + ColoredInterval (0.1, 1, 1e-3) +
				ColoredInterval (1, 10, 5e-2) + ColoredInterval (10, XS_EL_EN_MAXIMUM_, 0.1)+
				ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/2) + 	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/2) + 	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/30) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/30); 		//fine area
		break;
	}
	case (InelasticXS): {
		energy_range_ = ColoredInterval (11.5, EN_MAXIMUM_, 0.01);
		break;
	}
	case (ElasticResXS): {
		energy_range_ = ColoredInterval (XS_EL_EN_MINIMUM_*0.1, 0.1, 5e-4) + ColoredInterval (0.1, 1, 1e-3) +
				ColoredInterval (1, 10, 5e-2) + ColoredInterval (10, XS_EL_EN_MAXIMUM_, 0.1) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/2) + 		//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/2) + 		//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/30) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/30);	//fine area
		break;
	}
	case (XSIntegral): {
		energy_range_ = ColoredInterval (0, 0.1, 2e-4) + ColoredInterval (0.1, 1, 8e-4) +
				ColoredInterval (1, 10, 1e-2) + ColoredInterval (10, EN_MAXIMUM_, 0.02) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80) +	//fine area
				ColoredInterval (11.5, EN_MAXIMUM_, 0.003);
		break;
	}
	case (PlotElastic): {
		//from 1e-3 eV to 0.1 eV with step 5e-4 eV, etc.
		energy_range_ = ColoredInterval (XS_EL_EN_MINIMUM_, 0.1, 3e-4) + ColoredInterval (0.1, 1, 9e-4);
		break;
	}
	case (PlotResonance_3o2): {
		energy_range_ = ColoredInterval (En_3o2_ - 110*Width_3o2_, En_3o2_ + 110*Width_3o2_, Width_3o2_/3) +//coarse area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80); 		//fine area
		break;
	}
	case (PlotResonance_1o2): {
		energy_range_ = ColoredInterval (En_1o2_ - 110*Width_1o2_, En_1o2_ + 110*Width_1o2_, Width_1o2_/3) +//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80); 		//fine area
		break;
	}
	case (PlotResonances): { //used for plotting total and differential XS in resonances region
		energy_range_ = ColoredInterval (En_1o2_ - 110*Width_1o2_, En_1o2_ + 110*Width_1o2_, Width_1o2_/3) + 	//coarse area
				ColoredInterval (En_3o2_ - 110*Width_3o2_, En_3o2_ + 110*Width_3o2_, Width_3o2_/3) + 			//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) + 		//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80); 			//fine area
		break;
	}
	case (PlotDiffXS): {  //used for plotting test of differential to total cross sections
		energy_range_ = ColoredInterval (1e-4, 0.1, 7e-4) +
				ColoredInterval (0.1, 1, 7e-3) +
				ColoredInterval (1, EN_MAXIMUM_, 0.086);
			//ColoredInterval(En_1o2_ - 110 * Width_1o2_, En_1o2_ + 110 * Width_1o2_, Width_1o2_ / 3) + 	//coarse area
			//ColoredInterval(En_3o2_ - 110 * Width_3o2_, En_3o2_ + 110 * Width_3o2_, Width_3o2_ / 3) + 	//coarse area
			//ColoredInterval(En_1o2_ - 15 * Width_1o2_, En_1o2_ + 15 * Width_1o2_, Width_1o2_ / 80) + 	//fine area
			//ColoredInterval(En_3o2_ - 15 * Width_3o2_, En_3o2_ + 15 * Width_3o2_, Width_3o2_ / 80); 	//fine area
		break;
	}
	case (PlotInelasticXS): {
		energy_range_ = ColoredInterval (11.5, XS_EL_EN_MAXIMUM_, 0.007) +
				ColoredInterval (XS_EL_EN_MAXIMUM_, 100, 0.1);
		break;
	}
	case (PlotElasticResXS): {
		energy_range_ = ColoredInterval (XS_EL_EN_MINIMUM_, 0.1, 3e-4) + ColoredInterval (0.1, 1, 9e-4) +
				ColoredInterval (1, 10, 7e-3) + ColoredInterval (10, XS_EL_EN_MAXIMUM_, 0.016) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(XS_EL_EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(XS_EL_EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(XS_EL_EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) +//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(XS_EL_EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80);//fine area
		break;
	}
	case (PlotAllXS): {
		energy_range_ = ColoredInterval (1e-4, 0.1, 2e-4) + ColoredInterval (0.1, 1, 7e-4) +
				ColoredInterval (1, 10, 7e-3) + ColoredInterval (10, XS_EL_EN_MAXIMUM_, 0.016) +
				ColoredInterval (En_1o2_ - 110*Width_1o2_, std::min(XS_EL_EN_MAXIMUM_, En_1o2_ + 110*Width_1o2_), Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 110*Width_3o2_, std::min(XS_EL_EN_MAXIMUM_, En_3o2_ + 110*Width_3o2_), Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(XS_EL_EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) +//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(XS_EL_EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80) +//fine area
				ColoredInterval (11.5, XS_EL_EN_MAXIMUM_, 0.003);
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

double InelasticProcess::operator ()(double E) //returns cross section in 1e-16 cm^2
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

double InelasticProcess::BB_XS(double E)
{
	double gamma = (E+e_mass_eVconst)/e_mass_eVconst;
	double gamma2 = gamma*gamma;
	double beta2 = 1.0-1.0/gamma2;
	double BBconst = 8.0*M_PI*a_bohr_SIconst*a_bohr_SIconst*Ry_eVconst*Ry_eVconst/e_mass_eVconst;
	return (Oscillator_strength_/(En_threshold_*beta2))*log(beta2*gamma2*e_mass_eVconst/(2.0*En_threshold_)-beta2)*BBconst*E/(E+En_threshold_+ArExper_->E_Ionization);
}

double InelasticProcess::Exp_XS(double E)
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

ArExperimental::ArExperimental(void): total_elastic_cross(3, 5) /*fit by 3rd order polynomial*/, max_process_ID(0), E_Ionization(0)
{
	std::cout << "Reading Ar exprimental data..." << std::endl;
	std::ifstream inp;
	inp.open("data/ArScatteringCross.dat");
	std::string line, word;
	while (!inp.eof()) {
		std::getline(inp, line);
		if (line.size()>=2)
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double k = a_h_bar_2e_m_e_SIconst*sqrt(std::stod(word));
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double XS = std::stod(word); //Cross section is stored in 1e-20 m^2 in both files and this program.
		total_elastic_cross.push(k, XS);
	}
	inp.close();
	for (int l=0;l<6;++l)
		if (l<4) {
			phase_shifts_pos_.push_back(DataVector(3,4)); /*interpolation by 3rd order polynomial because there's data for 11eV*/
			phase_shifts_pos_.back().use_rightmost(true);
			phase_shifts_neg_.push_back(DataVector(3,4));
			phase_shifts_neg_.back().use_rightmost(true);
		} else {
			phase_shifts_pos_.push_back(DataVector(3,5)); /*fit by 3rd order polynomial TODO: tests, tests*/
			phase_shifts_pos_.back().use_rightmost(true);
			phase_shifts_neg_.push_back(DataVector(3,4));
			phase_shifts_neg_.back().use_rightmost(true);
		}
	inp.open("data/McEachranArPhaseShifts.dat");
	unsigned int MAX_L = 6; //L == (MAX_L-1) is the last one
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
					phase_shifts_neg_[l].push(k, vals[l]);
				phase_shifts_pos_[l].push(k, vals[l]);
			}
		}
		for (int l=0; l<(MAX_L-1); ++l) {
			if (is_d_1[l])
				phase_shifts_neg_[l+1].push(k, vals_1[l]);
		}
	}
	inp.close();
	//loading ionization and excitation data
	inp.open("data/ArIonizations_Magboltz.dat");
	read_inelastic(inp, ionizations);
	inp.close();
	if (ionizations.size()) {
		E_Ionization = DBL_MAX;
		for (int i =0, end_ = ionizations.size();i!=end_;++i)
			E_Ionization = std::min(ionizations[i].get_En_thresh(), E_Ionization);
	} 
	inp.open("data/ArExitations_Magboltz.dat");
	read_inelastic(inp, excitations);
	inp.close();
	if (excitations.size()) {
		First_excitation_En = DBL_MAX;
		for (int i = 0, end_ = excitations.size(); i != end_; ++i)
			First_excitation_En = std::min(excitations[i].get_En_thresh(), First_excitation_En);
	}
	std::cout << "Finished reading Ar exprimental data." << std::endl;
}

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

InelasticProcess * ArExperimental::FindInelastic(short ID)
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

unsigned int ArExperimental::max_L (long double k)
{
	return std::max(phase_shifts_pos_.size(), phase_shifts_neg_.size())-1;
}

void ArExperimental::phase_shift (long double k, unsigned int l, long double &ps_pos,  long double &ps_neg)
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
void ArAllData::argon_phase_values_exp(long double k, unsigned int l, long double &ps_p, long double &ps_n)
{
	long double angle_pos, angle_neg;
	ArExper_.phase_shift(k, l, angle_pos, angle_neg);
	ps_p = angle_pos;
	ps_n = angle_neg;
	if (1 == l) {
		long double E = std::pow(k / a_h_bar_2e_m_e_SIconst, 2.0);//recalculation from atomic units to energy
		long double cot3o2 = -2 * (E - En_3o2_) / Width_3o2_;
		long double cot1o2 = -2 * (E - En_1o2_) / Width_1o2_;
		long double cos3o2 = cot3o2 / std::sqrt(1 + cot3o2*cot3o2);
		long double cos1o2 = cot1o2 / std::sqrt(1 + cot1o2*cot1o2);
		ps_p += std::acos(cos3o2);
		ps_n += std::acos(cos1o2);
	}
}

//k is in atomic units
void ArAllData::argon_phase_values_MERT5(long double k, unsigned int l, long double &ps_p, long double &ps_n)
{
	//see Kurokawa Phys. Rev. A84 2011, MERT5+ fit http://dx.doi.org/10.1103/PhysRevA.84.062717
	double A = -1.365;
	double D = 80.5;
	double F = -153;
	double G = 31.0;
	double A1 = 8.8;
	double H = 29.7;
	double alpha_d = 11.08;
	double alpha_q = 0.0;
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
long double ArAllData::argon_time_delay_j(long double E, int J, int L)
{
	if (L!=2)
		return 0;
	if (((J != L - 1) && (J != L + 1)) || J < 0) {
		std::cerr << "ArAllData::argon_time_delay_j: Error: No momentum J=" << J << "/2 for e(s=1/2)-Ar(s=0) scattering with L =" << L / 2 << std::endl;
		return 0;
	}
	if (J == L-1) {
		return 2* h_bar_eVconst / (1 + std::pow(2 * (E - En_1o2_) / Width_1o2_, 2))/Width_1o2_;
	} else {
		return 2* h_bar_eVconst / (1 + std::pow(2 * (E - En_3o2_) / Width_3o2_, 2))/Width_3o2_;
	}
}

//Probability on scatter at certain angle going through total momentum J/2 and orbital momentum L/2. (J&L are expressed in halves)
//Can be negative.
//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
//E in eV
long double ArAllData::argon_scatter_probability_j(long double E, long double theta, int J, int L, int mode)
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
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) = 0;
	unsigned int L_MAX = 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = L_MAX_;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (PHASES_EN_MINIMUM_>E)
			E = PHASES_EN_MINIMUM_;
		k = a_h_bar_2e_m_e_SIconst*sqrt(E);
		if (E<THRESH_E_PHASES_) {
			L_MAX = L_MAX_;
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
long double ArAllData::argon_cross_elastic_diff(long double E, long double theta, int mode) {
	//different formulas are used for E<0.24eV and E>0.24eV!
	unsigned int L_MAX = 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	if (1!=mode&&2!=mode)
	{
		if (PHASES_EN_MINIMUM_>E)
			k = a_h_bar_2e_m_e_SIconst*sqrt(PHASES_EN_MINIMUM_);
	}
	long double cross = argon_scatter_spin_flip_amplitude_sq(E, theta, mode);
	cross += argon_scatter_spin_nonflip_amplitude_sq(E, theta, mode);
	cross *= M_PI / (2.0*pow(k, 2));
	cross = std::max(cross, (long double)0);
	return cross*a_bohr_SIconst*a_bohr_SIconst; //const is multiplied by 1e10
}

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
long double ArAllData::argon_scatter_spin_flip_amplitude_sq(long double E, long double theta, int mode)
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) = 0;
	unsigned int L_MAX = 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = L_MAX_;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (PHASES_EN_MINIMUM_>E)
			E = PHASES_EN_MINIMUM_;
		k = a_h_bar_2e_m_e_SIconst*sqrt(E);
		if (E<THRESH_E_PHASES_) {
			L_MAX = L_MAX_;
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

long double ArAllData::argon_scatter_spin_nonflip_amplitude_sq(long double E, long double theta, int mode)
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) = 0;
	unsigned int L_MAX = 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = L_MAX_;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (PHASES_EN_MINIMUM_>E)
			E = PHASES_EN_MINIMUM_;
		k = a_h_bar_2e_m_e_SIconst*sqrt(E);
		if (E<THRESH_E_PHASES_) {
			L_MAX = L_MAX_;
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
long double ArAllData::argon_cross_elastic(long double E, int mode) //Tabulation of this function must be done carefully. Do not forget near 0 case.
{
	long double cross = 0;
	switch (mode) {
	case 1: {
		long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
														// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
		unsigned int L_MAX = L_MAX_;
		for (unsigned int l = 0; l <= L_MAX; ++l) {
			long double phase_l_p = 0;
			long double phase_l_n = 0;
			argon_phase_values_MERT5(k, l, phase_l_p, phase_l_n);
			cross += (l + 1)*sin(phase_l_p)*sin(phase_l_p) + l*sin(phase_l_n)*sin(phase_l_n);
		}
		cross *= 4.0*M_PI / pow(k, 2);
		cross *= a_bohr_SIconst*a_bohr_SIconst;
		break;
	}
	case 2: {
		//this is extrapolation of experimental data
		cross = ArExper_.total_elastic_cross(a_h_bar_2e_m_e_SIconst*sqrt(E), a_h_bar_2e_m_e_SIconst*sqrt(E));
		//this part is calculated using phase shifts
		cross += argon_cross_resonance_1o2(E);
		cross += argon_cross_resonance_3o2(E);
		break;
	}
	case 3: {
		long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
		unsigned int L_MAX = ArExper_.max_L(k);;
		for (unsigned int l = 0; l <= L_MAX; ++l) {
			long double phase_l_p = 0;
			long double phase_l_n = 0;
			argon_phase_values_exp(k, l, phase_l_p, phase_l_n); //includes Feshbach resonances
			cross += (l + 1)*sin(phase_l_p)*sin(phase_l_p) + l*sin(phase_l_n)*sin(phase_l_n);
		}
		cross *= 4.0*M_PI / pow(k, 2);
		cross *= a_bohr_SIconst*a_bohr_SIconst;
		break;
	}
	case 4: {
		if (E<XS_EL_EN_MINIMUM_) {
			cross = 7.491 + E*(argon_cross_elastic(XS_EL_EN_MINIMUM_, 1) - 7.491) / XS_EL_EN_MINIMUM_; //linear behavior. MERT5 is used at XS_EL_EN_MINIMUM_
		}
		else {
			if (E <THRESH_E_XS_) {
				cross = argon_cross_elastic(E, 1); //MERT5
			}
			else {
				cross = argon_cross_elastic(E, 2); //Experiment
			}
		}
		break;
	}
	default: {
#ifndef D_EN_SMOOTH_
		cross = argon_cross_elastic(E, 4);
#else
		if (E<XS_EL_EN_MINIMUM_) {
			cross = 7.491 + E*(argon_cross_elastic(XS_EL_EN_MINIMUM_, 1) - 7.491) / XS_EL_EN_MINIMUM_; //linear behavior. MERT5 is used at XS_EL_EN_MINIMUM_
		}
		else {
			double E_l = THRESH_E_XS_ - D_EN_SMOOTH_;
			double E_r = THRESH_E_XS_ + D_EN_SMOOTH_;
			if (E <E_l) {
				cross = argon_cross_elastic(E, 1);//MERT5
				break;
			}
			if (E >E_r) {
				cross = argon_cross_elastic(E, 2);//Experiment
				break;
			}
			cross = ((E_r - E)*argon_cross_elastic(E, 1) + (E - E_l)*argon_cross_elastic(E, 2)) / (E_r - E_l); //linear smoothing between MERT5 and experiment
		}
#endif
		break;
	}
	}
	return cross;
}

/*long double ArAllData::argon_delay_1o2_probability(long double E, long double theta)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	long double prob = argon_scatter_probability_j(E, theta, 1, 2);
	prob = std::max(prob, (long double)0); //TODO: since p can be <0 the formula for the delay is wrong, because 1>=P>=0 for observables
	prob *= M_PI / (2.0*pow(k, 2)) *a_bohr_SIconst*a_bohr_SIconst;
	//if (prob<0)
	//	std::cout<<"ArAllData: Warning: P1/2 = "<<prob<<"\t at E = "<<E<<"\t th = "<<theta<<std::endl;
	return prob / argon_cross_elastic_diff(E, theta);
}

long double ArAllData::argon_delay_3o2_probability(long double E, long double theta)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	long double prob = argon_scatter_probability_j(E, theta, 3, 2);
	prob = std::max(prob, (long double)0); //TODO: since p can be <0 the formula for the delay is wrong, because 1>=P>=0 for observables
	prob *= M_PI / (2.0*pow(k, 2)) *a_bohr_SIconst*a_bohr_SIconst;
	//if (prob<0)
	//	std::cout<<"ArAllData: Warning: P1/2 = "<<prob<<"\t at E = "<<E<<"\t th = "<<theta<<std::endl;
	return prob / argon_cross_elastic_diff(E, theta);
}*/

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
long double ArAllData::argon_delay_spin_flip (long double E, long double theta, int mode)
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) = 0;
	unsigned int L_MAX = 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = L_MAX_;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (PHASES_EN_MINIMUM_>E)
			E = PHASES_EN_MINIMUM_;
		k = a_h_bar_2e_m_e_SIconst*sqrt(E);
		if (E<THRESH_E_PHASES_) {
			L_MAX = L_MAX_;
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
		l_p *= 2; //all cosines are taken from double angles
		l_n *= 2;
		cos_sum +=AP1(cos_th, l, 1)*cos(l_p);
		cos_sum +=-1*AP1(cos_th, l, 1)*cos(l_n);
		sin_sum +=AP1(cos_th, l, 1)*sin(l_p);
		sin_sum +=-1*AP1(cos_th, l, 1)*sin(l_n);
		if (l!=0) {
			sin_delayed_sum+= -2*argon_time_delay_j(E, 2*l -1, 2*l) * -1*AP1(cos_th, l, 1) * sin(l_n);
			cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l -1, 2*l) * -1*AP1(cos_th, l, 1) *cos(l_n);
		}
		sin_delayed_sum+= -2*argon_time_delay_j(E, 2*l +1, 2*l) * AP1(cos_th, l, 1) *sin(l_p);
		cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l +1, 2*l) * AP1(cos_th, l, 1) *cos(l_p);
	}
	long double cos_2 = cos_sum*cos_sum;
	long double sin_2 = sin_sum*sin_sum;
	return cos_2/(cos_2 + sin_2) * (sin_delayed_sum * sin_2/cos_2 + cos_delayed_sum /cos_sum);
}

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
long double ArAllData::argon_delay_spin_nonflip (long double E, long double theta, int mode)
{
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double, unsigned int, long double &, long double &) = 0;
	unsigned int L_MAX = 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = L_MAX_;
		phase_values = &ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values = &ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (PHASES_EN_MINIMUM_>E)
			E = PHASES_EN_MINIMUM_;
		k = a_h_bar_2e_m_e_SIconst*sqrt(E);
		if (E<THRESH_E_PHASES_) {
			L_MAX = L_MAX_;
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
			sin_delayed_sum+= -2*argon_time_delay_j(E, 2*l -1, 2*l) * P1(cos_th, l) * (l)*sin(l_n);
			cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l -1, 2*l) * P1(cos_th, l) * (l)*cos(l_n);
		}
		sin_delayed_sum+= -2*argon_time_delay_j(E, 2*l +1, 2*l) * P1(cos_th, l) * (l+1)*sin(l_p);
		cos_delayed_sum+= 2*argon_time_delay_j(E, 2*l +1, 2*l) * P1(cos_th, l) * (l+1)*cos(l_p);
	}
	long double cos_2 = cos_sum*cos_sum;
	long double sin_2 = sin_sum*sin_sum;
	return cos_2/(cos_2 + sin_2) * (sin_delayed_sum * sin_2/cos_2 + cos_delayed_sum /cos_sum);
}

long double ArAllData::argon_delay_spin_nonflip_prob (long double E, long double theta, int mode)
{
	long double B = argon_scatter_spin_flip_amplitude_sq(E, theta, mode);
	long double A = argon_scatter_spin_nonflip_amplitude_sq(E, theta, mode);
	return A/(A+B);
}

long double ArAllData::argon_ResNBrS_spectrum(long double W, long double E) //Normalization constant is arbitrary. Not dependant on E
{
	return std::exp(-0.5*std::pow((W - (ArExper_.First_excitation_En - En_1o2_)) / Width_1o2_, 2))
		+ std::exp(-0.5*std::pow((W - (ArExper_.First_excitation_En - En_3o2_)) / Width_3o2_, 2));
}

long double ArAllData::argon_ResNBrS_XS(long double E) //Normalization constant is taken from "global_definitions.h"
{
	double width_factor = std::pow(Width_3o2_ / 2.0, 2) / (std::pow(Width_3o2_ / 2.0, 2) + std::pow((E - En_3o2_) / 2.0, 2));
	width_factor += std::pow(Width_1o2_ / 2.0, 2) / (std::pow(Width_1o2_ / 2.0, 2) + std::pow((E - En_1o2_) / 2.0, 2));
	return RESONANCE_NBrS_XS_*width_factor;
}

long double ArAllData::argon_back_scatter_prob(long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
													// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross = 0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper_.max_L(k);
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double phase_l_p = 0;
		long double phase_l_n = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, phase_l_p, phase_l_n);
		else
			argon_phase_values_exp(k, l, phase_l_p, phase_l_n);
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double phase_f_p = phase_l_p;
			long double phase_f_n = phase_l_n;
			if (l != f) {
				if (E<THRESH_E_PHASES_)
					argon_phase_values_MERT5(k, f, phase_f_p, phase_f_n);
				else
					argon_phase_values_exp(k, f, phase_f_p, phase_f_n);
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
long double ArAllData::argon_TM_forward(long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
													// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross = 0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper_.max_L(k);
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double phase_l_p = 0;
		long double phase_l_n = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, phase_l_p, phase_l_n);
		else
			argon_phase_values_exp(k, l, phase_l_p, phase_l_n);
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double phase_f_p = phase_l_p;
			long double phase_f_n = phase_l_n;
			if (l != f) {
				if (E<THRESH_E_PHASES_)
					argon_phase_values_MERT5(k, f, phase_f_p, phase_f_n);
				else
					argon_phase_values_exp(k, f, phase_f_p, phase_f_n);
			}
			long double cos_l_f = cos(phase_l_p - phase_f_p);
			W += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl_transf(l, f, 0, 1, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl(l, f, 0, 1, 1e-5);
		}
	}
	return W / cross;
}

long double ArAllData::argon_TM_backward(long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
													// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross = 0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper_.max_L(k);
	for (unsigned int l = 0; l <= L_MAX; ++l) {
		long double phase_l_p = 0;
		long double phase_l_n = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, phase_l_p, phase_l_n);
		else
			argon_phase_values_exp(k, l, phase_l_p, phase_l_n);
		for (unsigned int f = l; f <= L_MAX; ++f) {
			long double phase_f_p = phase_l_p;
			long double phase_f_n = phase_l_n;
			if (l != f) {
				if (E<THRESH_E_PHASES_)
					argon_phase_values_MERT5(k, f, phase_f_p, phase_f_n);
				else
					argon_phase_values_exp(k, f, phase_f_p, phase_f_n);
			}
			long double cos_l_f = cos(phase_l_p - phase_f_p);
			W += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl_transf(l, f, -1, 0, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross += ((l == f) ? 1.0 : 2.0)*(2 * l + 1)*(2 * f + 1)*sin(phase_l_p)*sin(phase_f_p)*cos_l_f*Int_PlPl(l, f, -1, 0, 1e-5);
		}
	}
	return W / cross;
}

long double ArAllData::argon_cross_resonance_3o2(long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E);
	long double ps_p_0, ps_n_0;
	ArExper_.phase_shift(k, 1, ps_p_0, ps_n_0);
	long double ps_p = ps_p_0, ps_n = ps_n_0;
	long double cot3o2 = -2 * (E - En_3o2_) / Width_3o2_;
	long double cot1o2 = -2 * (E - En_1o2_) / Width_1o2_;
	long double cos3o2 = cot3o2 / std::sqrt(1 + cot3o2*cot3o2);
	long double cos1o2 = cot1o2 / std::sqrt(1 + cot1o2*cot1o2);
	ps_p += std::acos(cos3o2);
	ps_n += std::acos(cos1o2);

	long double cross = 2 * (pow(sin(ps_p), 2) - pow(sin(ps_p_0), 2));
	cross *= 4.0*M_PI / pow(k, 2);
	return cross * a_bohr_SIconst * a_bohr_SIconst;
}

long double ArAllData::argon_cross_resonance_1o2(long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E);
	long double ps_p_0, ps_n_0;
	ArExper_.phase_shift(k, 1, ps_p_0, ps_n_0);
	long double ps_p = ps_p_0, ps_n = ps_n_0;
	long double cot3o2 = -2 * (E - En_3o2_) / Width_3o2_;
	long double cot1o2 = -2 * (E - En_1o2_) / Width_1o2_;
	long double cos3o2 = cot3o2 / std::sqrt(1 + cot3o2*cot3o2);
	long double cos1o2 = cot1o2 / std::sqrt(1 + cot1o2*cot1o2);
	ps_p += std::acos(cos3o2);
	ps_n += std::acos(cos1o2);

	long double cross = (pow(sin(ps_n), 2) - pow(sin(ps_n_0), 2));
	cross *= 4.0*M_PI / pow(k, 2);
	return cross * a_bohr_SIconst * a_bohr_SIconst;
}

void ArDataTables::read_data (std::ifstream &inp, DataVector &data, long double y_factor)
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

ArDataTables::ArDataTables():
	total_cross_elastic_fname("data_derived/total_cross_section_elastic.dat"),
	integral_table_fname("data_derived/cross_integrals.dat"),
	theta_table_fname("data_derived/theta_probabilities.dat"),
	time_delay_spin_nonflip_prob_fname("data_derived/time_delay_spin_nonflip_probabilities.dat"),
	time_delay_spin_flip_fname("data_derived/time_delay_spin_flip.dat"),
	time_delay_spin_nonflip_fname("data_derived/time_delay_spin_nonflip.dat"),
	total_resonance_NBrS_spectrum_fname("data_derived/ResNBrS_spectrum.dat")
{
	std::cout<<"Loading Ar data tables"<<std::endl;
	integral_table_ = new FunctionTable; //shared between processes - read only
	theta_table_ = new FunctionTable; //shared between processes - read only
	time_delay_spin_flip_table_= new FunctionTable; //shared between processes - read only
	time_delay_spin_nonflip_table_= new FunctionTable; //shared between processes - read only
	time_delay_spin_nonflip_prob_table_= new FunctionTable; //shared between processes - read only

	ensure_file(total_cross_elastic_fname);
	ensure_file(integral_table_fname);
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
		generate_total_cross_elastic_table();
		str.open(total_cross_elastic_fname, std::ios_base::trunc);
		total_cross_elastic_.write(str, "E[eV]\tXS elastic [1e-20 m^2]");
		str.close();
	} else {
		total_cross_elastic_.read(inp);
		inp.close();
		if (total_cross_elastic_.size() < total_cross_elastic_.getNused()) {
			std::cout << "Failed to load \"" << total_cross_elastic_fname << "\"" << std::endl;
			generate_total_cross_elastic_table();
			str.open(total_cross_elastic_fname, std::ios_base::trunc);
			total_cross_elastic_.write(str, "E[eV]\tXS elastic [1e-20 m^2]");
			str.close();
		}
	}

	inp.open(integral_table_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout<<"Failed to load \""<<integral_table_fname<<"\""<<std::endl;
		generate_integral_table();
		str.open(integral_table_fname, std::ios_base::trunc|std::ios_base::binary);
		integral_table_->write(str);
		str.close();
	} else {
		integral_table_->read(inp);
		inp.close();
		if (integral_table_->is_empty()) {
			std::cout<<"Failed to load \""<<integral_table_fname<<"\""<<std::endl;
			generate_integral_table();
			str.open(integral_table_fname, std::ios_base::trunc|std::ios_base::binary);
			integral_table_->write(str);
			str.close();
		}
	}

	inp.open(theta_table_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout<<"Failed to load \""<<theta_table_fname<<"\""<<std::endl;
		generate_theta_table();
		str.open(theta_table_fname, std::ios_base::trunc|std::ios_base::binary);
		theta_table_->write(str);
		str.close();
	} else {
		theta_table_->read(inp);
		inp.close();
		if (theta_table_->is_empty()) {
			std::cout<<"Failed to load \""<<theta_table_fname<<"\""<<std::endl;
			generate_theta_table();
			str.open(theta_table_fname, std::ios_base::trunc|std::ios_base::binary);
			theta_table_->write(str);
			str.close();
		}
	}

	inp.open(time_delay_spin_nonflip_prob_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << time_delay_spin_nonflip_prob_fname << "\"" << std::endl;
		generate_time_delay_spin_nonflip_prob_table();
		str.open(time_delay_spin_nonflip_prob_fname, std::ios_base::trunc | std::ios_base::binary);
		time_delay_spin_nonflip_prob_table_->write(str);
		str.close();
	} else {
		time_delay_spin_nonflip_prob_table_->read(inp);
		inp.close();
		if (time_delay_spin_nonflip_prob_table_->is_empty()) {
			std::cout << "Failed to load \"" << time_delay_spin_nonflip_prob_fname << "\"" << std::endl;
			generate_time_delay_spin_nonflip_prob_table();
			str.open(time_delay_spin_nonflip_prob_fname, std::ios_base::trunc | std::ios_base::binary);
			time_delay_spin_nonflip_prob_table_->write(str);
			str.close();
		}
	}

	inp.open(time_delay_spin_nonflip_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << time_delay_spin_nonflip_fname << "\"" << std::endl;
		generate_time_delay_spin_nonflip_table();
		str.open(time_delay_spin_nonflip_fname, std::ios_base::trunc | std::ios_base::binary);
		time_delay_spin_nonflip_table_->write(str);
		str.close();
	} else {
		time_delay_spin_nonflip_table_->read(inp);
		inp.close();
		if (time_delay_spin_nonflip_table_->is_empty()) {
			std::cout << "Failed to load \"" << time_delay_spin_nonflip_fname << "\"" << std::endl;
			generate_time_delay_spin_nonflip_table();
			str.open(time_delay_spin_nonflip_fname, std::ios_base::trunc | std::ios_base::binary);
			time_delay_spin_nonflip_table_->write(str);
			str.close();
		}
	}

	inp.open(time_delay_spin_flip_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << time_delay_spin_flip_fname << "\"" << std::endl;
		generate_time_delay_spin_flip_table();
		str.open(time_delay_spin_flip_fname, std::ios_base::trunc | std::ios_base::binary);
		time_delay_spin_flip_table_->write(str);
		str.close();
	} else {
		time_delay_spin_flip_table_->read(inp);
		inp.close();
		if (time_delay_spin_flip_table_->is_empty()) {
			std::cout << "Failed to load \"" << time_delay_spin_flip_fname << "\"" << std::endl;
			generate_time_delay_spin_flip_table();
			str.open(time_delay_spin_flip_fname, std::ios_base::trunc | std::ios_base::binary);
			time_delay_spin_flip_table_->write(str);
			str.close();
		}
	}

	inp.open(total_resonance_NBrS_spectrum_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		std::cout << "Failed to load \"" << total_resonance_NBrS_spectrum_fname << "\"" << std::endl;
		generate_ResNBrS_spectrum_table();
		str.open(total_resonance_NBrS_spectrum_fname, std::ios_base::trunc);
		total_resonance_NBrS_spectrum_.write(str, "W[eV]\tFprob");
		str.close();
	} else {
		total_resonance_NBrS_spectrum_.read(inp);
		inp.close();
		if (total_resonance_NBrS_spectrum_.size()==0) {
			std::cout << "Failed to load \"" << total_resonance_NBrS_spectrum_fname << "\"" << std::endl;
			generate_ResNBrS_spectrum_table();
			str.open(total_resonance_NBrS_spectrum_fname, std::ios_base::trunc);
			total_resonance_NBrS_spectrum_.write(str, "W[eV]\tFprob");
			str.close();
		}
	}

	std::cout<<"Finished loading Ar data tables"<<std::endl;
}

void ArDataTables::DeleteData(void)
{
	delete integral_table_;
	delete theta_table_;
	delete time_delay_spin_flip_table_;
	delete time_delay_spin_nonflip_table_;
	delete time_delay_spin_nonflip_prob_table_;
	integral_table_ = NULL;
	theta_table_ = NULL;
	time_delay_spin_flip_table_ = NULL;
	time_delay_spin_nonflip_table_ = NULL;
	time_delay_spin_nonflip_prob_table_ = NULL;
}

void ArDataTables::generate_total_cross_elastic_table(void)
{
	std::cout << "Calculating total elastic cross section..." << std::endl;
	total_cross_elastic_.clear();
	total_cross_elastic_.setOrder(1); //interpolation with 1st order polynomial using 2 points
	total_cross_elastic_.setNused(2);

	int err;
	double E, cross;
	EnergyScanner EnRange(EnergyScanner::ElasticXS);
	while (true) {
		E = EnRange.Next(err);
		if (0 != err)
			break;
		cross = ArAllData_.argon_cross_elastic(E);
		total_cross_elastic_.push_back(E, cross);
	}
}

void ArDataTables::generate_integral_table(void) //uses tabulated total cross section
{
	std::cout<<"Generating cross section integrals..."<<std::endl;
	// v1.x-v3.x & v7.x
	ColoredRange energy_Y_range = ColoredInterval (0, 0.01, 1e-4) + ColoredInterval (0, EN_MAXIMUM_, 1e-3) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80);	//fine area

	/* v4.x
	ColoredRange energy_Y_range = ColoredInterval (0, 0.005, 1e-5) + ColoredInterval (0.005, 1, 1e-4) + ColoredInterval (0, EN_MAXIMUM_, 1e-3) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80);	//fine area
	*/
	// v5.x-v6.x
	/*ColoredRange energy_Y_range = ColoredInterval (0, 0.1, 1e-3) + ColoredInterval (0.1, 1, 2e-3) + ColoredInterval (0, EN_MAXIMUM_, 2.5e-3) +
			ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/2) +	//coarse area
			ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/2) +	//coarse area
			ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/20) + 	//fine area
			ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/20);	//fine area
	*/
	for (long int Ey_i=0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i!=Ey_ind_end_;++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		//v1.x-v4.x & v7.x
		ColoredRange energy_range = ColoredInterval (0, 0.1, 2e-4) + ColoredInterval (0.1, 1, 8e-4) +
					ColoredInterval (1, 10, 1e-2) + ColoredInterval (10, EN_MAXIMUM_, 0.02) +
					ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/5) +	//coarse area
					ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/5) +	//coarse area
					ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) + 	//fine area
					ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80) +	//fine area
					ColoredInterval (11.5, EN_MAXIMUM_, 0.003);
		//^dE defined by characteristics of cross section
		//dE defiend by sqrt(E/E-Ey) - finer dE when E is closer to Ey
		energy_range = energy_range + ColoredInterval(Ey, 1.1*Ey, Ey/500.0) + ColoredInterval(1.1*Ey, 10*Ey, Ey/100.0) + ColoredInterval(10*Ey, 50*Ey, Ey/10.0);

		// v5.x-v6.x
		/*ColoredRange energy_range = ColoredInterval (0, 0.1, 1e-3) + ColoredInterval (0.1, 1, 2e-3) +
					ColoredInterval (1, 10, 2e-2) + ColoredInterval (10, EN_MAXIMUM_, 0.02) +
					ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/2) +	//coarse area
					ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/2) +	//coarse area
					ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/20) + 	//fine area
					ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/20) +	//fine area
					ColoredInterval (11.5, EN_MAXIMUM_, 0.003);
		energy_range = energy_range + ColoredInterval(Ey, 1.1*Ey, Ey/500.0) + ColoredInterval(1.1*Ey, 10*Ey, Ey/10.0) + ColoredInterval(10*Ey, 50*Ey, Ey/1.0);
		*/
		energy_range.Trim(Ey, EN_MAXIMUM_);
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
				Int += 0.5*(TotalCrossSection(E) + TotalCrossSection(Ey))*sqrt(E*(E - Ey));
			} else {
				Int+=TotalCrossSection(E)*sqrt(E/(E-Ey))*(E-E_prev);
			}
			E_prev = E;
			integral_table_->push(E, Ey, Int);
		}
	}
}

void ArDataTables::generate_theta_table(void) //uses tabulated total cross section
{
	std::cout<<"Generating theta probability function tables..."<<std::endl;
	ColoredRange energy_Y_range = ColoredInterval (0, 0.01, 1e-4) + ColoredInterval (0, EN_MAXIMUM_, 1e-3) +
				ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/10) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/10) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/200) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/200);		//fine area
	std::vector<double> diff_XS, F, thetas;
	diff_XS.resize(ANGLE_POINTS_, 0);
	thetas.resize(ANGLE_POINTS_, 0);
	F.resize(ANGLE_POINTS_, 0);
	for (long int Ey_i=0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i!=Ey_ind_end_; ++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		long double (ArAllData::*diff_cross)(long double, long double, int) = 0;
		diff_cross = &ArAllData::argon_cross_elastic_diff;
		if (0 == diff_cross) {
			std::cout<<"Error: no differential cross section function"<<std::endl;
			break;
		}
		for (std::size_t i=0, i_end_=diff_XS.size(); i!=i_end_; ++i) {
			thetas[i] = i*M_PI/(double)(i_end_-1);
			F[i]=0;
			diff_XS[i]=(ArAllData_.*diff_cross)(Ey, thetas[i], 0);
			if (i!=0) {
				F[i] = F[i-1] + 0.5*(diff_XS[i]+diff_XS[i-1])*(thetas[i]-thetas[i-1]);//Integral
			}
		}
		for (std::size_t i=0, i_end_=diff_XS.size(); i!=i_end_; ++i)
			F[i] /= F[i_end_-1];//normalize probability function
		for (std::size_t th_i=0, th_i_end_=thetas.size(); th_i!=th_i_end_;++th_i)
			theta_table_->push(thetas[th_i], Ey, F[th_i]);
	}
}

void ArDataTables::generate_time_delay_spin_flip_table(void)
{
	std::cout << "Generating spin flip time delay tables..." << std::endl;
	ColoredRange energy_Y_range =
		ColoredInterval(En_3o2_ - 30 * Width_3o2_, En_3o2_ + 30 * Width_3o2_, Width_3o2_ / 10) +	//coarse area
		ColoredInterval(En_3o2_ - 15 * Width_3o2_, En_3o2_ + 15 * Width_3o2_, Width_3o2_ / 200)+		//fine area
		ColoredInterval(En_1o2_ - 30 * Width_1o2_, En_1o2_ + 30 * Width_1o2_, Width_1o2_ / 10) +	//coarse area
		ColoredInterval(En_1o2_ - 15 * Width_1o2_, En_1o2_ + 15 * Width_1o2_, Width_1o2_ / 200);		//fine area
	energy_Y_range.Trim(0, EN_MAXIMUM_);
	for (long int Ey_i = 0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i != Ey_ind_end_; ++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		for (std::size_t th_i = 0, th_i_end_ = ANGLE_POINTS_; th_i != th_i_end_; ++th_i) {
			double theta = th_i*M_PI / (ANGLE_POINTS_ - 1);
			time_delay_spin_flip_table_->push(theta, Ey, ArAllData_.argon_delay_spin_flip(Ey, theta));
		}
	}
}

void ArDataTables::generate_time_delay_spin_nonflip_table(void)
{
	std::cout << "Generating spin nonflip time delay tables..." << std::endl;
	ColoredRange energy_Y_range =
		ColoredInterval(En_3o2_ - 30 * Width_3o2_, En_3o2_ + 30 * Width_3o2_, Width_3o2_ / 10) +	//coarse area
		ColoredInterval(En_3o2_ - 15 * Width_3o2_, En_3o2_ + 15 * Width_3o2_, Width_3o2_ / 200)+		//fine area
		ColoredInterval(En_1o2_ - 30 * Width_1o2_, En_1o2_ + 30 * Width_1o2_, Width_1o2_ / 10) +	//coarse area
		ColoredInterval(En_1o2_ - 15 * Width_1o2_, En_1o2_ + 15 * Width_1o2_, Width_1o2_ / 200);		//fine area
	energy_Y_range.Trim(0, EN_MAXIMUM_);
	for (long int Ey_i = 0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i != Ey_ind_end_; ++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		for (std::size_t th_i = 0, th_i_end_ = ANGLE_POINTS_; th_i != th_i_end_; ++th_i) {
			double theta = th_i*M_PI / (ANGLE_POINTS_ - 1);
			time_delay_spin_nonflip_table_->push(theta, Ey, ArAllData_.argon_delay_spin_nonflip(Ey, theta));
		}
	}
}

void ArDataTables::generate_time_delay_spin_nonflip_prob_table(void)
{
	std::cout << "Generating spin nonflip time delay probability tables..." << std::endl;
	ColoredRange energy_Y_range =
		ColoredInterval(En_3o2_ - 30 * Width_3o2_, En_3o2_ + 30 * Width_3o2_, Width_3o2_ / 10) +	//coarse area
		ColoredInterval(En_3o2_ - 15 * Width_3o2_, En_3o2_ + 15 * Width_3o2_, Width_3o2_ / 200)+		//fine area
		ColoredInterval(En_1o2_ - 30 * Width_1o2_, En_1o2_ + 30 * Width_1o2_, Width_1o2_ / 10) +	//coarse area
		ColoredInterval(En_1o2_ - 15 * Width_1o2_, En_1o2_ + 15 * Width_1o2_, Width_1o2_ / 200);		//fine area
	energy_Y_range.Trim(0, EN_MAXIMUM_);
	for (long int Ey_i = 0, Ey_ind_end_ = energy_Y_range.NumOfIndices(); Ey_i != Ey_ind_end_; ++Ey_i) {
		double Ey = energy_Y_range.Value(Ey_i);
		for (std::size_t th_i = 0, th_i_end_ = ANGLE_POINTS_; th_i != th_i_end_; ++th_i) {
			double theta = th_i*M_PI / (ANGLE_POINTS_ - 1);
			time_delay_spin_nonflip_prob_table_->push(theta, Ey, ArAllData_.argon_delay_spin_nonflip_prob(Ey, theta));
		}
	}
}

void ArDataTables::generate_ResNBrS_spectrum_table(void)
{
	std::cout << "Generating resonance NBrS spectrum probability tables..." << std::endl;
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
}

long double ArDataTables::TotalCrossSection (double E)
{
	long double XS_total = 0;
	XS_total = XS_elastic(E);
	for (int i = 0, end_ = ArAllData_.ArExper_.ionizations.size(); i != end_; ++i) {
		XS_total += ArAllData_.ArExper_.ionizations[i](E);
	}
	for (int i = 0, end_ = ArAllData_.ArExper_.excitations.size(); i != end_; ++i) {
		XS_total += ArAllData_.ArExper_.excitations[i](E);
	}
	XS_total += ArAllData_.argon_ResNBrS_XS(E);
	return XS_total;
}

long double ArDataTables::CrossSection (double E, short type)
{
	if (Event::Elastic == type) {
		return XS_elastic(E);
	}
	if (Event::ResNBrS == type) {
		return ArAllData_.argon_ResNBrS_XS(E);
	}
	if (type>=Event::Ionization) {
		short ID = type - Event::Ionization;
		InelasticProcess *p = ArAllData_.ArExper_.FindInelastic(ID);
		if (NULL!=p)
			return (*p)(E);
	}
	return 0;
}

long double ArDataTables::XS_elastic(double E)
{
	return total_cross_elastic_(E, E);
}

//No Ramsauer minimum
/*long double ArDataTables::XS_elastic(double E)
{
	double E1 = 0.3286;
	double XS1 = 0.349195;
	if (E<E1) {
		return E*XS1/E1; //linear increase from (0,0). No Ramsauer minimum.
	}
	return total_cross_elastic_(E, E);
}*/

double ArDataTables::generate_theta (double E, short type, double Rand)
{
#ifdef	ANGLE_UNIFORM_
	Rand = Rand*2.0 - 1.0;
	return std::acos(Rand);
#endif
	if (type>=Event::Ionization) {//Considered uniform.
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	if (Event::Elastic == type) {
		return theta_table_->find_E(E, Rand);
	}
	if (Event::ResNBrS == type) {//Considered uniform.
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	//CODE below is for untabulated case
	std::vector<double> diff_XS, F, thetas;
	diff_XS.resize(ANGLE_POINTS_, 0);
	thetas.resize(ANGLE_POINTS_, 0);
	F.resize(ANGLE_POINTS_, 0);
	long double (ArAllData::*diff_cross)(long double, long double, int) = 0;
	if (Event::Elastic == type) {
		diff_cross = &ArAllData::argon_cross_elastic_diff;
	}
	if (0 == diff_cross) {
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	for (std::size_t i=0, i_end_=diff_XS.size(); i!=i_end_; ++i) {
		thetas[i] = i*M_PI/(double)(i_end_-1);
		diff_XS[i]=(ArAllData_.*diff_cross)(E, thetas[i], 0);
		if (i!=0) {
			F[i] = F[i-1] + 0.5*(diff_XS[i]+diff_XS[i-1])*(thetas[i]-thetas[i-1]);//Integral
		}
	}
	for (std::size_t i=0, i_end_=diff_XS.size(); i!=i_end_; ++i)
		F[i] /= F[i_end_-1];//normalize probability function
	double theta = 0;
	for (std::size_t i=1, i_end_=diff_XS.size(); i!=i_end_; ++i) {
		if (Rand<=F[i]) {
			theta = thetas[i-1] + (Rand-F[i-1])*(thetas[i]-thetas[i-1])/(F[i]-F[i-1]);
			break;
		}
	}
	return theta;
}

double ArDataTables::generate_time_delay(double E, double theta, short type, double Rand)
{
	if (Event::Elastic != type) {
		return 0;
	}
	double P = (*time_delay_spin_nonflip_prob_table_)(theta, E);
	P = std::max(P, 0.0);
	if (Rand < P) {
		return (*time_delay_spin_nonflip_table_)(theta, E);
	} else {
		return (*time_delay_spin_flip_table_)(theta, E);
	}
	return 0;
}

double ArDataTables::generate_time_delay_untabulated(double E, double theta, short type, double Rand)
{
	if (Event::Elastic != type) {
		return 0;
	}
	double P = ArAllData_.argon_delay_spin_nonflip_prob(E, theta);
	if (Rand < P) {
		return ArAllData_.argon_delay_spin_nonflip(E, theta);
	} else {
		return ArAllData_.argon_delay_spin_flip(E, theta);
	}
	return 0;
}

double ArDataTables::generate_ResNBrS_En_loss(double E, double theta, double Rand)
{
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

void ArDataTables::setOrder(int order)
{
	total_cross_elastic_.setOrder(order);
}

void ArDataTables::setNused(int N)
{
	total_cross_elastic_.setNused(N);
}

int ArDataTables::getOrder(void)
{
	return total_cross_elastic_.getOrder();
}

int ArDataTables::getNused(void)
{
	return total_cross_elastic_.getNused();
}
