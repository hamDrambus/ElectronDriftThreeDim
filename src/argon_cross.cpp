#include "argon_cross.h"

ArExperimental ArExper;
ArDataTables ArTables;

EnergyScanner::EnergyScanner(ScanType type): i(0), type_(type)
{
	/*ElasticXS, Resonance_3o2_XS, Resonance_1o2_XS, ResonancesXS,
		Resonance_3o2_DiffXS, Resonance_1o2_DiffXS, ResonancesDiffXS,
		DiffXS,	InelasticXS, ElasticResXS, XSIntegral, PlotElastic,
		PlotResonance_3o2, PlotResonance_1o2, PlotResonances,
		PlotDiffXS, PlotInelastic, PlotElasticResXS, PlotAllXS*/
	switch (type_) {
	case (ElasticXS): {
		//from 1e-3 eV to 0.1 eV with step 5e-4 eV, etc.
		energy_range_ = ColoredInterval (XS_EL_EN_MINIMUM_, 0.1, 5e-4) + ColoredInterval (0.1, 1, 1e-3) +
				ColoredInterval (1, 10, 5e-2) + ColoredInterval (10, XS_EL_EN_MAXIMUM_, 0.1);
		break;
	}
	case (Resonance_3o2_XS): {
		energy_range_ = ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/2) + 	//coarse area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/30); 		//fine area
		break;
	}
	case (Resonance_1o2_XS): {
		energy_range_ = ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/2) + 	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/30); 		//fine area
		break;
	}
	case (ResonancesXS): {
		energy_range_ = ColoredInterval (En_1o2_ - 100*Width_1o2_, En_1o2_ + 100*Width_1o2_, Width_1o2_/2) + 	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, En_3o2_ + 100*Width_3o2_, Width_3o2_/2) + 			//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/30) + 		//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/30); 		//fine area
		break;
	}
	case (Resonance_3o2_DiffXS): {
		energy_range_ =  ColoredInterval (En_3o2_ - 20*Width_3o2_, En_3o2_ + 20*Width_3o2_, 4*Width_3o2_);
		break;
	}
	case (Resonance_1o2_DiffXS): {
		energy_range_ = ColoredInterval (En_1o2_ - 20*Width_1o2_, En_1o2_ + 20*Width_1o2_, 4*Width_1o2_);
		break;
	}
	case (ResonancesDiffXS): {
		energy_range_ = ColoredInterval (En_1o2_ - 20*Width_1o2_, En_1o2_ + 20*Width_1o2_, 4*Width_1o2_) + 	//coarse area
				ColoredInterval (En_3o2_ - 20*Width_3o2_, En_3o2_ + 20*Width_3o2_, 4*Width_3o2_); 			//coarse area
		break;
	}
	case (DiffXS): {
		energy_range_ = ColoredInterval (PHASES_EN_MINIMUM_, 0.1, 1e-3) +
				ColoredInterval (0.1, 1, 1e-2) +
				ColoredInterval (1, PHASES_EN_MAXIMUM_, 0.1);
		break;
	}
	case (InelasticXS): {
		energy_range_ = ColoredInterval (11.5, EN_MAXIMUM_, 0.01);
		break;
	}
	case (ElasticResXS): {
		energy_range_ = ColoredInterval (XS_EL_EN_MINIMUM_, 0.1, 5e-4) + ColoredInterval (0.1, 1, 1e-3) +
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
		energy_range_ = ColoredInterval (XS_EL_EN_MINIMUM_, 0.1, 3e-4) + ColoredInterval (0.1, 1, 9e-4) +
				ColoredInterval (1, 10, 3e-2) + ColoredInterval (10, XS_EL_EN_MAXIMUM_, 0.086);
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
	case (PlotResonances): {
		energy_range_ = ColoredInterval (En_1o2_ - 110*Width_1o2_, En_1o2_ + 110*Width_1o2_, Width_1o2_/3) + 	//coarse area
				ColoredInterval (En_3o2_ - 110*Width_3o2_, En_3o2_ + 110*Width_3o2_, Width_3o2_/3) + 			//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, En_1o2_ + 15*Width_1o2_, Width_1o2_/80) + 		//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, En_3o2_ + 15*Width_3o2_, Width_3o2_/80); 			//fine area
		break;
	}
	case (PlotDiffXS): {
		energy_range_ = ColoredInterval (1e-4, 0.1, 7e-4) +
				ColoredInterval (0.1, 1, 7e-3) +
				ColoredInterval (1, EN_MAXIMUM_, 0.086);
		break;
	}
	case (PlotInelasticXS): {
		energy_range_ = ColoredInterval (11.5, XS_EL_EN_MAXIMUM_, 0.007) +
				ColoredInterval (XS_EL_EN_MAXIMUM_, 100, 0.1);
		break;
	}
	case (PlotElasticResXS): {
		energy_range_ = ColoredInterval (1e-4, 0.1, 2e-4) + ColoredInterval (0.1, 1, 7e-4) +
				ColoredInterval (1, 10, 7e-3) + ColoredInterval (10, XS_EL_EN_MAXIMUM_, 0.016) +
				ColoredInterval (En_1o2_ - 110*Width_1o2_, std::min(XS_EL_EN_MAXIMUM_, En_1o2_ + 110*Width_1o2_), Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 110*Width_3o2_, std::min(XS_EL_EN_MAXIMUM_, En_3o2_ + 110*Width_3o2_), Width_3o2_/5) +	//coarse area
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


InelasticProcess::InelasticProcess(std::string name, unsigned int ID, double En, double F, std::vector<double> &Ens, std::vector<double> &XSs):
		name_(name), ID_(ID), En_threshold_(En), Oscillator_strength_(F), exp_XS_(Ens, XSs, 3, 4)
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
	return (Oscillator_strength_/(En_threshold_*beta2))*log(beta2*gamma2*e_mass_eVconst/(2.0*En_threshold_)-beta2)*BBconst*E/(E+En_threshold_+ArExper.E_Ionization);
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
			phase_shifts_.push_back(DataVector(3,4)); /*interpolation by 3rd order polynomial because there's data for 11eV*/
			phase_shifts_.back().use_rightmost(true);
		} else {
			phase_shifts_.push_back(DataVector(3,5)); /*fit by 3rd order polynomial TODO: tests, tests*/
			phase_shifts_.back().use_rightmost(true);
		}
	inp.open("data/McEachranArPhaseShifts.dat");
	while (!inp.eof()) {
		std::getline(inp, line);
		if (line.size()>=2) {
			if ((line[0]=='/')&&(line[1]=='/'))
				continue;
			if ((line[0]=='\t')&&(line[1]=='\t')) //negative phase shifts are ignored
				continue;
		}
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double k = std::stod(word);
		//now process string "double\tdouble\tdouble\t..."
		//\t\t means that the double is omitted
		std::vector<double> vals;
		std::vector<bool> is_d;//is double presented in the table
		vals.resize(9,0);
		is_d.resize(9, false);
		for (int l=0;l<9;++l) {
			if (line.empty())
				break;
			if (line[0]=='\t') { //value is omitted in the table
				line.erase(line.begin());
				continue;
			}
			//now line does not start with \t
			word = strtoken(line,"\t"); //removes one \t after value
			vals[l] = std::stod(word);
			is_d[l] = true;
		}
		for (int l=0;l<6;++l) {
			if (is_d[l])
				phase_shifts_[l].push(k, vals[l]);
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
				if (E_vals.size()!=XS_vals.size()){
					std::size_t new_sz = std::min(E_vals.size(),XS_vals.size());
					std::cout<<"Line "<<line_count<<" warning! "<<proc_name<<" values size mismatch. Trimming data to size: "<<new_sz<<std::endl;
					if (E_vals.size()>new_sz)
						E_vals.erase(E_vals.begin()+new_sz, E_vals.end());
					if (XS_vals.size()>new_sz)
						XS_vals.erase(XS_vals.begin()+new_sz, XS_vals.end());
				}
				to.push_back(InelasticProcess(proc_name, max_process_ID, En_thresh, F, E_vals, XS_vals));
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
		to.push_back(InelasticProcess(proc_name, max_process_ID, En_thresh, F, E_vals, XS_vals));
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
	return phase_shifts_.size()-1;
}

long double ArExperimental::phase_shift (long double k, unsigned int l)
{
	if (l>max_L(k))
		return 0;
	return phase_shifts_[l](k, k);
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
	total_cross_resonance_3o2_fname("data_derived/total_cross_section_resonance3o2.dat"),
	total_cross_resonance_1o2_fname("data_derived/total_cross_section_resonance1o2.dat"),
	total_cross_elastic_(1,2), //interpolation with 1st order polynomial
	total_cross_resonance_3o2_(1,2),
	total_cross_resonance_1o2_(1,2)
{
	std::cout<<"Constructing Ar data tables"<<std::endl;
	ensure_file(total_cross_elastic_fname);
	ensure_file(total_cross_resonance_3o2_fname);
	ensure_file(total_cross_resonance_1o2_fname);
	total_cross_resonance_3o2_.set_out_value(0);
	total_cross_resonance_1o2_.set_out_value(0);

	std::ifstream inp;
	std::ofstream str;
	int err;
	inp.open(total_cross_elastic_fname);
	read_data(inp, total_cross_elastic_); //Cross section is stored in 1e-20 m^2 in both files and this program
	inp.close();
	if (total_cross_elastic_.size()<total_cross_elastic_.getNused()) {
		std::cout<<"Calculating total elastic cross section..."<<std::endl;
		total_cross_elastic_.clear();
		str.open(total_cross_elastic_fname, std::ios_base::trunc);
		str<<"//E[eV]\tXS elastic [1e-20 m^2]"<<std::endl;
		double E, cross;
		EnergyScanner EnRange(EnergyScanner::ElasticXS);
		while (true) {
			E = EnRange.Next(err);
			if (0!=err)
				break;
			cross = argon_cross_elastic(E);
			str<<E<<"\t"<<cross<<std::endl;
			total_cross_elastic_.push_back(E, cross);
		}
		str.close();
	}

	inp.open(total_cross_resonance_3o2_fname);
	read_data(inp, total_cross_resonance_3o2_);
	inp.close();
	if (total_cross_resonance_3o2_.size()<total_cross_resonance_3o2_.getNused()) {
		std::cout<<"Calculating total resonance 3/2 cross section..."<<std::endl;
		total_cross_resonance_3o2_.clear();
		str.open(total_cross_resonance_3o2_fname, std::ios_base::trunc);
		str<<"//E[eV]\tXS resonance 3/2 [1e-20 m^2]"<<std::endl;
		double E, cross;
		EnergyScanner EnRange(EnergyScanner::Resonance_3o2_XS);
		while (true) {
			E = EnRange.Next(err);
			if (0!=err)
				break;
			cross = argon_cross_resonance_3o2(E);
			str<<E<<"\t"<<cross<<std::endl;
			total_cross_resonance_3o2_.push_back(E, cross);
		}
		str.close();
	}

	inp.open(total_cross_resonance_1o2_fname);
	read_data(inp, total_cross_resonance_1o2_);
	inp.close();
	if (total_cross_resonance_1o2_.size()<total_cross_resonance_1o2_.getNused()) {
		std::cout<<"Calculating total resonance 1/2 cross section..."<<std::endl;
		total_cross_resonance_1o2_.clear();
		str.open(total_cross_resonance_1o2_fname, std::ios_base::trunc);
		str<<"//E[eV]\tXS resonance 1/2 [1e-20 m^2]"<<std::endl;
		double E, cross;
		EnergyScanner EnRange(EnergyScanner::Resonance_1o2_XS);
		while (true) {
			E = EnRange.Next(err);
			if (0!=err)
				break;
			cross = argon_cross_resonance_1o2(E);
			str<<E<<"\t"<<cross<<std::endl;
			total_cross_resonance_1o2_.push_back(E, cross);
		}
		str.close();
	}

	std::cout<<"Finished constructing Ar data tables"<<std::endl;
}

long double ArDataTables::TotalCrossSection (double E)
{
	long double XS_total = 0;
	for (int i = Event::Elastic, end_ =  Event::Elastic + 2 + ArExper.max_process_ID; i!=end_; ++i)
		XS_total+=CrossSection(E, i);
	return XS_total;
}

long double ArDataTables::CrossSection (double E, short type)
{
	if (Event::Elastic == type) {
		return XS_elastic(E);
	}
	if (Event::Resonance_3o2 == type) {
		return XS_resonance_3o2(E);
	}
	if (Event::Resonance_1o2 == type) {
		return XS_resonance_1o2(E);
	}
	if (type>=Event::Ionization) {
		short ID = type - Event::Ionization;
		InelasticProcess *p = ArExper.FindInelastic(ID);
		if (NULL!=p)
			return (*p)(E);
	}
	return 0;
}

long double ArDataTables::XS_elastic(double E)
{
	//if (E<EN_MINIMUM) //to avoid troubles. See Kurokawa for the choice of value.
	//	E = EN_MINIMUM;
	if (E<XS_EL_EN_MINIMUM_) {
		return 7.491 + E*(total_cross_elastic_(XS_EL_EN_MINIMUM_, XS_EL_EN_MINIMUM_) - 7.491)/XS_EL_EN_MINIMUM_; //linear behavior
	}
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

long double ArDataTables::XS_resonance_3o2(double E)
{
	return total_cross_resonance_3o2_(E, E);
}
long double ArDataTables::XS_resonance_1o2(double E)
{
	return total_cross_resonance_1o2_(E, E);
}

double ArDataTables::generate_Theta (double E, short type, double Rand) //TODO: tabulate
{
	if (type>=Event::Ionization) {//Considered uniform.
		return Rand*M_PI;
	}
	std::vector<double> diff_XS, F, thetas;
	diff_XS.resize(ANGLE_POINTS_, 0);
	thetas.resize(ANGLE_POINTS_, 0);
	F.resize(ANGLE_POINTS_, 0);
	long double (*diff_cross)(long double, long double);
	if (Event::Elastic == type) {
		diff_cross = argon_cross_elastic_diff;
	}
	if (Event::Resonance_3o2 == type) {
		diff_cross = argon_cross_resonance_3o2_diff;
	}
	if (Event::Resonance_1o2 == type) {
		diff_cross = argon_cross_resonance_1o2_diff;
	}
	for (std::size_t i=0, i_end_=diff_XS.size(); i!=i_end_; ++i) {
		thetas[i] = i*M_PI/(double)(i_end_-1);
		diff_XS[i]=diff_cross(E, thetas[i]);
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

void ArDataTables::setOrder(int order)
{
	total_cross_elastic_.setOrder(order);
	total_cross_resonance_3o2_.setOrder(order);
	total_cross_resonance_1o2_.setOrder(order);
}

void ArDataTables::setNused(int N)
{
	total_cross_elastic_.setNused(N);
	total_cross_resonance_3o2_.setNused(N);
	total_cross_resonance_1o2_.setNused(N);
}

int ArDataTables::getOrder(void)
{
	return total_cross_elastic_.getOrder();
}

int ArDataTables::getNused(void)
{
	return total_cross_elastic_.getNused();
}

//k is in atomic units
void argon_phase_values_exp(long double k, unsigned int l, long double &tan, long double &sin, long double &cos)
{
	double angle = ArExper.phase_shift(k, l);
	tan = std::tan(angle);
	sin = std::sin(angle);
	cos = std::cos(angle);
}

//k is in atomic units
void argon_phase_values_MERT5(long double k, unsigned int l, long double &tan, long double &sin, long double &cos)
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
	if (0==l) {
		tan = -A*k*(1+4*alpha_d*k*k*log(k)/3)-M_PI*alpha_d*k*k/3 + D*pow(k,3) + F* pow(k,4);
		tan/=(1+G*pow(k,3));
	} else {
		unsigned int l2 = 2*l;
		long double al = M_PI/((l2+3)*(l2+1)*(l2-1));
		long double bl = M_PI*(15*pow(l2+1,4)-140*pow(l2+1,2)+128)/(pow((l2+3)*(l2+1)*(l2-1),3)*(l2+5)*(l2-3));
		long double cl = 3*al/((l2+5)*(l2-3));
		tan = al*alpha_d*k*k+(bl*pow(alpha_d,2)+cl*alpha_q)*pow(k,4);
		if (1==l)
			tan+=H*pow(k,5) - A1*pow(k,3);
	}
	sin = fabs(tan)*sqrt(1/(1+tan*tan));
	cos = (tan>0?1.0:-1.0)*sqrt(1/(1+tan*tan));
}

//E in eV
long double argon_cross_elastic_diff (long double E, long double theta) {
	//different formulas are used for E<0.24eV and E>0.24eV!
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double cross = 0;
	long double cos_th = cos(theta);
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_PHASES_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			cross+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*P1(cos_th, l)*P2(cos_th, f);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			//P1 and P2 because each of them has separate cache.
		}
	}
	cross*=2*M_PI/pow(k,2);
	return cross*a_bohr_SIconst*a_bohr_SIconst; //const is multiplied bt 1e10
}

long double argon_cross_elastic (long double E)
{
	if (E <THRESH_E_XS_) {
		long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
		// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
		long double cross = 0;
		unsigned int L_MAX = L_MAX_;
		for (unsigned int l=0; l<=L_MAX; ++l) {
			long double sin_phase_l = 0;
			long double cos_phase_l = 1;
			long double tan_l = 0;
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
			cross+=(2*l+1)*sin_phase_l*sin_phase_l;
		}
		cross*=4*M_PI/pow(k,2);
		return cross*a_bohr_SIconst*a_bohr_SIconst;
	} else {
		return ArExper.total_elastic_cross(a_h_bar_2e_m_e_SIconst*sqrt(E), a_h_bar_2e_m_e_SIconst*sqrt(E));
	}
	return 0;
}

//for testing only
long double argon_cross_elastic_from_phases (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	long double cross = 0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		cross+=(2*l+1)*sin_phase_l*sin_phase_l;
	}
	cross*=4*M_PI/pow(k,2);
	return cross*a_bohr_SIconst*a_bohr_SIconst;
}

long double argon_back_scatter_prob (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_PHASES_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			W+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl(l,f, -1, 0, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
		}
		cross+=(2*l+1)*sin_phase_l*sin_phase_l;
	}
	return W/(2*cross);
}
//energy loss and input are in eV
long double argon_TM_forward (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_PHASES_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			W+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl_transf(l,f, 0, 1, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl(l,f, 0, 1, 1e-5);
		}
	}
	return W/cross;
}

long double argon_TM_backward (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper.max_L(k);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		if (E<THRESH_E_PHASES_)
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
		else
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				if (E<THRESH_E_PHASES_)
					argon_phase_values_MERT5(k, f, tan_f, sin_phase_f, cos_phase_f);
				else
					argon_phase_values_exp(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			W+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl_transf(l,f, -1, 0, 1e-5);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			cross+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*Int_PlPl(l,f, -1, 0, 1e-5);
		}
	}
	return W/cross;
}

long double argon_cross_resonance_3o2_diff (long double E, long double theta)
{
	return 0;
}

long double argon_cross_resonance_1o2_diff (long double E, long double theta)
{
	return 0;
}

long double argon_cross_resonance_3o2 (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E);
	long double sin_phase_1 = 0;
	long double cos_phase_1 = 1;
	long double tan_1 = 0;
	argon_phase_values_exp(k, 1, tan_1, sin_phase_1, cos_phase_1);
	long double cot = 2*(E - En_3o2_ )/Width_3o2_;
	//resonance phase changes from 0 to -Pi when energy changes from -inf to +inf, hence the choice of signs
	long double sin_3o2 = - sqrt(1.0/(1+pow( cot , 2)));
	long double cos_3o2 = - cot * sin_3o2;
	long double cross =8*M_PI*(pow(sin_phase_1*cos_3o2 + cos_phase_1*sin_3o2, 2) - sin_phase_1*sin_phase_1);
	cross /=pow (k, 2);
	return cross * a_bohr_SIconst * a_bohr_SIconst;
}

long double argon_cross_resonance_1o2 (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E);
	long double sin_phase_1 = 0;
	long double cos_phase_1 = 1;
	long double tan_1 = 0;
	argon_phase_values_exp(k, 1, tan_1, sin_phase_1, cos_phase_1);
	long double cot = 2*(E - En_1o2_ )/Width_1o2_;
	long double sin_1o2 = - sqrt(1.0/(1+pow( cot , 2)));
	long double cos_1o2 = - cot * sin_1o2;
	long double cross =4*M_PI*(pow(sin_phase_1*cos_1o2 + cos_phase_1*sin_1o2, 2) - sin_phase_1*sin_phase_1);
	cross /=pow (k, 2);
	return cross * a_bohr_SIconst * a_bohr_SIconst;
}

long double argon_back_resonance_3o2_prob (long double E)
{
	return argon_back_scatter_prob(E);
}

long double argon_back_resonance_1o2_prob (long double E)
{
	return argon_back_scatter_prob(E);
}

long double argon_TM_forward_resonance_3o2 (long double E)
{
	return argon_TM_forward(E);
}

long double argon_TM_forward_resonance_1o2 (long double E)
{
	return argon_TM_forward(E);
}

long double argon_TM_backward_resonance_3o2 (long double E)
{
	return argon_TM_backward(E);;
}

long double argon_TM_backward_resonance_1o2 (long double E)
{
	return argon_TM_backward(E);;
}
