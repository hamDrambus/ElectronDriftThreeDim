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
	case (DiffXS): { //not used
		energy_range_ = ColoredInterval (PHASES_EN_MINIMUM_, 0.1, 1e-3) +
				ColoredInterval (0.1, 1, 1e-2) + ColoredInterval (1, PHASES_EN_MAXIMUM_, 0.1);
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
	case (PlotDiffXS): {  //used for plotting test of differential to total cross sections
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

ArDataTables::ArDataTables(FunctionTable * int_table, FunctionTable * th_table):
	total_cross_elastic_fname("data_derived/total_cross_section_elastic.dat"),
	integral_table_fname("data_derived/cross_integrals.dat"),
	theta_table_fname("data_derived/theta_probabilities.dat"),
	total_cross_elastic_(1,2), //interpolation with 1st order polynomial
	total_cross_resonance_3o2_(1,2),
	total_cross_resonance_1o2_(1,2),
	integral_table_(int_table),
	theta_table_(th_table)
{
	std::cout<<"Constructing Ar data tables"<<std::endl;
	ensure_file(total_cross_elastic_fname);
	ensure_file(integral_table_fname);
	ensure_file(theta_table_fname);
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
			cross = ArAllData_.argon_cross_elastic(E);
			str<<E<<"\t"<<cross<<std::endl;
			total_cross_elastic_.push_back(E, cross);
		}
		str.close();
	}

	inp.open(integral_table_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		generate_integral_table();
		str.open(integral_table_fname, std::ios_base::trunc|std::ios_base::binary);
		integral_table_->write(str);
		str.close();
	} else {
		integral_table_->read(inp);
		inp.close();
		if (integral_table_->is_empty()) {
			generate_integral_table();
			str.open(integral_table_fname, std::ios_base::trunc|std::ios_base::binary);
			integral_table_->write(str);
			str.close();
		}
	}

	inp.open(theta_table_fname, std::ios_base::binary);
	if (!inp.is_open()) {
		generate_theta_table();
		str.open(theta_table_fname, std::ios_base::trunc|std::ios_base::binary);
		theta_table_->write(str);
		str.close();
	} else {
		theta_table_->read(inp);
		inp.close();
		if (theta_table_->is_empty()) {
			generate_theta_table();
			str.open(theta_table_fname, std::ios_base::trunc|std::ios_base::binary);
			theta_table_->write(str);
			str.close();
		}
	}

	std::cout<<"Finished constructing Ar data tables"<<std::endl;
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
				ColoredInterval (En_1o2_ - 100*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 100*Width_1o2_), Width_1o2_/5) +	//coarse area
				ColoredInterval (En_3o2_ - 100*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 100*Width_3o2_), Width_3o2_/5) +	//coarse area
				ColoredInterval (En_1o2_ - 15*Width_1o2_, std::min(EN_MAXIMUM_, En_1o2_ + 15*Width_1o2_), Width_1o2_/80) + 	//fine area
				ColoredInterval (En_3o2_ - 15*Width_3o2_, std::min(EN_MAXIMUM_, En_3o2_ + 15*Width_3o2_), Width_3o2_/80);	//fine area
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

long double ArDataTables::TotalCrossSection (double E)
{
	long double XS_total = 0;
	XS_total = XS_elastic(E) + XS_resonance_3o2(E) + XS_resonance_1o2(E);
	for (int i = 0, end_ = ArAllData_.ArExper_.ionizations.size(); i != end_; ++i) {
		XS_total += ArAllData_.ArExper_.ionizations[i](E);
	}
	for (int i = 0, end_ = ArAllData_.ArExper_.excitations.size(); i != end_; ++i) {
		XS_total += ArAllData_.ArExper_.excitations[i](E);
	}
	return XS_total;
}

long double ArDataTables::CrossSection (double E, short type)
{
	if (Event::Elastic == type) {
		return XS_elastic(E);
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

long double ArDataTables::XS_resonance_3o2(double E)
{
	return total_cross_resonance_3o2_(E, E);
}
long double ArDataTables::XS_resonance_1o2(double E)
{
	return total_cross_resonance_1o2_(E, E);
}

double ArDataTables::generate_Theta (double E, short type, double Rand) //DONE: tabulate
{
#ifdef	ANGLE_UNIFORM_
	Rand = Rand*2.0 - 1.0;
	return std::acos(Rand);
#endif
	/*if (Rand<0.4)
		return M_PI;
	return 0;*/
	if (type>=Event::Ionization) {//Considered uniform.
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	if (Event::Elastic == type) {
		return theta_table_->find_E(E, Rand);
	}
	std::vector<double> diff_XS, F, thetas;
	diff_XS.resize(ANGLE_POINTS_, 0);
	thetas.resize(ANGLE_POINTS_, 0);
	F.resize(ANGLE_POINTS_, 0);
	long double (ArAllData::*diff_cross)(long double, long double, int) = 0;
	if (Event::Elastic == type) {
		diff_cross = &ArAllData::argon_cross_elastic_diff;
	}
	if (0 == diff_cross)
		return Rand*M_PI;
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

ArAllData::ArAllData (void)
{}
//k is in atomic units
void ArAllData::argon_phase_values_exp(long double k, unsigned int l, long double &tan, long double &sin, long double &cos)
{
	double angle = ArExper_.phase_shift(k, l);
	tan = std::tan(angle);
	sin = std::sin(angle);
	cos = std::cos(angle);
}

//k is in atomic units
void ArAllData::argon_phase_values_MERT5(long double k, unsigned int l, long double &tan, long double &sin, long double &cos)
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

//mode = 0 or unexpected - standard formula used in simulation, called by default.
//mode = 1 - calculate MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_PHASES_)
//mode = 2 - calculate using extrapolation of experimental phase shifts (normally used only between THRESH_E_PHASES_ and EN_MAXIMUM_)
//E in eV
long double ArAllData::argon_cross_elastic_diff (long double E, long double theta, int mode) {
	//different formulas are used for E<0.24eV and E>0.24eV!
	void (ArAllData::*phase_values)(long double , unsigned int , long double &, long double &, long double &) = 0;
	unsigned int L_MAX = 0;
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
	switch (mode) {
	case 1: {
		L_MAX = L_MAX_;
		phase_values=&ArAllData::argon_phase_values_MERT5;
		break;
	}
	case 2: {
		L_MAX = ArExper_.max_L(k);
		phase_values=&ArAllData::argon_phase_values_exp;
		break;
	}
	default: {
		if (PHASES_EN_MINIMUM_>E)
			E = PHASES_EN_MINIMUM_;
		k = a_h_bar_2e_m_e_SIconst*sqrt(E);
		if (E<THRESH_E_PHASES_) {
			L_MAX = L_MAX_;
			phase_values=&ArAllData::argon_phase_values_MERT5;
		} else {
			L_MAX = ArExper_.max_L(k);
			phase_values=&ArAllData::argon_phase_values_exp;
		}
		break;
	}
	}
	long double cross = 0;
	LegendrePolynom P1, P2;
	long double cos_th = cos(theta);
	for (unsigned int l=0; l<=L_MAX; ++l) {
		long double sin_phase_l = 0;
		long double cos_phase_l = 1;
		long double tan_l = 0;
		((*this).*phase_values)(k, l, tan_l, sin_phase_l, cos_phase_l);
		for (unsigned int f=l; f<=L_MAX; ++f) {
			long double sin_phase_f = sin_phase_l;
			long double cos_phase_f = cos_phase_l;
			long double tan_f = tan_l;
			if (l!=f) {
				((*this).*phase_values)(k, f, tan_f, sin_phase_f, cos_phase_f);
			}
			long double cos_l_f = cos_phase_l*cos_phase_f + sin_phase_l*sin_phase_f;
			cross+=((l==f)?1.0:2.0)*(2*l+1)*(2*f+1)*sin_phase_l*sin_phase_f*cos_l_f*P1(cos_th, l)*P2(cos_th, f);
			//((l==f)?1.0:2.0) because sum by f starts not from 0 but from l
			//P1 and P2 because each of them has separate cache.
		}
	}
	cross*=2*M_PI/pow(k,2);
	return cross*a_bohr_SIconst*a_bohr_SIconst; //const is multiplied by 1e10
}

//mode = 0 or unexpected - standard formula used in simulation, called by default. Using smoothing between MERT5 and experimental XS
//mode = 1 - calculate using MERT5 phases (normally only between EN_MINIMUM_ and THRESH_E_XS_)
//mode = 2 - calculate using extrapolation of experimental total XS (normally used only between THRESH_E_XS_ and EN_MAXIMUM_) + resonances
//mode = 3 - calculate using experimental phases (for testing only) + resonances
//mode = 4 - same as 0, but no smoothing
//below XS_EL_EN_MINIMUM_ linear extrapolation to XS(0)=7.491 (1e-20 m^2) (in default mode)
long double ArAllData::argon_cross_elastic (long double E, int mode) //Tabulation of this function must be done carefully. Do not forget near 0 case.
{
	long double cross = 0;
	switch (mode) {
	case 1:{
		long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
		// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
		unsigned int L_MAX = L_MAX_;
		for (unsigned int l=0; l<=L_MAX; ++l) {
			long double sin_phase_l = 0;
			long double cos_phase_l = 1;
			long double tan_l = 0;
			argon_phase_values_MERT5(k, l, tan_l, sin_phase_l, cos_phase_l);
			cross+=(2*l+1)*sin_phase_l*sin_phase_l;
		}
		cross*=4*M_PI/pow(k,2);
		cross*=a_bohr_SIconst*a_bohr_SIconst;
		break;
	}
	case 2: {
		cross = ArExper_.total_elastic_cross(a_h_bar_2e_m_e_SIconst*sqrt(E), a_h_bar_2e_m_e_SIconst*sqrt(E));
		cross+=argon_cross_resonance_1o2(E);
		cross+=argon_cross_resonance_3o2(E);
		break;
	}
	case 3:{
		long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units
		unsigned int L_MAX = ArExper_.max_L(k);;
		for (unsigned int l=0; l<=L_MAX; ++l) {
			long double sin_phase_l = 0;
			long double cos_phase_l = 1;
			long double tan_l = 0;
			argon_phase_values_exp(k, l, tan_l, sin_phase_l, cos_phase_l);
			cross+=(2*l+1)*sin_phase_l*sin_phase_l;
		}
		cross*=4*M_PI/pow(k,2);
		cross*=a_bohr_SIconst*a_bohr_SIconst;
		cross+=argon_cross_resonance_1o2(E);
		cross+=argon_cross_resonance_3o2(E);
		break;
	}
	case 4: {
		if (E<XS_EL_EN_MINIMUM_) {
			cross = 7.491 + E*(argon_cross_elastic(XS_EL_EN_MINIMUM_, 1) - 7.491)/XS_EL_EN_MINIMUM_; //linear behavior. MERT5 is used at XS_EL_EN_MINIMUM_
		} else {
			if (E <THRESH_E_XS_) {
				cross  = argon_cross_elastic(E, 1); //MERT5
			} else {
				cross  = argon_cross_elastic(E, 2); //Experiment
			}
		}
		break;
	}
	default: {
#ifndef D_EN_SMOOTH_
		cross = argon_cross_elastic(E, 4);
#else
		if (E<XS_EL_EN_MINIMUM_) {
			cross = 7.491 + E*(argon_cross_elastic(XS_EL_EN_MINIMUM_, 1) - 7.491)/XS_EL_EN_MINIMUM_; //linear behavior. MERT5 is used at XS_EL_EN_MINIMUM_
		} else {
			double E_l = THRESH_E_XS_ - D_EN_SMOOTH_;
			double E_r = THRESH_E_XS_ + D_EN_SMOOTH_;
			if (E <E_l) {
				cross  = argon_cross_elastic(E, 1);//MERT5
				break;
			}
			if (E >E_r){
				cross  = argon_cross_elastic(E, 2);//Experiment
				break;
			}
			cross = ((E_r - E)*argon_cross_elastic(E, 1) + (E-E_l)*argon_cross_elastic(E, 2))/(E_r-E_l); //linear smoothing between MERT5 and experiment
		}
#endif
		break;
	}
	}
	return cross;
}

long double ArAllData::argon_back_scatter_prob (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper_.max_L(k);
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
long double ArAllData::argon_TM_forward (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper_.max_L(k);
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

long double ArAllData::argon_TM_backward (long double E)
{
	long double k = a_h_bar_2e_m_e_SIconst*sqrt(E); //recalculation from energy to atomic units is following:
	// k[atomic] = a_bohr * sqrt(2 * m_electron[SI] * q_electron[SI] * E[eV]) / h_bar(plank const)[SI].
	LegendrePolynom P1, P2;
	long double W = 0;
	long double cross =0;
	unsigned int L_MAX = (E<THRESH_E_PHASES_) ? L_MAX_ : ArExper_.max_L(k);
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

long double ArAllData::argon_cross_resonance_3o2_diff (long double E, long double theta, int mode)
{
	return argon_cross_elastic_diff(E, theta, mode);
}

long double ArAllData::argon_cross_resonance_1o2_diff (long double E, long double theta, int mode)
{
	return argon_cross_elastic_diff(E, theta, mode);
}

long double ArAllData::argon_cross_resonance_3o2 (long double E)
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

long double ArAllData::argon_cross_resonance_1o2 (long double E)
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

long double ArAllData::argon_back_resonance_3o2_prob (long double E)
{
	return argon_back_scatter_prob(E);
}

long double ArAllData::argon_back_resonance_1o2_prob (long double E)
{
	return argon_back_scatter_prob(E);
}

long double ArAllData::argon_TM_forward_resonance_3o2 (long double E)
{
	return argon_TM_forward(E);
}

long double ArAllData::argon_TM_forward_resonance_1o2 (long double E)
{
	return argon_TM_forward(E);
}

long double ArAllData::argon_TM_backward_resonance_3o2 (long double E)
{
	return argon_TM_backward(E);
}

long double ArAllData::argon_TM_backward_resonance_1o2 (long double E)
{
	return argon_TM_backward(E);
}
