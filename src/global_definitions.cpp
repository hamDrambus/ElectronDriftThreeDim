#include "global_definitions.h"

Settings gSettings;

std::string strtoken(std::string &in, std::string break_symbs)
{
	std::string out;
	while (!in.empty())
	{
		char a = in.front();
		in.erase(in.begin());
		bool break_ = false;
		for (auto h = break_symbs.begin(); h != break_symbs.end(); ++h)
			if (a == *h) {
				break_ = true;
				break;
			}
		if ((break_) && (out.empty()))
			continue;
		if (break_)
			return out;
		out.push_back(a);
	}
	return out;
}

void ensure_file(std::string fname)
{
	std::string folder = fname;
	while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
		folder.pop_back();
	if (!folder.empty())
		folder.pop_back();
	ensure_folder(folder);
}

void ensure_folder(std::string folder)
{
#if defined(__WIN32__)
	if (!folder.empty()) {
		DWORD ftyp = GetFileAttributesA(folder.c_str());
		if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES) {
			int code = system(("mkdir \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir error: " << GetLastError() << std::endl;
		}
	}
#else
	struct stat st;
	stat(folder.c_str(), &st);
	if (!S_ISDIR(st.st_mode)) {
		int code = system(("mkdir -p \"" + folder + "\"").c_str());
		if (code)
			std::cout << "mkdir -p error: " << code << std::endl;
	}
#endif //_WIN32__
}

char* c_str_cp (const std::string &str)
{
	std::size_t i_end_= str.size();
	char* ret = new char [i_end_+1];
	for (std::size_t i=0; i!=i_end_; ++i) {
		ret[i] = str[i];
	}
	ret[i_end_] = NULL;
	return ret;
}

Settings::Settings():is_valid_(false) {}

Settings::Settings(std::string fname):is_valid_(false)
{
	Load(fname);
}

bool Settings::Load(std::string fname)
{
	std::cout << "Settings::Load: Loading \"" << fname << "\"..." << std::endl;
	// Create an empty property tree object
	using boost::property_tree::ptree;
	using boost::property_tree::ptree_bad_data;
	using boost::property_tree::ptree_bad_path;
	using boost::property_tree::ptree_error;
	ptree pt;
	try {
		read_xml(fname, pt);
		ptree physics = pt.get_child("Settings.PhysicalConstants");
		phys_const_.e_charge_SI =	 physics.get<double>("electron_charge_SI");
		phys_const_.e_mass_SI =		 physics.get<double>("electron_mass_SI");
		phys_const_.e_mass_eV =		 physics.get<double>("electron_mass_eV");
		phys_const_.h_bar_SI =		 physics.get<double>("h_bar_SI");
		phys_const_.h_bar_eVs =		 physics.get<double>("h_bar_eVs");
		phys_const_.a_bohr_SI =		 physics.get<double>("a_bohr_SI"); //in 1e-10 m
		phys_const_.Ry_eV =			 physics.get<double>("Ry_energy_eV");
		phys_const_.boltzmann_SI =	 physics.get<double>("boltzmann_SI");
		phys_const_.light_speed_SI=	 physics.get<double>("light_speed_SI");
		phys_const_.a_h_bar_2eM_e_SI =
			phys_const_.a_bohr_SI*1e-10*sqrt(2* phys_const_.e_mass_SI * phys_const_.e_charge_SI)/ phys_const_.h_bar_SI;
		physics = physics.get_child("Argon");
		phys_const_.Ar_mass_eV =			 physics.get<double>("Ar_mass_eV");
		phys_const_.Ar_primal_line_nm =		 physics.get<double>("Ar_primal_line_nm");
		phys_const_.MERT5_Lmax =			 physics.get<unsigned int>("MERT5_Lmax");
		phys_const_.MERT5_A =				 physics.get<double>("MERT5_A");
		phys_const_.MERT5_D =				 physics.get<double>("MERT5_D");
		phys_const_.MERT5_F =				 physics.get<double>("MERT5_F");
		phys_const_.MERT5_G =				 physics.get<double>("MERT5_G");
		phys_const_.MERT5_A1 =				 physics.get<double>("MERT5_A1");
		phys_const_.MERT5_H =				 physics.get<double>("MERT5_H");
		phys_const_.MERT5_alpha_d =			 physics.get<double>("MERT5_alpha_d");
		phys_const_.MERT5_alpha_q =			 physics.get<double>("MERT5_alpha_q");
		
		phys_const_.phases_En_minimum =		 physics.get<double>("elastic_phases_En_min");
		phys_const_.phases_En_maximum =		 physics.get<double>("elastic_phases_En_max");
		phys_const_.phases_En_threshold =	 physics.get<double>("elastic_phases_En_threshold");
		phys_const_.XS_el_En_minimum =		 physics.get<double>("elastic_XS_En_min");
		phys_const_.XS_el_En_maximum =		 physics.get<double>("elastic_XS_En_max");
		phys_const_.XS_el_En_thresold =		 physics.get<double>("elastic_XS_En_threshold");
		phys_const_.XS_el_En_smooth_width =	 physics.get<double>("elastic_XS_En_smooth_width");
		phys_const_.no_ramsauer_minimum =	 physics.get<bool>("elastic_XS_no_Ramsauer_minimum", false);
		if (phys_const_.no_ramsauer_minimum)
			phys_const_.ramsauer_minimum_En = physics.get<double>("elastic_XS_Ramsauer_minimum_En");
		else
			phys_const_.ramsauer_minimum_En = physics.get<double>("elastic_XS_Ramsauer_minimum_En", 0);
		phys_const_.XS_el_at_0_En = physics.get<double>("elastic_XS_at_0_En");
		short model_candidate = physics.get<short>("scattering_angle_model", 0);
		try {
			phys_const_.scattering_angle_model = (PhysicalConstants::ScatterModel)model_candidate;
		} catch (std::exception &e) {
			phys_const_.scattering_angle_model = PhysicalConstants::ScatterModel::Normal;
		}

		phys_const_.En_3o2 = physics.get<double>("Feshbach_resonance_3o2_En");
		phys_const_.En_1o2 = physics.get<double>("Feshbach_resonance_1o2_En");
		phys_const_.Width_3o2 = physics.get<double>("Feshbach_resonance_3o2_Width");
		phys_const_.Width_1o2 = physics.get<double>("Feshbach_resonance_1o2_Width");
		phys_const_.resonance_NBrS_XS = physics.get<double>("Feshbach_resonance_NBrS_XS"); //in 1e-20 m^2
		phys_const_.resonance_En_loss = physics.get_optional<double>("Feshbach_resonance_NBrS_En_loss");
		
		ptree params = pt.get_child("Settings.ProgramConstants");
		prog_const_.is_test_version =	 params.get<bool>("is_test_version", false);
		prog_const_.thread_number =		 params.get<unsigned int>("thread_number");
		prog_const_.temperature =		 params.get<double>("temperature");
		prog_const_.pressure =			 params.get<double>("pressure");
		prog_const_.angle_discretization = params.get<unsigned int>("angle_discretization");
		prog_const_.maximal_energy =	 params.get<double>("maximal_energy");

		prog_const_.skip_history_rate = params.get_optional<int>("skip_history_rate");
		prog_const_.drift_distance_ignore_history = params.get_optional<double>("drift_distance_ignore_history");
		BOOST_FOREACH(ptree::value_type &w, params.get_child("RecordedValues")) {
			std::string key = w.first;
			bool value = params.get<bool>("RecordedValues." + key, false);
			if (value) {
				prog_const_.recorded_values[key] = value;
			}
		}

		prog_const_.data_folder = params.get<std::string>("data_location", "data");
		std::string prefix = prog_const_.data_folder +
			(prog_const_.data_folder.empty() ? "" : ((prog_const_.data_folder.back() == '/') ? "" : "/"));
		try {
			prog_const_.elastic_XS_phaseshift_fname = params.get<std::string>("DataFiles.elastic_XS_phaseshifts_data");
		} catch(ptree_error& e) {
			prog_const_.elastic_XS_phaseshift_fname = prefix + "McEachranArPhaseShifts.dat";
		}
		try {
			prog_const_.elastic_XS_fname = params.get<std::string>("DataFiles.elastic_XS_data");
		} catch (ptree_error& e) {
			prog_const_.elastic_XS_fname = prefix + "ArScatteringCross.dat";
		}
		try {
			prog_const_.ionization_XS_fname = params.get<std::string>("DataFiles.ionization_XS_data");
		} catch (ptree_error& e) {
			prog_const_.ionization_XS_fname = prefix + "ArIonizations_Magboltz.dat";
		}
		try {
			prog_const_.excitation_XS_fname = params.get<std::string>("DataFiles.excitation_XS_data");
		} catch (ptree_error& e) {
			prog_const_.excitation_XS_fname = prefix + "ArExcitations_Magboltz.dat";
		}
		prog_const_.tabulated_data_folder = params.get<std::string>("cache_data_folder");
		prog_const_.tabulated_data_folder = prog_const_.tabulated_data_folder +
			(prog_const_.tabulated_data_folder.empty() ? "" : ((prog_const_.tabulated_data_folder.back() == '/') ? "" : "/"));
		prog_const_.output_fname_pattern = params.get<std::string>("output_file");

		prog_const_.def_drift_distance = params.get<double>("drift_distance");
		prog_const_.def_n_electrons = params.get<unsigned int>("n_electrons");
		prog_const_.def_seed = params.get<ULong_t>("random_seed");
		std::string Ldrift = params.get<std::string>("drift_distance");
		std::string Ne = params.get<std::string>("n_electrons");
		std::string seed = params.get<std::string>("random_seed");

		std::string generator = params.get<std::string>("random_generator");
		prog_const_.random_generator = ProgramConstants::GeneratorClass::NONE;
		if (generator=="TRandom1")
			prog_const_.random_generator = ProgramConstants::GeneratorClass::TRand1;
		if (generator=="TRandom2")
			prog_const_.random_generator = ProgramConstants::GeneratorClass::TRand2;
		if (generator=="TRandom3")
			prog_const_.random_generator = ProgramConstants::GeneratorClass::TRand3;
		if (ProgramConstants::GeneratorClass::NONE==prog_const_.random_generator) {
			BOOST_PROPERTY_TREE_THROW(ptree_bad_data(
					std::string("conversion of type \"") + typeid(std::string).name() +
					"\" to ProgramConstants::GeneratorClass data failed", boost::any()));
		}
		BOOST_FOREACH(ptree::value_type &w,
			params.get_child("Runs"))
		{
			if (w.first == "Run") {
				RunParameters run_pars;
				ptree run_info = w.second;
				run_pars.field =			 run_info.get<double>("Td");
				run_pars.drift_distance =	 run_info.get<double>("drift_distance", prog_const_.def_drift_distance);
				run_pars.seed =				 run_info.get<ULong_t>("random_seed", prog_const_.def_seed);
				run_pars.n_electrons =		 run_info.get<unsigned int>("n_electrons", prog_const_.def_n_electrons);
				std::string Td_str =	run_info.get<std::string>("Td");
				std::string L_str =		(prog_const_.def_drift_distance == run_pars.drift_distance ? Ldrift : run_info.get<std::string>("drift_distance"));
				std::string Ne_str =	(prog_const_.def_n_electrons == run_pars.n_electrons ? Ne : run_info.get<std::string>("n_electrons"));
				std::string seed_str =	(prog_const_.def_seed == run_pars.seed ? seed : run_info.get<std::string>("random_seed"));
				std::string Td_key =	"($Td)";
				std::string L_key =		"($drift_distance)";
				std::string Ne_key =	"($n_electrons)";
				std::string seed_key =	"($random_seed)";
				std::string out_fname = prog_const_.output_fname_pattern;
				while (true) { //output filename in settings can contain references to run-sprecific parameters (usually to field in Td)
					std::size_t rep = out_fname.find(Td_key);
					if (std::string::npos != rep) {
						out_fname.replace(rep, Td_key.size(), Td_str);
						continue;
					}
					rep = out_fname.find(L_key);
					if (std::string::npos != rep) {
						out_fname.replace(rep, L_key.size(), L_str);
						continue;
					}
					rep = out_fname.find(Ne_key);
					if (std::string::npos != rep) {
						out_fname.replace(rep, Ne_key.size(), Ne_str);
						continue;
					}
					rep = out_fname.find(seed_key);
					if (std::string::npos != rep) {
						out_fname.replace(rep, seed_key.size(), seed_str);
						continue;
					}
					break;
				}
				run_pars.output_file = out_fname;
				prog_const_.run_specifics.push_back(run_pars);
			} else {
				std::cout << "Settings::Load: Unknown key \"" << w.first << "\" in \"Runs\". Ignored." << std::endl;
			}
		}
	} catch (ptree_bad_path& e) {
		std::cout << "Settings::Load: ptree_bad_path exception:" << std::endl;
		std::cout << e.what() << std::endl;
		goto fail_load;
	} catch (ptree_bad_data& e) {
		std::cout << "Settings::Load: ptree_bad_data exception:" << std::endl;
		std::cout << e.what() << std::endl;
		goto fail_load;
	} catch (ptree_error& e) {
		std::cout << "Settings::Load: ptree_error exception:" << std::endl;
		std::cout << e.what() << std::endl;
		goto fail_load;
	} catch (std::exception& e) {
		std::cout << "Settings::Load: std::exception:" << std::endl;
		std::cout << e.what() << std::endl;
		goto fail_load;
	}
	is_valid_ = true;
	std::cout << "Settings::Load: Successfully loaded \"" << fname << "\"" << std::endl;
	return true;
fail_load:
	is_valid_ = false;
	std::cout << "Settings::Load: Failed to load settings \"" << fname <<"\"" << std::endl;
	return false;
}

const PhysicalConstants* Settings::PhysConsts(void) const
{
	return &phys_const_;
}

const ProgramConstants* Settings::ProgConsts(void) const
{
	return &prog_const_;
}

bool Settings::isValid(void) const
{	return is_valid_;}
