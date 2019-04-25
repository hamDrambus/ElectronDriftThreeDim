#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "global_definitions.h"
#include "PolynomialFit.h"

//All energies in the program are in eV.
struct PhysicalConstants {
	long double e_charge_SI; //in coulombs (SI)
	long double e_mass_SI; //in kg (SI)
	long double e_mass_eV; //in eV
	long double h_bar_SI; //in SI
	long double h_bar_eVs; //in eV*s
	long double a_bohr_SI; //in meters (SI) multiplied by e10 for XS to be in 1e-20 m2
	long double Ry_eV; //Rydberg constant
	long double boltzmann_SI; //SI
							  //TODO: calculate after loading other constants; Calculated only once in order to increase the performance (the effect is not tested)
	long double a_h_bar_2eM_e_SI; //in SI a_bohr*sqrt(2*Me*e)/h_bar
	long double light_speed_SI;

	double Ar_mass_eV;
	double Ar_primal_line_nm;

	unsigned int MERT5_Lmax;
	double MERT5_A;
	double MERT5_D;
	double MERT5_F;
	double MERT5_G;
	double MERT5_A1;
	double MERT5_H;
	double MERT5_alpha_d;
	double MERT5_alpha_q;
	double phases_En_maximum;
	double phases_En_minimum;
	double phases_En_threshold; //MERT 5 is applied only below this energy for differential cross-section
	enum ScatterModel : short { Normal = 0, Uniform = 1 } scattering_angle_model;//(former #define ANGLE_UNIFORM_ for v9.x)
	bool no_ramsauer_minimum;
	double ramsauer_minimum_En;
	double XS_el_at_0_En;
	double XS_el_En_maximum;
	double XS_el_En_minimum;
	double XS_el_En_thresold;
	double XS_el_En_smooth_width; //MERT5 and experimental totalXS have small discontinuity at XS_el_En_thresold eV.

	double En_3o2;
	double En_1o2;
	double Width_3o2;
	double Width_1o2;

	double resonance_NBrS_XS;
	boost::optional<double> resonance_En_loss;
};

struct RunParameters {
	double field; //in Td
	double drift_distance; //in m
	ULong_t seed;
	unsigned int n_electrons;
	std::string output_file;
	boost::optional<PDF_routine> Ec_spectrum;
};

struct ProgramConstants {
	double temperature; //in SI (Kelvins)
	double pressure; //in SI (Pascals)
	unsigned int thread_number;
	double maximal_energy;
	unsigned int angle_discretization;
	boost::optional<double> drift_distance_ignore_history;
	boost::optional<int> skip_history_rate;
	bool is_test_version; //TODO: specify testing modules in settings.xml

	std::string data_folder;
	std::string elastic_XS_fname;
	std::string elastic_XS_phaseshift_fname;
	std::string excitation_XS_fname;
	std::string ionization_XS_fname;

	std::string output_fname_pattern;
	boost::optional<std::string> Ec_spectrum_fname_pattern;
	std::string tabulated_data_folder;

	double def_drift_distance;
	unsigned int def_n_electrons;
	ULong_t def_seed;
	enum GeneratorClass : short { NONE = 0, TRand1 = 1, TRand2 = 2, TRand3 = 3 } random_generator;

	std::map<std::string, bool> recorded_values;

	std::vector<RunParameters> run_specifics;
};

class Settings {
protected:
	bool is_valid_;
	PhysicalConstants phys_const_;
	ProgramConstants prog_const_;
public:
	bool isValid(void) const;
	Settings();
	Settings(std::string fname);
	bool Load(std::string fname);
	//bool Save(std::string fname) const;
	const PhysicalConstants* PhysConsts(void) const;
	const ProgramConstants* ProgConsts(void) const;
};

extern Settings gSettings;

#define En_1o2_ gSettings.PhysConsts()->En_1o2
#define En_3o2_ gSettings.PhysConsts()->En_3o2
#define Width_1o2_ gSettings.PhysConsts()->Width_1o2
#define Width_3o2_ gSettings.PhysConsts()->Width_3o2

#endif //SETTINGS_H_
