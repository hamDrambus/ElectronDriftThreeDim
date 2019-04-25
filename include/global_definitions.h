#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#if defined(__WIN32__)
#define NOMINMAX
#include "Windows4Root.h"
#include <direct.h>
#include <kbd.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif
#include <math.h>
#include <ctgmath>

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
	double Ar_mass_eV; //eV
	//TODO: calculate after loading other constants; Calculated only once in order to increase the performance (the effect is not tested)
	long double a_h_bar_2eM_e_SI; //in SI a_bohr*sqrt(2*Me*e)/h_bar

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
	enum ScatterModel : short {Normal = 0, Uniform = 1} scattering_angle_model;//(former #define ANGLE_UNIFORM_ for v9.x)
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
	unsigned int seed;
	unsigned int n_electrons;
	std::string output_file;
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
	std::string tabulated_data_folder;
	
	double def_drift_distance;
	unsigned int def_n_electrons;
	unsigned int def_seed;

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

#if defined(__WIN32__)
#define INVOKE_GNUPLOT(a) system(("start \"\" \"%GNUPLOT%\\gnuplot.exe\" --persist \"" + a + "\"").c_str())
#else
#define INVOKE_GNUPLOT(a) system(("gnome-terminal -- bash -c \"gnuplot \"" + a +"\"\"").c_str());
#endif //__WIN32__

std::string strtoken(std::string &in, std::string break_symbs);
void ensure_file(std::string fname); //makes sure file can be created later on
void ensure_folder(std::string folder);
char* c_str_cp (const std::string &str);

struct Event
{
	double En_start;
	double En_collision;
	double En_finish;
	double En_avr;

	double theta_start;
	double theta_collision;
	double theta_finish;
	double delta_theta;
	double delta_l;

	double pos_start;
	double pos_finish;
	double delta_x;
	double time_start;
	double delta_time;
	double delta_time_full; //with scattering (Wigner) time delay

	double photon_En;

	enum ProcessType : short {Overflow = -2, None = -1, Elastic = 0, ResNBrS =1, Ionization = 2}; //Manager and Ar XS functions depend on this.
	short process; //1 is ionization and from 2 to max are excitations (process = ID + 1).
	std::vector<double> CrossSections; //for each process starting from Elastic=1
	std::vector<double> CrossSectionsSum; //helper for random process selection
	//Debug info:
	double deb_log_rand; //-ln R * coef. which is equal to integral of XS
	double deb_solver_y_left;
	double deb_solver_y_right;
	double deb_solver_E_left;
	double deb_solver_E_right;
	double deb_solver_E_delta;
	double deb_N_integral;
};

#endif
