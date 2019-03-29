#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#if defined(__WIN32__)
#define NOMINMAX
#include "Windows4Root.h"
#include <direct.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif
#include <math.h>
#include <ctgmath>

#define THREADS_NUMBER_ 6
#define L_MAX_ 10
#define XS_EL_EN_MAXIMUM_ 20.0
#define XS_EL_EN_MINIMUM_ 1e-3
#define PHASES_EN_MAXIMUM_ 12
#define PHASES_EN_MINIMUM_ 5e-3
#define THRESH_E_PHASES_ 0.24
//^MERT 5 is applied only below this energy for backward probability and TM factors
#define THRESH_E_XS_ 1.0
//^MERT 5 is applied only below this energy for total cross section
/* These 2 thresholds are different because MERT5 excellently fits experimental total XS up until 1 eV (considering uncertainties).
 * MERT5 fit parameters are taken for the used experimental total XS (Kurokawa et al.)
 * On the other hand phase shift values obtained from experimental differential XS and MERT5 have discrepancy
 * from around 0.24 eV and above. Because of this P backward and TM have discontinuity.
 * Hence, even though MERT5 have right phase shifts around 1eV, situation is following:
 * From 5e-3 to 0.24 eV MERT5 PS are used for Pb and TM and diff. XS.
 * From 5e-3 to 1 eV MERT5 PS are used for total elastic XS.
 * Experimental PS are used otherwise.
 */
#define D_EN_SMOOTH_ 0.1
//^ MERT5 and experimental totalXS
#define EN_MAXIMUM_ 16.0
//^
#define ANGLE_POINTS_ 1001
//#define ANGLE_UNIFORM_
//^temporary for v9.x
#define En_3o2_ 11.103
#define En_1o2_ 11.270
#define Width_3o2_ 2.3e-3
#define Width_1o2_ 2.2e-3
//#define Width_3o2_ 2.0e-2
//#define Width_1o2_ 2.0e-2
//^in eV
#define RESONANCE_EN_LOSS_FACTOR_ 1.0
#define RESONANCE_EN_LOSS_ 2.5e-2
#define DRIFT_DISTANCE_ 3e-3
//^in m
#define DRIFT_DISTANCE_HISTORY 0
//^in m. Write info about electron only starting from this position
#define SKIP_HISTORY_ 0

const long double e_charge_SIconst = 1.60217662e-19; //in coulombs (SI)
const long double e_mass_SIconst = 9.10938356e-31; //in kg (SI)
const long double e_mass_eVconst = 5.109989461e5; //in eV
const long double h_bar_SIconst = 1.054571800e-34; //in SI
const long double h_bar_eVconst = 6.582119514e-16; //in eV*s
const long double a_bohr_SIconst = 5.2917721092e-1; //in meters (SI) multiplied by e10 for XS to be in 1e-20 m2
const long double Ry_eVconst = 13.605693; //Rydberg constant
const long double a_h_bar_2e_m_e_SIconst = 2.711063352e-1; //in SI a_bohr*sqrt(2*Me*e)/h_bar
const long double boltzmann_SIconst = 1.38064852e-23; //SI
//const double resonance_time_const = 2.9918725e-13;//in s. = h_bar /Width
const double Ar_mass_eVconst = 3.726e10; //eV

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
	double delta_time_full; //with resonance delay
	enum ProcessType : short {Overflow = -2, None = -1, Elastic = 0, Ionization = 1}; //Manager and Ar XS functions depend on this.
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
