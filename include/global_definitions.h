#ifndef GLOBAL_DEFINITIONS_H
#define GLOBAL_DEFINITIONS_H

#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <map>
#include <boost/optional.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lexical_cast.hpp>

#if defined(_WIN32)||defined(_WIN64)
#define NOMINMAX
#ifndef _NO_CERN_ROOT 
#include "Windows4Root.h"
#else //_NO_CERN_ROOT
#include <Windows.h>
#endif //_NO_CERN_ROOT
#include <direct.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif
#include <math.h>
#include <ctgmath>
#ifndef _NO_CERN_ROOT
#include <Rtypes.h>
#endif //_NO_CERN_ROOT

#if defined(_WIN32)||defined(_WIN64)
#define INVOKE_GNUPLOT(a) system(("start \"\" \"%GNUPLOT%\\gnuplot.exe\" --persist \"" + a + "\"").c_str())
#else //defined(_WIN32)||defined(_WIN64) 
#define INVOKE_GNUPLOT(a) system(("gnome-terminal -- bash -c \"gnuplot \"" + a +"\"\"").c_str());
#endif //defined(_WIN32)||defined(_WIN64)

std::string strtoken(std::string &in, std::string break_symbs);
void open_output_file(std::string name, std::ofstream &str, std::ios_base::openmode _mode);
void ensure_file(std::string fname); //makes sure file can be created later on
void ensure_folder(std::string folder);
char* c_str_cp (const std::string &str);

#define ELECTRON_NAME "electron"
#define ARGON_NAME "Argon"
#define ARGON_VAN_DER_WAALS_NAME "Argon Van der Waals Molecule"

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
	short process; //process meanings depend on interacting particle ID. ProcessType is now deprecated for positive values.
	std::size_t particle_ID; //At the moment electron interacts only with Argon and Argon Van der Waals molecule
	//particle_ID corresponds to the indexes of gParticleTable
	std::size_t particle_ID_finish;
	std::size_t particle_ID_target;
	//Debug info:
	double deb_log_rand; //-ln R * coef. which is equal to integral of XS
	double deb_solver_y_left;
	double deb_solver_y_right;
	double deb_solver_E_left;
	double deb_solver_E_right;
	double deb_solver_E_delta;
	double deb_N_integral;
};

#endif //GLOBAL_DEFINITIONS_H
