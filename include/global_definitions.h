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
#include <Rtypes.h>

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

#endif //GLOBAL_DEFINITIONS_H
