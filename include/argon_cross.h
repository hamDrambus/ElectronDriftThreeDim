#ifndef ARGON_CROSS_H
#define ARGON_CROSS_H

/*	TODO: add comprehensive explanation of used data and calculations of cross sections
 *
 */

#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <ctgmath>
#include "global_definitions.h"
#include "PolynomialFit.h"
#include "LegendrePolynomials.h"
#include "ColoredInterval.h"

class EnergyScanner
{
public:
	enum ScanType: short {ElasticXS, Resonance_3o2_XS, Resonance_1o2_XS, ResonancesXS,
			Resonance_3o2_DiffXS, Resonance_1o2_DiffXS, ResonancesDiffXS,
			DiffXS,	InelasticXS, ElasticResXS, XSIntegral, PlotElastic,
			PlotResonance_3o2, PlotResonance_1o2, PlotResonances,
			PlotDiffXS, PlotInelasticXS, PlotElasticResXS, PlotAllXS};
protected:
	unsigned int i;
	ScanType type_;
	ColoredRange energy_range_;
public:
	EnergyScanner(ScanType type);
	long double Next(int& err);
	void Reset(void);
};

class ArExperimental;
class InelasticProcess
{
protected:
	ArExperimental *ArExper_; //only for ionization energy for BB_XS(E);
	std::string name_;
	DataVector exp_XS_;
	unsigned int ID_;
	double En_threshold_;
	double Oscillator_strength_;
public:
	InelasticProcess(std::string name, unsigned int ID, double En, double F, std::vector<double> &Ens, std::vector<double> &XSs,
		ArExperimental *ArExper);
	double operator ()(double E); //returns cross section in 1e-16 cm^2
	double BB_XS(double E);
	double Exp_XS(double E);
	double get_En_thresh (void) const;
	std::string get_name (void) const;
	unsigned int get_ID (void) const;
};

class ArExperimental
{
protected:
	void read_inelastic(std::ifstream &inp, std::vector<InelasticProcess> &to);
public:
	std::vector<DataVector> phase_shifts_; //TODO: account for phase plus and phase minus
	DataVector total_elastic_cross;
	std::vector<InelasticProcess> excitations;
	std::vector<InelasticProcess> ionizations;
	short max_process_ID;
	double E_Ionization;
	ArExperimental(void);
	InelasticProcess * FindInelastic(short ID);
	unsigned int max_L (long double k);
	long double phase_shift (long double k, unsigned int l);
};

class ArAllData {
public:
	ArExperimental ArExper_;
	ArAllData(void);
	void argon_phase_values_exp(long double k, unsigned int l, long double &tan, long double &sin, long double &cos);
	void argon_phase_values_MERT5(long double k, unsigned int l, long double &tan, long double &sin, long double &cos);
	//E in eV, theta in radians, output is in m
	long double argon_cross_elastic_diff(long double E, long double theta);
	long double argon_cross_elastic(long double E);
	long double argon_cross_elastic_from_phases(long double E);
	long double argon_back_scatter_prob(long double E);
	long double argon_TM_forward(long double E);
	long double argon_TM_backward(long double E);

	long double argon_cross_resonance_3o2_diff(long double E, long double theta);
	long double argon_cross_resonance_1o2_diff(long double E, long double theta);
	long double argon_cross_resonance_3o2(long double E);
	long double argon_cross_resonance_1o2(long double E);
	long double argon_back_resonance_3o2_prob(long double E);
	long double argon_back_resonance_1o2_prob(long double E);
	long double argon_TM_forward_resonance_3o2(long double E);
	long double argon_TM_forward_resonance_1o2(long double E);
	long double argon_TM_backward_resonance_3o2(long double E);
	long double argon_TM_backward_resonance_1o2(long double E);
};

class ArDataTables //loads data from default files if presented. If not then values are calculated and files are created.
{
public:
	ArAllData ArAllData_;
protected:
	std::string total_cross_elastic_fname;
	std::string total_cross_resonance_3o2_fname;
	std::string total_cross_resonance_1o2_fname;

	DataVector total_cross_elastic_;
	DataVector total_cross_resonance_3o2_;
	DataVector total_cross_resonance_1o2_;
	void read_data (std::ifstream &inp, DataVector &data, long double y_factor = 1);
public:
	ArDataTables();
	long double CrossSection (double E, short type);
	long double TotalCrossSection (double E);
	long double XS_elastic(double E);
	long double XS_resonance_3o2(double E);
	long double XS_resonance_1o2(double E);

	double generate_Theta (double E, short type, double Rand); //TODO: tabulate

	void setOrder(int order);
	void setNused(int N);
	int getOrder(void);
	int getNused(void);
};

#endif

