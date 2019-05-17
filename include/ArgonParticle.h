#ifndef ARGON_PARTICLE_H_
#define ARGON_PARTICLE_H_

#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include "global_definitions.h"
#include "PolynomialFit.h"
#include "FunctionTable.h"
#include "Settings.h"
#include "LegendrePolynomials.h"
#include "ColoredInterval.h"
#include "Particle.h"

class EnergyScanner
{
public:
	enum ScanType: short {ElasticXS,
			InelasticXS, ElasticResXS, XSIntegral, PlotElastic,
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
	double operator ()(double E) const; //returns cross section in 1e-16 cm^2
	double BB_XS(double E) const;
	double Exp_XS(double E) const;
	double get_En_thresh (void) const;
	std::string get_name (void) const;
	unsigned int get_ID (void) const;
};

class ArExperimental
{
protected:
	void read_inelastic(std::ifstream &inp, std::vector<InelasticProcess> &to);
	bool is_valid_;
public:
	std::vector<DataVector> phase_shifts_pos_;
	std::vector<DataVector> phase_shifts_neg_;
	DataVector total_elastic_cross;
	std::vector<InelasticProcess> excitations;
	std::vector<InelasticProcess> ionizations;
	short max_process_ID;
	double E_Ionization;
	double First_excitation_En;
	ArExperimental(void);
	const InelasticProcess * FindInelastic(short ID) const;
	unsigned int max_L (long double k) const;
	void phase_shift (long double k, unsigned int l, long double &ps_pos, long double &ps_neg) const;
	bool isValid(void) const;
};

/*This class reads some experimental values from ./data
* and calculates Ar cross sections or factors for 1D simulation.
* MERT5 fit parameters are hard-coded here as well, not in the ./data
*/
class ArAllData {
protected:
	long double argon_cross_resonance_3o2(long double E) const;
	long double argon_cross_resonance_1o2(long double E) const;
	//^Separate functions for resonances are required because elastic XS is taken from experiment by extrapolating.
public:
	ArExperimental ArExper_;
	ArAllData(void);
	//returns positive (J=L+1/2) and negative (J=L-1/2) phaseshifts
	void argon_phase_values_exp(long double k, unsigned int l, long double &ps_p, long double &ps_n) const;
	void argon_phase_values_MERT5(long double k, unsigned int l, long double &ps_p, long double &ps_n) const;
	//E in eV, theta in radians, output is in 1e-20 m^2
	long double argon_cross_elastic_diff(long double E, long double theta, int mode = 0) const;
	long double argon_cross_elastic(long double E, int mode = 0) const;

	long double argon_back_scatter_prob(long double E);
	long double argon_TM_forward(long double E);
	long double argon_TM_backward(long double E);

	//have to make these 2 public for testing. UPD: these functions are useless. Time delay is calculated differently
	long double argon_scatter_probability_j(long double E, long double theta, int J, int L, int mode=0) const;
	long double argon_time_delay_j(long double E, int J, int L) const;
	long double argon_scatter_spin_flip_amplitude_sq(long double E, long double theta, int mode=0) const;
	long double argon_scatter_spin_nonflip_amplitude_sq(long double E, long double theta, int mode=0) const;

	//Wrong approach:
	//long double argon_delay_1o2_probability(long double E, long double theta);
	//long double argon_delay_3o2_probability(long double E, long double theta);
	long double argon_delay_spin_flip (long double E, long double theta, int mode =0) const;
	long double argon_delay_spin_nonflip (long double E, long double theta, int mode =0) const;
	//long double argon_delay_spin_flip_prob (long double E, long double theta, int mode =0);
	long double argon_delay_spin_nonflip_prob (long double E, long double theta, int mode =0) const;

	long double argon_ResNBrS_spectrum(long double W, long double E) const; //Normalization constant is arbitrary. Not dependant on E
	long double argon_ResNBrS_XS(long double E) const; //Normalization constant is taken from "global_definitions.h"
};

class ArgonParticle : public Particle {
public:
	ArAllData ArAllData_;
	//No deep copy
	FunctionTable *theta_table_; 	//shared between threads after it is initialized
	FunctionTable *time_delay_spin_flip_table_; 	//shared between threads after it is initialized
	FunctionTable *time_delay_spin_nonflip_table_; 	//shared between threads after it is initialized
	FunctionTable *time_delay_spin_nonflip_prob_table_; 	//shared between threads after it is initialized
protected:
	std::string total_cross_elastic_fname; //contains Feshbach resonances
	std::string theta_table_fname;
	std::string time_delay_spin_nonflip_prob_fname;
	std::string time_delay_spin_flip_fname;
	std::string time_delay_spin_nonflip_fname;
	std::string total_resonance_NBrS_spectrum_fname;

	DataVector total_cross_elastic_;
	DataVector total_resonance_NBrS_spectrum_; //Stores probability function, not the spectrum itself
	void read_data (std::ifstream &inp, DataVector &data, long double y_factor = 1);
	bool generate_total_cross_elastic_table(void);
	bool generate_theta_table(void);
	bool generate_time_delay_spin_flip_table(void);
	bool generate_time_delay_spin_nonflip_table(void);
	bool generate_time_delay_spin_nonflip_prob_table(void);
	bool generate_ResNBrS_spectrum_table(void);
public:

	void setOrder(int order)
	{
		total_cross_elastic_.setOrder(order);
	}
	void setNused(int N)
	{
		total_cross_elastic_.setNused(N);
	}
	int getOrder(void)
	{
		return total_cross_elastic_.getOrder();
	}
	int getNused(void)
	{
		return total_cross_elastic_.getNused();
	}

	void DeleteData(void);
public:
	ArgonParticle(void);
	virtual ~ArgonParticle();
	ArgonParticle (const ArgonParticle & ) = default;

	virtual bool isValid(void) const;

	virtual unsigned int GetQauntStateSize(const Particle *target, double E, double theta, unsigned int process) const;
	virtual double GetCrossSection(const Particle *target, double E, unsigned int process) const;
	virtual double GetCrossSection(const Particle *target, double E, double theta, unsigned int process) const;
	virtual std::vector<const Particle*> GetFinalStates(const Particle *target, double E, double theta, unsigned int process) const;

	virtual double GenerateScatterAngle(const Particle *target, double E, double theta, unsigned int process, double Rand) const;
	virtual double GenerateEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const;
	//For neutral bremsstrahlung or deexcitation
	virtual double GeneratePhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const;
	virtual double GenerateTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const;

	//Untabulated functions:
	virtual unsigned int GenerateUntabProcess(const Particle *target, double E, double theta, double Rand) const;
	virtual double GenerateUntabScatterAngle(const Particle *target, double E, double theta, unsigned int process, double Rand) const;
	virtual double GenerateUntabEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const;
	virtual double GenerateUntabPhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const;
	virtual double GenerateUntabTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const;
};

#endif //ARGON_PARTICLE_H_
