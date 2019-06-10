#ifndef MANAGER_H_
#define MANAGER_H_

#include <fstream>

#include <boost/random.hpp>
#ifndef _NO_CERN_ROOT
#include <TRandom3.h>
#include <TRandom2.h>
#include <TRandom1.h>
#include <TFile.h>
#include <TTree.h>
#endif //_NO_CERN_ROOT

#include "global_definitions.h"
#include "Settings.h"
#include "ParticleTable.h"
#include "Mixture.h"

class Manager
{
protected:
	//TODO: create random generator wrapper in order to place all #ifndefs into single place
#ifndef _NO_CERN_ROOT
	TRandom * random_generator_;
#endif //_NO_CERN_ROOT
	boost::random::hellekalek1995 *boost_random_generator;
	boost::optional<unsigned long int> initial_seed_;
	boost::optional<double> Concentration_;
	boost::optional<double> eField_;
	boost::optional<double> Coefficient_;
	boost::optional<double> Drift_distance_;
	boost::optional<std::size_t> run_index_;

	unsigned long int e_first_seed_;
#ifndef _NO_CERN_ROOT
	TTree * sim_data_;
	TTree * processes_data_;
#endif //_NO_CERN_ROOT
	std::vector<std::vector<unsigned long int> > processes_counters_; //particle_ID->process_ID
	std::vector<std::vector<short> > processes_IDs_; //particle_ID->process_ID
	std::vector<std::vector<std::string> > processes_legends_; //particle_ID->process_ID->NULL terminated String

	Event event_;
	const Mixture *material_;
	int skip_counter_;
	bool skipping_early_events;
	unsigned long int num_of_events;
	//Event current_event;

	void DoStepLength (Event &event);
	void Solve (long double LnR, Event &event);
	void Solve_table (long double LnR, Event &event);
	void Solve_test (long double LnR, Event &event);
	void DoScattering(Event &event);
	void PostStepAction(Event &event);
	void DoGotoNext(Event &event);
	void InitTree (void);
	//long double XS_integral(long double from, long double to, long double Eny, Event &event);
	long double XS_integral_for_test(long double from, long double to, long double Eny, long double dE);
	//long double XS_integral_table(long double from, long double to, long double Eny, Event &event);
	long double Path_integral (long double x, long double cos_th);
	double Uniform();
public:
	virtual bool isReady(void) const;
	virtual void Initialize(void);
	virtual void Initialize(Event &event);
	virtual void Clear(void); //fully clears manager state, subsequent resetting of parameters is required. 
	void DoStep(Event &event);
	bool IsFinished(Event &event);
	void setParameters(double Concetr /*in SI*/, double E /*in SI*/, double drift_distance /*in m*/);
	void setParameters(double T /*in K*/, double Pressure /*in SI*/, double E /*in SI*/, double drift_distance /*in m*/);
	bool setInitialSeed(unsigned long int seed);
	boost::optional<unsigned long int> getInitialSeed(void) const;
	bool setRunIndex(std::size_t index);
	boost::optional<std::size_t> getRunIndex(void) const;
	Manager(const Mixture *material);
	virtual ~Manager();
	void LoopSimulation(void);
	void WriteHistory(std::string root_fname);
	void Test(void);
};

#endif 
