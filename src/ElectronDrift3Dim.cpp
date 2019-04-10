//============================================================================
// Name        : ElectronDrift3Dim.cpp
// Author      : Frolov Egor geffdroid@gmail.com
// Version     :
// Copyright   : No copyrights
// Description : uses boost and ROOT libraries
//============================================================================

#include <iostream>
#include <TApplication.h>
#include "argon_cross.h"
#include "MTManager.h"
#include "tests.h"

void process_runs_in_thread(void* manager)
{
	((MTManager*)manager)->ProcessAll();
}

bool Process(void) {
	std::cout << "Starting simulation of electron drift:" << std::endl;
	std::cout << "  Material: pure gaseous argon" << std::endl;
	auto start = std::chrono::system_clock::now();
	unsigned int N_threads = gSettings.ProgConsts()->thread_number;
	if (N_threads < 1u)
		N_threads = 1u;
	std::vector<std::thread> pThreads;
	std::vector<MTManager*> _submanagers;
	std::vector<ArDataTables*> ar_data;
	//TODO: It is currently implied that calling FunctionTable methods is thread-safe because there is only reading and no changes in their internal states.
	//It would be better to fix this fact in the code explicitly (cost methods, locks, etc.)
	ArDataTables ArDataTables_;
	for (unsigned int n = 0u; n < N_threads; ++n) {
		ar_data.push_back(new ArDataTables(ArDataTables_)); //shallow copy, no repeated table constructions.
		_submanagers.push_back(new MTManager(ar_data[n], n));
		pThreads.push_back(std::thread());
	}
	if (true == gSettings.ProgConsts()->is_test_version) {
		//MTManager test_man(&ArDataTables_, -1, 1, seed);
		//ArDataTables_.ArAllData_.argon_cross_elastic_diff(0.1, 0.5);
		test_all(&ArDataTables_);
		//test_man.Test();
	}
	for (std::size_t run = 0, run_end_ = gSettings.ProgConsts()->run_specifics.size(); run != run_end_; ++run) {
		auto run_start_t = std::chrono::system_clock::now();
		const RunParameters *rp = &gSettings.ProgConsts()->run_specifics[run];
		double concentration = gSettings.ProgConsts()->pressure / (gSettings.ProgConsts()->temperature*gSettings.PhysConsts()->boltzmann_SI);
		unsigned int N_thread_active = (N_threads > rp->n_electrons ? rp->n_electrons: N_threads);
		std::cout << "Starting run #" << run << std::endl;
		std::cout << "  Electric field = \t" << rp->field << " Td" << std::endl;
		std::cout << "  Temperature = \t" << gSettings.ProgConsts()->temperature<< " K" << std::endl;
		std::cout << "  Pressure = \t" << gSettings.ProgConsts()->pressure << " Pa" << std::endl;
		std::cout << "  Concentration = \t" << concentration << " m^-3" << std::endl;
		std::cout << "  Number of electrons = \t" << rp->n_electrons << std::endl;
		std::cout << "  Drift length = \t" << rp->drift_distance << " m" << std::endl;
		std::cout << "  Thread number = \t" << N_thread_active << std::endl;
		std::cout << "  Output file = \t \"" << rp->output_file << "\"" << std::endl;
		std::cout << "  Random generator seed = \t" << rp->seed << std::endl;

		ensure_file(rp->output_file);
		int N_extra = rp->n_electrons % (N_thread_active>0 ? N_thread_active : 1);
		std::vector<int> N_e(N_thread_active, 0);
		for (unsigned int n = 0u; n < N_thread_active; ++n) { //distribute electrons among the processes as evenly as possible
			N_e[n] = rp->n_electrons / N_thread_active + ((N_extra>0) ? 1 : 0);
			--N_extra;
		}

		for (unsigned int n = 0u; n < N_thread_active; ++n) {
			_submanagers[n]->Clear();
			_submanagers[n]->setParameters(concentration, 1e-21*rp->field*concentration, rp->drift_distance);
			_submanagers[n]->setNelectons(N_e[n]);
			_submanagers[n]->setInitialSeed(rp->seed + n);
			_submanagers[n]->setRunIndex(run);
			pThreads[n]=std::thread(process_runs_in_thread, _submanagers[n]);
		}
		for (unsigned int n = 0u; n < N_thread_active; ++n) {
			std::this_thread::sleep_for(std::chrono::milliseconds(1000));
			pThreads[n].join();
			if (0 != n)
				_submanagers[0]->Merge(_submanagers[n]);
		}
		//Merged() into _submanagers[0];
		if (N_thread_active>0)
			_submanagers[0]->WriteHistory(rp->output_file);

		auto run_end_t = std::chrono::system_clock::now();
		std::chrono::duration<double> diff = run_end_t - run_start_t;
		std::cout << "Finished run #" << run << std::endl;
		std::cout << "  Electric field = \t" << rp->field << " Td" << std::endl;
		std::cout << "  Temperature = \t" << gSettings.ProgConsts()->temperature << " K" << std::endl;
		std::cout << "  Pressure = \t" << gSettings.ProgConsts()->pressure << " Pa" << std::endl;
		std::cout << "  Concentration = \t" << concentration << " m^-3" << std::endl;
		std::cout << "  Number of electrons = \t" << rp->n_electrons << std::endl;
		std::cout << "  Drift length = \t" << rp->drift_distance << " m" << std::endl;
		std::cout << "  Thread number = \t" << N_thread_active << std::endl;
		std::cout << "  Output file = \t \"" << rp->output_file << "\"" << std::endl;
		std::cout << "  Random generator seed = \t" << rp->seed << std::endl;
		std::cout << "  Elapsed time  = \t" << diff.count() << " s." << std::endl;
	}

	for (int n = 0; n < N_threads; ++n) {
		delete _submanagers[n];
		delete ar_data[n];
	}
	ArDataTables_.DeleteData(); //TODO: implement shared pointers on data tables properly.
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end - start;
	std::cout << "Finished simulation." << std::endl;
	std::cout << "Elapsed time  = \t" << diff.count() << " s." << std::endl;
	return true;
}

//input is settings.xml file. All general and specific info is stored there.
int main(int argn, char * argv[]) {
	std::string settings_fname;
	if (argn != 1) {
		if (argn>2)
			std::cout << "Warning! Only single parameter (settings filename) is used." << std::endl;
		settings_fname = argv[1];
	} else {
		std::cout << "Using default settings file \"settings.xml\"" << std::endl;
		settings_fname = "settings.xml";
	}
	gSettings.Load(settings_fname);
	if (!gSettings.isValid()) {
		return 1;
	}
	ensure_file("tests/t");
	TApplication* app = NULL;
	if (gSettings.ProgConsts()->is_test_version) {
		int n_par = 0;
		char **f = NULL;
		app = new TApplication("test_app",&n_par,f);
	}
	Process();
	if (gSettings.ProgConsts()->is_test_version) {
		if (NULL!=app)
			app->Run();
		std::string a;
		std::cout<<"Enter something: ";
		std::cin>>a;
	}
	return 0;
}
