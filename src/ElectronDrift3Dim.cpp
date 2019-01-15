//============================================================================
// Name        : ElectronDrift3Dim.cpp
// Author      : Frolov Egor geffdroid@gmail.com
// Version     :
// Copyright   : No copyrights
// Description :
//============================================================================

#include <iostream>
#include <TApplication.h>
#include "argon_cross.h"
#include "MTManager.h"
#include "tests.h"

void process_runs_in_thread(void* manager)
{
	((MTManager*)manager)->ProcessAll();
	TCondition* cond = ((MTManager*)manager)->getCondition();
	TMutex* mutex = ((MTManager*)manager)->getThreadMutex();
	if (0 != mutex->TryLock()) {//means that the main thread is waiting for the signal
		cond->Signal();
		//std::cout << "Signal()" << std::endl;
	}
	//std::cout << "Exiting thread" << std::endl;
}

void Process(int N_threads, unsigned int seed, unsigned int num_of_electrons, double concentration, double field, std::string root_fname) {
	if (N_threads < 1)
		N_threads = 1;
	std::vector<TThread*> pThreads;
	std::vector<MTManager*> _submanagers;
	std::vector<TMutex*> mutexes;
	std::vector<TMutex*> thread_mutexes;
	std::vector<TCondition*> conditions;
	std::vector<ArDataTables*> ar_data;

	ArDataTables ArDataTables_;
	if (N_threads > num_of_electrons)
		N_threads = num_of_electrons;
	int N_acc = 0;
	for (int n = 0; n < N_threads; ++n) {
		mutexes.push_back(new TMutex());
		conditions.push_back(new TCondition(mutexes[n]));
		thread_mutexes.push_back(new TMutex());
		ar_data.push_back(new ArDataTables(ArDataTables_)); //copy, no repeated readings and table constructions.
		if (n == (N_threads - 1)) {
			_submanagers.push_back(new MTManager(ar_data[n], n, num_of_electrons - N_acc, seed + n));
		} else {
			int N_e = num_of_electrons / N_threads;
			N_acc += N_e;
			_submanagers.push_back(new MTManager(ar_data[n], n, N_e, seed + n));
		}
		_submanagers[n]->setCondition(conditions[n]);
		_submanagers[n]->setThreadMutex(thread_mutexes[n]);
		pThreads.push_back(new TThread(("MTManager_" + std::to_string(n)).c_str(),
			&process_runs_in_thread, _submanagers[n]));
	}
	MTManager test_man(&ArDataTables_, -1, 1, seed);
	//test_all(ar_data[1]);
	//test_man.Test();
	
	for (int n = 0; n < N_threads; ++n) {
		_submanagers[n]->SetParameters(concentration, field);
		pThreads[n]->Run(); //if it is the last iteration, submanager (AnalysisManager) clears its data
	}
	//TThread::Ps();
	for (int n = 0; n < N_threads; ++n) {
		if (0 != thread_mutexes[n]->TryLock()) { //thread is already executed, so no wait required
		}
		else {
			conditions[n]->Wait();
		}
		thread_mutexes[n]->UnLock();
		if (0 != n)
			_submanagers[0]->Merge(_submanagers[n]);
	}
	//Merged() into _submanagers[0];
	_submanagers[0]->WriteHistory(root_fname);
	
	for (int n = 0; n < N_threads; ++n) {
		pThreads[n]->Delete();
		delete _submanagers[n];
		delete ar_data[n];
		conditions[n]->Delete();
		delete mutexes[n];
		delete thread_mutexes[n];
	}
}

int main(int argn, char * argv[]) {
	//int n_par = 0;
	//char **f = NULL;
	//TApplication app("test_app",&n_par,f);
	std::string a;
	std::string root_fname = "Output/eData_3Td.root";
	double Td = 7; //=E/N in 1e-21 in Si
	unsigned int seed = 42;
	unsigned int num_of_electrons = 1;
	double pressure = 1.015e5;
	double temperature = 87;
	if (argn!=1) {
		if (argn<5) {
			std::cout<<"Input error, the following parameters must be specified:"<<std::endl;
			std::cout<<"(int)<seed>\t(double)<Td>\t(int)<number of electrons>\t(string)<output file name>"<<std::endl;
			return 0;
		}
		if (argn>5)
			std::cout<<"Warning! Only first 4 parameters are used."<<std::endl;
		seed = std::stoi(argv[1]);
		Td = std::stod(argv[2]);
		num_of_electrons = std::abs(std::stoi(argv[3]));
		root_fname = argv[4];
	}
	ensure_file(root_fname);
	ensure_file("tests/t");

	double concentration = pressure / (temperature*boltzmann_SIconst);
	double field = Td*1e-21 * concentration;
	auto start = std::chrono::system_clock::now();
	Process(THREADS_NUMBER_, seed, num_of_electrons, concentration, field, root_fname);
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> diff = end-start;
	std::cout<<"Elapsed time:\t"<<diff.count()<<" s"<<std::endl;
	std::cout<<"Enter something: ";
	//std::cin>>a;
	//app.Run();
	return 0;
}
