#include "MTManager.h"

MTManager::MTManager(ArDataTables *Ar_tables, int instance, int N_electrons, UInt_t RandomSeed) :
	Manager(Ar_tables, RandomSeed), condition_(NULL), thread_mutex_(NULL), instance_(instance), N_electrons_(N_electrons)
{}

void MTManager::ProcessAll(void)
{
	unsigned int incr = N_electrons_>1000 ? 1 + N_electrons_ / 100 : 1;
	for (unsigned int i = 0; i<N_electrons_; ++i) {
		this->LoopSimulation();
		if (0 == (i + 1) % incr)
			std::cout <<"Thread "<<instance_<<": "<< i + 1 << "/" << N_electrons_ <<" Seed "<<start_seed_<<std::endl;
	}
	if (0 != N_electrons_%incr) {
		std::cout << "Thread " << instance_ << ": " << N_electrons_ << "/" << N_electrons_<<" Seed "<<start_seed_ << std::endl;
	}
}

void MTManager::setCondition(TCondition* cond)
{	condition_ = cond;}
TCondition* MTManager::getCondition(void) const
{	return condition_;}
void MTManager::setThreadMutex(TMutex* mutex)
{	thread_mutex_ = mutex;}
TMutex* MTManager::getThreadMutex(void) const
{	return thread_mutex_;}

void MTManager::Merge(MTManager *with)
{
	TList *list = new TList;
	list->Add(with->sim_data_);
	sim_data_->Merge(list);
	with->sim_data_->Reset();
	for (int i = 0, i_end_ = processes_size_; i!=i_end_; ++i) {
		processes_counters_[i]+=with->processes_counters_[i];
		with->processes_counters_[i] = 0;
	}
}
