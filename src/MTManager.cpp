#include "MTManager.h"

MTManager::MTManager(const Mixture *material, int instance) :
	Manager(material), instance_(instance)
{}

void MTManager::ProcessAll(void)
{
	if (!isReady()) {
		std::cout << "MTManager::ProcessAll: some of the parameters are not initiaized, exiting" << std::endl;
		return;
	}
	unsigned int incr = *N_electrons_> 1000 ? 1 + *N_electrons_ / 100 : 1;
	for (unsigned int i = 0; i<N_electrons_; ++i) {
		this->LoopSimulation();
		if (0 == (i + 1) % incr)
			std::cout <<"Thread "<<instance_<<": "<< i + 1 << "/" << *N_electrons_ <<" Seed "<<e_first_seed_<<std::endl;
	}
	if (0 != *N_electrons_%incr || 0u == N_electrons_) {
		std::cout << "Thread " << instance_ << ": " << *N_electrons_ << "/" << *N_electrons_<<" Seed "<< e_first_seed_ << std::endl;
	}
}

bool MTManager::setNelectons(unsigned int Ne)
{
	N_electrons_ = Ne;
	return true;
}

boost::optional<unsigned int> MTManager::getNelectons(void) const
{
	return N_electrons_;
}

bool MTManager::isReady(void) const
{
	return Manager::isReady() && boost::none != N_electrons_;
}

void MTManager::Clear(void)
{
	N_electrons_ = boost::none;
	run_index_ = boost::none;
	Manager::Clear();
}

void MTManager::Merge(MTManager *with)
{
	TList *list = new TList;
	list->Add(with->sim_data_);
	sim_data_->Merge(list);
	with->sim_data_->Reset();
	for (std::size_t i = 0, i_end_ = processes_counters_.size(); i != i_end_; ++i) {
		for (std::size_t pr = 0, pr_end_ = processes_counters_[i].size(); pr != pr_end_; ++pr) {
			processes_counters_[i][pr]+=with->processes_counters_[i][pr];
			with->processes_counters_[i][pr] = 0;
		}
	}
}
