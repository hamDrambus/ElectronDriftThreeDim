#include "ParticleTable.h"

ParticleTable gParticleTable;

ParticleTable::ParticleTable() : electron_(NULL), argon_(NULL), argon_van_der_waals_(NULL)
{}

ParticleTable::~ParticleTable()
{
	if (electron_)
		delete electron_;
	if (argon_)
		delete argon_;
	if (argon_van_der_waals_)
		delete argon_van_der_waals_;
}

void ParticleTable::Load(void)
{
	electron_ = new ElectronParticle();
	argon_ = new ArgonParticle;
	argon_van_der_waals_ = new ArgonVanDerWaalsParticle(argon_);
	list_.push_back(electron_);
	list_.push_back(argon_);
	list_.push_back(argon_van_der_waals_);
	for (auto p = list_.begin(), p_end_ = list_.end(); p!=p_end_; ++p)
		if (!(*p)->isValid())
			std::cerr<<"Error ParticleTable::ParticleTable(): invalid particle "<<(*p)->GetName()<<std::endl;
}

const Particle* ParticleTable::GetParticle(std::string name) const
{
	for (auto i = list_.begin(), i_end_ = list_.end(); i!=i_end_; ++i) {
		if ((*i)->GetName() == name)
			return *i;
	}
	return NULL;
}

int ParticleTable::GetParticleID(std::string name) const
{
	for (std::size_t i = 0, i_end_ = list_.size(); i!=i_end_; ++i) {
		if (list_[i]->GetName() == name)
			return i;
	}
	return -1;
}

int ParticleTable::GetParticleID(const Particle* part) const
{
	for (std::size_t i = 0, i_end_ = list_.size(); i!=i_end_; ++i) {
		if (list_[i] == part)
			return i;
	}
	return -1;
}
