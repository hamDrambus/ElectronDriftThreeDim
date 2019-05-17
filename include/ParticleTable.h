#ifndef PARTICLE_TABLE_H_
#define PARTICLE_TABLE_H_

#include "ArgonParticle.h"
#include "ArgonVanDerWaalsParticle.h"
#include "ElectronParticle.h"

class ParticleTable {
protected:
	std::vector<const Particle*> list_;
	ElectronParticle electron_;
	ArgonParticle argon_;
	ArgonVanDerWaalsParticle argon_van_der_waals_;
public:
	ParticleTable();
	~ParticleTable();
	const Particle* GetParticle(std::string name) const;
	const Particle* GetParticle(std::size_t index) const {
		if (index>=list_.size())
			return NULL;
		return list_[index];
	}
	std::size_t GetNParticle(void) const {
		return list_.size();
	}
	int GetParticleID(std::string name) const;
	int GetParticleID(const Particle*) const;
};

extern ParticleTable gParticleTable;

#endif //PARTICLE_TABLE_H_
