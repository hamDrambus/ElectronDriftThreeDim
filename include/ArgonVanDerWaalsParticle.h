#ifndef ARGON_VAN_DER_WAALS_PARTICLE_H_
#define ARGON_VAN_DER_WAALS_PARTICLE_H_

#include "global_definitions.h"
#include "Particle.h"

class ArgonVanDerWaalsParticle : public Particle {
public:
	ArgonVanDerWaalsParticle(void);
	virtual ~ArgonVanDerWaalsParticle();
	ArgonVanDerWaalsParticle (const ArgonVanDerWaalsParticle & ) = default;

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

#endif //ARGON_VAN_DER_WAALS_PARTICLE_H_
