#include "ArgonVanDerWaalsParticle.h"

ArgonVanDerWaalsParticle::ArgonVanDerWaalsParticle(ArgonParticle* argon): ArgonParticle(*argon)
{
	name_ = ARGON_VAN_DER_WAALS_NAME;
	mass_ = 2*gSettings.PhysConsts()->Ar_mass_eV;
	width_ = 0;
	is_valid_ = true;

	XS_En_sweeper_ += ColoredInterval (gSettings.PhysConsts()->Dissoc_attachment_En_thresh, gSettings.PhysConsts()->XS_el_En_maximum, 1e-3);

	auto entry = processes_.find(ELECTRON_NAME);
	if (processes_.end() != entry) { //Add new process to the argon
		std::size_t max_index = entry->second.size();
		(entry->second)[max_index] = "\"Dissociative attachment\"";
	}
}

ArgonVanDerWaalsParticle::~ArgonVanDerWaalsParticle() {};

double ArgonVanDerWaalsParticle::GetCrossSection(const Particle *target, double E) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GetCrossSection: Error: NULL target" << std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end() == procs) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported target particle \"" << target->GetName() << "\"" << std::endl;
		return 0;
	}
	if (target->GetName() == ELECTRON_NAME) {
		double output = 2 * ArgonParticle::GetCrossSection(target, E);
		if (E > gSettings.PhysConsts()->Dissoc_attachment_En_thresh)
			output += gSettings.PhysConsts()->Dissoc_attachment_XS;
		return output;
	}
	return 0;
}

double ArgonVanDerWaalsParticle::GetCrossSection(const Particle *target, double E, unsigned int process) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GetCrossSection: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GetCrossSection: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()==ELECTRON_NAME) {
		if ((procs->second.size() - 1) == process) { //Dissociative attachment is set as last process
			if (E > gSettings.PhysConsts()->Dissoc_attachment_En_thresh)
				return gSettings.PhysConsts()->Dissoc_attachment_XS;
			return 0;
		}
		return 2.0 * ArgonParticle::GetCrossSection(target, E, process);
	}
	return 0;
}

double ArgonVanDerWaalsParticle::GetCrossSection(const Particle *target, double E, double theta, unsigned int process) const
{
	return ArgonParticle::GetCrossSection(target, E, theta, process);
}

std::vector<const Particle*> ArgonVanDerWaalsParticle::GetFinalStates(const Particle *target, double E, double theta, unsigned int process) const
{
	return ArgonParticle::GetFinalStates(target, E, theta, process);
}

double ArgonVanDerWaalsParticle::GenerateScatterAngle(const Particle *target, double E, unsigned int process, double Rand) const
{
	return ArgonParticle::GenerateScatterAngle(target, E, process, Rand);
}

double ArgonVanDerWaalsParticle::GenerateEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateEnergyLoss: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateEnergyLoss: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateEnergyLoss: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()==ELECTRON_NAME) {
		if ((procs->second.size() - 1) == process) { //Dissociative attachment is set as last process, same energy loss as in elastic scaterring
			return ArgonParticle::GenerateEnergyLoss(target, E, theta, 0, Rand);
		}
		return ArgonParticle::GenerateEnergyLoss(target, E, theta, process, Rand);
	}
	return 0;
}

//For neutral bremsstrahlung or deexcitation. Different from energy loss
double ArgonVanDerWaalsParticle::GeneratePhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GeneratePhoton: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GeneratePhoton: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GeneratePhoton: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName()== ELECTRON_NAME) {
		if ((procs->second.size() - 1) == process) { //Dissociative attachment is set as last process, same photon as in elastic scaterring
			return ArgonParticle::GeneratePhoton(target, E, theta, 0, Rand);
		}
		return ArgonParticle::GeneratePhoton(target, E, theta, process, Rand);
	}
	return 0;
}

double ArgonVanDerWaalsParticle::GenerateTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateTimeDelay: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	if (target->GetName() == ELECTRON_NAME) {
		if ((procs->second.size() - 1) == process) { //Dissociative attachment is set as last process
			return -1 * log(1 - Rand)*gSettings.PhysConsts()->Argon_ion_decay_time;
		}
		return ArgonParticle::GenerateTimeDelay(target, E, theta, process, Rand);
	}
	return 0;
}

//Untabulated functions:
unsigned int ArgonVanDerWaalsParticle::GenerateUntabProcess(const Particle *target, double E, double Rand) const
{
	return GenerateProcess(target, E, Rand);
}

double ArgonVanDerWaalsParticle::GenerateUntabScatterAngle(const Particle *target, double E, unsigned int process, double Rand) const
{
	return GenerateScatterAngle(target, E, process, Rand);
}

double ArgonVanDerWaalsParticle::GenerateUntabEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GenerateEnergyLoss(target, E, theta, process, Rand);
}

double ArgonVanDerWaalsParticle::GenerateUntabPhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GeneratePhoton(target, E, theta, process, Rand);
}

double ArgonVanDerWaalsParticle::GenerateUntabTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GenerateTimeDelay(target, E, theta, process, Rand);
}
