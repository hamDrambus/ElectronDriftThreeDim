#include "ElectronParticle.h"

ElectronParticle::ElectronParticle(void) : Particle()
{
	name_ = "electron";
	mass_ = gSettings.PhysConsts()->e_mass_eV;
	width_ = 0;
	XS_En_sweeper_;
	processes_;
}

ElectronParticle::~ElectronParticle() {};

unsigned int ElectronParticle::GetQauntStateSize(const Particle *target, double E, double theta, unsigned int process) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GetQauntStateSize: Error: NULL target"<<std::endl;
		return 0;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GetQauntStateSize: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		return 0;
	}
	return 0;
}

double ElectronParticle::GetCrossSection(const Particle *target, double E, unsigned int process) const
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
	return 0;
}

double ElectronParticle::GetCrossSection(const Particle *target, double E, double theta, unsigned int process) const
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
	return 0;
}

std::vector<const Particle*> ElectronParticle::GetFinalStates(const Particle *target, double E, double theta, unsigned int process) const
{
	std::vector<const Particle*> out;
	if (NULL == target) {
		std::cerr << GetName() << "::GetFinalStates: Error: NULL target"<<std::endl;
		out.push_back(NULL);
		return out;
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GetFinalStates: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		out.push_back(NULL);
		return out;
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GetFinalStates: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		out.push_back(NULL);
		return out;
	}
	out.push_back(this);
	out.push_back(target);
	return out;
}

double ElectronParticle::GenerateScatterAngle(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	if (NULL == target) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: NULL target"<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	auto procs = processes_.find(target->GetName());
	if (processes_.end()==procs) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	if (procs->second.end()==procs->second.find(process)) {
		std::cerr << GetName() << "::GenerateScatterAngle: Error: unsupported process "<<process<<" for particle \""<<target->GetName()<<"\""<<std::endl;
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	Rand = Rand*2.0 - 1.0;
	return std::acos(Rand);
}

double ElectronParticle::GenerateEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const
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
	return 0;
}

//For neutral bremsstrahlung or deexcitation. Different from energy loss
double ElectronParticle::GeneratePhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const
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
	return 0;
}

double ElectronParticle::GenerateTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const
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
	return 0;
}

//Untabulated functions:
unsigned int ElectronParticle::GenerateUntabProcess(const Particle *target, double E, double theta, double Rand) const
{
	return GenerateProcess(target, E, theta, Rand);
}

double ElectronParticle::GenerateUntabScatterAngle(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GenerateScatterAngle(target, E, theta, process, Rand);
}

double ElectronParticle::GenerateUntabEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GenerateEnergyLoss(target, E, theta, process, Rand);
}

double ElectronParticle::GenerateUntabPhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GeneratePhoton(target, E, theta, process, Rand);
}

double ElectronParticle::GenerateUntabTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const
{
	return GenerateTimeDelay(target, E, theta, process, Rand);
}
