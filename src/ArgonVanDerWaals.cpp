#include "ArgonVanDerWaalsParticle.h"

ArgonVanDerWaalsParticle::ArgonVanDerWaalsParticle(void) : Particle()
{
	name_ = ARGON_VAN_DER_WAALS_NAME;
	mass_ = 2*gSettings.PhysConsts()->Ar_mass_eV;
	width_ = 0;
	is_valid_ = true;

	XS_En_sweeper_ = ColoredInterval (gSettings.PhysConsts()->Dissoc_attachment_En_thresh, gSettings.PhysConsts()->XS_el_En_maximum, 1e-3);

	std::map<unsigned int, std::string> entry;
	entry[0u] = "\"Elastic scattering\""; //has zero XS for now
	entry[1u] = "\"Dissociative attachment\"";
	processes_["electron"] = entry;
}

ArgonVanDerWaalsParticle::~ArgonVanDerWaalsParticle() {};

unsigned int ArgonVanDerWaalsParticle::GetQauntStateSize(const Particle *target, double E, double theta, unsigned int process) const
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
	if (target->GetName()=="electron") {
		if (0 == process) {
			return 0;
		}
		if (1 == process) {
			return gSettings.PhysConsts()->Dissoc_attachment_XS;
		}
	}
	return 0;
}

double ArgonVanDerWaalsParticle::GetCrossSection(const Particle *target, double E, double theta, unsigned int process) const
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
	if (target->GetName()=="electron") {
		if (0 == process) {
			return 0;
		}
		return GetCrossSection(target, E, process) / (4*M_PI); //considered isotropic for other processes
	}
	return 0;
}

std::vector<const Particle*> ArgonVanDerWaalsParticle::GetFinalStates(const Particle *target, double E, double theta, unsigned int process) const
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

double ArgonVanDerWaalsParticle::GenerateScatterAngle(const Particle *target, double E, unsigned int process, double Rand) const
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
	if (target->GetName()=="electron") {
		Rand = Rand*2.0 - 1.0;
		return std::acos(Rand);
	}
	Rand = Rand*2.0 - 1.0;
	return std::acos(Rand);
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
	if (target->GetName()=="electron") {
		if (0 == process) {
			long double gamma_f = target->GetMass() / GetMass();
			double EnergyLoss = 2*(1-cos(theta))*E*gamma_f /pow(1 + gamma_f, 2);
			if (gSettings.PhysConsts()->resonance_En_loss) {
				double width_factor = std::pow(Width_3o2_/2.0, 2)/(std::pow(Width_3o2_/2.0, 2) + std::pow((E-En_3o2_)/2.0, 2));
				width_factor += std::pow(Width_1o2_/2.0, 2)/(std::pow(Width_1o2_/2.0, 2) + std::pow((E-En_1o2_)/2.0, 2));
				width_factor *= 0.5**gSettings.PhysConsts()->resonance_En_loss;
				EnergyLoss += width_factor;
			}
			return EnergyLoss;
		}

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
	if (target->GetName()=="electron") {
		return 0;
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
	if (target->GetName()=="electron") {
		if (0 == process) {
			return 0;
		}
		if (1 == process) {
			return -1*log(1-Rand)*gSettings.PhysConsts()->Argon_ion_decay_time;
		}
		return 0;
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
