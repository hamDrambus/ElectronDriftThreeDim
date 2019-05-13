#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "global_definitions.h"
#include "Settings.h"

//TODO: this class and its daughters are required for generalization on scattering of e in mixtures and potentially other particles as well
class Particle {
protected:
	//Particle type --> list of processes
	std::map<std::string, std::map<unsigned int, std::string> > processes_;
	std::string name_; //must be unique for every derived class
	double width_; //in eV
	double mass_; //in eV
	bool is_valid_;
public:
	Particle(void) {
		is_valid_ = false;
		name_ = "GeneralParticle";
	}
	virtual ~Particle() = 0;

	std::string GetName(void) const { return name_; }
	unsigned int GetMass(void) const { return mass_ }
	double GetWidth(void) const { return width_; } //in eV
	double GetHalfLife(void) const {
		return width_ > 0 ? gSettings.PhysConsts()->h_bar_eVs/width_ : DBL_MAX;
	} //in seconds

	unsigned int GetProcessSize(const Partice *target) const {
		auto procs = processes_.find(target->GetName();
		if (processes_.end() != procs)
			return procs->second.size();
		return 0;
	}
	std::string GetProcessName(const Partice *target, unsigned int process) const {
		auto procs = processes_.find(target->GetName();
		if (processes_.end() != procs) {
			auto pr = procs->second.find(process);
			if (procs->second.end()!=pr)
				return pr->second;
		}
		return "";
	}
	bool isValid(void) const { return is_valid_; }
	
	virtual unsigned int GetQauntStateSize(const Partice *target, double E, double theta, unsigned int process) const = 0;
	virtual double GetCrossSection(const Partice *target, double E, unsigned int process) const = 0;
	virtual double GetCrossSection(const Partice *target, double E, double theta, unsigned int process) const = 0;
	virtual std::vector<Particle*> GetFinalStates(const Partice *target, double E, double theta, unsigned int process) const = 0;

	virtual unsigned int GenerateProcess(const Partice *target, double E, double theta, double Rand) const = 0;
	virtual double GenerateScatterAngle(const Partice *target, double E, double theta, unsigned int process, double Rand) const = 0;
	//For neutral bremsstrahlung or deexcitation
	virtual double GeneratePhoton(const Partice *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateTimeDelay(const Partice *target, double E, double theta, unsigned int process, double Rand) const = 0;

	//Untabulated functions:
	virtual unsigned int GenerateUntabProcess(const Partice *target, double E, double theta, double Rand) const = 0;
	virtual double GenerateUntabScatterAngle(const Partice *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateUntabPhoton(const Partice *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateUntabTimeDelay(const Partice *target, double E, double theta, unsigned int process, double Rand) const = 0;
}

//Call Normalize() after setting up compounds
class Mixture {
protected:
	bool is_valid_;
	std::map<std::string, double> components_;
	std::string name_;
public:
	Mixture(std::string mixture_name) {
		is_valid_ = false;
		name_ = mixture_name;
	}
	virtual ~Mixture() = 0;

	bool isValid(void) const { return is_valid_; }
	std::string GetName(void) const { return name_; }
	void SetName(std::string name) { name_ = name; }
	void Normalize(void) {
		while (true) {
			for (auto i = components_.begin(), i_end_ = components_.end(); i != i_end_; ++i)
				if (i->second <= 0) {
					std::cerr << "Mixture(\"" << name_ << "\")::Normalize: Warning! Negative or zero concentration of \"" \
						<< i->first << "\" is removed." << std::endl;
					components_.erase(i);
					continue;
				}
			break;
		}
		double Norm = 0;
		for (auto i = components_.begin(), i_end_ = components_.end(); i != i_end_; ++i)
			Norm += i->second;
		for (auto i = components_.begin(), i_end_ = components_.end(); i != i_end_; ++i)
			i->second /= Norm;
		is_valid_ = !components_.empty();
	}
	
	virtual double GetCrossSection(const Partice *incident, const Partice *incident, double E) const = 0;
	virtual long double GetXSIntegral(const Partice *incident, long double En_from, long double En_to, long double En_y) const = 0;
	virtual long double GetEnergyFromXSIntegral(const Partice *incident, long double En_y, long double integral_value) const = 0;
	virtual const Partice* GenerateScatteringParticle(const Partice *incident, double E, double Rand) const = 0;

	//Untabulated functions:
	virtual double GetUntabCrossSection(const Partice *incident, double E) const = 0;
	virtual long double GetUntabXSIntegral(const Partice *incident, long double En_from, long double En_to, long double En_y) const = 0;
	virtual long double GetUntabEnergyFromXSIntegral(const Partice *incident, long double En_y, long double integral_value) const = 0;
	virtual const Partice* GenerateUntabScatteringParticle(const Partice *incident, double E, double Rand) const = 0;
}

#endif //PARTICLE_H_
