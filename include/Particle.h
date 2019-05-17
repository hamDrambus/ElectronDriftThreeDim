#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "global_definitions.h"
#include "Settings.h"
#include "ColoredInterval.h"

//TODO: this class and its daughters are required for generalization of scattering e in mixtures and potentially other particles as well
class Particle {
protected:
	//Particle type --> list of processes
	std::map<std::string, std::map<unsigned int, std::string> > processes_;
	std::string name_; //must be unique for every derived class
	double width_; //in eV
	double mass_; //in eV
	bool is_valid_;
	ColoredRange XS_En_sweeper_; //TODO: make different for each incident particle
public:
	Particle(void) {
		is_valid_ = false;
		name_ = "GeneralParticle";
		mass_ = 0;
		width_ = 0;
	}
	virtual ~Particle() {}

	std::string GetName(void) const { return name_; }
	unsigned int GetMass(void) const { return mass_; }
	double GetWidth(void) const { return width_; } //in eV
	double GetHalfLife(void) const {
		return width_ > 0 ? gSettings.PhysConsts()->h_bar_eVs/width_ : DBL_MAX;
	} //in seconds

	unsigned int GetProcessSize(const Particle *target) const {
		if (NULL == target)
			return 0;
		auto procs = processes_.find(target->GetName());
		if (processes_.end() != procs)
			return procs->second.size();
		return 0;
	}
	std::string GetProcessName(const Particle *target, unsigned int process) const {
		if (NULL == target)
			return "";
		auto procs = processes_.find(target->GetName());
		if (processes_.end() != procs) {
			auto pr = procs->second.find(process);
			if (procs->second.end()!=pr)
				return pr->second;
		}
		return "";
	}
	ColoredRange GetEnergies (void) const { return XS_En_sweeper_; }
	virtual bool isValid(void) const { return is_valid_; }
	
	virtual unsigned int GetQauntStateSize(const Particle *target, double E, double theta, unsigned int process) const = 0;
	virtual double GetCrossSection(const Particle *target, double E, unsigned int process) const = 0;
	virtual double GetCrossSection(const Particle *target, double E, double theta, unsigned int process) const = 0;
	virtual std::vector<const Particle*> GetFinalStates(const Particle *target, double E, double theta, unsigned int process) const = 0;

	//returns negative value in case of error
	virtual int GenerateProcess(const Particle *target, double E, double theta, double Rand) const {
		if (NULL == target) {
			std::cerr << GetName() << "::GenerateProcess: Error: NULL target"<<std::endl;
			return 0;
		}
		auto procs = processes_.find(target->GetName());
		if (processes_.end()==procs) {
			std::cerr << GetName() << "::GenerateProcess: Error: unsupported target particle \""<<target->GetName()<<"\""<<std::endl;
			return 0;
		}
		std::size_t n_procs = procs->second.size();
		std::vector<double> CrossSections(n_procs, 0.0), CrossSectionsSum(n_procs, 0.0);
		//TODO: allocating memory each time is quite expensive. Similar issue for Mixture. Need to create cache
		//All sizes and particle interations are static (so far), so vectors can be resized for each incident particle
		for (std::size_t ind = 0, ind_end_ = procs->second.size(); ind!=ind_end_; ++ind) {
			CrossSections[ind] = GetCrossSection(target, E, theta, ind);
			CrossSectionsSum[ind] = std::max(CrossSections[ind], 0.0) + ((ind==0) ? 0.0 : CrossSectionsSum[ind - 1]);
		}
		for (unsigned int ind = 0, ind_end_ = procs->second.size(); ind!=ind_end_; ++ind) {
			CrossSectionsSum[ind] /= CrossSectionsSum[ind_end_ - 1];
			if (Rand < CrossSectionsSum[ind]) {
				return ind;
			}
		}
		std::cerr << GetName() << "::GenerateProcess: Error: failed to select process from "<<n_procs<<" candidates"<<std::endl;
		return -1;
	}
	virtual double GenerateScatterAngle(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;
	//For neutral bremsstrahlung or deexcitation
	virtual double GeneratePhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;

	//Untabulated functions:
	virtual unsigned int GenerateUntabProcess(const Particle *target, double E, double theta, double Rand) const = 0;
	virtual double GenerateUntabScatterAngle(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateUntabEnergyLoss(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateUntabPhoton(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;
	virtual double GenerateUntabTimeDelay(const Particle *target, double E, double theta, unsigned int process, double Rand) const = 0;
};

#endif //PARTICLE_H_
