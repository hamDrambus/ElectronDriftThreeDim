#ifndef MIXTURE_H_
#define MIXTURE_H_

#include "global_definitions.h"
#include "Particle.h"
#include "FunctionTable.h"
#include "ParticleTable.h"

class Mixture {
protected:
	bool is_valid_;
	std::vector<std::pair<const Particle*, double>> components_;
	std::string name_;
	std::string integral_table_fname;
	FunctionTable *integral_table_;
public:
	Mixture(std::string mixture_name);
	virtual ~Mixture();

	bool isValid(void) const {
		return is_valid_;
	}
	std::string GetName(void) const { return name_; }
	void SetName(std::string name) { name_ = name; }
	void AddComponent(const Particle* part, double fraction);
	void AddComponent(std::string name, double fraction);
	std::size_t GetNComponent(void) const {
		return components_.size();
	}
	const Particle* GetComponent(std::size_t index) const {
		if (index<0 || index >=components_.size())
			return NULL;
		return components_[index].first;
	}
	double GetComponentFraction(std::size_t index) const {
		if (index<0 || index >=components_.size())
			return 0;
		return components_[index].second;
	}
	const Particle* GetDominatingParticle(const Particle *incident, double E) const;
	void Prepare(void);

	virtual double GetCrossSection(const Particle *incident, double E) const;
	virtual long double GetXSIntegral(const Particle *incident, long double En_from, long double En_to, long double En_y) const;
	virtual long double GetEnergyFromXSIntegral(const Particle *incident, long double En_y, long double integral_value) const;
	virtual const Particle* GenerateScatteringParticle(const Particle *incident, double E, double Rand) const;
	//Untabulated functions:
	virtual double GetUntabCrossSection(const Particle *incident, double E) const;
	virtual long double GetUntabXSIntegral(const Particle *incident, long double En_from, long double En_to, long double En_y) const;
	virtual long double GetUntabEnergyFromXSIntegral(const Particle *incident, long double En_y, long double integral_value) const;
	virtual const Particle* GenerateUntabScatteringParticle(const Particle *incident, double E, double Rand) const;
protected:
	virtual bool generate_integral_table(void);
};

#endif //MIXTURE_H_
