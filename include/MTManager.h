#ifndef MTMANAGER_H_
#define MTMANAGER_H_

#include <thread>
#include "Manager.h"

class MTManager: public Manager
{
protected:
	int instance_;
	boost::optional<unsigned int> N_electrons_;
public:
	MTManager(const Mixture *material, int instance);
	void ProcessAll(void);
	bool setNelectons(unsigned int Ne);
	boost::optional<unsigned int> getNelectons(void) const;
	void Merge(MTManager *with);
	bool isReady(void) const;
	void Clear(void);
};

#endif 

