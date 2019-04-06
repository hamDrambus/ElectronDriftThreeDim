#ifndef MTMANAGER_H_
#define MTMANAGER_H_

#include <thread>
#include "Manager.h"

class MTManager: public Manager
{
protected:
	int instance_;
	boost::optional<unsigned int> N_electrons_;
	boost::optional<std::size_t> run_index_;
public:
	MTManager(ArDataTables *Ar_tables, int instance);
	void ProcessAll(void);
	bool setRunIndex(std::size_t index);
	boost::optional<std::size_t> getRunIndex(void) const;
	bool setNelectons(unsigned int Ne);
	boost::optional<unsigned int> getNelectons(void) const;
	void Merge(MTManager *with);
	bool isReady(void) const;
	void Clear(void);
};

#endif 

