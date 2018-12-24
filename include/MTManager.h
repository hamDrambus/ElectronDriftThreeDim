#ifndef MTMANAGER_H_
#define MTMANAGER_H_

#include <TMutex.h>
#include <TThread.h>
#include "Manager.h"

class MTManager: public Manager
{
protected:
	TCondition* condition_;
	TMutex* thread_mutex_;
	int instance_;
	int N_electrons_;
public:
	MTManager(ArDataTables *Ar_tables, int instance, int N_electrons, UInt_t RandomSeed = 42);
	void ProcessAll(void);
	void setCondition(TCondition* cond);
	TCondition* getCondition(void) const;
	void setThreadMutex(TMutex* mutex);
	TMutex* getThreadMutex(void) const;
	void Merge(MTManager *with);
};

#endif 

