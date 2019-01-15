#ifndef FUNCTION_TABLE_H
#define FUNCTION_TABLE_H

#include <iostream>
#include <vector>
#include <deque>
#include "PolynomialFit.h"

//parameters are [0]+[1]*x+[2]*x^2+...
class FunctionTable {
protected:
	std::deque<std::vector<double> > _Es;
	std::deque<std::vector<double> > _ys;
	std::deque<double> _Eys;
	std::pair<long int, long int> find_Ey_indexes (double Ey);
	std::pair<long int, long int> find_E_indexes (double E, std::size_t Ey_index);
public:
	FunctionTable(void);
	virtual double operator ()(double E, double Ey);
	double find_E (double Ey, double val);
	double find_Ey (double E, double val);
	void push (double E, double Ey, double val);
	void clear (void);
};

#endif
