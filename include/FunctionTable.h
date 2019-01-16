#ifndef FUNCTION_TABLE_H
#define FUNCTION_TABLE_H

#include <iostream>
#include <vector>
#include <deque>
#include "PolynomialFit.h"

//for monotone rising function only!
//F(E1,Ey)>F(E2,Ey) if E1>E2
//Ey:[0, Emax = 16eV]
//E: [Ey,Emax=16eV]
//F(Ey,Ey) = 0; - not necessary condition
class FunctionTable {
protected:
	std::deque<std::vector<double> > _Es;
	std::deque<std::vector<double> > _ys;
	std::deque<double> _Eys;
	std::pair<long int, long int> find_Ey_indexes (double Ey);
	std::pair<long int, long int> find_E_indexes (double E, std::size_t Ey_index);
	std::pair<long int, long int> find_E_indexes_by_value (double val, std::size_t Ey_index);
public:
	FunctionTable(void);
	virtual double operator ()(double E, double Ey);
	double find_E (double Ey, double val);
	void push (double E, double Ey, double val);
	void clear (void);

	void read (std::ifstream& str);
	void write (std::ofstream& str);
	bool is_empty(void);
};

#endif
