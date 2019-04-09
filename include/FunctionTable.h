#ifndef FUNCTION_TABLE_H
#define FUNCTION_TABLE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include "global_definitions.h"
#include "PolynomialFit.h"

//find_E_indexes_by_value works only for monotone rising function!
//F(E1,Ey)>F(E2,Ey) if E1>E2
//Ey:[0, Emax = 16eV]
//E: [Ey,Emax=16eV] - not necessary - for theta tables E===theta [0, M_PI] and F - probability function
//F(Ey,Ey) = 0; - not necessary condition
class FunctionTable {
protected:
	std::deque<std::vector<double> > _Es;
	std::deque<std::vector<double> > _ys;
	std::deque<double> _Eys;
	std::pair<long int, long int> find_Ey_indexes (double Ey) const;
	std::pair<long int, long int> find_E_indexes (double E, std::size_t Ey_index) const;
	std::pair<long int, long int> find_E_indexes_by_value (double val, std::size_t Ey_index) const;
public:
	FunctionTable(void);
	virtual double operator ()(double E, double Ey) const;
	double find_E (double Ey, double val) const;
	void push (double E, double Ey, double val);
	void clear (void);

	void read (std::ifstream& str);
	void write(std::string fname);
	void write (std::ofstream& str);
	void plot_E_Ey (void);
	bool is_empty(void) const;
};

#endif
