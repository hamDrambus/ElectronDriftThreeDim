#ifndef POLYNOMIAL_FIT_H
#define POLYNOMIAL_FIT_H

#include <iostream>
#include <fstream>
#include <string>
#include <cfloat>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include "global_definitions.h"

//TVectorD parameters are [0]+[1]*x+[2]*x^2+...
class PolynomialFit {
protected:
	int _order;
public:
	PolynomialFit(int order);
	virtual ~PolynomialFit();
	virtual void setOrder(int n); //TODO: actually is is a bad practice to call virtual method from the constructor
	//but it is ok here, since derivative class only limits setOrder() possible values to {2}

	int getOrder(void) const;

	//in_x0 - in what poInt_t set zero x (In the SG filter it is convenient to set x_in
	//to the poInt_t in which value is calculated
	virtual std::vector<double> operator ()(const std::vector<double> &xs_in, const std::vector<double> &ys_in, boost::optional<double> &in_x0) const;
	//Fit only part of a vector. offset+N_points-1 must in the range of the vector
	virtual std::vector<double> operator ()(const std::vector<double> &xs_in, const std::vector<double> &ys_in,
		int offset, int N_points, boost::optional<double> &in_x0) const;
};

//Wraps PolynomialFit: stores raw data, N points used in every fit and last region (cache_n_from, cache_n_to) in which fit/interpolation took place.
//The latter is required for optimization, because polynomial coefficients are not updated unless necessary (x moved from previous region)
class DataVector { //TODO: cache is useless in practice, but affects multithreading, there is need to store separate copies of DataVector for each thread.
	//TODO: either remove altogether, or event better determine automatically when its usage is correct (make this optional as well)
protected:
	std::vector<double> xs;
	std::vector<double> ys;
	PolynomialFit fitter;
	int N_used;

	bool use_left, use_right, is_set_left, is_set_right;
	double left_value, right_value;
public:
	DataVector(int fit_order = 1, int N_used = 2);
	DataVector(std::vector < double> &xx, std::vector<double> &yy, int fit_order, int N_used);
	virtual ~DataVector();

	void initialize(std::vector < double> &xx, std::vector<double> &yy,  int fit_order, int N_used);
	void setOrder(int ord);
	int getOrder(void) const;
	void setNused(int ord);
	int getNused(void) const;

	//precedence goes to use_left-/right-most methods.
	void use_leftmost(bool use);
	void use_rightmost(bool use);
	void set_leftmost(double val);
	void unset_leftmost(void);
	void set_rightmost(double val);
	void unset_rightmost(void);
	void set_out_value(double val); //out of range value
	void unset_out_value();

	double operator()(double point, boost::optional<double> x0 = boost::none) const; //x0 = point is recommended to use. At least x0 must be close to point, or there will be large errors otherwise
	void push (double x, double y);
	void push_back (double x, double y);
	void erase (std::size_t n);
	void clear (void);
	std::size_t size (void) const;
	double getX(std::size_t n) const;
	double getY(std::size_t n) const;

	//save/load full state except cache from file
	virtual void read(std::ifstream& str);
	void write(std::string fname, std::string comment = "") const;
	void write(std::ofstream& str, std::string comment="") const;
protected:
	void get_indices(double point, int &n_min, int &n_max) const; //[n_min, n_max] are used, not [n_min,n_max). N_used==n_max-n_min+1>=order+1
	double calculate(double x, double x0, const std::vector<double>& coefs) const;
};

//TODO: create common parent for PDF_routine and DataVector
//does not support indefinite domain.
//!TODO: use boost or some libraries like normal ppl.
class PDF_routine : public DataVector { //Probability Density Function routine - load not normalized distribution from file,
	// construct cumulative DF and use it to generate values.
	//TODO: calling of set/getOrder set/getNused is forbidden. PDF is fixed for order =1, but in can be generalized at the cost of performance.
	//TODO: use 3rd party code? CERN's ROOT or GSL has relevant methods.
protected:
	bool pdf_to_cdf(void);
public:
	PDF_routine();
	PDF_routine(std::vector < double> &pdf_xx, std::vector<double> &pdf_yy);
	virtual ~PDF_routine();
	virtual void read(std::ifstream& str);
	virtual bool read(std::string& fname);
	//std::vector<double> cdf_xs; == DataVector::xs
	std::vector<double> cdf_ys;
	double generate(double Rand) const; //has no internal random engine
};

#endif
