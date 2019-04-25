#ifndef POLYNOMIAL_FIT_H
#define POLYNOMIAL_FIT_H

#include <iostream>
#include <fstream>
#include <string>
#include <TMath.h>
#include "TMatrixD.h"
#include "TVectorD.h"
#include "global_definitions.h"

//parameters are [0]+[1]*x+[2]*x^2+...
class PolynomialFit { //TODO: storing _last_coefs is useless. (affects multithreading, there is need to store separate copies of DataVector for each thread)
protected:
	Int_t _order;
	TVectorD _last_coefs;
public:
	PolynomialFit(Int_t order);
	virtual ~PolynomialFit();
	virtual void setOrder(Int_t n); //TODO: actually is is a bad practice to call virtual method from the constructor
	//but is is ok here, since derivative class only limits setOrder() possible values to {2}

	Int_t getOrder(void) const;
	void getCoefs(TVectorD &pars) const;

	virtual Int_t operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in, double in_x0=0); //in_x0 - in what poInt_t set zero x (In the SG filter it is convenient to set x_in
	//to the poInt_t in which value is calculated
	virtual Int_t operator ()(std::vector<double> &xs_in, std::vector<double> &ys_in,
		int offset, int N_points, double in_x0=0); //only for a part of a vector. offset+N_points-1 must in the range of vector
};

//Wraps PolynomialFit: stores raw data, N points used in every fit and last region (cache_n_from, cache_n_to) in which fit/interpolation took place.
//The latter is required for optimization, because polynomial coefficients are not updated unless necessary (x moved from previous region)
class DataVector { //TODO: cache is useless in practice, but affects multithreading, there is need to store separate copies of DataVector for each thread.
	//TODO: either remove altogether, or event better determine automatically when its usage is correct (make this optional as well)
protected:
	std::vector<double> xs;
	std::vector<double> ys;
	PolynomialFit fitter;
	double x0_used;
	Int_t N_used;
	Int_t cache_n_from, cache_n_to;
	bool isCached;

	bool use_left, use_right, is_set_left, is_set_right;
	double left_value, right_value;
public:
	DataVector(Int_t fit_order = 1, Int_t N_used = 2);
	DataVector(std::vector < double> &xx, std::vector<double> &yy, Int_t fit_order, Int_t N_used);
	virtual ~DataVector();

	void initialize(std::vector < double> &xx, std::vector<double> &yy,  Int_t fit_order, Int_t N_used);
	void setOrder(Int_t ord);
	Int_t getOrder(void);
	void setNused(Int_t ord);
	Int_t getNused(void);
	std::vector<double> getCoefs(void);

	//precedence goes to use_left-/right-most methods.
	void use_leftmost(bool use);
	void use_rightmost(bool use);
	void set_leftmost(double val);
	void unset_leftmost(void);
	void set_rightmost(double val);
	void unset_rightmost(void);
	void set_out_value(double val); //out of range value
	void unset_out_value();

	double operator()(double point);
	double operator()(double point, double x0); //x0 = point is recommended to use. At least x0 must be close to point, or there will be large errors otherwise
	void push (double x, double y);
	void push_back (double x, double y);
	void erase (std::size_t n);
	void clear (void);
	std::size_t size (void);
	double getX(std::size_t n);
	double getY(std::size_t n);

	//save/load full state except cache from file
	virtual void read(std::ifstream& str);
	void write(std::string fname, std::string comment = "");
	void write(std::ofstream& str, std::string comment="");
protected:
	void get_indices(double point, int &n_min, int &n_max); //[n_min, n_max] are used, not [n_min,n_max). N_used==n_max-n_min+1>=order+1
	double calculate(double x);
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
