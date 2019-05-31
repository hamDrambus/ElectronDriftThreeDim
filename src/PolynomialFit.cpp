#include "PolynomialFit.h"

namespace uBLAS = boost::numeric::ublas;

PolynomialFit::PolynomialFit(int order)
{
	setOrder(order);
}

PolynomialFit::~PolynomialFit()
{}

void PolynomialFit::setOrder(int n)
{
	if (n < 0)
		_order = 2;
	else
		_order = n;
}

int PolynomialFit::getOrder(void) const
{
	return _order;
}

std::vector<double> PolynomialFit::operator ()(const std::vector<double> &xs_in, const std::vector<double> &ys_in, boost::optional<double> &in_x0) const {
	return (*this)(xs_in, ys_in, 0, xs_in.size(), in_x0);
}

std::vector<double> PolynomialFit::operator ()(const std::vector<double> &xs_in, const std::vector<double> &ys_in,
	int offset, int N_points, boost::optional<double> &in_x0) const//only for a part of a vector
{
	std::vector<double> out;
	if (xs_in.size() != ys_in.size()) {
		std::cout<<"PolynomialFit::operator(): Error: x-y data size mismatch"<<std::endl;
		return out;
	}
	if ((xs_in.size()-offset) < N_points) {
		std::cout<<"PolynomialFit::operator(): Error: N points is out of range:"<<std::endl;
		std::cout<<"\tx.size()="<<xs_in.size()<<" offset="<<offset<<" N_points="<<N_points<<std::endl;
		return out;
	}
	if (offset < 0) {
		std::cout<<"PolynomialFit::operator(): Error: offset is out of range:"<<std::endl;
		std::cout<<"\tx.size()="<<xs_in.size()<<" offset="<<offset<<" N_points="<<N_points<<std::endl;
		return out;
	}
	if (N_points < (_order + 1)) {
		std::cout<<"PolynomialFit::operator(): Error: no enough N points for fit:"<<std::endl;
		std::cout<<"\torder="<<_order<<" N_points="<<N_points<<std::endl;
		return out;
	}
	in_x0 = (in_x0 ? *in_x0 : xs_in[offset]); //It is bad to set x0 to some fixed value (e.g. 0) because
	//interpolating too far from it will result in unstable results due to limited precision.
	//Ideally x0 should be set to the point at which we interpolate the data.
	if (1 == _order) {
		out.resize(2);
		out[1] = (ys_in[offset + 1] - ys_in[offset]) / (xs_in[offset + 1] - xs_in[offset]);
		out[0] = ys_in[offset] + (*in_x0 - xs_in[offset])*out[1];
		//^value at in_x0 point
	} else {
		uBLAS::matrix<double> mat(N_points, _order + 1);
		for (int col = 0, col_end_ = mat.size2(); col < col_end_; ++col)
			for (int row = 0, row_end_ = mat.size1(); row < row_end_; ++row)
				mat(row, col) = pow(xs_in[offset + row] - *in_x0, col);
		uBLAS::vector<double> Y(N_points);
		for (int row = 0, row_end_ = Y.size(); row < row_end_; ++row)
			Y[row] = ys_in[offset + row];
		//Solve the equation mat^T*mat*X = mat^T*Y for X via LU decomposition (mat is generally not diagonal)
		Y = uBLAS::prod(uBLAS::trans(mat), Y);
		mat = uBLAS::prod(uBLAS::trans(mat), mat);
		int res = uBLAS::lu_factorize(mat);
		if (res != 0)
			return out;
		uBLAS::inplace_solve(mat, Y, uBLAS::unit_lower_tag());
		uBLAS::inplace_solve(mat, Y, uBLAS::upper_tag());
		out.resize(Y.size());
		std::copy(Y.begin(), Y.end(), out.begin());
	}
	if (out.size() != (_order+1)) {
		out.resize(0);
		return out;
	}
	return out;
}

//=========================================================

DataVector::DataVector(int fit_order, int N_used_) :
		fitter(fit_order), use_left(false), use_right(false), is_set_left(false),
		is_set_right (false), left_value(DBL_MAX), right_value(DBL_MAX)
{
	N_used = std::max(1, N_used_);
}

DataVector::DataVector(std::vector < double> &xx, std::vector<double> &yy, int fit_order, int N_used_): DataVector(fit_order, N_used_)
{
	initialize(xx, yy, fit_order, N_used_);
}
DataVector::~DataVector() {}

void DataVector::initialize(std::vector < double> &xx, std::vector<double> &yy, int fit_order, int N_used_)
{
	N_used = std::max(1, N_used_);
	fitter.setOrder(fit_order);
	if (xx.size()!=yy.size()) {
		std::cout<<"DataVector::initialize(): Error: x and y data size mismatch!"<<std::endl;
		return;
	}
	xs = xx;
	ys = yy;
	double mem;
	for (long unsigned int i=0, i_end = std::min(xs.size(),ys.size());i!=i_end;++i) { //sorting ys(xs) by xs
		auto min_index = i;
		double min_val = xs[i];
		for (auto j=i;j!=i_end;++j)
			if (xs[j]<min_val) {
				min_val = xs[j];
				min_index = j;
			}
		if (min_index!=i) {
			mem = xs[min_index];
			xs[min_index] = xs[i];
			xs[i] = mem;
			mem = ys[min_index];
			ys[min_index] = ys[i];
			ys[i] = mem;
		}
	}
}

void DataVector::setOrder(int ord)
{	fitter.setOrder(ord); }

int DataVector::getOrder(void) const
{	return fitter.getOrder(); }

void DataVector::setNused(int N)
{	N_used = std::max(1, N); }

int DataVector::getNused(void) const
{	return N_used; }

void DataVector::use_leftmost(bool use)
{
	use_left = use;
}
void DataVector::use_rightmost(bool use)
{
	use_right = use;
}
void DataVector::set_leftmost(double val)
{
	is_set_left = true;
	left_value = val;
}
void DataVector::unset_leftmost(void)
{
	is_set_left = false;
}
void DataVector::set_rightmost(double val)
{
	is_set_right = true;
	right_value = val;
}
void DataVector::unset_rightmost(void)
{
	is_set_right = false;
}
void DataVector::set_out_value(double val)
{
	set_leftmost(val);
	set_rightmost(val);
}
void DataVector::unset_out_value(void)
{
	unset_leftmost();
	unset_rightmost();
}

void DataVector::push (double x, double y) //do not disrupt order
{
	for (auto i=xs.begin(), i_end = xs.end(), j=ys.begin(), j_end = ys.end(); (i!=i_end)&&(j!=j_end);++i, ++j) {
		if (x==*i)
			return; //do not insert points with equal x
		if (*i>x) {
			xs.insert(i,x);
			ys.insert(j,y);
			return;
		}
	}
	xs.insert(xs.end(),x);
	ys.insert(ys.end(),y);
}

void DataVector::push_back (double x, double y)//does not check that the new array is ordered, but faster
{
	xs.push_back(x);
	ys.push_back(y);
}

void DataVector::erase (std::size_t n)
{
	xs.erase(xs.begin()+n);
	ys.erase(ys.begin()+n);
}

void DataVector::clear (void)
{
	xs.clear();
	ys.clear();
}

std::size_t DataVector::size (void) const
{	return std::min(xs.size(), ys.size());}

double DataVector::getX(std::size_t n) const
{	return xs[n]; }

double DataVector::getY(std::size_t n) const
{	return ys[n]; }

void DataVector::read(std::ifstream& str) //TODO: add try/catch for handling stoi and stod
{
	clear();
	std::string line, word;
	int line_n = 0;
	while (!str.eof() && str.is_open()) {
		std::getline(str, line);
		++line_n;
		if (1 == line_n) {
			//parse "//Order	N_used	use_left use_right is_set_left is_set_right left_value right_value"
			double dval;
			int ival;
			int word_n = 0;
			if (line.size() < 2)
				goto wrong_format;
			if ((line[0] != '/') || (line[1] != '/'))
				goto wrong_format;
			line.erase(line.begin(), line.begin() + 2);
			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			ival = std::stoi(word);
			setOrder(ival);

			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			ival = std::stoi(word);
			setNused(ival);
			
			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			ival = std::stoi(word);
			use_leftmost(ival);

			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			ival = std::stoi(word);
			use_rightmost(ival);

			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			ival = std::stoi(word);
			is_set_left = ival;

			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			ival = std::stoi(word);
			is_set_right = ival;

			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			dval = std::stoi(word);
			left_value = dval;

			word = strtoken(line, "\t ");
			++word_n;
			if (word.empty())
				goto wrong_format;
			dval = std::stoi(word);
			right_value = dval;
			continue;
		wrong_format:
			std::cerr << "DataVector::read: Error on line " << line_n << " word "<<word_n<<": wrong header format" << std::endl;
			return;
		}
		if (line.size() >= 2) //Ignore simple c style comment
			if ((line[0] == '/') && (line[1] == '/'))
				continue;
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double x = std::stod(word);
		word = strtoken(line, "\t ");
		if (word.empty())
			break;
		double val = std::stod(word);
		push_back(x, val);
	}
}

void DataVector::write(std::string fname, std::string comment) const
{
	std::ofstream str;
	str.open(fname, std::ios_base::trunc);
	write(str, comment);
	str.close();
}

void DataVector::write(std::ofstream& str, std::string comment) const
{
	//"//Order	N_used	use_left use_right is_set_left is_set_right left_value right_value"
	str << "//" << getOrder() << "\t" << N_used << "\t" << (use_left ? 1 : 0) << "\t" << (use_right ? 1 : 0)
		<< "\t" << (is_set_left ? 1 : 0) << "\t" << (is_set_right ? 1 : 0) << "\t" << left_value << "\t" << right_value << std::endl;
	str << "//" << comment << std::endl;
	for (std::size_t i = 0, i_end_ = size(); i != i_end_; ++i) {
		str << xs[i] << "\t" << ys[i] << std::endl;
	}
}

double DataVector::operator()(double point, boost::optional<double> x0) const
{
	if (xs.empty()||ys.empty())
		return DBL_MAX;
	if (point < xs.front()) {
		if (use_left)
			return ys.front(); //TODO: maybe add scaling
		if (is_set_left)
			return left_value;
	}
	if (point > xs.back()) {
		if (use_right)
			return ys.back();
		if (is_set_right)
			return right_value;
	}
	int n_min, n_max;
	get_indices(point, n_min, n_max);
	std::vector<double> coefs = fitter(xs, ys, n_min, n_max-n_min+1, x0); //n_max-n_min+1==N_used
	if (0!=coefs.size())
		return calculate(point, *x0, coefs);
	return DBL_MAX;
}

//[n_min, n_max] are used, not [n_min,n_max). N_used==n_max-n_min+1>=order+1
void DataVector::get_indices(double x, int &n_min, int &n_max) const
{
	if (x <= xs.front())
	{
		n_min = 0;
		n_max = N_used-1; //order is set to be matching the number of points
		return;
	}
	if (x > xs.back())
	{
		n_max = xs.size() - 1;
		n_min = n_max - N_used + 1;
		return;
	}
	std::size_t left = 0;
	std::size_t right = std::max(left, xs.size() - 2);
	std::size_t out_ = (left+right)/2;
	while (true) {
		if ((x > xs[out_]) && (x <= xs[out_ + 1]))
			break;
		if (x <= xs[out_]) {
			right = out_ - 1;
			out_ = (left + right) / 2;
		} else {
			left = out_ + 1;
			out_ = (left + right) / 2;
		}
	}
	//out_ is such, that x>xs[out] and x<=xs[out+1]
	n_min = out_ - (N_used-1) / 2; //asymmetrical interpolation in the case of odd order.
	n_max = n_min + (N_used-1);
	if (n_min < 0)
	{
		n_min = 0;
		n_max = n_min + (N_used-1);
	}
	if (n_max >= xs.size())
	{
		n_max = xs.size() - 1;
		n_min = n_max - N_used+1;
	}
}

double DataVector::calculate(double x, double x0, const std::vector<double>& coefs) const
{
	std::size_t order = coefs.size();
	double out_ = 0;
	for (std::size_t o_O = 0; o_O < order; ++o_O)
		out_ += coefs[o_O]*pow(x-x0, o_O);
	return out_;
}

PDF_routine::PDF_routine(): DataVector(1, 2)
{
	set_out_value(0);
	//TODO: same for cdf? left=0, right =1
}

PDF_routine::PDF_routine(std::vector < double> &pdf_xx, std::vector<double> &pdf_yy)
{
	initialize(pdf_xx, pdf_yy, 1, 2);
	pdf_to_cdf();
}

bool PDF_routine::pdf_to_cdf(void)
{
	cdf_ys.resize(xs.size());
	double X_prev = xs.front();
	double Y_prev = 0;
	for (long int i = 0, i_end_ = cdf_ys.size(); i != i_end_; ++i) {
		double X = xs[i];
		double Y = ys[i];
		cdf_ys[i] = 0;
		if (i != 0) {
			cdf_ys[i] = cdf_ys[i-1] + 0.5*(Y + Y_prev)*(X - X_prev);//Integral
		}
		X_prev = X;
		Y_prev = Y;
	}
	for (long int i = 0, i_end_ = cdf_ys.size(); i != i_end_; ++i)
		cdf_ys[i] /= cdf_ys[i_end_ - 1];//normalize probability function
	return true;
}

PDF_routine::~PDF_routine()
{}

bool PDF_routine::read(std::string& fname)
{
	std::ifstream str;
	str.open(fname);
	if (!str.is_open() || !str.good()) {
		std::cout << "PDF_routine::read:: Warning: Could not open file \"" << fname<<"\"" << std::endl;
		return false;
	}
	read(str);
	if (size() < getNused()) {
		std::cout << "PDF_routine::read:: Warning: file \"" << fname << "\" contains too little data" << std::endl;
		return false;
	}
	return true;
}

void PDF_routine::read(std::ifstream& str)
{
	clear();
	std::vector<double> ixx, iyy;//it's better to push whole read array instead of element-wise filling
	std::string line, word;
	int line_n = 0;
	while (!str.eof() && str.is_open()) {
		std::getline(str, line);
		++line_n;
		if (1 == line_n) {
			try { //parse "//Order	N_used	use_left use_right is_set_left is_set_right left_value right_value"
				double dval;
				int ival;
				int word_n = 0;
				if (line.size() < 2)
					throw std::invalid_argument("Header line expected \"//Order\\tNused\\t...\"");
				if ((line[0] != '/') || (line[1] != '/'))
					throw std::invalid_argument("Header line expected \"//Order\\tNused\\t...\"");
				line.erase(line.begin(), line.begin() + 2);
				word = strtoken(line, "\t ");
				++word_n;
				ival = std::stoi(word);
				setOrder(ival);

				word = strtoken(line, "\t ");
				++word_n;
				ival = std::stoi(word);
				setNused(ival);

				word = strtoken(line, "\t ");
				++word_n;
				ival = std::stoi(word);
				use_leftmost(ival);

				word = strtoken(line, "\t ");
				++word_n;
				ival = std::stoi(word);
				use_rightmost(ival);

				word = strtoken(line, "\t ");
				++word_n;
				ival = std::stoi(word);
				is_set_left = ival;

				word = strtoken(line, "\t ");
				++word_n;
				ival = std::stoi(word);
				is_set_right = ival;

				word = strtoken(line, "\t ");
				++word_n;
				dval = std::stoi(word);
				left_value = dval;

				word = strtoken(line, "\t ");
				++word_n;
				dval = std::stoi(word);
				right_value = dval;
				continue;
			} catch (std::exception &e) {
				setOrder(1);
				setNused(2);
				use_leftmost(false);
				use_rightmost(false);
				set_out_value(0); //out of range value
				continue;
			}
		}
		if (line.size() >= 2) //Ignore simple c style comment
			if ((line[0] == '/') && (line[1] == '/'))
				continue;
		double val, x;
		try {
			word = strtoken(line, "\t ");
			x = std::stod(word);
			word = strtoken(line, "\t ");
			val = std::stod(word);
		} catch (std::exception &e) {
			continue;
		}
		push_back(x, val);
		ixx.push_back(x);
		iyy.push_back(val);
	}
	initialize(ixx, iyy, 1, 2);
	pdf_to_cdf();
}

double PDF_routine::generate(double Rand) const //has no internal random engine
{
	for (std::size_t i = 0, i_end_ = cdf_ys.size(); i != i_end_; ++i) {
		if (Rand <= cdf_ys[i]) {
			if (0 == i)
				return xs[i];
			double x0 = xs[i-1], x1 = xs[i];
			double y0 = cdf_ys[i-1], y1 = cdf_ys[i];
			return x0 + (x1 - x0)*(Rand - y0) / (y1 - y0);
		}
	}
	return  (xs.size() ? xs.back() : 0);
}
