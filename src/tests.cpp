#include "tests.h"

void test_polynomial_fit (void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::string fname_raw = prefix + "sinx_to_x_raw.txt";
	std::string fname_fit1 = prefix + "sinx_to_x_fit1.txt";
	std::string fname_fit2 = prefix + "sinx_to_x_fit2.txt";
	std::string fname_fit3 = prefix + "sinx_to_x_fit3.txt";
	std::vector<double> xs, ys;
	xs.resize(50);
	ys.resize(50);
	std::ofstream str;
	open_output_file(fname_raw, str, std::ios_base::trunc);
	for (int i=0; i!=50;++i) { //not sorted
		double x = i*M_PI*3/100 + 0.02;
		xs[ i>=25 ? i-25 : i+25 ] = x;
		double y = sin(x)/x;
		ys[ i>=25 ? i-25 : i+25 ] = y;
		str<<boost::lexical_cast<std::string>(x)<<"\t"<<boost::lexical_cast<std::string>(y)<<std::endl;
	}
	str.close();
	std::vector<double> ys1 = ys;
	ys1.pop_back();
	//TESTING CONSTRUCTION:
	std::cout << "DataVector::insert test:" << std::endl;
	DataVector data0(1, 2);
	data0.insert(0, 1);
	data0.insert(0, 2);
	data0.insert(1, 2.5);
	data0.insert(2, 2.5);
	data0.insert(2, 2.5);
	std::cout << "DataVector must be [3]: "<<std::endl;
	std::cout << "X:\t0.0\t1.0\t2.0"<<std::endl;
	std::cout << "Y:\t2.0\t2.5\t2.5"<<std::endl;
	std::cout << "DataVector is ["<<data0.size()<<"]"<<std::endl;
	std::cout << "X:";
	for (std::size_t i = 0, i_end_ = data0.size(); i!=i_end_; ++i)
		std::cout<<"\t"<<data0.getX(i);
	std::cout << std::endl << "Y:";
	for (std::size_t i = 0, i_end_ = data0.size(); i!=i_end_; ++i)
		std::cout<<"\t"<<data0.getY(i);
	std::cout << std::endl;
	//TESTING PARAMETERS:
	std::cout<<"Testing DataVector parameters (order and n of used points)"<<std::endl;
	DataVector data1(3, 4);
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#1 Must throw error:"<<std::endl;
	data1.initialize(xs, ys1, 3, 4);
	std::cout<<"data(0.5)(xs, ys1, 3, 4): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#2 Must return DBL_MAX("<<boost::lexical_cast<std::string>(DBL_MAX)<<"):"<<std::endl;
	data1.initialize(xs, ys, 3, 3);
	std::cout<<"data(0.5)(xs, ys, 3, 3): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#3 Must return ~"<<boost::lexical_cast<std::string>(sin(0.5)/0.5)<<":"<<std::endl;
	data1.setNused(4);
	std::cout<<"data(0.5)(xs, ys, 3, 4): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#4 Must return DBL_MAX("<<boost::lexical_cast<std::string>(DBL_MAX)<<"):"<<std::endl;
	data1.initialize(xs, ys, 50, 51);
	std::cout<<"data(0.5)(xs, ys, 50, 51): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#5 Must return ~"<<boost::lexical_cast<std::string>(sin(0.5)/0.5)<<":"<<std::endl;
	data1.initialize(xs, ys, 49, 50);
	std::cout<<"data(0.5)(xs, ys, 49, 50): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#6 Must return ~"<<boost::lexical_cast<std::string>(sin(0.5)/0.5)<<":"<<std::endl;
	data1.setOrder(0);
	data1.setNused(1);
	std::cout<<"data(0.5)(xs, ys, 0, 1): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#7 Must return ~"<<boost::lexical_cast<std::string>(sin(0.5)/0.5)<<":"<<std::endl;
	data1.setOrder(0);
	data1.setNused(50);
	std::cout<<"data(0.5)(xs, ys, 0, 50): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;
	std::cout<<"#8 Must return DBL_MAX("<<boost::lexical_cast<std::string>(DBL_MAX)<<"):"<<std::endl;
	data1.setOrder(0);
	data1.setNused(51);
	std::cout<<"data(0.5)(xs, ys, 0, 51): "<<data1(0.5)<<std::endl;
	std::cout<<"-------------------------------------------------------"<<std::endl;

	//TESTING QUALITY:
	data1.setOrder(4);
	data1.setNused(5);
	open_output_file(fname_fit1, str, std::ios_base::trunc);
	for (int i = 0; i<150; ++i) {
		double x = i*M_PI*3/250;
		double y = data1(x, x);
		str<<boost::lexical_cast<std::string>(x)<<"\t"<<boost::lexical_cast<std::string>(y)<<std::endl;
	}
	str.close();
	data1.setOrder(6);
	data1.setNused(15);
	data1.set_out_value(0.5);
	open_output_file(fname_fit2, str, std::ios_base::trunc);
	for (int i = 0; i<150; ++i) {
		double x = i*M_PI*3/250;
		double y = data1(x, x);
		str<<boost::lexical_cast<std::string>(x)<<"\t"<<boost::lexical_cast<std::string>(y)<<std::endl;
	}
	str.close();
	str.close();
	data1.setOrder(2);
	data1.setNused(10);
	data1.unset_out_value();
	open_output_file(fname_fit3, str, std::ios_base::trunc);
	for (int i = 0; i<150; ++i) {
		double x = i*M_PI*3/250;
		double y = data1(x, x);
		str<<boost::lexical_cast<std::string>(x)<<"\t"<<boost::lexical_cast<std::string>(y)<<std::endl;
	}
	str.close();
	std::string name = prefix + "test_fit.sc";
	open_output_file(name, str, std::ios_base::trunc);
	str<<"plot \""<<fname_raw<<"\" u 1:2 title \"raw sin(x)/x\""<<std::endl;
	str<<"replot \""<<fname_fit1<<"\" u 1:2 title \"fit1\""<<std::endl;
	str<<"replot \""<<fname_fit2<<"\" u 1:2 title \"fit2\""<<std::endl;
	str<<"replot \""<<fname_fit3<<"\" u 1:2 title \"fit3\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_2_dim_table (void) //TODO: improve this test
{
	FunctionTable table;
	//F(x,y) = -0.5x + y. x:[0,2] y:[x,1]
	std::function<double(double, double)> F = [](double x, double y)->double  {
		double X = std::max(x, 0.0);
		X = std::min(X, 2.0);
		double Y = std::max(y, X);
		Y = std::min(Y, 1.0);
		return -0.5 * X + Y;
	};
	for (int ix = 0, ix_end_=11; ix!=ix_end_; ++ix) {
		for (int iy = 0, iy_end_ = ix_end_ - ix; iy!=iy_end_; ++iy) {
			double Ey = ix * 2.0/(ix_end_-1);
			double E = (iy_end_==1) ? Ey : (Ey + iy * (1.0-Ey)/(iy_end_ -1));
			table.push(Ey, E, F(Ey, E));
		}
	}


	std::cout<<"(X, Y):\t(0.0,0.0)\t(2.0,1,0)\t(0.0,1.0)\t(0.5,0.5)\t(0.2,0.7)"<<std::endl;
	std::cout<<"Table:\t"<<boost::lexical_cast<std::string>(table(0,0))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(table(1,2))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(table(1,0.0))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(table(0.5,0.5))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(table(0.7,0.2))<<std::endl;
	std::cout<<"Must be:"<<std::endl;
	std::cout<<"F(X,Y):\t"<<boost::lexical_cast<std::string>(F(0, 0))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(F(1, 2))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(F(1, 0.0))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(F(0.5, 0.5))<<"\t";
	std::cout<<boost::lexical_cast<std::string>(F(0.7, 0.2))<<std::endl;

	std::cout<<"Testing index search functions:"<<std::endl;
	FunctionTable table2;
	for (std::size_t i = 0; i!=155; ++i) {
		table2.push(0, i/10.0, std::pow(i/10.0, 2));
	}
	boost::optional<std::pair<std::size_t, std::size_t>> inds = table2.getY_data(0).getY_indices(1.01);
	std::cout<<"getY_indices(1.01) = ["<<inds->first<<", " << inds->second <<"]. Must be [10, 11]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(1.54);
	std::cout<<"getY_indices(1.54) = ["<<inds->first<<", " << inds->second <<"]. Must be [12, 13]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(-1);
	std::cout<<"getY_indices(-1) = ["<<inds->first<<", " << inds->second <<"]. Must be [0, 0]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(0.0);
	std::cout<<"getY_indices(0.0) = ["<<inds->first<<", " << inds->second <<"]. Must be [0, 0]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(0.005);
	std::cout<<"getY_indices(0.005) = ["<<inds->first<<", " << inds->second <<"]. Must be [0, 1]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(0.01);
	std::cout<<"getY_indices(0.01) = ["<<inds->first<<", " << inds->second <<"]. Must be [1, 1]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(237.16);
	std::cout<<"getY_indices(237.16) = ["<<inds->first<<", " << inds->second <<"]. Must be [154, 154]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(300);
	std::cout<<"getY_indices(300) = ["<<inds->first<<", " << inds->second <<"]. Must be [154, 154]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(1.54);
	std::cout<<"getY_indices(1.54) = ["<<inds->first<<", " << inds->second <<"]. Must be [12, 13]"<<std::endl;
	inds = table2.getY_data(0).getY_indices(200);
	std::cout<<"getY_indices(200) = ["<<inds->first<<", " << inds->second <<"]. Must be [141, 142]"<<std::endl;
}

void test_phase_shift_fit (void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::string fname_McEachran = prefix + "phase_shifts_McEachran_";
	std::string fname_MERT = prefix + "phase_shifts_MERT.txt";
	std::string fname_phase_fit = prefix + "phase_shifts_fit_exp.txt";

	std::ofstream str;
	const std::vector<DataVector> *ps_pos, *ps_neg;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	if (!argon->isValid()) {
		std::cerr<<"test_phase_shift_fit:: invalid particle, quitting test"<<std::endl;
		return;
	}
	ps_pos = &(argon->ArAllData_.ArExper_.phase_shifts_pos_);
	ps_neg = &(argon->ArAllData_.ArExper_.phase_shifts_neg_);
	for (unsigned int l=0, l_end_ = std::max(ps_pos->size(), ps_neg->size()); l!=l_end_; ++l) {
		std::string fname = fname_McEachran + std::to_string(l) + ".txt";
		open_output_file(fname, str, std::ios_base::trunc);
		str<<"E[eV]\tphase shift "<<l<<std::endl;
		if ((l<ps_pos->size())&&(l<ps_neg->size())) {
			double k_l_neg_written = -DBL_MAX;
			for (std::size_t i = 0, i_end_ = (*ps_pos)[l].size(); i!=i_end_; ++i) {//write in file either PS_l+ (E1) or PS_l- (E2) and when E1==E2, use their average.
				//ps_pos[l] and ps_neg[l] are supposed to be sorted by energy (k).
				double k = (*ps_pos)[l].getX(i);
				double PS_l = (*ps_pos)[l].getY(i);
				for (std::size_t j = 0, j_end_ = (*ps_neg)[l].size(); j!=j_end_; ++j) {
					if (((*ps_neg)[l].getX(j)<k)&&((*ps_neg)[l].getX(j)>k_l_neg_written)){
						str << boost::lexical_cast<std::string>(pow((*ps_neg)[l].getX(j)/gSettings.PhysConsts()->a_h_bar_2eM_e_SI, 2)) << "\t"
							<< boost::lexical_cast<std::string>((*ps_neg)[l].getY(j))<<std::endl;
						k_l_neg_written = (*ps_neg)[l].getX(j);
					}
					if ((*ps_neg)[l].getX(j)==k) {
						PS_l = (PS_l * (l+1) + l*(*ps_neg)[l].getY(j))/(2.0*l + 1); //weighted average!.
						break;
					}
				}
				str << boost::lexical_cast<std::string>(pow(k/ gSettings.PhysConsts()->a_h_bar_2eM_e_SI, 2))<< "\t"
					<< boost::lexical_cast<std::string>(PS_l)<<std::endl;
			}
		}
		if (l<ps_pos->size()) { //only ps_pos is present
			for (std::size_t i = 0, i_end_ = (*ps_pos)[l].size(); i!=i_end_; ++i) {
				str << boost::lexical_cast<std::string>(pow((*ps_pos)[l].getX(i)/ gSettings.PhysConsts()->a_h_bar_2eM_e_SI, 2)) << "\t"
					<< boost::lexical_cast<std::string>((*ps_pos)[l].getY(i))<<std::endl;
			}
		} else { //only ps_neg is present
			for (std::size_t i = 0, i_end_ = (*ps_neg)[l].size(); i!=i_end_; ++i) {
				str << boost::lexical_cast<std::string>(pow((*ps_neg)[l].getX(i)/ gSettings.PhysConsts()->a_h_bar_2eM_e_SI, 2)) << "\t"
					<< boost::lexical_cast<std::string>((*ps_neg)[l].getY(i)) <<std::endl;
			}
		}
		str.close();
	}

	int err;
	double E, k;
	{
		open_output_file(fname_phase_fit, str, std::ios_base::trunc);
		str<<"E[eV]\tphase shifts 0, 1, ... "<<std::endl;
		EnergyScanner EnRange(EnergyScanner::PlotDiffXS);
		while (true) {
			E = EnRange.Next(err);
			if (0!=err)
				break;
			k = sqrt(E)*gSettings.PhysConsts()->a_h_bar_2eM_e_SI;
			str << boost::lexical_cast<std::string>(E)<<"\t";
			for (std::size_t l = 0, l_end_ = std::max(ps_pos->size(), ps_neg->size()); l!=l_end_; ++l) {
				double PS_l;
				if ((l<ps_pos->size())&&(l<ps_neg->size())) {
					PS_l = ((*ps_pos)[l](k,k) * (l+1) + l*(*ps_neg)[l](k, k))/(2.0*l + 1); //weighted average!.;
				}
				if (l<ps_pos->size()) { //only ps_pos is present
					PS_l = (*ps_pos)[l](k,k);
				} else { 				//only ps_neg is present
					PS_l = (*ps_neg)[l](k,k);
				}
				str << boost::lexical_cast<std::string>(PS_l)<<"\t";
			}
			str << std::endl;
		}
		str.close();
	}

	{
		open_output_file(fname_MERT, str, std::ios_base::trunc);
		str<<"E[eV]\tphase shifts 0, 1, ... "<<std::endl;
		EnergyScanner EnRange(EnergyScanner::PlotDiffXS);
		while (true) {
			E = EnRange.Next(err);
			if (0!=err)
				break;
			if (E>gSettings.PhysConsts()->XS_el_En_thresold)
				break;
			k = sqrt(E)*gSettings.PhysConsts()->a_h_bar_2eM_e_SI;
			str << boost::lexical_cast<std::string>(E) << "\t";
			for (std::size_t l = 0, l_end = std::max(ps_pos->size(), ps_neg->size()); l!=l_end; ++l) {
				long double ps_p, ps_n;
				argon->ArAllData_.argon_phase_values_MERT5(k, l, ps_p, ps_n);
				str << boost::lexical_cast<std::string>(ps_p) << "\t";
			}
			str<<std::endl;
		}
		str.close();
	}

	for (std::size_t l = 0, l_end = std::max(ps_pos->size(), ps_neg->size()); l!=l_end; ++l) {
		std::string name = prefix + "test_phase_shift_" + std::to_string(l) + ".sc";
		open_output_file(name, str, std::ios_base::trunc);
		str<<"set logscale x"<<std::endl;
		str<<"plot '"<<fname_McEachran + std::to_string(l) + ".txt" <<"' u 1:2 title 'McEachran_"<<l<<"'"<<std::endl;
		str<<"replot '"<<fname_MERT <<"' u 1:"<<2+l<<" w line lc rgb \"#FF0000\" title 'MERT_"<<l<<"'"<<std::endl;
		str<<"replot '"<<fname_phase_fit <<"' u 1:"<<2+l<<" w line lc rgb \"#000000\" title 'fit_"<<l<<"'"<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}
}

void test_legendre_polynomial(void)
{
	LegendrePolynom Pl;
	std::cout<<"Pl(3,0.36) =\t "<<boost::lexical_cast<std::string>(Pl(0.36, 3))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.42336"<<std::endl;
	std::cout<<"Pl(4,0.36) =\t "<<boost::lexical_cast<std::string>(Pl(0.36, 4))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.0375168"<<std::endl;
	std::cout<<"Pl(20,0.36) =\t "<<boost::lexical_cast<std::string>(Pl(0.36, 20))<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.0542800664"<<std::endl;
	std::cout<<"Pl(50,0.5) =\t "<<boost::lexical_cast<std::string>(Pl(0.5, 50))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.031059099239"<<std::endl;
	std::cout<<"Pl(49,0.5) =\t "<<boost::lexical_cast<std::string>(Pl(0.5, 49))<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.086292778960940"<<std::endl;
	std::cout<<"Pl(48,0.5) =\t "<<boost::lexical_cast<std::string>(Pl(0.5, 48))<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.118866275929531"<<std::endl;
	std::cout<<"Pl(47,0.5) =\t "<<boost::lexical_cast<std::string>(Pl(0.5, 47))<<std::endl;
	std::cout<<"Wolfram alpha:\t 0.032013921114504"<<std::endl;
	std::cout<<"Pl(6,-1) =\t "<<boost::lexical_cast<std::string>(Pl(-1.0, 6))<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	std::cout<<"Pl(6,0) =\t "<<boost::lexical_cast<std::string>(Pl(0, 6))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.3125"<<std::endl;
	std::cout<<"Pl(6,1) =\t "<<boost::lexical_cast<std::string>(Pl(1, 6))<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	std::cout<<"Pl(7,-1) =\t "<<boost::lexical_cast<std::string>(Pl(-1.0, 7))<<std::endl;
	std::cout<<"Wolfram alpha:\t -1"<<std::endl;
	std::cout<<"Pl(7,0) =\t "<<boost::lexical_cast<std::string>(Pl(0, 7))<<std::endl;
	std::cout<<"Wolfram alpha:\t 0"<<std::endl;
	std::cout<<"Pl(7,1) =\t "<<boost::lexical_cast<std::string>(Pl(1, 7))<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	LegendrePolynom Pl2;
	std::cout<<"===Pl2==="<<std::endl;
	std::cout<<"Pl(7,-1) =\t "<<boost::lexical_cast<std::string>(Pl2(-1, 7))<<std::endl;
	std::cout<<"Wolfram alpha:\t -1"<<std::endl;
	std::cout<<"Pl(7,-1) =\t "<<boost::lexical_cast<std::string>(Pl2(-1, 7))<<std::endl;
	std::cout<<"Wolfram alpha:\t -1"<<std::endl;
	std::cout<<"Pl(8,-1) =\t "<<boost::lexical_cast<std::string>(Pl2(-1, 8))<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	std::cout<<"Pl(8,-1) =\t "<<boost::lexical_cast<std::string>(Pl2(-0.9, 8))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.409685903515625"<<std::endl;
	LegendrePolynom Pl3;
	std::cout<<"===Pl3==="<<std::endl;
	std::cout<<"Pl(0,-1) =\t "<<boost::lexical_cast<std::string>(Pl3(-1, 0))<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	std::cout<<"Pl(0,-1) =\t "<<boost::lexical_cast<std::string>(Pl3(-1, 0))<<std::endl;
	std::cout<<"Wolfram alpha:\t 1"<<std::endl;
	AssociatedLegendrePolynom APl;
	std::cout<<"Associated Legendre polynomials:"<<std::endl;
	std::cout<<"APl(3, 3, 0.36) =\t"<<boost::lexical_cast<std::string>(APl(0.36, 3, 3))<<std::endl;
	std::cout<<"Wolfram alpha:\t -12.18062527"<<std::endl;
	std::cout<<"APl(2, 1, -0.36) =\t"<<boost::lexical_cast<std::string>(APl(-0.36, 2, 1))<<std::endl;
	std::cout<<"Wolfram alpha:\t 1.0075884874"<<std::endl;
	std::cout<<"APl(1, 1, -0.36) =\t"<<boost::lexical_cast<std::string>(APl(-0.36, 1, 1))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.93295230317"<<std::endl;
	std::cout<<"APl(10, 1, 0.27) =\t"<<boost::lexical_cast<std::string>(APl(0.27, 10, 1))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.73106380723"<<std::endl;
	std::cout<<"APl(10, 1, 0.27) =\t"<<boost::lexical_cast<std::string>(APl(0.27, 10, 1))<<std::endl;
	std::cout<<"Wolfram alpha:\t -0.73106380723"<<std::endl;
}

void test_legendre_intregral (void)
{
	int Ncalls = 10;
	auto start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_1_0 (20, 13);
	}*/
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 13, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 13, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_high_dx = end-start;
	/*start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 13, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(13,x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_1_0 (20, 13)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<Int_PlPl (20, 13, -1, 0, 1e-5)<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<Int_PlPl (20, 13, -1, 0, 1e-7)<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (20, 13, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.00217109"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_1_0 (20, 18);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 18, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 18, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 18, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(18,x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_1_0 (20, 18)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<boost::lexical_cast<std::string>(Int_PlPl (20, 18, -1, 0, 1e-5))<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<boost::lexical_cast<std::string>(Int_PlPl (20, 18, -1, 0, 1e-7))<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (20, 18, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.0"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_1_0 (18, 18);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (18, 18, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (18, 18, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (18, 18, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(18,x)Pl(18,x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_1_0 (18, 18)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<boost::lexical_cast<std::string>(Int_PlPl (18, 18, -1, 0, 1e-5))<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<boost::lexical_cast<std::string>(Int_PlPl (18, 18, -1, 0, 1e-7))<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (18, 18, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.027027"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_0_1 (20, 15);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 15, 0, 1, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 15, 0, 1, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl (20, 15, 0, 1, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(15,x)dx [0.0 ; 1.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_0_1 (20, 15)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<boost::lexical_cast<std::string>(Int_PlPl (20, 15, 0, 1, 1e-5))<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<boost::lexical_cast<std::string>(Int_PlPl (20, 15, 0, 1, 1e-7))<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl (20, 15, 0, 1, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.00307571"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf_1_0 (20, 15);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, -1, 0, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, -1, 0, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, -1, 0, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(15,x)(1-x)dx [-1.0 ; 0.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_transf_1_0 (20, 15)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<boost::lexical_cast<std::string>(Int_PlPl_transf (20, 15, -1, 0, 1e-5))<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<boost::lexical_cast<std::string>(Int_PlPl_transf (20, 15, -1, 0, 1e-7))<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl_transf (20, 15, -1, 0, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: -0.00307571"<<std::endl;
	std::cout<<"========================"<<std::endl;

	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf_0_1 (20, 15);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_my = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, 0, 1, 1e-5);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_low_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, 0, 1, 1e-7);
	}
	end = std::chrono::high_resolution_clock::now();
	diff_high_dx = end-start;
	start = std::chrono::high_resolution_clock::now();
	/*for (int i= 0; i< Ncalls; ++i) {
		Int_PlPl_transf (20, 15, 0, 1, 1e-8);
	}*/
	end = std::chrono::high_resolution_clock::now();
	diff_highest_dx = end-start;

	std::cout<<"Pl(20,x)Pl(15,x)dx [0.0 ; 1.0]:"<<std::endl;
	std::cout<<"Method\tValue\ttime[s] for "<<Ncalls<<" calls"<<std::endl;
	//std::cout<<"Sum:\t"<<Int_PlPl_transf_0_1 (20, 15)<<"\t"<<diff_my.count()<<std::endl;
	std::cout<<"dx=1e-5:\t"<<boost::lexical_cast<std::string>(Int_PlPl_transf (20, 15, 0, 1, 1e-5))<<"\t"<<diff_low_dx.count()<<std::endl;
	std::cout<<"dx=1e-7:\t"<<boost::lexical_cast<std::string>(Int_PlPl_transf (20, 15, 0, 1, 1e-7))<<"\t"<<diff_high_dx.count()<<std::endl;
	//std::cout<<"dx=1e-8:\t"<<Int_PlPl_transf (20, 15, 0, 1, 1e-8)<<"\t"<<diff_highest_dx.count()<<std::endl;
	std::cout<<"Wolfram Alpha: 0.00307571"<<std::endl;
	std::cout<<"========================"<<std::endl;
}

void test_colored_interval (void)
{
	ColoredInterval ran1 (-1,0, 1);
	ColoredInterval ran2 (0,1, 1);
	ColoredInterval ran3 (1,4, 2);
	ColoredInterval ran4 (2,3, 0.5);
	ColoredInterval ran5 (2.8,3, 0.1);
	ColoredInterval ran6 (5,6, 0.1);
	ColoredInterval ran7 (-2,10, 5);
	ColoredRange range1, range2;
	std::cout<<"Sum of following intervals [from; color; to]:"<<std::endl;
	std::cout<<"[-1; 1; 0] + [0; 1; 1] + [2; 0.5; 3] + [-2; 5; 10]"<<std::endl;
	std::cout<<"Result: "<<std::endl;
	range1 = ran1 + ran2 + ran4 + ran7;
	range1.Print(std::cout);
	std::cout<<"========================"<<std::endl;
	std::cout<<"Sum of following intervals [from; color; to]:"<<std::endl;
	std::cout<<"[-1; 1; 0] + [0; 1; 1] + [2; 0.5; 3] + [-2; 5; 10] + "<<std::endl;
	std::cout<<"+ [1; 2; 4] + [2.8; 0.1; 3] + [5; 0.1; 6]"<<std::endl;
	std::cout<<"Result: "<<std::endl;
	range2 = range1;
	range2 = range2 + ran3 + ran5 + ran6;
	range2.Print(std::cout);
	std::cout<<"========================"<<std::endl;
	std::cout<<"Scan over following intervals [from; step; to]:"<<std::endl;
	std::cout<<"[0; 0.5; 1] + [1; 1; 3] + [6; 0.5; 7]"<<std::endl;
	ColoredInterval ran8 (0,1, 0.5);
	ColoredInterval ran9 (1,3, 1);
	ColoredInterval ran10 (6,7, 0.5);
	range1 = ran8 + ran9 +ran10;
	for (int i=0, end_=range1.NumOfIndices(); i!=end_; ++i) {
		std::cout<<boost::lexical_cast<std::string>(range1.Value(i))<<"; ";
	}
	std::cout<<std::endl;
	std::cout<<"========================"<<std::endl;
}

void test_factor_helper (void)
{
	LargeFactorHelper helper;
	std::cout<<"396*500*22/2^5/11"<<std::endl;
	helper.MultiplyBy(396);
	helper.MultiplyBy(500);
	helper.MultiplyBy(22);
	helper.DivideBy(2, 5);
	helper.DivideBy(11);
	helper.Print();
	std::cout<<"Direct = ";
	long double val = 396.0*500*22/pow(2,5)/11;
	std::cout<<boost::lexical_cast<std::string>(val)<<std::endl;
	std::cout<<"Helper = "<< boost::lexical_cast<std::string>(helper.Output())<<std::endl;
	std::cout<<"Wolfram = 12375"<<std::endl;
	helper.Clear();
	//=================================================
	std::cout<<"C(9,20)/80"<<std::endl;
	helper.MultiplyByFactorial(20);
	helper.DivideByFactorial(9);
	helper.DivideByFactorial(11);
	helper.DivideBy(80);
	helper.Print();
	std::cout<<"Direct = ";
	val = 1;
	for (int i=12;i<=20;++i)
		val*=i;
	for (int i=1;i<=9;++i)
		val/=i;
	val/=80;
	std::cout<<boost::lexical_cast<std::string>(val)<<std::endl;
	std::cout<<"Helper = "<< boost::lexical_cast<std::string>(helper.Output())<<std::endl;
	std::cout<<"Wolfram = 2099.5"<<std::endl;
	helper.Clear();
	//=================================================
	std::cout<<"C(9,20)*C(20, 60)/80/3^6/10!"<<std::endl;
	helper.MultiplyByFactorial(20);
	helper.DivideByFactorial(9);
	helper.DivideByFactorial(11);
	helper.MultiplyByFactorial(60);
	helper.DivideByFactorial(20);
	helper.DivideByFactorial(40);
	helper.DivideBy(80);
	helper.DivideBy(3,6);
	helper.DivideByFactorial(10);
	helper.Print();
	std::cout<<"Direct = ";
	val = 1;
	for (int i=12;i<=20;++i)
		val*=i;
	for (int i=2;i<=9;++i)
		val/=i;
	for (int i=41;i<=60;++i)
		val*=i;
	for (int i=2;i<=20;++i)
		val/=i;
	val/=80;
	val/=pow(3,6);
	for (int i=2;i<=10;++i)
		val/=i;
	std::cout<<val<<std::endl;
	std::cout<<"Helper = "<< boost::lexical_cast<std::string>(helper.Output())<<std::endl;
	std::cout<<"Wolfram = 3.32682902e+9"<<std::endl;
	helper.Clear();
	//=================================================
	LargeFactorHelper helper1;
	std::cout<<"30*27*17/125/121/8 + 3^4*13*5/17/11/1024:"<<std::endl;
	helper.MultiplyBy(30);
	helper.MultiplyBy(27);
	helper.MultiplyBy(17);
	helper.DivideBy(125);
	helper.DivideBy(121);
	helper.DivideBy(8);
	helper.Print();
	helper1.MultiplyBy(3,4);
	helper1.MultiplyBy(13);
	helper1.MultiplyBy(5);
	helper1.DivideBy(17);
	helper1.DivideBy(11);
	helper1.DivideBy(1024);
	helper1.Print();
	std::cout<<"Direct = ";
	val = 30.0*27*17/125/121/8 + 9.0*9*13*5/17/11/1024;
	std::cout<<boost::lexical_cast<std::string>(val)<<std::endl;
	helper+=helper1;
	std::cout<<"Helper = "<< boost::lexical_cast<std::string>(helper.Output())<<std::endl;
	helper.Clear();
	helper1.Clear();
	std::cout<<"27/81 - 6/18:"<<std::endl;
	helper.MultiplyBy(27);
	helper.DivideBy(81);
	helper.Print();
	helper1.MultiplyBy(6);
	helper1.DivideBy(18);
	helper1.Print();
	std::cout<<"Direct = ";
	val = 27.0/81 - 6.0/18;
	std::cout<<boost::lexical_cast<std::string>(val)<<std::endl;
	helper-=helper1;
	std::cout<<"Helper = "<< boost::lexical_cast<std::string>(helper.Output())<<std::endl;
	helper.Clear();
	helper1.Clear();
}

void test_diff_tot_cross (void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::string fname_diff = prefix + "diff_cross_elastic_10eV.txt";
	std::string fname_tot_MERT5 = prefix + "total_elastic_from_diff_MERT5.txt";
	std::string fname_tot_EXP = prefix + "total_elastic_from_diff_EXP.txt";
	std::string fname_tot = prefix + "total_elastic_from_diff.txt";
	std::string fname_angle_profiles = prefix + "diff_cross_elastic_profiles.txt";
	std::ofstream str;
	int err;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	const Particle* electron = gParticleTable.GetParticle(ELECTRON_NAME);
	if (!electron->isValid() || !argon->isValid()) {
		std::cerr<<"test_diff_tot_cross:: invalid particle, quitting test"<<std::endl;
		return;
	}
#ifndef _NO_CERN_ROOT
	{
		Double_t *ths, *XSs;
		ths = new Double_t[400];
		XSs = new Double_t[400];
		double Int = 0;
		double Energy = 10.0;//En_3o2_ - 0.5*Width_3o2_;
		for (int i = 0; i<400; ++i) {
			ths[i] = i*M_PI / 399.0;
			XSs[i] = argon->ArAllData_.argon_cross_elastic_diff(Energy, ths[i]);
			if (i != 0)
				Int += 0.5*(XSs[i] + XSs[i - 1])*(ths[i] - ths[i - 1]);
		}
		TCanvas *c1 = new TCanvas((std::string("diff. XS ") + std::to_string(Energy) + "eV").c_str(),
			(std::string("diff. XS ") + std::to_string(Energy) + "eV").c_str(), 900, 700);
		TGraph *gr = new TGraph(400, ths, XSs);
		TH1D * hist = new TH1D((std::string("generated thetas ") + std::to_string(Energy) + "eV").c_str(),
			(std::string("generated thetas ") + std::to_string(Energy) + "eV").c_str(), 300, 0, M_PI);
		TRandom *random_generator_ = new TRandom1(42);
		for (int h = 0; h<100000; ++h)
			hist->Fill(argon->GenerateScatterAngle(electron, Energy, 0, random_generator_->Uniform()));
		double Norm = 0;
		for (int bin = 0, bin_end = hist->GetNbinsX() + 1; bin != bin_end; ++bin)
			Norm += hist->GetBinContent(bin)*hist->GetBinWidth(bin);
		Norm /= Int;
		for (int bin = 0, bin_end = hist->GetNbinsX() + 1; bin != bin_end; ++bin)
			hist->SetBinContent(bin, hist->GetBinContent(bin) / (Norm));

		hist->Draw();
		gr->Draw("L");
	}
#else //_NO_CERN_ROOT
	std::cerr << "test_diff_tot_cross: Warning: ROOT is turned off, can't draw histogram of scattering angle generation" << std::endl;
#endif //_NO_CERN_ROOT
	{
		open_output_file(fname_tot_MERT5, str, std::ios_base::trunc);
		str<<"E[eV]\tXS from diff MERT5 [1e-20m^2]\tXS tot MERT5 PS [1e-20m^2]"<<std::endl;
		{
			EnergyScanner eScan(EnergyScanner::PlotDiffXS);
			while (true) {
				double E = eScan.Next(err);
				if ((0!=err)||(E>2.0))
					break;
				long double integral = 0;
				for (int j=0; j<10001; ++j)
					integral+=(M_PI/10000.0)*argon->ArAllData_.argon_cross_elastic_diff(E, j*M_PI/10000.0, 1)*sin(j*M_PI/10000.0);
				str << boost::lexical_cast<std::string>(E) << "\t"
					<< boost::lexical_cast<std::string>(integral) << "\t"
					<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic(E, 1)) << std::endl;
			}
		}
		str.close();

		open_output_file(fname_tot_EXP, str, std::ios_base::trunc);
		str<<"E[eV]\tXS from diff EXP [1e-20m^2]\tXS tot EXP [1e-20m^2]"<<std::endl;
		{
			EnergyScanner eScan(EnergyScanner::PlotDiffXS);
			while (true) {
				double E = eScan.Next(err);
				if ((0!=err))
					break;
				if (E<0.5)
					continue;
				long double integral = 0;
				for (int j=0;j<10001; ++j)
					integral+=(M_PI/10000.0)*argon->ArAllData_.argon_cross_elastic_diff(E, j*M_PI/10000.0, 2)*sin(j*M_PI/10000.0);
				str << boost::lexical_cast<std::string>(E) << "\t"
					<< boost::lexical_cast<std::string>(integral) << "\t"
					<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic(E, 2)) << std::endl;
			}
		}
		str.close();

		open_output_file(fname_tot, str, std::ios_base::trunc);
		str<<"E[eV]\tXS from diff [1e-20m^2]\tXS tot [1e-20m^2]"<<std::endl;
		{
			EnergyScanner eScan(EnergyScanner::PlotDiffXS);
			while (true) {
				double E = eScan.Next(err);
				if ((0!=err))
					break;
				long double integral = 0;
				for (int j=0;j<10001; ++j)
					integral+=(M_PI/10000.0)*argon->ArAllData_.argon_cross_elastic_diff(E, j*M_PI/10000.0)*sin(j*M_PI/10000.0);
				str << boost::lexical_cast<std::string>(E) << "\t"
					<< boost::lexical_cast<std::string>(integral) << "\t"
					<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic(E)) << std::endl;
			}
		}
		str.close();

		std::string name = prefix + "test_diff_XS_elastic_by_total.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str<<"set logscale x"<<std::endl;
		str<<"set logscale y"<<std::endl;
		str<<"set key top left"<<std::endl;
		str<<"plot \""<<fname_tot_MERT5<<"\" u 1:2 title \"total from MERT5 diff. XS\""<<std::endl;
		str<<"replot \""<<fname_tot_MERT5<<"\" u 1:3 title \"total from MERT5 phases\""<<std::endl;
		str<<"replot \""<<fname_tot_EXP<<"\" u 1:2 lc rgb \"#000000\" title \"total from EXP diff. XS\""<<std::endl;
		str<<"replot \""<<fname_tot_EXP<<"\" u 1:3 title \"total from EXP\""<<std::endl;
		str<<"replot \""<<fname_tot<<"\" u 1:2 w lines title \"total XS from diff. (mixed)\""<<std::endl;
		str<<"replot \""<<fname_tot<<"\" u 1:3 w lines title \"total XS (mixed)\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		open_output_file(fname_angle_profiles, str, std::ios_base::trunc);
		str << "E[eV]\tXS diff 22.5deg [1e-20m^2]\tXS diff 45deg\tXS diff 90deg\tXS diff 112.5deg\tXS diff 135deg" << std::endl;
		EnergyScanner eScan(EnergyScanner::PlotResonances);
		double th0 = 22.5*(M_PI) / 180;
		double th1 = 45*(M_PI) / 180;
		double th2 = 90*(M_PI) / 180;
		double th3 = 112.5*(M_PI) / 180;
		double th4 = 135*(M_PI) / 180;
		while (true) {
			double E = eScan.Next(err);
			if ((0 != err))
				break;
			str << boost::lexical_cast<std::string>(E) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th0)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th1)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th2)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th3)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th4)) <<std::endl;
			}
		str.close();
		std::string name = prefix + "test_diff_XS_elastic_profiles.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str << "set key top right" << std::endl;
		str << "set xlabel \"E [eV]\"" << std::endl;
		str << "plot \"" << fname_angle_profiles << "\" u 1:2 w l title \"diff. XS 22.5 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:3 w l title \"diff. XS 45 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:4 w l lc rgb \"#000000\" title \"diff. XS 90 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:5 w l title \"diff. XS 112.5 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:6 w l title \"diff. XS 135.5 deg.\"" << std::endl;
		str << "pause -1" << std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

}

void test_J_L_probabilities(void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::string fname_diff = prefix + "diff_cross_elastic_10eV_Pjl.txt";
	std::string fname_tot = prefix + "total_elastic_from_diff_Pjl.txt";
	std::string fname_angle_profiles = prefix + "diff_cross_elastic_profiles_Pjl.txt";
	std::ofstream str;
	int err;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	if (!argon->isValid()) {
		std::cerr<<"test_total_cross_all:: invalid particle, quitting test"<<std::endl;
		return;
	}

	{
		double E = gSettings.PhysConsts()->En_3o2-0.05;
		unsigned int N_th = 600;
		int N_lz = 0, N_gz = 0;
		double I_lz = 0;
		double I_gz = 0;
		for (unsigned int i = 0; i<=N_th; ++i) {
			double theta = i*M_PI / N_th;
			double k = gSettings.PhysConsts()->a_h_bar_2eM_e_SI*sqrt(E);
			for (unsigned int L = 0; L < argon->ArAllData_.ArExper_.max_L(k); ++L) {
				double prob = argon->ArAllData_.argon_scatter_probability_j(E, theta, 2 * L + 1, 2 * L);
				if (prob < 0) {
					++N_lz;
					I_lz += prob;
				} else {
					I_gz += prob;
					++N_gz;
				}
				if (L > 0) {
					prob = argon->ArAllData_.argon_scatter_probability_j(E, theta, 2 * L - 1, 2 * L);
					if (prob < 0) {
						++N_lz;
						I_lz += prob;
					} else {
						I_gz += prob;
						++N_gz;
					}
				}
			}
		}
		std::cout << "Test of probability of scattering through defined J L " << std::endl;
		std::cout << "E = " << boost::lexical_cast<std::string>(E) << "; N_theta = "<< N_th << std::endl;
		std::cout << "Num of negative probabilities: "<< N_lz << std::endl;
		std::cout << "Num of positive probabilities: "<< N_gz << std::endl;
		std::cout << "Sum of negative probabilities: " << I_lz << std::endl;
		std::cout << "Sum of positive probabilities: " << I_gz << std::endl;
	}

	{
		open_output_file(fname_diff, str, std::ios_base::trunc);
		str<<"theta[deg]\tdiff.XS_standard[1e-20m^2]\tdiff.XS_Pjl[1e-20m^2]"<<std::endl;
		for (int i=0; i<600; ++i) {
			double th = i*(M_PI)/599;
			str << boost::lexical_cast<std::string>(th*180/M_PI) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(10.0, th)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(10.0, th, 3)) << std::endl;
		}
		str.close();
		std::string name = prefix + "test_diff_XS_elastic_10eV_Pjl.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str<<"plot \""<<fname_diff<<"\" u 1:2 title \"Diff. XS at 10 eV\""<<std::endl;
		str << "replot \"" << fname_diff << "\" u 1:3 title \"Diff. XS at 10 eV using l-j probabilities\"" << std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}
	
	{
		open_output_file(fname_tot, str, std::ios_base::trunc);
		str << "E[eV]\tXS tot [1e-20m^2]\tXS from diff [1e-20m^2]\tXS from diff Pjl [1e-20m^2]" << std::endl;
		{
			EnergyScanner eScan(EnergyScanner::PlotDiffXS);
			while (true) {
				double E = eScan.Next(err);
				if ((0 != err))
					break;
				long double integral = 0;
				long double integral2 = 0;
				for (int j = 0; j < 10001; ++j) {
					integral += (M_PI / 10000.0)*argon->ArAllData_.argon_cross_elastic_diff(E, j*M_PI / 10000.0)*sin(j*M_PI / 10000.0);
					integral2 += (M_PI / 10000.0)*argon->ArAllData_.argon_cross_elastic_diff(E, j*M_PI / 10000.0, 3)*sin(j*M_PI / 10000.0);
				}
				str << boost::lexical_cast<std::string>(E) << "\t"
					<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic(E)) << "\t"
					<< boost::lexical_cast<std::string>(integral) << "\t"
					<< boost::lexical_cast<std::string>(integral2) << std::endl;
			}
		}
		str.close();

		std::string name = prefix + "test_diff_XS_elastic_by_total_Pjl.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str << "set logscale x" << std::endl;
		str << "set logscale y" << std::endl;
		str << "set key top left" << std::endl;
		str << "plot \"" << fname_tot << "\" u 1:2 w lines title \"total XS\"" << std::endl;
		str << "replot \"" << fname_tot << "\" u 1:3 w lines title \"total XS from diff.\"" << std::endl;
		str << "replot \"" << fname_tot << "\" u 1:4 w lines title \"total XS from diff. by Pjl\"" << std::endl;
		str << "pause -1" << std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		open_output_file(fname_angle_profiles, str, std::ios_base::trunc);
		str << "E[eV]\tXS diff 22.5deg [1e-20m^2]\tXS diff 45deg\tXS diff 90deg\tXS diff 112.5deg\tXS diff 135deg";
		str << "\tXS diff Pjl 22.5deg\tXS diff Pjl 45deg\tXS diff Pjl 90deg\tXS diff Pjl 112.5deg\tXS diff Pjl 135deg" << std::endl;
		EnergyScanner eScan(EnergyScanner::PlotResonances);
		double th0 = 22.5*(M_PI) / 180;
		double th1 = 45 * (M_PI) / 180;
		double th2 = 90 * (M_PI) / 180;
		double th3 = 112.5*(M_PI) / 180;
		double th4 = 135 * (M_PI) / 180;
		while (true) {
			double E = eScan.Next(err);
			if ((0 != err))
				break;
			str << boost::lexical_cast<std::string>(E) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th0)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th1)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th2)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th3)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th4)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th0, 3)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th1, 3)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th2, 3)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th3, 3)) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic_diff(E, th4, 3)) << std::endl;
		}
		str.close();
		std::string name = prefix + "test_diff_XS_elastic_profiles_Pjl.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str << "set key top right" << std::endl;
		str << "set xlabel \"E [eV]\"" << std::endl;
		str << "plot \"" << fname_angle_profiles << "\" u 1:2 w l title \"diff. XS 22.5 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:3 w l title \"diff. XS 45 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:4 w l lc rgb \"#000000\" title \"diff. XS 90 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:5 w l title \"diff. XS 112.5 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:6 w l title \"diff. XS 135.5 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:7 w l title \"diff. XS Pjl 22.5 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:8 w l title \"diff. XS Pjl 45 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:9 w l title \"diff. Pjl XS 90 deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:10 w l title \"diff. XS 112.5 Pjl deg.\"" << std::endl;
		str << "replot \"" << fname_angle_profiles << "\" u 1:11 w l title \"diff. XS 135.5 Pjl deg.\"" << std::endl;
		str << "pause -1" << std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}
}

void test_time_delay()
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	double max_delay = 1.05*6*gSettings.PhysConsts()->h_bar_eVs/std::min(gSettings.PhysConsts()->Width_1o2, gSettings.PhysConsts()->Width_3o2);
	double min_delay = -max_delay;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	const Particle* electron = gParticleTable.GetParticle(ELECTRON_NAME);
	if (!electron->isValid() || !argon->isValid()) {
		std::cerr<<"test_time_delay:: invalid particle, quitting test"<<std::endl;
		return;
	}
#ifndef _NO_CERN_ROOT
	{
		//11.075-11.125
		//11.075-11.30
		double Efr = 11.075, Eto = 11.30;
		std::stringstream En1, En2;
		En1.precision(5);
		En1.precision(5);
		En1<<Efr;
		En2<<Eto;
		std::string title0 = std::string("Time delay ") +En1.str()+"-"+En2.str()+" eV";
		std::string title1 = std::string("Time delay at ") +En1.str()+"-"+En2.str()+" eV";
		std::string title2 = std::string("Tabulated time delay at ") +En1.str()+"-"+En2.str()+" eV";

		TCanvas *c1 = new TCanvas(title0.c_str(), title0.c_str(), 1000, 700);
		c1->SetLogy();
		TLegend* legend = new TLegend(0.50, 0.75, 0.9, 0.9);
		//legend->SetHeader("");
		legend->SetMargin(0.25);
		TH1D * hist1 = new TH1D(title1.c_str(), title1.c_str(), 600, min_delay, max_delay);
		TH1D * hist2 = new TH1D(title2.c_str(), title2.c_str(), 600, min_delay, max_delay);
		hist2->SetLineColor(kRed);
		TRandom *random_generator_ = new TRandom1(42);
		//unsigned long int N0 = 0, Nn0 = 0;
		for (unsigned int h = 0; h < 2000000u; ++h) {
			double Energy = Efr + (Eto-Efr)*random_generator_->Uniform();
			double theta = argon->GenerateUntabScatterAngle(electron, Energy, 0, random_generator_->Uniform());
			double delay = argon->GenerateUntabTimeDelay(electron, Energy, theta, 0, random_generator_->Uniform());
			hist1->Fill(delay);
			Energy = Efr + (Eto-Efr)*random_generator_->Uniform();
			theta = argon->GenerateScatterAngle(electron, Energy, 0, random_generator_->Uniform());
			delay = argon->GenerateTimeDelay(electron, Energy, theta, 0, random_generator_->Uniform());
			hist2->Fill(delay);
		}
		legend->AddEntry(hist1, (std::string("N = 2000000, Mean = ") + boost::lexical_cast<std::string>(hist1->GetMean())).c_str(), "l");
		legend->AddEntry(hist2, (std::string("Tabulated. N = 2000000, Mean = ") + boost::lexical_cast<std::string>(hist2->GetMean())).c_str(), "l");
		hist1->Draw();
		hist2->Draw("same");
		legend->Draw("same");
	}
#else //_NO_CERN_ROOT
	std::cerr << "test_time_delay: Warning: ROOT is turned off, can't draw histogram of time delay over energy region." << std::endl;
#endif //_NO_CERN_ROOT
#ifndef _NO_CERN_ROOT
	{
		std::vector<double> Es = {11.093, gSettings.PhysConsts()->En_3o2 - gSettings.PhysConsts()->Width_3o2, gSettings.PhysConsts()->En_3o2, gSettings.PhysConsts()->En_1o2};
		std::vector<Color_t> colors = {kBlue, kBlack, kRed, kGreen};
		TLegend* legend = new TLegend(0.50, 0.75, 0.9, 0.9);
		//legend->SetHeader("");
		legend->SetMargin(0.25);
		std::string title0 = std::string("Time delay for different E");
		TCanvas *c1 = new TCanvas(title0.c_str(), title0.c_str(), 1000, 700);
		c1->SetLogy();
		for (std::size_t i=0; i!=Es.size(); ++i) {
			std::stringstream En;
			En.precision(6);
			En<<Es[i];
			std::string title1 = std::string("Time delay fixed E") + (i==0 ? std::string() : En.str() + " eV");
			TH1D * hist1 = new TH1D(title1.c_str(), title1.c_str(), 600, min_delay, max_delay);
			hist1->SetLineColor(colors[i]);
			TRandom *random_generator_ = new TRandom1(42);
			//unsigned long int N0 = 0, Nn0 = 0;
			for (int h = 0; h < 2000000; ++h) {
				double Energy = Es[i];
				double theta = argon->GenerateScatterAngle(electron, Energy, 0, random_generator_->Uniform());
				double delay = argon->GenerateTimeDelay(electron, Energy, theta, 0, random_generator_->Uniform());
				hist1->Fill(delay);
			}
			legend->AddEntry(hist1, (std::string("N = 2000000, ")+En.str()+" eV, Mean = " + boost::lexical_cast<std::string>(hist1->GetMean())).c_str(), "l");
			if (i==0)
				hist1->Draw();
			else
				hist1->Draw("same");
		}
		legend->Draw("same");
	}
#else //_NO_CERN_ROOT
	std::cerr << "test_time_delay: Warning: ROOT is turned off, can't draw histogram of time delay at fixed energies." << std::endl;
#endif //_NO_CERN_ROOT
	{
		std::vector<double> Es = { 11.093, gSettings.PhysConsts()->En_3o2 - gSettings.PhysConsts()->Width_3o2, gSettings.PhysConsts()->En_3o2, gSettings.PhysConsts()->En_1o2 };
		for (std::size_t i = 0; i != Es.size(); ++i) {
			std::stringstream En;
			En.precision(6);
			En << Es[i];
			std::string fname = prefix + "spin_flip-nonflip_prob_and_dt_" + En.str() + "eV.txt";
			std::ofstream str;
			open_output_file(fname, str, std::ios_base::trunc);
			str << "//E[eV] = " << Es[i] << std::endl;
			str << "//theta\tP_sflip\tT_sflip[s]\tT_snonflip[s]" << std::endl;
			for (unsigned int th = 0, th_end_ = 1301; th < th_end_; ++th) {
				double theta = M_PI * th / (th_end_ - 1);
				str << boost::lexical_cast<std::string>(theta) << "\t"
					<< boost::lexical_cast<std::string>(1 - std::max((*argon->time_delay_spin_nonflip_prob_table_)(theta, Es[i]), 0.0))
					<< "\t" << boost::lexical_cast<std::string>((*argon->time_delay_spin_flip_table_)(theta, Es[i]))
					<< "\t" << boost::lexical_cast<std::string>((*argon->time_delay_spin_nonflip_table_)(theta, Es[i])) << std::endl;
			}
			str.close();
			std::string name = prefix + "spin_flip-nonflip_prob_and_dt_" + En.str() + "eV.sc";
			open_output_file(name, str, std::ios_base::trunc);
			str << "set ytics nomirror" << std::endl;
			str << "set y2tics" << std::endl;
			str << "set xlabel \"theta [rad.]\"" << std::endl;
			str << "set y2label \"Probability\"" << std::endl;
			str << "set ylabel \"Time [s]\"" << std::endl;
			str << "plot \"" << fname << "\" u 1:3 axis x1y1 title \"T spin flip " << En.str() << " eV\"" << std::endl;
			str << "replot \"" << fname << "\" u 1:4 axis x1y1 title \"T spin nonflip " << En.str() << " eV\"" << std::endl;
			str << "replot \"" << fname << "\" u 1:2 axis x1y2 title \"P spin flip " << En.str() << " eV\"" << std::endl;
			str << "replot \"" << fname << "\" u 1:($2*$3+(1-$2)*$4) axis x1y1 w lines lc rgb \"#000000\" title \"<T> " << En.str() << " eV\"" << std::endl;
			//str<<"replot \""<<fname<<"\" u 1:(1-$2) title \"P spin nonflip "<<En.str()<<" eV\""<<std::endl;
			str << "pause -1" << std::endl;
			str.close();
			INVOKE_GNUPLOT(name);
		}
	}

}

void test_backward_scatter_prob (void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::string fname = prefix + "backward_scattering_prob.txt";
	std::ofstream str;
	int err;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	if (!argon->isValid()) {
		std::cerr<<"test_backward_scatter_prob:: invalid particle, quitting test"<<std::endl;
		return;
	}

	{
		open_output_file(fname, str, std::ios_base::trunc);
		str<<"E[eV]\tW_backward"<<std::endl;
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str << boost::lexical_cast<std::string>(E) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_back_scatter_prob(E)) << std::endl;
		}
		str.close();
	}
	std::string name = prefix + "test_backward_scatter.sc";
	open_output_file(name, str, std::ios_base::trunc);
	//str<<"set logscale x"<<std::endl;
	str<<"plot \""<<fname<<"\" u 1:2 title \"W_backward\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_TM_forward (void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::string fname = prefix + "TM_forward.txt";
	std::ofstream str;
	int err;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	if (!argon->isValid()) {
		std::cerr<<"test_TM_forward:: invalid particle, quitting test"<<std::endl;
		return;
	}

	{
		open_output_file(fname, str, std::ios_base::trunc);
		str<<"E[eV]\tTM_forward"<<std::endl;
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str << boost::lexical_cast<std::string>(E) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_TM_forward(E)) << std::endl;
		}
		str.close();
	}
	std::string name = prefix + "test_TM_forward.sc";
	open_output_file(name, str, std::ios_base::trunc);
	//str<<"set logscale x"<<std::endl;
	str<<"plot \""<<fname<<"\" u 1:2 title \"TM_forward\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_TM_backward (void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::string fname = prefix + "TM_backward.txt";
	std::ofstream str;
	int err;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	if (!argon->isValid()) {
		std::cerr<<"test_TM_backward:: invalid particle, quitting test"<<std::endl;
		return;
	}

	{
		open_output_file(fname, str, std::ios_base::trunc);
		str<<"E[eV]\tTM_backward"<<std::endl;
		EnergyScanner eScan(EnergyScanner::PlotDiffXS);
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str << boost::lexical_cast<std::string>(E) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_TM_backward(E)) << std::endl;
		}
		str.close();
	}
	std::string name = prefix + "test_TM_backward.sc";
	open_output_file(name, str, std::ios_base::trunc);
	//str<<"set logscale x"<<std::endl;
	str<<"plot \""<<fname<<"\" u 1:2 title \"TM_backward\""<<std::endl;
	str<<"pause -1"<<std::endl;
	str.close();
	INVOKE_GNUPLOT(name);
}

void test_total_cross_all (void)
{
	std::string prefix = gSettings.ProgConsts()->test_folder;
	std::ofstream str, str1;

	std::string fname_XS = prefix + "total_elastic_and_resonances_XS.txt";
	std::string fname_XS_res = prefix + "total_resonances_XS.txt";
	std::string fname_XS_ext = prefix + "excitation_XS.txt";
	std::string fname_XS_ext2 = prefix + "excitation_XS2.txt";
	int err;
	const ArgonParticle* argon = (ArgonParticle*) gParticleTable.GetParticle(ARGON_NAME);
	const Particle* electron = gParticleTable.GetParticle(ELECTRON_NAME);
	if (!electron->isValid() || !argon->isValid()) {
		std::cerr<<"test_total_cross_all:: invalid particle, quitting test"<<std::endl;
		return;
	}

	{
		EnergyScanner EnRange(EnergyScanner::PlotAllXS);
		open_output_file(fname_XS, str, std::ios_base::trunc);
		str<<std::scientific;
		str<<"E[eV]\tXS elastic+Feshbach resonances total [1e-20 m^2]\tXS elastic total from table"<<std::endl;
		open_output_file(fname_XS_ext, str1, std::ios_base::trunc);
		str1<<std::scientific;
		str1<<"E[eV]\tXS S ext.[1e-20 m^2]\tXS P ext.\tXS sum ext.\tXS ion."<<std::endl;
		while (true) {
			double E = EnRange.Next(err);
			if (0!=err)
				break;

			str<<E<<"\t"<<argon->ArAllData_.argon_cross_elastic(E) << "\t" << argon->GetCrossSection(electron, E, 0) <<std::endl;
			double XS_S=0, XS_P=0, XS_EXT=0, XS_ION=0;
			for (int j =0, end_ = argon->ArAllData_.ArExper_.ionizations.size(); j!=end_; ++j) {
				XS_ION+= argon->ArAllData_.ArExper_.ionizations[j](E);
			}
			for (int j =0, end_ = argon->ArAllData_.ArExper_.excitations.size(); j!=end_; ++j) {
				std::string name = argon->ArAllData_.ArExper_.excitations[j].get_name();
				if (name.find("S")!=std::string::npos)
					XS_S+= argon->ArAllData_.ArExper_.excitations[j](E);
				if (name.find("P")!=std::string::npos)
					XS_P+= argon->ArAllData_.ArExper_.excitations[j](E);
				XS_EXT+= argon->ArAllData_.ArExper_.excitations[j](E);
			}
			XS_S = std::max(1e-4, XS_S); //for logscale
			XS_P = std::max(1e-4, XS_P);
			XS_EXT = std::max(1e-4, XS_EXT);
			XS_ION = std::max(1e-4, XS_ION);
			str1 << boost::lexical_cast<std::string>(E) << "\t"
				 << boost::lexical_cast<std::string>(XS_S) << "\t"
				 << boost::lexical_cast<std::string>(XS_P) << "\t"
				 << boost::lexical_cast<std::string>(XS_EXT) << "\t"
				 << boost::lexical_cast<std::string>(XS_ION) << std::endl;
		}
		str.close();
		str1.close();
		std::string name = prefix + "test_table_XS.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str << "set logscale x" << std::endl;
		str << "plot \"" << fname_XS << "\" u 1:2 w lines title \"XS from function\"" << std::endl;
		str << "replot \"" << fname_XS << "\" u 1:3 title \"XS from table\"" << std::endl;
		str << "pause -1" << std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
		name = prefix + "test_XS_all.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str << "set logscale x" << std::endl;
		str << "set logscale y" << std::endl;
		str << "plot \"" << fname_XS << "\" u 1:2 w lines lc rgb \"#000000\" title \"XS elastic + Fishbach resonance\"" << std::endl;
		str << "replot \"" << fname_XS_ext << "\" u 1:2 w lines lc rgb \"#0000CC\" title \"XS S excitation\"" << std::endl;
		str << "replot \"" << fname_XS_ext << "\" u 1:3 w lines lc rgb \"#10CA73\" title \"XS P excitation\"" << std::endl;
		str << "replot \"" << fname_XS_ext << "\" u 1:4 w lines lc rgb \"#D90F2B\" title \"XS sum excitation\"" << std::endl;
		str << "replot \"" << fname_XS_ext << "\" u 1:5 w lines lc rgb \"#D90FD0\" title \"XS ionization\"" << std::endl;
		str << "pause -1" << std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		EnergyScanner EnRange(EnergyScanner::PlotInelasticXS);
		open_output_file(fname_XS_ext2, str, std::ios_base::trunc);
		str<<std::scientific;
		str<<"E[eV]\tXS S ext.[1e-20 m^2]\tXS P ext.\tXS sum ext.\tXS ion."<<std::endl;
		while (true) {
			double E = EnRange.Next(err);
			if (0!=err)
				break;
			double XS_S=0, XS_P=0, XS_EXT=0, XS_ION=0;
			for (int j =0, end_ = argon->ArAllData_.ArExper_.ionizations.size(); j!=end_; ++j) {
				XS_ION+= argon->ArAllData_.ArExper_.ionizations[j](E);
			}
			for (int j =0, end_ = argon->ArAllData_.ArExper_.excitations.size(); j!=end_; ++j) {
				std::string name = argon->ArAllData_.ArExper_.excitations[j].get_name();
				if (name.find("S")!=std::string::npos)
					XS_S+= argon->ArAllData_.ArExper_.excitations[j](E);
				if (name.find("P")!=std::string::npos)
					XS_P+= argon->ArAllData_.ArExper_.excitations[j](E);
				XS_EXT+= argon->ArAllData_.ArExper_.excitations[j](E);
			}
			XS_S = std::max(1e-4, XS_S); //for logscale
			XS_P = std::max(1e-4, XS_P);
			XS_EXT = std::max(1e-4, XS_EXT);
			XS_ION = std::max(1e-4, XS_ION);
			str << boost::lexical_cast<std::string>(E) << "\t"
				<< boost::lexical_cast<std::string>(XS_S) << "\t"
				<< boost::lexical_cast<std::string>(XS_P) << "\t"
				<< boost::lexical_cast<std::string>(XS_EXT) << "\t"
				<< boost::lexical_cast<std::string>(XS_ION) << std::endl;
		}
		str.close();
		std::string name = prefix + "test_XS_ext.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str << "set logscale x" << std::endl;
		str << "set logscale y" << std::endl;
		str << "plot \"" << fname_XS << "\" u 1:2 w lines lc rgb \"#000000\" title \"XS elastic + Fishbach resonance\"" << std::endl;
		str << "replot \"" << fname_XS_ext2 << "\" u 1:2 w lines lc rgb \"#0000CC\" title \"XS S excitation\"" << std::endl;
		str << "replot \"" << fname_XS_ext2 << "\" u 1:3 w lines lc rgb \"#10CA73\" title \"XS P excitation\"" << std::endl;
		str << "replot \"" << fname_XS_ext2 << "\" u 1:4 w lines lc rgb \"#D90F2B\" title \"XS sum excitation\"" << std::endl;
		str << "replot \"" << fname_XS_ext2 << "\" u 1:5 w lines lc rgb \"#D90FD0\" title \"XS ionization\"" << std::endl;
		str << "pause -1" << std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}

	{
		EnergyScanner eScan(EnergyScanner::PlotResonances);
		open_output_file(fname_XS_res, str, std::ios_base::trunc);
		str<<"E[eV]\tXS resonance[1e-20 m^2]\tXS resonance from table"<<std::endl;
		int err;
		while (true) {
			double E = eScan.Next(err);
			if (0!=err)
				break;
			str << boost::lexical_cast<std::string>(E) << "\t"
				<< boost::lexical_cast<std::string>(argon->ArAllData_.argon_cross_elastic(E)) << "\t"
				<< boost::lexical_cast<std::string>(argon->GetCrossSection(electron, E, 0)) << std::endl;
		}
		str.close();
		std::string name = prefix + "test_resonance_XS.sc";
		open_output_file(name, str, std::ios_base::trunc);
		str<<"plot \"" << fname_XS_res << "\" u 1:3 w l lc rgb \"#000000\" title \"Resonance XS from table\"" << std::endl;
		str<<"replot \""<< fname_XS_res <<"\" u 1:2 w p lc rgb \"#FF0000\" title \"Resonance XS from function\""<<std::endl;
		str<<"pause -1"<<std::endl;
		str.close();
		INVOKE_GNUPLOT(name);
	}
}

void test_all (void)
{
	std::cout<<"Testing polynomial fit:"<<std::endl;
	test_polynomial_fit ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing function table:"<<std::endl;
	test_2_dim_table ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	/*
	std::cout<<"Testing phase shifts fit:"<<std::endl;
	test_phase_shift_fit (ArTables);
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*/
	/*std::cout<<"Testing factor helping class:"<<std::endl;
	test_factor_helper ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*//*
	std::cout<<"Testing legendre polynomials:"<<std::endl;
	test_legendre_polynomial ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
*//*
	std::cout<<"Testing colored intervals:"<<std::endl;
	test_colored_interval ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
*//*
	std::cout<<"Testing integrals of legendre polynomials:"<<std::endl;
	test_legendre_intregral ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*//*
	std::cout<<"Testing differential cross section:"<<std::endl;
	test_diff_tot_cross ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*/
/*	std::cout<<"Testing J-L probabilities:"<<std::endl;
	test_J_L_probabilities ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
*/
/*	std::cout<<"Testing backward scatter probability:"<<std::endl;
	test_backward_scatter_prob ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing forward momentum transfer factor:"<<std::endl;
	test_TM_forward ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;

	std::cout<<"Testing backward momentum transfer factor:"<<std::endl;
	test_TM_backward ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*/
	/*
	std::cout<<"Testing total cross sections:"<<std::endl;
	test_total_cross_all ();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*//*
	std::cout<<"Testing time delay:"<<std::endl;
	test_time_delay();
	std::cout<<"==============================================="<<std::endl<<std::endl<<std::endl;
	*/
	std::cout<<"Testing finished."<<std::endl<<std::endl<<std::endl;
}
