#include "FunctionTable.h"

FunctionTable::FunctionTable(void)
{}

std::pair<long int, long int> FunctionTable::find_Ey_indexes (double Ey)
{
	std::pair<long int, long int> out(-1,-1);
	if (_Eys.empty()) {
		return out;
	} else {
		if (Ey >= _Eys.back()) {
			out.first = _Eys.size()-1;
			out.second = out.first;
			return out;
		}
		if (Ey < _Eys.front()) {
			out.first = 0;
			out.second = 0;
			return out;
		}
		if (Ey==_Eys.front()) {
			out.first = 0;
			out.second = _Eys.size()>1 ? 1 : 0;
			return out;
		}
		std::size_t left = 0;
		std::size_t right = std::max(left, _Eys.size() - 2);
		std::size_t Ey_index = (left+right)/2;
		while (true) {
			if ((Ey > _Eys[Ey_index]) && (Ey <= _Eys[Ey_index + 1]))
				break;
			if (Ey <= _Eys[Ey_index]) {
				right = Ey_index - 1;
				Ey_index = (left + right) / 2;
			} else {
				left = Ey_index + 1;
				Ey_index = (left + right) / 2;
			}
		}
		//Ey_index is such, that Ey>_Eys[Ey_index] and Ey<=_Eys[Ey_index+1]
		out.first = Ey_index;
		out.second = Ey_index + 1;
	}
	return out;
}

std::pair<long int, long int> FunctionTable::find_E_indexes (double E, std::size_t Ey_index)
{
	std::pair<long int, long int> out(-1,-1);
	if (_Es[Ey_index].empty()) {
		return out;
	} else {
		if (E >= _Es[Ey_index].back()) {
			out.first = _Es[Ey_index].size()-1;
			out.second = out.first;
			return out;
		}
		if (E < _Es[Ey_index].front()) {
			out.first = 0;
			out.second = 0;
			return out;
		}
		if (E==_Es[Ey_index].front()) {
			out.first = 0;
			out.second = _Es[Ey_index].size()>1 ? 1 : 0;
			return out;
		}
		std::size_t left = 0;
		std::size_t right = std::max(left, _Es[Ey_index].size() - 2);
		std::size_t E_index = (left+right)/2;
		while (true) {
			if ((E > _Es[Ey_index][E_index]) && (E <= _Es[Ey_index][E_index + 1]))
				break;
			if (E <= _Es[Ey_index][E_index]) {
				right = E_index - 1;
				E_index = (left + right) / 2;
			} else {
				left = E_index + 1;
				E_index = (left + right) / 2;
			}
		}
		//Ey_index is such, that Ey>_Eys[Ey_index] and Ey<=_Eys[Ey_index+1]
		out.first = E_index;
		out.second = E_index + 1;
	}
	return out;
}

std::pair<long int, long int> FunctionTable::find_E_indexes_by_value (double val, std::size_t Ey_index)
{
	std::pair<long int, long int> out(-1,-1);
	if (_Es[Ey_index].empty()) {
		return out;
	} else {
		if (val >= _ys[Ey_index].back()) {
			out.first = _ys[Ey_index].size()-1;
			out.second = out.first;
			return out;
		}
		if (val < _ys[Ey_index].front()) {
			out.first = 0;
			out.second = 0;
			return out;
		}
		if (val==_ys[Ey_index].front()) {
			out.first = 0;
			out.second = _ys[Ey_index].size()>1 ? 1 : 0;
			return out;
		}
		std::size_t left = 0;
		std::size_t right = std::max(left, _ys[Ey_index].size() - 2);
		std::size_t E_index = (left+right)/2;
		while (true) {
			if ((val > _ys[Ey_index][E_index]) && (val <= _ys[Ey_index][E_index + 1]))
				break;
			if (val <= _ys[Ey_index][E_index]) {
				right = E_index - 1;
				E_index = (left + right) / 2;
			} else {
				left = E_index + 1;
				E_index = (left + right) / 2;
			}
		}
		//Ey_index is such, that Ey>_Eys[Ey_index] and Ey<=_Eys[Ey_index+1]
		out.first = E_index;
		out.second = E_index + 1;
	}
	return out;
}

double FunctionTable::operator ()(double E, double Ey)
{
	std::pair<long int, long int> Ey_indexes = find_Ey_indexes(Ey);
	if (-1==Ey_indexes.first)
		return 0;
	if (Ey_indexes.first==Ey_indexes.second) {
		std::pair<long int, long int> E_inds = find_E_indexes(E, Ey_indexes.first);
		if (-1==E_inds.first)
			return 0;
		if (E_inds.first==E_inds.second) {
			return _ys[Ey_indexes.first][E_inds.first];
		} else {
			double y0 = _ys[Ey_indexes.first][E_inds.first], y1 = _ys[Ey_indexes.first][E_inds.second];
			double E0 = _Es[Ey_indexes.first][E_inds.first], E1 = _Es[Ey_indexes.first][E_inds.second];
			return y0 + (y1-y0)*(E-E0)/(E1-E0);
		}
	} else {
		std::pair<long int, long int> E_inds_left = find_E_indexes(E, Ey_indexes.first);
		std::pair<long int, long int> E_inds_right = find_E_indexes(E, Ey_indexes.second);
		if ((-1==E_inds_right.first)&&(-1==E_inds_left.first))
			return 0;
		if (-1==E_inds_right.first) {
			if (E_inds_left.first==E_inds_left.second) {
				return _ys[Ey_indexes.first][E_inds_left.first];
			} else {
				double y0 = _ys[Ey_indexes.first][E_inds_left.first], y1 = _ys[Ey_indexes.first][E_inds_left.second];
				double E0 = _Es[Ey_indexes.first][E_inds_left.first], E1 = _Es[Ey_indexes.first][E_inds_left.second];
				return y0 + (y1-y0)*(E-E0)/(E1-E0);
			}
		}
		if (-1==E_inds_left.first) {
			if (E_inds_right.first==E_inds_right.second) {
				return _ys[Ey_indexes.second][E_inds_right.first];
			} else {
				double y0 = _ys[Ey_indexes.second][E_inds_right.first], y1 = _ys[Ey_indexes.second][E_inds_right.second];
				double E0 = _Es[Ey_indexes.second][E_inds_right.first], E1 = _Es[Ey_indexes.second][E_inds_right.second];
				return y0 + (y1-y0)*(E-E0)/(E1-E0);
			}
		}
		if ((E_inds_left.first==E_inds_left.second)&&(E_inds_right.first==E_inds_right.second)) {
			double y0 = _ys[Ey_indexes.first][E_inds_left.first], y1 = _ys[Ey_indexes.second][E_inds_right.first];
			double E0 = _Eys[Ey_indexes.first], E1 = _Eys[Ey_indexes.second];
			return y0 + (y1-y0)*(Ey-E0)/(E1-E0);
		}
		if (E_inds_left.first==E_inds_left.second) {//three point linear interpolation
			double y0_0 = _ys[Ey_indexes.first][E_inds_left.first], y1_0 = _ys[Ey_indexes.second][E_inds_right.first], y1_1 = _ys[Ey_indexes.second][E_inds_right.second];
			double E0_0 = _Es[Ey_indexes.first][E_inds_left.first], E1_0 = _Es[Ey_indexes.second][E_inds_right.first], E1_1 = _Es[Ey_indexes.second][E_inds_right.second];
			double Ey0 = _Eys[Ey_indexes.first], Ey1 = _Eys[Ey_indexes.second];
			//first interpolate along Ey axis, then E axis
			double E0 = E0_0 + (Ey-Ey0)*(E1_0-E0_0)/(Ey1-Ey0);
			double E1 = E0_0 + (Ey-Ey0)*(E1_1-E0_0)/(Ey1-Ey0);
			double y0 = y0_0 + (Ey-Ey0)*(y1_0-y0_0)/(Ey1-Ey0);
			double y1 = y0_0 + (Ey-Ey0)*(y1_1-y0_0)/(Ey1-Ey0);
			return y0 + (y1-y0)*(E-E0)/(E1-E0);
		}
		if (E_inds_right.first==E_inds_right.second) {//three point linear interpolation
			double y0_0 = _ys[Ey_indexes.first][E_inds_left.first], y1_0 = _ys[Ey_indexes.second][E_inds_right.first], y0_1 = _ys[Ey_indexes.first][E_inds_left.second];
			double E0_0 = _Es[Ey_indexes.first][E_inds_left.first], E1_0 = _Es[Ey_indexes.second][E_inds_right.first], E0_1 = _Es[Ey_indexes.first][E_inds_left.second];
			double Ey0 = _Eys[Ey_indexes.first], Ey1 = _Eys[Ey_indexes.second];
			//first interpolate along Ey axis, then E axis
			double E0 = E0_0 + (Ey-Ey0)*(E1_0-E0_0)/(Ey1-Ey0);
			double E1 = E0_1 + (Ey-Ey0)*(E1_0-E0_1)/(Ey1-Ey0);
			double y0 = y0_0 + (Ey-Ey0)*(y1_0-y0_0)/(Ey1-Ey0);
			double y1 = y0_1 + (Ey-Ey0)*(y1_0-y0_1)/(Ey1-Ey0);
			return y0 + (y1-y0)*(E-E0)/(E1-E0);
		}
		//4 point linear interpolation
		double y0_0 = _ys[Ey_indexes.first][E_inds_left.first], y0_1 = _ys[Ey_indexes.first][E_inds_left.second],
				y1_0 = _ys[Ey_indexes.second][E_inds_right.first], y1_1 = _ys[Ey_indexes.second][E_inds_right.second];
		double E0_0 = _Es[Ey_indexes.first][E_inds_left.first], E0_1 = _Es[Ey_indexes.first][E_inds_left.second],
				E1_0 = _Es[Ey_indexes.second][E_inds_right.first], E1_1 = _Es[Ey_indexes.second][E_inds_right.second];
		double Ey0 = _Eys[Ey_indexes.first], Ey1 = _Eys[Ey_indexes.second];
		//first interpolate along Ey axis, then E axis
		double E0 = E0_0 + (Ey-Ey0)*(E1_0-E0_0)/(Ey1-Ey0);
		double E1 = E0_1 + (Ey-Ey0)*(E1_1-E0_1)/(Ey1-Ey0);
		double y0 = y0_0 + (Ey-Ey0)*(y1_0-y0_0)/(Ey1-Ey0);
		double y1 = y0_1 + (Ey-Ey0)*(y1_1-y0_1)/(Ey1-Ey0);
		return y0 + (y1-y0)*(E-E0)/(E1-E0);
	}
	return 0;
}

double FunctionTable::find_E (double Ey, double val)
{
	std::pair<long int, long int> Ey_indexes = find_Ey_indexes(Ey);
	if (-1==Ey_indexes.first)
		return -1;
	if (Ey_indexes.first==Ey_indexes.second) {
		std::pair<long int, long int> E_inds = find_E_indexes_by_value(val, Ey_indexes.first);
		if (-1==E_inds.first) {
			return -1;
		}
		if (E_inds.first==E_inds.second) {
			return _Es[Ey_indexes.first][E_inds.first];
		} else {
			double y0 = _ys[Ey_indexes.first][E_inds.first], y1 = _ys[Ey_indexes.first][E_inds.second];
			double E0 = _Es[Ey_indexes.first][E_inds.first], E1 = _Es[Ey_indexes.first][E_inds.second];
			return E0 + (E1-E0)*(val-y0)/(y1-y0);
		}
	} else {
		std::pair<long int, long int> E_inds_left = find_E_indexes_by_value(val, Ey_indexes.first);
		std::pair<long int, long int> E_inds_right = find_E_indexes_by_value(val, Ey_indexes.second);
		if ((-1==E_inds_right.first)&&(-1==E_inds_left.first))
			return -1;
		if (-1==E_inds_right.first) {
			if (E_inds_left.first==E_inds_left.second) {
				return _Es[Ey_indexes.first][E_inds_left.first];
			} else {
				double y0 = _ys[Ey_indexes.first][E_inds_left.first], y1 = _ys[Ey_indexes.first][E_inds_left.second];
				double E0 = _Es[Ey_indexes.first][E_inds_left.first], E1 = _Es[Ey_indexes.first][E_inds_left.second];
				return E0 + (E1-E0)*(val-y0)/(y1-y0);
			}
		}
		if (-1==E_inds_left.first) {
			if (E_inds_right.first==E_inds_right.second) {
				return _ys[Ey_indexes.second][E_inds_right.first];
			} else {
				double y0 = _ys[Ey_indexes.second][E_inds_right.first], y1 = _ys[Ey_indexes.second][E_inds_right.second];
				double E0 = _Es[Ey_indexes.second][E_inds_right.first], E1 = _Es[Ey_indexes.second][E_inds_right.second];
				return E0 + (E1-E0)*(val-y0)/(y1-y0);
			}
		}
		if ((E_inds_left.first==E_inds_left.second)&&(E_inds_right.first==E_inds_right.second)) {
			double E0 = _Es[Ey_indexes.first][E_inds_left.first], E1 = _Es[Ey_indexes.second][E_inds_right.first];
			double Ey0 = _Eys[Ey_indexes.first], Ey1 = _Eys[Ey_indexes.second];
			return E0 + (E1-E0)*(Ey-Ey0)/(Ey1-Ey0);
		}
		if (E_inds_left.first==E_inds_left.second) {//three point linear interpolation
			double y0_0 = _ys[Ey_indexes.first][E_inds_left.first], y1_0 = _ys[Ey_indexes.second][E_inds_right.first], y1_1 = _ys[Ey_indexes.second][E_inds_right.second];
			double E0_0 = _Es[Ey_indexes.first][E_inds_left.first], E1_0 = _Es[Ey_indexes.second][E_inds_right.first], E1_1 = _Es[Ey_indexes.second][E_inds_right.second];
			double Ey0 = _Eys[Ey_indexes.first], Ey1 = _Eys[Ey_indexes.second];
			//first interpolate along Ey axis, then E axis
			double E0 = E0_0 + (Ey-Ey0)*(E1_0-E0_0)/(Ey1-Ey0);
			double E1 = E0_0 + (Ey-Ey0)*(E1_1-E0_0)/(Ey1-Ey0);
			double y0 = y0_0 + (Ey-Ey0)*(y1_0-y0_0)/(Ey1-Ey0);
			double y1 = y0_0 + (Ey-Ey0)*(y1_1-y0_0)/(Ey1-Ey0);
			return E0 + (E1-E0)*(val-y0)/(y1-y0);
		}
		if (E_inds_right.first==E_inds_right.second) {//three point linear interpolation
			double y0_0 = _ys[Ey_indexes.first][E_inds_left.first], y1_0 = _ys[Ey_indexes.second][E_inds_right.first], y0_1 = _ys[Ey_indexes.first][E_inds_left.second];
			double E0_0 = _Es[Ey_indexes.first][E_inds_left.first], E1_0 = _Es[Ey_indexes.second][E_inds_right.first], E0_1 = _Es[Ey_indexes.first][E_inds_left.second];
			double Ey0 = _Eys[Ey_indexes.first], Ey1 = _Eys[Ey_indexes.second];
			//first interpolate along Ey axis, then E axis
			double E0 = E0_0 + (Ey-Ey0)*(E1_0-E0_0)/(Ey1-Ey0);
			double E1 = E0_1 + (Ey-Ey0)*(E1_0-E0_1)/(Ey1-Ey0);
			double y0 = y0_0 + (Ey-Ey0)*(y1_0-y0_0)/(Ey1-Ey0);
			double y1 = y0_1 + (Ey-Ey0)*(y1_0-y0_1)/(Ey1-Ey0);
			return E0 + (E1-E0)*(val-y0)/(y1-y0);
		}
		//4 point linear interpolation
		double y0_0 = _ys[Ey_indexes.first][E_inds_left.first], y0_1 = _ys[Ey_indexes.first][E_inds_left.second],
				y1_0 = _ys[Ey_indexes.second][E_inds_right.first], y1_1 = _ys[Ey_indexes.second][E_inds_right.second];
		double E0_0 = _Es[Ey_indexes.first][E_inds_left.first], E0_1 = _Es[Ey_indexes.first][E_inds_left.second],
				E1_0 = _Es[Ey_indexes.second][E_inds_right.first], E1_1 = _Es[Ey_indexes.second][E_inds_right.second];
		double Ey0 = _Eys[Ey_indexes.first], Ey1 = _Eys[Ey_indexes.second];
		//first interpolate along Ey axis, then E axis
		double E0 = E0_0 + (Ey-Ey0)*(E1_0-E0_0)/(Ey1-Ey0);
		double E1 = E0_1 + (Ey-Ey0)*(E1_1-E0_1)/(Ey1-Ey0);
		double y0 = y0_0 + (Ey-Ey0)*(y1_0-y0_0)/(Ey1-Ey0);
		double y1 = y0_1 + (Ey-Ey0)*(y1_1-y0_1)/(Ey1-Ey0);
		return E0 + (E1-E0)*(val-y0)/(y1-y0);
	}
}

void FunctionTable::push (double E, double Ey, double val)
{
	std::size_t Ey_index = 0;
	if (_Eys.empty()) {
		_Eys.push_back(Ey);
		_Es.push_back(std::vector<double>());
		_ys.push_back(std::vector<double>());
	} else {
		if (Ey > _Eys.back()) {
			_Eys.push_back(Ey);
			_Es.push_back(std::vector<double>());
			_ys.push_back(std::vector<double>());
			goto second;
		}
		if (Ey==_Eys.back()) {
			Ey_index = _Eys.size()-1;
			goto second;
		}
		if (Ey < _Eys.front()) {
			Ey_index = 0;
			_Eys.insert(_Eys.begin(), Ey);
			_Es.insert(_Es.begin(), std::vector<double>());
			_ys.insert(_ys.begin(), std::vector<double>());
			goto second;
		}
		if (Ey==_Eys.front()) {
			Ey_index =0;
			goto second;
		}
		std::size_t left = 0;
		std::size_t right = std::max(left, _Eys.size() - 2);
		Ey_index = (left+right)/2;
		while (true) {
			if ((Ey > _Eys[Ey_index]) && (Ey <= _Eys[Ey_index + 1]))
				break;
			if (Ey <= _Eys[Ey_index]) {
				right = Ey_index - 1;
				Ey_index = (left + right) / 2;
			} else {
				left = Ey_index + 1;
				Ey_index = (left + right) / 2;
			}
		}
		//Ey_index is such, that Ey>_Eys[Ey_index] and Ey<=_Eys[Ey_index+1]
		++Ey_index;
		if (Ey==_Eys[Ey_index]) {
			goto second;
		} else {
			_Eys.insert(_Eys.begin()+Ey_index,Ey);
			_Es.insert(_Es.begin()+Ey_index, std::vector<double>());
			_ys.insert(_ys.begin()+Ey_index, std::vector<double>());
			goto second;
		}
	}
	second:;
	std::size_t E_index = 0;
	if (_Es[Ey_index].empty()) {
		_Es[Ey_index].push_back(E);
		_ys[Ey_index].push_back(val);
	} else {
		if (E > _Es[Ey_index].back()) {
			_Es[Ey_index].push_back(E);
			_ys[Ey_index].push_back(val);
			return;
		}
		if (E==_Es[Ey_index].back()) {
			E_index =_Es[Ey_index].size()-1;
			_ys[Ey_index][E_index] = val;
			return;
		}
		if (E < _Es[Ey_index].front()) {
			E_index = 0;
			_Es[Ey_index].insert(_Es[Ey_index].begin(), E);
			_ys[Ey_index].insert(_ys[Ey_index].begin(), val);
			return;
		}
		if (E==_Es[Ey_index].front()) {
			E_index =0;
			_ys[Ey_index][E_index] = val;
			return;
		}
		std::size_t left = 0;
		std::size_t right = std::max(left, _Es[Ey_index].size() - 2);
		E_index = (left+right)/2;
		while (true) {
			if ((E > _Es[Ey_index][E_index]) && (E <= _Es[Ey_index][E_index + 1]))
				break;
			if (E <= _Es[Ey_index][E_index]) {
				right = E_index - 1;
				E_index = (left + right) / 2;
			} else {
				left = E_index + 1;
				E_index = (left + right) / 2;
			}
		}
		//E_index is such, that E>_Es[Ey_index][E_index] and E<=_Es[Ey_index][E_index+1]
		++E_index;
		if (E==_Es[Ey_index][E_index]) {
			_ys[Ey_index][E_index] = val;
		} else {
			_Es[Ey_index].insert(_Es[Ey_index].begin()+E_index, E);
			_ys[Ey_index].insert(_ys[Ey_index].begin()+E_index, val);
		}
	}
}

void FunctionTable::clear (void)
{
	_Es.clear();
	_ys.clear();
	_Eys.clear();
}

void FunctionTable::read (std::ifstream& str)
{
	clear();
	double val;
	size_t size, size_E;
	std::size_t counter, counter_E;
	bool valid = true;
	while (!str.eof()) {
		str.read((char*)&size,sizeof(std::size_t));
		_Eys.resize(size);
		_Es.resize(size);
		_ys.resize(size);
		counter = 0;
		while (!str.eof() && counter != size) {
			str.read((char*)&size_E, sizeof(std::size_t));
			_Es[counter].resize(size_E);
			_ys[counter].resize(size_E);
			counter_E = 0;
			while (!str.eof() && counter_E != size_E) {
				str.read((char*)&val,sizeof(double));
				_Es[counter][counter_E] = val;
				if (str.eof()) {
					valid = false;
					break;
				}
				str.read((char*)&val,sizeof(double));
				_ys[counter][counter_E] = val;
				++counter_E;
			}
			++counter;
		}
		break;
	}
	if ((counter_E!=size_E)||(counter!=size)||!valid) {
		std::cout<<"FunctionTable: Error while reading file"<<std::endl;
		clear();
	}
}

void FunctionTable::write (std::ofstream& str)
{
	std::size_t real_size= _Eys.size();
	str.write((char*)&real_size, sizeof(std::size_t));
	for (std::size_t Ey_ind = 0, Ey_ind_end_= real_size; Ey_ind != Ey_ind_end_; ++Ey_ind) {
		str.write((char*)&_Eys[Ey_ind], sizeof(double));
		std::size_t size= _Es[Ey_ind].size();
		str.write((char*)&size, sizeof(std::size_t));
		for (std::size_t E = 0, E_end_= real_size; E != E_end_; ++E) {
			str.write((char*)&_Es[Ey_ind][E], sizeof(double));
			str.write((char*)&_ys[Ey_ind][E], sizeof(double));
		}
	}
}

bool FunctionTable::is_empty(void)
{
	return _Eys.empty();
}
