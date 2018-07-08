#pragma once

#include <vector>

template <class T>
class RESIDUAL_CORE {
	T _h;
protected:
	RESIDUAL_CORE(T h) : _h(h) {}

	//GETTER-METHODS:
	T h() { return _h; }

	//CALCULATION-METHODS:
	virtual T calculateRESIDUAL(const std::vector<T>& stencil) = 0;
};

template <class T>
class RESIDUAL_3x3 : public RESIDUAL_CORE<T> {
	using RESIDUAL_CORE<T>::h;
	/*
	0	1	2
	3	4	5
	6	7	8
	*/
protected:
	RESIDUAL_3x3<T>(T h) : RESIDUAL_CORE<T>(h) {}

	//CALCULATION-METHODS:
	virtual T calculateRESIDUAL(const std::vector<T> & Z) {
		T tmp = (Z[5] - Z[3]) / (2. * h());
		T I = (1. + tmp * tmp)*((Z[7] - 2. * Z[4] + Z[1]) / (h() * h()));
		T II = -1. / (8. * h() * h() * h() * h())*(Z[5] - Z[3])*(Z[1] - Z[7])*(Z[2] - Z[0] - Z[8] + Z[6]);
		T tmp2 = (Z[1] - Z[7]) / (2. * h());
		T III = (1. + tmp2 * tmp2)*((Z[3] - 2. * Z[4] + Z[5]) / (h() * h()));
		return I + II + III;
	}
};

template <class T>
class RESIDUAL_STAR : public RESIDUAL_CORE<T> {
	using RESIDUAL_CORE<T>::h;
	/*
	x	0	x
	1	2	3
	x	4	x
	*/
protected:
	RESIDUAL_STAR<T>(T h) : RESIDUAL_CORE<T>(h) {}

	//CALCULATION-METHODS:
	virtual T calculateRESIDUAL(const std::vector<T> & Z) {
		T I = (Z[0] + Z[1] - 4 * Z[2] + Z[3] + Z[4]) / (h() * h());
		return I;
	}
};