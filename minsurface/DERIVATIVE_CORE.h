#pragma once

#include <vector>
#include <string>

template <class T>
class DERIVATIVE_CORE {
protected:
	//CALCULATION-METHODS:
	virtual T calculateDERIVATIVE(unsigned int x, unsigned int derivative_to) = 0;
};

template <class T>
class DERIVATIVE_DIFFERENCEQUOTIENT : public DERIVATIVE_CORE<T> {
	T _epsilon;
protected:
	DERIVATIVE_DIFFERENCEQUOTIENT(T epsilon) : _epsilon(epsilon) {}

	//GETTER-METHODS:
	virtual std::vector<T> getSTENCILS(unsigned int i) = 0;
	virtual T getRESIDUAL(unsigned int i) = 0;

	//CALCULATION-METHODS:
	virtual T calculateRESIDUAL(const std::vector<T> & stencil) = 0;
	virtual T calculateDERIVATIVE(unsigned int x, unsigned int derivative_to) {
		T R = getRESIDUAL(x);
		std::vector<T> stencil = getSTENCILS(x);
		stencil[derivative_to] += _epsilon;
		T Re = calculateRESIDUAL(stencil);
		return (Re - R) / _epsilon;
	}
};

template <class T>
class DERIVATIVE_HAND_3x3 : public DERIVATIVE_CORE<T> {
	T _h;
protected:
	DERIVATIVE_HAND_3x3(T h) : _h(h) {}

	//GETTER-METHODS:
	virtual std::vector<T> getSTENCILS(unsigned int i) = 0;

	//CALCULATION-METHODS:
	virtual T calculateDERIVATIVE(unsigned int x, unsigned int derivative_to) {
		std::vector<T> stencil = getSTENCILS(x);
		T a = stencil[0];
		T b = stencil[1];
		T c = stencil[2];
		T d = stencil[3];
		T e = stencil[4];
		T f = stencil[5];
		T g = stencil[6];
		T i = stencil[7];
		T j = stencil[8];
		if (derivative_to == 0 || derivative_to == 8) {
			return ((b - i) * (f - d)) / (8 * _h * _h * _h * _h);
		}
		else if (derivative_to == 1) {
			return (1 / (8 * _h * _h * _h * _h)) * (-a * d + a * f + 4 * b * d - 8 * b * e + 4 * b * f + c * d - c * f + 2 * d * d - 4 * d * f + d * g - 4 * d * i - d * j + 8 * e * i + 2 * f * f - f * g - 4 * f * i + f * j + 8 * _h * _h);
		}
		else if (derivative_to == 2 || derivative_to == 6) {
			return ((b - i) * (d - f)) / (8 * _h * _h * _h * _h);
		}
		else if (derivative_to == 3) {
			return (1 / (8 * _h * _h * _h * _h)) * (-a * b + a * i + 2 * b * b + b * c + 4 * b * d - 4 * b * f + b * g - 4 * b * i - b * j - c * i - 8 * d * e + 4 * d * i + 8 * e * f - 4 * f * i - g * i + 2 * i * i + i * j + 8 * _h * _h);
		}
		else if (derivative_to == 4) {
			return -(b * b - 2 * b* i + d * d - 2 * d * f + f * f + i * i + 8 * _h * _h) / (2 * _h * _h * _h * _h);
		}
		else if (derivative_to == 5) {
			return (1 / (8 * _h * _h * _h * _h)) * (a * b - a * i + 2 * b * b - b * c - 4 * b * d + 4 * b * f - b * g - 4 * b * i + b * j + c * i + 8 * d * e - 4 * d * i - 8 * e * f + 4 * f * i + g * i + 2 * i * i - i * j + 8 * _h * _h);
		}
		else if (derivative_to == 7) {
			return (1 / (8 * _h * _h * _h * _h)) * (a * d - a * f - 4 * b * d + 8 * b * e - 4 * b * f - c * d + c * f + 2 * d * d - 4 * d * f - d * g + 4 * d * i + d * j - 8 * e * i + 2 * f * f + f * g + 4 * f * i - f * j + 8 * _h * _h);
		}
		else {
			throw std::string("NOT REACHABLE!");
			return -1;
		}
	}
};