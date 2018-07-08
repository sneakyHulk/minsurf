#pragma once

#include <vector>

template<class T>
class SYSTEM_MATRIX_CORE {
	unsigned int _m;
protected:
	SYSTEM_MATRIX_CORE(unsigned int m) : _m(m) {}

	//GETTER-METHODS:
	unsigned int m() { return _m; }

	//SETTER-METHODS:
	virtual void insertSYSTEM_MATRIX(unsigned int i, unsigned int j, T value) = 0;

	//CALCULATION-METHODS:
	virtual void calcualteSYSTEM_MATRIX(std::vector<std::vector<T>>& surface) = 0;
};

template<class T>
class LAPLACIAN_OPERATOR : public SYSTEM_MATRIX_CORE<T> {
	using SYSTEM_MATRIX_CORE<T>::m;
	/*
	x	0	x
	1	2	3
	x	4	x
	*/
protected:
	LAPLACIAN_OPERATOR(unsigned int m) : SYSTEM_MATRIX_CORE<T>(m) {}

	//SETTER-METHODS:
	virtual void insertSYSTEM_MATRIX(unsigned int i, unsigned int j, T value) = 0;

	//CALCULATION-METHODS:
	void calcualteSYSTEM_MATRIX(std::vector<std::vector<T>>& surface) {
#pragma omp parallel for
		for (unsigned int i = 0; i < m() * m(); i++) {
			insertSYSTEM_MATRIX(i, i, -4);
		}
#pragma omp parallel for
		for (unsigned int i = 0; i < m() * m() - 1; i++) {
			if (!((i + 1) % m() == 0)) {
				insertSYSTEM_MATRIX(i, i + 1, 1);
			}
		}
#pragma omp parallel for
		for (unsigned int i = 1; i < m() * m(); i++) {
			if (!(i % m() == 0)) {
				insertSYSTEM_MATRIX(i, i - 1, 1);
			}
		}
#pragma omp parallel for
		for (unsigned int i = 0; i < m() * m() - m(); i++) {
			insertSYSTEM_MATRIX(i, i + m(), 1);
		}
#pragma omp parallel for
		for (unsigned int i = m(); i < m() * m(); i++) {
			insertSYSTEM_MATRIX(i, i - m(), 1);
		}
	}
};

template <class T>
class JACOBIAN_HARD_3x3 : public SYSTEM_MATRIX_CORE<T> {
	using SYSTEM_MATRIX_CORE<T>::m;
	/*
	0	1	2
	3	4	5
	6	7	8
	*/
public:
	JACOBIAN_HARD_3x3(unsigned int m) : SYSTEM_MATRIX_CORE<T>(m) {}

	//SETTER-METHODS:
	virtual void insertSYSTEM_MATRIX(unsigned int i, unsigned int j, T value) = 0;

	//CALCULATION-METHODS:
	virtual T calculateDERIVATIVE(unsigned int x, unsigned int derivative_to) = 0;
	void calcualteSYSTEM_MATRIX(std::vector<std::vector<T>>& surface) {
		insertSYSTEM_MATRIX(0, 0, calculateDERIVATIVE(0, 4));
		insertSYSTEM_MATRIX(0, 1, calculateDERIVATIVE(0, 5));
		insertSYSTEM_MATRIX(0, m(), calculateDERIVATIVE(0, 7));
		insertSYSTEM_MATRIX(0, m() + 1, calculateDERIVATIVE(0, 8));

#pragma omp parallel for
		for (int i = 1; i < m() - 1; i++) {
			insertSYSTEM_MATRIX(i, i - 1, calculateDERIVATIVE(i, 3));
			insertSYSTEM_MATRIX(i, i, calculateDERIVATIVE(i, 4));
			insertSYSTEM_MATRIX(i, i + 1, calculateDERIVATIVE(i, 5));
			insertSYSTEM_MATRIX(i, m() + i - 1, calculateDERIVATIVE(i, 6));
			insertSYSTEM_MATRIX(i, m() + i, calculateDERIVATIVE(i, 7));
			insertSYSTEM_MATRIX(i, m() + i + 1, calculateDERIVATIVE(i, 8));
		}

		insertSYSTEM_MATRIX(m() - 1, m() - 2, calculateDERIVATIVE(m() - 1, 3));
		insertSYSTEM_MATRIX(m() - 1, m() - 1, calculateDERIVATIVE(m() - 1, 4));
		insertSYSTEM_MATRIX(m() - 1, 2 * m() - 2, calculateDERIVATIVE(m() - 1, 6));
		insertSYSTEM_MATRIX(m() - 1, 2 * m() - 1, calculateDERIVATIVE(m() - 1, 7));

#pragma omp parallel for
		for (int j = 1; j < m() - 1; j++) {
			insertSYSTEM_MATRIX(m() * j, (j - 1) * m(), calculateDERIVATIVE(m() * j, 1));
			insertSYSTEM_MATRIX(m() * j, (j - 1) * m() + 1, calculateDERIVATIVE(m() * j, 2));
			insertSYSTEM_MATRIX(m() * j, j * m(), calculateDERIVATIVE(m() * j, 4));
			insertSYSTEM_MATRIX(m() * j, j * m() + 1, calculateDERIVATIVE(m() * j, 5));
			insertSYSTEM_MATRIX(m() * j, (j + 1) * m(), calculateDERIVATIVE(m() * j, 7));
			insertSYSTEM_MATRIX(m() * j, (j + 1) * m() + 1, calculateDERIVATIVE(m() * j, 8));

			for (unsigned int i = 1; i < m() - 1; i++) {
				insertSYSTEM_MATRIX(m() * j + i, (j - 1) * m() + (i - 1), calculateDERIVATIVE(m() * j + i, 0));
				insertSYSTEM_MATRIX(m() * j + i, (j - 1) * m() + i, calculateDERIVATIVE(m() * j + i, 1));
				insertSYSTEM_MATRIX(m() * j + i, (j - 1) * m() + (i + 1), calculateDERIVATIVE(m() * j + i, 2));
				insertSYSTEM_MATRIX(m() * j + i, j * m() + (i - 1), calculateDERIVATIVE(m() * j + i, 3));
				insertSYSTEM_MATRIX(m() * j + i, j * m() + i, calculateDERIVATIVE(m() * j + i, 4));
				insertSYSTEM_MATRIX(m() * j + i, j * m() + (i + 1), calculateDERIVATIVE(m() * j + i, 5));
				insertSYSTEM_MATRIX(m() * j + i, (j + 1) * m() + (i - 1), calculateDERIVATIVE(m() * j + i, 6));
				insertSYSTEM_MATRIX(m() * j + i, (j + 1) * m() + i, calculateDERIVATIVE(m() * j + i, 7));
				insertSYSTEM_MATRIX(m() * j + i, (j + 1) * m() + (i + 1), calculateDERIVATIVE(m() * j + i, 8));
			}

			insertSYSTEM_MATRIX(m() * j + m() - 1, (j - 1) * m() + (m() - 2), calculateDERIVATIVE(m() * j + m() - 1, 0));
			insertSYSTEM_MATRIX(m() * j + m() - 1, (j - 1) * m() + (m() - 1), calculateDERIVATIVE(m() * j + m() - 1, 1));
			insertSYSTEM_MATRIX(m() * j + m() - 1, j * m() + (m() - 2), calculateDERIVATIVE(m() * j + m() - 1, 3));
			insertSYSTEM_MATRIX(m() * j + m() - 1, j * m() + (m() - 1), calculateDERIVATIVE(m() * j + m() - 1, 4));
			insertSYSTEM_MATRIX(m() * j + m() - 1, (j + 1) * m() + (m() - 2), calculateDERIVATIVE(m() * j + m() - 1, 6));
			insertSYSTEM_MATRIX(m() * j + m() - 1, (j + 1) * m() + (m() - 1), calculateDERIVATIVE(m() * j + m() - 1, 7));
		}

		insertSYSTEM_MATRIX(m() * (m() - 1), m() * (m() - 2), calculateDERIVATIVE(m() * (m() - 1), 1));
		insertSYSTEM_MATRIX(m() * (m() - 1), m() * (m() - 2) + 1, calculateDERIVATIVE(m() * (m() - 1), 2));
		insertSYSTEM_MATRIX(m() * (m() - 1), m() * (m() - 1), calculateDERIVATIVE(m() * (m() - 1), 4));
		insertSYSTEM_MATRIX(m() * (m() - 1), m() * (m() - 1) + 1, calculateDERIVATIVE(m() * (m() - 1), 5));

#pragma omp parallel for
		for (int i = 1; i < m() - 1; i++) {
			insertSYSTEM_MATRIX(m() * (m() - 1) + i, m() * (m() - 2) + (i - 1), calculateDERIVATIVE(m() * (m() - 1) + i, 0));
			insertSYSTEM_MATRIX(m() * (m() - 1) + i, m() * (m() - 2) + i, calculateDERIVATIVE(m() * (m() - 1) + i, 1));
			insertSYSTEM_MATRIX(m() * (m() - 1) + i, m() * (m() - 2) + (i + 1), calculateDERIVATIVE(m() * (m() - 1) + i, 2));
			insertSYSTEM_MATRIX(m() * (m() - 1) + i, m() * (m() - 1) + (i - 1), calculateDERIVATIVE(m() * (m() - 1) + i, 3));
			insertSYSTEM_MATRIX(m() * (m() - 1) + i, m() * (m() - 1) + i, calculateDERIVATIVE(m() * (m() - 1) + i, 4));
			insertSYSTEM_MATRIX(m() * (m() - 1) + i, m() * (m() - 1) + (i + 1), calculateDERIVATIVE(m() * (m() - 1) + i, 5));
		}

		insertSYSTEM_MATRIX(m() * m() - 1, m() * (m() - 1) - 2, calculateDERIVATIVE(m() * m() - 1, 0));
		insertSYSTEM_MATRIX(m() * m() - 1, m() * (m() - 1) - 1, calculateDERIVATIVE(m() * m() - 1, 1));
		insertSYSTEM_MATRIX(m() * m() - 1, m() * m() - 2, calculateDERIVATIVE(m() * m() - 1, 3));
		insertSYSTEM_MATRIX(m() * m() - 1, m() * m() - 1, calculateDERIVATIVE(m() * m() - 1, 4));
	}
};