#pragma once

#include <Eigen/SparseCore>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

#include "DERIVATIVE_CORE.h"
#include "DISCRETIZATION_CORE.h"
#include "SYSTEM_MATRIX_CORE.h"
#include "RESIDUAL_CORE.h"

/*template <class T>
class DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT : public DISCRETIZATION_CORE<T>, public RESIDUAL_3x3<T>, public DERIVATIVE_DIFFERENCEQUOTIENT<T>, public JACOBIAN_HARD_3x3<T> {
	Eigen::SparseMatrix<T, Eigen::RowMajor> JR;
	Eigen::Matrix<T, Eigen::Dynamic, 1> R;
	Eigen::Matrix<T, Eigen::Dynamic, 1> Z;
	Eigen::Matrix<T, Eigen::Dynamic, 1> X;
public:
	DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT(std::vector<std::string>& input) : DISCRETIZATION_CORE<T>(input), RESIDUAL_3x3<T>(DISCRETIZATION_CORE<T>::h()), DERIVATIVE_DIFFERENCEQUOTIENT<T>(10e-5), JACOBIAN_HARD_3x3<T>(DISCRETIZATION_CORE<T>::m()) {
		Z = Eigen::Array<T, Eigen::Dynamic, 1>::Constant(DISCRETIZATION_CORE<T>::n(), DISCRETIZATION_CORE<T>::calculate_arithmetic_mean());

		calculate();

		DISCRETIZATION_CORE<T>::output();
	}
	DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT(DISCRETIZATION_CORE<T>* dis, Eigen::Matrix<T, Eigen::Dynamic, 1> Z) : DISCRETIZATION_CORE<T>(*dis), RESIDUAL_3x3<T>(DISCRETIZATION_CORE<T>::h()), DERIVATIVE_DIFFERENCEQUOTIENT<T>(10e-5), JACOBIAN_HARD_3x3<T>(DISCRETIZATION_CORE<T>::m()), Z(Z) {
		calculate();

		DISCRETIZATION_CORE<T>::output();
	}
	void calculate() {
		JR.resize(DISCRETIZATION_CORE<T>::n(), DISCRETIZATION_CORE<T>::n());
		JR.reserve(Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(DISCRETIZATION_CORE<T>::n(), 9));
		R.resize(DISCRETIZATION_CORE<T>::n());
		X.resize(DISCRETIZATION_CORE<T>::n());

		DISCRETIZATION_CORE<T>::solve_core();
	}
private:
	//GETTER-METHODS:
	virtual std::vector<T> getSTENCILS(unsigned int i) { return DISCRETIZATION_CORE<T>::stencils()[i]; }
	virtual T getRESIDUAL(unsigned int i) { return R(i); }
	virtual T getZ(unsigned int i) { return Z(i); }
	virtual bool stop(unsigned int count) {
		return (count <= 2000 && (X.norm() / DISCRETIZATION_CORE<T>::n()) > 10e-7);
	}

	//SETTER-METHODS:
	virtual void insertRESIDUAL(unsigned int i, T value) { R(i) = value; }
	virtual void insertSYSTEM_MATRIX(unsigned int i, unsigned int j, T value) { JR.coeffRef(i, j) = value; }

	//CALCULATION-METHODS:
	virtual std::vector<T> calculateSTENCIL(unsigned int i, unsigned int j) {
		return { DISCRETIZATION_CORE<T>::surface()[i][j],  DISCRETIZATION_CORE<T>::surface()[i][j + 1], DISCRETIZATION_CORE<T>::surface()[i][j + 2],
			DISCRETIZATION_CORE<T>::surface()[i + 1][j], DISCRETIZATION_CORE<T>::surface()[i + 1][j + 1], DISCRETIZATION_CORE<T>::surface()[i + 1][j + 2],
			DISCRETIZATION_CORE<T>::surface()[i + 2][j], DISCRETIZATION_CORE<T>::surface()[i + 2][j + 1], DISCRETIZATION_CORE<T>::surface()[i + 2][j + 2] };
	}
	virtual T calculateRESIDUAL(const std::vector<T> & stencil) {
		return RESIDUAL_3x3<T>::calculateRESIDUAL(stencil);
	}
	virtual T calculateDERIVATIVE(unsigned int x, unsigned int derivative_to) {
		return DERIVATIVE_DIFFERENCEQUOTIENT<T>::calculateDERIVATIVE(x, derivative_to);
	}
	virtual void solve() {
		JACOBIAN_HARD_3x3<T>::calcualteSYSTEM_MATRIX(DISCRETIZATION_CORE<T>::surface());

		double timestamp = omp_get_wtime();
		Eigen::BiCGSTAB<Eigen::SparseMatrix<T, Eigen::RowMajor>> solver;
		solver.setTolerance(10e-3);
		solver.compute(JR);
		X = solver.solve(R);
		//std::cout << "Anzahl der Iterationen  " << std::endl << "zur Loesung des LGS:    " << solver.iterations() << std::endl;
		//std::cout << "Geschaetzter Fehler:    " << solver.error() << std::endl;
		//std::cout << "Zeit fuer Loesung:      " << omp_get_wtime() - timestamp << std::endl;

		Z = Z - X;
		//std::cout << "Newton Fehlernorm:      " << (X.norm() / DISCRETIZATION_CORE<T>::n()) << std::endl;
	}
};*/

template <class T>
class NEWTON : public DISCRETIZATION_CORE<T>, public RESIDUAL_3x3<T>, public DERIVATIVE_HAND_3x3<T>, public JACOBIAN_HARD_3x3<T> {
protected:
	using DISCRETIZATION_CORE<T>::h; using DISCRETIZATION_CORE<T>::m; using DISCRETIZATION_CORE<T>::n; using DISCRETIZATION_CORE<T>::s;
	typedef Eigen::SparseMatrix<T, Eigen::RowMajor> MT;
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VT;
	MT JR;
	VT R, Z, X, DR;
	T norm_Fz_prev/*||F(z[k-1])||*/, norm_Fz_JRX_prev/*||F(z[k-1]) - F'(z[k-1]) * X[k-1]||*/;
	T ny, ny_prev, ny_max;
	T epsilon, delta;
	int max_newton_iterations, k_tilde, safety_output;

	NEWTON(std::vector<std::string>& input) : DISCRETIZATION_CORE<T>(input), RESIDUAL_3x3<T>(h()), DERIVATIVE_HAND_3x3<T>(h()), JACOBIAN_HARD_3x3<T>(m()), ny(0.1), norm_Fz_JRX_prev(0.0), norm_Fz_prev(0.0), ny_prev(0.0), ny_max(0.9) {
		Z = Eigen::Array<T, Eigen::Dynamic, 1>::Constant(n(), DISCRETIZATION_CORE<T>::calculate_arithmetic_mean());

		start();
	}
	NEWTON(std::vector<std::string>& input, char * vtk_filename) : DISCRETIZATION_CORE<T>(input, vtk_filename), RESIDUAL_3x3<T>(h()), DERIVATIVE_HAND_3x3<T>(h()), JACOBIAN_HARD_3x3<T>(m()), ny(0.1), norm_Fz_JRX_prev(0.0), norm_Fz_prev(0.0), ny_prev(0.0), ny_max(0.9) {
		Z.resize(n());

#pragma omp parallel for
		for (unsigned int i = 0; i < m(); i++) {
			for (unsigned int j = 0; j < m(); j++) {
				Z(m()*i + j) = DISCRETIZATION_CORE<T>::surface()[i][j];
			}
		}

		start();
	}
	NEWTON(DISCRETIZATION_CORE<T>* dis, VT Z) : DISCRETIZATION_CORE<T>(*dis), RESIDUAL_3x3<T>(h()), DERIVATIVE_HAND_3x3<T>(h()), JACOBIAN_HARD_3x3<T>(m()), Z(Z), ny(0.1), norm_Fz_JRX_prev(0.0), norm_Fz_prev(0.0), ny_prev(0.0), ny_max(0.9) {
		start();
	}

	void start() {
		max_newton_iterations = test_on_input_value<int>("MAX NEWTON ITERATIONS", 8, DISCRETIZATION_CORE<T>::input()); // TEST ON MAX NEWTON ITERATIONS
		safety_output = test_on_input_value<int>("NEWTON ITERATIONS PER OUTPUT", 9, DISCRETIZATION_CORE<T>::input()); // TEST ON NEWTON ITERATIONS PER OUTPUT
		epsilon = test_on_input_value<T>("EPSILON", 10, DISCRETIZATION_CORE<T>::input()); // TEST ON EPSILON
		delta = test_on_input_value<T>("DELTA", 11, DISCRETIZATION_CORE<T>::input()); // TEST ON DELTA
		k_tilde = test_on_input_value<int>("K_TILDE", 12, DISCRETIZATION_CORE<T>::input()); // TEST ON K_TILDE

		JR.resize(n(), n());
		JR.reserve(Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(n(), 9));
		R.resize(n());
		X.resize(n());
		DR.resize(n());
	}

	//GETTER-METHODS:
	virtual std::vector<T> getSTENCILS(unsigned int i) { return DISCRETIZATION_CORE<T>::stencils()[i]; }
	virtual T getZ(unsigned int i) { return Z(i); }
	virtual bool stop() { return (DISCRETIZATION_CORE<T>::count()++ <= max_newton_iterations && (DR.norm() / n()) > epsilon); }

	//SETTER-METHODS:
	virtual void insertRESIDUAL(unsigned int i, T value) { R(i) = value; }
	virtual void insertSYSTEM_MATRIX(unsigned int i, unsigned int j, T value) { JR.coeffRef(i, j) = value; }

	//CALCULATION-METHODS:
	virtual std::vector<T> calculateSTENCIL(unsigned int i, unsigned int j) {
		return { DISCRETIZATION_CORE<T>::surface()[i][j],  DISCRETIZATION_CORE<T>::surface()[i][j + 1], DISCRETIZATION_CORE<T>::surface()[i][j + 2],
			DISCRETIZATION_CORE<T>::surface()[i + 1][j], DISCRETIZATION_CORE<T>::surface()[i + 1][j + 1], DISCRETIZATION_CORE<T>::surface()[i + 1][j + 2],
			DISCRETIZATION_CORE<T>::surface()[i + 2][j], DISCRETIZATION_CORE<T>::surface()[i + 2][j + 1], DISCRETIZATION_CORE<T>::surface()[i + 2][j + 2] };
	}
	virtual T calculateRESIDUAL(const std::vector<T> & stencil) { return RESIDUAL_3x3<T>::calculateRESIDUAL(stencil); }
	virtual T calculateDERIVATIVE(unsigned int x, unsigned int derivative_to) { return DERIVATIVE_HAND_3x3<T>::calculateDERIVATIVE(x, derivative_to); }
	virtual void solve() = 0;
	void updateDRhelper() {
#pragma omp parallel for
		for (int i = 0; i < m(); i++) {
			for (unsigned int j = 0; j < m(); j++) {
				unsigned int pos = i * m() + j;
				std::vector<T> stencil = calculateSTENCIL(i, j);
				DISCRETIZATION_CORE<T>::stencils()[pos] = stencil;
				DR(pos) = calculateRESIDUAL(stencil);
			}
		}
	}
	void updateDR() {
#pragma omp parallel for
		for (int i = 1; i < s() - 1; i++) {
			for (unsigned int j = 1; j < s() - 1; j++) {
				unsigned int pos = (i - 1) * m() + (j - 1);
				DISCRETIZATION_CORE<T>::surface()[i][j] = Z(pos) + X(pos);
			}
		}
		updateDRhelper();
	}
	void output() {
		std::cout << "NEWTON ITERATION: " << DISCRETIZATION_CORE<T>::count() << ", " << "SCALED ERROR: '" << (DR.norm() / n()) << "', OUTPUTTING ." << std::flush;
		double timestamp = omp_get_wtime();
		DISCRETIZATION_CORE<T>::output();
		std::cout << ".. TOOK (" << omp_get_wtime() - timestamp << "s)" << std::endl;
	}
	void calculateNY() {
		//ny: Forcing term for linear system solver
		ny = (R.norm() - norm_Fz_JRX_prev) / norm_Fz_prev;
		ny < 0 ? ny = ny * (-1.) : 0;
		T temp = pow(ny_prev, (1. + sqrt(5.)) / 2.);
		(ny < temp) && (temp > 0.1) ? ny = temp : 0;
		ny > ny_max ? ny = ny_max : 0;
	}
	void updateNY_calculateK() {
		norm_Fz_prev = R.norm();
		norm_Fz_JRX_prev = (R + JR * X).norm();
		ny_prev = ny;

		static bool ignore_warnings = false;
		if (!ignore_warnings) {
			static int k = 0;
			if ((DR - R).norm() / DR.norm() < delta) {
				k = k + 1;
			}
			else {
				k = 0;
			}
			if (!(k < k_tilde)) {
				std::string tmp;
				std::cout << "SLOW CONVERGENCE!, CONTINUE? (y/n) ";
				std::cin >> tmp;
				while (tmp != "y" && tmp != "n") {
					std::cout << "SLOW CONVERGENCE!, CONTINUE? (y/n) ";
					std::cin >> tmp;
				}
				if (tmp == "n") {
					std::cout << "LAST "; NEWTON<T>::output();
					throw std::string("STOPPING NEWTON BECAUSE OF SLOW CONVERGENCE!");
				}

				ignore_warnings = true;
			}
		}
	}
	void lin_solve() {
		Eigen::BiCGSTAB<MT> solver;
		solver.setTolerance(ny);
		solver.compute(JR);
		X = solver.solve((-1)*R);
		if (X != X) {
			std::cout << "LAST "; NEWTON<T>::output();
			throw std::string("STOPPING NEWTON BECAUSE OF NON CONVERGENCE! TRY DIFFERENT BOUNDARY VALUES OR METHOD!");
		}
	}
};

template <class T>
class NEWTON_CLASSIC : public NEWTON<T> {
	using DISCRETIZATION_CORE<T>::m; using DISCRETIZATION_CORE<T>::s;
	using NEWTON<T>::JR; using NEWTON<T>::Z; using NEWTON<T>::X; using NEWTON<T>::R; using NEWTON<T>::DR;
	using NEWTON<T>::norm_Fz_prev; using NEWTON<T>::norm_Fz_JRX_prev;
	using NEWTON<T>::ny; using NEWTON<T>::ny_prev; using NEWTON<T>::ny_max;
	using NEWTON<T>::epsilon; using NEWTON<T>::delta;
	using NEWTON<T>::max_newton_iterations;	using NEWTON<T>::k_tilde; using NEWTON<T>::safety_output;
	typedef Eigen::SparseMatrix<T, Eigen::RowMajor> MT;
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VT;
public:
	NEWTON_CLASSIC(std::vector<std::string>& input) : NEWTON<T>(input) { start(); }
	NEWTON_CLASSIC(std::vector<std::string>& input, char * vtk_filename) : NEWTON<T>(input, vtk_filename) { start(); }
	NEWTON_CLASSIC(DISCRETIZATION_CORE<T>* dis, VT Z) : NEWTON<T>(dis, Z) { start(); }
private:
	void start() {
		DISCRETIZATION_CORE<T>::solve_core(); std::cout << "SUCCESS! - LAST "; NEWTON<T>::output(); 
	}
	virtual void solve() {
		static bool first_iter = true;

		if (DISCRETIZATION_CORE<T>::count() % safety_output == 0) {
			NEWTON<T>::output();
		}

		JACOBIAN_HARD_3x3<T>::calcualteSYSTEM_MATRIX(DISCRETIZATION_CORE<T>::surface());

		if (!first_iter) {
			NEWTON<T>::calculateNY();
		}

		NEWTON<T>::lin_solve();
		
		NEWTON<T>::updateDR();
		Z = Z + X;

		NEWTON<T>::updateNY_calculateK();

		first_iter = false;
	}
};

template <class T>
class NEWTON_DAMPED_TRIVIAL : public NEWTON<T> {
	using DISCRETIZATION_CORE<T>::m; using DISCRETIZATION_CORE<T>::s;
	using NEWTON<T>::JR; using NEWTON<T>::Z; using NEWTON<T>::X; using NEWTON<T>::R; using NEWTON<T>::DR;
	using NEWTON<T>::norm_Fz_prev; using NEWTON<T>::norm_Fz_JRX_prev;
	using NEWTON<T>::ny; using NEWTON<T>::ny_prev; using NEWTON<T>::ny_max;
	using NEWTON<T>::epsilon; using NEWTON<T>::delta;
	using NEWTON<T>::max_newton_iterations;	using NEWTON<T>::k_tilde; using NEWTON<T>::safety_output;
	typedef Eigen::SparseMatrix<T, Eigen::RowMajor> MT;
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VT;
	T lambda, min_lambda;
public:
	NEWTON_DAMPED_TRIVIAL(std::vector<std::string>& input) : NEWTON<T>(input), lambda(1) { start(); }
	NEWTON_DAMPED_TRIVIAL(std::vector<std::string>& input, char * vtk_filename) : NEWTON<T>(input, vtk_filename), lambda(1) { start(); }
	NEWTON_DAMPED_TRIVIAL(DISCRETIZATION_CORE<T>* dis, VT Z) : NEWTON<T>(dis, Z), lambda(1) { start(); }
private:
	void start() { 
		min_lambda = test_on_input_value<T>("MIN LAMBDA", 13, DISCRETIZATION_CORE<T>::input()); // TEST ON MIN LAMBDA
		DISCRETIZATION_CORE<T>::solve_core(); std::cout << "SUCCESS! - LAST "; NEWTON<T>::output();
	}
	virtual void solve() {
		static bool first_iter = true;

		if (DISCRETIZATION_CORE<T>::count() % safety_output == 0) {
			NEWTON<T>::output();
		}

		JACOBIAN_HARD_3x3<T>::calcualteSYSTEM_MATRIX(DISCRETIZATION_CORE<T>::surface());

		if (!first_iter) {
			NEWTON<T>::calculateNY();
		}

		NEWTON<T>::lin_solve();

		if (first_iter) {
			DAMPER();
		}
		else {
			if (R.norm() > DR.norm()) {
				Z = Z + lambda * X;
				if (lambda < 1) {
					lambda = 2 * lambda;
				}
			}
			else {
				DAMPER();
			}
		}

		NEWTON<T>::updateNY_calculateK();
		
		first_iter = false;
	}

	void DAMPER() {
		do {
			updateDR();

			lambda = lambda / 2;
		} while (R.norm() < DR.norm() && 2 * lambda > min_lambda);

		lambda = 2 * lambda;
		Z = Z + lambda * X;
	}

	void updateDR() {
#pragma omp parallel for
		for (int i = 1; i < s() - 1; i++) {
			for (unsigned int j = 1; j < s() - 1; j++) {
				unsigned int pos = (i - 1) * m() + (j - 1);
				DISCRETIZATION_CORE<T>::surface()[i][j] = Z(pos) + lambda * X(pos);
			}
		}
		NEWTON<T>::updateDRhelper();
	}
};

template <class T>
class NEWTON_DAMPED_BACKTRACKING : public NEWTON<T> {
	using DISCRETIZATION_CORE<T>::n; using DISCRETIZATION_CORE<T>::m; using DISCRETIZATION_CORE<T>::s;
	using NEWTON<T>::JR; using NEWTON<T>::Z; using NEWTON<T>::X; using NEWTON<T>::R; using NEWTON<T>::DR;
	using NEWTON<T>::norm_Fz_prev; using NEWTON<T>::norm_Fz_JRX_prev;
	using NEWTON<T>::ny; using NEWTON<T>::ny_prev; using NEWTON<T>::ny_max;
	using NEWTON<T>::epsilon; using NEWTON<T>::delta;
	using NEWTON<T>::max_newton_iterations;	using NEWTON<T>::k_tilde; using NEWTON<T>::safety_output;
	typedef Eigen::SparseMatrix<T, Eigen::RowMajor> MT;
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VT;
	VT DReps, DR0;
	T t, thetha, thetha_max, thetha_min;
public:
	NEWTON_DAMPED_BACKTRACKING(std::vector<std::string>& input) : NEWTON<T>(input), t(0.0001), thetha(1.), thetha_min(0.1), thetha_max(0.5) { start(); }
	NEWTON_DAMPED_BACKTRACKING(std::vector<std::string>& input, char * vtk_filename) : NEWTON<T>(input, vtk_filename), t(0.0001), thetha(1.), thetha_min(0.1), thetha_max(0.5) { start(); }
	NEWTON_DAMPED_BACKTRACKING(DISCRETIZATION_CORE<T>* dis, VT Z) : NEWTON<T>(dis, Z), t(0.0001), thetha(1.), thetha_min(0.1), thetha_max(0.5) { start(); }
private:
	void start() {
		DReps.resize(n());
		DR0.resize(n());
		DISCRETIZATION_CORE<T>::solve_core(); std::cout << "SUCCESS! - LAST "; NEWTON<T>::output();
	}
	virtual void solve() {
		static bool first_iter = true;

		if (DISCRETIZATION_CORE<T>::count() % safety_output == 0) {
			NEWTON<T>::output();
		}

		JACOBIAN_HARD_3x3<T>::calcualteSYSTEM_MATRIX(DISCRETIZATION_CORE<T>::surface());

		if (!first_iter) {
			NEWTON<T>::calculateNY();
		}

		NEWTON<T>::lin_solve();

		NEWTON<T>::updateDR();
		T p_0 = R.squaredNorm() / 2.;
		while (DR.norm() > ((1. - t * (1. - ny)) * R.norm())) {
			updateDReps();
			updateDR0();
			T p_prime_0 = (DReps.squaredNorm() - DR0.squaredNorm()) / (2. * 0.000001);
			T p_1 = DR.squaredNorm() / 2.;
			thetha = -p_prime_0 / (2. * (p_1 - p_0 - p_prime_0));
			if (thetha < thetha_min) { thetha = thetha_min; }
			if (thetha > thetha_max) { thetha = thetha_max; }
			X = thetha * X;
			ny = 1. - thetha * (1. - ny);
			NEWTON<T>::updateDR();
		}
		Z = Z + X;

		NEWTON<T>::updateNY_calculateK();

		first_iter = false;
	}

	void updateDReps() {
#pragma omp parallel for
		for (int i = 1; i < s() - 1; i++) {
			for (unsigned int j = 1; j < s() - 1; j++) {
				unsigned int pos = (i - 1) * m() + (j - 1);
				DISCRETIZATION_CORE<T>::surface()[i][j] = Z(pos) + 0.000001 * X(pos);
			}
		}
#pragma omp parallel for
		for (int i = 0; i < m(); i++) {
			for (unsigned int j = 0; j < m(); j++) {
				unsigned int pos = i * m() + j;
				std::vector<T> stencil = NEWTON<T>::calculateSTENCIL(i, j);
				DISCRETIZATION_CORE<T>::stencils()[pos] = stencil;
				DReps(pos) = NEWTON<T>::calculateRESIDUAL(stencil);
			}
		}
	}
	void updateDR0() {
#pragma omp parallel for
		for (int i = 1; i < s() - 1; i++) {
			for (unsigned int j = 1; j < s() - 1; j++) {
				unsigned int pos = (i - 1) * m() + (j - 1);
				DISCRETIZATION_CORE<T>::surface()[i][j] = Z(pos);
			}
		}
#pragma omp parallel for
		for (int i = 0; i < m(); i++) {
			for (unsigned int j = 0; j < m(); j++) {
				unsigned int pos = i * m() + j;
				std::vector<T> stencil = NEWTON<T>::calculateSTENCIL(i, j);
				DISCRETIZATION_CORE<T>::stencils()[pos] = stencil;
				DR0(pos) = NEWTON<T>::calculateRESIDUAL(stencil);
			}
		}
	}
};

template <class T>
class LAPLACE : public DISCRETIZATION_CORE<T>, public RESIDUAL_STAR<T>, public LAPLACIAN_OPERATOR<T> {
	using DISCRETIZATION_CORE<T>::h; using DISCRETIZATION_CORE<T>::m; using DISCRETIZATION_CORE<T>::n; using DISCRETIZATION_CORE<T>::s;
	typedef Eigen::SparseMatrix<T, Eigen::RowMajor> MT;
	typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VT;
	MT A;
	VT R, Z, X;
public:
	LAPLACE(std::vector<std::string>& input) : DISCRETIZATION_CORE<T>(input), RESIDUAL_STAR<T>(h()), LAPLACIAN_OPERATOR<T>(m()) {
		A.resize(n(), n());
		A.reserve(Eigen::Matrix<T, Eigen::Dynamic, 1>::Constant(n(), 5));
		R.resize(n());
		X.resize(n());
		Z.resize(n());

		DISCRETIZATION_CORE<T>::solve_core();

		if (DISCRETIZATION_CORE<T>::input()[0].size() == 1) {
			std::cout << "OUTPUTTING ." << std::flush;
			double timestamp = omp_get_wtime();
			DISCRETIZATION_CORE<T>::output();
			std::cout << ".. TOOK (" << omp_get_wtime() - timestamp << "s)" << std::endl;
			return;
		}
		if (!std::isdigit(static_cast<unsigned char>(DISCRETIZATION_CORE<T>::input()[0][1]))) {
			throw std::string("SELECTION DIGIT FOR SOLVER IS NON DIGIT!");
		}
		int D = std::stoi(DISCRETIZATION_CORE<T>::input()[0].substr(1, DISCRETIZATION_CORE<T>::input()[0].size() - 1));
		if (D == 2) {
			NEWTON_CLASSIC<T>(static_cast<DISCRETIZATION_CORE<T>*>(this), Z);
		}
		/*else if (D == 1) {
		DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT<T>(static_cast<DISCRETIZATION_CORE<T>*>(this), Z);
		}*/
		else if (D == 3) {
			NEWTON_DAMPED_TRIVIAL<T>(static_cast<DISCRETIZATION_CORE<T>*>(this), Z);
		}
		else if (D == 4) {
			NEWTON_DAMPED_BACKTRACKING<T>(static_cast<DISCRETIZATION_CORE<T>*>(this), Z);
		}
		else {
			std::stringstream tmp;
			tmp << "DISCRETIZATION CORE '" << D << "' IS NOT SUPPORTED!";
			throw tmp.str();
		}
	}
private:
	//GETTER-METHODS:
	virtual std::vector<T> getSTENCILS(unsigned int i) { return DISCRETIZATION_CORE<T>::stencils()[i]; }
	virtual T getZ(unsigned int i) { return Z(i); }
	virtual bool stop() { return false; }

	//SETTER-METHODS:
	virtual void insertRESIDUAL(unsigned int i, T value) { R(i) = value; }
	virtual void insertSYSTEM_MATRIX(unsigned int i, unsigned int j, T value) { A.coeffRef(i, j) = value; }

	//CALCULATION-METHODS:
	virtual std::vector<T> calculateSTENCIL(unsigned int i, unsigned int j) {
		return { DISCRETIZATION_CORE<T>::surface()[i][j + 1],
			DISCRETIZATION_CORE<T>::surface()[i + 1][j], DISCRETIZATION_CORE<T>::surface()[i + 1][j + 1], DISCRETIZATION_CORE<T>::surface()[i + 1][j + 2],
			DISCRETIZATION_CORE<T>::surface()[i + 2][j + 1] };
	}
	virtual T calculateRESIDUAL(const std::vector<T> & stencil) { return RESIDUAL_STAR<T>::calculateRESIDUAL(stencil); }
	virtual void solve() {
		LAPLACIAN_OPERATOR<T>::calcualteSYSTEM_MATRIX(DISCRETIZATION_CORE<T>::surface());

		A *= 1 / (h() * h());

		Eigen::BiCGSTAB<MT> solver;
		solver.setTolerance(10e-7);
		solver.compute((-1)*A);
		X = solver.solve(R);
		if (X != X) {
			DISCRETIZATION_CORE<T>::output();
			throw std::string("STOPPING LAPLACE BECAUSE OF NON CONVERGENCE! TRY DIFFERENT BOUNDARY VALUES OR METHOD!");
		}

		Z = X;
	}
};