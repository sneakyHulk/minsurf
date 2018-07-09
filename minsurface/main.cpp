#include <cstdio>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <exception>



#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <omp.h>
#include <Eigen/SparseCore>
#include "DISCRETIZATIONS.h"
#include <unistd.h>

using namespace Eigen;

void prepare_input(std::vector<std::string>& input, char * filename) {
	try {
		std::ifstream ifs(filename, std::ifstream::in);
		while (ifs.good()) {
			std::string tmp;
			char c = ifs.get();
			while (c != EOF && c != '\n') {
				tmp.push_back(c);
				c = ifs.get();
			}
			input.push_back(tmp);
		}
	}
	catch (std::exception e) {
		std::cout << e.what() << std::endl;
	}

	//KOMMERNTARE + LEERE ZEILEN ABFANGEN:
	for (unsigned int i = 0; i < input.size();) {
		for (unsigned int j = 0; j < input[i].size(); j++) {
			if (input[i][j] == '/') {
				if (++j < input[i].size()) {
					if (input[i][j] == '/') {
						input[i].erase(input[i].begin() + (j - 1), input[i].end());
						continue;
					}
				}
			}
			if (input[i][j] == ' ' || input[i][j] == '\r') {
				input[i].erase(input[i].begin() + j);
				j--;
			}
		}
		if (input[i].empty()) {
			input.erase(input.begin() + i);
			continue;
		}
		i++;
	}
}




int main(int argc, char *argv[]) {
	std::vector<std::string> input;
	double timestamp = omp_get_wtime();
	if (argc == 2) {
		prepare_input(input, argv[1]);
		//DISCRETIZATION_CORE:
		try {
			if (input.empty()) {
				throw std::string("INPUTFILE IS EMPTY!");
			}
			if (!std::isdigit(static_cast<unsigned char>(input[0][0]))) {
				throw std::string("SELECTION DIGIT FOR SOLVER IS NON DIGIT!");
			}
			int D = input[0][0] - '0';

			if (input.size() <= 1) {
				throw std::string("DATA TYPE NOT ENTERED!");
			}
			std::string type = input[1];

			if (type == "float") {
				if (D == 2) {
					NEWTON_CLASSIC<float> x(input);
				}
				/*else if (D == 1) {
				DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT<double> x(input);
				}*/
				else if (D == 0) {
					LAPLACE<float> x(input);
				}
				else if (D == 3) {
					NEWTON_DAMPED_TRIVIAL<float> x(input);
				}
				else if (D == 4) {
					NEWTON_DAMPED_BACKTRACKING<float> x(input);
				}
				else {
					std::stringstream tmp;
					tmp << "DISCRETIZATION CORE '" << D << "' IS NOT SUPPORTED!" << std::endl;
					throw tmp.str();
				}
			}
			else if (type == "double") {
				if (D == 2) {
					NEWTON_CLASSIC<double> x(input);
				}
				/*else if (D == 1) {
				DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT<double> x(input);
				}*/
				else if (D == 0) {
					LAPLACE<double> x(input);
				}
				else if (D == 3) {
					NEWTON_DAMPED_TRIVIAL<double> x(input);
				}
				else if (D == 4) {
					NEWTON_DAMPED_BACKTRACKING<double> x(input);
				}
				else {
					std::stringstream tmp;
					tmp << "DISCRETIZATION CORE '" << D << "' IS NOT SUPPORTED!" << std::endl;
					throw tmp.str();
				}
			}
			else if (type == "longdouble") {
				if (D == 2) {
					NEWTON_CLASSIC<long double> x(input);
				}
				/*else if (D == 1) {
				DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT<double> x(input);
				}*/
				else if (D == 0) {
					LAPLACE<long double> x(input);
				}
				else if (D == 3) {
					NEWTON_DAMPED_TRIVIAL<long double> x(input);
				}
				else if (D == 4) {
					NEWTON_DAMPED_BACKTRACKING<long double> x(input);
				}
				else {
					std::stringstream tmp;
					tmp << "DISCRETIZATION CORE '" << D << "' IS NOT SUPPORTED!" << std::endl;
					throw tmp.str();
				}
			}
			else {
				std::stringstream tmp;
				tmp << "DATA TYPE '" << type << "' IS NOT SUPPORTED!" << std::endl;
				throw tmp.str();
			}
		}
		catch (std::string e) {
			std::cout << e << std::endl;
		}
	}
	else if (argc == 3) {
		prepare_input(input, argv[1]);
		//DISCRETIZATION_CORE:
		try {
			if (input.empty()) {
				throw std::string("INPUTFILE IS EMPTY!");
			}
			if (!std::isdigit(static_cast<unsigned char>(input[0][0]))) {
				throw std::string("SELECTION DIGIT FOR SOLVER IS NON DIGIT!");
			}
			int D = input[0][0] - '0';

			if (input.size() <= 1) {
				throw std::string("DATA TYPE NOT ENTERED!");
			}
			std::string type = input[1];

			if (type == "float") {
				if (D == 2) {
					NEWTON_CLASSIC<float> x(input, argv[2]);
				}
				/*else if (D == 1) {
				DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT<double> x(input);
				}*/
				else if (D == 3) {
					NEWTON_DAMPED_TRIVIAL<float> x(input, argv[2]);
				}
				else if (D == 4) {
					NEWTON_DAMPED_BACKTRACKING<float> x(input, argv[2]);
				}
				else {
					std::stringstream tmp;
					tmp << "DISCRETIZATION CORE '" << D << "' IS NOT SUPPORTED!" << std::endl;
					throw tmp.str();
				}
			}
			else if (type == "double") {
				if (D == 2) {
					NEWTON_CLASSIC<double> x(input, argv[2]);
				}
				/*else if (D == 1) {
				DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT<double> x(input);
				}*/
				else if (D == 3) {
					NEWTON_DAMPED_TRIVIAL<double> x(input, argv[2]);
				}
				else if (D == 4) {
					NEWTON_DAMPED_BACKTRACKING<double> x(input, argv[2]);
				}
				else {
					std::stringstream tmp;
					tmp << "DISCRETIZATION CORE '" << D << "' IS NOT SUPPORTED!" << std::endl;
					throw tmp.str();
				}
			}
			else if (type == "longdouble") {
				if (D == 2) {
					NEWTON_CLASSIC<long double> x(input, argv[2]);
				}
				/*else if (D == 1) {
				DISCRETIZATION_3x3_SOLVER_EIGEN_JACOBIAN_HARD_DERIVATIVE_DIFFERENCEQUOTIENT<double> x(input);
				}*/
				else if (D == 3) {
					NEWTON_DAMPED_TRIVIAL<long double> x(input, argv[2]);
				}
				else if (D == 4) {
					NEWTON_DAMPED_BACKTRACKING<long double> x(input, argv[2]);
				}
				else {
					std::stringstream tmp;
					tmp << "DISCRETIZATION CORE '" << D << "' IS NOT SUPPORTED!" << std::endl;
					throw tmp.str();
				}
			}
			else {
				std::stringstream tmp;
				tmp << "DATA TYPE '" << type << "' IS NOT SUPPORTED!" << std::endl;
				throw tmp.str();
			}
		}
		catch (std::string e) {
			std::cout << e << std::endl;
		}
	}
	else {
		std::cout << "MISSING INPUT FILE, TRY: ./main <config-FILE> OR sudo ./main <config-FILE> <output-FILE>" << std::endl;
	}

	std::cout << "TOTAL TIME: \t" << omp_get_wtime() - timestamp << std::endl;
	return 0;
}