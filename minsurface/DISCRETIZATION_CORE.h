#pragma once

#define VTK

#ifdef VTK
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPointData.h>
#include <vtkPolyDataNormals.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#endif // VTK

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <typeinfo>

#include "INTERPRETER.h"

template<class D>
D test_on_input_value(const char * what, unsigned int pos_in_input, const std::vector<std::string>& input) {
	D result;
	try {
		if (input.size() <= pos_in_input) {
			std::stringstream tmp;
			tmp << what << " NOT ENTERED!";
			throw tmp.str();
		}
		if (std::is_same<D, int>::value) {
			result = std::stoi(input[pos_in_input]);
		}
		else if (std::is_same<D, float>::value) {
			result = std::stof(input[pos_in_input]);
		}
		else if (std::is_same<D, double>::value) {
			result = std::stod(input[pos_in_input]);
		}
		else if (std::is_same<D, long double>::value) {
			result = std::stold(input[pos_in_input]);
		}
		else {
			std::stringstream tmp;
			tmp << "DATA TYPE '" << typeid(D).name() << "' IS NOT SUPPORTED!" << std::endl;
			throw tmp.str();
		}
		if (!(result > 0)) {
			std::stringstream tmp;
			tmp << what << " HAS TO BE LARGER THAN 0, IS: '" << result << "'!";
			throw tmp.str();
		}
	}
	catch (std::invalid_argument e) {
		std::stringstream tmp;
		tmp << "ERROR AT CONVERTING TO '" << typeid(D).name() << "' FOR " << what << " '" << e.what() << "'!";
		throw tmp.str();
	}
	catch (std::out_of_range e) {
		std::stringstream tmp;
		tmp << "ERROR AT CONVERTING TO '" << typeid(D).name() << "' FOR " << what << " '" << e.what() << "'!";
		throw tmp.str();
	}
	return result;
}

template <class T>
class  DISCRETIZATION_CORE {
	int _s;
	unsigned int _n;
	unsigned int _m;
	T _h;
	std::vector<std::vector<T>> _surface;
	std::vector<std::vector<T>> _stencils;
	std::string _path;
	std::vector<std::string> _input;
	unsigned int _count;

protected:
	DISCRETIZATION_CORE(std::vector<std::string>& input) : _input(input), _count(1) {
		_s = test_on_input_value<int>("GRID POINTS PER ROW", 2, _input); //GRID-POINTS PER ROW
		_h = 1. / (_s - 1);
		_m = _s - 2;
		_n = _m * _m;

		//OUTPUT-PATH:
		if (_input.size() <= 3) {
			std::stringstream tmp;
			tmp << "OUTPUT-PATH NOT ENTERED!";
			throw tmp.str();
		}
		_path = _input[3];
		//NORTH:	4*x^4-6*x^3+3*x^2-x
		std::cout << "NORTH-";
		if (_input.size() <= 4) {
			std::stringstream tmp;
			tmp << "NORTH FUNCTION NOT ENTERED!";
			throw tmp.str();
		}
		INTERPRETER north = INTERPRETER(_input[4]);
		//EAST:		sin(3.1415*x^(1/2))
		std::cout << "EAST-";
		if (_input.size() <= 5) {
			std::stringstream tmp;
			tmp << "EAST FUNCTION NOT ENTERED!";
			throw tmp.str();
		}
		INTERPRETER east = INTERPRETER(_input[5]);
		//SOUTH:	x^2
		std::cout << "SOUTH-";
		if (_input.size() <= 6) {
			std::stringstream tmp;
			tmp << "SOUTH FUNCTION NOT ENTERED!";
			throw tmp.str();
		}
		INTERPRETER south = INTERPRETER(_input[6]);
		//WEST:		arctan(0.9*x)+1-1.733*x^3
		std::cout << "WEST-";
		if (_input.size() <= 7) {
			std::stringstream tmp;
			tmp << "WEST FUNCTION NOT ENTERED!";
			throw tmp.str();
		}
		INTERPRETER west = INTERPRETER(_input[7]);
		
		_surface.resize(_s);
		for (unsigned int i = 0; i < _s; i++) {
			_surface[i].resize(_s);
		}
		_stencils.resize(_n);
		for (unsigned int i = 0; i < (_s /*- 1*/); i++) {
			_surface[0][i] = north.calculate({ i * 1.0 / (_s - 1) });
			_surface[i][_s - 1] = east.calculate({ i * 1.0 / (_s - 1) });
			_surface[_s - 1][_s - 1 - i] = south.calculate({ i * 1.0 / (_s - 1) });
			_surface[_s - 1 - i][0] = west.calculate({ i * 1.0 / (_s - 1) });
		}

		double north_west = std::abs(north.calculate({ 0 }) - west.calculate({ 1 }));
		double north_east = std::abs(north.calculate({ 1 }) - east.calculate({ 0 }));
		double south_east = std::abs(south.calculate({ 0 }) - east.calculate({ 1 }));
		double south_west = std::abs(south.calculate({ 1 }) - west.calculate({ 0 }));

		if (north_west > 0.1) { std::cout << "BOUNDARY VALUES IN THE NORTH-WEST CORNER DO NOT MATCH! TOTAL DIFFERENCE: " << north_west << std::endl; }
		if (north_east > 0.1) { std::cout << "BOUNDARY VALUES IN THE NORTH-EAST CORNER DO NOT MATCH! TOTAL DIFFERENCE: " << north_east << std::endl; }
		if (south_east > 0.1) { std::cout << "BOUNDARY VALUES IN THE SOUTH-EAST CORNER DO NOT MATCH! TOTAL DIFFERENCE: " << south_east << std::endl; }
		if (south_west > 0.1) { std::cout << "BOUNDARY VALUES IN THE SOUTH-WEST CORNER DO NOT MATCH! TOTAL DIFFERENCE: " << south_west << std::endl; }

		if (north_west > 0.1 || north_east > 0.1 || south_east > 0.1 || south_west > 0.1) {
			std::string tmp;
			std::cout << "CONTINUE? (y/n) ";
			std::cin >> tmp;
			while (tmp != "y" && tmp != "n") {
				std::cout << "CONTINUE? (y/n) ";
				std::cin >> tmp;
			}
			if (tmp == "n") {
				throw std::string("BREAK!");
			}
		}
		/*
		double max_slope;
		double change_in_slope;
		std::vector<double> slopeNORTH(_s - 2); slopeNORTH.reserve(_s - 2);
		std::vector<double> slopeEAST(_s - 2); slopeEAST.reserve(_s - 2);
		std::vector<double> slopeSOUTH(_s - 2); slopeSOUTH.reserve(_s - 2);
		std::vector<double> slopeWEST(_s - 2); slopeWEST.reserve(_s - 2);

		for (unsigned int i = 1; i < (_s - 1); i++) {
		slopeNORTH[i] = (_surface[0][i] - _surface[0][i - 1]) / _h;
		slopeEAST[i] = (_surface[i][_s - 1] - _surface[i - 1][_s - 1]) / _h;
		slopeSOUTH[i] = (_surface[_s - 1][i] - _surface[_s - 1][i - 1]) / _h;
		slopeWEST[i] = (_surface[i][0] - _surface[i - 1][0]) / _h;
		}

		double max_change_in_slope = 0;

		for (unsigned i = 1; i < slopeNORTH.size(); i++) {
		change_in_slope = slopeNORTH[i - 1] - slopeNORTH[i];
		if (max_change_in_slope < change_in_slope) {
		max_change_in_slope = change_in_slope;
		}
		//std::cout << change_in_slope << std::endl;
		change_in_slope = slopeEAST[i - 1] - slopeEAST[i];
		if (max_change_in_slope < change_in_slope) {
		max_change_in_slope = change_in_slope;
		}
		//std::cout << change_in_slope << std::endl;
		change_in_slope = slopeSOUTH[i - 1] - slopeSOUTH[i];
		if (max_change_in_slope < change_in_slope) {
		max_change_in_slope = change_in_slope;
		}
		//std::cout << change_in_slope << std::endl;
		change_in_slope = slopeWEST[i - 1] - slopeWEST[i];
		if (max_change_in_slope < change_in_slope) {
		max_change_in_slope = change_in_slope;
		}
		//std::cout << change_in_slope << std::endl;
		}

		max_slope = 0;

		max_slope = std::max<double>(std::abs(*std::max_element(slopeNORTH.begin(), slopeNORTH.end())), std::abs(*std::max_element(slopeNORTH.begin(), slopeNORTH.end())));
		if (std::max<double>(std::abs(*std::max_element(slopeEAST.begin(), slopeEAST.end())), std::abs(*std::max_element(slopeEAST.begin(), slopeEAST.end()))) > max_slope) {
		max_slope = std::max<double>(std::abs(*std::max_element(slopeEAST.begin(), slopeEAST.end())), std::abs(*std::max_element(slopeEAST.begin(), slopeEAST.end())));
		}
		if (max_slope < std::max<double>(std::abs(*std::max_element(slopeSOUTH.begin(), slopeSOUTH.end())), std::abs(*std::max_element(slopeSOUTH.begin(), slopeSOUTH.end())))) {
		max_slope = std::max<double>(std::abs(*std::max_element(slopeSOUTH.begin(), slopeSOUTH.end())), std::abs(*std::max_element(slopeSOUTH.begin(), slopeSOUTH.end())));
		}
		if (max_slope <std::max<double>(std::abs(*std::max_element(slopeWEST.begin(), slopeWEST.end())), std::abs(*std::max_element(slopeWEST.begin(), slopeWEST.end())))) {
		max_slope = std::max<double>(std::abs(*std::max_element(slopeWEST.begin(), slopeWEST.end())), std::abs(*std::max_element(slopeWEST.begin(), slopeWEST.end())));
		}

		std::cout << "MAXIMUM SLOPE: " << max_slope << std::endl;
		std::cout << "MAXIMUM CHANGE IN SLOPE: " << max_change_in_slope << std::endl;
		std::cout << "-------------------------------------------------------------------------" << std::endl;

		std::string tmp;
		std::cout << "CONTINUE? (y/n) ";
		std::cin >> tmp;
		while (tmp != "y" && tmp != "n") {
		std::cout << "CONTINUE? (y/n) ";
		std::cin >> tmp;
		}
		if (tmp == "n") {
		throw std::string("BREAK!");
		}
		*/
	}

	DISCRETIZATION_CORE(std::vector<std::string>& input, char * vtk_filename) : _input(input), _count(1) {
		//OUTPUT-PATH:
		if (_input.size() <= 3) {
			std::stringstream tmp;
			tmp << "OUTPUT-PATH NOT ENTERED!";
			throw tmp.str();
		}
		_path = _input[3];
		read_vtp(vtk_filename);
		_h = 1. / (_s - 1);
		_m = _s - 2;
		_n = _m * _m;
		_stencils.resize(_n);
	}

	//GETTER-METHODS:
	unsigned int s() { return _s; }
	unsigned int n() { return _n; }
	unsigned int m() { return _m; }
	T h() { return _h; }
	std::vector<std::vector<T>>& surface() { return _surface; }
	std::vector<std::vector<T>>& stencils() { return _stencils; }
	std::vector<std::string>& input() { return _input; }
	unsigned int& count() { return _count; }
	virtual T getZ(unsigned int i) = 0;
	virtual bool stop() = 0;

	//SETTER-METHODS:
	virtual void insertRESIDUAL(unsigned int i, T value) = 0;

	//CALCULATION-METHODS:
	virtual std::vector<T> calculateSTENCIL(unsigned int i, unsigned int j) = 0;
	virtual T calculateRESIDUAL(const std::vector<T>& stencil) = 0;
	virtual void solve() = 0;
	void solve_core() {
		do {
			updateSurface();

#pragma omp parallel for
			for (int i = 0; i < _m; i++) {
				for (unsigned int j = 0; j < _m; j++) {
					unsigned int pos = i * _m + j;
					std::vector<T> stencil = calculateSTENCIL(i, j);
					_stencils[pos] = stencil;
					insertRESIDUAL(pos, calculateRESIDUAL(stencil));
				}
			}

			solve();

		} while (stop());

		updateSurface();
	}
	//Calculates the arithmetic mean over the boundary values. The result are the values of the start vector for the newton iteration:
	T calculate_arithmetic_mean() {
		T result = 0;
#pragma omp parallel for reduction(+:result)
		for (int i = 0; i < (_s - 1); i++) {
			result += _surface[0][i] + _surface[i][0] + _surface[_s - 1][i] + _surface[i][_s - 1];
		}
		result /= 4 * (_s - 1);
		return result;
	}

	void output() {
#ifdef VTK
		vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::Take(custom_reader());

		vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
#if VTK_MAJOR_VERSION <= 5
		writer->SetInput(polydata);
#else
		writer->SetInputData(polydata);
#endif
		writer->SetFileName(_path.c_str());
		writer->Write();
#else
		std::ofstream out(_path);
		for (unsigned int i = 0; i < _surface.size(); i++) {
			for (unsigned int j = 0; j < _surface[i].size(); j++) {
				out << _surface[i][j] << " ";
			}
			out << std::endl;
		}
#endif // VTK
	}

private:
	void updateSurface() {
#pragma omp parallel for
		for (int i = 1; i < _s - 1; i++) {
			for (unsigned int j = 1; j < _s - 1; j++) {
				unsigned int pos = (i - 1) * _m + (j - 1);
				_surface[i][j] = getZ(pos);
			}
		}
	}

#ifdef VTK
	vtkPolyData * custom_reader() {
		vtkIdType number_of_points, number_of_triangles, N;
		N = _s;
		number_of_points = N * N;
		number_of_triangles = 2 * (N - 1)*(N - 1);
		vtkSmartPointer<vtkPoints> points
			= vtkSmartPointer<vtkPoints>::New();
		points->SetNumberOfPoints(number_of_points);

		//punkte einlesen
		for (vtkIdType i = 0; i < N; i++) {
			for (vtkIdType j = 0; j < N; j++) {
				double x, y, z;
				x = i * _h;
				y = j * _h;
				z = _surface[i][j];
				points->SetPoint(i*N + j, x, y, z);
			}
		}

		//Dreiecke erstellen
		vtkSmartPointer<vtkCellArray> polys
			= vtkSmartPointer<vtkCellArray>::New();
		for (vtkIdType i = 0; i < N - 1; i++)
		{
			for (vtkIdType j = 0; j < N - 1; j++)
			{
				polys->InsertNextCell(3);
				polys->InsertCellPoint(i*N + j);
				polys->InsertCellPoint(i*N + j + 1);
				polys->InsertCellPoint((i + 1)*N + j + 1);

				polys->InsertNextCell(3);
				polys->InsertCellPoint(i*N + j);
				polys->InsertCellPoint((i + 1)*N + j);
				polys->InsertCellPoint((i + 1)*N + j + 1);
			}
		}

		vtkPolyData * polydata = vtkPolyData::New();
		polydata->SetPoints(points);
		polydata->SetPolys(polys);
		return polydata;
	}
#endif // VTK

	void read_vtp(char * filename) {
#ifdef VTK
		// Read the file
		double timestamp = omp_get_wtime();
		
		vtkSmartPointer<vtkXMLPolyDataReader> reader =
			vtkSmartPointer<vtkXMLPolyDataReader>::New();

		
		reader->SetFileName(filename);		
		reader->Update();

		std::cout << "READING '" << filename << "' ." << std::flush;

		// Extract the polydata
		vtkSmartPointer<vtkPolyData> polydata =
			reader->GetOutput();
		std::cout << "." << std::flush;

		// Get the number of points in the polydata
		unsigned int idNumPointsInFile = polydata->GetNumberOfPoints();
		_s = std::sqrt(idNumPointsInFile);

		_surface.resize(_s);

		for (unsigned int i = 0; i < _s; i++) {
			_surface[i].resize(_s);
			for (unsigned int j = 0; j < _s; j++) {
				double p[3];
				polydata->GetPoint(i*_s + j, p);
				_surface[i][j] = p[2];
			}
		}
		std::cout << ". TOOK (" << omp_get_wtime() - timestamp << "s)" << std::endl;
#else
		throw std:string("NOT IMPLEMENTED!");
		return;
#endif // VTK
	}
};