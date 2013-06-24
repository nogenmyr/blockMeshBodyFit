/**
 * \file src/amgpp.cpp
 * \brief The amg command-line interface.
 * 
 * Copyright 2008, 2009 Tobias Preclik
 * 
 * This file is part of amgpp.
 * 
 * Amgpp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Amgpp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with amgpp.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <string.h>
#include "AmgSolver.h"
#include "GaussianElimination.h"
#include "Log.h"

unsigned int loglevel_ = WARNING;


/**
 * \brief Print out usage information.
 */
void usage() {
	std::cout << "amgtest\n";
	std::cout << "--matrix <A.dat> --rhs <b.dat>\n";
	std::cout << "--poisson <n>\n";
	std::cout << "--stencil <nw> <n> <ne> <w> <c> <e> <sw> <s> <se> <gridsize> --rhs <b.dat>\n";
	std::cout << "[--loglevel info|warning|critical|error] [--nu <presmooth> <postsmooth>]\n";
	std::cout << "[--sdep <thetapos> <thetaneg>] [--fmg] [--coarsening amg1r5|direct|standard]\n";
	std::cout << std::flush;
}


/**
 * \brief Setup a matrix from a stencil on an \f$n \times n\f$ grid.
 */
void setupMatrixFromStencil(double stencil[3][3], MatrixCRS& A, std::size_t n) {
	unsigned int N = n * n;

	std::vector<std::size_t> capacities(N, 9);
	A = MatrixCRS(N, N, capacities);

	for (std::size_t i = 0; i < N; ++i) {
		if (i >= n) {
			if (i % n != 0 && stencil[2][0] != 0.0)
				A.appendRowElement(i, i - n - 1, stencil[2][0]);
			if (stencil[2][1] != 0.0)
				A.appendRowElement(i, i - n, stencil[2][1]);
			if (i % n != n - 1 && stencil[2][2] != 0.0)
				A.appendRowElement(i, i - n + 1, stencil[2][2]);
		}

		if (i % n != 0 && stencil[1][0] != 0.0) {
			A.appendRowElement(i, i - 1, stencil[1][0]);
		}

		A.appendRowElement(i, i, stencil[1][1]);

		if (i % n != n - 1 && stencil[1][2] != 0.0) {
			A.appendRowElement(i, i + 1, stencil[1][2]);
		}

		if (i < N - n) {
			if (i % n != 0 && stencil[0][0] != 0.0)
				A.appendRowElement(i, i + n - 1, stencil[0][0]);
			if (stencil[0][1] != 0.0)
				A.appendRowElement(i, i + n, stencil[0][1]);
			if (i % n != n - 1 && stencil[0][2] != 0.0)
				A.appendRowElement(i, i + n + 1, stencil[0][2]);
		}
	}
}


/**
 * \brief Setup a Poisson problem on the unit square.
 *
 * Sets up the Poisson problem on the unit square. The right-hand-side is chosen
 * so that \f$u(x,y) = (x^2-x^4)(y^4-y^2)\f$ is the solution.
 *
 * \param system On return the system matrix and right-hand-side will be set up.
 * \param n The number of intervals in each direction.
 */
void setupPoisson(Lse& system, int n) {
	double stencil[3][3] = {
		{ 0.0,  1.0,  0.0},
		{ 1.0, -4.0,  1.0},
		{ 0.0,  1.0,  0.0}
	};

	unsigned int N = (n-1) * (n-1);
	double h = 1.0 / n;

	setupMatrixFromStencil(stencil, system.A_, n - 1);

	// setup right-hand-side
	system.b_.resize(N, false);
	for (unsigned int i = 0; i < N; ++i) {
		double x = (1 + i % (n - 1)) * h;
		double y = (1 + i / (n - 1)) * h;

		system.b_[i] = -2 * ((1 - 6 * x * x) * y * y * (1 - y * y) + (1 - 6 * y * y) * x * x * (1 - x * x)) * h * h;
	}

	// note: u(x,y) is zero on the boundary so we do not need to correct right-hand-side

	// write reference solution
	std::ofstream fout("reference-solution.dat");
	for (unsigned int i = 0; i < N; ++i) {
		double x = (1 + i % (n - 1)) * h;
		double y = (1 + i / (n - 1)) * h;

		if (i % (n-1) == 0)
			fout << std::endl;
		fout << x << " " << y << " " << (x*x - x*x*x*x)*(y*y*y*y-y*y) << std::endl;
	}
	fout.close();
}


void outputGrid(Lse& system, int nDim, double hx = 1.0, double hy = 1.0) {
	unsigned int N = system.size();
	int mDim = N / nDim;

	std::ofstream fout("solution.dat");
	for (std::size_t i = 1; i <= nDim; ++i) {
		fout << '\n';
		for (std::size_t j = 1; j <= mDim; ++j)
			fout << j * hx << " " << i * hy << " " << system.x_[(i - 1) * mDim + (j - 1)] << std::endl;
	}
	fout.close();
}


int main(int argc, char* argv[]) {
	Lse system;
	bool initialRand(false);
	int nDim = 0;
	int type = 0;
	int nu1 = 1, nu2 = 1;
	double thetaPos = 0.0, thetaNeg = 0.4;
	int coarsening = 0;
	bool rhs = false;
	bool fmg = false;

	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "--poisson") == 0 && i + 1 < argc && type == 0) {
			setupPoisson(system, atoi(argv[i + 1]));

			rhs = true;
			nDim = atoi(argv[i + 1]) - 1;
			type = 1;
			i += 1;
		}
		else if (strcmp(argv[i], "--stencil") == 0 && i + 10 < argc && type == 0) {
			double stencil[3][3];

			for (int j = 0; j < 3; ++j)
				for (int k = 0; k < 3; ++k)
					stencil[j][k] = atof(argv[i + 1 + j*3 + k]);

			int n = atoi(argv[i + 10]);
			setupMatrixFromStencil(stencil, system.A_, n);

			nDim = n;
			type = 2;
			i += 10;
		}
		else if (strcmp(argv[i], "--matrix") == 0 && i + 1 < argc && type == 0) {
			std::ifstream matrixIn(argv[i + 1]);
			if (!matrixIn) {
				std::cerr << "Error: Cannot open matrix input file " << argv[i + 1] << "." << std::endl;
				exit(EXIT_FAILURE);
			}
			matrixIn >> system.A_;
			
			nDim = 1;
			type = 3;
			i += 1;
		}
		else if (strcmp(argv[i], "--rhs") == 0 && i + 1 < argc) {
			std::ifstream rhsIn(argv[i + 1]);
			if (!rhsIn) {
				std::cerr << "Error: Cannot open right-hand-side input file " << argv[i + 1] << "." << std::endl;
				exit(EXIT_FAILURE);
			}
			rhsIn >> system.b_;

			rhs = true;
			i += 1;
		}
		else if ((strcmp(argv[i], "--log") == 0 || strcmp(argv[i], "-l") == 0) && i + 1 < argc) {
			if (strcmp(argv[i+1], "info") == 0)
				loglevel_ = INFO;
			else if (strcmp(argv[i+1], "warning") == 0)
				loglevel_ = WARNING;
			else if (strcmp(argv[i+1], "critical") == 0)
				loglevel_ = CRITICAL;
			else if (strcmp(argv[i+1], "error") == 0)
				loglevel_ = ERROR;
			else {
				std::cerr << "Error: Invalid log level  " << argv[i+1] << "." << std::endl;
				exit(EXIT_FAILURE);
			}

			i += 1;
		}
		else if (strcmp(argv[i], "--nu") == 0 && i + 2 < argc) {
			nu1 = atoi(argv[i+1]);
			nu2 = atoi(argv[i+2]);

			i += 2;
		}
		else if (strcmp(argv[i], "--sdep") == 0 && i + 2 < argc) {
			thetaPos = strtod(argv[i+1], NULL);
			thetaNeg = strtod(argv[i+2], NULL);

			i += 2;
		}
		else if (strcmp(argv[i], "--fmg") == 0) {
			fmg = true;
		}
		else if (strcmp(argv[i], "--rand") == 0) {
			initialRand = true;
		}
		else if ((strcmp(argv[i], "--coarsening") == 0 || strcmp(argv[i], "-c") == 0) && i + 1 < argc) {
			// choose coarsening strategy from coarsening_amg1r5, coarsening_direct and coarsening_standard
			if (strcmp(argv[i+1], "amg1r5") == 0)
				coarsening = 0;
			else if (strcmp(argv[i+1], "direct") == 0)
				coarsening = 1;
			else if (strcmp(argv[i+1], "standard") == 0)
				coarsening = 2;
			else {
				std::cerr << "Error: Invalid coarsening parameter " << argv[i+1] << std::endl;
				exit(EXIT_FAILURE);
			}

			i += 1;
		}
		else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
			usage();
			exit(EXIT_SUCCESS);
		}
		else {
			std::cerr << "Error: Invalid parameter " << argv[i] << "." << std::endl;
			usage();
			exit(EXIT_FAILURE);
		}
	}

	if (type == 0) {
		std::cerr << "Error: No system matrix set up." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (!rhs) {
		std::cerr << "Error: No right-hand side set up." << std::endl;
		exit(EXIT_FAILURE);
	}

	// verify that system dimensions match
	std::size_t n(system.b_.size());
	if (system.A_.rows() != n || system.A_.columns() != n) {
		std::cerr << "Error: Matrix must be square and must match right-hand-side vector." << std::endl;
		exit(EXIT_FAILURE);
	}
	
	std::cout << "System size: " << n << std::endl;

	// setup initial solution
	system.x_.resize(n, false);

	if (!initialRand)
		system.x_ = 0.0;
	else
		for (std::size_t i = 0; i < n; ++i)
			system.x_[i] = 1.0 * rand() / (RAND_MAX + 1.0);

	// setup AMG solver
	AmgSolver<GaussianElimination> amg;

	// no Full Multigrid
	amg.setStartupExecution(fmg);

	// V(1,1)-cycle
	amg.setSmoothingParameters(nu1, nu2);

	// set strong dependency thresholds for positive and negative entries
	amg.setCoarseningParameters(thetaPos, thetaNeg);

	// choose coarsening strategy from coarsening_amg1r5, coarsening_direct and coarsening_standard
	amg.setCoarseningStrategy(AmgSolver<GaussianElimination>::coarsening_amg1r5);

	// solve system
	bool solved = amg.solve(system);

	if (loglevel_ > INFO) {
		if (solved)
			std::cout << "Solved the LSE in " << amg.getLastIterations() << " multigrid cycles." << std::endl;
		else
			std::cout << "WARNING: Did not solve the LSE within accuracy. (" << amg.getLastPrecision() << ")" << std::endl;
	}

	double hx = 1.0, hy = 1.0;
	if (type == 1)
		hx = hy = 1.0 / (nDim + 1);

	if (nDim != 0)
		outputGrid(system, nDim, hx, hy);

	return 0;
}
