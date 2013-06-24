/**
 * \file src/GaussianElimination.cpp
 * \brief A direct Gaussian elimination for arbitrary linear systems.
 * 
 * Copyright 2006, 2007, 2008, 2009 Klaus Iglberger
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

#include "GaussianElimination.h"
#include "Log.h"
#include <vector>


/**
 * \brief The default constructor for the GaussianElimination class.
 */
GaussianElimination::GaussianElimination() :
	A_(),
	b_()
{}


/**
 * \brief TODO
 *
 * \param contacts TODO
 * \return void
 *
 * TODO
 */
bool GaussianElimination::solve(const MatrixCRS& A, const VectorN& b, VectorN& x)
{
	const std::size_t n(b.size());

	assert(A.rows() == A.columns() && A.rows() == n && "LSE data dimensions do not match.");

	// allocate helper data
	A_ = A;
	b_ = b;
	x.resize(n, false);
	
	std::size_t pi, pj;
	lastPrecision_ = 0.0;

	// Initializing the pivot vector
	std::vector<std::size_t> p(n);
	for (std::size_t j=0; j<n; ++j) {
		p[j] = j;
	}

	// Gaussian elimination
	for (std::size_t j=0; j<n; ++j)
	{
		std::size_t max(j);
		double max_val(std::fabs(A_[p[max]*n+j]));

		// Partial search for pivot
		for (std::size_t i=j+1; i<n; ++i) {
			if (std::fabs(A_[p[i]*n+j]) > max_val) {
				max = i;
				max_val = std::fabs(A_[p[max]*n+j]);
			}
		}

		// Swap rows so that pivot is on diagonal
		std::swap(p[max], p[j]);
		pj = p[j];

		if (A_[pj*n+j] != 0.0)
		{
			// Eliminate column below diagonal
			for (std::size_t i=j+1; i<n; ++i)
			{
				pi = p[i];
				double f = A_[pi*n+j] / A_[pj*n+j];

				A_[pi*n+j] = 0.0;
				for (std::size_t k=j+1; k<n; ++k) {
					A_[pi*n+k] -= A_[pj*n+k] * f;
				}

				b_[pi] -= b_[pj] * f;
			}
		}
		else {
			// Assert column is zero below diagonal: which should be true due to the pivot search
			for (std::size_t i=j+1; i<n; ++i) {
				assert(A_[p[i]*n+j] == 0.0 && "Fatal error in Gaussian elimination.");
			}
		}
	}

	// Backward substitution
	for (std::size_t i=n-1; i<n; --i)
	{
		pi = p[i];
		double rhs = b_[pi];

		for (std::size_t j=i+1; j<n; ++j) {
			rhs -= x[j] * A_[pi*n+j];
		}

		if (A_[pi*n+i] != 0.0) {
			x[i] = rhs / A_[pi*n+i];
		}
		else {
			// this will introduce errors in the solution
			x[i] = 0.0;
			lastPrecision_ = std::max(lastPrecision_, std::fabs(rhs));
		}
	}

	if (lastPrecision_ < threshold_) {
		if (loglevel_ <= INFO)
			std::cout << "Solved the linear system using Gaussian elimination." << std::endl;
	}
	else if (loglevel_ <= WARNING)
		std::cout << "WARNING: Did not solve the linear system within accuracy. (" << lastPrecision_ << ")" << std::endl;

	lastIterations_ = 1;

	return lastPrecision_ < threshold_;
}
