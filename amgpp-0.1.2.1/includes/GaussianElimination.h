/**
 * \file includes/GaussianElimination.h
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

#ifndef GAUSSIANELIMINATION_H_
#define GAUSSIANELIMINATION_H_

#include "Solver.h"
#include "MatrixCRS.h"
#include "MatrixNxN.h"
#include "VectorN.h"
#include "Lse.h"


/**
 * \brief A Gaussian elimination for solving arbitrary linear systems.
 *
 * TODO
 */
class GaussianElimination : public Solver
{
public:
	explicit GaussianElimination();

	bool solve(Lse& lse);
	bool solve(const MatrixCRS& A, const VectorN& b, VectorN& x);

private:
	MatrixNxN A_;
	VectorN b_;
};


/**
 * \brief TODO
 *
 * \param lse TODO
 * \return bool
 *
 * TODO
 */
inline bool GaussianElimination::solve(Lse& lse) {
	return solve(lse.A_, lse.b_, lse.x_);
}

#endif
