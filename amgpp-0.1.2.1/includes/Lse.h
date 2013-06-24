/**
 * \file includes/Lse.h
 * \brief A data structure for linear systems of equations.
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

#ifndef LSE_H_
#define LSE_H_

#include <limits.h>
#include "MatrixCRS.h"
#include "VectorN.h"

/**
 * \brief A linear system of equations data-structure.
 *
 * Represents a linear system of equations \f$ A x + b = 0 \f$.
 */
class Lse {
public:
	MatrixCRS A_;
	VectorN   b_;
	VectorN   x_;

	double residual(std::size_t i) const;
	double residual() const;
	std::size_t size() const;
};
	
inline double Lse::residual(std::size_t i) const {
	return b_[i] - A_.multiplyRow(i, x_);
}

inline double Lse::residual() const {
	double r_max(0);
	for (std::size_t i = 0; i < size(); ++i)
		r_max = std::max(r_max, std::fabs(residual(i)));

	return r_max;
}

inline std::size_t Lse::size() const {
	return x_.size();
}

#endif
