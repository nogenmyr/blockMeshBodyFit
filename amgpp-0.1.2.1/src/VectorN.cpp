/**
 * \file src/VectorN.cpp
 * \brief Implementation of an arbitrary sized vector.
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

#include "VectorN.h"


/**
 * \brief Changing the size of the vector.
 *
 * \param n The new size of the vector.
 * \param preserve \a true if the old values of the vector should be preserved, \a false if not.
 * \return void
 *
 * This function resizes the vector using the given size to \a n. During this operation, new
 * dynamic memory may be allocated in case the capacity of the vector is too small. Therefore
 * this function potentially changes all vector elements. In order to preserve the old vector
 * values, the \a preserve flag can be set to \a true. However, new vector elements are not
 * initialized!\n
 * The following example illustrates the resize operation of a vector of size 2 to a vector of
 * size 4. The new, uninitialized elements are marked with \a x:

 \f[
 \left(\begin{array}{*{2}{c}}
 1 & 2 \\
 \end{array}\right)

 \Longrightarrow

 \left(\begin{array}{*{4}{c}}
 1 & 2 & x & x \\
 \end{array}\right)
 \f]
 */
void VectorN::resize(std::size_t n, bool preserve) {
	if (n == size_) return;

	if (preserve) {
		double* v = new double[n];
		const std::size_t minsize(std::min(n, size_));

		for (std::size_t i=0; i<minsize; ++i)
			v[i] = v_[i];

		std::swap(v_, v);
		delete [] v;
		capacity_ = n;
	}
	else if (n > capacity_) {
		double* v = new double[n];
		std::swap(v_, v);
		delete [] v;
		capacity_ = n;
	}

	size_ = n;
}
