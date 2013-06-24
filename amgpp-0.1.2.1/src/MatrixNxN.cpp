/**
 * \file src/MatrixNxN.cpp
 * \brief Implementation of a quadratic NxN matrix.
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

#include "MatrixNxN.h"


/**
 * \brief Reading the matrix from a file.
 *
 * \param file The name of the input file.
 * \return void
 * \exception std::runtime_error Input error.
 *
 * This function reads a matrix from the specified file \a file. The file has to contain the
 * matrix data in the following format:
 *
 *	\code
 *	#size
 *	m(0,0)
 *	m(0,1)
 *	m(0,2)
 *	...
 *	m(1,0)
 *	m(1,1)
 *	...
 *	\endcode
 *
 * where \f$ m(i,j), i,j \in [0..N-1], \f$ specify the matrix elements. In case the output file
 * could not be opened or an input error occured, a \a std::runtime_error is thrown.
 */
void MatrixNxN::read(const char* file)
{
	std::ifstream in(file, std::ifstream::in);
	if (!in.is_open()) {
		throw std::runtime_error("File could not be opened!");
	}

	std::size_t msize(0);
	if (!(in >> msize) || msize == 0) {
		throw std::runtime_error("Matrix size could not be extracted!");
	}

	MatrixNxN tmp(msize);

	const std::size_t sqrsize(msize*msize);
	for (std::size_t i=0; i<sqrsize; ++i) {
		if (!(in >> tmp.v_[i])) {
			throw std::runtime_error("Error during matrix extraction!");
		}
	}

	swap(tmp);

	in.close();
}


/**
 * \brief Writing the matrix to a file.
 *
 * \param file The name of the output file.
 * \param prec The number of non-zero digits displayed in the output file.
 * \return void
 * \exception std::runtime_error Output error.
 *
 * This function writes the matrix to the specified file \a file using the following format:
 *
 *	\code
 *	#size
 *	m(0,0)
 *	m(0,1)
 *	m(0,2)
 *	...
 *	m(1,0)
 *	m(1,1)
 *	...
 *	\endcode
 *
 * where \f$ m(i,j), i,j \in [0..N-1], \f$ specify the matrix elements. In case the output file
 * could not be opened, a \a std::runtime_error is thrown.
 *
 * \b Note: All previous data is replaced by the matrix data!
 */
void MatrixNxN::write(const char* file, std::streamsize prec) const
{
	std::ofstream out(file, std::ofstream::out | std::ostream::trunc);
	if (!out.is_open()) {
		throw std::runtime_error("File could not be opened!");
	}

	out << size_ << "\n";

	out.precision(prec);

	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		out << v_[i] << "\n";

	out.close();
}
