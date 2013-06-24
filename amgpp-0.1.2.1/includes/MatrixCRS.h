/**
 * \file includes/MatrixCRS.h
 * \brief Implementation of a compressed row storage datastructure for sparse
 *        matrices.
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

#ifndef MATRIXCRS_H_
#define MATRIXCRS_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>
#include <limits>
#include <assert.h>
#include "VectorN.h"


struct IndexValuePair {
	double value_;
	std::size_t index_;

	IndexValuePair();
	IndexValuePair(const double& value, std::size_t index);
	bool operator<(const IndexValuePair& ivp) const;
};


inline IndexValuePair::IndexValuePair() {
}


inline IndexValuePair::IndexValuePair(const double& value, std::size_t index) : value_(value), index_(index) {
}




/**
 * \brief Efficient implementation of a compressed row storage datastructure
 *        for \f$ M \times N \f$ sparse matrices.
 */
class MatrixCRS {
public:
	explicit inline MatrixCRS();
	explicit inline MatrixCRS(std::size_t m, std::size_t n);
	explicit inline MatrixCRS(std::size_t m, std::size_t n, std::size_t capacity);
	explicit        MatrixCRS(std::size_t m, std::size_t n, const std::vector<std::size_t>& capacities);

	inline const IndexValuePair& operator[](std::size_t index)                const;
	inline const IndexValuePair& operator()(std::size_t i, std::size_t j_rel) const;
	inline MatrixCRS&            operator*=(double rhs);
	inline MatrixCRS&            operator*=(const MatrixCRS& rhs);
	       bool                  operator==(const MatrixCRS& rhs)             const;
	inline const MatrixCRS       operator-()                                  const;

	inline std::size_t     rows()                                             const;
	inline std::size_t     columns()                                          const;
	inline std::size_t     nonzeros(std::size_t m)                            const;
	inline std::size_t     capacity()                                         const;
	inline std::size_t     capacity(std::size_t m)                            const;
	       void            reserveRowElements(std::size_t m, std::size_t capacity);
	       void            freeReservedRowElements(std::size_t m, std::size_t limit = std::numeric_limits<std::size_t>::max());
	       void            appendRowElement(std::size_t m, std::size_t n, double value);
	       void            insertRowElement(std::size_t m, std::size_t n, double value);
	inline double          multiplyRow(std::size_t i, const VectorN& v)       const;
	inline std::size_t     find(std::size_t i, std::size_t j)                 const;
	inline std::size_t     begin(std::size_t m)                               const;
	inline std::size_t     end(std::size_t m)                                 const;
	inline const double&   getValue(std::size_t index)                        const;
	inline double&         getValue(std::size_t index);
	inline void            clear();
	       const MatrixCRS getTranspose()                                     const;
	inline void            swap(MatrixCRS& m) throw();
	inline void            scaleDiagonal(double factor);

private:
	std::size_t cols_;                     //!< The current number of columns of the matrix.
	std::vector<IndexValuePair> values_;   //!< The dynamically allocated matrix elements.
	std::vector<size_t> row_indices_;      //!< The indices of the first element of each row in the values_ field.
	std::vector<size_t> row_lengths_;      //!< The number of nonzero entries for each row.

	friend const MatrixCRS operator+(const MatrixCRS& rhs, const MatrixCRS& lhs);
	friend const MatrixCRS operator-(const MatrixCRS& rhs, const MatrixCRS& lhs);
	friend const MatrixCRS operator*(const MatrixCRS& rhs, const MatrixCRS& lhs);
	friend const MatrixCRS operator*(const MatrixCRS& mat, double scalar);
	friend const MatrixCRS operator*(double scalar, const MatrixCRS& mat);
	friend const VectorN   operator*(const MatrixCRS& mat, const VectorN& vec);
};


/**
 * \brief The default constructor for MatrixCRS.
 */
inline MatrixCRS::MatrixCRS() : cols_(0), values_(), row_indices_(), row_lengths_() {
}


/**
 * \brief Constructor for a matrix of size \f$ M \times N \f$.
 *
 * \param m The number of rows of the matrix.
 * \param n The number of columns of the matrix.
 *
 * The matrix is initialized to the zero matrix and has no free capacity.
 */
inline MatrixCRS::MatrixCRS(std::size_t m, std::size_t n) : cols_(n), values_(), row_indices_(m + 1, 0), row_lengths_(m, 0) {
}


/**
 * \brief Constructor for a matrix of size \f$ M \times N \f$.
 *
 * \param m The number of rows of the matrix.
 * \param n The number of columns of the matrix.
 * \param cap The number of elements allocated in advance.
 *
 * The matrix is initialized to the zero matrix.
 */
inline MatrixCRS::MatrixCRS(std::size_t m, std::size_t n, std::size_t cap) : cols_(n), values_(cap), row_indices_(m + 1, 0), row_lengths_(m, 0) {
}


/**
 * \brief Direct read-access to the matrix elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..capacity-1]\f$.
 * \return Constant reference to the accessed value.
 */
inline const IndexValuePair& MatrixCRS::operator[](std::size_t index) const {
	return values_[index];
}


/**
 * \brief Direct read-access to the matrix elements.
 *
 * \param i Absolute row index. The index has to be in the range \f$[0..M-1]\f$.
 * \param j Relative column index. The index has to be in the range \f$[0..nonzeros(i)-1]\f$.
 * \return Constant reference to the accessed value.
 */
inline const IndexValuePair& MatrixCRS::operator()(std::size_t i, std::size_t j_rel) const {
	assert(i<rows() && j_rel<nonzeros(i) && "Invalid matrix access index");
	return values_[row_indices_[i] + j_rel];
}


/**
 * \brief Unary minus operator for the inversion of a matrix (\f$ A = -B \f$).
 *
 * \return The inverse of the matrix.
 */
inline const MatrixCRS MatrixCRS::operator-() const {
	MatrixCRS tmp(*this);

	for (std::size_t i = 0; i < tmp.rows(); ++i) {
		for (std::size_t j = tmp.begin(i); j < tmp.end(i); ++j)
			tmp.values_[j].value_ = -tmp.values_[j].value_;
	}

	return tmp;
}


/**
 * \brief Multiplication assignment operator for the multiplication between a matrix and
 * \brief a scalar value (\f$ A*=s \f$).
 *
 * \param rhs The right-hand-side scalar value for the multiplication.
 * \return Reference to the matrix.
 */
inline MatrixCRS& MatrixCRS::operator*=(double rhs) {
	for (std::size_t i = 0; i < rows(); ++i) {
		for (std::size_t j = begin(i); j < end(i); ++j)
			values_[j].value_ *= rhs;
	}

	return *this;
}


/**
 * \brief Multiplication assignment operator for the multiplication between two matrices
 * \brief (\f$ A*=B \f$).
 *
 * \param rhs The right-hand-side matrix for the multiplication.
 * \return Reference to the matrix.
 *
 * This operator works only for square matrices
 */
inline MatrixCRS& MatrixCRS::operator*=(const MatrixCRS& rhs) {
	MatrixCRS tmp((*this) * rhs);
	return this->operator=(tmp);
}


/**
 * \brief Returns the number of rows of the matrix.
 *
 * \return The number of rows of the matrix.
 */
inline std::size_t MatrixCRS::rows() const {
	return row_lengths_.size();
}


/**
 * \brief Returns the number of columns of the matrix.
 *
 * \return The number of columns of the matrix.
 */
inline std::size_t MatrixCRS::columns() const {
	return cols_;
}


/**
 * \brief Returns the number of entries in the row m.
 *
 * \param m The index of the row.
 * \return The number of entries in the row m.
 */
inline std::size_t MatrixCRS::nonzeros(std::size_t m) const {
	return row_lengths_[m];
}


/**
 * \brief Returns the total capacity of the matrix.
 *
 * \return The total capacity of the matrix.
 */
inline std::size_t MatrixCRS::capacity() const {
	return values_.size();
}


/**
 * \brief Returns the current capacity of the row m.
 *
 * \param m The index of the row.
 * \return The number of allocated entries for the row m.
 */
inline std::size_t MatrixCRS::capacity(std::size_t m) const {
	return row_indices_[m+1] - row_indices_[m];
}


/**
 * \brief Locates the entry (i, j) in the matrix.
 *
 * \param i The index of the row.
 * \param j The index of the column.
 * \return Returns the direct access index of the element if present or
 *         std::numeric_limits<std::size_t>::max() otherwise.
 *
 * The method requires that the matrix entries are ordered by their column indices. Then a binary
 * search for the element can be performed which takes logarithmic time.
 */
inline std::size_t MatrixCRS::find(std::size_t i, std::size_t j) const {
	std::size_t direct_index = std::lower_bound(values_.begin() + begin(i), values_.begin() + end(i), IndexValuePair(0, j)) - values_.begin();

	if (direct_index != end(i) && values_[direct_index].index_ == j)
		return direct_index;

	return std::numeric_limits<std::size_t>::max();
}


/**
 * \brief Returns the access index of the first entry in the m-th row.
 *
 * \param m The index of the row.
 * \return The access index of the first entry in the m-th row.
 */
inline std::size_t MatrixCRS::begin(std::size_t m) const {
	return row_indices_[m];
}


/**
 * \brief Returns the access index of the entry after the last entry in the m-th row.
 *
 * \param m The index of the row.
 * \return The access index of the entry after the last entry in the m-th row.
 */
inline std::size_t MatrixCRS::end(std::size_t m) const {
	return row_indices_[m] + row_lengths_[m];
}


/**
 * \brief Returns a read-only reference to the value of the specified element.
 *
 * \param index The direct access index of the element.
 * \return A constant reference to the value of the specified element.
 */
inline const double& MatrixCRS::getValue(std::size_t index) const {
	return values_[index].value_;
}


/**
 * \brief Returns a read-write reference to the value of the specified element.
 *
 * \param index The direct access index of the element.
 * \return A reference to the value of the specified element.
 */
inline double& MatrixCRS::getValue(std::size_t index) {
	return values_[index].value_;
}


/**
 * \brief Clearing the \f$ M \times N \f$ matrix.
 *
 * \return void
 */
inline void MatrixCRS::clear() {
	std::fill(row_indices_.begin(), row_indices_.end(), 0);
	std::fill(row_lengths_.begin(), row_lengths_.end(), 0);
}


/**
 * \brief Multiplies a matrix row with a vector.
 * \param i The index of the row.
 * \param v The vector.
 *
 * \return Returns the product of the specified matrix row with the given vector.
 */
inline double MatrixCRS::multiplyRow(std::size_t i, const VectorN& v) const {
	double accu = 0.0;
	const std::size_t j_end = end(i);

	for (std::size_t j = begin(i); j < j_end; ++j) {
		accu += values_[j].value_ * v[values_[j].index_];
	}

	return accu;
}


/**
 * \brief Swapping the contents of two matrices.
 *
 * \param m The matrix to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
inline void MatrixCRS::swap(MatrixCRS& m) throw() {
	std::swap(cols_, m.cols_);
	std::swap(values_, m.values_);
	std::swap(row_indices_, m.row_indices_);
	std::swap(row_lengths_, m.row_lengths_);
}


/**
 * \brief Scales the diagonal of the matrix.
 *
 * \param factor The multiplication factor.
 * \return void
 */
inline void MatrixCRS::scaleDiagonal(double factor) {
	const std::size_t n(std::min(rows(), columns()));
	std::size_t j;

	for (std::size_t i = 0; i < n; ++i) {
		j = find(i, i);
		if (j != std::numeric_limits<std::size_t>::max())
			getValue(j) *= factor;
	}
}


/**
 * \brief Global output operator for MxN matrices.
 *
 * The output starts with the number of rows, columns and total non-zero entries
 * each seperated by a single space. Each non-zero entry is listed on a line
 * of its own, which consists of its row and column index followed by the value
 * each seperated by a single space.
 *
 * \param os Reference to the output stream.
 * \param m Reference to a constant matrix object.
 * \return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, const MatrixCRS& m) {
	std::size_t c(0);
	for (std::size_t i = 0; i < m.rows(); ++i)
		c += m.nonzeros(i);

	os << m.rows() << " " << m.columns() << " " << c << '\n';
	for (std::size_t i = 0; i < m.rows(); ++i) {
		for (std::size_t j = m.begin(i); j < m.end(i); ++j)
			os << i << " " << m[j].index_ << " " << m[j].value_ << '\n';
	}

	return os;
}


/**
 * \brief Global input operator for MxN matrices.
 *
 * The input format begins with the number of rows, the number of columns and
 * the total number of non-zero entries. All non-zero entries are listed in
 * row-major order. Each entry consists of its row and column index followed
 * by the entry's value. All quantities must be separated by white-space of any
 * kind.
 *
 * \param is Reference to the output stream.
 * \param m Reference to a matrix object.
 * \return Reference to the output stream.
 */
inline std::istream& operator>>(std::istream& is, MatrixCRS& A) {
	std::size_t n, m, c;
	is >> n >> m >> c;

	A = MatrixCRS(n, m, c);

	int i = 0, j = -1, i_last, j_last;
	double x;
	A.reserveRowElements(0, c);

	for (int k = 0; k < c; ++k) {
		i_last = i; j_last = j;
		is >> i >> j >> x;
		
		assert(((i == i_last && j > j_last) || i > i_last) && "Input must be ordered row-wise.");
		assert(i >= 0 && i < n && j >= 0 && j < m && "Entry indices out of range.");
		
		if (i > i_last) {
			for (int l = i_last; l < i; ++l)
				A.freeReservedRowElements(l);
		}

		A.appendRowElement(i, j, x);
	}

	return is;
}


/**
 * \brief Checks the given matrix for not-a-number elements.
 *
 * \param m The matrix to be checked for not-a-number elements.
 * \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
 */
inline bool isnan(const MatrixCRS& m)
{
	for (std::size_t i = 0; i < m.rows(); ++i) {
		for (std::size_t j = m.begin(i); j < m.end(i); ++j)
			if (isnan(m[j].value_)) return true;
	}

	return false;
}


/**
 * \brief Returns a matrix containing the absolute values of each single element of \a m.
 *
 * \param m The floating point input matrix.
 * \return The absolute value of each single element of \a m.
 *
 * The \a fabs function calculates the absolute value of each element of the input matrix
 * \a m. This function can only be applied to  floating point matrices. For matrices of
 * integral data type, the pe::abs(const MatrixCRS&) function can be used.
 */
inline const MatrixCRS fabs(const MatrixCRS& m) {
	MatrixCRS tmp(m);

	for (std::size_t i = 0; i < tmp.rows(); ++i) {
		for (std::size_t j = tmp.begin(i); j < tmp.end(i); ++j)
			tmp.getValue(j) = std::fabs(tmp[j].value_);
	}

	return tmp;
}


/**
 * \brief Swapping the contents of two matrices.
 *
 * \param a The first matrix to be swapped.
 * \param b The second matrix to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
inline void swap(MatrixCRS& a, MatrixCRS& b) throw() {
	a.swap(b);
}


/**
 * \brief Multiplication operator for the multiplication of a matrix and a scalar value
 * \brief (\f$ A=B*s \f$).
 *
 * \param mat The left-hand side matrix for the multiplication.
 * \param scalar The right-hand side scalar value for the multiplication.
 * \return The scaled result matrix.
 */
inline const MatrixCRS operator*(const MatrixCRS& mat, double scalar) {
	MatrixCRS tmp(mat);

	for (std::size_t i=0; i < tmp.rows(); ++i) {
		for (std::size_t j = tmp.begin(i); j < tmp.end(i); ++j)
			tmp.getValue(j) = tmp[j].value_ * scalar;
	}

	return tmp;
}


/**
 * \brief Multiplication operator for the multiplication of a scalar value and a matrix
 * \brief (\f$ A=s*B \f$).
 * \ingroup MatrixCRS
 *
 * \param scalar The left-hand side scalar value for the multiplication.
 * \param mat The right-hand side matrix for the multiplication.
 * \return The scaled result matrix.
 */
inline const MatrixCRS operator*(double scalar, const MatrixCRS& mat) {
	MatrixCRS tmp(mat);

	for (std::size_t i=0; i < tmp.rows(); ++i) {
		for (std::size_t j = tmp.begin(i); j < tmp.end(i); ++j)
			tmp.getValue(j) = scalar * tmp[j].value_;
	}

	return tmp;
}


/**
 * \brief Multiplication operator for the multiplication of a matrix and a vector
 * \brief (\f$ \vec{a}=B*\vec{c} \f$).
 *
 * \param mat The left-hand side matrix for the multiplication.
 * \param vec The right-hand-side vector for the multiplication.
 * \return The resulting vector.
 */
inline const VectorN operator*(const MatrixCRS& mat, const VectorN& vec) {
	assert(vec.size() == mat.columns() && "Matrix and vector sizes do not match");

	VectorN tmp(mat.rows());

	for (std::size_t i=0; i < mat.rows(); ++i) {
		tmp[i] = 0;
		for (std::size_t j = mat.begin(i); j < mat.end(i); ++j) {
			tmp[i] += mat[j].value_ * vec[mat[j].index_];
		}
	}

	return tmp;
}

#endif
