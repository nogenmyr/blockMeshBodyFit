/**
 * \file src/MatrixCRS.cpp
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

#include "MatrixCRS.h"


bool IndexValuePair::operator<(const IndexValuePair& ivp) const {
	return index_ < ivp.index_;
}


/**
 * \brief Constructor for a matrix of size \f$ M \times N \f$.
 *
 * \param m The number of rows of the matrix.
 * \param n The number of columns of the matrix.
 * \param capacities The number of elements which should be reserved for each row.
 *
 * The matrix is initialized to the zero matrix.
 */
MatrixCRS::MatrixCRS(std::size_t m, std::size_t n, const std::vector<std::size_t>& capacities) : cols_(n), values_(), row_indices_(m + 1), row_lengths_(m, 0) {
	assert(capacities.size() == m && "The number of rows does not match the size of the capacity vector.");

	std::size_t cap = 0;
	row_indices_[0] = 0;
	for (std::size_t i = 0; i < m; ++i) {
		row_indices_[i+1] = row_indices_[i] + capacities[i];
		cap += capacities[i]; 
	}

	values_.resize(cap);
}


/**
 * \brief Tests if the two matrices are semantically equal.
 *
 * \return Returns true if the two matrices are semantically equal, otherwise false.
 */
bool MatrixCRS::operator==(const MatrixCRS& rhs) const {
	if (rows() != rhs.rows() || columns() != rhs.columns())
		return false;

	const MatrixCRS& lhs(*this);

	for (std::size_t i = 0; i < rows(); ++i) {
		std::size_t j_lhs(lhs.begin(i)), j_rhs(rhs.begin(i));
		std::size_t j_lhs_end(lhs.end(i)), j_rhs_end(rhs.end(i));

		while (j_lhs < j_lhs_end && j_rhs < j_rhs_end) {
			if (lhs[j_lhs].index_ < rhs[j_rhs].index_) {
				if (lhs[j_lhs++].value_ != 0.0)
					return false;
			}
			else if (lhs[j_lhs].index_ > rhs[j_rhs].index_) {
				if (rhs[j_rhs++].value_ != 0.0)
					return false;
			}
			else {
				if (lhs[j_lhs++].value_ != rhs[j_rhs++].value_)
					return false;
			}
		}

		while (j_lhs < j_lhs_end) {
			if (lhs[j_lhs++].value_ != 0.0)
				return false;
		}

		while (j_rhs < j_rhs_end) {
			if (rhs[j_rhs++].value_ != 0.0)
				return false;
		}
	}

	return true;
}


/**
 * \brief Reserves at least capacity elements for the m-th row.
 *
 * \param m The index of the row.
 * \param cap The minimum target capacity of the row.
 */
void MatrixCRS::reserveRowElements(std::size_t m, std::size_t cap) {
	if (capacity(m) >= cap)
		return;

	std::size_t additional = cap - capacity(m);

	// enlarge entry storage if necessary
	if (capacity() - row_indices_[rows()] < additional) {
		values_.resize(row_indices_[rows()] + additional);
	}

	// move successive rows
	std::copy_backward(&values_[row_indices_[m+1]], &values_[row_indices_[rows()]], &values_[row_indices_[rows()] + additional]);
	for (std::size_t i = m+1; i <= rows(); ++i)
		row_indices_[i] += additional;
}


/**
 * \brief Frees reserved row elements for the m-th row.
 * 
 * \param m The index of the row.
 * \param limit The maximum number of elements freed, by default INT_MAX.
 * 
 * At most limit unused reserved row elements from the m-th row will be transfered to the
 * next row. This will take constant time unless the next row is already filled. In this case
 * the process needs time in the order of the number of entries in the next row.
 */ 
void MatrixCRS::freeReservedRowElements(std::size_t m, std::size_t limit) {
	std::size_t additional = std::min(capacity(m) - nonzeros(m), limit);
	if (additional == 0)
		return;
	
	// shrink reserved storage
	std::copy(&values_[row_indices_[m+1]], &values_[row_indices_[m+1] + row_lengths_[m+1]], &values_[row_indices_[m+1] - additional]);
	row_indices_[m+1] -= additional;
}


/**
 * \brief Appends an element to the m-th row.
 *
 * \param m The absolute index of the entry's row.
 * \param n The absolute index of the entry's column.
 * \param value The value of the entry.
 *
 * Appends an element to the m-th row under the condition that the column index is larger than
 * the index of the previous element in the row. Furthermore there must be enough reserved space
 * in the row. This allows an insertion in constant time.
 */
void MatrixCRS::appendRowElement(std::size_t m, std::size_t n, double value) {
	assert(end(m) < begin(m+1) && "Not enough reserved space in the current row.");
	assert((nonzeros(m) == 0 || n > (*this)(m, nonzeros(m)-1).index_) && "Index is not strictly increasing.");

	values_[end(m)].value_ = value;
	values_[end(m)].index_ = n;
	++row_lengths_[m];
}


/**
 * \brief Inserts an element into the m-th row.
 *
 * \param m The absolute index of the entry's row.
 * \param n The absolute index of the entry's column.
 * \param value The value of the entry.
 *
 * Inserts an element to the m-th row and reallocates space if necessary.
 */
void MatrixCRS::insertRowElement(std::size_t m, std::size_t n, double value) {
	if (end(m) >= begin(m+1))
		reserveRowElements(m, capacity(m) + 1);

	std::size_t j = end(m);
	++row_lengths_[m];

	values_[j].value_ = value;
	values_[j].index_ = n;

	while (j != begin(m) && values_[j].index_ < values_[j-1].index_) {
		std::swap(values_[j], values_[j-1]);
		--j;
	}

	assert((j == begin(m) || values_[j-1].index_ < values_[j].index_) && "Double entry detected.");
}


/**
 * \brief Calculation of the transpose of the matrix.
 *
 * \return The transpose of the matrix.
 */
const MatrixCRS MatrixCRS::getTranspose() const {
	// count the number of entries per column
	std::vector<std::size_t> column_lengths(columns(), 0);
	for (std::size_t i = 0; i < rows(); ++i)
		for (std::size_t j = begin(i); j < end(i); ++j)
			++column_lengths[values_[j].index_];

	// setup tranpose and reserve the correct number of entries per row
	MatrixCRS tmp(columns(), rows(), column_lengths);

	// append elements to rows of transpose
	for (std::size_t i = 0; i < rows(); ++i)
		for (std::size_t j = begin(i); j < end(i); ++j)
			tmp.appendRowElement(values_[j].index_, i, values_[j].value_);

	return tmp;
}


/**
 * \brief Addition operator for the addition of two matrices (\f$ A=B+C \f$).
 *
 * \param lhs The left-hand side matrix for the matrix addition.
 * \param rhs The right-hand side matrix to be added to the left-hand side matrix.
 * \return The sum of the two matrices.
 */
const MatrixCRS operator+(const MatrixCRS& lhs, const MatrixCRS& rhs) {
	assert((lhs.rows() == rhs.rows() && rhs.columns() == rhs.columns()) && "Matrix sizes do not match");

	// analyze matrices to predict the required row storage capacities
	std::vector<std::size_t> capacities(lhs.rows(), 0);
	for (std::size_t i = 0; i < lhs.rows(); ++i) {
		std::size_t j_l(0), j_r(0);
		std::size_t accu(lhs.nonzeros(i) + rhs.nonzeros(i));

		while (j_l < lhs.nonzeros(i) && j_r < rhs.nonzeros(i)) {
			if (lhs(i, j_l).index_ < rhs(i, j_r).index_)
				++j_l;
			else if (lhs(i, j_l).index_ > rhs(i, j_r).index_)
				++j_r;
			else {
				--accu; ++j_l; ++j_r;
			}
		}

		capacities[i] = accu;
	}

	// merge matrices 
	MatrixCRS tmp(lhs.rows(), lhs.columns(), capacities);
	for (std::size_t i = 0; i < lhs.rows(); ++i) {
		std::size_t j_l(0), j_r(0);

		while (j_l < lhs.nonzeros(i) && j_r < rhs.nonzeros(i)) {
			if (lhs(i, j_l).index_ < rhs(i, j_r).index_) {
				tmp.appendRowElement(i, lhs(i, j_l).index_, lhs(i, j_l).value_);
				++j_l;
			}
			else if (lhs(i, j_l).index_ > rhs(i, j_r).index_) {
				tmp.appendRowElement(i, rhs(i, j_r).index_, rhs(i, j_r).value_);
				++j_r;
			}
			else {
				tmp.appendRowElement(i, lhs(i, j_l).index_, lhs(i, j_l).value_ + rhs(i, j_r).value_);
				++j_l; ++j_r;
			}
		}

		while (j_l < lhs.nonzeros(i)) {
			tmp.appendRowElement(i, lhs(i, j_l).index_, lhs(i, j_l).value_);
			++j_l;
		}

		while (j_r < rhs.nonzeros(i)) {
			tmp.appendRowElement(i, rhs(i, j_r).index_, rhs(i, j_r).value_);
			++j_r;
		}
	}

	return tmp;
}


/**
 * \brief Subtraction operator for the subtraction of two matrices (\f$ A=B-C \f$).
 *
 * \param lhs The left-hand side matrix for the matrix subtraction.
 * \param rhs The right-hand-side matrix to be subtracted from the left-hand side matrix.
 * \return The difference of the two matrices.
 */
const MatrixCRS operator-(const MatrixCRS& lhs, const MatrixCRS& rhs) {
	assert((lhs.rows() == rhs.rows() && rhs.columns() == rhs.columns()) && "Matrix sizes do not match");

	// analyze matrices to predict the required row storage capacities
	std::vector<std::size_t> capacities(lhs.rows(), 0);
	for (std::size_t i = 0; i < lhs.rows(); ++i) {
		std::size_t j_l(0), j_r(0);
		std::size_t accu(lhs.nonzeros(i) + rhs.nonzeros(i));

		while (j_l < lhs.nonzeros(i) && j_r < rhs.nonzeros(i)) {
			if (lhs(i, j_l).index_ < rhs(i, j_r).index_)
				++j_l;
			else if (lhs(i, j_l).index_ > rhs(i, j_r).index_)
				++j_r;
			else {
				--accu; ++j_l; ++j_r;
			}
		}

		capacities[i] = accu;
	}

	// merge matrices 
	MatrixCRS tmp(lhs.rows(), lhs.columns(), capacities);
	for (std::size_t i = 0; i < lhs.rows(); ++i) {
		std::size_t j_l(0), j_r(0);

		while (j_l < lhs.nonzeros(i) && j_r < rhs.nonzeros(i)) {
			if (lhs(i, j_l).index_ < rhs(i, j_r).index_) {
				tmp.appendRowElement(i, lhs(i, j_l).index_, lhs(i, j_l).value_);
				++j_l;
			}
			else if (lhs(i, j_l).index_ > rhs(i, j_r).index_) {
				tmp.appendRowElement(i, rhs(i, j_r).index_, rhs(i, j_r).value_);
				++j_r;
			}
			else {
				tmp.appendRowElement(i, lhs(i, j_l).index_, lhs(i, j_l).value_ - rhs(i, j_r).value_);
				++j_l; ++j_r;
			}
		}

		while (j_l < lhs.nonzeros(i)) {
			tmp.appendRowElement(i, lhs(i, j_l).index_, lhs(i, j_l).value_);
			++j_l;
		}

		while (j_r < rhs.nonzeros(i)) {
			tmp.appendRowElement(i, rhs(i, j_r).index_, rhs(i, j_r).value_);
			++j_r;
		}
	}

	return tmp;
}


/**
 * \brief Multiplication operator for the multiplication of two matrices (\f$ A=B*C \f$).
 *
 * \param lhs The left-hand side matrix for the multiplication.
 * \param rhs The right-hand-side matrix for the multiplication.
 * \return The resulting matrix.
 */
const MatrixCRS operator*(const MatrixCRS& lhs, const MatrixCRS& rhs) {
	assert(lhs.columns() == rhs.rows() && "Matrix sizes do not match");

	std::vector<IndexValuePair> row(rhs.columns());
	MatrixCRS                   rhs_trans = rhs.getTranspose();
	std::size_t                 row_len;

	MatrixCRS tmp(lhs.rows(), rhs.columns());
	for (std::size_t i = 0; i < lhs.rows(); ++i) {
		row_len = 0;

		for (std::size_t j = 0; j < row.size(); ++j) {
			std::size_t k_l(0), k_r(0);
			double accu = 0.0;

			while (k_l < lhs.nonzeros(i) && k_r < rhs_trans.nonzeros(j)) {
				if (lhs(i, k_l).index_ < rhs_trans(j, k_r).index_)
					++k_l;
				else if (lhs(i, k_l).index_ > rhs_trans(j, k_r).index_)
					++k_r;
				else {
					accu += lhs(i, k_l).value_ * rhs_trans(j, k_r).value_;
					++k_l; ++k_r;
				}
			}

			if (accu != 0)
				row[row_len++] = IndexValuePair(accu, j);
		}

		// merge sparse vector into matrix row
		tmp.reserveRowElements(i, row_len);
		for (std::size_t j = 0; j < row_len; ++j)
			tmp.appendRowElement(i, row[j].index_, row[j].value_);
	}

	return tmp;
}

