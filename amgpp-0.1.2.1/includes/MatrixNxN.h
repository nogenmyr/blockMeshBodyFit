/**
 * \file includes/MatrixNxN.h
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

#ifndef MATRIXNXN_H_
#define MATRIXNXN_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "VectorN.h"
#include "MatrixCRS.h"


/**
 * \brief Efficient implementation of a \f$ N \times N \f$ matrix.
 *
 * The MatrixNxN class is the representation of a quadratic \f$ N \times N \f$ matrix with a
 * total of \f$ N^2 \f$ dynamically allocated elements. These elements can be directly accessed
 * with the 1D subscript operator or with the 2D function operator. The matrix is stored in a
 * row-wise fashion:
 *
 *	\f[\left (\begin{array}{*{5}{c}}
 *	0      & 1       & 2       & \cdots & N-1         \\
 *	N      & N+1     & N+2     & \cdots & 2 \cdot N-1 \\
 *	\vdots & \vdots  & \vdots  & \ddots & \vdots      \\
 *	N^2-N  & N^2-N+1 & N^2-N+2 & \cdots & N^2-1       \\
 *	\end{array}\right)\f]
 *
 * MatrixNxN can be used with any element type. The arithmetic operators for matrix/matrix,
 * matrix/vector and matrix/element operations with the same element type work for any element
 * type as long as the element type supports the arithmetic operation. Arithmetic operations
 * between matrices, vectors and elements of different element types are only supported for
 * all data types supported by the MathTrait class template (for details see the MathTrait
 * class description).
 *
 *	\code
 *	MatrixNxN< double > a, b, c;
 *	MatrixNxN< float  > d;
 *	MatrixNxN< std::complex<double> > e, f, g;
 *	MatrixNxN< std::complex<float>  > h;
 *	
 *	...         // Appropriate resizing
 *	
 *	c = a + b;  // OK: Same element type, supported
 *	c = a + d;  // OK: Different element types, supported by the MathTrait class template
 *	
 *	g = e + f;  // OK: Same element type, supported
 *	g = e + h;  // Error: Different element types, not supported by the MathTrait class template
 *	\endcode
 */
class MatrixNxN
{
public:
	explicit inline        MatrixNxN();
	explicit inline        MatrixNxN(std::size_t n);
	explicit inline        MatrixNxN(std::size_t n, double init);
	inline                 MatrixNxN(const MatrixNxN& m);
	template<typename Other, std::size_t N>
	inline                 MatrixNxN(const Other (&rhs)[N][N]);
	inline                 ~MatrixNxN();

	template<typename Other, std::size_t N>
	inline MatrixNxN&      operator=(const Other (&rhs)[N][N]);
	inline MatrixNxN&      operator=(double rhs);
	inline MatrixNxN&      operator=(const MatrixNxN& rhs);
	inline MatrixNxN&      operator=(const MatrixCRS& rhs);
	inline double&         operator[](std::size_t i);
	inline const double&   operator[](std::size_t i)                  const;
	inline double&         operator()(std::size_t i, std::size_t j);
	inline const double&   operator()(std::size_t i, std::size_t j)   const;
	inline const MatrixNxN operator-()                                const;
	inline MatrixNxN&      operator+=(const MatrixNxN& rhs);
	inline MatrixNxN&      operator-=(const MatrixNxN& rhs);
	template<typename Other>
	inline MatrixNxN&      operator*=(Other rhs);
	inline MatrixNxN&      operator*=(const MatrixNxN& rhs);
	inline const MatrixNxN operator+(const MatrixNxN& rhs);
	inline const MatrixNxN operator-(const MatrixNxN& rhs);
	inline const MatrixNxN operator*(const MatrixNxN& rhs);
	template<typename Other>
	inline const MatrixNxN operator*(Other scalar);
	inline const VectorN   operator*(const VectorN& vec);
	inline bool            operator==(const MatrixNxN& rhs);
	template<typename Other>
	inline bool            operator==(Other scalar);
	inline bool            operator!=(const MatrixNxN& rhs);
	template<typename Other>
	inline bool            operator!=(Other scalar);

	inline std::size_t     size()                                     const;
	inline std::size_t     rows()                                     const;
	inline std::size_t     columns()                                  const;
	inline void            reset();
	inline void            clear();
	inline void            resize(std::size_t n, bool preserve=false);
	inline void            extend(std::size_t n, bool preserve=false);
	inline const MatrixNxN getTranspose()                             const;
	inline bool            isSymmetric()                              const;
	template<typename Other>
	inline MatrixNxN&      scale(Other scalar);
	inline bool            equal(double scalar)                       const;
	inline bool            equal(const MatrixNxN& mat)                const;
	inline void            swap(MatrixNxN& m) throw();

	void read(const char* file);
	void write(const char* file, std::streamsize prec=6)              const;

private:
	std::size_t size_;        //!< The current size/dimension of the matrix.
	std::size_t capacity_;    //!< The maximum capacity of the matrix.
	double* v_;               //!< The dynamically allocated matrix elements.
	/**<
	 * Access to the matrix elements is gained via the subscript or
	 * function call operator. The order of the elements is
	 * 	\f[\left (\begin{array}{*{5}{c}}
	 * 	0      & 1       & 2       & \cdots & N-1         \\
	 * 	N      & N+1     & N+2     & \cdots & 2 \cdot N-1 \\
	 * 	\vdots & \vdots  & \vdots  & \ddots & \vdots      \\
	 * 	N^2-N  & N^2-N+1 & N^2-N+2 & \cdots & N^2-1       \\
	 * 	\end{array}\right)\f]
	 */

	template<typename Other>
	friend const MatrixNxN operator*(Other scalar, const MatrixNxN& lhs);
};


template<typename Other>
inline bool operator==(Other scalar, const MatrixNxN& mat);

template<typename Other>
inline bool operator!=(Other scalar, const MatrixNxN& mat);

inline std::ostream& operator<<(std::ostream& os, const MatrixNxN& m);

inline bool isnan(const MatrixNxN& v);

inline const MatrixNxN fabs(const MatrixNxN& v);

inline void swap(MatrixNxN& a, MatrixNxN& b) throw();

template<typename Other>
inline const MatrixNxN operator*(Other scalar, const MatrixNxN& lhs);


/**
 * \brief The default constructor for MatrixNxN.
 */
inline MatrixNxN::MatrixNxN() :
	size_(0), capacity_(0), v_(0)
{}


/**
 * \brief Constructor for a matrix of size \f$ n \times n \f$. No element initialization is performed!
 *
 * \param n The size of the matrix in both dimensions.
 *
 * \b Note: This constructor is only responsible to allocate the required dynamic memory. No
 *          element initialization is performed!
 */
inline MatrixNxN::MatrixNxN(std::size_t n) : size_(n), capacity_(n), v_(new double[n*n])
{}


/**
 * \brief Constructor for a homogenous initialization of all \f$ n \times n \f$ matrix elements.
 *
 * \param n The size of the matrix in both dimensions.
 * \param init The initial value of the matrix elements.
 *
 * All matrix elements are initialized with the specified value.
 */
inline MatrixNxN::MatrixNxN(std::size_t n, double init) : size_(n), capacity_(n), v_(new double[n*n])
{
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] = init;
}


/**
 * \brief The copy constructor for MatrixNxN.
 *
 * \param m Matrix to be copied.
 *
 * The copy constructor is explicitly defined due to the required dynamic memory management
 * and in order to enable/facilitate NRV optimization.
 */
inline MatrixNxN::MatrixNxN(const MatrixNxN& m) : size_(m.size_), capacity_(m.size_), v_(new double[size_*size_])
{
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] = m.v_[i];
}


/**
 * \brief Array initialization of all matrix elements.
 *
 * \param rhs \f$ N \times N \f$ dimensional array for the initialization.
 * \return Reference to the assigned matrix.
 *
 * This constructor offers the option to directly initialize the elements of the matrix:
 *
 *	\code
 *	const real init[3][3] = { { 1, 2, 3 },
 *	                          { 4, 5 },
 *	                          { 7, 8, 9 } };
 *	MatrixNxN<real> A = init;
 *	\endcode
 *
 * The matrix is sized accoring to the size of the array and initialized with the given values.
 * Missing values are initialized with zero (as e.g. the value 6 in the example).
 */
template<typename Other, std::size_t N>
inline MatrixNxN::MatrixNxN(const Other (&rhs)[N][N]) : size_(N), capacity_(N), v_(new double[N*N])
{
	for (std::size_t i=0; i<N; ++i)
		for (std::size_t j=0; j<N; ++j)
			v_[i*N+j] = rhs[i][j];
}


/**
 * \brief The destructor for MatrixNxN.
 */
inline MatrixNxN::~MatrixNxN()
{
	delete [] v_;
}


/**
 * \brief Array assignment to all matrix elements.
 *
 * \param rhs \f$ N \times N \f$ dimensional array for the asignment.
 * \return Reference to the assigned matrix.
 *
 * This assignment operator offers the option to directly set all elements of the matrix:
 *
 *	\code
 *	const real init[3][3] = { { 1, 2, 3 },
 *	                          { 4, 5 },
 *	                          { 7, 8, 9 } };
 *	MatrixNxN<real> A;
 *	A = init;
 *	\endcode
 *
 * The matrix is resized accoring to the size of the array and initialized with the given values.
 * Missing values are initialized with zero (as e.g. the value 6 in the example).
 */
template<typename Other, std::size_t N>
inline MatrixNxN& MatrixNxN::operator=(const Other (&rhs)[N][N])
{
	resize(N, false);

	for (std::size_t i=0; i<N; ++i)
		for (std::size_t j=0; j<N; ++j)
			v_[i*N+j] = rhs[i][j];

	return *this;
}


/**
 * \brief Homogenous assignment to all matrix elements.
 *
 * \param rhs Scalar value to be assigned to all matrix elements.
 * \return Reference to the assigned matrix.
 */
inline MatrixNxN& MatrixNxN::operator=(double rhs)
{
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] = rhs;
	return *this;
}




/**
 * \brief Copy assignment operator for MatrixNxN.
 *
 * \param rhs Matrix to be copied.
 * \return Reference to the assigned matrix.
 *
 * The matrix is resized according to the given \f$ N \times N \f$ matrix and initialized as a
 * copy of this matrix.
 */
inline MatrixNxN& MatrixNxN::operator=(const MatrixNxN& rhs)
{
	if (&rhs == this) return *this;

	resize(rhs.size_, false);

	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] = rhs.v_[i];

	return *this;
}


/**
 * \brief Assignment operator for sparse matrix instances.
 *
 * \param rhs Sparse matrix to be copied.
 * \return Reference to the assigned matrix.
 *
 * The matrix is resized according to the given sparse matrix and initialized as a
 * copy of this matrix.
 */
inline MatrixNxN& MatrixNxN::operator=(const MatrixCRS& rhs)
{
	assert(rhs.rows() == rhs.columns () && "Number of rows and columns do not match.");

	resize (rhs.rows(), false);

	(*this) = 0.0;
	for (std::size_t i = 0; i < rhs.rows(); ++i) {
		for (std::size_t j = rhs.begin(i); j < rhs.end(i); ++j) {
			(*this)(i, rhs[j].index_) = rhs[j].value_; 
		}
	}

	return *this;
}


/**
 * \brief 1D-access to the matrix elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..N^2-1]\f$.
 * \return Reference to the accessed value.
 *
 * In case assert() is active, this operator performs an index check.
 */
inline double& MatrixNxN::operator[](std::size_t index)
{
	assert(index < size_*size_ && "Invalid matrix access index");
	return v_[index];
}


/**
 * \brief 1D-access to the matrix elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..N^2-1]\f$.
 * \return Reference to the accessed value.
 *
 * In case assert() is active, this operator performs an index check.
 */
inline const double& MatrixNxN::operator[](std::size_t index) const
{
	assert(index < size_*size_ && "Invalid matrix access index");
	return v_[index];
}


/**
 * \brief 2D-access to the matrix elements.
 *
 * \param i Access index for the row. The index has to be in the range \f$[0..N-1]\f$.
 * \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
 * \return Reference to the accessed value.
 */
inline double& MatrixNxN::operator()(std::size_t i, std::size_t j)
{
	assert(i<size_ && j<size_ && "Invalid matrix access index");
	return v_[i*capacity_+j];
}


/**
 * \brief 2D-access to the matrix elements.
 *
 * \param i Access index for the row. The index has to be in the range \f$[0..N-1]\f$.
 * \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
 * \return Reference to the accessed value.
 */
inline const double& MatrixNxN::operator()(std::size_t i, std::size_t j) const
{
	assert(i<size_ && j<size_ && "Invalid matrix access index");
	return v_[i*capacity_+j];
}


/**
 * \brief Unary minus operator for the inversion of a matrix (\f$ A = -B \f$).
 *
 * \return The inverse of the matrix.
 */
inline const MatrixNxN MatrixNxN::operator-() const
{
	MatrixNxN tmp(size_);

	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		tmp[i] = -v_[i];
	return tmp;
}


/**
 * \brief Addition assignment operator for the addition of two matrices (\f$ A+=B \f$).
 *
 * \param rhs The right-hand side matrix to be added to the matrix.
 * \return Reference to the matrix.
 */
inline MatrixNxN& MatrixNxN::operator+=(const MatrixNxN& rhs)
{
	assert(rhs.size_ == size_ && "Matrix sizes do not match");

	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] += rhs.v_[i];

	return *this;
}


/**
 * \brief Subtraction assignment operator for the subtraction of two matrices (\f$ A-=B \f$).
 *
 * \param rhs The right-hand side matrix to be subtracted from the matrix.
 * \return Reference to the matrix.
 */
inline MatrixNxN& MatrixNxN::operator-=(const MatrixNxN& rhs)
{
	assert(rhs.size_ == size_ && "Matrix sizes do not match");

	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] -= rhs.v_[i];

	return *this;
}


/**
 * \brief Multiplication assignment operator for the multiplication between a matrix and
 * \brief a scalar value (\f$ A*=s \f$).
 *
 * \param rhs The right-hand side scalar value for the multiplication.
 * \return Reference to the matrix.
 */
template<typename Other>
inline MatrixNxN& MatrixNxN::operator*=(Other rhs)
{
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] *= rhs;

	return *this;
}


/**
 * \brief Multiplication assignment operator for the multiplication between two matrices
 * \brief (\f$ A*=B \f$).
 *
 * \param rhs The right-hand side matrix for the multiplication.
 * \return Reference to the matrix.
 */
inline MatrixNxN& MatrixNxN::operator*=(const MatrixNxN& rhs)
{
	assert(rhs.size_ == size_ && "Matrix sizes do not match");

	MatrixNxN tmp(size_, 0);

	for (std::size_t i=0; i<size_; ++i) {
		for (std::size_t j=0; j<size_; ++j) {
			for (std::size_t k=0; k<size_; ++k) {
				tmp.v_[i*size_+k] += v_[i*size_+j] * rhs.v_[j*size_+k];
			}
		}
	}

	return this->operator=(tmp);
}


/**
 * \brief Addition operator for the addition of two matrices (\f$ A=B+C \f$).
 *
 * \param lhs The left-hand side matrix for the matrix addition.
 * \param rhs The right-hand side matrix to be added to the left-hand side matrix.
 * \return The sum of the two matrices.
 */
inline const MatrixNxN MatrixNxN::operator+(const MatrixNxN& rhs)
{
	assert(size_ == rhs.size_ && "Matrix sizes do not match");

	MatrixNxN tmp(size_);
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		tmp.v_[i] = v_[i] + rhs.v_[i];

	return tmp;
}


/**
 * \brief Subtraction operator for the subtraction of two matrices (\f$ A=B-C \f$).
 *
 * \param lhs The left-hand side matrix for the matrix subtraction.
 * \param rhs The right-hand side matrix to be subtracted from the matrix.
 * \return The difference of the two matrices.
 */
inline const MatrixNxN MatrixNxN::operator-(const MatrixNxN& rhs)
{
	assert(size_ == rhs.size_ && "Matrix sizes do not match");

	MatrixNxN tmp(size_);
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		tmp.v_[i] = v_[i] - rhs.v_[i];

	return tmp;
}


/**
 * \brief Multiplication operator for the multiplication of two matrices (\f$ A=B*C \f$).
 *
 * \param lhs The left-hand side matirx for the multiplication.
 * \param rhs The right-hand side matrix for the multiplication.
 * \return The resulting matrix.
 */
inline const MatrixNxN MatrixNxN::operator*(const MatrixNxN& rhs)
{
	assert(size_ == rhs.size_ && "Matrix sizes do not match");

	MatrixNxN tmp(size_, 0);

	for (std::size_t i=0; i<size_; ++i) {
		for (std::size_t j=0; j<size_; ++j) {
			for (std::size_t k=0; k<size_; ++k) {
				tmp.v_[i*size_+k] += v_[i*size_+j] * rhs.v_[j*size_+k];
			}
		}
	}

	return tmp;
}


/**
 * \brief Multiplication operator for the multiplication of a matrix and a scalar value
 * \brief (\f$ A=B*s \f$).
 *
 * \param mat The left-hand side matrix for the multiplication.
 * \param scalar The right-hand side scalar value for the multiplication.
 * \return The scaled result matrix.
 */
template<typename Other>
inline const MatrixNxN MatrixNxN::operator*(Other scalar)
{
	MatrixNxN tmp(size_);

	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		tmp.v_[i] = scalar * v_[i];

	return tmp;
}


/**
 * \brief Multiplication operator for the multiplication of a matrix and a vector
 * \brief (\f$ \vec{a}=B*\vec{c} \f$).
 *
 * \param mat The left-hand side matrix for the multiplication.
 * \param vec The right-hand side vector for the multiplication.
 * \return The resulting vector.
 */
inline const VectorN MatrixNxN::operator*(const VectorN& vec)
{
	assert(size_ == vec.size() && "Matrix and vector sizes do not match");

	VectorN tmp(size_, 0);

	for (std::size_t i=0; i<size_; ++i) {
		for (std::size_t j=0; j<size_; ++j) {
			tmp[i] += v_[i*size_+j] * vec[j];
		}
	}

	return tmp;
}


/**
 * \brief Equality operator for the comparison of two matrices.
 *
 * \param lhs The left-hand side matrix for the comparison.
 * \param rhs The right-hand side matrix for the comparison.
 * \return \a true if the two matrices are equal, \a false if not.
 */
inline bool MatrixNxN::operator==(const MatrixNxN& rhs)
{
	return equal(rhs);
}


/**
 * \brief Equality operator for the comparison of a matrix and a scalar value.
 *
 * \param mat The left-hand side matrix for the comparison.
 * \param scalar The right-hand side scalar value for the comparison.
 * \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
 *
 * If all values of the matrix are equal to the scalar value, the equality test returns true,
 * otherwise false.
 */
template<typename Other>
inline bool MatrixNxN::operator==(Other scalar)
{
	return equal(scalar);
}


/**
 * \brief Inequality operator for the comparison of two matrices.
 *
 * \param lhs The left-hand side matrix for the comparison.
 * \param rhs The right-hand side matrix for the comparison.
 * \return \a true if the two matrices are not equal, \a false if they are equal.
 */
inline bool MatrixNxN::operator!=(const MatrixNxN& rhs)
{
	return !equal(rhs);
}


/**
 * \brief Inequality operator for the comparison of a matrix and a scalar value.
 *
 * \param mat The left-hand side matrix for the comparison.
 * \param scalar The right-hand side scalar value for the comparison.
 * \return \a true if at least one element of the matrix is different from the scalar, \a false if not.
 *
 * If one value of the matrix is inequal to the scalar value, the inequality test returns true,
 * otherwise false.
 */
template<typename Other>
inline bool MatrixNxN::operator!=(Other scalar)
{
	return !equal(scalar);
}


/**
 * \brief Returns the current size/dimension of the matrix.
 *
 * \return The size of the matrix.
 *
 * The function returns the current size N of the \f$ N^2 \f$ matrix.
 */
inline std::size_t MatrixNxN::size() const
{
	return size_;
}


/**
 * \brief Returns the current number of rows of the matrix.
 *
 * \return The number of rows of the matrix.
 */
inline std::size_t MatrixNxN::rows() const
{
	return size_;
}


/**
 * \brief Returns the current number of columns of the matrix.
 *
 * \return The number of columns of the matrix.
 */
inline std::size_t MatrixNxN::columns() const
{
	return size_;
}


/**
 * \brief Reset to the default initial values.
 *
 * \return void
 */
inline void MatrixNxN::reset()
{
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] = 0.0;
}


/**
 * \brief Clearing the \f$ N \times N \f$ matrix.
 *
 * \return void
 *
 * After the clear() function, the size of the matrix is 0.
 */
inline void MatrixNxN::clear()
{
	size_ = 0;
}


/**
 * \brief Changing the size of the matrix.
 *
 * \param n The new size of the matrix.
 * \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
 * \return void
 *
 * This function resizes the matrix using the given size to \f$ n \times n \f$. During this
 * operation, new dynamic memory may be allocated in case the capacity of the matrix is too
 * small. Therefore this function potentially changes all matrix elements. In order to preserve
 * the old matrix values, the \a preserve flag can be set to \a true. However, new matrix
 * elements are not initialized!\n
 * The following example illustrates the resize operation of a \f$ 2 \times 2 \f$ matrix to a
 * \f$ 3 \times 3 \f$ matrix. The new, uninitialized elements are marked with \a x:
 *
 *	\f[
 *	\left(\begin{array}{*{2}{c}}
 *	1 & 2 \\
 *	3 & 4 \\
 *	\end{array}\right)
 *	
 *	\Longrightarrow
 *	
 *	\left(\begin{array}{*{3}{c}}
 *	1 & 2 & x \\
 *	3 & 4 & x \\
 *	x & x & x \\
 *	\end{array}\right)
 *	\f]
 */
inline void MatrixNxN::resize(std::size_t n, bool preserve)
{
	if (n == size_) return;

	if (preserve)
	{
		double* v = new double[n*n];
		const std::size_t minsize(std::min(n, size_));

		for (std::size_t i=0; i<minsize; ++i)
			for (std::size_t j=0; j<minsize; ++j)
				v[i*n+j] = v_[i*size_+j];

		std::swap(v_, v);
		delete [] v;
		capacity_ = n;
	}
	else if (n > capacity_) {
		double* v = new double[n*n];
		std::swap(v_, v);
		delete [] v;
		capacity_ = n;
	}

	size_ = n;
}


/**
 * \brief Extending the size of the matrix.
 *
 * \param n Number of additional rows and columns.
 * \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
 * \return void
 *
 * This function increases the matrix size by \a n rows and \a n columns. During this operation,
 * new dynamic memory may be allocated in case the capacity of the matrix is too small. Therefore
 * this function potentially changes all matrix elements. In order to preserve the old matrix
 * values, the \a preserve flag can be set to \a true. However, new matrix elements are not
 * initialized!
 */
inline void MatrixNxN::extend(std::size_t n, bool preserve)
{
	resize(size_+n, preserve);
}


/**
 * \brief Calculation of the transpose of the matrix.
 *
 * \return The transpose of the matrix.
 */
inline const MatrixNxN MatrixNxN::getTranspose() const
{
	MatrixNxN tmp(size_);

	for (std::size_t i=0; i<size_; ++i) {
		for (std::size_t j=0; j<size_; ++j) {
			tmp.v_[j*size_+i] = v_[i*size_+j];
		}
	}

	return tmp;
}


/**
 * \brief Calculation of the transpose of the matrix.
 *
 * \return The transpose of the matrix.
 */
inline bool MatrixNxN::isSymmetric() const
{
	const std::size_t iend(size_-1);

	for (std::size_t i=0; i<iend; ++i) {
		for (std::size_t j=i+1; j<size_; ++j) {
			if (v_[i*size_+j] != v_[j*size_+i])
				return false;
		}
	}

	return true;
}


/**
 * \brief Scaling of the matrix by the scalar value \a scalar (\f$ A=B*s \f$).
 *
 * \param scalar The scalar value for the matrix scaling.
 * \return Reference to the matrix.
 */
template<typename Other>
inline MatrixNxN& MatrixNxN::scale(Other scalar)
{
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		v_[i] *= scalar;

	return *this;
}


/**
 * \brief Equality operator for the comparison of a matrix and a scalar value.
 *
 * \param scalar The scalar value for the comparison.
 * \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
 *
 * If all values of the matrix are equal to the scalar value, the equality test returns true,
 * otherwise false.
 */
inline bool MatrixNxN::equal(double scalar) const
{
	// In order to compare the vector and the scalar value, the data values of the lower-order
	// data type are converted to the higher-order data type within the equal function.
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		if (v_[i] != scalar) return false;
	return true;
}


/**
 * \brief Equality operator for the comparison of two vectors.
 *
 * \param mat The right-hand side matrix for the comparison.
 * \return \a true if the two vectors are equal, \a false if not.
 */
inline bool MatrixNxN::equal(const MatrixNxN& mat) const
{
	assert(mat.size_ == size_ && "Matrix sizes do not match");

	// In order to compare the vector and the scalar value, the data values of the lower-order
	// data type are converted to the higher-order data type within the equal function.
	const std::size_t sqrsize(size_*size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		if (v_[i] != mat.v_[i]) return false;
	return true;
}


/**
 * \brief Swapping the contents of two matrices.
 *
 * \param m The matrix to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
inline void MatrixNxN::swap(MatrixNxN& m) throw()
{
	std::swap(size_, m.size_);
	std::swap(capacity_, m.capacity_);
	std::swap(v_, m.v_);
}


/**
 * \brief Equality operator for the comparison of a scalar value and a matrix.
 *
 * \param scalar The left-hand side scalar value for the comparison.
 * \param mat The right-hand side matrix for the comparison.
 * \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
 *
 * If all values of the matrix are equal to the scalar value, the equality test returns true,
 * otherwise false.
 */
template<typename Other>
inline bool operator==(Other scalar, const MatrixNxN& mat)
{
	return mat.equal(scalar);
}


/**
 * \brief Inequality operator for the comparison of a scalar value and a matrix.
 *
 * \param scalar The left-hand side scalar value for the comparison.
 * \param mat The right-hand side matrix for the comparison.
 * \return \a true if at least one element of the matrix is different from the scalar, \a false if not.
 *
 * If one value of the matrix is inequal to the scalar value, the inequality test returns true,
 * otherwise false.
 */
template<typename Other>
inline bool operator!=(Other scalar, const MatrixNxN& mat)
{
	return !mat.equal(scalar);
}


/**
 * \brief Global output operator for NxN matrices.
 *
 * \param os Reference to the output stream.
 * \param m Reference to a constant matrix object.
 * \return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, const MatrixNxN& m)
{
	for (std::size_t i=0; i<m.rows(); ++i) {
		for (std::size_t j=0; j<m.columns(); ++j) {
			os << std::setw(14) << m(i,j);
		}
		os << "\n";
	}

	return os;
}


/**
 * \brief Checks the given matrix for not-a-number elements.
 *
 * \param m The matrix to be checked for not-a-number elements.
 * \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
 */
inline bool isnan(const MatrixNxN& m)
{
	const std::size_t sqrsize(m.size()*m.size());
	for (std::size_t i=0; i<sqrsize; ++i)
		if (isnan(m[i])) return true;
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
 * integral data type, the pe::abs(const MatrixN&) function can be used.
 */
inline const MatrixNxN fabs(const MatrixNxN& m)
{
	MatrixNxN tmp(m.size());

	const std::size_t sqrsize(m.size()*m.size());
	for (std::size_t i=0; i<sqrsize; ++i)
		tmp[i] = std::fabs(m[i]);

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
inline void swap(MatrixNxN& a, MatrixNxN& b) throw()
{
	a.swap(b);
}


/**
 * \brief Multiplication operator for the multiplication of a scalar value and a matrix
 * \brief (\f$ A=s*B \f$).
 *
 * \param scalar The left-hand side scalar value for the multiplication.
 * \param mat The right-hand side matrix for the multiplication.
 * \return The scaled result matrix.
 */
template<typename Other>
inline const MatrixNxN operator*(Other scalar, const MatrixNxN& mat)
{
	MatrixNxN tmp(mat.size_);

	const std::size_t sqrsize(mat.size_*mat.size_);
	for (std::size_t i=0; i<sqrsize; ++i)
		tmp.v_[i] = scalar * mat.v_[i];

	return tmp;
}

#endif
