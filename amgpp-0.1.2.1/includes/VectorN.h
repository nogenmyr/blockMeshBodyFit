/**
 * \file includes/VectorN.h
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

#ifndef VECTORN_H_
#define VECTORN_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <assert.h>


/**
 * \brief Efficient implementation of an arbitrary sized vector.
 * \ingroup vectorN
 *
 * The VectorN class is the representation of a vector with an arbitrary number N of dynamically
 * allocated elements. The elements can be accessed directly with the subscript operator. The
 * order of the elements is as following:

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right)\f]

 * VectorN can be used with any element type. The arithmetic operators for vector/vector and
 * vector/element operations with the same element type work for any element type as long as
 * the element type supports the arithmetic operation. Arithmetic operations between vectors
 * and elements of different element types are only supported for all data types supported by
 * the MathTrait class template (for details see the MathTrait class description).

   \code
   VectorN< double > a, b, c;
   VectorN< float  > d;
   VectorN< std::complex<double> > e, f, g;
   VectorN< std::complex<float>  > h;

   ...         // Appropriate resizing

   c = a + b;  // OK: Same element type, supported
   c = a + d;  // OK: Different element types, supported by the MathTrait class template

   g = e + f;  // OK: Same element type, supported
   g = e + h;  // Error: Different element types, not supported by the MathTrait class template
   \endcode
 */
class VectorN {
public:
	explicit inline      VectorN();
	explicit inline      VectorN(std::size_t n);
	explicit inline      VectorN(std::size_t n, double init);
	         inline      VectorN(const VectorN& v);
	template<std::size_t N>
	         inline      VectorN(const double (&rhs)[N]);
	         inline      ~VectorN();

	template<std::size_t N>
	inline VectorN&      operator= (const double (&rhs)[N]);
	inline VectorN&      operator= (double rhs);
	inline VectorN&      operator= (const VectorN& rhs);
	inline double&       operator[](std::size_t index);
	inline const double& operator[](std::size_t index)      const;
	inline const VectorN operator- ()                       const;
	inline VectorN&      operator+=(const VectorN& rhs);
	inline VectorN&      operator-=(const VectorN& rhs);
	inline VectorN&      operator*=(double rhs);
	inline VectorN&      operator/=(double rhs);

	inline std::size_t   size()                             const;
	inline std::size_t   capacity()                         const;
	inline void          reset();
	inline void          clear();
	       void          resize(std::size_t n, bool preserve=true);
	inline void          extend(std::size_t n, bool preserve=true);
	inline void          reserve(std::size_t n);
	inline double        length()                           const;
	inline double        sqrLength()                        const;
	inline VectorN&      normalize();
	inline const VectorN getNormalized()                    const;
	inline VectorN&      scale(double scalar);
	inline bool          equal(double scalar)               const;
	inline bool          equal(const VectorN& vec)          const;
	inline double        min()                              const;
	inline double        max()                              const;
	inline double        absmin()                           const;
	inline double        absmax()                           const;
	inline void          swap(VectorN& v) throw();

private:
	std::size_t size_;      //!< The current size/dimension of the vector.
	std::size_t capacity_;  //!< The maximum capacity of the vector.
	double* v_;             //!< The dynamically allocated vector elements.

	friend const VectorN operator+(const VectorN& lhs, const VectorN& rhs);
	friend const VectorN operator-(const VectorN& lhs, const VectorN& rhs);
	friend double        operator*(const VectorN& lhs, const VectorN& rhs);
	friend const VectorN operator*(const VectorN& vec, double scalar);
	friend const VectorN operator*(double scalar, const VectorN& vec);
	friend const VectorN operator/(const VectorN& vec, double scalar);
};


/**
 * \brief The default constructor for VectorN.
 */
inline VectorN::VectorN() : size_(0), capacity_(0), v_(0) {
}


/**
 * \brief Constructor for a vector of size \a n. No element initialization is performed!
 *
 * \param n The size of the vector.
 *
 * \b Note: This constructor is only responsible to allocate the required dynamic memory. No
 *          element initialization is performed!
 */
inline VectorN::VectorN(std::size_t n) : size_(n), capacity_(n), v_(new double[n]) {
}


/**
 * \brief Constructor for a homogenous initialization of all \a n vector elements.
 *
 * \param n The size of the vector.
 * \param init The initial value of the vector elements.
 *
 * All vector elements are initialized with the specified value.
 */
inline VectorN::VectorN(std::size_t n, double init) : size_(n), capacity_(n), v_(new double[n]) {
	for (std::size_t i=0; i<size_; ++i)
		v_[i] = init;
}


/**
 * \brief The copy constructor for VectorN.
 *
 * \param v Vector to be copied.
 *
 * The copy constructor is explicitly defined due to the required dynamic memory management
 * and in order to enable/facilitate NRV optimization.
 */
inline VectorN::VectorN(const VectorN& v) : size_(v.size_), capacity_(v.size_), v_(new double[size_]) {
	for (std::size_t i=0; i<size_; ++i)
		v_[i] = v.v_[i];
}


/**
 * \brief Array initialization of all vector elements.
 *
 * \param rhs N-dimensional array for the initialization.
 * \return Reference to the assigned vector.
 *
 * This assignment operator offers the option to directly initialize the elements of the vector:

 \code
 const real init[4] = { 1, 2, 3 };
 VectorN<real> v = init;
 \endcode

 * The vector is sized accoring to the size of the array and initialized with the given values.
 * Missing values are initialized with zero (as e.g. the fourth element in the example).
 */
template <std::size_t N>
inline VectorN::VectorN(const double (&rhs)[N]) : size_(N), capacity_(N), v_(new double[N]) {
	for (std::size_t i=0; i<N; ++i)
		v_[i] = rhs[i];
}


/**
 * \brief The destructor for VectorN.
 */
inline VectorN::~VectorN() {
	delete [] v_;
}


/**
 * \brief Array assignment to all vector elements.
 *
 * \param rhs N-dimensional array for the assignment.
 * \return Reference to the assigned vector.
 *
 * This assignment operator offers the option to directly set all elements of the vector:

 \code
 const real init[4] = { 1, 2, 3 };
 VectorN<real> v;
 v = init;
 \endcode

 * The vector is resized accoring to the size of the array and initialized with the given values.
 * Missing values are initialized with zero (as e.g. the fourth element in the example).
 */
template<std::size_t N>
inline VectorN& VectorN::operator=(const double (&rhs)[N]) {
	resize(N, false);

	for (std::size_t i=0; i<N; ++i)
		v_[i] = rhs[i];

	return *this;
}


/**
 * \brief Homogenous assignment to all vector elements.
 *
 * \param rhs Scalar value to be assigned to all vector elements.
 * \return Reference to the assigned vector.
 */
inline VectorN& VectorN::operator=(double rhs) {
	for (std::size_t i=0; i<size_; ++i)
		v_[i] = rhs;

	return *this;
}


/**
 * \brief Copy assignment operator for VectorN.
 *
 * \param rhs Vector to be copied.
 * \return Reference to the assigned vector.
 *
 * The vector is resized according to the given N-dimensional vector and initialized as a
 * copy of this vector.
 */
inline VectorN& VectorN::operator=(const VectorN& rhs)
{
	if (&rhs == this) return *this;

	resize(rhs.size_, false);

	for (std::size_t i=0; i<size_; ++i)
		v_[i] = rhs.v_[i];

	return *this;
}


/**
 * \brief Subscript operator for the direct access to the vector elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..N]\f$.
 * \return Reference to the accessed value.
 */
inline double& VectorN::operator[](std::size_t index) {
	assert(index < size_ && "Invalid vector access index");
	return v_[index];
}


/**
 * \brief Subscript operator for the direct access to the vector elements.
 *
 * \param index Access index. The index has to be in the range \f$[0..N]\f$.
 * \return Reference to the accessed value.
 */
inline const double& VectorN::operator[](std::size_t index) const {
	assert(index < size_ && "Invalid vector access index");
	return v_[index];
}


/**
 * \brief Unary minus operator for the inversion of a vector (\f$ \vec{a} = -\vec{b} \f$).
 *
 * \return The inverse of the vector.
 */
inline const VectorN VectorN::operator-() const {
	VectorN tmp(size_);
	for (std::size_t i=0; i<size_; ++i)
		tmp[i] = -v_[i];
	return tmp;
}


/**
 * \brief Addition assignment operator for the addition of two vectors (\f$ \vec{a}+=\vec{b} \f$).
 *
 * \param rhs The right-hand side vector to be added to the vector.
 * \return Reference to the vector.
 */
inline VectorN& VectorN::operator+=(const VectorN& rhs) {
	assert(rhs.size_ == size_ && "Vector sizes do not match");
	for (std::size_t i=0; i<size_; ++i)
		v_[i] += rhs.v_[i];
	return *this;
}


/**
 * \brief Subtraction assignment operator for the subtraction of two vectors
 * \brief (\f$ \vec{a}-=\vec{b} \f$).
 *
 * \param rhs The right-hand side vector to be subtracted from the vector.
 * \return Reference to the vector.
 */
inline VectorN& VectorN::operator-=(const VectorN& rhs) {
	assert(rhs.size_ == size_ && "Vector sizes do not match");
	for (std::size_t i=0; i<size_; ++i)
		v_[i] -= rhs.v_[i];
	return *this;
}


/**
 * \brief Multiplication assignment operator for the multiplication between a vector and
 * \brief a scalar value (\f$ \vec{a}*=s \f$).
 *
 * \param rhs The right-hand side scalar value for the multiplication.
 * \return Reference to the vector.
 */
inline VectorN& VectorN::operator*=(double rhs) {
	for (std::size_t i=0; i<size_; ++i)
		v_[i] *= rhs;
	return *this;
}


/**
 * \brief Division assignment operator for the division of a vector by a scalar value
 * \brief (\f$ \vec{a}/=s \f$).
 *
 * \param rhs The right-hand side scalar value for the division.
 * \return Reference to the vector.
 *
 * \b Note: A division by zero is only checked by an user assert.
 */
inline VectorN& VectorN::operator/=(double rhs) {
	assert(rhs != 0.0 && "Division by zero detected");

	const double tmp = 1.0 / rhs;
	for (std::size_t i=0; i<size_; ++i)
		v_[i] *= tmp;
	
	return *this;
}


/**
 * \brief Returns the current size/dimension of the vector.
 *
 * \return The size of the vector.
 */
inline std::size_t VectorN::size() const {
	return size_;
}


/**
 * \brief Returns the maximum capacity of the vector.
 *
 * \return The capacity of the vector.
 */
inline std::size_t VectorN::capacity() const {
	return capacity_;
}


/**
 * \brief Reset to the default initial values.
 *
 * \return void
 */
inline void VectorN::reset() {
	for (std::size_t i=0; i<size_; ++i)
		v_[i] = 0.0;
}


/**
 * \brief Clearing the vector.
 *
 * \return void
 *
 * After the clear() function, the size of the vector is 0.
 */
inline void VectorN::clear() {
	size_ = 0;
}


/**
 * \brief Extending the size of the vector.
 *
 * \param n Number of additional vector elements.
 * \param preserve \a true if the old values of the vector should be preserved, \a false if not.
 * \return void
 *
 * This function increases the vector size by \a n elements. During this operation, new dynamic
 * memory may be allocated in case the capacity of the vector is too small. Therefore this
 * function potentially changes all vector elements. In order to preserve the old vector values,
 * the \a preserve flag can be set to \a true. However, new vector elements are not initialized!
 */
inline void VectorN::extend(std::size_t n, bool preserve) {
	resize(size_+n, preserve);
}


/**
 * \brief Setting the minimum capacity of the vector.
 *
 * \param n The new minimum capacity of the vector.
 * \return void
 *
 * This function increases the capacity of the vector to at least \a n elements. The current
 * values of the vector elements are preserved.
 */
inline void VectorN::reserve(std::size_t n) {
	if (n > capacity_) {
		// Allocating a new array
		double* tmp = new double[n];

		// Replacing the old array
		std::copy(v_, v_+size_, tmp);
		std::swap(tmp, v_);
		capacity_ = n;
		delete [] tmp;
	}
}


/**
 * \brief Calculation of the vector length \f$|\vec{a}|\f$.
 *
 * \return The length of the vector.
 *
 * The return type of the length function depends on the actual type of the vector instance:
 *
 * <table border="0" cellspacing="0" cellpadding="1">
 *    <tr>
 *       <td width="250px"> \b double </td>
 *       <td width="100px"> \b Length </td>
 *    </tr>
 *    <tr>
 *       <td>float</td>
 *       <td>float</td>
 *    </tr>
 *    <tr>
 *       <td>integral data types and double</td>
 *       <td>double</td>
 *    </tr>
 *    <tr>
 *       <td>long double</td>
 *       <td>long double</td>
 *    </tr>
 * </table>
 */
inline double VectorN::length() const
{
	double sum = 0.0;
	for (std::size_t i=0; i<size_; ++i)
		sum += v_[i] * v_[i];

	return std::sqrt(sum);
}




/**
 * \brief Calculation of the vector square length \f$|\vec{a}|^2\f$.
 *
 * \return The square length of the vector.
 */
inline double VectorN::sqrLength() const {
	double sum = 0.0;
	for (std::size_t i=0; i<size_; ++i)
		sum += v_[i] * v_[i];

	return sum;
}


/**
 * \brief Normalization of the vector (\f$|\vec{a}|=1\f$).
 *
 * \return Reference to the vector.
 *
 * Normalization of the vector to a length of 1. This operation is only defined for floating
 * point vectors. The attempt to use this function for an integral vector results in a compile
 * time error.
 */
inline VectorN& VectorN::normalize()
{
	const double len(length());

	if (len == 0.0)
		return *this;

	const double ilen(1.0 / len);

	for (std::size_t i=0; i<size_; ++i)
		v_[i] *= ilen;

	return *this;
}


/**
 * \brief Calculation of the normalized vector (\f$|\vec{a}|=1\f$).
 *
 * \return The normalized vector.
 *
 * The function returns the normalized vector. This operation is only defined for floating
 * point vectors. The attempt to use this function for an integral vector results in a compile
 * time error.
 */
inline const VectorN VectorN::getNormalized() const {
	const double len(length());

	if (len == 0.0)
		return *this;

	const double ilen(1.0 / len);
	VectorN tmp(size_);

	for (std::size_t i=0; i<size_; ++i)
		tmp[i] = v_[i] * ilen;

	return tmp;
}


/**
 * \brief Scaling of the vector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
 *
 * \param scalar The scalar value for the vector scaling.
 * \return Reference to the vector.
 */
inline VectorN& VectorN::scale(double scalar) {
	for (std::size_t i=0; i<size_; ++i)
		v_[i] *= scalar;

	return *this;
}


/**
 * \brief Equality operator for the comparison of a vector and a scalar value.
 *
 * \param scalar The scalar value for the comparison.
 * \return \a true if all elements of the vector are equal to the scalar, \a false if not.
 *
 * If all values of the vector are equal to the scalar value, the equality test returns true,
 * otherwise false.
 */
inline bool VectorN::equal(double scalar) const {
	for (std::size_t i=0; i<size_; ++i)
		if (v_[i] != scalar) return false;

	return true;
}


/**
 * \brief Equality operator for the comparison of two vectors.
 *
 * \param vec The right-hand side vector for the comparison.
 * \return \a true if the two vectors are equal, \a false if not.
 */
inline bool VectorN::equal(const VectorN& vec) const {
	assert(vec.size_ == size_ && "Vector sizes do not match");

	for (std::size_t i=0; i<size_; ++i)
		if (v_[i] != vec.v_[i]) return false;

	return true;
}


/**
 * \brief Returns the smallest element of the vector.
 *
 * \return The smallest vector element.
 *
 * In case the vector currently has a size of 0, the returned value is 0.
 */
inline double VectorN::min() const {
	if (size_ == 0) return 0.0;

	double minimum(v_[0]);
	for (std::size_t i=1; i<size_; ++i)
		minimum = std::min(minimum, v_[i]);

	return minimum;
}


/**
 * \brief Returns the largest element of the vector.
 *
 * \return The largest vector element.
 *
 * In case the vector currently has a size of 0, the returned value is 0.
 */
inline double VectorN::max() const {
	if (size_ == 0) return 0.0;

	double maximum(v_[0]);
	for (std::size_t i=1; i<size_; ++i)
		maximum = std::max(maximum, v_[i]);

	return maximum;
}


/**
 * \brief Returns the minimum absolute value of all vector components.
 *
 * \return The minimum absolute value of all vector components.
 *
 * In case the vector currently has a size of 0, the returned value is 0.
 */
inline double VectorN::absmin() const {
	if (size_ == 0) return 0.0;

	double minimum(std::abs(v_[0]));
	for (std::size_t i=1; i<size_; ++i)
		minimum = std::min(minimum, std::abs(v_[i]));

	return minimum;
}


/**\brief Returns the maximum absolute value of all vector components.
 *
 * \return The maximum absolute value of all vector components.
 *
 * In case the vector currently has a size of 0, the returned value is 0.
 */
inline double VectorN::absmax() const {
	if (size_ == 0) return 0.0;

	double maximum(std::abs(v_[0]));
	for (std::size_t i=1; i<size_; ++i)
		maximum = std::max(maximum, std::abs(v_[i]));

	return maximum;
}




/**
 * \brief Swapping the contents of two vectors.
 *
 * \param v The vector to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
inline void VectorN::swap(VectorN& v) throw() {
	std::swap(size_, v.size_);
	std::swap(capacity_, v.capacity_);
	std::swap(v_, v.v_);
}


/**
 * \brief Equality operator for the comparison of two vectors.
 *
 * \param lhs The left-hand side vector for the comparison.
 * \param rhs The right-hand side vector for the comparison.
 * \return \a true if the two vectors are equal, \a false if not.
 */
inline bool operator==(const VectorN& lhs, const VectorN& rhs) {
	return lhs.equal(rhs);
}


/**
 * \brief Equality operator for the comparison of a vector and a scalar value.
 *
 * \param vec The left-hand side vector for the comparison.
 * \param scalar The right-hand side scalar value for the comparison.
 * \return \a true if all elements of the vector are equal to the scalar, \a false if not.
 *
 * If all values of the vector are equal to the scalar value, the equality test returns true,
 * otherwise false.
 */
inline bool operator==(const VectorN& vec, double scalar) {
	return vec.equal(scalar);
}


/**
 * \brief Equality operator for the comparison of a scalar value and a vector.
 * \ingroup vectorN
 *
 * \param scalar The left-hand side scalar value for the comparison.
 * \param vec The right-hand side vector for the comparison.
 * \return \a true if all elements of the vector are equal to the scalar, \a false if not.
 *
 * If all values of the vector are equal to the scalar value, the equality test returns true,
 * otherwise false.
 */
inline bool operator==(double scalar, const VectorN& vec) {
	return vec.equal(scalar);
}


/**
 * \brief Inequality operator for the comparison of two vectors.
 *
 * \param lhs The left-hand side vector for the comparison.
 * \param rhs The right-hand side vector for the comparison.
 * \return \a true if the two vectors are not equal, \a false if they are equal.
 */
inline bool operator!=(const VectorN& lhs, const VectorN& rhs) {
	return !lhs.equal(rhs);
}


/**
 * \brief Inequality operator for the comparison of a vector and a scalar value.
 *
 * \param vec The left-hand side vector for the comparison.
 * \param scalar The right-hand side scalar value for the comparison.
 * \return \a true if at least one element of the vector is different from the scalar, \a false if not.
 *
 * If one value of the vector is inequal to the scalar value, the inequality test returns true,
 * otherwise false.
 */
inline bool operator!=(const VectorN& vec, double scalar) {
	return !vec.equal(scalar);
}


/**
 * \brief Inequality operator for the comparison of a scalar value and a vector.
 *
 * \param scalar The left-hand side scalar value for the comparison.
 * \param vec The right-hand side vector for the comparison.
 * \return \a true if at least one element of the vector is different from the scalar, \a false if not.
 *
 * If one value of the vector is inequal to the scalar value, the inequality test returns true,
 * otherwise false.
 */
inline bool operator!=(double scalar, const VectorN& vec) {
	return !vec.equal(scalar);
}


/**
 * \brief Global output operator for arbitrary sized vectors.
 *
 * The output starts with the size of the vector followed by
 * each entry on a line of its own.
 *
 * \param os Reference to the output stream.
 * \param v Reference to a constant vector object.
 * \return Reference to the output stream.
 */
inline std::ostream& operator<<(std::ostream& os, const VectorN& v) {
	os << v.size() << '\n';
	for (std::size_t i=0; i<v.size(); ++i)
		os << v[i] << '\n';
	return os;
}


/**
 * \brief Global input operator for arbitrary sized vectors.
 *
 * The input format begins with the size of the vector followed by
 * the entries of the vector. The numbers are seperated by white space.
 *
 * \param is Reference to the input stream.
 * \param v Reference to a vector object.
 * \return Reference to the input stream.
 */
inline std::istream& operator>>(std::istream& is, VectorN& v) {
	std::size_t n;
	is >> n;

	v.resize(n, false);

	for (std::size_t i=0; i<n; ++i)
		is >> v[i];

	return is;
}


/**
 * \brief Checks the given vector for not-a-number elements.
 *
 * \param v The vector to be checked for not-a-number elements.
 * \return \a true if at least one element of the vector is not-a-number, \a false otherwise.
 */
inline bool isnan(const VectorN& v) {
	for (std::size_t i=0; i<v.size(); ++i) {
		if (isnan(v[i])) return true;
	}

	return false;
}


/**
 * \brief Returns a vector containing the absolute values of each single element of \a v.
 *
 * \param v The floating point input vector.
 * \return The absolute value of each single element of \a v.
 *
 * The \a fabs function calculates the absolute value of each element of the input vector \a v.
 * This function can only be applied to floating point vectors. For vectors of integral data
 * type, the pe::abs(const VectorN&) function can be used.
 */
inline const VectorN fabs(const VectorN& v) {
	VectorN tmp(v.size());
	for (std::size_t i=0; i<v.size(); ++i)
		tmp[i] = std::fabs(v[i]);

	return tmp;
}


/**
 * \brief Returns the minimum component of the vector.
 *
 * \return The minimum component of the vector.
 */
inline double min(const VectorN& v) {
	return v.min();
}


/**
 * \brief Returns the maximum component of the vector.
 *
 * \return The maximum component of the vector.
 */
inline double max(const VectorN& v) {
	return v.max();
}


/**
 * \brief Returns the minimum absolute value of all vector components.
 *
 * \return The minimum absolute value of all vector components.
 */
inline double absmin(const VectorN& v) {
	return v.absmin();
}


/**
 * \brief Returns the maximum absolute value of all vector components.
 *
 * \return The maximum absolute value of all vector components.
 */
inline double absmax(const VectorN& v) {
	return v.absmax();
}


/**
 * \brief Swapping the contents of two vectors.
 *
 * \param a The first vector to be swapped.
 * \param b The second vector to be swapped.
 * \return void
 * \exception no-throw guarantee.
 */
inline void swap(VectorN& a, VectorN& b) throw() {
	a.swap(b);
}



/**
 * \brief Addition operator for the addition of two vectors (\f$ \vec{a}=\vec{b}+\vec{c} \f$).
 *
 * \param lhs The left-hand side vector for the vector addition.
 * \param rhs The right-hand side vector for the vector addition.
 * \return The sum of the two vectors.
 *
 * The operator returns a vector of the higher-order data type of the two involved vector data
 * types \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data type
 * as long as the data type has an addition operator. In case \a T1 and \a T2 differ, the
 * operator only works for data types supported by the MathTrait class template.
 */
inline const VectorN operator+(const VectorN& lhs, const VectorN& rhs) {
	assert(lhs.size_ == rhs.size_ && "Vector sizes do not match");

	VectorN tmp(lhs.size_);
	for (std::size_t i=0; i<lhs.size_; ++i)
		tmp[i] = lhs.v_[i] + rhs.v_[i];

	return tmp;
}


/**
 * \brief Subtraction operator for the subtraction of two vectors (\f$ \vec{a}=\vec{b}-\vec{c} \f$).
 *
 * \param lhs The left-hand side vector for the vector subtraction.
 * \param rhs The right-hand side vector to be subtracted from the vector.
 * \return The difference of the two vectors.
 *
 * The operator returns a vector of the higher-order data type of the two involved vector data
 * types \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data type
 * as long as the data type has a subtraction operator. In case \a T1 and \a T2 differ, the
 * operator only works for data types supported by the MathTrait class template.
 */
inline const VectorN operator-(const VectorN& lhs, const VectorN& rhs) {
	assert(lhs.size_ == rhs.size_ && "Vector sizes do not match");

	VectorN tmp(lhs.size_);
	for (std::size_t i=0; i<lhs.size_; ++i)
		tmp[i] = lhs.v_[i] - rhs.v_[i];

	return tmp;
}


/**
 * \brief Multiplication operator for the scalar product (inner product) of two vectors
 * \brief (\f$ s=\vec{a}*\vec{b} \f$).
 *
 * \param lhs The left-hand side vector for the inner product.
 * \param rhs The right-hand side vector for the inner product.
 * \return The scalar product.
 *
 * The operator returns a scalar value of the higher-order data type of the two involved vector
 * data types \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data
 * type as long as the data type has a multiplication operator. In case \a T1 and \a T2 differ,
 * the operator only works for data types supported by the MathTrait class template.
 */
inline double operator*(const VectorN& lhs, const VectorN& rhs) {
	assert(lhs.size_ == rhs.size_ && "Vector sizes do not match");

	double sp(0);
	for (std::size_t i=0; i<lhs.size_; ++i)
		sp += lhs.v_[i] * rhs.v_[i];

	return sp;
}


/**
 * \brief Multiplication operator for the multiplication of a vector and a scalar value
 * \brief (\f$ \vec{a}=\vec{b}*s \f$).
 *
 * \param vec The left-hand side vector for the multiplication.
 * \param scalar The right-hand side scalar value for the multiplication.
 * \return The scaled result vector.
 *
 * The operator returns a vector of the higher-order data type of the two involved vector data
 * types \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data type
 * as long as the data type has a multiplication operator. In case \a T1 and \a T2 differ, the
 * operator only works for data types supported by the MathTrait class template.
 */
inline const VectorN operator*(const VectorN& vec, double scalar) {
	VectorN tmp(vec.size_);

	for (std::size_t i=0; i<vec.size_; ++i)
		tmp[i] = scalar * vec.v_[i];

	return tmp;
}


/**
 * \brief Multiplication operator for the multiplication of a scalar value and a vector
 * \brief (\f$ \vec{a}=s*\vec{b} \f$).
 *
 * \param scalar The left-hand side scalar value for the multiplication.
 * \param vec The right-hand side vector for the multiplication.
 * \return The scaled result vector.
 *
 * The operator returns a vector of the higher-order data type of the two involved vector data
 * types \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data type
 * as long as the data type has a multiplication operator. In case \a T1 and \a T2 differ, the
 * operator only works for data types supported by the MathTrait class template.
 */
inline const VectorN operator*(double scalar, const VectorN& vec)
{
	VectorN tmp(vec.size_);

	for (std::size_t i=0; i<vec.size_; ++i)
		tmp[i] = scalar * vec.v_[i];

	return tmp;
}


/**
 * \brief Division operator for the divison of a vector by a scalar value
 * \brief (\f$ \vec{a}=\vec{b}/s \f$).
 *
 * \param vec The left-hand side vector for the division.
 * \param scalar The right-hand side scalar value for the division.
 * \return The scaled result vector.
 *
 * The operator returns a vector of the higher-order data type of the two involved vector data
 * types \a T1 and \a T2. In case \a T1 and \a T2 match, the operator works for any data type
 * as long as the data type has a division operator. In case \a T1 and \a T2 differ, the
 * operator only works for data types supported by the MathTrait class template.
 *
 * \b Note: A division by zero is only checked by an user assert.
 */
inline const VectorN operator/(const VectorN& vec, double scalar) {
	assert(scalar != 0.0 && "Division by zero detected");

	VectorN tmp(vec.size_);

	const double idiv(1.0 / scalar);
	for (std::size_t i=0; i<vec.size_; ++i)
		tmp[i] = vec.v_[i] * idiv;

	return tmp;
}

#endif
