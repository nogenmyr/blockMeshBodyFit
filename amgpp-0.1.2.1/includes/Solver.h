/**
 * \file includes/Solver.h
 * \brief Header for the base class for all solver classes.
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

#ifndef SOLVER_H_
#define SOLVER_H_

#include <limits>
#include <algorithm>


/**
 * \brief Base class for all solver classes.
 */
class Solver {
public:
   explicit Solver();

   inline std::size_t getMaxIterations()  const;
   inline std::size_t getLastIterations() const;
   inline double      getLastPrecision()  const;
   inline double      getThreshold()      const;
   inline void        setMaxIterations(std::size_t maxIterations);
   inline void        setThreshold(double threshold);

protected:
   std::size_t maxIterations_;
   std::size_t lastIterations_;
   double lastPrecision_;
   double threshold_;
};


/**
 * \brief The default constructor.
 */
inline Solver::Solver() :
   maxIterations_(10000),
   lastIterations_(0),
   lastPrecision_(std::numeric_limits<double>::infinity()),
   threshold_(1.0e-7)
{
}


/**
 * \brief Returns the maximum number of iterations the solver may spend solving the problem.
 *
 * \return The maximum number of iterations spent in the solver.
 */
inline std::size_t Solver::getMaxIterations() const {
   return maxIterations_;
}


/**
 * \brief Returns the number of iterations spent in the last solution process.
 *
 * \return The number of iterations spent in the last solution process.
 */
inline std::size_t Solver::getLastIterations() const {
   return lastIterations_;
}


/**
 * \brief Returns the precision of the solution after the solution process.
 *
 * \return The precision of the solution after the solution process.
 * 
 * The solver is not enforced to compute the precision after the solution. Instead it can just
 * report infinity as the last precision.
 */
inline double Solver::getLastPrecision() const {
   return lastPrecision_;
}


/**
 * \brief Returns the threshold which classifies a solution as good enough.
 *
 * \return The threshold for the solution quality.
 */
inline double Solver::getThreshold() const {
   return threshold_;
}


/**
 * \brief Sets the maximum number of iterations the solver may spend solving the problem.
 *
 * \param maxIterations The maximum number of iterations spent in the solver.
 */
inline void Solver::setMaxIterations(std::size_t maxIterations) {
   maxIterations_ = maxIterations;
}


/**
 * \brief Sets the threshold which classifies a solution as good enough.
 *
 * \param threshold The threshold for the solution quality.
 */
inline void Solver::setThreshold(double threshold) {
   threshold_ = threshold;
}

#endif
