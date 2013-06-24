/**
 * \file includes/AmgSolver.h
 * \brief An algebraic multigrid black box solver.
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

#ifndef AMGSOLVER_H_
#define AMGSOLVER_H_

#include <vector>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include "Lse.h"
#include "Log.h"
#include "Solver.h"


/**
 * \brief An algebraic multigrid black box solver.
 *
 * TODO
 */
template<typename DirectSolver>
class AmgSolver : public Solver {
public:
	enum CoarseningStrategy {
		coarsening_amg1r5, coarsening_direct, coarsening_standard
	};

private:
	static const std::size_t undecided_point;
	static const std::size_t fine_point;
	static const std::size_t coarse_point;

public:
	explicit AmgSolver();

	bool solve(Lse& system);

	DirectSolver& getDirectSolver();
	const DirectSolver& getDirectSolver() const;

	void setCoarseningStrategy(CoarseningStrategy strategy);
	void setStartupExecution(bool fmg);
	void setSmoothingParameters(std::size_t nu1, std::size_t nu2);
	void setCoarseningParameters(double thetaPos, double thetaNeg);

private:
	void splitPoints(std::size_t level, double thetaPos, double thetaNeg);
	bool coarsen(std::size_t level, CoarseningStrategy strategy = coarsening_amg1r5);
	bool coarsenByAmg1r5(std::size_t level);
	bool coarsenByDirectInterpolation(std::size_t level);
	bool coarsenByStandardInterpolation(std::size_t level);
	bool startup();
	bool cycle(std::size_t level);
	bool solveDirect(std::size_t level);
	void smooth(std::size_t level);
	void initSmoother(std::size_t level);
	void computeGalerkinProduct(const MatrixCRS& P, const MatrixCRS& A, MatrixCRS& R, MatrixCRS& Ac);
	void outputSquareGrid(std::ostream& out, std::size_t level) const;
	void outputCoarseningGraph(std::ostream& out, const std::string& name, std::size_t level) const;
	std::string formatString(const std::string& format, std::size_t n);

	struct Grid {
		Lse system_;

		~Grid();

		std::vector<std::size_t> diag_index_;
		VectorN diag_inv_;

		VectorN                  tmp_;
		VectorN                  r_;
		std::vector<std::size_t> children_;
		MatrixCRS                P_;         //!< Prolongation operator
		MatrixCRS                R_;         //!< Restriction operator

		std::vector< std::vector<std::size_t> > S_;
		std::vector< std::vector<std::size_t> > St_;
	};

	std::vector<Grid> grids_;        //!< The full hierarchy of grids.
	std::size_t nu0_;                //!< The number of recursive descents in the full multigrid startup.
	std::size_t nu1_;                //!< The number of presmoothing sweeps.
	std::size_t nu2_;                //!< The number of postsmoothing sweeps.
	std::size_t mu_;                 //!< The number of recursive descents in the multigrid cycle. 
	bool fmg_;                       //!< Flag that indicates if a full multigrid startup is to be done.
	double thetaPos_;                  //!< Threshold for positive off-diagonal entries determining strong dependency.
	double thetaNeg_;                  //!< Threshold for positive off-diagonal entries determining strong dependency.
	std::size_t directSolveLimit_;   //!< The system size limit up to which the direct solver should be employed. 
	DirectSolver directSolver_;      //!< Solver instance to solve the coarsest grid.
	CoarseningStrategy strategy_;    //!< Coarsening strategy to be applied.

	std::vector<std::size_t> tmp_;
};

template<typename DirectSolver>
const std::size_t AmgSolver<DirectSolver>::undecided_point(std::numeric_limits<std::size_t>::max() - 2);

template<typename DirectSolver>
const std::size_t AmgSolver<DirectSolver>::coarse_point(std::numeric_limits<std::size_t>::max() - 1);

template<typename DirectSolver>
const std::size_t AmgSolver<DirectSolver>::fine_point(std::numeric_limits<std::size_t>::max());


/**
 * \brief The default constructor for the Amg class.
 */
template<typename DirectSolver>
AmgSolver<DirectSolver>::AmgSolver() :
	nu0_  (1),
	nu1_  (1),
	nu2_  (1),
	mu_   (1),
	fmg_  (true),
	thetaPos_(0.2),
	thetaNeg_(0.2),
	directSolveLimit_(128),
	directSolver_(),
	strategy_(coarsening_amg1r5) {
}


/**
 * \brief The destructor for the grid structure.
 *
 * This is essentially the default destructor and is only necessary due to compiler warnings
 * issued because of inlining failures.
 */
template<typename DirectSolver>
AmgSolver<DirectSolver>::Grid::~Grid() {
}


/**
 * \brief Solves the provided system using the black-box AMG solver.
 *
 * \param system The system to solve.
 * \return void
 *
 * The solution process first sets up the grid hierarchy, initializes the smoothers and
 * optionally starts with a full-multigrid startup. Afterwards the process proceeds in multigrid
 * cycles until the residual is below the threshold.
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::solve(Lse& system) {
	std::size_t n(system.size());
	bool converged(false);

	// setup grid hierarchy
	grids_.clear();
	grids_.push_back(Grid());
	grids_[0].system_ = system;
	grids_[0].r_.resize(n, false);

	//directSolveLimit_ = std::numeric_limits<std::size_t>::max();   //TODO
	//fmg_ = false;
	while (grids_.back().system_.size() > directSolveLimit_) {
		if (!coarsen(grids_.size() - 1, strategy_)) {
			if (loglevel_ <= WARNING) {
				std::cout << "WARNING: Failed to coarsen below direct solver limit." << std::endl;
			}
			break;
		}

		//std::ostringstream name;
		//name << "G" << grids_.size() - 2;
		//outputCoarseningGraph(std::cout, name.str(), grids_.size() - 2);
	}

	if (loglevel_ <= INFO)
		std::cout << "Created " << grids_.size() << " grid(s)." << std::endl;

	// prepare smoother on each grid
	for (std::size_t i = 0; i < grids_.size(); ++i)
		initSmoother(i);

	// check initial solution
	// [nothing to be done with LSEs]

	// compute initial residual
	lastPrecision_ = grids_[0].system_.residual();
	if (lastPrecision_ < threshold_)
		converged = true;

	// full multigrid startup
	if (!converged && fmg_)
		converged = startup();

	// multigrid cycles
	std::size_t it(0);
	for (; !converged && it<maxIterations_; ++it) {
		converged = cycle(0);
		//for (std::size_t i = 0; i < nu1_ + nu2_; ++i)
		//	smooth(0);	//TODO

		lastPrecision_ = grids_[0].system_.residual();
		if (lastPrecision_ < threshold_)
			converged = true;
	}

	if (converged) {
		if (loglevel_ <= INFO)
			std::cout << "Solved the LSE in " << it << " multigrid cycles." << std::endl;
	}
	else if (loglevel_ <= WARNING)
		std::cout << "WARNING: Did not solve the LSE within accuracy. (" << lastPrecision_ << ")" << std::endl;

	lastIterations_ = it;
	system.x_ = grids_[0].system_.x_;

	return converged;
}


/**
 * \brief Specifies the coarsening strategy used in the solution process.
 *
 * \param strategy The coarsening strategy used in the solution process. 
 * \return void
 *
 * There are several possibilities to select coarse-grid points and setup the interpolation
 * operator. The most common strategies are \c coarsening_direct and \c coarsening_standard.
 * An older strategy implemented in the Amg1R5 code, is sort of a trade-off between the two and
 * is identified by \c coarsening_amg1r5. Cf. section A.7.2.1 in "Multigrid"
 * [Trottenberg, Oosterlee, Schüller] and the documentation of the corresponding coarsening
 * utility functions.
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::setCoarseningStrategy(CoarseningStrategy strategy) {
	strategy_ = strategy;
}


/**\brief Sets if a full multigrid is performed or not.
 *
 * \param fmg True if a full multigrid startup should be performed. 
 * \return void
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::setStartupExecution(bool fmg) {
	fmg_ = fmg;
}


/**
 * \brief Sets the number of pre- and postsmoothing sweeps.
 *
 * \param nu1 The number of presmoothing sweeps. 
 * \param nu2 The number of postsmoothing sweeps. 
 * \return void
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::setSmoothingParameters(std::size_t nu1, std::size_t nu2) {
	nu1_ = nu1;
	nu2_ = nu2;
}


/**
 * \brief Sets the coarsening parameters.
 *
 * \param thetaPos The parameter tweaking the fraction of the positive off-diagonal entries strongly influencing the diagonal entry. 
 * \param thetaNeg The parameter tweaking the fraction of the negative off-diagonal entries strongly influencing the diagonal entry. 
 * \return void
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::setCoarseningParameters(double thetaPos, double thetaNeg) {
	thetaPos_ = thetaPos;
	thetaNeg_ = thetaNeg;
}


/**
 * \brief Provides access to the direct solver.
 *
 * \return Returns a reference to the direct solver.
 */
template<typename DirectSolver>
DirectSolver& AmgSolver<DirectSolver>::getDirectSolver() {
	return directSolver_;
}


/**
 * \brief Provides read-only access to the direct solver.
 *
 * \return Returns a read-only reference to the direct solver.
 */
template<typename DirectSolver>
const DirectSolver& AmgSolver<DirectSolver>::getDirectSolver() const {
	return directSolver_;
}


/**
 * \brief Initializes the smoother on a particular level.
 *
 * \param level The level on which the smoother is to be initialized.
 * \return void
 *
 * After setting up the system on a level the smoother can initialize its helper structures. E.g.
 * Gauss-Seidel like smoothers usually preinvert diagonal entries.
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::initSmoother(std::size_t level) {
	Grid&            grid(grids_[level]);
	std::size_t      n   (grid.system_.size());
	const MatrixCRS& A   (grid.system_.A_);

	// allocate helper data
	grid.diag_index_.resize(n, false);
	grid.diag_inv_.resize(n, false);

	// locate diagonal entries in system matrix and precompute inverses
	for (std::size_t i = 0; i < n; ++i) {
		grid.diag_index_[i] = A.find(i, i);

		assert(grid.diag_index_[i] != std::numeric_limits<double>::infinity() && "Diagonal entry missing.");
		assert(A[grid.diag_index_[i]].value_ != 0.0 && "Invalid zero diagonal element.");

		grid.diag_inv_[i] = 1.0 / A[grid.diag_index_[i]].value_;
	}
}


/**
 * \brief Splits the points on grid \a level into coarse- and fine-grid points.
 *
 * \param level The level on which the splitting should be performed.
 * \param thetaPos The parameter controlling the strong dependency threshold for positive off-diagonal entries.
 * \param thetaNeg The parameter controlling the strong dependency threshold for negative off-diagonal entries.
 * \return void
 *
 * The splitting algorithm is refered to as "standard coarsening" in "Multigrid"
 * [Trottenberg, Oosterlee, Schüller]. The resulting types of the points (\c fine_point or
 * \c coarse_point) are stored in the grid's \a children_ vector on return.
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::splitPoints(std::size_t level, double thetaPos, double thetaNeg) {
	Grid&                    grid     (grids_[level]);
	std::size_t              n        (grid.system_.size());
	const MatrixCRS&         A        (grid.system_.A_);
	std::vector<std::size_t> qualities(n, 0);

	if (loglevel_ <= INFO)
		std::cout << "Constructing strong dependency subgraph..." << std::endl;

	// mark all points as undecided
	grid.children_.clear();
	grid.children_.resize(n, undecided_point);

	// allocate enough memory for strong dependency graph adjacency lists
	grid.S_.clear();
	grid.S_.resize(n);
	grid.St_.clear();
	grid.St_.resize(n);
	for (std::size_t i = 0; i < n; ++i) {
		grid.S_ [i].clear();
		grid.S_ [i].resize(A.nonzeros(i) - 1);
		grid.St_[i].clear();
		grid.St_[i].resize(A.nonzeros(i) - 1);
	}

	// construct adjacency lists of strong dependency subgraph
	for (std::size_t i = 0; i < n; ++i) {

		// scan row for absolute maximum off-diagonal entry (for strong dependency threshold)
		double sdepThresholdPos(0), sdepThresholdNeg(0);
		for (std::size_t j = A.begin(i); j < A.end(i); ++j) {
			if (A[j].index_ != i) {
				sdepThresholdPos = std::max(sdepThresholdPos, A[j].value_);
				sdepThresholdNeg = std::min(sdepThresholdNeg, A[j].value_);
			}
		}

		if (sdepThresholdPos == 0.0 && sdepThresholdNeg == 0.0) {
			// mark as fine-grid point because it has no non-zero off-diagonal elements and does not need to be interpolated
			grid.children_[i] = fine_point;
		}

		sdepThresholdPos *= thetaPos;
		sdepThresholdNeg *= thetaNeg;

		std::size_t S_size(0);
		for (std::size_t j = A.begin(i); j < A.end(i); ++j) {
			if (A[j].index_ != i && A[j].value_ != 0 && (A[j].value_ >= sdepThresholdPos || A[j].value_ <= sdepThresholdNeg)) {
				grid.S_ [i          ][S_size++                ] = A[j].index_;
				grid.St_[A[j].index_][qualities[A[j].index_]++] = i;
			}
		}

		// shorten S_ adjacency list to its filled part
		grid.S_[i].resize(S_size);
	}

	// shorten St_ adjacency lists to their filled part
	for (std::size_t i = 0; i < n; ++i)
		grid.St_[i].resize(qualities[i]);

	if (loglevel_ <= INFO)
		std::cout << "Selecting coarse-grid points..." << std::endl;

	// setup priority queue with helper structure to accelerate modification of the item priorities
	typedef std::pair<std::size_t, std::size_t> PqueueItem;
	typedef std::set<PqueueItem> Pqueue;

	Pqueue candidates;
	std::vector<Pqueue::iterator> candidate_positions(n);

	// insert all undecided grid points with their qualities as priorities into the priority queue
	for (std::size_t i = 0; i < n; ++i)
		if (grid.children_[i] == undecided_point)
			candidate_positions[i] = candidates.insert(PqueueItem(qualities[i], i)).first;

	// add candidates to coarse grid points until no more are available
	while (!candidates.empty()) {
		// remove grid point with highest quality from priority queue
		Pqueue::iterator last(--candidates.end());
		std::size_t i(last->second);
		assert(grid.children_[i] == undecided_point && "Coarse contact priority queue corrupted.");
		if (loglevel_ <= INFO)
			std::cout << "Classified point " << i << " as coarse grid point with quality " << last->first << "." << std::endl;
		candidates.erase(last);

		// grid point i becomes new coarse grid point
		grid.children_[i] = coarse_point;

		// all points in S^T_i become fine grid points
		for (std::size_t j = 0; j < grid.St_[i].size(); ++j) {
			std::size_t p(grid.St_[i][j]);
			if (grid.children_[p] == undecided_point) {
				// p becomes fine grid point
				grid.children_[p] = fine_point;
				assert(candidate_positions[p]->second == p && "Priority queue helper structure corrupted.");
				candidates.erase(candidate_positions[p]);

				// the qualities of all undecided points on which p strongly depends are incremented 
				for (std::size_t k = 0; k < grid.S_[p].size(); ++k) {
					std::size_t q(grid.S_[p][k]);

					if (grid.children_[q] == undecided_point) {
						// increment qualities of yet undecided points by removing and readding them
						Pqueue::iterator& it(candidate_positions[q]);
						assert(it->second == q && "Priority queue helper structure corrupted.");
						std::size_t quality(it->first + 1);

						//Pqueue::iterator predecessor(it);
						//--predecessor;
						candidates.erase(it);
						//candidates_positions[q] = candidates.insert(predecessor, PqueueItem(quality, q));
						candidates.insert(PqueueItem(quality, q));
						candidate_positions[q] = candidates.find(PqueueItem(quality, q));
						if (loglevel_ <= INFO)
							std::cout << "Incremented quality of point " << q << " to " << quality << "." << std::endl;
					}
				}
			}
		}
	}
}


/**
 * \brief Coarsens a particular grid level employing the given strategy.
 *
 * \param level The grid level to coarsen.
 * \param strategy The strategy to use.
 * \return Returns true if coarsening was successful.
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::coarsen(std::size_t level, CoarseningStrategy strategy) {
	switch (strategy) {
		case coarsening_direct:
			return coarsenByDirectInterpolation(level);
		case coarsening_standard:
			return coarsenByStandardInterpolation(level);
		case coarsening_amg1r5:
			return coarsenByAmg1r5(level);
		default:
			return false;
	}
}


/**
 * \brief Coarsens a particular grid level.
 *
 * \param level The grid level to coarsen.
 * \return Returns true if coarsening was successful.
 *
 * In the coarsening process the system matrix of the given level is analyzed and coarse-grid
 * point candidates are selected. An interpolation matrix \f$P\f$ is constructed and the coarse system
 * computed using the Galerkin product \f$P^T A P\f$. This coarsening process is well suited for
 * matrices similar to M-matrices. That is the off-diagonal entries should be essentially
 * non-positive. The approach is the same as presented in "A Multigrid Tutorial (2nd edition)"
 * [Briggs, Henson, McCormick] and as implemented in AMG1R5. It is a compromise between the direct
 * and standard interpolation presented in "Multigrid" [Trottenberg, Oosterlee, Schüller]. Refer
 * to remark A.7.8 for a comparison to the other interpolations.
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::coarsenByAmg1r5(std::size_t level) {
	assert(level < grids_.size() && "Invalid grid level parameter for coarsening.");

	if (level < grids_.size() - 1)
		return true;

	if (grids_[level].system_.size() <= 1)
		return false;

	if (loglevel_ <= INFO)
		std::cout << "Coarsening grid " << level << "..." << std::endl;

	grids_.push_back(Grid());

	Grid&       grid      (grids_[level]);
	Grid&       coarsegrid(grids_[level+1]);
	std::size_t n         (grid.system_.size());
	const MatrixCRS& A         (grid.system_.A_);

	// standard C/F splitting without strong positive couplings
	splitPoints(level, 0.0, thetaNeg_);

	if (loglevel_ <= INFO) {
		std::size_t cpoints(0);
		for (std::size_t i = 0; i < n; ++i) {
			if (grid.children_[i] == coarse_point)
				cpoints++;
		}
		std::cout << cpoints << "/" << n << " coarse-grid points after first pass..." << std::endl;
	}

	if (loglevel_ <= INFO)
		std::cout << "Reclassifying points in second phase..." << std::endl;

	// 2nd phase: reconsider fine-grid points which depend strongly on other fine grid points but
	// do not have a common coarse-grid point
	for (std::size_t i = 0; i < n; ++i) {
		assert(grid.children_[i] != undecided_point && "A point was left unclassified in the coarsening process.");
		if (grid.children_[i] == coarse_point)
			continue;

		std::size_t tentative_cpoint(fine_point);

		for (std::size_t j = 0; j < grid.S_[i].size(); ++j) {
			std::size_t S_ij(grid.S_[i][j]);

			if (grid.children_[S_ij] == coarse_point)
				continue;

			// We need to check if the intersection of S_{S_ij} and C_i is empty. Note that
			// C_i is a subset of S_i and S_i as well as S_{S_ij} are sorted. Thus we can
			// compute the intersection in linear time.

			std::size_t pi(0), pj(0);
			while (pi < grid.S_[i].size() && pj < grid.S_[S_ij].size()) {
				if (grid.S_[i][pi] < grid.S_[S_ij][pj])
					++pi;
				else if (grid.S_[i][pi] > grid.S_[S_ij][pj])
					++pj;
				else if (grid.children_[grid.S_[i][pi]] != fine_point)
					// fine grid points i and j have a common C point
					break;
				else {
					// fine grid points i and j have a common F point
					++pi; ++pj;
				}
			}

			if (pi >= grid.S_[i].size() || pj >= grid.S_[S_ij].size()) {
				// no common strongly influencing coarse-grid point
				if (tentative_cpoint == fine_point) {
					// S_ij becomes tentative coarse-grid point
					if (loglevel_ <= INFO)
						std::cout << "Point " << S_ij << " tentatively becomes coarse-grid point." << std::endl;
					tentative_cpoint = S_ij;
					grid.children_[S_ij] = coarse_point;
				}
				else {
					// replace tentative coarse-grid point by i
					if (loglevel_ <= INFO)
						std::cout << "Tentative coarse-grid point " << tentative_cpoint << " replaced by " << i << "." << std::endl;
					assert(tentative_cpoint != i && "Invalid tentative coarse-grid point.");
					grid.children_[i]                = coarse_point;
					grid.children_[tentative_cpoint] = fine_point;
					tentative_cpoint = i;
					break;
				}
			}
		}
	}

	// assign ascending indices to coarse-grid points
	std::size_t cpoints(0);
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] == coarse_point)
			grid.children_[i] = cpoints++;
	}

	if (loglevel_ <= INFO)
		std::cout << "Selected " << cpoints << " out of " << n << " points as coarse-grid points." << std::endl;

	if (cpoints == n) {
		// if all points happen to be classified as coarse-grid points, we give up to prevent infinite coarsening loop
		grids_.pop_back();
		return false;
	}

	if (loglevel_ <= INFO)
		std::cout << "Setting up interpolation operator..." << std::endl;

	// allocate prolongation/interpolation operator P
	std::vector<std::size_t>& capacities(tmp_);
	capacities.resize(n);
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] != fine_point)
			capacities[i] = 1;
		else
			capacities[i] = grid.S_[i].size();  // C_i^s \subset S_i
	}

	MatrixCRS& P(coarsegrid.P_);
	P = MatrixCRS(n, cpoints, capacities);

	// construct prolongation/interpolation operator P
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] != fine_point) {
			// coarse-grid points are injected
			P.appendRowElement(i, grid.children_[i], 1.0);
			continue;
		}

		// touch all j in C_i^s and remember matrix access index in temporary array
		tmp_.clear();
		tmp_.resize(n, std::numeric_limits<std::size_t>::max());
		for (std::size_t j = 0; j < grid.S_[i].size(); ++j) {
			std::size_t S_ij(grid.S_[i][j]);

			if (grid.children_[S_ij] == fine_point)
				continue;

			P.appendRowElement(i, grid.children_[S_ij], 0.0);
			tmp_[S_ij] = P.end(i) - 1;
		}

		// distribute weakly influencing coarse- and fine-grid points to diagonal and add strongly
		// influencing coarse-grid points to prolongation operator
		double accu(0);   // sum weakly influencing neighbors (D_i^w) and diagonal entry in accu

		std::size_t jS(0);
		for (std::size_t jA = A.begin(i); jA < A.end(i);) {
			std::size_t j(A[jA].index_);

			if (jS < grid.S_[i].size() && grid.S_[i][jS] < j) {
				++jS;
			}
			else if (jS < grid.S_[i].size() && grid.S_[i][jS] == j) {
				if (grid.children_[j] != fine_point) {
					// j is in the strongly influencing neighboring coarse interpolatory set C_i^s
					P.getValue(tmp_[j]) -= A[jA].value_;
				}
				else {
					// j is in the strongly influencing neighboring fine-grid point set F_i^s
					double sum(0);
					VectorN v(P.nonzeros(i), 0.0);

					std::size_t pi(0), pj(0), kA(A.begin(j));
					while (pi < grid.S_[i].size() && pj < grid.S_[j].size()) {
						if (grid.S_[i][pi] < grid.S_[j][pj])
							++pi;
						else if (grid.S_[i][pi] > grid.S_[j][pj])
							++pj;
						else {
							if (grid.children_[grid.S_[i][pi]] != fine_point) {
								// common strongly influencing C point
								while (A[kA].index_ < grid.S_[i][pi]) {
									++kA;
									assert(kA < A.end(j) && "Strong dependency subgraph is inconsistent with matrix.");
								}
								assert(A[kA].index_ == grid.S_[i][pi] && "Strong dependency subgraph is inconsistent with matrix.");

								sum += A[kA].value_;
								v[tmp_[grid.S_[i][pi]] - P.begin(i)] = A[kA].value_;
							}
							++pi; ++pj;
						}
					}

					assert(sum != 0 && "Strongly influencing coarse grid off-diagonal entries sum to zero.");

					v *= A[jA].value_ / sum;

					for (std::size_t k = 0; k < v.size(); ++k)
						P.getValue(P.begin(i) + k) -= v[k];
				}

				++jA;
			}
			else {
				// j is a weakly influencing neighbor (D_i^w)
				accu += A[jA].value_;
				++jA;
			}
		}

		accu = 1.0 / accu;
		for (std::size_t jP = P.begin(i); jP < P.end(i); ++jP)
			P.getValue(jP) *= accu;
	}

	computeGalerkinProduct(P, A, grid.R_, coarsegrid.system_.A_);

	// allocating data on coarse grid
	coarsegrid.r_.clear();
	coarsegrid.r_.resize(cpoints);
	coarsegrid.tmp_.clear();
	coarsegrid.tmp_.resize(cpoints);
	coarsegrid.system_.x_.clear();
	coarsegrid.system_.x_.resize(cpoints);
	coarsegrid.system_.b_.clear();
	coarsegrid.system_.b_.resize(cpoints);

	if (loglevel_ <= INFO)
		std::cout << "Coarsened grid level " << level << "." << std::endl;

	return true;
}


/**
 * \brief Coarsens a particular grid level.
 *
 * \param level The grid level to coarsen.
 * \return Returns true if coarsening was successful.
 *
 * In the coarsening process the system matrix of the given level is analyzed and coarse-grid
 * point candidates are selected. An interpolation matrix \f$P\f$ is constructed and the coarse system
 * computed using the Galerkin product \f$P^T A P\f$. This coarsening process is well suited for
 * symmetric matrices with large negative \em and positive off-diagonal entries.
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::coarsenByDirectInterpolation(std::size_t level) {
	assert(level < grids_.size() && "Invalid grid level parameter for coarsening.");

	if (level < grids_.size() - 1)
		return true;

	if (grids_[level].system_.size() <= 1)
		return false;

	if (loglevel_ <= INFO)
		std::cout << "Coarsening grid " << level << "..." << std::endl;

	grids_.push_back(Grid());

	Grid&       grid      (grids_[level]);
	Grid&       coarsegrid(grids_[level+1]);
	std::size_t n         (grid.system_.size());
	const MatrixCRS& A    (grid.system_.A_);

	// standard C/F splitting
	splitPoints(level, thetaPos_, thetaNeg_);

	// assign ascending indices to coarse-grid points
	std::size_t cpoints(0);
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] == coarse_point)
			grid.children_[i] = cpoints++;
	}

	if (loglevel_ <= INFO)
		std::cout << "Selected " << cpoints << " out of " << n << " points as coarse-grid points." << std::endl;

	if (cpoints == n) {
		// if all points happen to be classified as coarse-grid points, we give up to prevent infinite coarsening loop
		grids_.pop_back();
		return false;
	}

	if (loglevel_ <= INFO)
		std::cout << "Setting up interpolation operator..." << std::endl;

	// allocate prolongation/interpolation operator P
	std::vector<std::size_t>& capacities(tmp_);
	capacities.resize(n);
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] != fine_point)
			capacities[i] = 1;
		else
			capacities[i] = grid.S_[i].size();  // C_i^s \subset S_i
	}

	MatrixCRS& P(coarsegrid.P_);
	P = MatrixCRS(n, cpoints, capacities);

	// construct prolongation/interpolation operator P using direct interpolation
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] != fine_point) {
			// coarse-grid points are injected
			P.appendRowElement(i, grid.children_[i], 1.0);
			continue;
		}

		double diag(0);
		double alpha_nom(0), alpha_denom(0);
		double beta_nom(0), beta_denom(0);

		std::size_t jS(0);
		for (std::size_t jA = A.begin(i); jA < A.end(i);) {
			std::size_t j(A[jA].index_);

			if (jS < grid.S_[i].size() && grid.S_[i][jS] < j) {
				++jS;
			}
			else if (j != i) {
				if (jS < grid.S_[i].size() && grid.S_[i][jS] == j && grid.children_[j] != fine_point) {
					// strongly influencing coarse-grid point (= interpolatory point)

					if (A[jA].value_ >= 0)
						beta_denom += A[jA].value_;
					else
						alpha_denom += A[jA].value_;

					// add matrix entry for interpolatory point to interpolation operator
					P.appendRowElement(i, grid.children_[j], A[jA].value_);
					++jS;
				}
				else {
					// either a weakly influencing neighbor (D_i^w) or a strongly influencing fine-grid point (= non-interpolatory point)
					if (A[jA].value_ >= 0)
						beta_nom += A[jA].value_;
					else
						alpha_nom += A[jA].value_;
				}

				++jA;
			}
			else {
				// diagonal entry
				diag = A[jA++].value_;
			}
		}

		alpha_nom += alpha_denom;
		beta_nom += beta_denom;

		// no interpolatory point at all
		if (alpha_denom == 0.0 && beta_denom == 0.0) {
			// TODO use standard instead of direct interpolation
			assert(false && "No interpolatory point available for interpolation. Not implemented yet: Use standard interpolation as fallback.");
		}

		// no positive interpolatory point
		if (beta_denom == 0.0)
			diag += beta_nom;

		// no negative interpolatory point
		if (alpha_denom == 0.0)
			diag += alpha_nom;

		assert(diag != 0.0 && "Modified diagonal entry became zero.");

		double coeffPos(-beta_nom / (beta_denom * diag));
		double coeffNeg(-alpha_nom / (alpha_denom * diag));

		// correct row entries of interpolation operator
		for (std::size_t jP = P.begin(i); jP != P.end(i); ++jP) {
			if (P[jP].value_ >= 0)
				P.getValue(jP) *= coeffPos;
			else
				P.getValue(jP) *= coeffNeg;
		}
	}

	computeGalerkinProduct(P, A, grid.R_, coarsegrid.system_.A_);

	// allocating data on coarse grid
	coarsegrid.r_.clear();
	coarsegrid.r_.resize(cpoints);
	coarsegrid.tmp_.clear();
	coarsegrid.tmp_.resize(cpoints);
	coarsegrid.system_.x_.clear();
	coarsegrid.system_.x_.resize(cpoints);
	coarsegrid.system_.b_.clear();
	coarsegrid.system_.b_.resize(cpoints);

	if (loglevel_ <= INFO)
		std::cout << "Coarsened grid level " << level << "." << std::endl;

	return true;
}




/**\brief Coarsens a particular grid level.
 *
 * \param level The grid level to coarsen.
 * \return Returns true if coarsening was successful.
 *
 * In the coarsening process the system matrix of the given level is analyzed and coarse-grid
 * point candidates are selected. An interpolation matrix \f$P\f$ is constructed and the coarse system
 * computed using the Galerkin product \f$P^T A P\f$. This coarsening process is well suited for
 * symmetric matrices with large negative \em and positive off-diagonal entries.
 * TODO
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::coarsenByStandardInterpolation(std::size_t level) {
	assert(level < grids_.size() && "Invalid grid level parameter for coarsening.");

	if (level < grids_.size() - 1)
		return true;

	if (grids_[level].system_.size() <= 1)
		return false;

	if (loglevel_ <= INFO)
		std::cout << "Coarsening grid " << level << "..." << std::endl;

	grids_.push_back(Grid());

	Grid&       grid      (grids_[level]);
	Grid&       coarsegrid(grids_[level+1]);
	std::size_t n         (grid.system_.size());
	const MatrixCRS& A    (grid.system_.A_);

	// standard C/F splitting
	splitPoints(level, thetaPos_, thetaNeg_);

	// assign ascending indices to coarse-grid points
	std::size_t cpoints(0);
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] == coarse_point)
			grid.children_[i] = cpoints++;
	}

	if (loglevel_ <= INFO)
		std::cout << "Selected " << cpoints << " out of " << n << " points as coarse-grid points." << std::endl;

	if (cpoints == n) {
		// if all points happen to be classified as coarse-grid points, we give up to prevent infinite coarsening loop
		grids_.pop_back();
		return false;
	}

	if (loglevel_ <= INFO)
		std::cout << "Setting up interpolation operator..." << std::endl;

	// allocate prolongation/interpolation operator P
	std::vector<std::size_t>& capacities(tmp_);
	capacities.resize(n);
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] != fine_point)
			capacities[i] = 1;
		else
			capacities[i] = grid.S_[i].size();  // C_i^s \subset S_i
	}

	MatrixCRS& P(coarsegrid.P_);
	P = MatrixCRS(n, cpoints, capacities);

	// construct prolongation/interpolation operator P using direct interpolation
	std::vector<bool> rowCstrong(n);   // Warning: do not initialize all vector elements in each iteration because it would change the complexity!
	for (std::size_t i = 0; i < n; ++i) {
		if (grid.children_[i] != fine_point) {
			// coarse-grid points are injected
			P.appendRowElement(i, grid.children_[i], 1.0);
			continue;
		}

		std::size_t interpolatoryPoints(0);

		// copy i-th row of A
		capacities.clear();
		capacities.resize(1);
		capacities[0] = n;
		MatrixCRS row(1, n, capacities);
		for (std::size_t jA = A.begin(i); jA < A.end(i); ++jA) {
			row.appendRowElement(0, A[jA].index_, A[jA].value_);
			rowCstrong[ A[jA].index_ ] = false;   // initialize
		}

		// eliminate strongly influencing fine-grid points
		for (std::size_t jS = 0; jS < grid.S_[i].size(); ++jS) {
			std::size_t j(grid.S_[i][jS]);

			if (grid.children_[j] == fine_point) {
				// eliminate e_j from row
				double coeff(-row[ row.find(0, j) ].value_ / A[ A.find(j, j) ].value_);

				MatrixCRS tmpRow(1, n, capacities);
				std::size_t jRow(row.begin(0)), jA(A.begin(j));
				while (jRow < row.end(0) && jA < A.end(j)) {
					if (row[jRow].index_ < A[jA].index_) {
						tmpRow.appendRowElement(0, row[jRow].index_, row[jRow].value_);
						++jRow;
					}
					else if (row[jRow].index_ > A[jA].index_) {
						tmpRow.appendRowElement(0, A[jA].index_, coeff * A[jA].value_);
						rowCstrong[A[jA].index_] = false;   // initialize
						++jA;
					}
					else {
						tmpRow.appendRowElement(0, A[jA].index_, row[jRow].value_ + coeff * A[jA].value_);
						++jRow; ++jA;
					}
				}

				while (jRow < row.end(0)) {
					tmpRow.appendRowElement(0, row[jRow].index_, row[jRow].value_);
					++jRow;
				}


				while (jA < A.end(j)) {
					tmpRow.appendRowElement(0, A[jA].index_, coeff * A[jA].value_);
					rowCstrong[A[jA].index_] = false;   // initialize
					++jA;
				}

				// construct union of C_i^s and C_j^s
				for (std::size_t jS2 = 0; jS2 < grid.S_[j].size(); ++jS2) {
					if (grid.children_[grid.S_[j][jS2]] != fine_point) {
						if (!rowCstrong[grid.S_[j][jS2]]) {
							rowCstrong[grid.S_[j][jS2]] = true;
							++interpolatoryPoints;
						}
					}
				}

				row = tmpRow;
			}
			else {
				if (!rowCstrong[j]) {
					rowCstrong[j] = true;
					++interpolatoryPoints;
				}
			}
		}

		assert(rowCstrong[i] == false && "Fine-grid point's interpolatory set contains itself.");

		P.reserveRowElements(i, interpolatoryPoints);

		double diag(0);
		double alpha_nom(0), alpha_denom(0);
		double beta_nom(0), beta_denom(0);

		for (std::size_t jRow = row.begin(0); jRow < row.end(0);) {
			std::size_t j(row[jRow].index_);

			if (j != i) {
				if (rowCstrong[j]) {
					// strongly influencing coarse-grid point (= interpolatory point)
					if (row[jRow].value_ >= 0)
						beta_denom += row[jRow].value_;
					else
						alpha_denom += row[jRow].value_;

					// add matrix entry for interpolatory point to interpolation operator
					P.appendRowElement(i, grid.children_[j], row[jRow].value_);
				}
				else {
					// either a weakly influencing neighbor (D_i^w) or a strongly influencing fine-grid point (= non-interpolatory point)
					if (row[jRow].value_ >= 0)
						beta_nom += row[jRow].value_;
					else
						alpha_nom += row[jRow].value_;
				}

				++jRow;
			}
			else {
				// diagonal entry
				diag = row[jRow++].value_;
			}
		}

		alpha_nom += alpha_denom;
		beta_nom += beta_denom;

		// no interpolatory point at all
		if (alpha_denom == 0.0 && beta_denom == 0.0) {
			assert(false && "No interpolatory point available for interpolation.");
		}

		// no positive interpolatory point
		if (beta_denom == 0.0)
			diag += beta_nom;

		// no negative interpolatory point
		if (alpha_denom == 0.0)
			diag += alpha_nom;

		assert(diag != 0.0 && "Modified diagonal entry became zero.");

		double coeffPos(-beta_nom / (beta_denom * diag));
		double coeffNeg(-alpha_nom / (alpha_denom * diag));

		// correct row entries of interpolation operator
		for (std::size_t jP = P.begin(i); jP != P.end(i); ++jP) {
			if (P[jP].value_ >= 0)
				P.getValue(jP) *= coeffPos;
			else
				P.getValue(jP) *= coeffNeg;
		}

		// TODO use truncation of interpolation operator to keep growth of the interpolatory set under control
	}

	computeGalerkinProduct(P, A, grid.R_, coarsegrid.system_.A_);

	// allocating data on coarse grid
	coarsegrid.r_.clear();
	coarsegrid.r_.resize(cpoints);
	coarsegrid.tmp_.clear();
	coarsegrid.tmp_.resize(cpoints);
	coarsegrid.system_.x_.clear();
	coarsegrid.system_.x_.resize(cpoints);
	coarsegrid.system_.b_.clear();
	coarsegrid.system_.b_.resize(cpoints);

	if (loglevel_ <= INFO)
		std::cout << "Coarsened grid level " << level << "." << std::endl;

	return true;
}


/**
 * \brief Computes the Galerkin product \f$ P^T A P \f$.
 *
 * \param P The interpolation operator.
 * \param A The system matrix.
 * \param R The sparse matrix which will contain the restriction operator \f$ P^T \f$ on return.
 * \param Ac The sparse matrix which will contain the Gelerkin product \f$ P^T A P \f$ on return.
 * \return void
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::computeGalerkinProduct(const MatrixCRS& P, const MatrixCRS& A, MatrixCRS& R, MatrixCRS& Ac) {
	assert(A.rows() == P.rows() && A.columns() == P.rows() && "Size of interpolation operator does not match size of system matrix.");

	if (loglevel_ <= INFO)
		std::cout << "Setting up restriction operator..." << std::endl;

	R = P.getTranspose();

	if (loglevel_ <= INFO)
		std::cout << "Computing Galerkin product..." << std::endl;
	
	Ac = R * A * P;
}


/**
 * \brief Executes a full-multigrid startup. 
 *
 * \return Returns if the solution is already converged.
 *
 * The full-multigrid startup starts on the coarsest level. The initial solution is assumed to
 * be zero. The coarsest system is solved and the error (which equals the solution in this case)
 * is prolongated to the next finer level. The procedure is repeated until the finest grid is
 * reached. To solve a system the regular multigrid scheme is applied \f$\nu_0\f$ times unless
 * it is the coarsest system, where the usual direct solver is applied.
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::startup() {
	assert(grids_.size() > 0 && "No grid initialized.");

	bool converged(false);

	// restrict right-hand-sides
	for (std::size_t level = 1; level < grids_.size(); ++level) {
		Grid& finegrid(grids_[level - 1]);
		Grid& coarsegrid(grids_[level]);

		for (std::size_t j = 0; j < coarsegrid.system_.size(); ++j)
			coarsegrid.system_.b_[j] = finegrid.R_.multiplyRow(j, finegrid.system_.b_);
	}

	// solve coarsest level
	grids_.back().system_.x_ = 0.0;
	cycle(grids_.size() - 1);

	// start-up
	for (std::size_t level = grids_.size() - 2; level < grids_.size(); --level) {
		Grid& finegrid(grids_[level]);
		Grid& coarsegrid(grids_[level + 1]);

		// prolongate error
		for (std::size_t i = 0; i < finegrid.system_.x_.size(); ++i)
			finegrid.system_.x_[i] = coarsegrid.P_.multiplyRow(i, coarsegrid.system_.x_);

		for (unsigned int i = 0; i < nu0_; ++i)
			converged = cycle(level);
	}

	if (loglevel_ <= INFO) {
		std::cout << grids_[0].system_.residual() << " (residual on finest level after FMG startup)" << std::endl;
		std::cout << "Full Multigrid Start-up finished." << std::endl;
	}

	return converged;
}




/**\brief Executes a multigrid cycle starting vom grid \a level .
 *
 * \param level The grid level from where to start the multigrid cycle.
 * \return Returns \c true if the residual drops below the desired threshold.
 *
 * TODO
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::cycle(std::size_t level)
{
	assert(level < grids_.size() && "Invalid grid level.");

	double residual;
	double r0(0.0);

	if (level == grids_.size() - 1)
		return solveDirect(level);

	if (loglevel_ <= INFO) {
		r0 = grids_[level].system_.residual();
		std::cout << std::string(level, '\t') << r0 << " (residual on grid " << level << " before presmoothing)" << std::endl;
	}

	// presmoothing
	for (std::size_t i = 0; i < nu1_; ++i)
		smooth(level);

	Grid& finegrid(grids_[level]);
	Grid& coarsegrid(grids_[level + 1]);

	if (loglevel_ <= INFO)
		std::cout << std::string(level, '\t') << finegrid.system_.residual() << " (residual on grid " << level << " before descent)" << std::endl;

	// compute LSE residual
	for (std::size_t j = 0; j < finegrid.system_.size(); ++j)
		finegrid.r_[j] = finegrid.system_.b_[j] - finegrid.system_.A_.multiplyRow(j, finegrid.system_.x_);

	// set initial solution to 0
	coarsegrid.system_.x_ = 0.0;

	// compute right-hand-side
	for (std::size_t j = 0; j < coarsegrid.system_.size(); ++j)
		coarsegrid.system_.b_[j] = finegrid.R_.multiplyRow(j, finegrid.r_);

	// recursive descent
	for (std::size_t i = 0; i < mu_; ++i) {
		// solve the coarse-grid residual equation A^2h*e^2h = r^2h
		if (!cycle(level + 1)) {
			if (loglevel_ <= INFO)
				std::cout << std::string(level, '\t') << "WARNING: Coarse-grid was not solved accurately!" << std::endl;
		}
	}

	// prolongate error
	for (std::size_t j = 0; j < finegrid.system_.size(); ++j)
		finegrid.system_.x_[j] += coarsegrid.P_.multiplyRow(j, coarsegrid.system_.x_);

	if (loglevel_ <= INFO)
		std::cout << std::string(level, '\t') << grids_[level].system_.residual() << " (residual on grid " << level << " after descent)" << std::endl;

	// postsmoothing
	for (std::size_t i = 0; i < nu2_; ++i)
		smooth(level);

	if (loglevel_ <= INFO) {
		residual = grids_[level].system_.residual();
		std::cout << std::string(level, '\t') << residual << " (residual on grid " << level << " after postsmoothing)" << std::endl;

		double convergence(residual / r0);
		if (convergence >= 1.0)
			std::cout << std::string(level, '\t') << convergence << " (convergence)" << std::endl;
		else if (convergence >= 0.9)
			std::cout << std::string(level, '\t') << convergence << " (convergence)" << std::endl;
		else
			std::cout << std::string(level, '\t') << convergence << " (convergence)" << std::endl;
	}

		return grids_[level].system_.residual() < threshold_;
}


/**
 * \brief Uses the direct solver to solve grid \a level .
 *
 * \return Returns true if the direct solver was able to solve the given system.
 */
template<typename DirectSolver>
bool AmgSolver<DirectSolver>::solveDirect(std::size_t level) {
	bool converged(false);

	converged = directSolver_.solve(grids_[level].system_);

	if (loglevel_ <= INFO) {
		std::cout << std::string(level, '\t') << grids_[level].system_.residual() << " (residual on level " << level << " after execution of direct solver)" << std::endl;
		if (converged)
			std::cout << std::string(level, '\t') << "Solved the LSE on grid " << level << "." << std::endl;
		else
			std::cout << std::string(level, '\t') << "WARNING: Did not solve the LSE on grid " << level << "." << std::endl;
	}

	return converged;
}


/**
 * \brief A Gauß-Seidel smoother.
 *
 * \param level The grid level to smooth.
 *
 * TODO
 */
template<typename DirectSolver>
void AmgSolver<DirectSolver>::smooth(std::size_t level) {
	Grid&       grid    (grids_[level]      );
	MatrixCRS&  A       (grid.system_.A_    );
	VectorN&    x       (grid.system_.x_    );
	VectorN&    b       (grid.system_.b_    );
	VectorN&    diag_inv(grid.diag_inv_     );
	std::size_t n       (grid.system_.size());

	for (std::size_t i = 0; i < n; ++i) {
		const double residual(b[i] - A.multiplyRow(i, x));
		x[i] += diag_inv[i] * residual;
	}
}


template<typename DirectSolver>
void AmgSolver<DirectSolver>::outputSquareGrid(std::ostream& out, std::size_t level) const {
	const Grid& l0(grids_[0]);
	const Grid& l(grids_[level]);

	std::size_t n(static_cast<std::size_t>(std::floor(std::sqrt(l0.system_.size()) + 0.5)) + 1);
	assert((n - 1) * (n - 1) == l0.system_.size() && "Invalid size.");
	double h = 1.0 / n;

	for (std::size_t i = 0; i < l0.system_.size(); ++i) {
		double x = (1 + i % (n - 1)) * h;
		double y = (1 + i / (n - 1)) * h;

		std::size_t p(i);
		for (std::size_t j = 0; j < level && p != fine_point; ++j)
			p = grids_[j].children_[p];

		if (p != fine_point) {
			out << x << " " << y << " " << l.system_.x_[p] << '\n';
		}
	}
}


template<typename DirectSolver>
void AmgSolver<DirectSolver>::outputCoarseningGraph(std::ostream& out, const std::string& name, std::size_t level) const {
	const Grid& grid(grids_[level]);
	const std::size_t n(grid.system_.size());

	out << "graph " << name << " {" << '\n';
	out << std::string(1, '\t') << "overlap = scale;\n";
	out << std::string(1, '\t') << "{" << '\n';
	out << std::string(2, '\t') << "node [shape = octagon]\n";
	for (std::size_t i = 0; i < n; ++i)
		if (grid.children_[i] == undecided_point)
			out << std::string(2, '\t') << i << '\n';
	out << std::string(1, '\t') << "}\n";

	out << std::string(1, '\t') << "{" << '\n';
	out << std::string(2, '\t') << "node [shape = ellipse]";
	for (std::size_t i = 0; i < n; ++i)
		if (grid.children_[i] == fine_point)
			out << " " << i;
	out << '\n';
	out << std::string(1, '\t') << "}\n";

	out << std::string(1, '\t') << "{" << '\n';
	out << std::string(2, '\t') << "node [shape = box]";
	for (std::size_t i = 0; i < n; ++i)
		if (grid.children_[i] != fine_point && grid.children_[i] != undecided_point)
			out << " " << i;
	out << '\n';
	out << std::string(1, '\t') << "}\n";

	typedef std::pair<std::size_t, std::size_t> edge_t;
	std::set<edge_t> edges;
	std::set<edge_t> dedges;

	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t jS = 0; jS < grid.S_[i].size(); ++jS) {
			std::size_t j(grid.S_[i][jS]);
			if (j < i) {
				edge_t flipedge(j, i);
				if (edges.find(flipedge) != edges.end()) {
					edges.erase(flipedge);
					dedges.insert(flipedge);
				}
				else {
					edges.insert(flipedge);
				}
			}
			else {
				edges.insert(edge_t(i, j));
			}
		}
	}

	out << std::string(1, '\t') << "{\n";
	out << std::string(2, '\t') << "edge [ style = solid ]\n";
	for (std::set<edge_t>::iterator it = edges.begin(); it != edges.end(); ++it)
		out << std::string(2, '\t') << it->first << " -- " << it->second << ";\n";
	out << std::string(1, '\t') << "}\n";

	out << std::string(1, '\t') << "{\n";
	out << std::string(2, '\t') << "edge [ style = bold ]\n";
	for (std::set<edge_t>::iterator it = dedges.begin(); it != dedges.end(); ++it)
		out << std::string(2, '\t') << it->first << " -- " << it->second << ";\n";
	out << std::string(1, '\t') << "}\n";

	out << std::string(1, '\t') << "{\n";
	out << std::string(2, '\t') << "edge [style = dotted]\n";
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t jA = grid.system_.A_.find(i, i) + 1; jA != grid.system_.A_.end(i); ++jA) {
			if (edges.find(edge_t(i, grid.system_.A_[jA].index_))  == edges.end() &&
					dedges.find(edge_t(i, grid.system_.A_[jA].index_)) == dedges.end()) {
				out << std::string(2, '\t') << i << " -- " << grid.system_.A_[jA].index_ << ";\n";
			}
		}
	}
	out << std::string(1, '\t') << "}\n";


	out << "}\n";
}


template<typename DirectSolver>
std::string AmgSolver<DirectSolver>::formatString(const std::string& format, std::size_t n) {
	std::size_t first = format.find('%');
	if (first == std::string::npos)
		first = format.size();

	std::size_t space = 0;
	while (first + space < format.size() && format[first + space] == '%') {
		++space;
	}
	//if (first == format.size())
	//   return format;

	std::ostringstream sstr;
	sstr << format.substr(0, first);
	sstr << std::setw(space) << std::setfill('0') << std::right << n;
	sstr << format.substr(first + space, format.size() - first - space);

	return sstr.str();
}

#endif
