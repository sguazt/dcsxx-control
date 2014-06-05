/**
 * \file src/dcs/control/solver/dare_slicot5.hpp
 *
 * \brief Discrete-time Algebraic Riccati Equations solver.
 *
 * Computes the unique stabilizing solution \f$X\f$ of the discrete-time
 * algebraic Riccati equation
 * \f[
 *   A^{T}XA - E^{T}XE - (A^{T}XB+S)(B^{T}XB+R)^{-1}(B^{T}X A+S)^{T} + Q = 0
 * \f]
 * or, equivalently (if \f$R\f$ is nonsingular):
 * \f[
 *   E^{T}XE = F^{T}XF-F^{T}XB(B^{T}XB+R)^{-1}B^{T}XF + Q - SR^{-1}S^{T}
 * \f]
 * where \f$F=A-BR^{-1}S^{T}\f$.
 *
 * Beside the solution \f$X\f$, the solver also computes the \e gain matrix
 * \f[                       -1
 *   G = (B^{T}XB + R)^{-1}(B^{T}XA + S^{T})
 * \f]
 * and the vector \f$L\f$ of <em>closed-loop eigenvalues</em> associated to the
 * generalized Schur decomposition of matrix pencil \f$(A-BG,E)\f$.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright 2009 Marco Guazzone (marco.guazzone@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef DCS_CONTROL_SOLVER_DARE_SLICOT5_HPP
#define DCS_CONTROL_SOLVER_DARE_SLICOT5_HPP


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublasx/container/generalized_diagonal_matrix.hpp>
#include <boost/numeric/ublasx/operation/cat.hpp>
#include <boost/numeric/ublasx/operation/diag.hpp>
#include <boost/numeric/ublasx/operation/lu.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/ql.hpp>
#include <boost/numeric/ublasx/operation/qz.hpp>
#include <boost/numeric/ublasx/operation/rcond.hpp>
#include <boost/numeric/ublasx/operation/seq.hpp>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <stdexcept>


namespace dcs { namespace control {

namespace ublas = ::boost::numeric::ublas;
namespace ublasx = ::boost::numeric::ublasx;

template <typename ValueT>
class dare_solver
{
	public: typedef ValueT value_type;
	public: typedef typename ublas::type_traits<value_type>::real_type real_type;
	private: typedef ublas::matrix<value_type,ublas::column_major> work_matrix_type;
	private: typedef ublas::vector< ::std::complex<value_type> > work_vector_type;
	public: typedef work_matrix_type matrix_type;
	public: typedef work_vector_type vector_type;
	private: typedef typename ublas::matrix_traits<work_matrix_type>::size_type size_type;


	// Default constructor
	public: dare_solver()
	{
		// empty
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT
			>
		dare_solver(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q)
	{
		solve(A, B, Q);
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT,
				typename RMatrixExprT
			>
		dare_solver(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q, ublas::matrix_expression<RMatrixExprT> const& R)
	{
		solve(A, B, Q, R);
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT,
				typename RMatrixExprT,
				typename SMatrixExprT
			>
		dare_solver(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q, ublas::matrix_expression<RMatrixExprT> const& R, ublas::matrix_expression<SMatrixExprT> const& S)
	{
		solve(A, B, Q, R, S);
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT,
				typename RMatrixExprT,
				typename SMatrixExprT,
				typename EMatrixExprT
			>
		dare_solver(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q, ublas::matrix_expression<RMatrixExprT> const& R, ublas::matrix_expression<SMatrixExprT> const& S, ublas::matrix_expression<EMatrixExprT> const& E)
	{
		solve(A, B, Q, R, S, E);
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT
			>
		void solve(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q)
	{
		work_matrix_type R = ublas::identity_matrix<value_type>(ublas::num_columns(B));

		solve(A, B, Q, R);
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT,
				typename RMatrixExprT
			>
		void solve(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q, ublas::matrix_expression<RMatrixExprT> const& R)
	{
		work_matrix_type S = ublas::zero_matrix<value_type>(ublas::num_rows(B), ublas::num_columns(B));

		solve(A, B, Q, R, S);
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT,
				typename RMatrixExprT,
				typename SMatrixExprT
			>
		void solve(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q, ublas::matrix_expression<RMatrixExprT> const& R, ublas::matrix_expression<SMatrixExprT> const& S)
	{
		work_matrix_type E = ublas::identity_matrix<value_type>(ublas::num_rows(B));

		solve(A, B, Q, R, S, E);
	}


	public: template <
				typename AMatrixExprT,
				typename BMatrixExprT,
				typename QMatrixExprT,
				typename RMatrixExprT,
				typename SMatrixExprT,
				typename EMatrixExprT
			>
		void solve(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q, ublas::matrix_expression<RMatrixExprT> const& R, ublas::matrix_expression<SMatrixExprT> const& S, ublas::matrix_expression<EMatrixExprT> const& E)
	{
		size_type n = ublasx::num_rows(B);
		size_type m = ublasx::num_columns(B);

		// precondition: n > 0 && m > 0
		DCS_ASSERT(
			n > 0 && m > 0,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] Empty matrix.")
		);
		// precondition: A is square && num_rows[A] == num_rows[B]
		DCS_ASSERT(
			ublasx::num_rows(A) == n && ublasx::num_columns(A) == n,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] A must be a square matrix and must have  the same number of rows of B.")
		);
		// precondition: Q is square && size[Q] == size[A]
		DCS_ASSERT(
			ublasx::num_rows(Q) == n && ublasx::num_columns(Q) == n,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] A and Q must be the same size.")
		);
		// precondition: Q is symmetric
		DCS_ASSERT(
			ublas::norm_1(Q-ublas::trans(Q)) <= 100*eps*ublas::norm_1(Q),
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] Q must be symmetric.")
		);
		// precondition: R is square && num_rows[R] == num_columns[B]
		DCS_ASSERT(
			ublasx::num_rows(R) == m && ublasx::num_columns(R) == m,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be a square matrix and must have  the same number of columns of B.")
		);
		// precondition: R is symmetric
		DCS_ASSERT(
			ublas::norm_1(R-ublas::trans(R)) <= 100*eps*ublas::norm_1(R),
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be symmetric.")
		);
		// precondition: size[S] == size[B]
		DCS_ASSERT(
			ublasx::num_rows(S) == n && ublasx::num_columns(S) == m,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] B and S must be the same size.")
		);
		// precondition: E is square && size[E] == size[A]
		DCS_ASSERT(
			ublasx::num_rows(E) == n && ublasx::num_columns(E) == n,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] A and Q must be the same size.")
		);
		// precondition: E != I_n && E is not singular
		DCS_ASSERT(
			/*E != ublas::identity_matrix(n) && */ublasx::rcond(E) >= eps,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be nonsingular.")
		);


		// reset state
		X_.resize(0, 0, false);
		G_.resize(0, 0, false);

		// Create temporary copies of Q, R, S matrices and make sure that Q and R
		// are symmetric matrices.
		work_matrix_type work_Q((Q + ublas::trans(Q))/real_type(2));
		work_matrix_type work_R((R + ublas::trans(R))/real_type(2));
		work_matrix_type work_S(S);

		// Scale Q,R,S so that norm_1(Q)+norm_1(R)+norm_1(S) = 1
		real_type scale_factor(0);
		scale_factor += ublas::norm_1(work_Q)+ublas::norm_1(work_R)+ublas::norm_1(S);
		// Chek if the scaling factor is too near to zero.
		if (scale_factor <= eps)
		{
			scale_factor = real_type(1);
		}
		work_Q /= scale_factor;
		work_R /= scale_factor;
		work_S /= scale_factor;

		// Set up extended Symplectic pencil (H,J) of the form:
		//             |A   0   B|     |E   0   0|
		//   H - z J = |Q  -E'  S| - z |0  -A'  0|
		//             |S'  0   R|     |0  -B'  0|
		// of dimension (2N+M)x(2N+M), where N and M are the number of rows and of
		// of columns of matrix B, respectively.
		work_matrix_type H; // (2*n+m,2*n+m)
		work_matrix_type J; // (2*n+m,2*n+m)
		H = ublasx::cat_columns(
				ublasx::cat_columns(
					ublasx::cat_rows(
						ublasx::cat_rows(
							A,
							ublas::zero_matrix<value_type>(n,n)
						),
						B
					),
					ublasx::cat_rows(
						ublasx::cat_rows(
							work_Q,
							-ublas::trans(E)
						),
						work_S
					)
				),
				ublasx::cat_rows(
					ublasx::cat_rows(
						ublas::trans(work_S),
						ublas::zero_matrix<value_type>(m, n)
					),
					work_R
				)
			);
		J = ublasx::cat_columns(
				ublasx::cat_columns(
					ublasx::cat_rows(
						E,
						ublas::zero_matrix<value_type>(n, n+m)
					),
					ublasx::cat_rows(
						ublas::zero_matrix<value_type>(n, n),
						ublasx::cat_rows(
							-ublas::trans(A),
							ublas::zero_matrix<value_type>(n, m)
						)
					)
				),
				ublasx::cat_rows(
					ublasx::cat_rows(
						ublas::zero_matrix<value_type>(m, n),
						-ublas::trans(B)
					),
					ublas::zero_matrix<value_type>(m)
				)
			);

		// Undo scaling of the data matrices.
		//work_Q *= scale_factor;
		//work_R *= scale_factor;
		//work_S *= scale_factor;

		// Temporary copies of Q, R, and S not needed any more.
		work_Q.resize(0, 0, false);
		work_R.resize(0, 0, false);
		work_S.resize(0, 0, false);

		// Compress the pencils using QL factorization.
		// The compressed form is (see [1]):
		//   \lambda H - J
		// where H and J are 2n-by-2n matrices.
		// [1] Van Dooren, P.
		//     A Generalized Eigenvalue Approach for Solving Riccati Equations.
		//     SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981.
		//
		size_type n2 = n*2;
		size_type n2pm = n2 + m;
		ublasx::ql_decomposition<value_type> ql(ublas::subrange(H, 0, n2pm, n2, n2pm));
		work_matrix_type tmp_HJ;
		tmp_HJ = ublas::subrange(H, 0, n2pm, 0, n2);
		H = ublas::subrange(ql.tlprod(tmp_HJ), 0, n2, 0, n2);
		tmp_HJ = ublas::subrange(J, 0, n2pm, 0, n2);
		J = ublas::subrange(ql.tlprod(tmp_HJ), 0, n2, 0, n2);

		// Check the singularity of the L factor in the QL factorization:
		// if singular, then the extended matrix pencil is also singular.
		if (ublasx::rcond(ql.L()) <= eps)
		{
			DCS_DEBUG_TRACE("The extended matrix pencil is singular");
			throw ::std::runtime_error("[dcs::control::dare_solver::solve] Cannot compute DARE solution: the extended matrix pencil is singular.");
		}

		// The natural tendency of the QZ algorithm to get the largest
		// eigenvalues in the leading part of the matrix pair is
		// exploited, by computing the unstable eigenvalues (i.e., the ones
		// outside the unit circle) of the permuted matrix pair (J,H).
		// This is equivalent to take the eigenvalues inside the unit circle of
		// the matrix pair (H,J).

		//ublasx::qz_decomposition<work_matrix_type,work_matrix_type> qz;
		//qz.decompose(H, J, ublasx::udi_qz_eigenvalues);
		work_matrix_type Z;
		//ublasx::qz_decompose_inplace(J, H, qzQ, qzZ, ublasx::udo_qz_eigenvalues);//FIXME:test me
		//TODO: clear qzQ since it seems it is not used
		ublasx::qz_decomposition<value_type> qz;
		qz.decompose(J, H);
		qz.reorder(ublasx::udo_qz_eigenvalues);
		Z = qz.Z();
		// Compute the n stable closed-loop eigenvalues of the system matrix
		// A-BG, where G is the optimal feedback matrix computed based on the
		// solution matrix X.
		// These eigenvalues correspond to the the trailing n generalized
		// eigenvalues of the QZ decomposition
		l_ = ublas::subrange(qz.eigenvalues(), n, n2);

		// Select submatrices X1 and X2 out of the matrix Z which define the
		// solution X = X2 * inv(X1).
		// Since X = X' we may obtain X as the solution of the system of
		// linear equations X1' X = X2', where
		//    X1 = Z(1:n, 1:n),
		//    X2 = Z(n+1:2n, 1:n).
		work_matrix_type X1 = ublas::subrange(Z, 0, n, 0, n);
		work_matrix_type X2 = ublas::subrange(Z, n, n2, 0, n);

		// Z is not needed any more
		Z.resize(0, 0, false);

		//// Check if X1 is singular
		//X1_norm = ublas::norm_1(X1);

		// Solve the system X1' X = X2'
		//work_matrix_type X;
		size_type sing;
		sing = ublasx::lu_solve(ublas::trans(X1), ublas::trans(X2), X_);
	//	ublas::permutation_matrix<size_type> P(n);
	//	X1 = ublas::trans(X1);
	//	ublasx::lu_decompose_inplace(X1, P);
	//	if (ublasx::rcond(X1) < eps)
	//	{
	//		// Nearly singular matrix
	//		//FIXME: what to do?
	//	}
	//	X = ublasx::lu_apply(X1, P, ublas::trans(X2));
		if (sing)
		{
			// Nearly singular matrix
			DCS_DEBUG_TRACE("Cannot compute DARE solution: nearly singular matrix.");
			throw ::std::runtime_error("[dcs::control::dare_solver::solve] Cannot compute DARE solution: nearly singular matrix.");
		}

		// Make sure the solution X is symmetrix
		X_ = scale_factor * (X_ + ublas::trans(X_))/real_type(2);

		// Compute the Gain matrix G = (B'*X*B+R)\(B'*X*A+S');
		work_matrix_type BTX(ublas::prod(ublas::trans(B), X_));
		sing = ublasx::lu_solve(
				ublas::prod(BTX, B) + R,
				ublas::prod(BTX, A) + ublas::trans(S),
				G_
		);
		if (sing)
		{
			// Nearly singular matrix
			DCS_DEBUG_TRACE("Cannot compute DARE gain matrix: nearly singular matrix.");
			throw ::std::runtime_error("[dcs::control::dare_solver::solve] Cannot compute DARE gain matrix: nearly singular matrix.");
		}
	}


	public: matrix_type solution() const
	{
		return X_;
	}


	/// Return the gain matrix: \f$(B^T X B + R)^{-1}( B^T X A + S^T)\f$.
	public: matrix_type gain() const
	{
		return G_;
	}


	/**
	 * \brief Return the closed-loop eigenvalues vector.
	 *
	 * The closed-loop eigenvalues vector \f$\lambda\f$ is computed as:
	 * \f[
	 *   \lambda = \operatorname{eig}(A-BG, E)
	 * \f]
	 * where \f$G\f$ is the gain matrix.
	 */
	public: vector_type eigenvalues() const
	{
		return l_;
	}


	private: static const real_type eps;


	private: work_matrix_type X_; /// The solution matrix
	private: work_matrix_type G_; /// The gain matrix
	private: work_vector_type l_; /// The closed-loop eigenvalues vector.
};


template <typename ValueT>
const typename ublas::type_traits<ValueT>::real_type dare_solver<ValueT>::eps = ::std::numeric_limits<typename ublas::type_traits<ValueT>::real_type>::epsilon();


/**
 * \brief Solve the discrete-time algebraic Ricciati equation (DARE).
 *
 * Computes the unique stabilizing solution X of the discrete-time algebraic
 * Riccati equation
 * \f[
 *   A^{T}XA-E^{T}XE-(A^{T}XB+S)(B^{T}XB+R)^{-1}(B^{T}XA+S)+Q=0
 * \f]
 * The dare function also returns the gain matrix
 * \f[
 *   G=(B^{T}XB+R)^{-1}(B^{T}XA+S)
 * \f]
 * and the vector L of closed loop eigenvalues, where
 * \f[
 *   L=\operatorname{eigen}(A-B*G,E)
 * \f]
 *
 * \author Marco Guazzone, marco.guazzone@mfn.unipmn.it
 */
template <
	typename AMatrixExprT,
	typename BMatrixExprT,
	typename QMatrixExprT,
	typename RMatrixExprT,
	typename SMatrixExprT,
	typename EMatrixExprT,
	typename XMatrixT
>
void dare(ublas::matrix_expression<AMatrixExprT> const& A, ublas::matrix_expression<BMatrixExprT> const& B, ublas::matrix_expression<QMatrixExprT> const& Q, ublas::matrix_expression<RMatrixExprT> const& R, ublas::matrix_expression<SMatrixExprT> const& S, ublas::matrix_expression<EMatrixExprT> const& E, XMatrixT& X)
{
#if 0
	typedef AMatrixExprT A_matrix_type;
	typedef BMatrixExprT B_matrix_type;
	typedef QMatrixExprT Q_matrix_type;
	typedef RMatrixExprT R_matrix_type;
	typedef SMatrixExprT S_matrix_type;
	typedef EMatrixExprT E_matrix_type;
	typedef XMatrixT X_matrix_type;
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<A_matrix_type>::value_type,
				typename ublas::promote_traits<
					typename ublas::matrix_traits<B_matrix_type>::value_type,
					typename ublas::promote_traits<
						typename ublas::matrix_traits<Q_matrix_type>::value_type,
						typename ublas::promote_traits<
							typename ublas::matrix_traits<R_matrix_type>::value_type,
								typename ublas::promote_traits<
									typename ublas::matrix_traits<S_matrix_type>::value_type,
									typename ublas::promote_traits<
										typename ublas::matrix_traits<E_matrix_type>::value_type,
										typename ublas::matrix_traits<X_matrix_type>::value_type
									>::promote_type
								>::promote_type
							>::promote_type
						>::promote_type
					>::promote_type
				>::promote_type value_type;
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<A_matrix_type>::size_type,
				typename ublas::promote_traits<
					typename ublas::matrix_traits<B_matrix_type>::size_type,
					typename ublas::promote_traits<
						typename ublas::matrix_traits<Q_matrix_type>::size_type,
						typename ublas::promote_traits<
							typename ublas::matrix_traits<R_matrix_type>::size_type,
								typename ublas::promote_traits<
									typename ublas::matrix_traits<S_matrix_type>::size_type,
									typename ublas::promote_traits<
										typename ublas::matrix_traits<E_matrix_type>::size_type,
										typename ublas::matrix_traits<X_matrix_type>::size_type
									>::promote_type
								>::promote_type
							>::promote_type
						>::promote_type
					>::promote_type
				>::promote_type size_type;
	typedef typename ublas::type_traits<value_type>::real_type real_type;
	typedef ublas::matrix<value_type> work_matrix_type; //FIXME: adjust orientation?

	size_type n = ublasx::num_rows(B);
	size_type m = ublasx::num_columns(B);

	const double eps = ::std::numeric_limits<double>::epsilon();

	// precondition: n > 0 && m > 0
	DCS_ASSERT(
		n > 0 && m > 0,
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] Empty matrix.")
	);
	// precondition: A is square && num_rows[A] == num_rows[B]
	DCS_ASSERT(
		ublasx::num_rows(A) == n && ublasx::num_columns(A) == n,
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] A must be a square matrix and must have  the same number of rows of B.")
	);
	// precondition: Q is square && size[Q] == size[A]
	DCS_ASSERT(
		ublasx::num_rows(Q) == n && ublasx::num_columns(Q) == n,
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] A and Q must be the same size.")
	);
	// precondition: Q is symmetric
	DCS_ASSERT(
		ublas::norm_1(Q-ublas::trans(Q)) <= 100*eps*ublas::norm_1(Q),
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] Q must be symmetric.")
	);
	// precondition: R is square && num_rows[R] == num_columns[B]
	DCS_ASSERT(
		ublasx::num_rows(R) == m && ublasx::num_columns(R) == m,
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be a square matrix and must have  the same number of columns of B.")
	);
	// precondition: R is symmetric
	DCS_ASSERT(
		ublas::norm_1(R-ublas::trans(R)) <= 100*eps*ublas::norm_1(R),
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be symmetric.")
	);
	// precondition: size[S] == size[B]
	DCS_ASSERT(
		ublasx::num_rows(S) == n && ublasx::num_columns(S) == m,
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] B and S must be the same size.")
	);
	// precondition: E is square && size[E] == size[A]
	DCS_ASSERT(
		ublasx::num_rows(E) == n && ublasx::num_columns(E) == n,
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] A and Q must be the same size.")
	);
	// precondition: E != I_n && E is not singular
	DCS_ASSERT(
		/*E != ublas::identity_matrix(n) && */ublasx::rcond(E) >= eps,
		throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be nonsingular.")
	);

	// Create temporary copies of Q, R, S matrices and make sure that Q and R
	// are symmetric matrices.
	work_matrix_type work_Q((Q + ublas::trans(Q))/real_type(2));
	work_matrix_type work_R((R + ublas::trans(R))/real_type(2));
	work_matrix_type work_S(S);

	// Scale Q,R,S so that norm_1(Q)+norm_1(R)+norm_1(S) = 1
	real_type scale_factor(0);
	scale_factor += ublas::norm_1(work_Q)+ublas::norm_1(work_R)+ublas::norm_1(S);
	// Chek if the scaling factor is too near to zero.
	if (scale_factor <= eps)
	{
		scale_factor = real_type(1);
	}
	work_Q /= scale_factor;
	work_R /= scale_factor;
	work_S /= scale_factor;

	// Set up extended Symplectic pencil (H,J) of the form:
	//             |A   0   B|     |E   0   0|
	//   H - z J = |Q  -E'  S| - z |0  -A'  0|
	//             |S'  0   R|     |0  -B'  0|
	// of dimension (2N+M)x(2N+M), where N and M are the number of rows and of
	// of columns of matrix B, respectively.
	work_matrix_type H; // (2*n+m,2*n+m)
	work_matrix_type J; // (2*n+m,2*n+m)
	H = ublasx::cat_columns(
			ublasx::cat_columns(
				ublasx::cat_rows(
					ublasx::cat_rows(
						A,
						ublas::zero_matrix<value_type>(n,n)
					),
					B
				),
				ublasx::cat_rows(
					ublasx::cat_rows(
						work_Q,
						-ublas::trans(E)
					),
					work_S
				)
			),
			ublasx::cat_rows(
				ublasx::cat_rows(
					ublas::trans(work_S),
					ublas::zero_matrix<value_type>(m, n)
				),
				work_R
			)
		);
	J = ublasx::cat_columns(
			ublasx::cat_columns(
				ublasx::cat_rows(
					E,
					ublas::zero_matrix<value_type>(n, n+m)
				),
				ublasx::cat_rows(
					ublas::zero_matrix<value_type>(n, n),
					ublasx::cat_rows(
						-ublas::trans(A),
						ublas::zero_matrix<value_type>(n, m)
					)
				)
			),
			ublasx::cat_rows(
				ublasx::cat_rows(
					ublas::zero_matrix<value_type>(m, n),
					-ublas::trans(B)
				),
				ublas::zero_matrix<value_type>(m)
			)
		);

	// Undo scaling of the data matrices.
	//work_Q *= scale_factor;
	//work_R *= scale_factor;
	//work_S *= scale_factor;

	// Temporary copies of Q, R, and S not needed any more.
	work_Q.resize(0, 0, false);
	work_R.resize(0, 0, false);
	work_S.resize(0, 0, false);

	// Compress the pencils using QL factorization.
	// The compressed form is (see [1]):
	//   \lambda H - J
	// where H and J are 2n-by-2n matrices.
	// [1] Van Dooren, P.
	//     A Generalized Eigenvalue Approach for Solving Riccati Equations.
	//     SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981.
	//
	size_type n2 = n*2;
	size_type n2pm = n2 + m;
	ublasx::ql_decomposition<value_type> ql(ublas::subrange(H, 0, n2pm, n2, n2pm));
	work_matrix_type tmp_HJ;
	tmp_HJ = ublas::subrange(H, 0, n2pm, 0, n2);
	H = ublas::subrange(ql.tlprod(tmp_HJ), 0, n2, 0, n2);
	tmp_HJ = ublas::subrange(J, 0, n2pm, 0, n2);
	J = ublas::subrange(ql.tlprod(tmp_HJ), 0, n2, 0, n2);

	// Check the singularity of the L factor in the QL factorization:
	// if singular, then the extended matrix pencil is also singular.
	if (ublasx::rcond(ql.L()) <= eps)
	{
		DCS_DEBUG_TRACE("The extended matrix pencil is singular");
		throw ::std::runtime_error("[dcs::control::dare] Cannot compute DARE solution: the extended matrix pencil is singular.");
	}

	// The natural tendency of the QZ algorithm to get the largest
	// eigenvalues in the leading part of the matrix pair is
	// exploited, by computing the unstable eigenvalues (i.e., the ones
	// outside the unit circle) of the permuted matrix pair (J,H).
	// This is equivalent to take the eigenvalues inside the unit circle of
	// the matrix pair (H,J).

	//ublasx::qz_decomposition<work_matrix_type,work_matrix_type> qz;
	//qz.decompose(H, J, ublasx::udi_qz_eigenvalues);
	work_matrix_type Z;
	//ublasx::qz_decompose_inplace(J, H, qzQ, qzZ, ublasx::udo_qz_eigenvalues);//FIXME:test me
	//TODO: clear qzQ since it seems it is not used
	ublasx::qz_decomposition<value_type> qz;
	qz.decompose(J, H);
	qz.reorder(ublasx::udo_qz_eigenvalues);
	Z = qz.Z();

	// Select submatrices X1 and X2 out of the matrix Z which define the
	// solution X = X2 * inv(X1).
	// Since X = X' we may obtain X as the solution of the system of
	// linear equations X1' X = X2', where
	//    X1 = Z(1:n, 1:n),
	//    X2 = Z(n+1:2n, 1:n).
	work_matrix_type X1 = ublas::subrange(Z, 0, n, 0, n);
	work_matrix_type X2 = ublas::subrange(Z, n, n2, 0, n);

	// Z is not needed any more
	Z.resize(0, 0, false);

	//// Check if X1 is singular
	//X1_norm = ublas::norm_1(X1);

	// Solve the system X1' X = X2'
	work_matrix_type X;
	size_type sing;
	sing = ublasx::lu_solve(ublas::trans(X1), ublas::trans(X2), X);
//	ublas::permutation_matrix<size_type> P(n);
//	X1 = ublas::trans(X1);
//	ublasx::lu_decompose_inplace(X1, P);
//	if (ublasx::rcond(X1) < eps)
//	{
//		// Nearly singular matrix
//		//FIXME: what to do?
//	}
//	X = ublasx::lu_apply(X1, P, ublas::trans(X2));
	if (sing)
	{
		// Nearly singular matrix
		DCS_DEBUG_TRACE("Cannot compute DARE solution: nearly singular matrix.");
		throw ::std::runtime_error("[dcs::control::dare] Cannot compute DARE solution: nearly singular matrix.");
	}

	// Make sure the solution X is symmetrix
	X = scale_factor * (X + ublas::trans(X))/real_type(2);
#endif

	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<AMatrixExprT>::value_type,
				typename ublas::promote_traits<
					typename ublas::matrix_traits<BMatrixExprT>::value_type,
					typename ublas::promote_traits<
						typename ublas::matrix_traits<QMatrixExprT>::value_type,
						typename ublas::promote_traits<
							typename ublas::matrix_traits<RMatrixExprT>::value_type,
								typename ublas::promote_traits<
									typename ublas::matrix_traits<SMatrixExprT>::value_type,
									typename ublas::promote_traits<
										typename ublas::matrix_traits<EMatrixExprT>::value_type,
										typename ublas::matrix_traits<XMatrixT>::value_type
									>::promote_type
								>::promote_type
							>::promote_type
						>::promote_type
					>::promote_type
				>::promote_type value_type;

	dare_solver<value_type> solver;
	solver.solve(A, B, Q, R, S, E);

	X = solver.solution();
}


}} // Namespace dcs::control


#endif // DCS_CONTROL_SOLVER_DARE_SLICOT5_HPP
