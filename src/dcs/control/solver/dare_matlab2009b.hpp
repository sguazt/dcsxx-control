/**
 * \file src/dcs/control/solver/dare.hpp
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
 * <hr/>
 *
 * Copyright (C) 2009-2010  Distributed Computing System (DCS) Group, Computer
 * Science Department - University of Piemonte Orientale, Alessandria (Italy).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */

#ifndef DCS_CONTROL_DARE_HPP
#define DCS_CONTROL_DARE_HPP


#include <algorithm>
#include <boost/mpl/if.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublasx/container/generalized_diagonal_matrix.hpp>
#include <boost/numeric/ublasx/container/sequence_vector.hpp>
#include <boost/numeric/ublasx/operation/abs.hpp>
#include <boost/numeric/ublasx/operation/any.hpp>
#include <boost/numeric/ublasx/operation/arithmetic_ops.hpp>
#include <boost/numeric/ublasx/operation/balance.hpp>
#include <boost/numeric/ublasx/operation/cat.hpp>
#include <boost/numeric/ublasx/operation/diag.hpp>
#include <boost/numeric/ublasx/operation/empty.hpp>
#include <boost/numeric/ublasx/operation/find.hpp>
#include <boost/numeric/ublasx/operation/log2.hpp>
#include <boost/numeric/ublasx/operation/lu.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/pow2.hpp>
#include <boost/numeric/ublasx/operation/qr.hpp>
#include <boost/numeric/ublasx/operation/qz.hpp>
#include <boost/numeric/ublasx/operation/rcond.hpp>
#include <boost/numeric/ublasx/operation/relational_ops.hpp>
#include <boost/numeric/ublasx/operation/rep.hpp>
#include <boost/numeric/ublasx/operation/reshape.hpp>
#include <boost/numeric/ublasx/operation/round.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <boost/numeric/ublasx/operation/tril.hpp>
#include <boost/numeric/ublasx/operation/triu.hpp>
#include <boost/type_traits/is_complex.hpp>
#include <cmath>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <functional>
#include <iostream>
#include <stdexcept>


namespace dcs { namespace control {

namespace ublas = ::boost::numeric::ublas;
namespace ublasx = ::boost::numeric::ublasx;


namespace detail { namespace /*<unnamed>*/ {

template <typename MatrixExprT, typename SVectorT, typename PVectorT>
typename ublasx::balance_traits<MatrixExprT>::balanced_matrix_type mat_scale(ublas::matrix_expression<MatrixExprT> const& A,
																			 SVectorT& s,
																			 PVectorT& p,
																			 bool full_balance,
																			 bool permute)
{
	typedef typename ublasx::balance_traits<MatrixExprT>::balanced_matrix_type balanced_matrix_type;
	typedef typename ublas::matrix_traits<MatrixExprT>::size_type size_type;

	balanced_matrix_type B;

	// Compute scaling
	if (full_balance)
	{
		B = ublasx::balance(A, s, true, false);
	}
	else
	{
		//TODO: see MATLAB mscale.m (use of matscale function)
		throw ::std::runtime_error("[dcs::control::detail::mat_scale] Regularized row/column balancing not yet supported.");
	}

	// Compute permutation
	size_type n(ublasx::num_rows(A));
	if (ublasx::size(p) != n)
	{
		p.resize(n, false);
	}
	if (permute)
	{
		if (ublas::norm_1(ublasx::tril(A,-2)) == 0)
		{
			// Real upper Hessenberg: no permutation
			p = ublasx::sequence_vector<size_type>(0, n);
		}
		else if (ublas::norm_1(ublasx::triu(A,2)) == 0)
		{
			// Real lower Hessenberg
			p = ublasx::sequence_vector<size_type>(n-1, -1, n);
		}
		else
		{
			SVectorT tmp_s;
			ublasx::balance(A, tmp_s, p, true, true);
		}

		// Permute rows and columns of B according to p
		balanced_matrix_type tmp_B(ublasx::num_rows(B), ublasx::num_columns(B));
		for (size_type i = 0; i < n; ++i)
		{
			size_type k(p[i]);

			ublas::row(tmp_B, i) = ublas::row(B, k);
		}
		for (size_type i = 0; i < n; ++i)
		{
			size_type k(p[i]);

			ublas::column(B, i) = ublas::column(tmp_B, k);
		}
	}
	else
	{
		// Sequential order
		p = ublasx::sequence_vector<size_type>(0, n);
	}

	return B;
}


template <typename AMatrixExprT, typename LVectorExprT, typename RVectorExprT>
ublas::matrix<
	typename ublas::matrix_traits<AMatrixExprT>::value_type
> left_right_mat_scale(ublas::matrix_expression<AMatrixExprT> const& A,
					   ublas::vector_expression<LVectorExprT> const& l,
					   ublas::vector_expression<RVectorExprT> const& r)
{
	typedef typename ublas::matrix_traits<AMatrixExprT>::value_type value_type;
	typedef typename ublas::matrix_traits<AMatrixExprT>::size_type size_type;

	ublas::matrix<value_type> X(A);

	size_type nr(ublasx::num_rows(A));
	size_type nc(ublasx::num_columns(A));

	if (!ublasx::empty(l))
	{
		X = ublas::element_prod(
					ublasx::rep(l, 1, nc),
					X
				);
	}
	if (!ublasx::empty(r))
	{
		X = ublas::element_prod(
				X,
				ublasx::rep(
						ublasx::reshape(r, 1, nc),
						nr,
						1
					)
			);
	}

	return X;
}


template <typename AMatrixExprT, typename LMatrixExprT, typename RMatrixExprT>
ublas::matrix<
	typename ublas::matrix_traits<AMatrixExprT>::value_type
> left_right_mat_scale(ublas::matrix_expression<AMatrixExprT> const& A,
					   ublas::matrix_expression<LMatrixExprT> const& L,
					   ublas::matrix_expression<RMatrixExprT> const& R)
{
	typedef typename ublas::matrix_traits<AMatrixExprT>::value_type value_type;
	typedef typename ublas::matrix_traits<AMatrixExprT>::size_type size_type;

	ublas::matrix<value_type> X(A);

	size_type nr(ublasx::num_rows(A));
	size_type nc(ublasx::num_columns(A));

	if (!ublasx::empty(L))
	{
		X = ublas::element_prod(
					ublasx::rep(
							ublasx::reshape(L, nr, 1),
							1,
							nc
						),
					X
				);
	}
	if (!ublasx::empty(R))
	{
		X = ublas::element_prod(
				X,
				ublasx::rep(
						ublasx::reshape(R, nc, 1),
						nr,
						1
					)
			);
	}
}


template <typename HMatrixT, typename JMatrixT, typename LVectorT, typename X1MatrixT, typename X2MatrixT, typename XMatrixT>
void gdare(HMatrixT& H, JMatrixT& J, ::std::size_t n, ::std::size_t m, LVectorT& l, X1MatrixT& X1, X2MatrixT& X2, XMatrixT& X, bool balance, bool factorize)
{
	//   - inputs: H, J, n, balance, factorize
	//   - outputs:
	//      - X1, X2, sx, L (if 'factorize' == true)
	//      - X, L (otherwise)
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<HMatrixT>::size_type,
				typename ublas::matrix_traits<JMatrixT>::size_type>::promote_type size_type;
	typedef typename ublas::promote_traits<
				typename ublas::matrix_traits<HMatrixT>::value_type,
				typename ublas::matrix_traits<JMatrixT>::value_type>::promote_type value_type;
	typedef typename ublas::type_traits<value_type>::real_type real_type;
	typedef ublas::matrix<value_type> work_matrix_type;
	typedef ublas::vector<value_type> work_vector_type;

	size_type n2 = n*2;
	size_type n2pm = n2+m;
//	size_type n2p1 = n2+1;
//	size_type mp1 = m+1;

	const real_type eps = ::std::numeric_limits<real_type>::epsilon();

	// Scale Symplectic matrix/pencil before compression to preserve
	// Symplectic structure (D = state matrix scaling).

	ublas::vector<real_type> sx; // scaling vector
	ublas::vector<size_type> p; // permutation vector
	if (balance)
	{
		// BEGIN of MATLAB arescale

		// Form matrix M to be balanced
		work_matrix_type M;
		real_type tol_zero(1e2*eps);
		if (ublasx::empty(J))
		{
			M = H;
		}
		else
		{
			work_matrix_type Hd;
			work_matrix_type Jd;
			real_type nh;
			real_type nj;

			//Hd = H - ublasx::diag(ublasx::diag(H));
			Hd = H - ublasx::generalized_diagonal_adaptor<HMatrixT>(H);
			nh = ublas::norm_1(
						ublas::subrange(Hd, 0, n, 0, n)
					)
					+
					ublasx::norm_1(
						ublas::subrange(Hd, n, n2, n, n2)
					);
			//Jd = J-ublasx::diag(ublasx::diag(J));
			Jd = J - ublasx::generalized_diagonal_adaptor<JMatrixT>(J);
			nj = ublas::norm_1(
						ublas::subrange(Jd, 0, n, 0, n)
					)
					+
					ublasx::norm_1(
						ublas::subrange(Jd, n, n2, n, n2)
					);
			if (nh > 0 && nj > 0)
			{
				M = nj*ublasx::abs(H) + nh*ublasx::abs(J);
			}
			else
			{
				M = ublasx::abs(H) + ublasx::abs(J);
			}
		}

		// Small parasitic entries can trick scaling into making X1 nearly
		// singular (see tare:lvlTwo_Hinf1).
		// Zero out such entries in the magnitude matrix M
		// Also helps identifying near-triangularizing permutation, e.g.,
		//  [1 eps;1/eps 1] -> (ignore eps) -> [1 1/eps;eps 1]
		if (m == 0)
		{
			work_matrix_type Mu(n2,2);
			ublas::column(Mu, 1) = ublasx::abs(ublasx::diag(M));
			ublas::column(Mu, 2) = ublas::scalar_vector<real_type>(n2, 1);
			Mu = tol_zero * (Mu+ublas::trans(Mu));
			size_type nr_M(ublasx::num_rows(M));
			size_type nc_M(ublasx::num_columns(M));
			for (size_type r = 0; r < nr_M; ++r)
			{
				for (size_type c = 0; c < nc_M; ++c)
				{
					if (::std::abs(M(r,c)) < Mu(r,c))
					{
						M(r,c) = real_type/*zero*/();
					}
				}
			}
		}

		// Rescale magnitude matrix.
		// 1. Use two-step approach to acquire permutation, see g162709
		// 2. Acquire permutation making H(perm,perm) more upper triangular.
		//    This enhances numerics when H can be permuted to nearly
		//    upper-triangular (the Schur method may underperform the eigenvalue
		//    one without this permutation).
		// 3. Full balancing needed when F,G<<1, but can hurt when A has small
		//    entries (roundoff) and F=0 or G=0.
		ublas::vector<real_type> s;
		ublas::vector<size_type> pp;
		detail::mat_scale(M, s, pp, true, true);
		pp = ublasx::find(pp, ::std::bind2nd(::std::less<size_type>(), n2));
		p.resize(ublasx::size(pp), false);
		for (size_type i = 0; i < n2; ++i)
		{
			p[pp[i]] = i;
		}
		// unconstrained balancing diag(D1,D2,DR)
		s = ublasx::log2(s);

		// Impose the constraint that diagonal scalings must be of the form 
		// diag(D,1./D,DR). 
		sx = ublasx::round(
					(ublas::subrange(s, n, n2)-ublas::subrange(s, 0, n))/real_type(2)
				);  // D=sqrt(D1/D2)
		ublas::subrange(s, 0, n) = ublasx::pow2(sx);
		ublas::subrange(s, n, n2) = ublasx::pow2(-sx);
		ublas::subrange(s, n2, n2+m) = ublasx::pow2(-ublas::subrange(s, n2, n2+m));
		sx = ublas::subrange(s, 0, n); // n-by-1 vector

		// Rescale H,J and return diagonal scaling of state matrix
		H = detail::left_right_mat_scale(H, s, real_type(1)/s);
		if (!ublasx::empty(J))
		{
		   J = detail::left_right_mat_scale(J, s, real_type(1)/s);
		}

		// END of MATLAB arescale
	}
	else
	{
		// No balance

		sx = ublas::scalar_vector<value_type>(n, 1);
		p = ublasx::sequence_vector<size_type>(0, n);
	}

	work_matrix_type E;

	E = ublas::subrange(J, 0, n, 0, n);

	// Compression step on H(:,n2+1:n2+m) = [S1;-S2;R]
	if (m > 0)
	{
		ublasx::qr_decomposition<value_type> qr;

		qr.decompose(ublas::subrange(H, 0, n2pm, n2, n2pm));
		work_matrix_type tmp_Qt;
		tmp_Qt = ublas::trans(
						ublas::subrange(qr.Q(), 0, n2pm, m, n2pm)
					),
		H = ublas::prod(
				tmp_Qt,
				ublas::subrange(H, 0, n2pm, 0, n2)
			);
		J = ublas::prod(
				tmp_Qt,
				ublas::subrange(J, 0, n2pm, 0, n2)
			);
	}

	// QZ algorithm
	//
	// NOTE: usual formulation is in term of the pencil (H,J), but
	//       generalized eigenvalues have a tendency to deflate out in the
	//       "desired" order, so work with (J.H) instead.
	//       Specifically, the natural tendency of the QZ algorithm to get
	//       the largest eigenvalues in the leading part of the matrix pair
	//       is exploited, by computing the unstable eigenvalues (i.e., the
	//       ones outside the unit circle) of the permuted matrix pair
	//       (J,H).
	//       This is equivalent to take the eigenvalues inside the unit
	//       circle of the matrix pair (H,J).

	size_type n_p(ublasx::size(p));
	work_matrix_type HH(n_p,n_p);
	work_matrix_type JJ(n_p,n_p);
	// permute H and J
	for (size_type r = 0; r < n_p; ++r)
	{
		for (size_type c = 0; c < n_p; ++c)
		{
			HH(r,c) = H(p[r],p[c]);
			JJ(r,c) = J(p[r],p[c]);
		}
	}
	// QZ decompose
	ublasx::qz_decomposition<value_type> qz;
	qz.decompose(JJ, HH);
	// Reorder eigenvalues to push eigenvalues outside the unit circle
	// (inside for (H,J)) to the top
	qz.reorder(ublasx::udo_qz_eigenvalues);
	JJ = qz.S();
	HH = qz.T();
	work_matrix_type ZZ;
	work_matrix_type Z(n_p,n_p);
	ZZ = qz.Z();
	for (size_type i = 0; i < n_p; ++i)
	{
		ublas::row(Z, i) = ublas::row(ZZ, p[i]);
	}
	// Compute the n stable closed-loop eigenvalues of the system matrix
	// A-BG, where G is the optimal feedback matrix computed based on the
	// solution matrix X.
	// These eigenvalues correspond to the the trailing n generalized
	// eigenvalues of the QZ decomposition
	l = qz.eigenvalues();

	// Account for non-identity E matrix and orthonormalize basis
	if (E != ublas::identity_matrix<value_type>(n,n))
	{
		work_matrix_type EZ;
		EZ = ublas::prod(E, ublas::subrange(Z, 0, n, 0, n));
		EZ.resize(n2, n);
		ublas::subrange(EZ, n, n2, 0, n) = ublas::subrange(Z, n, n2, 0, n);
		ublasx::qr_decomposition<value_type> qr;
		qr.decompose(EZ);
		Z = ublas::subrange(qr.Q(), 0, n2, 0, n);
	}

//	work_matrix_type X1;
//	work_matrix_type X2;
	X1 = ublas::subrange(Z, 0, n, 0, n);
	X2 = ublas::subrange(Z, n, n2, 0, n);

	// Check that stable invariant subspace was properly extracted
	{
		// BEGIN of MATLAB arecheckout

		work_matrix_type X12;
		X12 = ublas::prod(ublas::trans(X1), X2);
		// Solution asymmetry
		real_type asym(ublas::norm_1(X12-ublas::trans(X12)));

		work_vector_type al(ublasx::abs(l));
//			l = ublasx::find(
//						ublasx::abs(l),
//						::std::bind2nd(::std::greater<real_type>(), real_type(1))
//					);

		bool has_abs_le_1;
		bool has_abs_gt_1;
		// any(!abs(l) > 1)
		has_abs_le_1 = ublasx::any(
								ublas::subrange(al, 0, n),
								::std::bind2nd(
										::std::less_equal<real_type>(),
										real_type(1)
									)
							);
		has_abs_gt_1 = ublasx::any(
								ublas::subrange(al, n, n2),
								::std::bind2nd(
										::std::greater<real_type>(),
										real_type(1)
									)
							);
		if (has_abs_le_1
			|| has_abs_gt_1
			|| asym > ::std::max(1.0e3*eps, 0.1*ublas::norm_1(X12)))
		{
			// Could not (reliably) isolate stable invariant subspace of
			// dimension n
			//report = -1;
			X1 = X2
			   = work_matrix_type();
			sx = work_vector_type();

			::std::clog << "[Warning] Unable to solve the specified Riccati equation because the Symplectic spectrum is too near the imaginary axis." << ::std::endl;

			return;
		}
		else
		{
			//report = 0;
			if (asym > ::std::sqrt(eps))
			{
				::std::clog << "[Warning] Solution may be inaccurate due to poor scaling or eigenvalues near the stability boundary." << ::std::endl;
			}
		}
		// END of MATLAB arecheckout

		// Last n eigenvalues are inside the unit circle
		//l = ublas::conj(ublas::subrange(l, n, n2));
		if (::boost::is_complex<value_type>::value)
		{
			// eig(H,J)=[t,1/conj(t)], |t|<1 => L = eig(J,H)=[1/t,conj(t)] 
			//                               => t = conj(L(n+1:2*n))

			l = ublas::conj(ublas::subrange(l, n, n2));
		}
		else
		{
			l = ublas::subrange(l, n, n2);
		}
	}

	if (factorize)
	{
		//TODO
	}
	else
	{
		// Compute the Riccati solution X=sx*(X2/X1)*sx
		// Solve the system X*X1=X2 (== X1' X = X2')
		size_type sing;
		sing = ublasx::lu_solve(ublas::trans(X1), ublas::trans(X2), X);
		if (sing)
		{
			// X1 is singular

			X = work_matrix_type();

			::std::clog << "[Warning] The E matrix must be nonsingular." << ::std::endl;
			return;
		}
		else
		{
			// Symmetrize
			X = (X+ublas::trans(X))/real_type(2);
			// Factor in scaling sx (X -> sx X sx)
			X = detail::left_right_mat_scale(X,sx,sx);
		}
	}
}

}} // Namespace detail::<unnamed>


/**
 * \brief Solver for the discrete-time algebraic Ricciati equation (DARE).
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
 * where \f$S\f$ is the cross-coupling matrix.
 * Usually, \f$R\f$, \f$S\f$ and \f$E\f$ are set to the values
 * \f{align}{
 *  R &= I,\\
 *  S &= 0,\\
 *  E &= I
 * \f}
 * resulting in the following Riccati equation form:
 * \f[
 *   A^{T}XA-X-A^{T}XB(B^{T}XB+I)^{-1}B^{T}XA+Q=0
 * \f]
 *
 * \note
 *  Inspired by the \c dare MATLAB function.
 *
 * \author Marco Guazzone, marco.guazzone@mfn.unipmn.it
 */
template <typename ValueT>
class dare_solver
{
	public: typedef ValueT value_type;
	public: typedef typename ublas::type_traits<value_type>::real_type real_type;
	private: typedef ublas::matrix<value_type,ublas::column_major> work_matrix_type;
	private: typedef ublas::vector<
						typename ::boost::mpl::if_<
							::boost::is_complex<value_type>,
							value_type,
							::std::complex<value_type>
						>::type
					> work_vector_type;
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
		// BEGIN of MATLAB dare

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
			ublas::norm_1(Q-ublas::trans(Q)) <= 100*eps_*ublas::norm_1(Q),
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] Q must be symmetric.")
		);
		// precondition: R is square && num_rows[R] == num_columns[B]
		DCS_ASSERT(
			ublasx::num_rows(R) == m && ublasx::num_columns(R) == m,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be a square matrix and must have  the same number of columns of B.")
		);
		// precondition: R is symmetric
		DCS_ASSERT(
			ublas::norm_1(R-ublas::trans(R)) <= 100*eps_*ublas::norm_1(R),
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
			/*E != ublas::identity_matrix(n) && */ublasx::rcond(E) >= eps_,
			throw ::std::invalid_argument("[dcsxx::control::operation::dare] R must be nonsingular.")
		);


		// reset state
		X_.resize(0, 0, false);
		G_.resize(0, 0, false);


		// Set up extended Symplectic pencil (H,J) of the form:
		//             | A  0   B|     |E  0  0|
		//   H - z J = |-Q  E' -S| - z |0  A' 0|
		//             | S' 0   R|     |0 -B' 0|
		// of dimension (2N+M)x(2N+M), where N and M are the number of rows and
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
							-Q,
							ublas::trans(E)
						),
						-S
					)
				),
				ublasx::cat_rows(
					ublasx::cat_rows(
						ublas::trans(S),
						ublas::zero_matrix<value_type>(m,n)
					),
					R
				)
			);
		J = ublasx::cat_columns(
				ublasx::cat_columns(
					ublasx::cat_rows(
						E,
						ublas::zero_matrix<value_type>(n,n+m)
					),
					ublasx::cat_rows(
						ublas::zero_matrix<value_type>(n,n),
						ublasx::cat_rows(
							ublas::trans(A),
							ublas::zero_matrix<value_type>(n,m)
						)
					)
				),
				ublasx::cat_rows(
					ublasx::cat_rows(
						ublas::zero_matrix<value_type>(m,n),
						-ublas::trans(B)
					),
					ublas::zero_matrix<value_type>(m)
				)
			);

		work_matrix_type X1;
		work_matrix_type X2;
		detail::gdare(H, J, n, m, l_, X1, X2, X_, true, false);

		// Compute gain matrix G=(B'*X*B+R)\(B'*X*A+S')
		work_matrix_type BtX(ublas::prod(ublas::trans(B), X_));
		size_type sing;
		sing = ublasx::lu_solve(
				ublas::prod(BtX, B) + R,
				ublas::prod(BtX, A) + ublas::trans(S),
				G_
		);
		if (sing)
		{
			// Nearly singular matrix
			DCS_DEBUG_TRACE("Cannot compute DARE gain matrix: nearly singular matrix.");
			throw ::std::runtime_error("[dcs::control::dare_solver::solve] Cannot compute DARE gain matrix: nearly singular matrix.");
		}

//		// Compute relative residual
//		// - T1 = A'*X*A - E'*X*E
//		// - T2 = (A'*X*B + S)*G
//		work_matrix_type AtX(ublas::prod(ublas::trans(A), X_));
//		T1_ = ublas::prod(AtX, A) - ublas::prod(ublas::prod(ublas::trans(E), X), E);
//		T2_ = ublas::prod((ublas::prod(AtX, B) + S), G_);
//		work_matrix_type Res = T1 - T2 + Q;
//		report = ublas::norm_1(Res)/(1+ublas::norm_1(T1_)+ublas::norm_1(T2_)+ublas::norm_1(Q));

		// END of MATLAB dare
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


	private: static const real_type eps_;


	private: work_matrix_type X_; /// The solution matrix
	private: work_matrix_type G_; /// The gain matrix
	private: work_vector_type l_; /// The closed-loop eigenvalues vector.
};


template <typename ValueT>
const typename ublas::type_traits<ValueT>::real_type dare_solver<ValueT>::eps_ = ::std::numeric_limits<typename ublas::type_traits<ValueT>::real_type>::epsilon();


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
 * where \f$S\f$ is the cross-coupling matrix.
 * Usually, \f$R\f$, \f$S\f$ and \f$E\f$ are set to the values
 * \f{align}{
 *  R &= I,\\
 *  S &= 0,\\
 *  E &= I
 * \f}
 * resulting in the following Riccati equation form:
 * \f[
 *   A^{T}XA-X-A^{T}XB(B^{T}XB+I)^{-1}B^{T}XA+Q=0
 * \f]
 *
 * \note
 *  Inspired by the \c dare MATLAB function.
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


#endif // DCS_CONTROL_DARE_HPP
