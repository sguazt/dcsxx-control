/**
 * \file dcs/control/dlqr.hpp
 *
 * \brief Discrete Linear Quadratic Regulator.
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

#ifndef DCS_CONTROL_DLQR_HPP
#define DCS_CONTROL_DLQR_HPP


#include <algorithm>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublasx/operation/cat.hpp>
#include <boost/numeric/ublasx/operation/eigen.hpp>
#include <boost/numeric/ublasx/operation/max.hpp>
#include <boost/numeric/ublasx/operation/min.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <dcs/debug.hpp>
#include <dcs/control/dare.hpp>
#include <stdexcept>


namespace dcs { namespace control {

/*
template <typename MatrixT>
void dlqr(MatrixT const& A, MatrixT& const B)
{
	MatrixT Q = ::dcs::math::la::identity_matrix<value_type>(n);
	MatrixT R = ::dcs::math::la::zero_matrix<value_type>(n);
	MatrixT S = ::dcs::math::la::identity_matrix<value_type>(n);

	dlqr(A, B, Q, R, S);
}

template <typename MatrixT>
void dlqr(MatrixT const& A, MatrixT& const B, MatrixT& const Q)
{
	typedef typename ::dcs::math::la::matrix_traits<MatrixT>::size_type size_type;

	size_type n = ::dcs::math::la::num_rows(B);
	size_type m = ::dcs::math::la::num_columns(B);
}
*/

/**
 * \brief Discrete Linear Quadratic Controller design.
 *
 * \param A The state matrix with size \f$n \times n\f$.
 * \param B The input matrix with size \f$n \times r\f$.
 * \param C The output matrix with size \f$m \times n\f$.
 * \param Q The error weighted matrix with size \f$n \times n\f$; must be a
 *  positive semidefinite symmetric matrix (i.e., \f$Q = Q^T \ge 0\f$).
 * \param R The control weighted matrix with size \f$r \times r\f$; must be a
 *  positive definite symmetric matrix (i.e., \f$R = R^T > 0\f$).
 * \param N The cross-term weighted matrix with size \f$n \times r\f$; must be a
 *  positive semidefinite symmetric matrix (i.e, \f$N = N^T \ge 0\f$).
 *
 * Given a dicrete time-invariant state-space system model:
 * \f{align*}{
 *   x(k+1) &= Ax(k)+Bu(k) \\
 *   y(k+1) &= Cx(k+1)
 * \f}
 * the Linear Quadratic Controller problem calculates the optimal gain matrix
 * K such that the state-feedback law
 * \[
 *   u(k) = -Kx(k)
 * \]
 * minimizes the following quadratic cost function:
 * \[
 *   J(u) = \sum_{k=1}^{\infty}{x^T(k)Qx(k) + u^T(k)Ru(k) + 2x^T(k)Nu(k)}
 * \]
 */
template <
	typename Matrix1T,
	typename Matrix2T,
	typename Matrix3T,
	typename Matrix4T,
	typename Matrix5T
>
void dlqr(boost::numeric::ublas::matrix_expression<Matrix1T> const& A, boost::numeric::ublas::matrix_expression<Matrix2T> const& B, boost::numeric::ublas::matrix_expression<Matrix3T> const& Q, boost::numeric::ublas::matrix_expression<Matrix4T> const& R, boost::numeric::ublas::matrix_expression<Matrix5T> const& N)
{
	namespace ublasx = ::boost::numeric::ublasx;
	namespace ublas = ::boost::numeric::ublas;

	typedef Matrix1T A_matrix_type;
	typedef Matrix2T B_matrix_type;
	typedef Matrix3T Q_matrix_type;
	typedef Matrix4T R_matrix_type;
	typedef Matrix5T N_matrix_type;
	typedef typename ublasx::promote_traits<
				typename ublasx::matrix_traits<A_matrix_type>::value_type,
				typename ublasx::promote_traits<
					typename ublasx::matrix_traits<B_matrix_type>::value_type,
					typename ublasx::promote_traits<
						typename ublasx::matrix_traits<Q_matrix_type>::value_type,
						typename ublasx::promote_traits<
							typename ublasx::matrix_traits<R_matrix_type>::value_type,
							typename ublasx::matrix_traits<N_matrix_type>::value_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type value_type;
	typedef typename ublasx::promote_traits<
				typename ublasx::matrix_traits<A_matrix_type>::size_type,
				typename ublasx::promote_traits<
					typename ublasx::matrix_traits<B_matrix_type>::size_type,
					typename ublasx::promote_traits<
						typename ublasx::matrix_traits<Q_matrix_type>::size_type,
						typename ublasx::promote_traits<
							typename ublasx::matrix_traits<R_matrix_type>::size_type,
							typename ublasx::matrix_traits<N_matrix_type>::size_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type size_type;

	size_type A_nr = ublasx::num_rows(A);
	size_type A_nc = ublasx::num_columns(A);
	size_type B_nr = ublasx::num_rows(B);
	size_type B_nc = ublasx::num_columns(B);
	size_type Q_nr = ublasx::num_rows(Q);
	size_type Q_nc = ublasx::num_columns(Q);
	size_type R_nr = ublasx::num_rows(Q);
	size_type R_nc = ublasx::num_columns(Q);
	size_type N_nr = ublasx::num_rows(N);
	size_type N_nc = ublasx::num_columns(N);

	// precondition: A is square
	DCS_ASSERT(
		A_nr == A_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] State matrix A must be a square matrix.")
	);
	// precondition: num_rows(B) == num_rows(A)
	DCS_ASSERT(
		B_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqr] The number of rows of the input matrix B must be the same of the state matrix A.")
	);
	// precondition: Q is square
	DCS_ASSERT(
		Q_nr == Q_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] Error weighted matrix Q must be a square matrix.")
	);
	// precondition: num_rows(Q) == num_rows(A)
	DCS_ASSERT(
		Q_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqr] The number of rows of the error weighted matrix Q must be the same of the state matrix A.")
	);
	// precondition: R is square
	DCS_ASSERT(
		R_nr == R_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] Control weighted matrix Q must be a square matrix.")
	);
	// precondition: num_rows(R) == num_rows(B)
	DCS_ASSERT(
		R_nr == A_nr,
		throw ::std::invalid_argument("[dcs::control::dlqr] The number of rows of the control weighted matrix Q must be the same of the control matrix B.")
	);
	// precondition: num_rows(N) == num_rows(A) && num_columns(N) == num_columns(B)
	DCS_ASSERT(
		N_nr == A_nr && N_nc == B_nc,
		throw ::std::invalid_argument("[dcs::control::dlqr] The cross-term weighted matrix N must have the same number of rows of the state matrix A and be the same number of columns of the control matrix B.")
	);


	const value_type eps = ::std::numeric_limits<value_type>::epsilon();

	Q_matrix_type tmp_Q;
	R_matrix_type tmp_R;

	// Check/enforce symmetry and check positivity

	if (ublas::norm_1(ublas::trans(Q)-Q) > (100*eps*ublas::norm_1(Q)))
	{
		DCS_DEBUG_TRACE("The error weighted matrix Q is not symmetric and has been replaced by (Q+Q')/2.");
	//	tmp_Q = (Q+ublas::trans(Q))/static_cast<value_type>(2);
	}
	//else
	//{
	//	tmp_Q = Q;
	//}
	tmp_Q = (Q+ublas::trans(Q))/static_cast<value_type>(2);
	if (ublas::norm_1(ublas::trans(R)-R) > (100*eps*ublas::norm_1(R)))
	{
		DCS_DEBUG_TRACE("The control weighted matrix R is not symmetric and has been replaced by (R+R')/2.");
	//	tmp_R = (R+ublas::trans(R))/static_cast<value_type>(2);
	}
	//else
	//{
	//	tmp_R = R;
	//}
	tmp_R = (R+ublas::trans(R))/static_cast<value_type>(2);


	ublas::vector< std::complex<value_type> > v;

	ublas::vector<value_type> v_R;
	ublasx::eigenvalues(tmp_R, v);
	v_R = ublas::real(v);

	//FIXME: why using real(...)?
	//ublas::vector< std::complex<value_type> > v_QNR;
	ublas::vector<value_type> v_QNR;
	ublasx::eigenvalues(
		ublasx::cat_columns(
			ublasx::cat_rows(tmp_Q, N), 
			ublasx::cat_rows(ublas::trans(N), tmp_R) 
		),
		v
	);
	v_QNR = ublas::real(v);

	if (ublasx::min(v_R) <= value_type()/*zero*/)
	{
		throw ::std::runtime_error("[dcs::control::dlqr] The control weighted matrix R is not positive definite.");
	}
	else if (ublasx::min(v_QNR) < (-1.0e+2*eps*::std::max(value_type()/*zero*/, ublasx::max(v_QNR))))
	{
		DCS_DEBUG_TRACE("The matrix [Q N;N' R] is not positive semi-definite.");
	}

	ublas::matrix<value_type> E = ublas::identity_matrix<value_type>(A_nr);

	dare_solver<value_type> solver;
	solver.solve(A, B, Q, R, N, E);
std::cerr << "Solution = " << solver.solution() << std::endl;
std::cerr << "Gain = " << solver.gain() << std::endl;
std::cerr << "Eigenvalues = " << solver.eigenvalues() << std::endl;
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_DLQR_HPP
