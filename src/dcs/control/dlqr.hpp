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


#include <dcs/debug.hpp>
#include <dcs/math/la/container/dense_matrix.hpp>
#include <dcs/math/la/matrix_basic_operations.hpp>
#include <dcs/math/la/operation/cat.hpp>
#include <dcs/math/la/operation/eigen.hpp>
#include <dcs/math/la/operation/min.hpp>
#include <dcs/math/la/operation/num_columns.hpp>
#include <dcs/math/la/operation/num_rows.hpp>
#include <dcs/math/la/traits/matrix.hpp>
#include <dcs/math/la/traits/promote.hpp>
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
void dlqr(Matrix1T const& A, Matrix2T& const B, Matrix3T& const Q, Matrix4T& const R, Matrix5T& const N)
{
	namespace dcs_la = ::dcs::math::la;

	typedef Matrix1T A_matrix1_type;
	typedef Matrix2T B_matrix2_type;
	typedef Matrix3T Q_matrix3_type;
	typedef Matrix4T R_matrix4_type;
	typedef Matrix5T N_matrix5_type;
	typedef typename dcs_la::promote_traits<
				typename dcs_la::matrix_traitsA_<matrix_type>::value_type,
				typename dcs_la::promote_traits<
					typename dcs_la::matrix_traits<B_matrix_type>::value_type,
					typename dcs_la::promote_traits<
						typename dcs_la::matrix_traits<Q_matrix_type>::value_type,
						typename dcs_la::promote_traits<
							typename dcs_la::matrix_traits<R_matrix_type>::value_type,
							typename dcs_la::matrix_traits<N_matrix_type>::value_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type value_type;
	typedef typename dcs_la::promote_traits<
				typename dcs_la::matrix_traits<A_matrix_type>::size_type,
				typename dcs_la::promote_traits<
					typename dcs_la::matrix_traits<B_matrix_type>::size_type,
					typename dcs_la::promote_traits<
						typename dcs_la::matrix_traits<Q_matrix_type>::size_type,
						typename dcs_la::promote_traits<
							typename dcs_la::matrix_traits<R_matrix_type>::size_type,
							typename dcs_la::matrix_traits<N_matrix_type>::size_type
						>::promote_type
					>::promote_type
				>::promote_type
			>::promote_type size_type;
	typedef dcs_la::dense_matrix<value_type> tmp_matrix_type;

	size_type A_nr = dcs_la::num_rows(A);
	size_type A_nc = dcs_la::num_columns(A);
	size_type B_nr = dcs_la::num_rows(B);
	size_type B_nc = dcs_la::num_columns(B);
	size_type Q_nr = dcs_la::num_rows(Q);
	size_type Q_nc = dcs_la::num_columns(Q);
	size_type R_nr = dcs_la::num_rows(Q);
	size_type R_nc = dcs_la::num_columns(Q);
	size_type N_nr = dcs_la::num_rows(S);
	size_type N_nc = dcs_la::num_columns(S);

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

	if (dcs_la::norm_1(dcs_la::trans(Q)-Q) > (100*eps*dcs_la::norm_1(Q)))
	{
		DCS_DEBUG_TRACE("The error weighted matrix Q is not symmetric and has been replaced by (Q+Q')/2.");
		tmp_Q = (Q+dcs_la::trans(Q))/static_cast<value_type>(2);
	}
	else
	{
		tmp_Q = Q;
	}
	if (dcs_la::norm_1(dcs_la::trans(R)-R) > (100*eps*dcs_la::norm_1(R)))
	{
		DCS_DEBUG_TRACE("The control weighted matrix R is not symmetric and has been replaced by (R+R')/2.");
		tmp_R = (R+dcs_la::trans(R))/static_cast<value_type>(2);
	}
	else
	{
		tmp_R = R;
	}


	v_R = dcs_la::eigenvalues(tmp_R);
	v_QNR = dcs_la::real(
				dcs_la::eigenvalues(
					dcs_la::cat_rows<tmp_matrix_type>(
						dcs_la::cat_columns<tmp_matrix_type>(tmp_Q, NN), 
						dcs_la::cat_columns<tmp_matrix_type>(dcs_la::trans(NN), tmp_R) 
					)
				)
		);
	if (dcs_la::min(v_R) <= 0)
	{
		throw ::std::runtime_error("[dcs::control::dlqr] The control weighted matrix R is not positive definite.");
	}
	else if (dcs_la::min(v_QNR) < (1e2*eps*::std::max(0, dcs_la::max(v_QNR))))
	{
		DCS_DEBUG_TRACE("The matrix [Q N;N' R] is not positive semi-definite.");
	}

	dare(A, B, Q, R, NN, E);
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_DLQR_HPP
