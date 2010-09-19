/**
 * \file src/dcs/control/analysis/controllability.hpp
 *
 * \brief Controllability for a state-space system.
 *
 * A system with internal state vector \f$\mathbf{x}\f$ is called
 * \e controllable if and only if the system states can be changed by changing
 * the system input.
 *
 * More specifically, given a continuous-time state-space system model:
 * \f[
 *   \dot{\mathbf{x}}(t) = \mathbf{A}\mathbf{x}(t)+\mathbf{B}\mathbf{u}(t)
 *   \mathbf{y}(t) = \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f]
 * or a discrete-time state-space system model:
 * \f[
 *   \mathbf{x}(t) = \mathbf{A}\mathbf{x}(t-1)+\mathbf{B}\mathbf{u}(t-1)
 *   \mathbf{y}(t) = \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f]
 * the following definitions apply [1]:
 * - The system is called <em>state controllable</em> at time \f$t_0\f$ if
 *   it there exists an unconstrained input signal \f$\mathbf{u}(t)\f$ that
 *   willl transfer an initial state to any final state in a finite time
 *   interval \f$t_0 \le t \le t_1\f$.
 * - The system is called <em>complete state controllable</em> (or simply
 *   \e controllable) at time \f$t_0\f$ if every state \f$x_0\f$ in the
 *   state-space model is controllable.
 * - The system is called <em>complete output controllable</em> if there exists
 *   an unconstrained input signal \f$\mathbf{u}(t)\f$ that will transfer any
 *   given initial output \f$\mathbf{y}(t_0)\f$ to any final output
 *   \f$\mathbf{y}(t_1)\f$ in a finite time interval \f$t_0 \le t \le t_1\f$.
 * .
 *
 * Formally, system controllability is defined in terms of the
 * <em>controllability matrix</em>.
 *
 * Given a continuous-time state-space system model:
 * \f{align*}
 *   \dot{\mathbf{x}}(t) &= \mathbf{A}\mathbf{x}(t)+\mathbf{B}\mathbf{u}(t),\\
 *   \mathbf{y}(t) &= \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f}
 * or a discrete-time state-space system model:
 * \f{align*}
 *   \mathbf{x}(t) &= \mathbf{A}\mathbf{x}(t-1)+\mathbf{B}\mathbf{u}(t-1),\\
 *   \mathbf{y}(t) &= \mathbf{C}\mathbf{x}(t)+\mathbf{D}\mathbf{u}(t)
 * \f}
 * where:
 * - \f$\mathbf{A}\f$ is a \f$n \times n}\f$ matrix,
 * - \f$\mathbf{B}\f$ is a \f$n \times m}\f$ matrix,
 * - \f$\mathbf{C}\f$ is a \f$p \times n}\f$ matrix,
 * - \f$\mathbf{D}\f$ is a \f$p \times m}\f$ matrix.
 * .
 * We define:
 * - the <em>(state) controllability matrix</em> as the following
 *   \f$n \times nm\f$ matrix:
 *   \f[
 *     \begin{pmatrix}
 *       \mathbf{B} & \mathbf{A}\mathbf{B} & \cdots & \mathbf{A}^{n-1}\mathbf{B}
 *     \end{pmatrix}
 *   \f]
 * - the <em>output controllability matrix</em> as the following
 *   \f$p \times (n+1)m\f$ matrix:
 *   \f[
 *     \begin{pmatrix}
 *       \mathbf{C}\mathbf{B} & \mathbf{C}\mathbf{A}\mathbf{B} & \cdots & \mathbf{C}\mathbf{A}^{n-1}\mathbf{B} & \mathbf{D}
 *     \end{pmatrix}
 *   \f]
 * .
 * Then, we say that:
 * - A system is <em>complete state controllable</em> (or simply
 *   \e controllable) if and only if the associated controllability matrix is of
 *   rank \f$n\f$.
 * - A system is <em>complete output controllable</em> (or simply
 *   \e output controllable) if and only if the associated output
 *   controllability matrix is of rank \f$p\f$.
 * .
 *
 * References:
 * -# K. Ogata, <em>Modern Control Engineering</em>, 3rd edition, Prentice Hall, 2002.
 * .
 *
 * \todo
 *  It might be a good idea to add an additional parameter to the function for
 *  checking controllability, representing the tolerance with which estimating
 *  the rank of the controllability matrix. However, this may be a problem if we
 *  change in the future the way we check for controllability. As a matter of
 *  fact, there are several ways for checking the controllability of a system.
 *  For instance, in Octave, controllability is determined by applying Arnoldi
 *  iteration with complete re-orthogonalization to obtain an orthogonal basis
 *  of the Krylov subspace.
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

#ifndef DCS_CONTROL_ANALYSIS_HPP
#define DCS_CONTROL_ANALYSIS_HPP


#include <algorithm>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/rank.hpp>


namespace dcs { namespace control {

/**
 * \brief Type traits for building the type of the controllability matrix given
 *  the type of the state and the input matrices.
 *
 * \tparam AMatrixT The type of the state matrix for a state-space system model.
 * \tparam BMatrixT The type of the input matrix for a state-space system model.
 */
template <
	typename AMatrixT,
	typename BMatrixT
>
struct controllability_matrix_traits
{
	/// The type of the controllability matrix.
	typedef typename ::boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type type;
};


/**
 * \brief Type traits for building the type of the output controllability matrix
 *  given the type of the state, the input, the output, and the direct
 *  transmission matrices.
 *
 * \tparam AMatrixT The type of the state matrix for a state-space system model.
 * \tparam BMatrixT The type of the input matrix for a state-space system model.
 * \tparam CMatrixT The type of the control matrix for a state-space system
 *  model.
 * \tparam DMatrixT The type of the direct transmission matrix for a state-space
 *  system model.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT
>
struct output_controllability_matrix_traits
{
	/// The type of the output controllability matrix.
	typedef typename ::boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type type;
};


/**
 * \brief Compute the (state) controllability matrix.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \return The state controllability matrix.
 */
template <
	typename AMatrixT,
	typename BMatrixT
>
typename controllability_matrix_traits<
	AMatrixT,
	BMatrixT
>::type make_controllability_matrix(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A, ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef typename controllability_matrix_traits<AMatrixT,BMatrixT>::type matrix_type;
	typedef typename ublas::matrix_traits<matrix_type>::value_type value_type;
	typedef typename ublas::matrix_traits<matrix_type>::size_type size_type;

	size_type n = ublasx::num_rows(B);
	size_type m = ublasx::num_columns(B);

	// [B AB A^2*B ... A^{n-1}*B]
	matrix_type C = ublas::zero_matrix<value_type>(n, n*m);
	ublas::subrange(C, 0, n, 0, m) = B;
	for (size_type k = 1; k < n; ++k)
	{
		size_type km = k*m;

		// C[:,k*m+1:(k+1)*m] = A*C[:,(k-1)*m+1:k*m]
		ublas::subrange(C, 0, n, km, km+m) = ublas::prod(A, ublas::subrange(C, 0, n, km-m, km));
	}

	return C;
}


/**
 * \brief Compute the (state) controllability matrix.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \return The state controllability matrix.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT
>
typename output_controllability_matrix_traits<
	AMatrixT,
	BMatrixT,
	CMatrixT,
	DMatrixT
>::type make_output_controllability_matrix(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A, ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
										   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C, ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef typename output_controllability_matrix_traits<AMatrixT,BMatrixT,CMatrixT,DMatrixT>::type matrix_type;
	typedef typename ublas::matrix_traits<matrix_type>::value_type value_type;
	typedef typename ublas::matrix_traits<matrix_type>::size_type size_type;

	size_type n = ublasx::num_rows(B);
	size_type m = ublasx::num_columns(B);
	size_type p = ublasx::num_rows(C);
	size_type nr = p;
	size_type nc = (n+1)*m;

	// [C*B C*A*B C*A^2*B ... C*A^{n-1}*B D]
	// Compute: [B A*B A^2*B ... A^{n-1}*B]
	matrix_type P = make_controllability_matrix(A, B);
	size_type xnr = ublasx::num_rows(P);
	size_type xnc = ublasx::num_columns(P);

	// Make enough space for the output controllability matrix
	P.resize(::std::max(nr, xnr), nc, true);

	// Compute: C*[B A*B A^2*B ... A^{n-1}*B] == [C*B C*A*B C*A^{2}*B C*A^{n-1}*B]
	ublas::subrange(P, 0, nr, 0, xnc) = ublas::prod(C, ublas::subrange(P, 0, xnr, 0, xnc));

	// Add D to the last m columns
	ublas::subrange(P, 0, nr, xnc, nc) = D;

	// Make P of the right size
	P.resize(nr, nc, true);

	return P;
}


/**
 * \brief Test if a system is complete state controllable.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \return \c true if the state-space system \f$(\mathbf{A},\mathbf{B})\f$ is
 *  complete state controllable; \c false otherwise.
 */
template <
	typename AMatrixT,
	typename BMatrixT
>
bool is_controllable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A, ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
{
	return	::boost::numeric::ublasx::rank(
				make_controllability_matrix(A, B)
			)
			==
			::boost::numeric::ublasx::num_rows(B);
}


/**
 * \brief Test if a system is complete output controllable.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \tparam CMatrixT The type of the output matrix for the state-space system
 *  model.
 * \tparam DMatrixT The type of the direct transmission matrix for the
 *  state-space system  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \param C The output matrix for the state-space system model.
 * \param D The direct transmission matrix for the state-space system model.
 * \return \c true if the state-space system \f$(\mathbf{A},\mathbf{B},
 *  \mathbf{C},\mathbf{D})\f$ is complete output controllable; \c false
 *  otherwise.
 */
template <
	typename AMatrixT,
	typename BMatrixT,
	typename CMatrixT,
	typename DMatrixT
>
bool is_output_controllable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A, ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
							::boost::numeric::ublas::matrix_expression<CMatrixT> const& C, ::boost::numeric::ublas::matrix_expression<DMatrixT> const& D)
{
	return	::boost::numeric::ublasx::rank(
				make_output_controllability_matrix(A, B, C, D)
			)
			==
			::boost::numeric::ublasx::num_rows(C);
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_ANALYSIS_HPP
