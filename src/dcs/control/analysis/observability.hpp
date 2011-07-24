/**
 * \file src/dcs/control/analysis/observability.hpp
 *
 * \brief Observability for a state-space system.
 *
 * A system is said to be observable if any initial state \f$x(k_0)\f$ can be
 * estimated from the control sequence \f$u(k)\f$, with
 * \f$k=k_0,k_0+1,\ldots,k_f-1\f$, and the measurements \f$y(k)\f$, with
 * \f$k=k_0,k_0+1,\ldots,k_f\f$.
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
 * - The system is said to be <em>completely observable</em> if every state
 *   \f$\mathbf{x}(t_0)\f$ can be determined from the observation of the output
 * .
 *   vector \f$\mathbf{y}(t)\f$ over a finite time interval \f$f_0 \le t \le t_1\f$.
 * The system is therefore completely observable if every transition of the
 * state eventually affects every element of the output vectory.
 * The concept of observability is useful in solving the problem of
 * reconstructing unmeasurable state variables from measurable variables in the
 * minimun possible length of time.
 *
 * Formally, system observability is defined in terms of the
 * <em>observability matrix</em>.
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
 * - \f$\mathbf{B}\f$ is a \f$n \times p}\f$ matrix,
 * - \f$\mathbf{C}\f$ is a \f$m \times n}\f$ matrix,
 * - \f$\mathbf{D}\f$ is a \f$m \times p}\f$ matrix.
 * .
 * We define the <em>observability matrix</em> as the following
 * \f$nm \times n\f$ matrix:
 * \f[
 *   \begin{pmatrix}
 *     \mathbf{C} \\
 *     \mathbf{C}\mathbf{A} \\
 *     \cdots \\
 *     \mathbf{B}\mathbf{A}^{n-1}
 *   \end{pmatrix}
 * \f]
 * Then, we say that a system is <em>complete observable</em> (or simply
 *   \e observable) if and only if the associated observability matrix is of
 *   rank \f$n\f$.
 *
 * A weaker concept of observability is \e detectability (a concept which is
 * dual to the one of stabilizability): for a partially observable system, if
 * the unobservable modes are stable and the observable modes are unstable, the
 * system is said to be detectable.
 *
 * \note
 *  Computing the rank of the observability matrix is not recommended for
 *  observability testing.
 *  Indeed, the resulting observability matrix will be numerically singular for
 *  most systems with more than a handful of states.
 *  This fact is well documented in the control literature.
 *  For example, see section III in [2].
 *
 * References:
 * -# K. Ogata, <em>Modern Control Engineering</em>, 4th edition, Prentice Hall, 2010.
 * -# C.C. Paige, <em>Properties of Numerival Algorithms Related to Computing Controllability</em>, IEEE Transactions on Automatic Control, 26(1), 1981.
 *    [http://lawww.epfl.ch/webdav/site/la/users/105941/public/NumCompCtrl.pdf]
 * .
 *
 * \todo
 *  It might be a good idea to add an additional parameter to the function for
 *  checking observability, representing the tolerance with which estimating
 *  the rank of the observability matrix. However, this may be a problem if we
 *  change in the future the way we check for observability. As a matter of
 *  fact, there are several ways for checking the observability of a system.
 *  For instance, in Octave, observability is determined by applying Arnoldi
 *  iteration with complete re-orthogonalization to obtain an orthogonal basis
 *  of the Krylov subspace.
 *
 * <hr/>
 *
 * Copyright (C) 2009-2011  Distributed Computing System (DCS) Group, Computer
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

#ifndef DCS_CONTROL_ANALYSIS_OBSERVABILITY_HPP
#define DCS_CONTROL_ANALYSIS_OBSERVABILITY_HPP


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
 * \brief Type traits for building the type of the observability matrix given
 *  the type of the state and the input matrices.
 *
 * \tparam AMatrixT The type of the state matrix for a state-space system model.
 * \tparam BMatrixT The type of the input matrix for a state-space system model.
 */
template <
	typename AMatrixT,
	typename CMatrixT
>
struct observability_matrix_traits
{
	/// The type of the observability matrix.
	typedef typename ::boost::numeric::ublas::matrix_temporary_traits<AMatrixT>::type type;
};


/**
 * \brief Compute the observability matrix.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam BMatrixT The type of the input matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param B The input matrix for the state-space system model.
 * \return The observability matrix.
 */
template <
	typename AMatrixT,
	typename CMatrixT
>
typename observability_matrix_traits<
	AMatrixT,
	CMatrixT
>::type make_observability_matrix(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
								  ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef typename observability_matrix_traits<AMatrixT,CMatrixT>::type matrix_type;
	typedef typename ublas::matrix_traits<matrix_type>::value_type value_type;
	typedef typename ublas::matrix_traits<matrix_type>::size_type size_type;

	size_type m = ublasx::num_rows(C);
	size_type n = ublasx::num_columns(C);

	// [C C*A C*A^2 ... C*A^{n-1}]
	matrix_type O = ublas::zero_matrix<value_type>(n*m, n);
	ublas::subrange(O, 0, m, 0, n) = C;
	for (size_type k = 1; k < n; ++k)
	{
		size_type km = k*m;

		// O[k*m+1:(k+1)*m,:] = A*O[(k-1)*m+1:k*m,:]
		ublas::subrange(O, km, km+m, 0, n) = ublas::prod(ublas::subrange(O, km-m, km, 0, n), A);
	}

	return O;
}


/**
 * \brief Test if a system is complete observable.
 *
 * \tparam AMatrixT The type of the state matrix for the state-space system
 *  model.
 * \tparam CMatrixT The type of the output matrix for the state-space system
 *  model.
 * \param A The state matrix for the state-space system model.
 * \param C The output matrix for the state-space system model.
 * \return \c true if the state-space system \f$(\mathbf{A},\mathbf{C})\f$ is
 *  complete observable; \c false otherwise.
 */
template <
	typename AMatrixT,
	typename CMatrixT
>
bool is_observable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C)
{
	return	::boost::numeric::ublasx::rank(
				make_observability_matrix(A, C)
			)
			==
			::boost::numeric::ublasx::num_columns(C);
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_ANALYSIS_OBSERVABILITY_HPP
