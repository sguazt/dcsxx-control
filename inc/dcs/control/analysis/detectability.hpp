/**
 * \file dcs/control/analysis/detectability.hpp
 *
 * \brief Detectability analysis for linear systems.
 *
 * A system:
 * - continuous-time case: \f$\dot{x}=Ax,\quad y=Cxf$
 * - discrete-time case: \f$x(k+1)=Ax(k)\quad y(k)=Cx(k)\f$
 * is detectable if all its unobservable modes decay to zero asymptotically.
 *
 * The detectability test can be performed in different ways:
 * - Eigenvector test [1]:
 *   - The continuous-time LTI system is detectable if and only if every
 *     eigenvector of A^T corresponding to an eigenvalue with a positive or zero
 *     real part is not in the kernel of C.
 *   - The discrete-time LTI system is detectable if and only if every
 *     eigenvector of A^T corresponding to an eigenvalue with magnitude larger
 *     or equal to 1 is not in the kernel of C.
 *   .
 * - Popov-Belevitch-Hautus [PBH] test [1,2]:
 *   - The continuous-time LTI system is detectable if and only if
 *     \f[
 *       \operatorname{rank}(\begin{pmatrix}A - \lambda I \\ C\end{pmatrix}) = n, \quad \forall\lambda\in\mathbb{C}: \Re[\lambda] \ge 0
 *     \f]
 *   - The discrete-time LTI system is detectable if and only if
 *     \f[
 *       \operatorname{rank}(\begin{pmatrix}A - \lambda I \\ C\end{pmatrix}) = n, \quad \forall\lambda\in\mathbb{C}: |\lambda| \ge 1
 *     \f]
 *   .
 * - Lyapunov test [1]:
 *   The LTI system (AB-LTI) is detectable if and only if there is a
 *   positive-definite solution \f$P\f$ to the following Lyapunov matrix
 *   inequality:
 *   - Continuous-time case: \f$A^{T}P + PA - C^{T}C < 0\f$
 *   - Discrete-time case: \f$A^{T}PA - P - C^{T}C < 0\f$
 *   .
 * - Stabilization test [2]:
 *   The matrix pair \f$(C,A)\f$ is detectable if and only if the dual par \f$(A^T,C^T)\f$ is stabilizable.
 * .
 *
 * \note
 *  Estimating the rank of a matrix is ill-conditioned; that is, it is very
 *  sensitive to roundoff errors and errors in the data.
 *
 * \todo
 *  Implement other (more numerically robust) tests (e.g., see Octave and SCILAB).
 *
 *
 * References:
 * -# J.P. Hesphana, "Linear Systems Theory", Princeton University Press, 2009.
 * -# W.J. Terrel, "Stability and Stabilization: An Introduction", Princeton University Press, 2009
 * .
 *
 * <hr/>
 *
 * Copyright (C) 2011  Distributed Computing System (DCS) Group, Computer
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
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 */

#ifndef DCS_CONTROL_ANALYSIS_DETECTABILITY_HPP
#define DCS_CONTROL_ANALYSIS_DETECTABILITY_HPP


#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublasx/operation/eigen.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/rank.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <cmath>
#include <complex>
#include <dcs/assert.hpp>
#include <dcs/macro.hpp>
#include <stdexcept>


namespace dcs { namespace control {

template <
    typename AMatrixT,
    typename CMatrixT,
	typename RealT
>
bool is_detectable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				   ::boost::numeric::ublas::matrix_expression<CMatrixT> const& C,
				   bool discrete,
				   RealT tol)
{
	DCS_MACRO_SUPPRESS_UNUSED_VARIABLE_WARNING( tol );

	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	// pre: A must be square
	DCS_ASSERT(
			ublasx::num_rows(A) == ublasx::num_columns(A),
			throw ::std::invalid_argument("[dcs::control::is_detectable] The A matrix must be square.")
		);
	// pre: num_rows(A) == num_rows(B)
	DCS_ASSERT(
			ublasx::num_rows(A) == ublasx::num_columns(C),
			throw ::std::invalid_argument("[dcs::control::is_detectable] The size of A and C matrices are not conformant.")
		);

	typedef typename ublas::promote_traits<typename AMatrixT::value_type,
										   typename CMatrixT::value_type>::promote_type value_type;
	typedef typename ublas::promote_traits<typename AMatrixT::size_type,
										   typename CMatrixT::size_type>::promote_type size_type;
	typedef typename ublas::type_traits<value_type>::real_type real_type;
	typedef ::std::complex<real_type> complex_type;


	// Use the Popov-Belevitch-Hautus [PBH] test for detectability.
	// E.g., see Hespanha, "Linear System Theory", 2009, pp.152-153.

	ublas::vector<complex_type> v;
	ublasx::eigenvalues(A, v);
	size_type nv(ublasx::size(v));
	size_type na(ublas::num_rows(A));
	size_type nc(ublas::num_rows(C));
	size_type ns(na+nc);
	ublas::matrix<complex_type> cA(A);
	ublas::identity_matrix<complex_type> I(na);
	for (size_type i = 0; i < nv; ++i)
	{
		if (discrete)
		{
			if (::std::abs(v(i)) < 1)
			{
				continue;
			}
		}
		else
		{
			if (v(i).real() < 0)
			{
				continue;
			}
		}

		//ublas::matrix<complex_type> AA(ublas::matrix<complex_type>(A)-v(i)*I);
		ublas::matrix<complex_type> S(ns,na);
		ublas::subrange(S, 0, na, 0, na) = cA-v(i)*I;
		ublas::subrange(S, na, ns, 0, na) = C;

		if (ublasx::rank(S) != na)
		{
			return false;
		}
	}

	return true;
}


template <
    typename AMatrixT,
    typename BMatrixT
>
inline
bool is_detectable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				 	 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					 bool discrete)
{
	return is_detectable(A, B, discrete, 0.0);
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_ANALYSIS_DETECTABILITY_HPP
