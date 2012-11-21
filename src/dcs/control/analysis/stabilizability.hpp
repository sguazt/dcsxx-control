/**
 * \file src/dcs/control/analysis/stabilizability.hpp
 *
 * \brief Stabilizability analysis for linear systems.
 *
 * A system:
 * - continuous-time case: \f$\dot{x}=Ax+Bu\f$
 * - discrete-time case: \f$x(k+1)=Ax(k)+Bu(k)\f$
 * is stabilizable if all its uncontrollable modes decay to zero asymptotically.
 *
 * The stabilizability test can be performed in different ways [1]:
 * - Eigenvector test:
 *   - The continuous-time LTI system is stabilizable if and only if every
 *     eigenvector of A^T corresponding to an eigenvalue with a positive or zero
 *     real part is not in the kernel of B^T.
 *   - The discrete-time LTI system is stabilizable if and only if every
 *     eigenvector of A^T corresponding to an eigenvalue with magnitude larger
 *     or equal to 1 is not in the kernel of B^T.
 *   .
 * - Popov-Belevitch-Hautus [PBH] test:
 *   - The continuous-time LTI system is stabilizable if and only if
 *     \f[
 *       \operatorname{rank}(\begin{pmatrix}A - \lambda I & B\end{pmatrix}) = n, \quad \forall\lambda\in\mathbb{C}: \Re[\lambda] \ge 0
 *     \f]
 *   - The discrete-time LTI system is stabilizable if and only if
 *     \f[
 *       \operatorname{rank}(\begin{pmatrix}A - \lambda I & B\end{pmatrix}) = n, \quad \forall\lambda\in\mathbb{C}: |\lambda| \ge 1
 *     \f]
 *   .
 * - Lyapunov test:
 *   The LTI system (AB-LTI) is sta bilizable if and only if there is a
 *   positive-definite solution \f$P\f$ to the following Lyapunov matrix
 *   inequality:
 *   - Continuous-time case: \f$AP + PA^{T} - BB^{T} < 0\f$
 *   - Discrete-time case: \f$APA^{T} - P - B B^{T} < 0\f$
 *   .
 * .
 *
 * \note
 *  Estimating the rank of a matrix is ill-conditioned; that is, it is very
 *  sensitive to roundoff errors and errors in the data.
 *
 * \todo
 *  Implement other (more numerically robust) tests (e.g., see Octave and SCILAB).
 *
 * References:
 * -# J.P. Hesphana, "Linear Systems Theory", Princeton University Press, 2009.
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
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */

#ifndef DCS_CONTROL_ANALYSIS_STABILIZABILITY_HPP
#define DCS_CONTROL_ANALYSIS_STABILIZABILITY_HPP


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
    typename BMatrixT,
	typename RealT
>
bool is_stabilizable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				 	 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					 bool discrete,
					 RealT tol)
{
	DCS_MACRO_SUPPRESS_UNUSED_VARIABLE_WARNING( tol );

	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	// pre: A must be square
	DCS_ASSERT(
			ublasx::num_rows(A) == ublasx::num_columns(A),
			throw ::std::invalid_argument("[dcs::control::is_stabilizable] The A matrix must be square.")
		);
	// pre: num_rows(A) == num_rows(B)
	DCS_ASSERT(
			ublasx::num_rows(A) == ublasx::num_rows(B),
			throw ::std::invalid_argument("[dcs::control::is_stabilizable] The size of A and B matrices are not conformant.")
		);

	typedef typename ublas::promote_traits<typename AMatrixT::value_type,
										   typename BMatrixT::value_type>::promote_type value_type;
	typedef typename ublas::promote_traits<typename AMatrixT::size_type,
										   typename BMatrixT::size_type>::promote_type size_type;
	typedef typename ublas::type_traits<value_type>::real_type real_type;
	typedef ::std::complex<real_type> complex_type;


	// Use the Popov-Belevitch-Hautus [PBH] test for stabilizability.
	// E.g., see Hespanha, "Linear System Theory", 2009, pp.125-126.

	ublas::vector<complex_type> v;
	ublasx::eigenvalues(A, v);
	size_type nv(ublasx::size(v));
	size_type na(ublas::num_rows(A));
	size_type nb(ublas::num_columns(B));
	size_type ns(na+nb);
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
		ublas::matrix<complex_type> S(na,ns);
		ublas::subrange(S, 0, na, 0, na) = cA-v(i)*I;
		ublas::subrange(S, 0, na, na, ns) = B;

		if (ublasx::rank(S) != na)
		{
			return false;
		}
	}

	return true;
/*
-- octave
	// Controllability staircase form
	[ac, ~, ~, ncont] = slab01od(A, B, tol);

	// Extract uncontrollable part of staircase form
    uncont_idx = ncont+1 : rows (a);
    auncont = ac(uncont_idx, uncont_idx);

    ## calculate poles of uncontrollable part
    eigw = eig (auncont);

  ## check whether uncontrollable poles are stable
  bool = __is_stable__ (eigw, ! dflg, tol);

-- scilab
    [n,u,ind,V,a,b]=contr(a,b,tol);
  n=sum(n);nc=n;
  if lhs==4 then c=c*u;x0=u'*x0;end
  if n<>na then
    //order evals uncont. part
    nn=n+1:na;
    [v,n1]=schur(a(nn,nn),part(typ,1))
    n=n+n1
    //new realization
    if lhs>2 then
    u(:,nn)=u(:,nn)*v
      if lhs==4 then
    a(:,nn)=a(:,nn)*v;a(nn,nn)=v'*a(nn,nn)
    b(nn,:)=v'*b(nn,:)
    c(:,nn)=c(:,nn)*v
    x0(nn)=v'*x0(nn)
      end;
    end;
  end;
  if lhs==4 then sl=syslin(dom,a,b,c,d,x0),end
  if lhs==5 then v=sl.B;sl=sl.A;end
*/

}


template <
    typename AMatrixT,
    typename BMatrixT
>
inline
bool is_stabilizable(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				 	 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					 bool discrete)
{
	return is_stabilizable(A, B, discrete, 0.0);
}

}} // Namespace dcs::control


#endif // DCS_CONTROL_ANALYSIS_STABILIZABILITY_HPP
