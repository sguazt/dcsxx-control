/**
 * \file dcs/control/bindings/slicot/sb02od.hpp
 *
 * \brief Binding for the SLICOT SB02OD Fortran routine.
 *
 * The purpose of the SB02OD routine is to solve for \f$X\f$ either the
 * continuous-time algebraic Riccati equation
 * \f[
 *  Q + A'X + XA - (L+XB)R^{-1}(L+XB)' = 0
 * \f]
 * or the discrete-time algebraic Riccati equation
 * \f[
 *  X = A'XA - (L+A'XB)(R + B'XB)^{-1}(L+A'XB)' + Q
 * \f]
 * where \f$A\f$, \f$B\f$, \f$Q\f$, \f$R\f$, and \f$L\f$ are
 * \f$N\f$-by-\f$N\f$, \f$N\f$-by-\f$M\f$, \f$N\f$-by-\f$N\f$,
 * \f$M\f$-by-\f$M\f$ and \f$N\f$-by-\f$M\f$ matrices,
 * respectively, such that \f$Q = C'C\f$, \f$R = D'D\f$ and \f$L = C'D\f$;
 * \f$X\f$ is an \f$N\f$-by-\f$N\f$ symmetric matrix.
 *
 * The routine also returns the computed values of the closed-loop
 * spectrum of the system, i.e., the stable eigenvalues
 * \f$\lambda(1),\ldots,\lambda(N)\f$ of the corresponding Hamiltonian or
 * symplectic pencil, in the continuous-time case or discrte-time case,
 * respectively.
 *
 * Optionally, matrix \f$G = BR^{-1}B'\f$ may be given instead of \f$B\f$ and
 * \f$R\f$.
 * Other options include the case with \f$Q\f$ and/or \f$R\f$ given in a
 * factored form, \f$Q = C'C\f$, \f$R = D'D\f$, and with L a zero matrix.
 *
 * The routine uses the method of deflating subspaces, based on
 * reordering the eigenvalues in a generalized Schur matrix pair.
 * A standard eigenproblem is solved in the continuous-time case if \f$G\f$ is
 * given.
 *
 * The routine uses a variant of the method of deflating subspaces
 * proposed by van Dooren [1]. See also [2], [3].
 * It is assumed that \f$(A,B)\f$ is stabilizable and \f$(C,A)\f$ is detectable.
 * Under these assumptions the algebraic Riccati equation is known to
 * have a unique non-negative definite solution.
 *
 * REFERENCES
 * -# Van Dooren, P.
 *    A Generalized Eigenvalue Approach for Solving Riccati Equations.
 *    SIAM J. Sci. Stat. Comp., 2, pp. 121-135, 1981.
 * -# Mehrmann, V.
 *    The Autonomous Linear Quadratic Control Problem. Theory and Numerical
 *    Solution.
 *    Lecture Notes in Control and Information Sciences, vol. 163,
 *    Springer-Verlag, Berlin, 1991.
 * -# Sima, V.
 *    Algorithms for Linear-Quadratic Optimization.
 *    Pure and Applied Mathematics: A Series of Monographs and
 *    Textbooks, vol. 200, Marcel Dekker, Inc., New York, 1996.
 * .
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright 2012 Marco Guazzone (marco.guazzone@gmail.com)
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

#ifndef DCS_CONTROL_BINDINGS_SLICOT_SB02OD_HPP
#define DCS_CONTROL_BINDINGS_SLICOT_SB02OD_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation/num_columns.hpp>
#include <boost/numeric/ublas/operation/num_rows.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <dcs/debug.hpp>
#include <dcs/exception.hpp>
#include <dcs/control/bindings/fortran.hpp>
#include <dcs/control/bindings/slicot/detail/slicot.hpp>
#include <stdexcept>


namespace dcs { namespace control { namespace bindings { namespace slicot {

inline
void sb02od(fortran_char const& dico,
			fortran_char const& jobb,
			fortran_char const& fact,
			fortran_char const& uplo,
			fortran_char const& jobl,
			fortran_char const& sort,
			fortran_int const& n,
			fortran_int const& m,
			fortran_int const& p,
			fortran_double_real const* a,
			fortran_int const& lda,
			fortran_double_real const* b,
			fortran_int const& ldb,
			fortran_double_real const* q,
			fortran_int const& ldq,
			fortran_double_real const* r,
			fortran_int const& ldr,
			fortran_double_real const* l,
			fortran_int const& ldl,
			fortran_double_real& rcond,
			fortran_double_real* x,
			fortran_int const& ldx,
			fortran_double_real* alfar,
			fortran_double_real* alfai,
			fortran_double_real* beta,
			fortran_double_real* s,
			fortran_int const& lds,
			fortran_double_real* t,
			fortran_int const& ldt,
			fortran_double_real* u,
			fortran_int const& ldu,
			fortran_double_real const& tol,
			fortran_int* iwork,
			fortran_double_real* dwork,
			fortran_int const& ldwork,
			fortran_logical* bwork)
{
	fortran_int info(0);

	DCS_CONTROL_BINDINGS_SLICOT_SB02OD(&dico,
									   &jobb,
									   &fact,
									   &uplo,
									   &jobl,
									   &sort,
									   &n,
									   &m,
									   &p,
									   a,
									   &lda,
									   b,
									   &ldb,
									   q,
									   &ldq,
									   r,
									   &ldr,
									   l,
									   &ldl,
									   &rcond,
									   x,
									   &ldx,
									   alfar,
									   alfai,
									   beta,
									   s,
									   &lds,
									   t,
									   &ldt,
									   u,
									   &ldu,
									   &tol,
									   iwork,
									   dwork,
									   &ldwork,
									   bwork,
									   &info);

	if (info != 0)
	{
		::std::ostringstream oss;

		switch (info)
		{
			case 1:
				oss << "The computed extended matrix pencil is singular, possibly due to rounding errors";
				break;
			case 2:
				oss << "The QZ/QR algorithm failed";
				break;
			case 3:
				oss << "Reordering of the (generalized) eigenvalues failed";
				break;
			case 4:
				oss << "After reordering, roundoff changed values of some complex eigenvalues so that leading eigenvalues in the (generalized) Schur form no longer satisfy the stability condition; this could also be caused due to scaling";
				break;
			case 5:
				oss << "The computed dimension of the solution does not equal the number of states";
				break;
			case 6:
				oss << "A singular matrix was encountered during the computation of the solution matrix";
				break;
			default:
				// OK
				break;
		}

		DCS_EXCEPTION_THROW(::std::runtime_error, oss.str());
	}
}

namespace detail { namespace /*<unnamed>*/ {

template <typename AMatrixT,
		  typename BMatrixT,
		  typename QMatrixT,
		  typename RMatrixT,
		  typename LMatrixT,
		  typename XMatrixT,
		  typename LambdaVectorT,
		  typename SMatrixT,
		  typename TMatrixT,
		  typename UMatrixT>
void sb02od_impl(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
				 ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
				 ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
				 ::boost::numeric::ublas::matrix_expression<LMatrixT> const& L,
				 bool discrete,
				 ::boost::numeric::ublas::matrix_container<XMatrixT>& X,
				 ::boost::numeric::ublas::vector_container<LambdaVectorT>& lambda,
				 ::boost::numeric::ublas::matrix_container<SMatrixT>& S,
				 ::boost::numeric::ublas::matrix_container<TMatrixT>& T,
				 ::boost::numeric::ublas::matrix_container<UMatrixT>& U,
				 ::boost::numeric::ublas::column_major_tag)
{
	namespace ublas = ::boost::numeric::ublas;

	typedef ublas::matrix<fortran_double_real, ublas::column_major> matrix_type;
	typedef ublas::vector<fortran_double_real> real_vector_type;
	typedef ublas::vector<fortran_int> int_vector_type;
	typedef ublas::vector<fortran_logical> logical_vector_type;

	const fortran_char dico(discrete ? 'D' : 'C');
	const fortran_char jobb('B');
	const fortran_char fact('N');
	const fortran_char uplo('U');
	const fortran_char jobl((ublas::num_rows(L) > 0 && ublas::num_columns(L) > 0) ? 'N' : 'Z');
	const fortran_char sort('S');

	const fortran_int n(ublas::num_rows(A)); // Number of states
	const fortran_int m(ublas::num_columns(B)); // Number of inputs
	const fortran_int p(0); // Number of outputs (not used because FACT == 'N')
	const fortran_int lda(::std::max(1, n));
	const fortran_int ldb(::std::max(1, n));
	const fortran_int ldq(::std::max(1, n));
	const fortran_int ldr(::std::max(1, m));
	const fortran_int ldl((jobl == 'Z') ? 1 : ::std::max(1, n));

	matrix_type tmp_L;
	if (jobl == 'Z')
	{
		tmp_L = matrix_type(ldl,m,0);
	}
	else
	{
		tmp_L = L;
		tmp_L.resize(ldl, m, false);
	}

	const fortran_int ldx(::std::max(1, n));
	X().resize(ldx, n, false);
	const fortran_int n2(2*n);
	real_vector_type tmp_alfar(n2,0);
	real_vector_type tmp_alfai(n2,0);
	real_vector_type tmp_beta(n2, 0);
	const fortran_int lds((jobb == 'B') ? ::std::max(1, n2+m) : ::std::max(1,n2));
	S().resize(lds, (jobb == 'B') ? (n2+m) : n2, false);
	const fortran_int ldt((jobb == 'G') ? ((dico == 'D') ? ::std::max(1,n2) : 1) : ::std::max(1, n2+m));
	T().resize(ldt, n2, false);
	const fortran_int ldu(::std::max(1, n2));
	U().resize(ldu, n2, false);

	const fortran_double_real tol(0); // use default value (i.e., eps - the machine precision)
	fortran_double_real rcond(0);

	const fortran_int liwork((jobb == 'B') ? ::std::max(::std::max(1, m), n2) : ::std::max(1,n2));
	int_vector_type iwork(liwork,0);
	const fortran_int ldwork((jobb == 'G') ? ((dico == 'C') ? ::std::max(3,6*n) : ::std::max(7*(n2+1)+16,16*n)) : ::std::max(::std::max(::std::max(7*(n2+1)+16, 16*n), n2+m), 3*m));
	real_vector_type dwork(ldwork,0);
	logical_vector_type bwork(n2,0);

	sb02od(dico,
		   jobb,
		   fact,
		   uplo,
		   jobl,
		   sort,
		   n,
		   m,
		   p,
		   A().data().begin(),
		   lda,
		   B().data().begin(),
		   ldb,
		   Q().data().begin(),
		   ldq,
		   R().data().begin(),
		   ldr,
		   tmp_L.data().begin(),
		   ldl,
		   rcond,
		   X().data().begin(),
		   ldx,
		   tmp_alfar.data().begin(),
		   tmp_alfai.data().begin(),
		   tmp_beta.data().begin(),
		   S().data().begin(),
		   lds,
		   T().data().begin(),
		   ldt,
		   U().data().begin(),
		   ldu,
		   tol,
		   iwork.data().begin(),
		   dwork.data().begin(),
		   ldwork,
		   bwork.data().begin());

	// Fill the lambda vector:
	// \lambda(k) = \frac{\alfa_r(k)+j\alfa_i(k)}{\beta(k)}
	lambda().resize(n2, false);
	typedef typename ublas::vector_traits<LambdaVectorT>::value_type complex_type;
	for (fortran_int k = 0; k < n2; ++k)
	{
		lambda()(k) = complex_type(tmp_alfar(k), tmp_alfai(k))/tmp_beta(k);
	}
}

template <typename AMatrixT,
		  typename BMatrixT,
		  typename QMatrixT,
		  typename RMatrixT,
		  typename LMatrixT,
		  typename XMatrixT,
		  typename LambdaVectorT,
		  typename SMatrixT,
		  typename TMatrixT,
		  typename UMatrixT>
inline
void sb02od_impl(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
				 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
				 ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
				 ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
				 ::boost::numeric::ublas::matrix_expression<LMatrixT> const& L,
				 bool discrete,
				 ::boost::numeric::ublas::matrix_container<XMatrixT>& X,
				 ::boost::numeric::ublas::vector_container<LambdaVectorT>& lambda,
				 ::boost::numeric::ublas::matrix_container<SMatrixT>& S,
				 ::boost::numeric::ublas::matrix_container<TMatrixT>& T,
				 ::boost::numeric::ublas::matrix_container<UMatrixT>& U,
				 ::boost::numeric::ublas::row_major_tag)
{
	namespace ublas = ::boost::numeric::ublas;

	typedef ublas::matrix<typename ublas::matrix_traits<AMatrixT>::value_type,ublas::column_major> A_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<BMatrixT>::value_type,ublas::column_major> B_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<QMatrixT>::value_type,ublas::column_major> Q_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<RMatrixT>::value_type,ublas::column_major> R_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<LMatrixT>::value_type,ublas::column_major> L_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<XMatrixT>::value_type,ublas::column_major> X_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<SMatrixT>::value_type,ublas::column_major> S_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<TMatrixT>::value_type,ublas::column_major> T_colmaj_matrix_type;
	typedef ublas::matrix<typename ublas::matrix_traits<UMatrixT>::value_type,ublas::column_major> U_colmaj_matrix_type;

	A_colmaj_matrix_type tmp_A(A);
	B_colmaj_matrix_type tmp_B(B);
	Q_colmaj_matrix_type tmp_Q(Q);
	R_colmaj_matrix_type tmp_R(R);
	L_colmaj_matrix_type tmp_L(L);
	X_colmaj_matrix_type tmp_X;
	S_colmaj_matrix_type tmp_S;
	T_colmaj_matrix_type tmp_T;
	U_colmaj_matrix_type tmp_U;

	sb02od_impl(tmp_A,
				tmp_B,
				tmp_Q,
				tmp_R,
				tmp_L,
				discrete,
				tmp_X,
				lambda,
				tmp_S,
				tmp_T,
				tmp_U,
				ublas::column_major_tag());

	X() = tmp_X;
	S() = tmp_S;
	T() = tmp_T;
	U() = tmp_U;
}

}} // Namespace detail::<unnamed>

template <typename AMatrixT,
		  typename BMatrixT,
		  typename QMatrixT,
		  typename RMatrixT,
		  typename LMatrixT,
		  typename XMatrixT,
		  typename LambdaVectorT,
		  typename SMatrixT,
		  typename TMatrixT,
		  typename UMatrixT>
inline
void sb02od(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
			::boost::numeric::ublas::matrix_expression<LMatrixT> const& L,
			bool discrete,
			::boost::numeric::ublas::matrix_container<XMatrixT>& X,
			::boost::numeric::ublas::vector_container<LambdaVectorT>& lambda,
			::boost::numeric::ublas::matrix_container<SMatrixT>& S,
			::boost::numeric::ublas::matrix_container<TMatrixT>& T,
			::boost::numeric::ublas::matrix_container<UMatrixT>& U)
{
	namespace ublas = ::boost::numeric::ublas;

	typedef typename ublas::matrix_traits<AMatrixT>::orientation_category A_orientation_category;
	typedef typename ublas::matrix_traits<BMatrixT>::orientation_category B_orientation_category;
	typedef typename ublas::matrix_traits<QMatrixT>::orientation_category Q_orientation_category;
	typedef typename ublas::matrix_traits<RMatrixT>::orientation_category R_orientation_category;
	typedef typename ublas::matrix_traits<LMatrixT>::orientation_category L_orientation_category;
	typedef typename ublas::matrix_traits<XMatrixT>::orientation_category X_orientation_category;
	typedef typename ublas::matrix_traits<SMatrixT>::orientation_category S_orientation_category;
	typedef typename ublas::matrix_traits<TMatrixT>::orientation_category T_orientation_category;
	typedef typename ublas::matrix_traits<UMatrixT>::orientation_category U_orientation_category;

	// pre: same orientation category
	BOOST_STATIC_ASSERT((
		::boost::mpl::and_<
			::boost::is_same<A_orientation_category,B_orientation_category>,
			::boost::mpl::and_<
				::boost::is_same<A_orientation_category,Q_orientation_category>,
				::boost::mpl::and_<
					::boost::is_same<A_orientation_category,R_orientation_category>,
					::boost::mpl::and_<
						::boost::is_same<A_orientation_category,L_orientation_category>,
						::boost::mpl::and_<
							::boost::is_same<A_orientation_category,X_orientation_category>,
							::boost::mpl::and_<
								::boost::is_same<A_orientation_category,S_orientation_category>,
								::boost::mpl::and_<
									::boost::is_same<A_orientation_category,T_orientation_category>,
									::boost::is_same<A_orientation_category,U_orientation_category>
								>
							>
						>
					>
				>
			>
		>::value
	));

	// Dispatch to the right implementation
	detail::sb02od_impl(A,
						B,
						Q,
						R,
						L,
						discrete,
						X,
						lambda,
						S,
						T,
						U,
						A_orientation_category());
}

}}}} // Namespace dcs::control::bindings::slicot

#endif // DCS_CONTROL_BINDINGS_SLICOT_SB02OD_HPP
