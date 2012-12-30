/**
 * \file dcs/control/design/fmpc.hpp
 *
 * \brief Fast Model Predictive Control using Online Optimization.
 *
 *
 * <hr/>
 *
 * Copyright (C) 2009-2012  Distributed Computing System (DCS) Group, Computer
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

#ifndef DCS_CONTROL_FMPC_HPP
#define DCS_CONTROL_FMPC_HPP


#include <algorithm>
#include <boost/static_assert.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/lapack.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/views.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublasx/operation/cat.hpp>
#include <boost/numeric/ublasx/operation/num_columns.hpp>
#include <boost/numeric/ublasx/operation/num_rows.hpp>
#include <boost/numeric/ublasx/operation/size.hpp>
#include <boost/numeric/ublasx/traits/layout_type.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <cmath>
#include <complex>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
//#include <dcs/control/detail/blas.hpp>
//#include <dcs/control/detail/lapack.hpp>
#include <dcs/macro.hpp>
#include <stdexcept>


namespace dcs { namespace control {

namespace detail { namespace /*<unnamed>*/ {


const int ione = 1;
const int itwo = 2;
const int ithree = 3;
const int iseven = 7;
const double fone = 1;
const double ftwo = 2;
const double fzero = 0;
const double fmone = -1;
const int quiet = 0;

///@{ Prototypes

template <typename RealT, typename SizeT>
void dnudz(RealT* A, RealT* B, RealT* At, RealT* Bt, RealT* eyen, RealT* eyem, RealT* Q, RealT* R, RealT* Qf, RealT* hp, RealT* rd, RealT* rp, SizeT T, SizeT n, SizeT m, SizeT nz, RealT kappa, RealT* dnu, RealT* dz);

template <typename RealT, typename SizeT>
void rdrp(RealT* A, RealT* B, RealT* z, RealT* nu, RealT* gf, RealT* gp, RealT* b, SizeT T, SizeT n, SizeT m, SizeT nz, RealT kappa, RealT* rd, RealT* rp, RealT* Ctnu);

template <typename RealT, typename SizeT>
void gfgphp(RealT* Q, RealT* R, RealT* Qf, RealT* zmax, RealT* zmin, RealT* z, SizeT T, SizeT n, SizeT m, SizeT nz, RealT* gf, RealT* gp, RealT* hp);

template <typename RealT, typename SizeT>
void resdresp(RealT* rd, RealT* rp, SizeT T, SizeT n, SizeT nz, RealT* resd, RealT* resp, RealT* res);

template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename QfMatrixT,
	typename ZMinVectorT,
	typename ZMaxVectorT,
	typename XVectorT,
	typename Z0VectorT,
	typename SizeT,
	typename RealT
>
void fmpc_solve_impl(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					 ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
					 ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
					 ::boost::numeric::ublas::matrix_expression<QfMatrixT> const& Qf,
					 ::boost::numeric::ublas::vector_expression<ZMinVectorT> const& zmin, 
					 ::boost::numeric::ublas::vector_expression<ZMaxVectorT> const& zmax,
					 ::boost::numeric::ublas::vector_expression<XVectorT> const& x,
					 ::boost::numeric::ublas::vector_container<Z0VectorT>& z0,
					 SizeT T,
					 SizeT niters,
					 RealT kappa);

template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename XMinVectorT,
	typename XMaxVectorT,
	typename UMinVectorT,
	typename UMaxVectorT,
	typename QfMatrixT,
	typename SizeT,
	typename RealT,
	typename X0MatrixT,
	typename U0MatrixT,
	typename X0VectorT,
	typename XMatrixT,
	typename UMatrixT,
	typename XVectorT,
	typename KMatrixT
>
void fmpc_step(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			   ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			   ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
			   ::boost::numeric::ublas::vector_expression<XMinVectorT> const& xmin, 
			   ::boost::numeric::ublas::vector_expression<XMaxVectorT> const& xmax,
			   ::boost::numeric::ublas::vector_expression<UMinVectorT> const& umin, 
			   ::boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
			   ::boost::numeric::ublas::matrix_expression<QfMatrixT> const& Qf,
			   SizeT T,
			   SizeT niters,
			   RealT kappa,
			   ::boost::numeric::ublas::matrix_expression<X0MatrixT> const& X0,
			   ::boost::numeric::ublas::matrix_expression<U0MatrixT> const& U0,
			   ::boost::numeric::ublas::vector_expression<X0VectorT> const& x0,
			   bool want_X,
			   ::boost::numeric::ublas::matrix_container<XMatrixT>& X,
			   bool want_U,
			   ::boost::numeric::ublas::matrix_container<UMatrixT>& U,
			   bool want_x,
			   ::boost::numeric::ublas::vector_container<XVectorT>& x,
			   bool want_K,
			   ::boost::numeric::ublas::matrix_container<KMatrixT>& K);

///@} Prototypes


#if 0
void printmat(double *A, int m, int n)
{
    double *dptr;
    int j, i;
    dptr = A;
	::std::cerr.precision(4);
    for (j = 0; j < m; j++)
    {
        for (i = 0; i < n; i++)
        {
            /*printf("%5.4f\t", *(dptr+m*i+j));*/
            ::std::cerr << ::std::fixed << *(dptr+m*i+j) << "\t";
        }
        /*printf("\n");*/
		::std::cerr << ::std::endl;
    }
    return;
}
#endif // 0

/* computes the search directions dz and dnu */
template <typename RealT, typename SizeT>
void dnudz(RealT* A, RealT* B, RealT* At, RealT* Bt, RealT* eyen, RealT* eyem, RealT* Q, RealT* R, RealT* Qf, RealT* hp, RealT* rd, RealT* rp, SizeT T, SizeT n, SizeT m, SizeT nz, RealT kappa, RealT* dnu, RealT* dz)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace bindings = ::boost::numeric::bindings;

	typedef RealT real_type;
	typedef SizeT size_type;
#if 0
	typedef ublas::array_adaptor<real_type> array_type;
	typedef ublas::matrix<real_type,ublas::column_major,array_type> matrix_type;
	typedef ublas::vector<real_type,array_type> vector_type;
#endif // 0

	size_type info;
    real_type* dptr(0);
	real_type* dptr1(0);
	real_type* dptr2(0);
	real_type* dptr3(0);
    size_type nT = n*T;
	size_type nn = n*n;
	size_type mm = m*m;
	size_type nnT = nn*T;
	size_type nnTm1 = nn*(T-1);
	size_type mmT = mm*T;
	size_type nm = n*m;
	size_type nmT = nm*T;

    /* allocate memory */
    real_type* PhiQ = new real_type[nnT];
    real_type* PhiR = new real_type[mmT];
    real_type* PhiinvQAt = new real_type[nnT];
    real_type* PhiinvRBt = new real_type[nmT];
    real_type* PhiinvQeye = new real_type[nnT];
    real_type* PhiinvReye = new real_type[mmT];
    real_type* CPhiinvrd = new real_type[nT];
    real_type* Yd = new real_type[nnT];
    real_type* Yud = new real_type[nnTm1];
    real_type* Ld = new real_type[nnT];
    real_type* Lld = new real_type[nnTm1];
    real_type* gam = new real_type[nT];
    real_type* v = new real_type[nT];
    real_type* be = new real_type[nT];
    real_type* temp = new real_type[n];
    real_type* tempmatn = new real_type[nn];
    real_type* tempmatm = new real_type[mm];
    real_type* Ctdnu = new real_type[nz];
    real_type* rdmCtdnu = new real_type[nz];

	// auxiliary matrices and vectors
#if 0
	matrix_type A1;
	matrix_type A2;
	matrix_type A3;
	vector_type v1;
	vector_type v2;
#endif // 0

    /* form PhiQ and PhiR */
    for (size_type i = 0; i < T-1; ++i)
    {
        dptr = PhiQ+nn*i; dptr1 = Q;
        for (size_type j = 0; j < nn; ++j)
        {
            *dptr = 2*(*dptr1);
            ++dptr; ++dptr1;
        }
        dptr = PhiQ+nn*i; dptr1 = hp+m*(i+1)+n*i;
        for (size_type j = 0; j < n; ++j)
        {
            *dptr = *dptr+kappa*(*dptr1);
            dptr = dptr+n+1; ++dptr1;
        }
        dptr = PhiR+mm*i; dptr1 = R;
        for (size_type j = 0; j < mm; ++j)
        {
            *dptr = 2*(*dptr1);
            ++dptr; ++dptr1;
        }
        dptr = PhiR+m*m*i; dptr1 = hp+i*(n+m);
        for (size_type j = 0; j < m; ++j)
        {
            *dptr = *dptr+kappa*(*dptr1);
            dptr = dptr+m+1; ++dptr1;
        }
    }
    
    dptr = PhiR+mm*(T-1); dptr1 = R;
    for (size_type j = 0; j < mm; ++j)
    {
        *dptr = 2*(*dptr1);
        ++dptr; ++dptr1;
    }
    dptr = PhiR+mm*(T-1); dptr1 = hp+(T-1)*(n+m);
    for (size_type j = 0; j < m; ++j)
    {
        *dptr = *dptr+kappa*(*dptr1);
        dptr = dptr+m+1; ++dptr1;
    }
    dptr = PhiQ+nn*(T-1); dptr1 = Qf;
    for (size_type j = 0; j < nn; ++j)
    {
        *dptr = 2*(*dptr1);
        ++dptr; ++dptr1;
    }
    dptr = PhiQ+nn*(T-1); dptr1 = hp+m*T+n*(T-1);
    for (size_type j = 0; j < n; ++j)
    {
        *dptr = *dptr+kappa*(*dptr1);
        dptr = dptr+n+1; ++dptr1;
    }

    /* compute PhiinvQAt, PhiinvRBt, PhiinvQeye, PhiinvReye */
    for (size_type i = 0; i < T; ++i)
    {
        dptr = PhiinvQAt+nn*i; dptr1 = At;
        for (size_type j = 0; j < nn; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr = dptr-nn; dptr1 = PhiQ+nn*i;
//		F77_CALL(dposv)("l",&n,&n,dptr1,&n,dptr,&n,&info);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr1));
		A2 = matrix_type(n, n, array_type(nn, dptr));
		bindings::lapack::posv(bindings::lower(A1), A2);
		::std::copy(A1.data().begin(), A1.data().end(), dptr1);
		::std::copy(A2.data().begin(), A2.data().end(), dptr);
#else
		LAPACK_DPOSV("l",&n,&n,dptr1,&n,dptr,&n, &info);
#endif // 0
        dptr = PhiinvQeye+nn*i; dptr1 = eyen;
        for (size_type j = 0; j < nn; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr = dptr-nn; dptr1 = PhiQ+nn*i;
//		F77_CALL(dtrtrs)("l","n","n",&n,&n,dptr1,&n,dptr,&n,&info);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr1));
		A2 = matrix_type(n, n, array_type(nn, dptr));
		bindings::lapack::trtrs(bindings::lower(A1), A2);
		::std::copy(A2.data().begin(), A2.data().end(), dptr);
#else
		LAPACK_DTRTRS("l","n","n",&n,&n,dptr1,&n,dptr,&n,&info);
#endif // 0
//		F77_CALL(dtrtrs)("l","t","n",&n,&n,dptr1,&n,dptr,&n,&info);
#if 0
		bindings::lapack::trtrs(bindings::trans(bindings::lower(A1)), A2);
		::std::copy(A2.data().begin(), A2.data().end(), dptr);
#else
		LAPACK_DTRTRS("l","t","n",&n,&n,dptr1,&n,dptr,&n,&info);
#endif // 0
    }
    for (size_type i = 0; i < T; ++i)
    {
        dptr = PhiinvRBt+nm*i; dptr1 = Bt;
        for (size_type j = 0; j < nm; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr = dptr-nm; dptr1 = PhiR+mm*i;
//		F77_CALL(dposv)("l",&m,&n,dptr1,&m,dptr,&m,&info);
#if 0
		A1 = matrix_type(m, m, array_type(mm, dptr1));
		A2 = matrix_type(m, n, array_type(nm, dptr));
		bindings::lapack::posv(bindings::lower(A1), A2);
		::std::copy(A1.data().begin(), A1.data().end(), dptr1);
		::std::copy(A2.data().begin(), A2.data().end(), dptr);
#else
		LAPACK_DPOSV("l",&m,&n,dptr1,&m,dptr,&m,&info);
#endif // 0
        dptr = PhiinvReye+mm*i; dptr1 = eyem;
        for (size_type j = 0; j < mm; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr = dptr-mm; dptr1 = PhiR+mm*i;
//		F77_CALL(dtrtrs)("l","n","n",&m,&m,dptr1,&m,dptr,&m,&info);
#if 0
		A1 = matrix_type(m, m, array_type(mm, dptr1));
		A2 = matrix_type(m, m, array_type(mm, dptr));
		bindings::lapack::trtrs(bindings::lower(A1), A2);
#else
		LAPACK_DTRTRS("l","n","n",&m,&m,dptr1,&m,dptr,&m,&info);
#endif // 0
//		F77_CALL(dtrtrs)("l","t","n",&m,&m,dptr1,&m,dptr,&m,&info);
#if 0
		bindings::lapack::trtrs(bindings::trans(bindings::lower(A1)), A2);
		::std::copy(A2.data().begin(), A2.data().end(), dptr);
#else
		LAPACK_DTRTRS("l","t","n",&m,&m,dptr1,&m,dptr,&m,&info);
#endif // 0
    }
    
    /* form Yd and Yud */
    dptr = Yud; dptr1 = PhiinvQAt; 
    for (size_type i = 0; i < nnTm1; ++i)
    {
        *dptr = -(*dptr1);
        ++dptr; ++dptr1;
    }
    dptr2 = Yd; dptr3 = PhiinvQeye;
    for (size_type i = 0; i < nn; ++i)
    {
        *dptr2 = *dptr3;
        ++dptr2; ++dptr3;
    }
    dptr2 = dptr2-nn;
//	F77_CALL(dgemm)("n","n",&n,&n,&m,&fone,B,&n,PhiinvRBt,&m,&fone,dptr2,&n);
#if 0
	A1 = matrix_type(n, m, array_type(nm, B));
	A2 = matrix_type(m, n, array_type(nm, PhiinvRBt));
	A3 = matrix_type(n, n, array_type(nn, dptr2));
	bindings::blas::gemm(fone, A1, A2, fone, A3);
	::std::copy(A3.data().begin(), A3.data().end(), dptr2);
#else
	BLAS_DGEMM("n","n",&n,&n,&m,&fone,B,&n,PhiinvRBt,&m,&fone,dptr2,&n);
#endif // 0
    for (size_type i = 1; i < T; ++i)
    {
        dptr = Yd+nn*i; dptr1 = PhiinvQeye+nn*i;
        for (size_type j = 0; j < nn; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr1 = PhiinvRBt+nm*i; dptr = dptr-nn;
//		F77_CALL(dgemm)("n","n",&n,&n,&m,&fone,B,&n,dptr1,&m,&fone,dptr,&n); 
#if 0
		A1 = matrix_type(n, m, array_type(nm, B));
		A2 = matrix_type(m, n, array_type(nm, dptr1));
		A3 = matrix_type(n, n, array_type(nn, dptr));
		bindings::blas::gemm(fone, A1, A2, fone, A3);
		::std::copy(A3.data().begin(), A3.data().end(), dptr);
#else
		BLAS_DGEMM("n","n",&n,&n,&m,&fone,B,&n,dptr1,&m,&fone,dptr,&n); 
#endif // 0
        dptr1 = PhiinvQAt+nn*(i-1);
//		F77_CALL(dgemm)("n","n",&n,&n,&n,&fone,A,&n,dptr1,&n,&fone,dptr,&n);
#if 0
		A1 = matrix_type(n, n, array_type(nn, A));
		A2 = matrix_type(n, n, array_type(nn, dptr1));
		A3 = matrix_type(n, n, array_type(nn, dptr));
		bindings::blas::gemm(fone, A1, A2, fone, A3);
		::std::copy(A3.data().begin(), A3.data().end(), dptr);
#else
		BLAS_DGEMM("n","n",&n,&n,&n,&fone,A,&n,dptr1,&n,&fone,dptr,&n);
#endif // 0
    }

    /* compute Lii */
    dptr = Ld; dptr1 = Yd; 
    for (size_type i = 0; i < nn; ++i)
    {
        *dptr = *dptr1;
        ++dptr; ++dptr1; 
    }
    dptr = dptr-nn; 
//	F77_CALL(dposv)("l",&n,&ione,dptr,&n,temp,&n,&info);
#if 0
	A1 = matrix_type(n, n, array_type(nn, dptr));
	A2 = matrix_type(n, ione, array_type(n, temp));
	bindings::lapack::posv(bindings::lower(A1), A2);
	::std::copy(A1.data().begin(), A1.data().end(), dptr);
	::std::copy(A2.data().begin(), A2.data().end(), temp);
#else
	LAPACK_DPOSV("l",&n,&ione,dptr,&n,temp,&n,&info);
#endif // 0
    for (size_type i = 1; i < T; ++i)
    {
        dptr = Ld+nn*(i-1); dptr1 = Yud+nn*(i-1); dptr2 = Lld+nn*(i-1);
        for (size_type j = 0; j < nn; ++j)
        {
            *dptr2 = *dptr1;
            ++dptr2; ++dptr1;
        }
        dptr2 = dptr2-nn;
//		F77_CALL(dtrtrs)("l","n","n",&n,&n,dptr,&n,dptr2,&n,&info);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr));
		A2 = matrix_type(n, n, array_type(nn, dptr2));
		bindings::lapack::trtrs(bindings::lower(A1), A2);
		::std::copy(A2.data().begin(), A2.data().end(), dptr2);
#else
		LAPACK_DTRTRS("l","n","n",&n,&n,dptr,&n,dptr2,&n,&info);
#endif // 0
        dptr1 = tempmatn;
        for (size_type j = 0; j < nn; ++j)
        {
            *dptr1 = *dptr2;
            ++dptr1; ++dptr2;
        }
        dptr1 = dptr1-nn; dptr2 = dptr2-nn;
//		F77_CALL(dgemm)("t","n",&n,&n,&n,&fone,dptr1,&n,eyen,&n,&fzero,dptr2,&n);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr1));
		A2 = matrix_type(n, n, array_type(nn, eyen));
		A3 = matrix_type(n, n, array_type(nn, dptr2));
		bindings::blas::gemm(fone, bindings::trans(A1), A2, fzero, A3);
		::std::copy(A3.data().begin(), A3.data().end(), dptr2);
#else
		BLAS_DGEMM("t","n",&n,&n,&n,&fone,dptr1,&n,eyen,&n,&fzero,dptr2,&n);
#endif // 0
        dptr = Ld+nn*i; dptr1 = Yd+nn*i;
        for (size_type j = 0; j < nn; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr = dptr-nn;
//		F77_CALL(dgemm)("n","t",&n,&n,&n,&fmone,dptr2,&n,dptr2,&n,&fone,dptr,&n);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr2));
		A2 = matrix_type(n, n, array_type(nn, dptr2));
		A3 = matrix_type(n, n, array_type(nn, dptr));
		bindings::blas::gemm(fmone, A1, bindings::trans(A2), fone, A3);
		::std::copy(A3.data().begin(), A3.data().end(), dptr);
#else
		BLAS_DGEMM("n","t",&n,&n,&n,&fmone,dptr2,&n,dptr2,&n,&fone,dptr,&n);
#endif // 0
//		F77_CALL(dposv)("l",&n,&ione,dptr,&n,temp,&n,&info);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr));
		A2 = matrix_type(n, ione, array_type(n, temp));
		bindings::lapack::posv(bindings::lower(A1), A2);
		::std::copy(A1.data().begin(), A1.data().end(), dptr);
		::std::copy(A2.data().begin(), A2.data().end(), temp);
#else
		LAPACK_DPOSV("l",&n,&ione,dptr,&n,temp,&n,&info);
#endif // 0
    }

    /* compute CPhiinvrd */
    dptr = CPhiinvrd; dptr1 = rd+m;
    for (size_type i = 0; i < n; ++i)
    {
        *dptr = *dptr1;
        ++dptr; ++dptr1;
    }
    dptr = dptr-n;
//	F77_CALL(dtrsv)("l","n","n",&n,PhiQ,&n,dptr,&ione);
#if 0
	A1 = matrix_type(n, n, array_type(nn, PhiQ));
	v1 = vector_type(n, array_type(n, dptr));
	bindings::blas::trsv(bindings::lower(A1), v1);
#else
	BLAS_DTRSV("l","n","n",&n,PhiQ,&n,dptr,&ione);
#endif // 0
//	F77_CALL(dtrsv)("l","t","n",&n,PhiQ,&n,dptr,&ione);
#if 0
//	A1 = matrix_type(n, n, array_type(nn, PhiQ));
//	v1 = vector_type(n, array_type(n, dptr));
	bindings::blas::trsv(bindings::trans(bindings::lower(A1)), v1);
	::std::copy(v1.data().begin(), v1.data().end(), dptr);
#else
	BLAS_DTRSV("l","t","n",&n,PhiQ,&n,dptr,&ione);
#endif // 0
    dptr2 = temp; dptr1 = rd;
    for (size_type i = 0; i < m; ++i)
    {
        *dptr2 = *dptr1;
        ++dptr2; ++dptr1;
    }
    dptr2 = dptr2-m;
//	F77_CALL(dtrsv)("l","n","n",&m,PhiR,&m,dptr2,&ione);
#if 0
	A1 = matrix_type(m, m, array_type(mm, PhiR));
	v1 = vector_type(m, array_type(m, dptr2));
	bindings::blas::trsv(bindings::lower(A1), v1);
#else
	BLAS_DTRSV("l","n","n",&m,PhiR,&m,dptr2,&ione);
#endif // 0
//	F77_CALL(dtrsv)("l","t","n",&m,PhiR,&m,dptr2,&ione);
#if 0
//	A1 = matrix_type(m, m, array_type(mm, PhiR));
//	v1 = vector_type(m, array_type(m, dptr2));
	bindings::blas::trsv(bindings::trans(bindings::lower(A1)), v1);
	::std::copy(v1.data().begin(), v1.data().end(), dptr2);
#else
	BLAS_DTRSV("l","t","n",&m,PhiR,&m,dptr2,&ione);
#endif // 0
//	F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,temp,&ione,&fone,dptr,&ione);
#if 0
	A1 = matrix_type(n, m, array_type(nm, B));
	v1 = vector_type(m, array_type(m, temp));
	v2 = vector_type(n, array_type(n, dptr));
	bindings::blas::gemv(fmone, A1, v1, fone, v2);
	::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
	BLAS_DGEMV("n",&n,&m,&fmone,B,&n,temp,&ione,&fone,dptr,&ione);
#endif // 0
    
    for (size_type i = 1; i < T; ++i)
    {
        dptr = CPhiinvrd+n*i; dptr1 = rd+m+i*(n+m);
        for (size_type j = 0; j < n; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr = dptr-n; dptr3 = PhiQ+nn*i;
//		F77_CALL(dtrsv)("l","n","n",&n,dptr3,&n,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr3));
		v1 = vector_type(n, array_type(n, dptr));
		bindings::blas::trsv(bindings::lower(A1), v1);
#else
		BLAS_DTRSV("l","n","n",&n,dptr3,&n,dptr,&ione);
#endif // 0
//		F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr,&ione);
#if 0
//		A1 = matrix_type(n, n, array_type(nn, dptr3));
//		v1 = vector_type(n, array_type(n, dptr));
		bindings::blas::trsv(bindings::trans(bindings::lower(A1)), v1);
		::std::copy(v1.data().begin(), v1.data().end(), dptr);
#else
		BLAS_DTRSV("l","t","n",&n,dptr3,&n,dptr,&ione);
#endif // 0
        dptr2 = temp; dptr1 = rd+i*(m+n);
        for (size_type j = 0; j < m; ++j)
        {
            *dptr2 = *dptr1;
            ++dptr2; ++dptr1;
        }
        dptr3 = PhiR+mm*i; dptr2 = dptr2-m;
//		F77_CALL(dtrsv)("l","n","n",&m,dptr3,&m,dptr2,&ione);
#if 0
		A1 = matrix_type(m, m, array_type(mm, dptr3));
		v1 = vector_type(m, array_type(m, dptr2));
		bindings::blas::trsv(bindings::lower(A1), v1);
#else
		BLAS_DTRSV("l","n","n",&m,dptr3,&m,dptr2,&ione);
#endif // 0
//		F77_CALL(dtrsv)("l","t","n",&m,dptr3,&m,dptr2,&ione);
#if 0
//		A1 = matrix_type(m, m, array_type(mm, dptr3));
//		v1 = vector_type(m, array_type(m, dptr2));
		bindings::blas::trsv(bindings::trans(bindings::lower(A1)), v1);
		::std::copy(v1.data().begin(), v1.data().end(), dptr2);
#else
		BLAS_DTRSV("l","t","n",&m,dptr3,&m,dptr2,&ione);
#endif // 0
//		F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,temp,&ione,&fone,dptr,&ione);
#if 0
		A1 = matrix_type(n, m, array_type(nm, B));
		v1 = vector_type(m, array_type(m, temp));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fmone, A1, v1, fone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&m,&fmone,B,&n,temp,&ione,&fone,dptr,&ione);
#endif // 0
        dptr2 = temp; dptr1 = rd+(i-1)*(n+m)+m;
        for (size_type j = 0; j < n; ++j)
        {
            *dptr2 = *dptr1;
            ++dptr2; ++dptr1;
        }
        dptr3 = PhiQ+nn*(i-1); dptr2 = dptr2-n;
//		F77_CALL(dtrsv)("l","n","n",&n,dptr3,&n,dptr2,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr3));
		v1 = vector_type(n, array_type(n, dptr2));
		bindings::blas::trsv(bindings::lower(A1), v1);
#else
		BLAS_DTRSV("l","n","n",&n,dptr3,&n,dptr2,&ione);
#endif // 0
//		F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr2,&ione);
#if 0
//		A1 = matrix_type(n, n, array_type(nn, dptr3));
//		v1 = vector_type(n, array_type(n, dptr2));
		bindings::blas::trsv(bindings::trans(bindings::lower(A1)), v1);
		::std::copy(v1.data().begin(), v1.data().end(), dptr2);
#else
		BLAS_DTRSV("l","t","n",&n,dptr3,&n,dptr2,&ione);
#endif // 0
//		F77_CALL(dgemv)("n",&n,&n,&fmone,A,&n,temp,&ione,&fone,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, A));
		v1 = vector_type(n, array_type(n, temp));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fmone, A1, v1, fone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&n,&fmone,A,&n,temp,&ione,&fone,dptr,&ione);
#endif // 0
    }

    /* form be */
    dptr = be; dptr1 = rp; dptr2 = CPhiinvrd;
    for (size_type i = 0; i < nT; ++i)
    {
        *dptr = (*dptr2)-(*dptr1);
        ++dptr; ++dptr1; ++dptr2;
    }

    /* solve for dnu */
    dptr = v; dptr1 = be;
    for (size_type i = 0; i < n; ++i)
    {
        *dptr = -(*dptr1);
        ++dptr; ++dptr1;
    }
    dptr = dptr-n;
//	F77_CALL(dtrsv)("l","n","n",&n,Ld,&n,dptr,&ione);
#if 0
	A1 = matrix_type(n, n, array_type(nn, Ld));
	v1 = vector_type(n, array_type(n, dptr));
	bindings::blas::trsv(bindings::lower(A1), v1);
	::std::copy(v1.data().begin(), v1.data().end(), dptr);
#else
	BLAS_DTRSV("l","n","n",&n,Ld,&n,dptr,&ione);
#endif // 0
    for (size_type i = 1; i < T; ++i)
    {
        dptr = v+i*n; dptr1 = v+(i-1)*n; dptr2 = be+i*n; 
        for (size_type j = 0; j < n; ++j)
        {
            *dptr = *dptr2;
            ++dptr; ++dptr2;
        }
        dptr = dptr-n; dptr3 = Lld+nn*(i-1);
//		F77_CALL(dgemv)("n",&n,&n,&fmone,dptr3,&n,dptr1,&ione,&fmone,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr3));
		v1 = vector_type(n, array_type(n, dptr1));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fmone, A1, v1, fmone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&n,&fmone,dptr3,&n,dptr1,&ione,&fmone,dptr,&ione);
#endif // 0
        dptr3 = Ld+nn*i;
//		F77_CALL(dtrsv)("l","n","n",&n,dptr3,&n,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr3));
		v1 = vector_type(n, array_type(n, dptr));
		bindings::blas::trsv(bindings::lower(A1), v1);
		::std::copy(v1.data().begin(), v1.data().end(), dptr);
#else
		BLAS_DTRSV("l","n","n",&n,dptr3,&n,dptr,&ione);
#endif // 0
    }
    dptr = dnu+n*(T-1); dptr1 = v+n*(T-1);
    for (size_type i = 0; i < n; ++i)
    {
        *dptr = *dptr1;
        ++dptr; ++dptr1;
    }
    dptr = dptr-n; dptr3 = Ld+nn*(T-1);
//	F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr,&ione);
#if 0
	A1 = matrix_type(n, n, array_type(nn, dptr3));
	v1 = vector_type(n, array_type(n, dptr));
	bindings::blas::trsv(bindings::trans(bindings::lower(A1)), v1);
	::std::copy(v1.data().begin(), v1.data().end(), dptr);
#else
	BLAS_DTRSV("l","t","n",&n,dptr3,&n,dptr,&ione);
#endif // 0
    for (size_type i = T-1; i > 0; --i)
    {
        dptr = dnu+n*(i-1); dptr1 = dnu+n*i; dptr2 = v+n*(i-1); 
        for (size_type j = 0; j < n; ++j)
        {
            *dptr = *dptr2;
            ++dptr; ++dptr2;
        }
        dptr = dptr-n; dptr3 = Lld+nn*(i-1);
//		F77_CALL(dgemv)("t",&n,&n,&fmone,dptr3,&n,dptr1,&ione,&fone,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr3));
		v1 = vector_type(n, array_type(n, dptr1));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fmone, A1, v1, fone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("t",&n,&n,&fmone,dptr3,&n,dptr1,&ione,&fone,dptr,&ione);
#endif // 0
        dptr3 = Ld+nn*(i-1);
//		F77_CALL(dtrsv)("l","t","n",&n,dptr3,&n,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr3));
		v1 = vector_type(n, array_type(n, dptr));
		bindings::blas::trsv(bindings::trans(bindings::lower(A1)), v1);
		::std::copy(v1.data().begin(), v1.data().end(), dptr);
#else
		BLAS_DTRSV("l","t","n",&n,dptr3,&n,dptr,&ione);
#endif // 0
    }

    /* form Ctdnu */
    for (size_type i = 0; i < T-1; ++i)
    {
        dptr = Ctdnu+i*(n+m); dptr1 = dnu+i*n;
//		F77_CALL(dgemv)("n",&m,&n,&fmone,Bt,&m,dptr1,&ione,&fzero,dptr,&ione);
#if 0
		A1 = matrix_type(m, n, array_type(nm, Bt));
		v1 = vector_type(n, array_type(n, dptr1));
		v2 = vector_type(m, array_type(m, dptr));
		bindings::blas::gemv(fmone, A1, v1, fzero, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&m,&n,&fmone,Bt,&m,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0
        dptr = Ctdnu+i*(n+m)+m; dptr2 = dnu+(i+1)*n;
        for (size_type j = 0; j < n; ++j)
        {
            *dptr = *dptr1;
            ++dptr; ++dptr1;
        }
        dptr = dptr-n;
//		F77_CALL(dgemv)("n",&n,&n,&fmone,At,&n,dptr2,&ione,&fone,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, At));
		v1 = vector_type(n, array_type(n, dptr2));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fmone, A1, v1, fone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&n,&fmone,At,&n,dptr2,&ione,&fone,dptr,&ione);
#endif // 0
    }
    
    dptr = Ctdnu+(T-1)*(n+m); dptr1 = dnu+(T-1)*n;
//	F77_CALL(dgemv)("n",&m,&n,&fmone,Bt,&m,dptr1,&ione,&fzero,dptr,&ione);
#if 0
	A1 = matrix_type(m, n, array_type(nm, Bt));
	v1 = vector_type(n, array_type(n, dptr1));
	v2 = vector_type(m, array_type(m, dptr));
	bindings::blas::gemv(fmone, A1, v1, fzero, v2);
	::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
	BLAS_DGEMV("n",&m,&n,&fmone,Bt,&m,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0
    dptr = dptr+m; 
    for (size_type i = 0; i < n; ++i)
    {
        *dptr = *dptr1;
        ++dptr; ++dptr1;
    }
    dptr = rdmCtdnu; dptr1 = Ctdnu; dptr2 = rd;
    for (size_type i = 0; i < nz; ++i)
    {
        *dptr = -(*dptr1)-(*dptr2);
        ++dptr; ++dptr1; ++dptr2;
    }

    /* solve for dz */
    for (size_type i = 0; i < T; ++i)
    {
        dptr = dz+(i+1)*m+i*n; dptr1 = rdmCtdnu+(i+1)*m+i*n;
        dptr2 = PhiinvQeye+nn*i;
//		F77_CALL(dgemv)("n",&n,&n,&fone,dptr2,&n,dptr1,&ione,&fzero,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, dptr2));
		v1 = vector_type(n, array_type(n, dptr1));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fone, A1, v1, fzero, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&n,&fone,dptr2,&n,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0
    }
    for (size_type i = 0; i < T; ++i)
    {
        dptr = dz+i*(m+n); dptr1 = rdmCtdnu+i*(m+n);
        dptr2 = PhiinvReye+mm*i;
//		F77_CALL(dgemv)("n",&m,&m,&fone,dptr2,&m,dptr1,&ione,&fzero,dptr,&ione);
#if 0
		A1 = matrix_type(m, m, array_type(mm, dptr2));
		v1 = vector_type(m, array_type(m, dptr1));
		v2 = vector_type(m, array_type(m, dptr));
		bindings::blas::gemv(fone, A1, v1, fzero, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&m,&m,&fone,dptr2,&m,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0
    }

    delete[] PhiQ;
    delete[] PhiR;
    delete[] PhiinvQAt;
    delete[] PhiinvRBt;
    delete[] PhiinvQeye;
    delete[] PhiinvReye;
    delete[] CPhiinvrd;
    delete[] Yd;
    delete[] Yud;
    delete[] Ld;
    delete[] Lld;
    delete[] gam;
    delete[] v;
    delete[] be;
    delete[] temp;
    delete[] tempmatn;
    delete[] tempmatm;
    delete[] Ctdnu;
    delete[] rdmCtdnu;
}

/* computes rd and rp */
template <typename RealT, typename SizeT>
void rdrp(RealT* A, RealT* B, RealT* z, RealT* nu, RealT* gf, RealT* gp, RealT* b, SizeT T, SizeT n, SizeT m, SizeT nz, RealT kappa, RealT* rd, RealT* rp, RealT* Ctnu)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace bindings = ::boost::numeric::bindings;

	typedef RealT real_type;
	typedef SizeT size_type;
#if 0
	typedef ublas::array_adaptor<real_type> array_type;
	typedef ublas::matrix<real_type,ublas::column_major,array_type> matrix_type;
	typedef ublas::vector<real_type,array_type> vector_type;
#endif // 0

    real_type* dptr(0);
	real_type* dptr1(0);
	real_type* dptr2(0);
	size_type nT = n*T;
#if 0
	size_type nn = n*n;
	size_type nm = n*m;

	// auxiliary matrices and vectors
	matrix_type A1;
	vector_type v1;
	vector_type v2;
#endif // 0

    real_type* Cz = new real_type[nT];
    
    /* compute Cz */
    dptr = Cz; dptr1 = z+m;
    for (size_type i = 0; i < n; ++i)
    {
        *dptr = *dptr1;
        ++dptr; ++dptr1;
    }
//	F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,z,&ione,&fone,Cz,&ione);
#if 0
	A1 = matrix_type(n, m, array_type(nm, B));
	v1 = vector_type(m, array_type(m, z));
	v2 = vector_type(n, array_type(n, Cz));
	bindings::blas::gemv(fmone, A1, v1, fone, v2);
	::std::copy(v2.data().begin(), v2.data().end(), Cz);
#else
	BLAS_DGEMV("n",&n,&m,&fmone,B,&n,z,&ione,&fone,Cz,&ione);
#endif // 0
    for (size_type i = 2; i <= T; ++i)
    {
        dptr = Cz+(i-1)*n; dptr1 = z+m+(i-2)*(n+m); 
        dptr2 = z+m+(i-1)*(m+n);
        for (size_type j = 0; j < n; ++j)
        {
            *dptr = *dptr2;
            ++dptr; ++dptr2;
        }
        dptr = dptr-n; 
//		F77_CALL(dgemv)("n",&n,&n,&fmone,A,&n,dptr1,&ione,&fone,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, A));
		v1 = vector_type(n, array_type(n, dptr1));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fmone, A1, v1, fone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&n,&fmone,A,&n,dptr1,&ione,&fone,dptr,&ione);
#endif // 0
        dptr1 = dptr1+n;
//		F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,dptr1,&ione,&fone,dptr,&ione);
#if 0
		A1 = matrix_type(n, m, array_type(nm, B));
		v1 = vector_type(m, array_type(m, dptr1));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(fmone, A1, v1, fone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&m,&fmone,B,&n,dptr1,&ione,&fone,dptr,&ione);
#endif // 0
    }
    /*
    dptr = Cz+(T-1)*n; dptr1 = z+m+(T-2)*(n+m);
    F77_CALL(dgemv)("n",&n,&n,&fmone,A,&n,dptr1,&ione,&fzero,dptr,&ione);
    dptr1 = dptr1+n;
    F77_CALL(dgemv)("n",&n,&m,&fmone,B,&n,dptr1,&ione,&fone,dptr,&ione);
    dptr1 = z+nz-n;
    for (i = 0; i < n; i++)
    {
        *dptr = *dptr+*dptr1;
        dptr++; dptr1++;
    }
    */

    /* compute Ctnu */
    dptr = Ctnu; dptr1 = Ctnu+m; dptr2 = nu;
    for (size_type i = 1; i <= T-1; ++i)
    {
//		F77_CALL(dgemv)("t",&n,&m,&fmone,B,&n,dptr2,&ione,&fzero,dptr,&ione);
#if 0
		A1 = matrix_type(n, m, array_type(nm, B));
		v1 = vector_type(n, array_type(n, dptr2));
		v2 = vector_type(m, array_type(m, dptr));
		bindings::blas::gemv(fmone, bindings::trans(A1), v1, fzero, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("t",&n,&m,&fmone,B,&n,dptr2,&ione,&fzero,dptr,&ione);
#endif // 0
        dptr = dptr+n+m;
        for (size_type j = 0; j < n; ++j)
        {
            *dptr1 = *dptr2;
            ++dptr1; ++dptr2;
        }
        dptr1 = Ctnu+m+(i-1)*(n+m);
//		F77_CALL(dgemv)("t",&n,&n,&fmone,A,&n,dptr2,&ione,&fone,dptr1,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, A));
		v1 = vector_type(n, array_type(n, dptr2));
		v2 = vector_type(n, array_type(n, dptr1));
		bindings::blas::gemv(fmone, bindings::trans(A1), v1, fone, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr1);
#else
		BLAS_DGEMV("t",&n,&n,&fmone,A,&n,dptr2,&ione,&fone,dptr1,&ione);
#endif // 0
        dptr1 = dptr1+n+m;
    }
//	F77_CALL(dgemv)("t",&n,&m,&fmone,B,&n,dptr2,&ione,&fzero,dptr,&ione);
#if 0
	A1 = matrix_type(n, m, array_type(nm, B));
	v1 = vector_type(n, array_type(n, dptr2));
	v2 = vector_type(m, array_type(m, dptr));
	bindings::blas::gemv(fmone, bindings::trans(A1), v1, fzero, v2);
	::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
	BLAS_DGEMV("t",&n,&m,&fmone,B,&n,dptr2,&ione,&fzero,dptr,&ione);
#endif // 0
    dptr = Ctnu+nz-n; dptr1 = nu+(T-1)*n;
    for (size_type i = 0; i < n; ++i)
    {
        *dptr = *dptr1;
        ++dptr; ++dptr1;
    }

    dptr = rp; dptr1 = Cz; dptr2 = b;
    for (size_type i = 0; i < nT; ++i)
    {
        *dptr = *dptr1-*dptr2;
        ++dptr; ++dptr1; ++dptr2;
    }
    dptr = rd; dptr1 = Ctnu; dptr2 = gf;
    for (size_type i = 0; i < nz; ++i)
    {
        *dptr = *dptr1+*dptr2;
        ++dptr; ++dptr1; ++dptr2;
    }
//	F77_CALL(daxpy)(&nz,&kappa,gp,&ione,rd,&ione);
#if 0
	v1 = vector_type(nz, array_type(nz, gp));
	v2 = vector_type(nz, array_type(nz, rd));
	bindings::blas::axpy(kappa, v1, v2);
	::std::copy(v2.data().begin(), v2.data().end(), rd);
#else
	BLAS_DAXPY(&nz,&kappa,gp,&ione,rd,&ione);
#endif // 0

    delete[] Cz;
}

/* computes gf, gp and hp */
template <typename RealT, typename SizeT>
void gfgphp(RealT* Q, RealT* R, RealT* Qf, RealT* zmax, RealT* zmin, RealT* z, SizeT T, SizeT n, SizeT m, SizeT nz, RealT* gf, RealT* gp, RealT* hp)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace bindings = ::boost::numeric::bindings;

	typedef RealT real_type;
	typedef SizeT size_type;
#if 0
	typedef ublas::array_adaptor<real_type> array_type;
	typedef ublas::matrix<real_type,ublas::column_major,array_type> matrix_type;
	typedef ublas::vector<real_type,array_type> vector_type;
#endif // 0

#if 0
	size_type nn = n*n;
	size_type mm = m*m;

	// auxiliary matrices and vectors
	matrix_type A1;
	vector_type v1;
	vector_type v2;
#endif // 0

    real_type* gp1 = new real_type[nz];
    real_type* gp2 = new real_type[nz];

    real_type* dptr(0);
	real_type* dptr1(0);
	real_type* dptr2(0);

    dptr = gp1; dptr1 = zmax; dptr2 = z;
    for (size_type i = 0; i < nz; ++i)
    {
        *dptr = 1.0/(*dptr1-*dptr2);
        ++dptr; ++dptr1; ++dptr2;
    }
    dptr = gp2; dptr1 = zmin; dptr2 = z;
    for (size_type i = 0; i < nz; ++i)
    {
        *dptr = 1.0/(*dptr2-*dptr1);
        ++dptr; ++dptr1; ++dptr2;
    }
    dptr = hp; dptr1 = gp1; dptr2 = gp2;
    for (size_type i = 0; i < nz; ++i)
    {
        *dptr = (*dptr1)*(*dptr1) + (*dptr2)*(*dptr2);
        ++dptr; ++dptr1; ++dptr2;
    }
    dptr = gp; dptr1 = gp1; dptr2 = gp2;
    for (size_type i = 0; i < nz; ++i)
    {
        *dptr = *dptr1-*dptr2;
        ++dptr; ++dptr1; ++dptr2;
    }
    
    dptr = gf; dptr1 = z; 
	size_type Tm1 = T-1;
    for (size_type i = 0; i < Tm1; ++i)
    {
//		F77_CALL(dgemv)("n",&m,&m,&ftwo,R,&m,dptr1,&ione,&fzero,dptr,&ione);
#if 0
		A1 = matrix_type(m, m, array_type(mm, R));
		v1 = vector_type(m, array_type(n, dptr1));
		v2 = vector_type(m, array_type(m, dptr));
		bindings::blas::gemv(ftwo, A1, v1, fzero, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&m,&m,&ftwo,R,&m,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0
        dptr = dptr+m; dptr1 = dptr1+m;
//		F77_CALL(dgemv)("n",&n,&n,&ftwo,Q,&n,dptr1,&ione,&fzero,dptr,&ione);
#if 0
		A1 = matrix_type(n, n, array_type(nn, Q));
		v1 = vector_type(n, array_type(n, dptr1));
		v2 = vector_type(n, array_type(n, dptr));
		bindings::blas::gemv(ftwo, A1, v1, fzero, v2);
		::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
		BLAS_DGEMV("n",&n,&n,&ftwo,Q,&n,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0
        dptr = dptr+n; dptr1 = dptr1+n;
    }
//	F77_CALL(dgemv)("n",&m,&m,&ftwo,R,&m,dptr1,&ione,&fzero,dptr,&ione);
#if 0
	A1 = matrix_type(m, m, array_type(mm, R));
	v1 = vector_type(m, array_type(m, dptr1));
	v2 = vector_type(m, array_type(m, dptr));
	bindings::blas::gemv(ftwo, A1, v1, fzero, v2);
	::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
	BLAS_DGEMV("n",&m,&m,&ftwo,R,&m,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0
    dptr = dptr+m; dptr1 = dptr1+m;
//	F77_CALL(dgemv)("n",&n,&n,&ftwo,Qf,&n,dptr1,&ione,&fzero,dptr,&ione);
#if 0
	A1 = matrix_type(n, n, array_type(nn, Qf));
	v1 = vector_type(n, array_type(n, dptr1));
	v2 = vector_type(n, array_type(n, dptr));
	bindings::blas::gemv(ftwo, A1, v1, fzero, v2);
	::std::copy(v2.data().begin(), v2.data().end(), dptr);
#else
	BLAS_DGEMV("n",&n,&n,&ftwo,Qf,&n,dptr1,&ione,&fzero,dptr,&ione);
#endif // 0

    delete[] gp1;
	delete[] gp2;
}

/* computes resd, resp, and res */
template <typename RealT, typename SizeT>
void resdresp(RealT* rd, RealT* rp, SizeT T, SizeT n, SizeT nz, RealT* resd, RealT* resp, RealT* res)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace bindings = ::boost::numeric::bindings;

	typedef RealT real_type;
	typedef SizeT size_type;
#if 0
	typedef ublas::array_adaptor<real_type> array_type;
	typedef ublas::vector<real_type,array_type> vector_type;
#endif // 0

    size_type nnu = T*n;

#if 0
	// auxiliary vector
	vector_type v;
#endif // 0

//	*resp = F77_CALL(dnrm2)(&nnu,rp,&ione);
#if 0
	v = vector_type(nnu, array_type(nnu, rp));
	*resp = bindings::blas::nrm2(v);
#else
	*resp = BLAS_DNRM2(&nnu,rp,&ione);
#endif // 0
//	*resd = F77_CALL(dnrm2)(&nz,rd,&ione);
#if 0
	v = vector_type(nz, array_type(nz, rd));
	*resd = bindings::blas::nrm2(v);
#else
	*resd = BLAS_DNRM2(&nz,rd,&ione);
#endif // 0
    *res = ::std::sqrt((*resp)*(*resp)+(*resd)*(*resd));
}

template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename QfMatrixT,
	typename ZMinVectorT,
	typename ZMaxVectorT,
	typename XVectorT,
	typename Z0VectorT,
	typename SizeT,
	typename RealT
>
void fmpc_solve_impl(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					 ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
					 ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
					 ::boost::numeric::ublas::matrix_expression<QfMatrixT> const& Qf,
					 ::boost::numeric::ublas::vector_expression<ZMinVectorT> const& zmin, 
					 ::boost::numeric::ublas::vector_expression<ZMaxVectorT> const& zmax,
					 ::boost::numeric::ublas::vector_expression<XVectorT> const& x,
					 ::boost::numeric::ublas::vector_container<Z0VectorT>& z0,
					 SizeT T,
					 SizeT niters,
					 RealT kappa)
{
	namespace ublas = ::boost::numeric::ublas;
	namespace ublasx = ::boost::numeric::ublasx;

	typedef AMatrixT A_matrix_type;
	typedef BMatrixT B_matrix_type;
	typedef QMatrixT Q_matrix_type;
	typedef RMatrixT R_matrix_type;
	typedef QfMatrixT Qf_matrix_type;
	typedef ZMinVectorT zmin_vector_type;
	typedef ZMaxVectorT zmax_vector_type;
	typedef XVectorT x_vector_type;
	typedef Z0VectorT z0_vector_type;
	typedef SizeT size_type;
	typedef RealT real_type;
	//TODO: type promotion among all matrices for the real_type
	typedef ublas::vector<real_type> work_vector_type;
	typedef ublas::matrix<real_type, ublas::column_major> work_matrix_type;
	typedef typename work_matrix_type::array_type matrix_array_type;
	typedef typename work_vector_type::array_type vector_array_type;

	size_type n(ublasx::num_columns(A));
	size_type m(ublasx::num_columns(B));
	size_type nz(ublasx::size(z0));
//	size_type maxiter(niters);
	real_type alpha(0.01);
	real_type beta(0.9);
	real_type tol(0.1);
	size_type Tn(T*n);

	work_vector_type nu(n*T, 0);

	work_vector_type b(Tn, 0);
	ublas::subrange(b, 0, n) = ublas::prod(A, x);

	work_vector_type z(z0); // (nz x 1)

	work_matrix_type At(ublas::trans(A));
	work_matrix_type Bt(ublas::trans(B));
	work_matrix_type eyen(ublas::identity_matrix<real_type>(n, n));
	work_matrix_type eyem(ublas::identity_matrix<real_type>(m, m));

	DCS_DEBUG_TRACE(std::endl << "iteration \t step \t\t rd \t\t\t rp");

    for (size_type k = 0; k < niters; ++k)
    {
		work_vector_type gf(nz);
		work_vector_type gp(nz);
		work_vector_type hp(nz);
		work_vector_type rd(nz);
		work_vector_type rp(Tn);
		work_vector_type Ctnu(nz);
		real_type resd;
		real_type resp;
		real_type res;

		gfgphp(const_cast<matrix_array_type&>(Q().data()).begin(),
			   const_cast<matrix_array_type&>(R().data()).begin(),
			   const_cast<matrix_array_type&>(Qf().data()).begin(),
			   const_cast<vector_array_type&>(zmax().data()).begin(),
			   const_cast<vector_array_type&>(zmin().data()).begin(),
			   z.data().begin(),
			   T,
			   n,
			   m,
			   nz,
			   gf.data().begin(),
			   gp.data().begin(),
			   hp.data().begin());
		rdrp(const_cast<matrix_array_type&>(A().data()).begin(),
			 const_cast<matrix_array_type&>(B().data()).begin(),
//			 const_cast<matrix_array_type&>(Q().data()).begin(),
//			 const_cast<matrix_array_type&>(R().data()).begin(),
//			 const_cast<matrix_array_type&>(Qf().data()).begin(),
			 z.data().begin(),
			 nu.data().begin(),
			 gf.data().begin(),
			 gp.data().begin(),
			 b.data().begin(),
			 T,
			 n,
			 m,
			 nz,
			 kappa,
			 rd.data().begin(),
			 rp.data().begin(),
			 Ctnu.data().begin());
		resdresp(rd.data().begin(),
				 rp.data().begin(),
				 T,
				 n,
				 nz,
				 &resd,
				 &resp,
				 &res);

        if (res < tol) break;

		work_vector_type dz(nz);
		work_vector_type dnu(Tn);

		dnudz(const_cast<matrix_array_type&>(A().data()).begin(),
			  const_cast<matrix_array_type&>(B().data()).begin(),
			  At.data().begin(),
			  Bt.data().begin(),
			  eyen.data().begin(),
			  eyem.data().begin(),
			  const_cast<matrix_array_type&>(Q().data()).begin(),
			  const_cast<matrix_array_type&>(R().data()).begin(),
			  const_cast<matrix_array_type&>(Qf().data()).begin(),
			  hp.data().begin(),
			  rd.data().begin(),
			  rp.data().begin(),
			  T,
			  n,
			  m,
			  nz,
			  kappa,
			  dnu.data().begin(),
			  dz.data().begin()); 

    	real_type s(1); 

        /* feasibility search */
        while (1)
        {
            bool cont(false);
            for (size_type i = 0; i < nz && !cont; ++i)
            {
				if (z(i)+s*dz(i) >= zmax()(i)
					||
					z(i)+s*dz(i) <= zmin()(i)
				)
				{
					cont = true;
				}
            }
            if (cont)
            {
                s *= beta;
                continue;
            }
            else
			{
                break;
			}
        }

		work_vector_type newnu(nu+s*dnu);
		work_vector_type newz(z+s*dz);
		real_type newresd;
		real_type newresp;
		real_type newres;

        /* insert backtracking line search */
        while (1)
        {
			work_vector_type newgf(nz);
			work_vector_type newgp(nz);
			work_vector_type newhp(nz);
			work_vector_type newrd(nz);
			work_vector_type newrp(Tn);
			work_vector_type newCtnu(nz);

			gfgphp(const_cast<matrix_array_type&>(Q().data()).begin(),
				   const_cast<matrix_array_type&>(R().data()).begin(),
				   const_cast<matrix_array_type&>(Qf().data()).begin(),
				   const_cast<vector_array_type&>(zmax().data()).begin(),
				   const_cast<vector_array_type&>(zmin().data()).begin(),
				   newz.data().begin(),
				   T,
				   n,
				   m,
				   nz,
				   newgf.data().begin(),
				   newgp.data().begin(),
				   newhp.data().begin());
			rdrp(const_cast<matrix_array_type&>(A().data()).begin(),
				 const_cast<matrix_array_type&>(B().data()).begin(),
//				 const_cast<matrix_array_type&>(Q().data()).begin(),
//				 const_cast<matrix_array_type&>(R().data()).begin(),
//				 const_cast<matrix_array_type&>(Qf().data()).begin(),
				 newz.data().begin(),
				 newnu.data().begin(),
				 newgf.data().begin(),
				 newgp.data().begin(),
				 b.data().begin(),
				 T,
				 n,
				 m,
				 nz,
				 kappa,
				 newrd.data().begin(),
				 newrp.data().begin(),
				 newCtnu.data().begin());
			resdresp(newrd.data().begin(),
					 newrp.data().begin(),
					 T,
					 n,
					 nz,
					 &newresd,
					 &newresp,
					 &newres);

            if (newres <= (1-alpha*s)*res) break;

            s *= beta;

			newnu = nu+s*dnu;
			newz = z+s*dz;
        }

		nu = newnu;
		z = newz;

		DCS_DEBUG_TRACE("    " << k << " \t\t " << s << " \t " << newresd << " \t\t " << newresp);
    }

	z0() = z;
}


template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename XMinVectorT,
	typename XMaxVectorT,
	typename UMinVectorT,
	typename UMaxVectorT,
	typename QfMatrixT,
	typename SizeT,
	typename RealT,
	typename X0MatrixT,
	typename U0MatrixT,
	typename X0VectorT,
	typename XMatrixT,
	typename UMatrixT,
	typename XVectorT,
	typename KMatrixT
>
void fmpc_step(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			   ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			   ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
			   ::boost::numeric::ublas::vector_expression<XMinVectorT> const& xmin, 
			   ::boost::numeric::ublas::vector_expression<XMaxVectorT> const& xmax,
			   ::boost::numeric::ublas::vector_expression<UMinVectorT> const& umin, 
			   ::boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
			   ::boost::numeric::ublas::matrix_expression<QfMatrixT> const& Qf,
			   SizeT T,
			   SizeT niters,
			   RealT kappa,
			   ::boost::numeric::ublas::matrix_expression<X0MatrixT> const& X0,
			   ::boost::numeric::ublas::matrix_expression<U0MatrixT> const& U0,
			   ::boost::numeric::ublas::vector_expression<X0VectorT> const& x0,
			   bool want_X,
			   ::boost::numeric::ublas::matrix_container<XMatrixT>& X,
			   bool want_U,
			   ::boost::numeric::ublas::matrix_container<UMatrixT>& U,
			   bool want_x,
			   ::boost::numeric::ublas::vector_container<XVectorT>& x,
			   bool want_K,
			   ::boost::numeric::ublas::matrix_container<KMatrixT>& K)
{
	DCS_MACRO_SUPPRESS_UNUSED_VARIABLE_WARNING(want_x);

	namespace ublas = ::boost::numeric::ublas;
	namespace bindings = ::boost::numeric::bindings;

#if 0
	typedef SizeT size_type;
	//typedef typename ublas::promote_traits<SizeT,fortran_int_t>::promote_type size_type; //TODO: type promotion
#else
	typedef fortran_int_t size_type;
#endif
	typedef RealT real_type; //TODO: type promotion
	typedef ublas::matrix<real_type, ublas::column_major> work_matrix_type;
	typedef ublas::vector<real_type> work_vector_type;

#if 0
#else
	size_type TT(T);
	size_type nn(niters);
#endif // 0
	size_type nx(ublas::num_columns(A));
	size_type nu(ublas::num_columns(B));
    size_type nz(T*(nx+nu));

	work_vector_type z(nz, 0);
    for (size_type i = 0; i < TT; ++i)
    {
		size_type i1(i*(nu+nx));
		size_type i2(i1+nu);
		size_type i3(i2+nx);
		ublas::subrange(z, i1, i2) = ublas::column(U0(), i);
		ublas::subrange(z, i2, i3) = ublas::column(X0(), i);
    }

	work_vector_type zmin(nz, 0);
	work_vector_type zmax(nz, 0);
    for (size_type i = 0; i < TT; ++i)
    {
		size_type i1(i*(nu+nx));
		size_type i2(i1+nu);
		size_type i3(i2+nx);
		ublas::subrange(zmin, i1, i2) = umin;
		ublas::subrange(zmin, i2, i3) = xmin;
		ublas::subrange(zmax, i1, i2) = umax;
		ublas::subrange(zmax, i2, i3) = xmax;
    }

	work_vector_type zrang(zmax-zmin);
	work_vector_type zminp(zmin+0.01*zrang);
	work_vector_type zmaxp(zmax-0.01*zrang);


    /* project z */
    for (size_type i = 0; i < nz; ++i) z(i) = z(i) > zmaxp(i) ? zmaxp(i) : z(i);
    for (size_type i = 0; i < nz; ++i) z(i) = z(i) < zminp(i) ? zminp(i) : z(i);

	x() = x0;

#if 0
    fmpc_solve_impl(A, B, Q, R, Qf, zmin, zmax, x, z, T, niters, kappa);
#else
    fmpc_solve_impl(A, B, Q, R, Qf, zmin, zmax, x, z, TT, nn, kappa);
#endif // 0

	if (want_X || want_U)
	{
		for (size_type i = 0; i < TT; ++i)
		{
			size_type i1(i*(nu+nx));
			size_type i2(i1+nu);
			size_type i3(i2+nx);
			if (want_U)
			{
				ublas::column(U(), i) = ublas::subrange(z, i1, i2);
			}
			if (want_X)
			{
				ublas::column(X(), i) = ublas::subrange(z, i2, i3);
			}
		}
	}

	// Compute the terminal control gain K = (R+B'Q_fB)^{-1}B'Q_fA
	if (want_K)
	{
//		int info;

		work_matrix_type BtQf(ublas::prod(ublas::trans(B), Qf)); // B'Q_f
		K() = ublas::prod(BtQf, A); // B'Q_fA
		BtQf = ublas::prod(BtQf, B) + R; // B'Q_fB+R
//		F77_CALL(dposv)("l",&nu,&nx,BtQf.data().begin(),&nu,K().data().begin(),&nu,&info); // Find K s.t. (B'Q_fB+R)K=B'Q_fA
		bindings::lapack::posv(bindings::lower(BtQf), K());
	}
}

}} // Namespace detail::<unnamed>


template <typename RealT, typename SizeT>
class fmpc_controller
{
	public: typedef RealT real_type;
	public: typedef SizeT size_type;
	public: typedef ::boost::numeric::ublas::matrix<real_type, ::boost::numeric::ublas::column_major> matrix_type;
	public: typedef ::boost::numeric::ublas::vector<real_type> vector_type;


	public: fmpc_controller()
	{
	}


	public: template <
				typename QMatrixT,
				typename RMatrixT,
				typename QfMatrixT,
				typename XMinVectorT,
				typename XMaxVectorT,
				typename UMinVectorT,
				typename UMaxVectorT>
		fmpc_controller(::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
						::boost::numeric::ublas::matrix_expression<QfMatrixT> const& Qf,
						::boost::numeric::ublas::vector_expression<XMinVectorT> const& xmin,
						::boost::numeric::ublas::vector_expression<XMaxVectorT> const& xmax,
						::boost::numeric::ublas::vector_expression<UMinVectorT> const& umin,
						::boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
						size_type horiz,
						real_type barrier,
						size_type niters)
	: Q_(Q),
	  R_(R),
	  Qf_(Qf),
	  xmin_(xmin),
	  xmax_(xmax),
	  umin_(umin),
	  umax_(umax),
	  T_(horiz),
	  kappa_(barrier),
	  niters_(niters),
	  X_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(Q), horiz)),
	  U_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(R), horiz)),
	  K_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(R), ::boost::numeric::ublasx::num_rows(Q))),
	  J_(0)
	{
	}


	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename QMatrixT,
				typename RMatrixT,
				typename XMinVectorT,
				typename XMaxVectorT,
				typename UMinVectorT,
				typename UMaxVectorT,
				typename QfMatrixT,
				typename X0MatrixT,
				typename U0MatrixT>
		fmpc_controller(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
						::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
						::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
						::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
						::boost::numeric::ublas::matrix_expression<QfMatrixT> const& Qf,
						::boost::numeric::ublas::vector_expression<XMinVectorT> const& xmin,
						::boost::numeric::ublas::vector_expression<XMaxVectorT> const& xmax,
						::boost::numeric::ublas::vector_expression<UMinVectorT> const& umin,
						::boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
						size_type horiz,
						real_type barrier,
						size_type niters,
						::boost::numeric::ublas::matrix_expression<X0MatrixT> const& X0,
						::boost::numeric::ublas::matrix_expression<U0MatrixT> const& U0)
	: Q_(Q),
	  R_(R),
	  Qf_(Qf),
	  xmin_(xmin),
	  xmax_(xmax),
	  umin_(umin),
	  umax_(umax),
	  T_(horiz),
	  kappa_(barrier),
	  niters_(niters),
	  X_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(Q), horiz)),
	  U_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(R), horiz)),
	  K_(::boost::numeric::ublas::zero_matrix<real_type>(::boost::numeric::ublasx::num_rows(R), ::boost::numeric::ublasx::num_rows(Q))),
	  J_(0)
	{
		solve(A, B, X0, U0/*, x0*/);
	}


	public: template <
				typename AMatrixT,
				typename BMatrixT,
				typename X0MatrixT,
				typename U0MatrixT>
		void solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
					 ::boost::numeric::ublas::matrix_expression<X0MatrixT> const& X0,
					 ::boost::numeric::ublas::matrix_expression<U0MatrixT> const& U0)
	{
		A_ = A;
		B_ = B;
		X_ = X0;
		U_ = U0;
//		x_ = x0;
	}


	public: template <
				typename AMatrixT,
				typename BMatrixT>
		void solve(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
					 ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B)
	{
		A_ = A;
		B_ = B;
	}


	public: template <typename XVectorT>
		vector_type control(::boost::numeric::ublas::vector_expression<XVectorT> const& x)
	{
		namespace ublas = ::boost::numeric::ublas;
		namespace ublasx = ::boost::numeric::ublasx;

		vector_type dummy_x;

		detail::fmpc_step(A_, B_, Q_, R_, xmin_, xmax_, umin_, umax_, Qf_, T_, niters_, kappa_, X_, U_, x, true, X_, true, U_, false, dummy_x, true, K_);

		// Get the next optimal action
		vector_type u(ublas::column(U_, 0));

		// compute stage cost
        vector_type xtQ(ublas::prod(x, Q_));
        vector_type utR(ublas::prod(u, R_));
        J_ = ublas::inner_prod(xtQ, x) + ublas::inner_prod(utR, u);

		// shift previous state and input trajectories
		// for warm start in next step
		X_ = ublasx::cat_rows(ublas::subrange(X_, 0, ublasx::size(x), 1, T_), ublas::zero_matrix<real_type>(ublasx::size(x), 1));
		U_ = ublasx::cat_rows(ublas::subrange(U_, 0, ublasx::size(u), 1, T_), ublas::zero_matrix<real_type>(ublasx::size(u), 1));

		return u;
	}


	public: matrix_type Q() const
	{
		return Q_;
	}


	public: matrix_type R() const
	{
		return R_;
	}


	public: matrix_type Qf() const
	{
		return Qf_;
	}


	public: vector_type xmin() const
	{
		return xmin_;
	}


	public: vector_type xmax() const
	{
		return xmax_;
	}


	public: vector_type umin() const
	{
		return umin_;
	}


	public: vector_type umax() const
	{
		return umax_;
	}


	public: size_type horizon() const
	{
		return T_;
	}


	public: size_type max_num_iterations() const
	{
		return niters_;
	}


	public: real_type barrier() const
	{
		return kappa_;
	}


	public: matrix_type predicted_states() const
	{
		return X_;
	}


	public: matrix_type optimal_inputs() const
	{
		return U_;
	}


	public: matrix_type gain() const
	{
		return K_;
	}


	private: real_type cost() const
	{
		return J_;
	}


	private: matrix_type Q_;
	private: matrix_type R_;
	private: matrix_type Qf_;
	private: vector_type xmin_;
	private: vector_type xmax_;
	private: vector_type umin_;
	private: vector_type umax_;
	private: size_type T_;
	private: real_type kappa_;
	private: size_type niters_;
	private: matrix_type A_;
	private: matrix_type B_;
	private: matrix_type X_;
	private: matrix_type U_;
//	private: matrix_type x_;
	private: matrix_type K_;
	private: real_type J_;
}; // fmpc_controller


template <
	typename AMatrixT,
	typename BMatrixT,
	typename QMatrixT,
	typename RMatrixT,
	typename XMinVectorT,
	typename XMaxVectorT,
	typename UMinVectorT,
	typename UMaxVectorT,
	typename QfMatrixT,
	typename SizeT,
	typename RealT,
	typename X0MatrixT,
	typename U0MatrixT,
	typename X0VectorT
>
inline
::boost::numeric::ublas::vector<RealT> fmpc_step(::boost::numeric::ublas::matrix_expression<AMatrixT> const& A,
			   ::boost::numeric::ublas::matrix_expression<BMatrixT> const& B,
			   ::boost::numeric::ublas::matrix_expression<QMatrixT> const& Q,
			   ::boost::numeric::ublas::matrix_expression<RMatrixT> const& R,
			   ::boost::numeric::ublas::matrix_expression<QfMatrixT> const& Qf,
			   ::boost::numeric::ublas::vector_expression<XMinVectorT> const& xmin, 
			   ::boost::numeric::ublas::vector_expression<XMaxVectorT> const& xmax,
			   ::boost::numeric::ublas::vector_expression<UMinVectorT> const& umin, 
			   ::boost::numeric::ublas::vector_expression<UMaxVectorT> const& umax,
			   SizeT T,
			   RealT kappa,
			   SizeT niters,
			   ::boost::numeric::ublas::matrix_container<X0MatrixT>& X0,
			   ::boost::numeric::ublas::matrix_container<U0MatrixT>& U0,
			   ::boost::numeric::ublas::vector_container<X0VectorT>& x0)
{
	namespace ublas = ::boost::numeric::ublas;

	typedef RealT real_type;
	typedef SizeT size_type;
	typedef ublas::matrix<real_type, ublas::column_major> work_matrix_type;
	typedef ublas::vector<real_type, ublas::column_major> work_vector_type;

	work_matrix_type dummyK;

	detail::fmpc_step(A, B, Q, R, xmin, xmax, umin, umax, Qf, T, niters, kappa, X0, U0, x0, true, X0, true, U0, true, x0, false, dummyK);

	return ublas::column(U0(), 0); 
}
}} // Namespace dcs::control


#endif // DCS_CONTROL_FMPC_HPP
