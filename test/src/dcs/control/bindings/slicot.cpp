/**
 * \file test/src/dcs/control/bindings/slicot.cpp
 *
 * \brief Test suite for bindings to the SLICOT Fortran library.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright (C) 2012       Marco Guazzone (marco.guazzone@gmail.com)
 *                          [Distributed Computing System (DCS) Group,
 *                           Computer Science Institute,
 *                           Department of Science and Technological Innovation,
 *                           University of Piemonte Orientale,
 *                           Alessandria (Italy)]
 *
 * This file is part of Foobar (below referred to as "this program").
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cstddef>
#include <dcs/assert.hpp>
#include <dcs/debug.hpp>
#include <dcs/exception.hpp>
#include <dcs/math/type/matrix.hpp>
#include <dcs/math/type/matrix_properties.hpp>
#include <dcs/math/type/vector.hpp>
#include <dcs/control/bindings/fortran.hpp>
#include <dcs/control/bindings/slicot.hpp>
#include <dcs/test.hpp>
#include <stdexcept>

namespace /*<unnamed>*/ {

const double tol(1.0e-5);

} // Namespace <unnamed>


DCS_TEST_DEF( sb02od_c_arrays )
{
	DCS_TEST_CASE("SB02OD - C Arrays");

	namespace math = dcs::math;
	namespace bindings = dcs::control::bindings;
	namespace slicot = dcs::control::bindings::slicot;

	typedef double real_type;
	typedef math::matrix<real_type,math::matrix_properties<math::column_major_storage_layout> > matrix_type;
	typedef math::vector<real_type> real_vector_type;
	typedef math::vector<bindings::fortran_int> int_vector_type;
	typedef math::vector<bindings::fortran_logical> logical_vector_type;

	const bindings::fortran_int nx(3);
	const bindings::fortran_int nu(1);

	const bindings::fortran_int lda = std::max(bindings::fortran_int(1),nx);
	matrix_type A(lda,nx);
	A(0,0) = 0;        A(0,1) =  0;        A(0,2) =  1;
	A(1,0) = 0;        A(1,1) =  0;        A(1,2) =  0;
	A(2,0) = 0.243365; A(2,1) = -0.563735; A(2,2) = -0.932137;

	const bindings::fortran_int ldb = std::max(bindings::fortran_int(1),nx);
	matrix_type B(ldb,nu);
	B(0,0) = 1;
	B(1,0) = 0;
	B(2,0) = 0.811669;

	const bindings::fortran_int ldq = std::max(bindings::fortran_int(1),nx);
	matrix_type Q(ldq,nx,0);
	Q(nx-1,nx-1) = 1;

	const bindings::fortran_int ldr = std::max(bindings::fortran_int(1),nu);
	matrix_type R(ldr,nu,1);

	const bindings::fortran_int ldl(1);
	matrix_type L(ldl,nx,0);

	double rcond(0);

	const bindings::fortran_int ldx = std::max(bindings::fortran_int(1),nx);
	matrix_type X(ldx,nx,0);
	real_vector_type alfar(2*nx);
	real_vector_type alfai(2*nx);
	real_vector_type beta(2*nx);
	const bindings::fortran_int lds = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type S(lds,2*nx+nu,0);
	const bindings::fortran_int ldt = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type T(ldt,2*nx,0);
	const bindings::fortran_int ldu = std::max(bindings::fortran_int(1),2*nx);
	matrix_type U(ldu,2*nx,0);
	const bindings::fortran_int liwork = std::max(std::max(bindings::fortran_int(1),nu),2*nx);
	int_vector_type iwork(liwork,0);
	const bindings::fortran_int ldwork = std::max(std::max(std::max(7*(2*nx+1)+16,16*nx),2*nx+nu),3*nu);
	real_vector_type dwork(ldwork);
	logical_vector_type bwork(2*nx, 0);

	const char dico('D');
	const char jobb('B');
	const char fact('N');
	const char uplo('U');
	const char jobl('Z');
	const char sort('S');

	matrix_type ok_X(nx,nx);
	ok_X(0,0) =  0.0671891998477877; ok_X(0,1) = -0.1556382535540959; ok_X(0,2) = -0.2977628887069395;
	ok_X(1,0) = -0.1556382535540959; ok_X(1,1) =  0.3605232094480238; ok_X(1,2) =  0.6897432336827665;
	ok_X(2,0) = -0.2977628887069395; ok_X(2,1) =  0.6897432336827665; ok_X(2,2) =  2.3481878868207247;

	slicot::sb02od(dico,
				   jobb,
				   fact,
				   uplo,
				   jobl,
				   sort,
				   nx,
				   nu,
				   0,
				   A.begin_data(),
				   lda,
				   B.begin_data(),
				   ldb,
				   Q.begin_data(),
				   ldq,
				   R.begin_data(),
				   ldr,
				   L.begin_data(),
				   ldl,
				   rcond,
				   X.begin_data(),
				   ldx,
				   alfar.begin_data(),
				   alfai.begin_data(),
				   beta.begin_data(),
				   S.begin_data(),
				   lds,
				   T.begin_data(),
				   ldt,
				   U.begin_data(),
				   ldu,
				   0,
				   iwork.begin_data(),
				   dwork.begin_data(),
				   ldwork,
				   bwork.begin_data());

	DCS_TEST_CHECK_VECTOR_CLOSE(X.begin_data(), ok_X.begin_data(), nx*nx, tol);
}

DCS_TEST_DEF( sb02od_ublas_arrays )
{
	DCS_TEST_CASE("SB02OD - Boost.uBLAS Arrays");

	namespace bindings = dcs::control::bindings;
	namespace slicot = dcs::control::bindings::slicot;
	namespace ublas = boost::numeric::ublas;

	typedef double real_type;
	typedef ublas::matrix<real_type,ublas::column_major> matrix_type;
	typedef ublas::vector<real_type> real_vector_type;
	typedef ublas::vector<bindings::fortran_int> int_vector_type;
	typedef ublas::vector<bindings::fortran_logical> logical_vector_type;

	const bindings::fortran_int nx(3);
	const bindings::fortran_int nu(1);

	const bindings::fortran_int lda = std::max(bindings::fortran_int(1),nx);
	matrix_type A(lda,nx);
	A(0,0) = 0;        A(0,1) =  0;        A(0,2) =  1;
	A(1,0) = 0;        A(1,1) =  0;        A(1,2) =  0;
	A(2,0) = 0.243365; A(2,1) = -0.563735; A(2,2) = -0.932137;

	const bindings::fortran_int ldb = std::max(bindings::fortran_int(1),nx);
	matrix_type B(ldb,nu);
	B(0,0) = 1;
	B(1,0) = 0;
	B(2,0) = 0.811669;

	const bindings::fortran_int ldq = std::max(bindings::fortran_int(1),nx);
	matrix_type Q(ldq,nx,0);
	Q(nx-1,nx-1) = 1;

	const bindings::fortran_int ldr = std::max(bindings::fortran_int(1),nu);
	matrix_type R(ldr,nu,1);

	const bindings::fortran_int ldl(1);
	matrix_type L(ldl,nx,0);

	double rcond(0);

	const bindings::fortran_int ldx = std::max(bindings::fortran_int(1),nx);
	matrix_type X(ldx,nx,0);
	real_vector_type alfar(2*nx);
	real_vector_type alfai(2*nx);
	real_vector_type beta(2*nx);
	const bindings::fortran_int lds = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type S(lds,2*nx+nu,0);
	const bindings::fortran_int ldt = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type T(ldt,2*nx,0);
	const bindings::fortran_int ldu = std::max(bindings::fortran_int(1),2*nx);
	matrix_type U(ldu,2*nx,0);
	const bindings::fortran_int liwork = std::max(std::max(bindings::fortran_int(1),nu),2*nx);
	int_vector_type iwork(liwork,0);
	const bindings::fortran_int ldwork = std::max(std::max(std::max(7*(2*nx+1)+16,16*nx),2*nx+nu),3*nu);
	real_vector_type dwork(ldwork);
	logical_vector_type bwork(2*nx, 0);

	const char dico('D');
	const char jobb('B');
	const char fact('N');
	const char uplo('U');
	const char jobl('Z');
	const char sort('S');

	matrix_type ok_X(nx,nx);
	ok_X(0,0) =  0.0671891998477877; ok_X(0,1) = -0.1556382535540959; ok_X(0,2) = -0.2977628887069395;
	ok_X(1,0) = -0.1556382535540959; ok_X(1,1) =  0.3605232094480238; ok_X(1,2) =  0.6897432336827665;
	ok_X(2,0) = -0.2977628887069395; ok_X(2,1) =  0.6897432336827665; ok_X(2,2) =  2.3481878868207247;

	slicot::sb02od(dico,
				   jobb,
				   fact,
				   uplo,
				   jobl,
				   sort,
				   nx,
				   nu,
				   0,
				   A.data().begin(),
				   lda,
				   B.data().begin(),
				   ldb,
				   Q.data().begin(),
				   ldq,
				   R.data().begin(),
				   ldr,
				   L.data().begin(),
				   ldl,
				   rcond,
				   X.data().begin(),
				   ldx,
				   alfar.data().begin(),
				   alfai.data().begin(),
				   beta.data().begin(),
				   S.data().begin(),
				   lds,
				   T.data().begin(),
				   ldt,
				   U.data().begin(),
				   ldu,
				   0,
				   iwork.data().begin(),
				   dwork.data().begin(),
				   ldwork,
				   bwork.data().begin());

	DCS_TEST_CHECK_VECTOR_CLOSE(X.data().begin(), ok_X.data().begin(), nx*nx, tol);
}

DCS_TEST_DEF( sb02od_ublas_arrays_2 )
{
	DCS_TEST_CASE("SB02OD - Boost.uBLAS Arrays #2");

	namespace bindings = dcs::control::bindings;
	namespace slicot = dcs::control::bindings::slicot;
	namespace ublas = boost::numeric::ublas;

	typedef double real_type;
	typedef ublas::matrix<real_type,ublas::column_major> matrix_type;
	typedef ublas::vector<real_type> real_vector_type;
	typedef ublas::vector<bindings::fortran_int> int_vector_type;
	typedef ublas::vector<bindings::fortran_logical> logical_vector_type;
	typedef ublas::vector< std::complex<real_type> > complex_vector_type;

	const bindings::fortran_int nx(3);
	const bindings::fortran_int nu(1);

	const bindings::fortran_int lda = std::max(bindings::fortran_int(1),nx);
	matrix_type A(lda,nx);
	A(0,0) = 0;        A(0,1) =  0;        A(0,2) =  1;
	A(1,0) = 0;        A(1,1) =  0;        A(1,2) =  0;
	A(2,0) = 0.243365; A(2,1) = -0.563735; A(2,2) = -0.932137;

	const bindings::fortran_int ldb = std::max(bindings::fortran_int(1),nx);
	matrix_type B(ldb,nu);
	B(0,0) = 1;
	B(1,0) = 0;
	B(2,0) = 0.811669;

	const bindings::fortran_int ldq = std::max(bindings::fortran_int(1),nx);
	matrix_type Q(ldq,nx,0);
	Q(nx-1,nx-1) = 1;

	const bindings::fortran_int ldr = std::max(bindings::fortran_int(1),nu);
	matrix_type R(ldr,nu,1);

	matrix_type L;

	//const bindings::fortran_int ldx = std::max(bindings::fortran_int(1),nx);
	matrix_type X;
	complex_vector_type lambda;
	//const bindings::fortran_int lds = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type S;
	//const bindings::fortran_int ldt = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type T;
	//const bindings::fortran_int ldu = std::max(bindings::fortran_int(1),2*nx);
	matrix_type U;

	matrix_type ok_X(nx,nx);
	ok_X(0,0) =  0.0671891998477877; ok_X(0,1) = -0.1556382535540959; ok_X(0,2) = -0.2977628887069395;
	ok_X(1,0) = -0.1556382535540959; ok_X(1,1) =  0.3605232094480238; ok_X(1,2) =  0.6897432336827665;
	ok_X(2,0) = -0.2977628887069395; ok_X(2,1) =  0.6897432336827665; ok_X(2,2) =  2.3481878868207247;

	slicot::sb02od(A, B, Q, R, L, true, X, lambda, S, T, U);

	//DCS_DEBUG_TRACE("X=" << X);
	//DCS_DEBUG_TRACE("lambda=" << lambda);
	//DCS_DEBUG_TRACE("S=" << S);
	//DCS_DEBUG_TRACE("T=" << T);
	//DCS_DEBUG_TRACE("U=" << U);

	DCS_TEST_CHECK_VECTOR_CLOSE(X.data().begin(), ok_X.data().begin(), nx*nx, tol);
}

DCS_TEST_DEF( sg02ad_c_arrays )
{
	DCS_TEST_CASE("SG02AD - C Arrays");

	namespace math = dcs::math;
	namespace bindings = dcs::control::bindings;
	namespace slicot = dcs::control::bindings::slicot;

	typedef double real_type;
	typedef math::matrix<real_type,math::matrix_properties<math::column_major_storage_layout> > matrix_type;
	typedef math::vector<real_type> real_vector_type;
	typedef math::vector<bindings::fortran_int> int_vector_type;
	typedef math::vector<bindings::fortran_logical> logical_vector_type;

	const bindings::fortran_int nx(3);
	const bindings::fortran_int nu(1);

	const bindings::fortran_int lda = std::max(bindings::fortran_int(1),nx);
	matrix_type A(lda,nx);
	A(0,0) = 0;        A(0,1) =  0;        A(0,2) =  1;
	A(1,0) = 0;        A(1,1) =  0;        A(1,2) =  0;
	A(2,0) = 0.243365; A(2,1) = -0.563735; A(2,2) = -0.932137;

	const bindings::fortran_int lde = std::max(bindings::fortran_int(1),nx);
	matrix_type E(lde,nx);
	E(0,0) = 1; E(0,1) = 0; E(0,2) = 0;
	E(1,0) = 0; E(1,1) = 1; E(1,2) = 0;
	E(2,0) = 0; E(2,1) = 0; E(2,2) = 1;

	const bindings::fortran_int ldb = std::max(bindings::fortran_int(1),nx);
	matrix_type B(ldb,nu);
	B(0,0) = 1;
	B(1,0) = 0;
	B(2,0) = 0.811669;

	const bindings::fortran_int ldq = std::max(bindings::fortran_int(1),nx);
	matrix_type Q(ldq,nx,0);
	Q(nx-1,nx-1) = 1;

	const bindings::fortran_int ldr = std::max(bindings::fortran_int(1),nu);
	matrix_type R(ldr,nu,1);

	const bindings::fortran_int ldl(1);
	matrix_type L(ldl,nx,0);

	double rcondu(0);

	const bindings::fortran_int ldx = std::max(bindings::fortran_int(1),nx);
	matrix_type X(ldx,nx,0);
	real_vector_type alfar(2*nx);
	real_vector_type alfai(2*nx);
	real_vector_type beta(2*nx);
	const bindings::fortran_int lds = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type S(lds,2*nx+nu,0);
	const bindings::fortran_int ldt = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type T(ldt,2*nx,0);
	const bindings::fortran_int ldu = std::max(bindings::fortran_int(1),2*nx);
	matrix_type U(ldu,2*nx,0);
	const bindings::fortran_int liwork = std::max(std::max(bindings::fortran_int(1),nu),2*nx);
	int_vector_type iwork(liwork,0);
	const bindings::fortran_int ldwork = std::max(std::max(std::max(7*(2*nx+1)+16,16*nx),2*nx+nu),3*nu);
	real_vector_type dwork(ldwork);
	logical_vector_type bwork(2*nx, 0);

	const char dico('D');
	const char jobb('B');
	const char fact('N');
	const char uplo('U');
	const char jobl('Z');
	const char scal('N');
	const char sort('S');
	const char acc('N');

	matrix_type ok_X(nx,nx);
	ok_X(0,0) =  0.0671891998477877; ok_X(0,1) = -0.1556382535540959; ok_X(0,2) = -0.2977628887069395;
	ok_X(1,0) = -0.1556382535540959; ok_X(1,1) =  0.3605232094480238; ok_X(1,2) =  0.6897432336827665;
	ok_X(2,0) = -0.2977628887069395; ok_X(2,1) =  0.6897432336827665; ok_X(2,2) =  2.3481878868207247;

	bindings::fortran_int iwarn(0);

	slicot::sg02ad(dico,
				   jobb,
				   fact,
				   uplo,
				   jobl,
				   scal,
				   sort,
				   acc,
				   nx,
				   nu,
				   0,
				   A.begin_data(),
				   lda,
				   E.begin_data(),
				   lde,
				   B.begin_data(),
				   ldb,
				   Q.begin_data(),
				   ldq,
				   R.begin_data(),
				   ldr,
				   L.begin_data(),
				   ldl,
				   rcondu,
				   X.begin_data(),
				   ldx,
				   alfar.begin_data(),
				   alfai.begin_data(),
				   beta.begin_data(),
				   S.begin_data(),
				   lds,
				   T.begin_data(),
				   ldt,
				   U.begin_data(),
				   ldu,
				   0,
				   iwork.begin_data(),
				   dwork.begin_data(),
				   ldwork,
				   bwork.begin_data(),
				   iwarn);

	DCS_TEST_CHECK_VECTOR_CLOSE(X.begin_data(), ok_X.begin_data(), nx*nx, tol);
}

DCS_TEST_DEF( sg02ad_ublas_arrays )
{
	DCS_TEST_CASE("SG02AD - Boost.uBLAS Arrays");

	namespace bindings = dcs::control::bindings;
	namespace slicot = dcs::control::bindings::slicot;
	namespace ublas = boost::numeric::ublas;

	typedef double real_type;
	typedef ublas::matrix<real_type,ublas::column_major> matrix_type;
	typedef ublas::vector<real_type> real_vector_type;
	typedef ublas::vector<bindings::fortran_int> int_vector_type;
	typedef ublas::vector<bindings::fortran_logical> logical_vector_type;

	const bindings::fortran_int nx(3);
	const bindings::fortran_int nu(1);

	const bindings::fortran_int lda = std::max(bindings::fortran_int(1),nx);
	matrix_type A(lda,nx);
	A(0,0) = 0;        A(0,1) =  0;        A(0,2) =  1;
	A(1,0) = 0;        A(1,1) =  0;        A(1,2) =  0;
	A(2,0) = 0.243365; A(2,1) = -0.563735; A(2,2) = -0.932137;

	const bindings::fortran_int lde = std::max(bindings::fortran_int(1),nx);
	matrix_type E(lde,nx);
	E(0,0) = 1; E(0,1) = 0; E(0,2) = 0;
	E(1,0) = 0; E(1,1) = 1; E(1,2) = 0;
	E(2,0) = 0; E(2,1) = 0; E(2,2) = 1;

	const bindings::fortran_int ldb = std::max(bindings::fortran_int(1),nx);
	matrix_type B(ldb,nu);
	B(0,0) = 1;
	B(1,0) = 0;
	B(2,0) = 0.811669;

	const bindings::fortran_int ldq = std::max(bindings::fortran_int(1),nx);
	matrix_type Q(ldq,nx,0);
	Q(nx-1,nx-1) = 1;

	const bindings::fortran_int ldr = std::max(bindings::fortran_int(1),nu);
	matrix_type R(ldr,nu,1);

	const bindings::fortran_int ldl(1);
	matrix_type L(ldl,nx,0);

	double rcondu(0);

	const bindings::fortran_int ldx = std::max(bindings::fortran_int(1),nx);
	matrix_type X(ldx,nx,0);
	real_vector_type alfar(2*nx);
	real_vector_type alfai(2*nx);
	real_vector_type beta(2*nx);
	const bindings::fortran_int lds = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type S(lds,2*nx+nu,0);
	const bindings::fortran_int ldt = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type T(ldt,2*nx,0);
	const bindings::fortran_int ldu = std::max(bindings::fortran_int(1),2*nx);
	matrix_type U(ldu,2*nx,0);
	const bindings::fortran_int liwork = std::max(std::max(bindings::fortran_int(1),nu),2*nx);
	int_vector_type iwork(liwork,0);
	const bindings::fortran_int ldwork = std::max(std::max(std::max(7*(2*nx+1)+16,16*nx),2*nx+nu),3*nu);
	real_vector_type dwork(ldwork);
	logical_vector_type bwork(2*nx, 0);

	const char dico('D');
	const char jobb('B');
	const char fact('N');
	const char uplo('U');
	const char jobl('Z');
	const char scal('N');
	const char sort('S');
	const char acc('N');

	matrix_type ok_X(nx,nx);
	ok_X(0,0) =  0.0671891998477877; ok_X(0,1) = -0.1556382535540959; ok_X(0,2) = -0.2977628887069395;
	ok_X(1,0) = -0.1556382535540959; ok_X(1,1) =  0.3605232094480238; ok_X(1,2) =  0.6897432336827665;
	ok_X(2,0) = -0.2977628887069395; ok_X(2,1) =  0.6897432336827665; ok_X(2,2) =  2.3481878868207247;

	bindings::fortran_int iwarn(0);

	slicot::sg02ad(dico,
				   jobb,
				   fact,
				   uplo,
				   jobl,
				   scal,
				   sort,
				   acc,
				   nx,
				   nu,
				   0,
				   A.data().begin(),
				   lda,
				   E.data().begin(),
				   lde,
				   B.data().begin(),
				   ldb,
				   Q.data().begin(),
				   ldq,
				   R.data().begin(),
				   ldr,
				   L.data().begin(),
				   ldl,
				   rcondu,
				   X.data().begin(),
				   ldx,
				   alfar.data().begin(),
				   alfai.data().begin(),
				   beta.data().begin(),
				   S.data().begin(),
				   lds,
				   T.data().begin(),
				   ldt,
				   U.data().begin(),
				   ldu,
				   0,
				   iwork.data().begin(),
				   dwork.data().begin(),
				   ldwork,
				   bwork.data().begin(),
				   iwarn);

	DCS_TEST_CHECK_VECTOR_CLOSE(X.data().begin(), ok_X.data().begin(), nx*nx, tol);
}

DCS_TEST_DEF( sg02ad_ublas_arrays_2 )
{
	DCS_TEST_CASE("SG02AD - Boost.uBLAS Arrays #2");

	namespace bindings = dcs::control::bindings;
	namespace slicot = dcs::control::bindings::slicot;
	namespace ublas = boost::numeric::ublas;

	typedef double real_type;
	typedef ublas::matrix<real_type,ublas::column_major> matrix_type;
	typedef ublas::vector<real_type> real_vector_type;
	typedef ublas::vector<bindings::fortran_int> int_vector_type;
	typedef ublas::vector<bindings::fortran_logical> logical_vector_type;
	typedef ublas::vector< std::complex<real_type> > complex_vector_type;

	const bindings::fortran_int nx(3);
	const bindings::fortran_int nu(1);

	const bindings::fortran_int lda = std::max(bindings::fortran_int(1),nx);
	matrix_type A(lda,nx);
	A(0,0) = 0;        A(0,1) =  0;        A(0,2) =  1;
	A(1,0) = 0;        A(1,1) =  0;        A(1,2) =  0;
	A(2,0) = 0.243365; A(2,1) = -0.563735; A(2,2) = -0.932137;

	const bindings::fortran_int lde = std::max(bindings::fortran_int(1),nx);
	matrix_type E(lde,nx);
	E(0,0) = 1; E(0,1) = 0; E(0,2) = 0;
	E(1,0) = 0; E(1,1) = 1; E(1,2) = 0;
	E(2,0) = 0; E(2,1) = 0; E(2,2) = 1;

	const bindings::fortran_int ldb = std::max(bindings::fortran_int(1),nx);
	matrix_type B(ldb,nu);
	B(0,0) = 1;
	B(1,0) = 0;
	B(2,0) = 0.811669;

	const bindings::fortran_int ldq = std::max(bindings::fortran_int(1),nx);
	matrix_type Q(ldq,nx,0);
	Q(nx-1,nx-1) = 1;

	const bindings::fortran_int ldr = std::max(bindings::fortran_int(1),nu);
	matrix_type R(ldr,nu,1);

	matrix_type L;

	//const bindings::fortran_int ldx = std::max(bindings::fortran_int(1),nx);
	matrix_type X;
	complex_vector_type lambda;
	//const bindings::fortran_int lds = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type S;
	//const bindings::fortran_int ldt = std::max(bindings::fortran_int(1),2*nx+nu);
	matrix_type T;
	//const bindings::fortran_int ldu = std::max(bindings::fortran_int(1),2*nx);
	matrix_type U;

	matrix_type ok_X(nx,nx);
	ok_X(0,0) =  0.0671891998477877; ok_X(0,1) = -0.1556382535540959; ok_X(0,2) = -0.2977628887069395;
	ok_X(1,0) = -0.1556382535540959; ok_X(1,1) =  0.3605232094480238; ok_X(1,2) =  0.6897432336827665;
	ok_X(2,0) = -0.2977628887069395; ok_X(2,1) =  0.6897432336827665; ok_X(2,2) =  2.3481878868207247;

	slicot::sg02ad(A, E, B, Q, R, L, true, X, lambda, S, T, U);

	//DCS_DEBUG_TRACE("X=" << X);
	//DCS_DEBUG_TRACE("lambda=" << lambda);
	//DCS_DEBUG_TRACE("S=" << S);
	//DCS_DEBUG_TRACE("T=" << T);
	//DCS_DEBUG_TRACE("U=" << U);

	DCS_TEST_CHECK_VECTOR_CLOSE(X.data().begin(), ok_X.data().begin(), nx*nx, tol);
}

int main()
{
	DCS_TEST_SUITE("DCS Control SLICOT C++ Interface");

	DCS_TEST_BEGIN();
		DCS_TEST_DO( sb02od_c_arrays );
		DCS_TEST_DO( sb02od_ublas_arrays );
		DCS_TEST_DO( sb02od_ublas_arrays_2 );
		DCS_TEST_DO( sg02ad_c_arrays );
		DCS_TEST_DO( sg02ad_ublas_arrays );
		DCS_TEST_DO( sg02ad_ublas_arrays_2 );
	DCS_TEST_END();
}
