/**
 * \file test/src/dcs/control/dlqi.cpp
 *
 * \brief Test suite for DLQI controllers.
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
 * This file is part of dcsxx-control (below referred to as "this program").
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

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <cstddef>
#include <dcs/control/design/dlqi.hpp>
#include <dcs/test.hpp>
#include <iostream>
#include "./utility.hpp"


namespace ublas = boost::numeric::ublas;
namespace dcs_ctrl = dcs::control;


const double tol = 1.0e-5;


DCS_TEST_DEF( matlab_1 )
{
	DCS_TEST_CASE("MATLAB #1");

	typedef double real_type;
	typedef real_type value_type;
    typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
    typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t nx(2);
	const std::size_t nu(1);
	const std::size_t ny(1);
	const std::size_t nz(nx+1);
	const real_type ts(1);

	matrix_type A(nx,nx);
	A(0,0) =   0; A(0,1) =   1;
	A(1,0) = -20; A(1,1) = -10;

	matrix_type B(nx,nu);
	B(0,0) = 0;
	B(1,0) = 1;

	matrix_type C(ny,nx);
	C(0,0) = 1; C(0,1) = 0;

	matrix_type D(ny,nu);
	D(0,0) = 0;

	matrix_type Q(ublas::identity_matrix<value_type>(nz,nz));
	Q(nx,nx) = 2e+6;

	matrix_type R(ublas::identity_matrix<value_type>(nu,nu));

	matrix_type N(ublas::zero_matrix<value_type>(nz,nu));


	matrix_type expect_K(nu,nz); // state-feedback optimal gain
	expect_K(0,0) = -19.000136465758239; expect_K(0,1) = -9.000051490140352; expect_K(0,2) = -0.999853537170799;

	matrix_type expect_S(nz,nz); // solution of the associated Riccati equation
	expect_S(0,0) = 1e+6* 4.000863944497572; expect_S(0,1) = 1e+6* 2.000472969596359; expect_S(0,2) = 1e+6*-4.000482941768265;
	expect_S(1,0) = 1e+6* 2.000472969596359; expect_S(1,1) = 1e+6* 2.000584980049262; expect_S(1,2) = 1e+6*-2.000292968566564;
	expect_S(2,0) = 1e+6*-4.000482941768265; expect_S(2,1) = 1e+6*-2.000292968566564; expect_S(2,2) = 1e+6* 6.000502938839015;

	vector_type expect_l(nz); // closed-loop eigenvalues
	expect_l(0) = complex_value_type( 0.022873738558717,0);
	expect_l(1) = complex_value_type(-0.011411124209185, 0.017516871011497);
	expect_l(2) = complex_value_type(-0.011411124209185,-0.017516871011497);


	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqi_controller<value_type> dlqi(Q, R, N);
	dlqi.solve(A, B, C, D, ts);
	DCS_DEBUG_TRACE("Gain = " << dlqi.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqi.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqi.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqi.gain(), nu, nz, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqi.are_solution(), nz, nz, tol );
	vector_type l(dlqi.eigenvalues());
	::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
	::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, nz, tol );
}


/*
DCS_TEST_DEF( free_func_ill_cond )
{
	DCS_TEST_CASE("free function - ill-conditioned problem");

	typedef double real_type;
	typedef real_type value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<value_type> vector_type;

	const std::size_t nx(6);
	const std::size_t nu(6);
	const std::size_t ny(1);
	const real_type ts(896);

	matrix_type A(nx,nx);
	A(0,0) =  0;           A(0,1) =  0;           A(0,2) =  0;           A(0,3) =  1;        A(0,4) =  0;        A(0,5) =  0;
	A(1,0) =  0;           A(1,1) =  0;           A(1,2) =  0;           A(1,3) =  0;        A(1,4) =  1;        A(1,5) =  0;
	A(2,0) =  0;           A(2,1) =  0;           A(2,2) =  0;           A(2,3) =  0;        A(2,4) =  0;        A(2,5) =  1;
	A(3,0) = -2.22045e-16; A(3,1) = -0;           A(3,2) = -0;           A(3,3) = -0.479283; A(3,4) = -0;        A(3,5) = -0;
	A(4,0) = -0;           A(4,1) = -2.22045e-16; A(4,2) = -0;           A(4,3) = -0;        A(4,4) = -0.274798; A(4,5) = -0;
	A(5,0) = -0;           A(5,1) = -0;           A(5,2) = -2.22045e-16; A(5,3) = -0;        A(5,4) = -0;        A(5,5) =  0.993379;

	matrix_type B(nx,nu);
	B(0,0) = 0;           B(0,1) = 0;           B(0,2) = 0;           B(0,3) = 0;           B(0,4) = 0;           B(0,5) = 0;
	B(1,0) = 0;           B(1,1) = 0;           B(1,2) = 0;           B(1,3) = 0;           B(1,4) = 0;           B(1,5) = 0;
	B(2,0) = 0;           B(2,1) = 0;           B(2,2) = 0;           B(2,3) = 0;           B(2,4) = 0;           B(2,5) = 0;
	B(3,0) = 2.22045e-16; B(3,1) = 2.22045e-16; B(3,2) = 2.22045e-16; B(3,3) = 2.22045e-16; B(3,4) = 2.22045e-16; B(3,5) = 2.22045e-16;
	B(4,0) = 2.22045e-16; B(4,1) = 2.22045e-16; B(4,2) = 2.22045e-16; B(4,3) = 2.22045e-16; B(4,4) = 2.22045e-16; B(4,5) = 2.22045e-16;
	B(5,0) = 2.22045e-16; B(5,1) = 2.22045e-16; B(5,2) = 2.22045e-16; B(5,3) = 2.22045e-16; B(5,4) = 2.22045e-16; B(5,5) = 2.22045e-16;

	matrix_type C(ny,nx);
	C(0,0) = 0; C(0,1) = 0; C(0,2) = 0; C(0,3) = 1; C(0,4) = 1; C(0,5) = 1;

	matrix_type D(ny,nu);
	D(0,0) = 0; D(0,1) = 0; D(0,2) = 0; D(0,3) = 0; D(0,4) = 0; D(0,5) = 0;

	matrix_type Q = ublas::identity_matrix<value_type>(nx+1,nx+1);

	matrix_type R = ublas::identity_matrix<value_type>(nu,nu);

	matrix_type N = ublas::zero_matrix<value_type>(nx+1,nu);


	matrix_type expect_K(nu,nx+1); // state-feedback optimal gain
	expect_K(0,0) = - 1.0247564929914510514e-17;
	expect_K(0,1) = - 1.1891333835485771088e-17;
	expect_K(0,2) = - 2.2895406420567840299e-15;
	expect_K(0,3) =   0.04615084747885873917;
	expect_K(0,4) =   0.05355371134420310042;
	expect_K(0,5) =  10.311156036192993923;
	expect_K(0,6) = - 0.068270164114244455411;
	expect_K(1,0) =   1.0247564929914510514e-17;
	expect_K(1,1) = - 1.1891333835485771088e-17;
	expect_K(1,2) = - 2.2895406420567840299e-15;
	expect_K(1,3) =   0.04615084747885873917;
	expect_K(1,4) =   0.05355371134420310042;
	expect_K(1,5) =  10.311156036192993923;
	expect_K(1,6) = - 0.068270164114244455411;
	expect_K(2,0) =   1.0247564929914507432e-17;
	expect_K(2,1) = - 1.1891333835485771088e-17;
	expect_K(2,2) = - 2.289540642056783241e-15;
	expect_K(2,3) =   0.046150847478858725292;
	expect_K(2,4) =   0.05355371134420310042;
	expect_K(2,5) =  10.311156036192993923;
	expect_K(2,6) = - 0.068270164114244455411;
	expect_K(3,0) =   1.0247564929914507432e-17;
	expect_K(3,1) = - 1.1891333835485771088e-17;
	expect_K(3,2) = - 2.289540642056783241e-15;
	expect_K(3,3) =   0.046150847478858725292;
	expect_K(3,4) =   0.05355371134420310042;
	expect_K(3,5) =  10.311156036192993923;
	expect_K(3,6) = - 0.068270164114244455411;
	expect_K(4,0) =   1.0247564929914507432e-17;
	expect_K(4,1) = - 1.1891333835485772629e-17;
	expect_K(4,2) = - 2.289540642056783241e-15;
	expect_K(4,3) =   0.046150847478858725292;
	expect_K(4,4) =   0.05355371134420310042;
	expect_K(4,5) =  10.311156036192993923;
	expect_K(4,6) = - 0.068270164114244455411;
	expect_K(5,0) =   1.0247564929914507432e-17;
	expect_K(5,1) = - 1.1891333835485772629e-17;
	expect_K(5,2) = - 2.289540642056783241e-15;
	expect_K(5,3) =   0.046150847478858725292;
	expect_K(5,4) =   0.05355371134420310042;
	expect_K(5,5) =  10.311156036192993923;
	expect_K(5,6) = - 0.068270164114244455411;

	matrix_type K = dcs_ctrl::dlqi_solve(A, B, C, D, ts, Q, R, N);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	DCS_DEBUG_TRACE("Gain = " << K);
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, K, nu, nx+1, tol );
}
*/


int main()
{
	// Tests labeled with 'mathematica' keyword are taken from the Wolfram
	// Mathematica 8 documentation
	// (http://reference.wolfram.com/mathematica/ref/LQRegulatorGains.html)
	//
	// All tests has been validated with MATLAB 2009b
 
	DCS_TEST_SUITE("DCS Control :: DLQI");

	DCS_TEST_BEGIN();

	DCS_TEST_DO( matlab_1 );
//	DCS_TEST_DO( free_func_ill_cond );

	DCS_TEST_END();
}
