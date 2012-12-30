/**
 * \file test/src/dcs/control/dlqr.cpp
 *
 * \brief Test suite for DLQR controllers.
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

#include <algorithm>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <cstddef>
#include <dcs/control/design/dlqr.hpp>
#include <dcs/test.hpp>
#include <iostream>
#include "./utility.hpp"


namespace ublas = boost::numeric::ublas;
namespace dcs_ctrl = dcs::control;


const double tol = 1.0e-5;


DCS_TEST_DEF( free_func )
{
	DCS_TEST_CASE("test_free_func");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<value_type> vector_type;

	const std::size_t n = 2;
	const std::size_t m = 1;

	matrix_type A(n,n);
	A(0,0) =  0;                    A(0,1) = 1;
	A(1,0) = -0.988011751199999956; A(1,1) = 1.98799160895100258;

	matrix_type B(n,m);
	B(0,0) = 0;
	B(1,0) = 1;

	matrix_type Q(n,n);
	Q(0,0) = 4.04792870679057178e-12; Q(0,1) = 0;
	Q(1,0) = 0;                       Q(1,1) = 0;

	matrix_type R(m,m);
	R(0,0) = 0.00810385427000000040;

	matrix_type N(n,m);
	N(0,0) = 0;
	N(1,0) = 0;


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -1.0e-3*0.792840452340859; expect_K(0,1) = 1.0e-3*0.802772572963461;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 1.0e-5*0.634804227345193; expect_S(0,1) =-1.0e-5*0.642756176807922;
	expect_S(1,0) =-1.0e-5*0.642756176807922; expect_S(1,1) = 1.0e-5*0.650824595759136;

	vector_type expect_e(n); // closed-loop eigenvalues
	expect_e(0) = 0.996904570704525;
	expect_e(1) = 0.990284265674475;

	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	matrix_type K = dcs_ctrl::dlqr_solve(A, B, Q, R, N);
	DCS_DEBUG_TRACE("Gain = " << K);
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, K, m, n, tol );
}


DCS_TEST_DEF( oo )
{
	DCS_TEST_CASE("test_oo");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n = 2;
	const std::size_t m = 1;

	matrix_type A(n,n);
	A(0,0) =  0;                    A(0,1) = 1;
	A(1,0) = -0.988011751199999956; A(1,1) = 1.98799160895100258;

	matrix_type B(n,m);
	B(0,0) = 0;
	B(1,0) = 1;

	matrix_type Q(n,n);
	Q(0,0) = 4.04792870679057178e-12; Q(0,1) = 0;
	Q(1,0) = 0;                       Q(1,1) = 0;

	matrix_type R(m,m);
	R(0,0) = 0.00810385427000000040;

	matrix_type N(n,m);
	N(0,0) = 0;
	N(1,0) = 0;


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -1.0e-3*0.792840452340859; expect_K(0,1) = 1.0e-3*0.802772572963461;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 1.0e-5*0.634804227345193; expect_S(0,1) =-1.0e-5*0.642756176807922;
	expect_S(1,0) =-1.0e-5*0.642756176807922; expect_S(1,1) = 1.0e-5*0.650824595759136;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_value_type(0.996904570704525,0);
	expect_l(1) = complex_value_type(0.990284265674475,0);

	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqr_controller<value_type> dlqr(Q, R, N);
	dlqr.solve(A, B);
	DCS_DEBUG_TRACE("Gain = " << dlqr.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqr.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqr.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqr.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqr.are_solution(), n, n, tol );
	vector_type l(dlqr.eigenvalues());
    ::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
    ::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );

	matrix_type X(16,n);

	X( 0,0) = -0.08; X( 0,1) =  0.04;
	X( 1,0) = -0.06; X( 1,1) =  0.02;
	X( 2,0) = -0.05; X( 2,1) = -0.02;
	X( 3,0) =  0.04; X( 3,1) = -0.01;
	X( 4,0) = -0.38; X( 4,1) =  0.09;
	X( 5,0) =  0.48; X( 5,1) =  0.12;
	X( 6,0) =  0.34; X( 6,1) =  0.04;
	X( 7,0) =  0.76; X( 7,1) =  0.50;
	X( 8,0) = -0.07; X( 8,1) =  0.33;
	X( 9,0) =  0.18; X( 9,1) =  0.20;
	X(10,0) = -0.07; X(10,1) =  0.02;
	X(11,0) = -0.31; X(11,1) = -0.04;
	X(12,0) =  0.06; X(12,1) = -0.09;
	X(13,0) =  0.69; X(13,1) = -1.00;
	X(14,0) = -0.32; X(14,1) = -0.02;
	X(15,0) =  0.36; X(15,1) =  0.40;

	matrix_type expect_U(m,16);

  expect_U(0,0) = -0.000095538139106; expect_U(0,1) = -0.000063625878600; expect_U(0,2) = -0.000023586571158; expect_U(0,3) =  0.000039741343823; expect_U(0,4) = -0.000373528903456; expect_U(0,5) =  0.000284230708368; expect_U(0,6) =  0.000237454850877; expect_U(0,7) =  0.000201172457297; expect_U(0,8) = -0.000320413780742; expect_U(0,9) = -0.000017843233171; expect_U(0,10) = -0.000071554283123; expect_U(0,11) = -0.000213669637307; expect_U(0,12) =  0.000119819958707; expect_U(0,13) =  0.001349832485079; expect_U(0,14) = -0.000237653493290; expect_U(0,15) = -0.000035686466343;

	matrix_type U;
	U = dlqr.control(X);
	DCS_DEBUG_TRACE("Input X = " << X);
	DCS_DEBUG_TRACE("Control U = " << U);
	DCS_DEBUG_TRACE("Control expect U = " << expect_U);
	DCS_TEST_CHECK_MATRIX_CLOSE(U, expect_U, m, 16, tol);
}


DCS_TEST_DEF( mathematica_1 )
{
	// This is the first example for discrete-time systems found in the
	// Mathematica 8 doc.
	// See:
	//   http://reference.wolfram.com/mathematica/ref/LQRegulatorGains.html

	DCS_TEST_CASE("Mathematica #1");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n = 4;
	const std::size_t m = 1;

	matrix_type A(n,n);
	A(0,0) = 0.4725; A(0,1) = 0.2376; A(0,2) = 0.0589; A(0,3) = 0.1971;
	A(1,0) = 0.1451; A(1,1) = 0.5669; A(1,2) = 0.2311; A(1,3) = 0.0439;
	A(2,0) = 0.0932; A(2,1) = 0.1190; A(2,2) = 0.5752; A(2,3) = 0.2319;
	A(3,0) = 0.2628; A(3,1) = 0.0757; A(3,2) = 0.1406; A(3,3) = 0.4465;

	matrix_type B(n,m);
	B(0,0) =  0.5711;
	B(1,0) = -0.3999;
	B(2,0) =  0.6899;
	B(3,0) =  0.8156;

	matrix_type Q(n,n);
	Q(0,0) = 0.5; Q(0,1) = 0.0; Q(0,2) = 0; Q(0,3) = 0;
	Q(1,0) = 0.0; Q(1,1) = 0.8; Q(1,2) = 0; Q(1,3) = 0;
	Q(2,0) = 0.0; Q(2,1) = 0.0; Q(2,2) = 2; Q(2,3) = 0;
	Q(3,0) = 0.0; Q(3,1) = 0.0; Q(3,2) = 0; Q(3,3) = 4;

	matrix_type R(m,m);
	R(0,0) = 1;

	matrix_type N(ublas::zero_matrix<value_type>(n,m));


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = 0.258263408785324; expect_K(0,1) = 0.130076037449530; expect_K(0,2) = 0.301712976754135; expect_K(0,3) = 0.393239174520131;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 0.820688030019406; expect_S(0,1) = 0.366157393574175; expect_S(0,2) = 0.142357177991738; expect_S(0,3) = 0.230214495860850;
	expect_S(1,0) = 0.366157393574175; expect_S(1,1) = 1.502565253736899; expect_S(1,2) = 0.491675367390775; expect_S(1,3) = 0.233049076800882;
	expect_S(2,0) = 0.142357177991738; expect_S(2,1) = 0.491675367390775; expect_S(2,2) = 2.770266171822138; expect_S(2,3) = 0.158115947330550;
	expect_S(3,0) = 0.230214495860850; expect_S(3,1) = 0.233049076800882; expect_S(3,2) = 0.158115947330550; expect_S(3,3) = 4.285246483933100;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_value_type(0.715676446034849, 0.000000000000000);
	expect_l(1) = complex_value_type(0.460342919154304, 0.000000000000000);
	expect_l(2) = complex_value_type(0.130363078014160, 0.048280928782719);
	expect_l(3) = complex_value_type(0.130363078014160,-0.048280928782719);

	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqr_controller<value_type> dlqr(Q, R, N);
	dlqr.solve(A, B);
	DCS_DEBUG_TRACE("Gain = " << dlqr.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqr.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqr.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqr.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqr.are_solution(), n, n, tol );
	vector_type l(dlqr.eigenvalues());
    ::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
    ::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );
}


DCS_TEST_DEF( mathematica_2 )
{
	// This is the second example for discrete-time systems found in the
	// Mathematica 8 doc.
	// See:
	//   http://reference.wolfram.com/mathematica/ref/LQRegulatorGains.html

	DCS_TEST_CASE("Mathematica #2");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n = 3;
	const std::size_t m = 1;

	matrix_type A(n,n);
	A(0,0) = 1; A(0,1) = 1; A(0,2) = -2;
	A(1,0) = 0; A(1,1) = 1; A(1,2) =  1;
	A(2,0) = 0; A(2,1) = 0; A(2,2) =  1;

	matrix_type B(n,m);
	B(0,0) = 1;
	B(1,0) = 0;
	B(2,0) = 1;

	matrix_type Q(ublas::identity_matrix<value_type>(n,n));

	matrix_type R(m,m);
	R(0,0) = 10;

	matrix_type N(ublas::zero_matrix<value_type>(n,m));


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = 0.131575602686720; expect_K(0,1) = 1.015926918639910; expect_K(0,2) = 1.316503115669289;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) =  1.0e+2*0.077212408523700; expect_S(0,1) = 1.0e+2*0.254481597239685; expect_S(0,2) = -1.0e+2*0.001210476586442;
	expect_S(1,0) =  1.0e+2*0.254481597239685; expect_S(1,1) = 1.0e+2*1.182679389650807; expect_S(1,2) =  1.0e+2*0.256345692556058;
	expect_S(2,0) = -1.0e+2*0.001210476586442; expect_S(2,1) = 1.0e+2*0.256345692556058; expect_S(2,2) =  1.0e+2*0.402837910468753;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_value_type(0.464520830811466, 0.249314173526600);
	expect_l(1) = complex_value_type(0.464520830811466,-0.249314173526600);
	expect_l(2) = complex_value_type(0.622879620021054, 0.000000000000000);


	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqr_controller<value_type> dlqr(Q, R, N);
	dlqr.solve(A, B);
	DCS_DEBUG_TRACE("Gain = " << dlqr.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqr.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqr.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqr.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqr.are_solution(), n, n, tol );
	vector_type l(dlqr.eigenvalues());
    ::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
    ::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );
}


DCS_TEST_DEF( mathematica_3 )
{
	// This is the third example for discrete-time systems found in the
	// Mathematica 8 doc.
	// See:
	//   http://reference.wolfram.com/mathematica/ref/LQRegulatorGains.html

	DCS_TEST_CASE("Mathematica #3");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n = 2;
	const std::size_t m = 2;

	matrix_type A(n,n);
	A(0,0) = 0.951230000000000; A(0,1) = 0.000000000000000;
	A(1,0) = 0.000000000000000; A(1,1) = 0.904838000000000;

	matrix_type B(n,m);
	B(0,0) =  4.900000000000000; B(0,1) = 4.900000000000000;
	B(1,0) = -0.019000000000000; B(1,1) = 0.009000000000000;


	matrix_type Q(n,n);
	Q(0,0) = 1.0e+2*0.000100000000000; Q(0,1) = 1.0e+2*0.000000000000000;
	Q(1,0) = 1.0e+2*0.000000000000000; Q(1,1) = 1.0e+2*1.000000000000000;


	matrix_type R(m,m);
	R(0,0) = 2.000000000000000; R(0,1) = 0.000000000000000;
	R(1,0) = 0.000000000000000; R(1,1) = 0.500000000000000;

	matrix_type N(ublas::zero_matrix<value_type>(n,m));


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = 0.022043168004973; expect_K(0,1) = -3.296571964433571;
	expect_K(1,0) = 0.078004554558866; expect_K(1,1) =  3.751723555757965;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) =  1.0e+2*0.000178886938538; expect_S(0,1) = -1.0e+2*0.001642946262000;
	expect_S(1,0) = -1.0e+2*0.001642946262000; expect_S(1,1) =  1.0e+2*3.736813637743126;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_value_type(0.459187566597895,0);
	expect_l(1) = complex_value_type(0.810246213513234,0);


	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqr_controller<value_type> dlqr(Q, R, N);
	dlqr.solve(A, B);
	DCS_DEBUG_TRACE("Gain = " << dlqr.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqr.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqr.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqr.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqr.are_solution(), n, n, tol );
	vector_type l(dlqr.eigenvalues());
    ::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
    ::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );
}


DCS_TEST_DEF( mathematica_4 )
{
	// This is the third example for discrete-time systems found in the
	// Mathematica 8 doc.
	// See:
	//   http://reference.wolfram.com/mathematica/ref/LQRegulatorGains.html

	DCS_TEST_CASE("Mathematica #4");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n = 3;
	const std::size_t m = 2;

	matrix_type A(n,n);
	A(0,0) = -0.5000; A(0,1) = 0.0000; A(0,2) = 0.0000;
	A(1,0) =  0.0000; A(1,1) = 0.9999; A(1,2) = 0.0000;
	A(2,0) =  0.0000; A(2,1) = 0.0000; A(2,2) = 0.8187;

	matrix_type B(n,m);
	B(0,0) =  0.1850; B(0,1) = 0.1974;
	B(1,0) =  0.1000; B(1,1) = 0.1390;
	B(2,0) =  0.1813; B(2,1) = 0.0000;


	matrix_type Q(n,n);
	Q(0,0) = 10; Q(0,1) = 0; Q(0,2) = 0;
	Q(1,0) =  0; Q(1,1) = 1; Q(1,2) = 0;
	Q(2,0) =  0; Q(2,1) = 0; Q(2,2) = 1;

	matrix_type R(m,m);
	R(0,0) = 0.4; R(0,1) = 0.1;
	R(1,0) = 0.1; R(1,1) = 0.1;

	matrix_type N(n,m);
	N(0,0) = 1; N(0,1) = 1;
	N(1,0) = 2; N(1,1) = 1;
	N(2,0) = 2; N(2,1) = 1;


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -0.10259451927409935; expect_K(0,1) = 1.0550552116729919;  expect_K(0,2) =  3.3699401712231705;
	expect_K(1,0) = -0.4714396669727725;  expect_K(1,1) = 0.4206832287681541; expect_K(1,2) = -1.796589906583708;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 13.104637764024423;   expect_S(0,1) =  0.2570942162727292; expect_S(0,2) =  0.18034745994439821;
	expect_S(1,0) =  0.25709421627274576; expect_S(1,1) = -2.741491441343835;  expect_S(1,2) = -6.024514042762959;
	expect_S(2,0) =  0.18034745994441229; expect_S(2,1) = -6.024514042762947;  expect_S(2,2) = -4.280531570924603;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_value_type( 0.850245380354058,0);
	expect_l(1) = complex_value_type(-0.35398128255121586,0);
	expect_l(2) = complex_value_type( 0.15942743551445882,0);


	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqr_controller<value_type> dlqr(Q, R, N);
	dlqr.solve(A, B);
	DCS_DEBUG_TRACE("Gain = " << dlqr.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqr.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqr.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqr.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqr.are_solution(), n, n, tol );
	vector_type l(dlqr.eigenvalues());
    ::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
    ::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );
}


int main()
{
	// Tests labeled with 'mathematica' keyword are taken from the Wolfram
	// Mathematica 8 documentation
	// (http://reference.wolfram.com/mathematica/ref/LQRegulatorGains.html)
	//
	// All tests has been validated with MATLAB 2009b
 
	DCS_TEST_SUITE("DCS Control :: DLQR");

	DCS_TEST_BEGIN();

	DCS_TEST_DO( free_func );
	DCS_TEST_DO( oo );
	DCS_TEST_DO( mathematica_1 );
	DCS_TEST_DO( mathematica_2 );
	DCS_TEST_DO( mathematica_3 );
	DCS_TEST_DO( mathematica_4 );

	DCS_TEST_END();
}
