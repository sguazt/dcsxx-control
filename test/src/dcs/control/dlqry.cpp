/**
 * \file test/src/dcs/control/dlqry.cpp
 *
 * \brief Test suite for the DLQRY controller.
 *
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <cstddef>
#include <dcs/control/design/dlqry.hpp>
#include <dcs/test.hpp>
#include <iostream>
#include "./utility.hpp"


namespace ublas = boost::numeric::ublas;
namespace dcs_ctrl = dcs::control;


const double tol = 1.0e-5;


DCS_TEST_DEF( free_func )
{
	DCS_DEBUG_TRACE("Test Case: test_free_func");

	typedef double value_type;
	typedef std::complex<double> complex_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_type> vector_type;

	const std::size_t n = 2; // number of states
	const std::size_t m = 1; // number of inputs
	const std::size_t p = 2; // number of outputs

	matrix_type A(n,n);
	A(0,0) =  0.250000000000000000; A(0,1) = 0.80000000000000000;
	A(1,0) = -0.988011751199999956; A(1,1) = 1.98799160895100258;

	matrix_type B(n,m);
	B(0,0) = 0;
	B(1,0) = 1;

	matrix_type C(p,n);
	C(0,0) = 2.011946497e-6; C(0,1) = 0;
	C(1,0) = 0.000000000000; C(1,1) = 1;

	matrix_type D(p,m);
	D(0,0) = 0.5;
	D(1,0) = 0.5;

	matrix_type Q(n,n);
	Q(0,0) = 4.04792870679057178e-12; Q(0,1) = 0;
	Q(1,0) = 0;                       Q(1,1) = 2;

	matrix_type R(m,m);
	R(0,0) = 0.00810385427000000040;

	matrix_type N(ublas::zero_matrix<value_type>(p,m));


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -0.120072092580768; expect_K(0,1) = 1.979799862835259;

	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	matrix_type K = dcs_ctrl::dlqry_solve(A, B, C, D, Q, R, N);
	DCS_DEBUG_TRACE("Gain = " << K);
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, K, m, n, tol );
}


DCS_TEST_DEF( oo )
{
	DCS_DEBUG_TRACE("Test Case: test_oo");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_type> vector_type;

	const std::size_t n = 2; // number of states
	const std::size_t m = 1; // number of inputs
	const std::size_t p = 2; // number of outputs

	matrix_type A(n,n);
	A(0,0) =  0.250000000000000000; A(0,1) = 0.80000000000000000;
	A(1,0) = -0.988011751199999956; A(1,1) = 1.98799160895100258;

	matrix_type B(n,m);
	B(0,0) = 0;
	B(1,0) = 1;

	matrix_type C(p,n);
	C(0,0) = 2.011946497e-6; C(0,1) = 0;
	C(1,0) = 0.000000000000; C(1,1) = 1;

	matrix_type D(p,m);
	D(0,0) = 0.5;
	D(1,0) = 0.5;

	matrix_type Q(n,n);
	Q(0,0) = 4.04792870679057178e-12; Q(0,1) = 0;
	Q(1,0) = 0;                       Q(1,1) = 2;

	matrix_type R(m,m);
	R(0,0) = 0.00810385427000000040;

	matrix_type N(ublas::zero_matrix<value_type>(p,m));


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -0.120072092580768; expect_K(0,1) = 1.979799862835259;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 0.062747695129094; expect_S(0,1) = 0.006690561505124;
	expect_S(1,0) = 0.006690561505124; expect_S(1,1) = 0.072218977786870;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_type(0.129095068582370, 0.824459776759475);
	expect_l(1) = complex_type(0.129095068582370,-0.824459776759475);


	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqry_controller<value_type> dlqry(Q, R, N);
	dlqry.solve(A, B, C, D);
	DCS_DEBUG_TRACE("Gain = " << dlqry.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqry.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqry.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqry.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqry.are_solution(), n, n, tol );
	vector_type l(dlqry.eigenvalues());
	::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
	::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );
}


DCS_TEST_DEF( mathematica_1 )
{
	DCS_DEBUG_TRACE("Test Case: Mathematica #1");

	// This is the third example for discrete-time systems found in the
	// Mathematica 8 doc.
	// See:
	//   http://reference.wolfram.com/mathematica/ref/LQOutputRegulatorGains.html
 
	typedef double value_type;
	typedef ::std::complex<value_type> complex_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_type> vector_type;

	const std::size_t n = 2; // number of states
	const std::size_t m = 1; // number of inputs
	const std::size_t p = 2; // number of outputs

	matrix_type A(n,n);
	A(0,0) = -0.05156; A(0,1) = -0.058770;
	A(1,0) =  0.01503; A(1,1) = -0.005887;

	matrix_type B(n,m);
	B(0,0) =  0.00000;
	B(1,0) = -0.03384;

	matrix_type C(p,n);
	C(0,0) =  1.1; C(0,1) = -0.5;
	C(1,0) = -0.7; C(1,1) =  1.0;

	matrix_type D(p,m);
	D(0,0) = 0.1;
	D(1,0) = 0.2;

	matrix_type Q(p,p);
	Q(0,0) = 10; Q(0,1) = 0.0;
	Q(1,0) =  0; Q(1,1) = 0.1;

	matrix_type R(m,m);
	R(0,0) = 1;

	matrix_type N(ublas::zero_matrix<value_type>(p,m));


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = 0.972149539638538; expect_K(0,1) = -0.442347371188903;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) =  11.141022805252103; expect_S(0,1) = -5.057580732711166;
	expect_S(1,0) = - 5.057580732711166; expect_S(1,1) =  2.418501386750889;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_type(-0.036208017520516, 0.050803722152156);
	expect_l(1) = complex_type(-0.036208017520516,-0.050803722152156);

	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqry_controller<value_type> dlqry(Q, R, N);
	dlqry.solve(A, B, C, D);
	DCS_DEBUG_TRACE("Gain = " << dlqry.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqry.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqry.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqry.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqry.are_solution(), n, n, tol );
	vector_type l(dlqry.eigenvalues());
	::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
	::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );
}


DCS_TEST_DEF( mathematica_2 )
{
	DCS_DEBUG_TRACE("Test Case: Mathematica #2");

	// This is the third example for discrete-time systems found in the
	// Mathematica 8 doc.
	// See:
	//   http://reference.wolfram.com/mathematica/ref/LQOutputRegulatorGains.html
 
	typedef double value_type;
	typedef ::std::complex<value_type> complex_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_type> vector_type;

	const std::size_t n = 4; // number of states
	const std::size_t m = 2; // number of inputs
	const std::size_t p = 2; // number of outputs

	matrix_type A(n,n);
	A(0,0) =  0.00; A(0,1) =  0.00; A(0,2) =  1; A(0,3) =  0;
	A(1,0) =  0.00; A(1,1) =  0.00; A(1,2) =  0; A(1,3) =  1;
	A(2,0) = -0.24; A(2,1) =  0.00; A(2,2) = -1; A(2,3) =  0;
	A(3,0) =  0.00; A(3,1) = -0.24; A(3,2) =  0; A(3,3) = -1;

	matrix_type B(n,m);
	B(0,0) = 0; B(0,1) = 0;
	B(1,0) = 0; B(1,1) = 0;
	B(2,0) = 1; B(2,1) = 0;
	B(3,0) = 0; B(3,1) = 1;

	matrix_type C(p,n);
	C(0,0) = -0.3; C(0,1) = -0.24; C(0,2) =  1; C(0,3) = -0.4;
	C(1,0) =  6.0; C(1,1) = -0.20; C(1,2) = 10; C(1,3) = -0.5;

	matrix_type D(p,m);
	D(0,0) = 0; D(0,1) = 1;
	D(1,0) = 0; D(1,1) = 1;

	matrix_type Q(ublas::identity_matrix<value_type>(p,p));

	matrix_type R(ublas::identity_matrix<value_type>(m,m));

	matrix_type N(ublas::identity_matrix<value_type>(p,m));


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -0.5727425667891235; expect_K(0,1) = -0.015333480978689376; expect_K(0,2) = -0.6509281335164934; expect_K(0,3) = -0.005825529959808867;
	expect_K(1,0) =  2.4734966016228954; expect_K(1,1) = -0.12406716180484925;  expect_K(1,2) =  4.353453877275993;  expect_K(1,3) = -0.27783321717671866;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 6.450086022615284;  expect_S(0,1) = 0.29484190309117786;  expect_S(0,2) = 7.92070501731365;    expect_S(0,3) = 0.3409790812924333;
	expect_S(1,0) = 0.2948419030911792; expect_S(1,1) = 0.013477610610715578; expect_S(1,2) = 0.36206582872544446; expect_S(1,3) = 0.015586601618962546;
	expect_S(2,0) = 7.920705017313646;  expect_S(2,1) = 0.3620658287254433;   expect_S(2,2) = 9.726625001794925;   expect_S(2,3) = 0.4187222791327824;
	expect_S(3,0) = 0.3409790812924328; expect_S(3,1) = 0.015586601618962477; expect_S(3,2) = 0.41872227913278426; expect_S(3,3) = 0.01802560980913717;

	vector_type expect_l(n); // closed-loop eigenvalues
	expect_l(0) = complex_type(-0.8240762158433446,0);
	expect_l(1) = complex_type(-0.5195595835434045,0);
	expect_l(2) = complex_type( 0.27784825088463877,0);
	expect_l(3) = complex_type(-0.005451100804677363,0);

	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqry_controller<value_type> dlqry(Q, R, N);
	dlqry.solve(A, B, C, D);
	DCS_DEBUG_TRACE("Gain = " << dlqry.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqry.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqry.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqry.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqry.are_solution(), n, n, tol );
	vector_type l(dlqry.eigenvalues());
	::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
	::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_l, l, n, tol );
}


DCS_TEST_DEF( eesim )
{
	DCS_DEBUG_TRACE("Test Case: Mathematica #2");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_type> vector_type;

	const std::size_t nx = 6; // number of states
	const std::size_t nu = 6; // number of inputs
	const std::size_t ny = 1; // number of outputs

	matrix_type A(nx,nx);
	A(0,0) =  0;                  A(0,1) =  0;                  A(0,2) =  0;                 A(0,3) =  1;                 A(0,4) =  0;                 A(0,5) =  0;
	A(1,0) =  0;                  A(1,1) =  0;                  A(1,2) =  0;                 A(1,3) =  0;                 A(1,4) =  1;                 A(1,5) =  0;
	A(2,0) =  0;                  A(2,1) =  0;                  A(2,2) =  0;                 A(2,3) =  0;                 A(2,4) =  0;                 A(2,5) =  1;
	A(3,0) = -0.0639186388233641; A(3,1) = -0;                  A(3,2) = -0;                 A(3,3) =  0.794311427335589; A(3,4) = -0;                 A(3,5) = -0;
	A(4,0) = -0;                  A(4,1) = -0.0635751341222181; A(4,2) = -0;                 A(4,3) = -0;                 A(4,4) =  0.794589468149401; A(4,5) = -0;
	A(5,0) = -0;                  A(5,1) = -0;                  A(5,2) = -0.063972851249695; A(5,3) = -0;                 A(5,4) = -0;                 A(5,5) =  0.794267523469565;

	matrix_type B(nx,nu);
	B(0,0) =  0;                  B(0,1) =  0;                  B(0,2) =  0;                  B(0,3) = 0;                 B(0,4) = 0;                 B(0,5) = 0;
	B(1,0) =  0;                  B(1,1) =  0;                  B(1,2) =  0;                  B(1,3) = 0;                 B(1,4) = 0;                 B(1,5) = 0;
	B(2,0) =  0;                  B(2,1) =  0;                  B(2,2) =  0;                  B(2,3) = 0;                 B(2,4) = 0;                 B(2,5) = 0;
	B(3,0) = -0.0187739831892682; B(3,1) = -0.0187739831892682; B(3,2) = -0.0187739831892682; B(3,3) = 0.233302674436667; B(3,4) = 0.233302674436667; B(3,5) = 0.233302674436667;
	B(4,0) = -0.018657202160102;  B(4,1) = -0.018657202160102;  B(4,2) = -0.018657202160102;  B(4,3) = 0.233185765885323; B(4,4) = 0.233185765885323; B(4,5) = 0.233185765885323;
	B(5,0) = -0.0187924305365045; B(5,1) = -0.0187924305365045; B(5,2) = -0.0187924305365045; B(5,3) = 0.233321119359592; B(5,4) = 0.233321119359592; B(5,5) = 0.233321119359592;

	matrix_type C(ny,nx);
	C(0,0) = 0; C(0,1) = 0; C(0,2) = 0; C(0,3) = 1; C(0,4) = 1; C(0,5) = 1;

	matrix_type D(ny,nu);
	D(0,0) = 0; D(0,1) = 0; D(0,2) = 0; D(0,3); D(0,4) = 0; D(0,5) = 0;

	matrix_type Q(ublas::identity_matrix<value_type>(ny,ny));

	matrix_type R(ublas::identity_matrix<value_type>(nu,nu));

	//matrix_type N(ublas::identity_matrix<value_type>(ny,nu));
	matrix_type N(ublas::zero_matrix<value_type>(ny,nu));


	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	dcs_ctrl::dlqry_controller<value_type> dlqry(Q, R, N);
	dlqry.solve(A, B, C, D);
	DCS_DEBUG_TRACE("Gain = " << dlqry.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqry.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqry.eigenvalues());
}


int main()
{
	// Tests labeled with 'mathematica' keyword are taken from the Wolfram
	// Mathematica 8 documentation
	// (http://reference.wolfram.com/mathematica/ref/LQOutputRegulatorGains.html)
	//
	// All tests has been validated with MATLAB 2009b
 
	DCS_TEST_BEGIN();

	DCS_TEST_DO( free_func );
	DCS_TEST_DO( oo );
	DCS_TEST_DO( mathematica_1 );
	DCS_TEST_DO( mathematica_2 );
	DCS_TEST_DO( eesim );

	DCS_TEST_END();
}
