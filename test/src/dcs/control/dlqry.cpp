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

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 0.062747695129094; expect_S(0,1) = 0.006690561505124;
	expect_S(1,0) = 0.006690561505124; expect_S(1,1) = 0.072218977786870;

	vector_type expect_e(n); // closed-loop eigenvalues
	expect_e(0) = complex_type(0.129095068582370, 0.824459776759475);
	expect_e(1) = complex_type(0.129095068582370,-0.824459776759475);

	matrix_type K = dcs_ctrl::dlqry_solve(A, B, C, D, Q, R, N);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
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

	vector_type expect_e(n); // closed-loop eigenvalues
	expect_e(0) = complex_type(0.129095068582370, 0.824459776759475);
	expect_e(1) = complex_type(0.129095068582370,-0.824459776759475);


	dcs_ctrl::dlqry_controller<value_type> dlqry(Q, R, N);
	dlqry.solve(A, B, C, D);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	DCS_DEBUG_TRACE("Gain = " << dlqry.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqry.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqry.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqry.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqry.are_solution(), n, n, tol );
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_e, dlqry.eigenvalues(), n, tol );
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

	matrix_type Q(n,n);
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

	vector_type expect_e(n); // closed-loop eigenvalues
	expect_e(0) = complex_type(-0.036208017520516, 0.050803722152156);
	expect_e(1) = complex_type(-0.036208017520516,-0.050803722152156);

	dcs_ctrl::dlqry_controller<value_type> dlqry(Q, R, N);
	dlqry.solve(A, B, C, D);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	DCS_DEBUG_TRACE("Gain = " << dlqry.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqry.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqry.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqry.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqry.are_solution(), n, n, tol );
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_e, dlqry.eigenvalues(), n, tol );
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

	DCS_TEST_END();
}
