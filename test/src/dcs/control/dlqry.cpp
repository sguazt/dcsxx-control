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


DCS_TEST_DEF( test_free_func )
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

	matrix_type N(p,m);
	N(0,0) = 1;
	N(1,0) = 0;


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -0.715925028553135; expect_K(0,1) = 1.381105915011857;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) =  1.208699380958362; expect_S(0,1) = -0.976250563077556;
	expect_S(1,0) = -0.976250563077556; expect_S(1,1) =  3.071182198140876;

	vector_type expect_e(n); // closed-loop eigenvalues
	expect_e(0) = complex_type(0.428442042494072, 0.431077736177637);
	expect_e(1) = complex_type(0.428442042494072,-0.431077736177637);

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


DCS_TEST_DEF( test_oo )
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

	matrix_type N(p,m);
	N(0,0) = 1;
	N(1,0) = 0;


	matrix_type expect_K(m,n); // state-feedback optimal gain
	expect_K(0,0) = -0.715925028553135; expect_K(0,1) = 1.381105915011857;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) =  1.208699380958362; expect_S(0,1) = -0.976250563077556;
	expect_S(1,0) = -0.976250563077556; expect_S(1,1) =  3.071182198140876;

	vector_type expect_e(n); // closed-loop eigenvalues
	expect_e(0) = complex_type(0.428442042494072, 0.431077736177637);
	expect_e(1) = complex_type(0.428442042494072,-0.431077736177637);

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
	DCS_TEST_BEGIN();

	DCS_TEST_DO( test_free_func );
	DCS_TEST_DO( test_oo );

	DCS_TEST_END();
}
