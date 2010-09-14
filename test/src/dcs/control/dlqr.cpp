#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cstddef>
#include <dcs/control/design/dlqr.hpp>
#include <dcs/test.hpp>
#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace dcs_ctrl = dcs::control;


const double tol = 1.0e-5;


DCS_TEST_DEF( test_free_func )
{
	DCS_DEBUG_TRACE("Test Case: test_free_func");

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

	matrix_type C(m,n);
	C(0,0) = 2.011946497e-6; C(0,1) = 0;

	matrix_type D(m,m);
	D(0,0) = 0;

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

	matrix_type K = dcs_ctrl::dlqr_solve(A, B, Q, R, N);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
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

	matrix_type C(m,n);
	C(0,0) = 2.011946497e-6; C(0,1) = 0;

	matrix_type D(m,m);
	D(0,0) = 0;

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

	dcs_ctrl::dlqr_controller<value_type> dlqr(Q, R, N);
	dlqr.solve(A, B);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("Q = " << Q);
	DCS_DEBUG_TRACE("R = " << R);
	DCS_DEBUG_TRACE("N = " << N);
	DCS_DEBUG_TRACE("Gain = " << dlqr.gain());
	DCS_DEBUG_TRACE("Riccati's Solution = " << dlqr.are_solution());
	DCS_DEBUG_TRACE("Closed-loop Eigenvalues = " << dlqr.eigenvalues());
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, dlqr.gain(), m, n, tol );
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_S, dlqr.are_solution(), n, n, tol );
	DCS_TEST_CHECK_VECTOR_CLOSE( expect_e, dlqr.eigenvalues(), n, tol );
}


int main()
{
	DCS_TEST_BEGIN();

	DCS_TEST_DO( test_free_func );
	DCS_TEST_DO( test_oo );

	DCS_TEST_END();
}
