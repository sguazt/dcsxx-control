#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cstddef>
#include <dcs/control/dlqr.hpp>
#include <dcs/test.hpp>
#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace dcs_ctrl = dcs::control;

DCS_TEST_DEF( test1 )
{
	DCS_DEBUG_TRACE("Test Case: test1");

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


	matrix_type expect_K(n,m); // state-feedback optimal gain
	expect_K(0,0) = -1.0e-3*0.792840452340859; expect_K(0,1) = 1.0e-3*0.802772572963461;

	matrix_type expect_S(n,n); // solution of the associated Riccati equation
	expect_S(0,0) = 1.0e-5*0.634804227345193; expect_S(0,1) =-1.0e-5*0.642756176807922;
	expect_S(1,0) =-1.0e-5*0.642756176807922; expect_S(1,1) = 1.0e-5*0.650824595759136;

	vector_type expect_E(n); // closed-loop eigenvalues
	expect_E(0) = 0.996904570704525;
	expect_E(1) = 0.990284265674475;

	dcs_ctrl::dlqr(A, B, Q, R, N);
}


int main()
{
	DCS_TEST_BEGIN();

	DCS_TEST_DO( test1 );

	DCS_TEST_END();
}
