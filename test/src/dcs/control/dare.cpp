#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <complex>
#include <dcs/control/solver/dare.hpp>
#include <dcs/debug.hpp>
#include <dcs/test.hpp>
#include <iostream>
#include "./utility.hpp"

namespace ublas = boost::numeric::ublas;
namespace ublasx = boost::numeric::ublasx;
namespace dcs_control = dcs::control;


const double tol = 1e-5;

//namespace detail { namespace /*<unnamed>*/ {
//
//template <typename T>
//struct complex_cmp
//{
//	bool operator()(::std::complex<T> const& a, ::std::complex<T> const& b) const
//	{
//		return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag());
//	}
//};
//
//}} // Namespace detail::<unnamed>*/


DCS_TEST_DEF( test_1 )
{
	DCS_DEBUG_TRACE("Test Case: Test #1");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	const std::size_t n(2);
	const std::size_t m(1);

	matrix_type A(n,n);
	A(0,0) = 1; A(0,1) = 0.786939;
	A(1,0) = 0; A(1,1) = 0.606531;

	matrix_type B(n,m);
	B(0,0) = 0.426123;
	B(1,0) = 0.786939;

	matrix_type Q(n,n);
	Q(0,0) = 1; Q(0,1) = 0;
	Q(1,0) = 0; Q(1,1) = 1;

	matrix_type R;
	R = ublas::identity_matrix<value_type>(m);

	matrix_type S;
	S = ublas::zero_matrix<value_type>(n,m);

	matrix_type E;
	E = ublas::identity_matrix<value_type>(n);

	matrix_type X;

	matrix_type expect_X(n,n);
	expect_X(0,0) = 2.401533210017093; expect_X(0,1) = 1.057915121181177;
	expect_X(1,0) = 1.057915121181177; expect_X(1,1) = 2.097053494085709;

	dcs_control::dare(A, B, Q, R, S, E, X);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "S = " << S );
	DCS_DEBUG_TRACE( "E = " << E );
	DCS_DEBUG_TRACE( "Solution X = " << X );
	DCS_TEST_CHECK_MATRIX_CLOSE(X, expect_X, n, n, tol);
}


DCS_TEST_DEF( test_2 )
{
	DCS_DEBUG_TRACE("Test Case: Test #2");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(2);
	const std::size_t m(1);

	matrix_type A(n,n);
	A(0,0) = 1; A(0,1) = 0.786939;
	A(1,0) = 0; A(1,1) = 0.606531;

	matrix_type B(n,m);
	B(0,0) = 0.426123;
	B(1,0) = 0.786939;

	matrix_type Q(n,n);
	Q(0,0) = 1; Q(0,1) = 0;
	Q(1,0) = 0; Q(1,1) = 1;

	matrix_type R;
	R = ublas::identity_matrix<value_type>(m);

	matrix_type S;
	S = ublas::zero_matrix<value_type>(n,m);

	matrix_type E;
	E = ublas::identity_matrix<value_type>(n);

	matrix_type expect_X(n,n);
	expect_X(0,0) = 2.401533210017093; expect_X(0,1) = 1.057915121181177;
	expect_X(1,0) = 1.057915121181177; expect_X(1,1) = 2.097053494085709;

//	vector_type expect_L(n);
//	expect_L(0) = std::complex<value_type>(0.376036021533870, 0.186272944433102);
//	expect_L(1) = std::complex<value_type>(0.376036021533870,-0.186272944433102);

	matrix_type expect_G(m,n);
	expect_G(0,0) = 0.538832818098109; expect_G(0,1) = 0.794025839341854;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R, S, E);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "S = " << S );
	DCS_DEBUG_TRACE( "E = " << E );
	DCS_DEBUG_TRACE( "Solution X = " << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
}


DCS_TEST_DEF( test_3 )
{
	DCS_DEBUG_TRACE("Test Case: Test #3");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< ::std::complex<value_type> > vector_type;

	const std::size_t n(2);
	const std::size_t m(1);

	matrix_type A(n,n);
	A(0,0) = 0.4; A(0,1) = 1.7;
	A(1,0) = 0.9; A(1,1) = 3.8;

	matrix_type B(n,m);
	B(0,0) = 0.8;
	B(1,0) = 2.1;

	matrix_type Q(n,n);
	Q(0,0) =  1; Q(0,1) = -1;
	Q(1,0) = -1; Q(1,1) =  1;

	matrix_type R;
	R = ublas::scalar_matrix<value_type>(m, m, 3);

	matrix_type expect_X(n,n);
	expect_X(0,0) = 1.53536988885652; expect_X(0,1) = 1.26226249343861;
	expect_X(1,0) = 1.26226249343861; expect_X(1,1) = 10.5595806602029;

	vector_type expect_l(n);
	expect_l(0) = ::std::complex<value_type>(-0.0022307284532799,0);
	expect_l(1) = ::std::complex<value_type>( 0.245448650750199,0);

	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.409151396044144; expect_G(0,1) = 1.72831474327036;


	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );

	vector_type l(solver.eigenvalues());
	::std::sort(l.begin(), l.end(), detail::complex_cmp<value_type>());
	::std::sort(expect_l.begin(), expect_l.end(), detail::complex_cmp<value_type>());
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(l, expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_1 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.1");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(2);
	const std::size_t m(1);

	matrix_type A(n,n);
	A(0,0) = 2; A(0,1) = -1;
	A(1,0) = 1; A(1,1) =  0;

	matrix_type B(n,m);
	B(0,0) = 1;
	B(1,0) = 0;

	matrix_type Q(n,n);
	Q(0,0) = 0; Q(0,1) = 0;
	Q(1,0) = 0; Q(1,1) = 1;

	matrix_type R(m,m);
	R(0,0) = 0;

	matrix_type expect_X(n,n);
	expect_X(0,0) = 1.0; expect_X(0,1) = 0.0;
	expect_X(1,0) = 0.0; expect_X(1,1) = 1.0;

	vector_type expect_l(n);
	expect_l(0) = 0.0;
	expect_l(1) = 0.0;

	matrix_type expect_G(m,n);
	expect_G(0,0) = 2.0; expect_G(0,1) = -1.0;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_2 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.2");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(2);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) = 0; A(0,1) =  1;
	A(1,0) = 0; A(1,1) = -1;

	matrix_type B(n,m);
	B(0,0) = 1; B(0,1) = 0;
	B(1,0) = 2; B(1,1) = 1;

	matrix_type Q(n,n);
	Q(0,0) = -4.0/11.0; Q(0,1) = -4.0/11.0;
	Q(1,0) = -4.0/11.0; Q(1,1) =  7.0/11.0;

	matrix_type R(m,m);
	R(0,0) = 9; R(0,1) = 3;
	R(1,0) = 3; R(1,1) = 1;

	matrix_type S(n,m);
	S(0,0) =  3; S(0,1) = 1;
	S(1,0) = -1; S(1,1) = 7;

	matrix_type expect_X(n,n);
	expect_X(0,0) = -0.014021341244239*1e+2; expect_X(0,1) =  0.130568663991581*1e+2;
	expect_X(1,0) =  0.130568663991581*1e+2; expect_X(1,1) = -1.256364927952908*1e+2;

	vector_type expect_l(n);
	expect_l(0) = -0.217058149756746;
	expect_l(1) =  0.687271691663810;

	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.940453958559404; expect_G(0,1) = -11.009835262327730;
	expect_G(1,0) = -1.782864114890681; expect_G(1,1) =  19.609003024189001;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R, S);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "S = " << S );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_3 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.3");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(2);
	const std::size_t m(1);

	matrix_type A(n,n);
	A(0,0) = 0; A(0,1) = 1;
	A(1,0) = 0; A(1,1) = 0;

	matrix_type B(n,m);
	B(0,0) = 0;
	B(1,0) = 1;

	matrix_type Q(n,n);
	Q(0,0) = 1; Q(0,1) = 2;
	Q(1,0) = 2; Q(1,1) = 4;

	matrix_type R(m,m);
	R(0,0) = 1;

	matrix_type expect_X(n,n);
	expect_X(0,0) = 1; expect_X(0,1) = 2;
	expect_X(1,0) = 2; expect_X(1,1) = 2+std::sqrt(5);

	vector_type expect_l(n);
	expect_l(0) =   0.0;
	expect_l(1) =  -0.381966011250105;

	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.0; expect_G(0,1) = 0.381966011250105;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_4 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.4");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(3);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) = 0.0; A(0,1) = 0.1; A(0,2) = 0.0;
	A(1,0) = 0.0; A(1,1) = 0.0; A(1,2) = 0.1;
	A(2,0) = 0.0; A(2,1) = 0.0; A(2,2) = 0.0;

	matrix_type B(n,m);
	B(0,0) = 1; B(0,1) = 0;
	B(1,0) = 0; B(1,1) = 0;
	B(2,0) = 0; B(2,1) = 1;

	matrix_type Q(n,n);
	Q(0,0) = 1e+5; Q(0,1) =    0; Q(0,2) =   0;
	Q(1,0) =    0; Q(1,1) = 1e+3; Q(1,2) =   0;
	Q(2,0) =    0; Q(2,1) =    0; Q(2,2) = -10;

	matrix_type R(m,m);
	R(0,0) = 0; R(0,1) = 0;
	R(1,0) = 0; R(1,1) = 1;


	matrix_type expect_X(n,n);
	expect_X(0,0) = 1e+5; expect_X(0,1) =    0; expect_X(0,2) = 0;
	expect_X(1,0) =    0; expect_X(1,1) = 1e+3; expect_X(1,2) = 0;
	expect_X(2,0) =    0; expect_X(2,1) =    0; expect_X(2,2) = 0;

	vector_type expect_l(n);
	//This is the result we got with MATLAB 2008a
	//expect_l(0) =  -1.0e-14*0.195676808090185;
	//expect_l(1) =   1.0e-14*0.0;
	//expect_l(2) =   1.0e-14*0.0;
	//Instead, this is the one we got with SLICOT 5.0
	expect_l(0) = 0;
	expect_l(1) = 0;
	expect_l(2) = 0;

	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.0; expect_G(0,1) = 0.100000000000000; expect_G(0,2) =  0;
	expect_G(1,0) =  0.0; expect_G(1,1) = 0.0;               expect_G(1,2) = -0.000000000000000;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_5 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.5");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(4);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) =  0.998; A(0,1) = 0.067; A(0,2) =  0.000; A(0,3) = 0.000;
	A(1,0) = -0.067; A(1,1) = 0.998; A(1,2) =  0.100; A(1,3) = 0.000;
	A(2,0) =  0.000; A(2,1) = 0.000; A(2,2) =  0.998; A(2,3) = 0.153;
	A(3,0) =  0.000; A(3,1) = 0.000; A(3,2) = -0.153; A(3,3) = 0.998;

	matrix_type B(n,m);
	B(0,0) =  0.0033; B(0,1) =  0.0200;
	B(1,0) =  0.1000; B(1,1) = -0.0007;
	B(2,0) =  0.0400; B(2,1) =  0.0073;
	B(3,0) = -0.0028; B(3,1) =  0.1000;

	matrix_type Q(n,n);
	Q(0,0) =  1.870; Q(0,1) = 0.000; Q(0,2) = 0.000; Q(0,3) = -0.244;
	Q(1,0) =  0.000; Q(1,1) = 0.744; Q(1,2) = 0.205; Q(1,3) =  0.000;
	Q(2,0) =  0.000; Q(2,1) = 0.205; Q(2,2) = 0.589; Q(2,3) =  0.000;
	Q(3,0) = -0.244; Q(3,1) = 0.000; Q(3,2) = 0.000; Q(3,3) =  1.048;

	matrix_type R(m,m);
	R(0,0) = 1; R(0,1) = 0;
	R(1,0) = 0; R(1,1) = 1;


	matrix_type expect_X(n,n);
	expect_X(0,0) =  30.707390002658851; expect_X(0,1) =  7.731389771619265; expect_X(0,2) =  3.966329567211186; expect_X(0,3) = - 4.901197596654572;
	expect_X(1,0) =   7.731389771619265; expect_X(1,1) = 11.829796382196321; expect_X(1,2) =  5.164569890757084; expect_X(1,3) =   0.278956010969016;
	expect_X(2,0) =   3.966329567211186; expect_X(2,1) =  5.164569890757084; expect_X(2,2) = 17.132194857924887; expect_X(2,3) =   1.573172972387152;
	expect_X(3,0) = - 4.901197596654572; expect_X(3,1) =  0.278956010969016; expect_X(3,2) =  1.573172972387152; expect_X(3,3) =  14.880017305642824;

	vector_type expect_l(n);
	expect_l(0) = std::complex<value_type>(0.924483957358999, 0.065175187407696);
	expect_l(1) = std::complex<value_type>(0.924483957358999,-0.065175187407696);
	expect_l(2) = std::complex<value_type>(0.921554823634630, 0.141844900635441);
	expect_l(3) = std::complex<value_type>(0.921554823634630,-0.141844900635441);

	matrix_type expect_G(m,n);
	expect_G(0,0) = 0.793645328788861; expect_G(0,1) = 1.237433329574766; expect_G(0,2) = 1.123694684785059; expect_G(0,3) = 0.148799363279981;
	expect_G(1,0) = 0.093940974504011; expect_G(1,1) = 0.158621967953197; expect_G(1,2) = 0.111849254879984; expect_G(1,3) = 1.264446426229076;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_6 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.6");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(4);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) =  1e-3* 984.750; A(0,1) = -1e-3* 79.903; A(0,2) =  1e-3*  0.9054; A(0,3) = -1e-3*  1.0765;
	A(1,0) =  1e-3*  41.588; A(1,1) =  1e-3*998.990; A(1,2) = -1e-3* 35.8550; A(1,3) =  1e-3* 12.6840;
	A(2,0) = -1e-3* 546.620; A(2,1) =  1e-3* 44.916; A(2,2) = -1e-3*329.9100; A(2,3) =  1e-3*193.1800;
	A(3,0) =  1e-3*2662.400; A(3,1) = -1e-3*100.450; A(3,2) = -1e-3*924.5500; A(3,3) = -1e-3*263.2500;

	matrix_type B(n,m);
	B(0,0) =  1e-4*   37.112; B(0,1) = 1e-4*7.361000;
	B(1,0) = -1e-4*  870.510; B(1,1) = 1e-4*0.093411;
	B(2,0) = -1e-4*11984.000; B(2,1) = 1e-4*4.137800;
	B(3,0) = -1e-4*31927.000; B(3,1) = 1e-4*9.253500;

	matrix_type Q(n,n);
	Q(0,0) = 0.01; Q(0,1) = 0.00; Q(0,2) = 0.00; Q(0,3) = 0.00;
	Q(1,0) = 0.00; Q(1,1) = 0.01; Q(1,2) = 0.00; Q(1,3) = 0.00;
	Q(2,0) = 0.00; Q(2,1) = 0.00; Q(2,2) = 0.01; Q(2,3) = 0.00;
	Q(3,0) = 0.00; Q(3,1) = 0.00; Q(3,2) = 0.00; Q(3,3) = 0.01;

	matrix_type R(m,m);
	R(0,0) = 1; R(0,1) = 0;
	R(1,0) = 0; R(1,1) = 1;


	matrix_type expect_X(n,n);
	expect_X(0,0) =  1.845979749417746; expect_X(0,1) = -0.056190418825482; expect_X(0,2) = -0.011141163569509; expect_X(0,3) = -0.010550944100757;
	expect_X(1,0) = -0.056190418825482; expect_X(1,1) =  2.047915333938012; expect_X(1,2) = -0.059430745710932; expect_X(1,3) =  0.012336403378050;
	expect_X(2,0) = -0.011141163569509; expect_X(2,1) = -0.059430745710932; expect_X(2,2) =  0.022799689365452; expect_X(2,3) =  0.000808498679231;
	expect_X(3,0) = -0.010550944100757; expect_X(3,1) =  0.012336403378050; expect_X(3,2) =  0.000808498679231; expect_X(3,3) =  0.011514128775705;


	vector_type expect_l(n);
	expect_l(0) = std::complex<value_type>(-0.265900128787536,+0.397166041362223);
	expect_l(1) = std::complex<value_type>(-0.265900128787536,-0.397166041362223);
	expect_l(2) = std::complex<value_type>( 0.985894853633392,+0.074733641936961);
	expect_l(3) = std::complex<value_type>( 0.985894853633392,-0.074733641936961);


	matrix_type expect_G(m,n);
	expect_G(0,0) = -0.032648216905453; expect_G(0,1) = -0.127202840052208; expect_G(0,2) =  0.042338962051211; expect_G(0,3) =   0.003014183257250;
	expect_G(1,0) =  0.001329925628570; expect_G(1,1) = -0.000143573351369; expect_G(1,2) = -0.000000888234976; expect_G(1,3) =  -0.000002486099077;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_7 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.7");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(4);
	const std::size_t m(4);

	matrix_type A(n,n);
	A(0,0) = -0.6; A(0,1) = -2.2; A(0,2) = -3.6; A(0,3) = -5.400018;
	A(1,0) =  1.0; A(1,1) =  0.6; A(1,2) =  0.8; A(1,3) =  3.399982;
	A(2,0) =  0.0; A(2,1) =  1.0; A(2,2) =  1.8; A(2,3) =  3.799982;
	A(3,0) =  0.0; A(3,1) =  0.0; A(3,2) =  0.0; A(3,3) = -0.999982;

	matrix_type B(n,m);
	B(0,0) =  1; B(0,1) = -1; B(0,2) = -1; B(0,3) = -1;
	B(1,0) =  0; B(1,1) =  1; B(1,2) = -1; B(1,3) = -1;
	B(2,0) =  0; B(2,1) =  0; B(2,2) =  1; B(2,3) = -1;
	B(3,0) =  0; B(3,1) =  0; B(3,2) =  0; B(3,3) =  1;

	matrix_type Q(n,n);
	Q(0,0) = 2; Q(0,1) = 1; Q(0,2) =  3; Q(0,3) =  6;
	Q(1,0) = 1; Q(1,1) = 2; Q(1,2) =  2; Q(1,3) =  5;
	Q(2,0) = 3; Q(2,1) = 2; Q(2,2) =  6; Q(2,3) = 11;
	Q(3,0) = 6; Q(3,1) = 5; Q(3,2) = 11; Q(3,3) = 22;

	matrix_type R(m,m);
	R(0,0) = 1; R(0,1) = 0; R(0,2) = 0; R(0,3) = 0;
	R(1,0) = 0; R(1,1) = 1; R(1,2) = 0; R(1,3) = 0;
	R(2,0) = 0; R(2,1) = 0; R(2,2) = 1; R(2,3) = 0;
	R(3,0) = 0; R(3,1) = 0; R(3,2) = 0; R(3,3) = 1;


	matrix_type expect_X(n,n);
	expect_X(0,0) =  2.817800285785478; expect_X(0,1) =  2.210160220269767; expect_X(0,2) =  4.997305238829639; expect_X(0,3) = 10.025262087702128;
	expect_X(1,0) =  2.210160220269767; expect_X(1,1) =  4.530484712640063; expect_X(1,2) =  6.257605864699351; expect_X(1,3) = 12.998245151137297;
	expect_X(2,0) =  4.997305238829639; expect_X(2,1) =  6.257605864699351; expect_X(2,2) = 13.192802180652688; expect_X(2,3) = 24.447703263082389;
	expect_X(3,0) = 10.025262087702128; expect_X(3,1) = 12.998245151137297; expect_X(3,2) = 24.447703263082389; expect_X(3,3) = 47.471191943996843;

	vector_type expect_l(n);
	expect_l(0) = 0.067912325238694;
	expect_l(1) = 0.163207174929605;
	expect_l(2) = 0.348653414955625;
	expect_l(3) = -0.999981999984721;


	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.250004803979136; expect_G(0,1) =  0.216553352201969; expect_G(0,2) =  0.460193182636258; expect_G(0,3) =  0.926746747874340;
	expect_G(1,0) =  0.717798364193832; expect_G(1,1) =  1.123538879388990; expect_G(1,2) =  1.813227965775154; expect_G(1,3) =  3.654563388552431;
	expect_G(2,0) = -0.038319084032009; expect_G(2,1) =  0.646201164736901; expect_G(2,2) =  1.172363846404621; expect_G(2,3) =  1.780245030303713;
	expect_G(3,0) = -0.000001137431440; expect_G(3,1) = -0.000001872546490; expect_G(3,2) = -0.000003258144388; expect_G(3,3) = -0.000007034919287;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	// NOTE: With a tolerance of 10^{-5} some check on the X matrix fails
	//       This may be caused by the difference in the algorithms used by
	//       our library and MATLAB.
	//       So we decided to lower the precision to 10^{-4}
	const double tol2 = 1e-4;

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol2);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_8 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.8");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(5);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) = 1e-4*9540.70; A(0,1) = 1e-4* 196.43; A(0,2) = 1e-4*  35.97; A(0,3) = 1e-4*   6.73; A(0,4) = 1e-4*   1.90;
	A(1,0) = 1e-4*4084.90; A(1,1) = 1e-4*4131.70; A(1,2) = 1e-4*1608.40; A(1,3) = 1e-4* 446.79; A(1,4) = 1e-4* 119.71;
	A(2,0) = 1e-4*1221.70; A(2,1) = 1e-4*2632.60; A(2,2) = 1e-4*3614.90; A(2,3) = 1e-4*1593.00; A(2,4) = 1e-4*1238.30;
	A(3,0) = 1e-4* 411.18; A(3,1) = 1e-4*1285.80; A(3,2) = 1e-4*2720.90; A(3,3) = 1e-4*2144.20; A(3,4) = 1e-4*4097.60;
	A(4,0) = 1e-4*  13.05; A(4,1) = 1e-4*  58.08; A(4,2) = 1e-4* 187.50; A(4,3) = 1e-4* 361.62; A(4,4) = 1e-4*9428.00;

	matrix_type B(n,m);
	B(0,0) = 1e-4*  4.34; B(0,1) = -1e-4*  1.22;
	B(1,0) = 1e-4*266.06; B(1,1) = -1e-4*104.53;
	B(2,0) = 1e-4*375.30; B(2,1) = -1e-4*551.00;
	B(3,0) = 1e-4*360.76; B(3,1) = -1e-4*660.00;
	B(4,0) = 1e-4* 46.17; B(4,1) = -1e-4* 91.48;

	matrix_type Q(n,n);
	Q(0,0) = 1; Q(0,1) = 0; Q(0,2) = 0; Q(0,3) = 0; Q(0,4) = 0;
	Q(1,0) = 0; Q(1,1) = 1; Q(1,2) = 0; Q(1,3) = 0; Q(1,4) = 0;
	Q(2,0) = 0; Q(2,1) = 0; Q(2,2) = 1; Q(2,3) = 0; Q(2,4) = 0;
	Q(3,0) = 0; Q(3,1) = 0; Q(3,2) = 0; Q(3,3) = 1; Q(3,4) = 0;
	Q(4,0) = 0; Q(4,1) = 0; Q(4,2) = 0; Q(4,3) = 0; Q(4,4) = 1;

	matrix_type R(m,m);
	R(0,0) = 1; R(0,1) = 0;
	R(1,0) = 0; R(1,1) = 1;


	matrix_type expect_X(n,n);
	expect_X(0,0) = 60.456378667853016; expect_X(0,1) = 4.672669320499056; expect_X(0,2) = 3.391010077549930; expect_X(0,3) = 2.183442069546651; expect_X(0,4) = 24.049867925873606;
	expect_X(1,0) =  4.672669320499056; expect_X(1,1) = 1.800604768836491; expect_X(1,2) = 0.643306876604517; expect_X(1,3) = 0.368582157291515; expect_X(1,4) =  2.810035381375605;
	expect_X(2,0) =  3.391010077549930; expect_X(2,1) = 0.643306876604517; expect_X(2,2) = 1.617206256697997; expect_X(2,3) = 0.373785641698017; expect_X(2,4) =  2.685608599483135;
	expect_X(3,0) =  2.183442069546651; expect_X(3,1) = 0.368582157291515; expect_X(3,2) = 0.373785641698017; expect_X(3,3) = 1.246083057616368; expect_X(3,4) =  2.056868888863928;
	expect_X(4,0) = 24.049867925873606; expect_X(4,1) = 2.810035381375605; expect_X(4,2) = 2.685608599483135; expect_X(4,3) = 2.056868888863928; expect_X(4,4) = 27.429360377605320;

	vector_type expect_l(n);
	expect_l(0) = 0.240702660917670;
	expect_l(1) = 0.050274226759839;
	expect_l(2) = 0.632461076685746;
	expect_l(3) = 0.976994439625731;
	expect_l(4) = 0.958597232843032;

	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.488379469793936; expect_G(0,1) =  0.088047339574689; expect_G(0,2) =  0.081945427902266; expect_G(0,3) =  0.050798306322876; expect_G(0,4) =  0.397748565018305;
	expect_G(1,0) = -0.624013969453445; expect_G(1,1) = -0.109166256049604; expect_G(1,2) = -0.112656060538442; expect_G(1,3) = -0.073282793503359; expect_G(1,4) = -0.585936198398728;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_9 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.9");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(6);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) = 0; A(0,1) = 1; A(0,2) = 0; A(0,3) = 0; A(0,4) = 0; A(0,5) = 0;
	A(1,0) = 0; A(1,1) = 0; A(1,2) = 1; A(1,3) = 0; A(1,4) = 0; A(1,5) = 0;
	A(2,0) = 0; A(2,1) = 0; A(2,2) = 0; A(2,3) = 0; A(2,4) = 0; A(2,5) = 0;
	A(3,0) = 0; A(3,1) = 0; A(3,2) = 0; A(3,3) = 0; A(3,4) = 1; A(3,5) = 0;
	A(4,0) = 0; A(4,1) = 0; A(4,2) = 0; A(4,3) = 0; A(4,4) = 0; A(4,5) = 1;
	A(5,0) = 0; A(5,1) = 0; A(5,2) = 0; A(5,3) = 0; A(5,4) = 0; A(5,5) = 0;

	matrix_type B(n,m);
	B(0,0) = 0; B(0,1) = 0;
	B(1,0) = 0; B(1,1) = 0;
	B(2,0) = 1; B(2,1) = 0;
	B(3,0) = 0; B(3,1) = 0;
	B(4,0) = 0; B(4,1) = 0;
	B(5,0) = 0; B(5,1) = 1;

	matrix_type Q(n,n);
	Q(0,0) = 1; Q(0,1) = 1; Q(0,2) = 0; Q(0,3) =  0; Q(0,4) =  0; Q(0,5) = 0;
	Q(1,0) = 1; Q(1,1) = 1; Q(1,2) = 0; Q(1,3) =  0; Q(1,4) =  0; Q(1,5) = 0;
	Q(2,0) = 0; Q(2,1) = 0; Q(2,2) = 0; Q(2,3) =  0; Q(2,4) =  0; Q(2,5) = 0;
	Q(3,0) = 0; Q(3,1) = 0; Q(3,2) = 0; Q(3,3) =  1; Q(3,4) = -1; Q(3,5) = 0;
	Q(4,0) = 0; Q(4,1) = 0; Q(4,2) = 0; Q(4,3) = -1; Q(4,4) =  1; Q(4,5) = 0;
	Q(5,0) = 0; Q(5,1) = 0; Q(5,2) = 0; Q(5,3) =  0; Q(5,4) =  0; Q(5,5) = 0;

	matrix_type R(m,m);
	R(0,0) = 3; R(0,1) = 0;
	R(1,0) = 0; R(1,1) = 1;

	matrix_type S(n,m);
	S(0,0) =  1; S(0,1) = 0;
	S(1,0) =  1; S(1,1) = 0;
	S(2,0) =  0; S(2,1) = 0;
	S(3,0) =  1; S(3,1) = 0;
	S(4,0) = -1; S(4,1) = 0;
	S(5,0) =  0; S(5,1) = 0;


	matrix_type expect_X(n,n);
	expect_X(0,0) =  0.776931379296968; expect_X(0,1) =  0.810457127580939; expect_X(0,2) = -0.150367062101950; expect_X(0,3) = -0.223068620703031; expect_X(0,4) =  0.256594368987001; expect_X(0,5) = -0.002115177331576;
	expect_X(1,0) =  0.810457127580939; expect_X(1,1) =  1.615873554290608; expect_X(1,2) =  0.682929643886609; expect_X(1,3) = -0.189542872419061; expect_X(1,4) = -0.005040700871270; expect_X(1,5) =  0.254095399301336;
	expect_X(2,0) = -0.150367062101950; expect_X(2,1) =  0.682929643886609; expect_X(2,2) =  1.485634705307614; expect_X(2,3) = -0.150367062101950; expect_X(2,4) = -0.016336231909491; expect_X(2,5) =  0.077827181137999;
	expect_X(3,0) = -0.223068620703031; expect_X(3,1) = -0.189542872419061; expect_X(3,2) = -0.150367062101950; expect_X(3,3) =  0.776931379296968; expect_X(3,4) = -0.743405631012998; expect_X(3,5) = -0.002115177331577;
	expect_X(4,0) =  0.256594368987001; expect_X(4,1) = -0.005040700871270; expect_X(4,2) = -0.016336231909491; expect_X(4,3) = -0.743405631012998; expect_X(4,4) =  1.481770561154726; expect_X(4,5) = -0.741674246035512;
	expect_X(5,0) = -0.002115177331576; expect_X(5,1) =  0.254095399301336; expect_X(5,2) =  0.077827181137999; expect_X(5,3) = -0.002115177331577; expect_X(5,4) = -0.741674246035512; expect_X(5,5) =  1.235707250512302;

	vector_type expect_l(n);
	expect_l(0) = std::complex<value_type>( 0.195825301202028, 0.642361399464209);
	expect_l(1) = std::complex<value_type>( 0.195825301202028,-0.642361399464209);
	expect_l(2) = std::complex<value_type>( 0.376862701311676, 0.000000000000000);
	expect_l(3) = std::complex<value_type>(-0.000000000000001, 0.000000013433172);
	expect_l(4) = std::complex<value_type>(-0.000000000000001,-0.000000013433172);
	expect_l(5) = std::complex<value_type>(-0.587066411126228, 0.000000000000000);

	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.223068620703031; expect_G(0,1) =  0.189542872419061; expect_G(0,2) =  0.150367062101950; expect_G(0,3) =  0.223068620703031; expect_G(0,4) = -0.256594368987002; expect_G(0,5) =  0.002115177331576;
	expect_G(1,0) = -0.007765239364716; expect_G(1,1) = -0.007544263584996; expect_G(1,2) =  0.108418825705555; expect_G(1,3) = -0.007765239364716; expect_G(1,4) =  0.007986215144435; expect_G(1,5) = -0.331813954691450;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R, S);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "S = " << S );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_darex_ex_1_10 )
{
	DCS_DEBUG_TRACE("Test Case: DAREX Ex. 1.10");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector< std::complex<value_type> > vector_type;

	const std::size_t n(9);
	const std::size_t m(3);

	matrix_type A(n,n);
	A(0,0) =  1e-2*87.010; A(0,1) = 1e-2* 13.500; A(0,2) =  1e-2* 1.159; A(0,3) =  1e-2* 0.05014; A(0,4) = -1e-2* 3.722; A(0,5) = 1e-2*0.03484; A(0,6) = 1e-2*0.00000; A(0,7) = 1e-2* 0.4242; A(0,8) = 1e-2* 0.7249;
	A(1,0) =  1e-2* 7.655; A(1,1) = 1e-2* 89.740; A(1,2) =  1e-2* 1.272; A(1,3) =  1e-2* 0.05504; A(1,4) = -1e-2* 4.016; A(1,5) = 1e-2*0.03743; A(1,6) = 1e-2*0.00000; A(1,7) = 1e-2* 0.4530; A(1,8) = 1e-2* 0.7499;
	A(2,0) = -1e-2*12.720; A(2,1) = 1e-2* 35.750; A(2,2) =  1e-2*81.700; A(2,3) =  1e-2* 0.14550; A(2,4) = -1e-2*10.280; A(2,5) = 1e-2*0.09870; A(2,6) = 1e-2*0.00000; A(2,7) = 1e-2* 1.1850; A(2,8) = 1e-2* 1.8720;
	A(3,0) = -1e-2*36.350; A(3,1) = 1e-2* 63.390; A(3,2) =  1e-2* 7.491; A(3,3) =  1e-2*79.66000; A(3,4) = -1e-2*27.350; A(3,5) = 1e-2*0.26530; A(3,6) = 1e-2*0.00000; A(3,7) = 1e-2* 3.1720; A(3,8) = 1e-2* 4.8820;
	A(4,0) = -1e-2*96.000; A(4,1) = 1e-2*164.590; A(4,2) = -1e-2*12.890; A(4,3) = -1e-2* 0.55970; A(4,4) =  1e-2* 7.142; A(4,5) = 1e-2*0.71080; A(4,6) = 1e-2*0.00000; A(4,7) = 1e-2* 8.4520; A(4,8) = 1e-2*12.5900;
	A(5,0) = -1e-2*66.440; A(5,1) = 1e-2* 11.296; A(5,2) = -1e-2* 8.889; A(5,3) = -1e-2* 0.38540; A(5,4) =  1e-2* 8.447; A(5,5) = 1e-2*1.36000; A(5,6) = 1e-2*0.00000; A(5,7) = 1e-2*14.4300; A(5,8) = 1e-2*10.1600;
	A(6,0) =  1e-2*41.020; A(6,1) = 1e-2* 69.300; A(6,2) = -1e-2* 5.471; A(6,3) = -1e-2* 0.23710; A(6,4) =  1e-2* 6.649; A(6,5) = 1e-2*1.24900; A(6,6) = 1e-2*0.01063; A(6,7) = 1e-2* 9.9970; A(6,8) = 1e-2* 6.9670;
	A(7,0) = -1e-2*17.990; A(7,1) = 1e-2* 30.170; A(7,2) = -1e-2* 2.393; A(7,3) = -1e-2* 0.10350; A(7,4) =  1e-2* 6.059; A(7,5) = 1e-2*2.21600; A(7,6) = 1e-2*0.00000; A(7,7) = 1e-2*21.3900; A(7,8) = 1e-2* 3.5540;
	A(8,0) = -1e-2*34.510; A(8,1) = 1e-2* 58.040; A(8,2) = -1e-2* 4.596; A(8,3) = -1e-2* 0.19890; A(8,4) =  1e-2*10.560; A(8,5) = 1e-2*1.98600; A(8,6) = 1e-2*0.00000; A(8,7) = 1e-2*21.9100; A(8,8) = 1e-2*21.5200;

	matrix_type B(n,m);
	B(0,0) = 1e-4* 4.76000; B(0,1) = -1e-4* 0.57010; B(0,2) = -1e-4*83.68000;
	B(1,0) = 1e-4* 0.87900; B(1,1) = -1e-4* 4.77300; B(1,2) = -1e-4* 2.73000;
	B(2,0) = 1e-4* 1.48200; B(2,1) = -1e-4*13.12000; B(2,2) =  1e-4* 8.87600;
	B(3,0) = 1e-4* 3.89200; B(3,1) = -1e-4*35.13000; B(3,2) =  1e-4*24.80000;
	B(4,0) = 1e-4*10.34000; B(4,1) = -1e-4*92.75000; B(4,2) =  1e-4*66.80000;
	B(5,0) = 1e-4* 7.20300; B(5,1) = -1e-4*61.59000; B(5,2) =  1e-4*38.34000;
	B(6,0) = 1e-4* 4.45400; B(6,1) = -1e-4*36.83000; B(6,2) =  1e-4*20.29000;
	B(7,0) = 1e-4* 1.97100; B(7,1) = -1e-4*15.54000; B(7,2) =  1e-4* 6.93700;
	B(8,0) = 1e-4* 3.77300; B(8,1) = -1e-4*30.28000; B(8,2) =  1e-4*14.69000;

	matrix_type Q(n,n);
    Q(0,0) = 50; Q(0,1) = 0; Q(0,2) = 0; Q(0,3) = 0; Q(0,4) =  0; Q(0,5) = 0; Q(0,6) = 0; Q(0,7) = 0; Q(0,8) = 0;
    Q(1,0) =  0; Q(1,1) = 0; Q(1,2) = 0; Q(1,3) = 0; Q(1,4) =  0; Q(1,5) = 0; Q(1,6) = 0; Q(1,7) = 0; Q(1,8) = 0;
    Q(2,0) =  0; Q(2,1) = 0; Q(2,2) = 0; Q(2,3) = 0; Q(2,4) =  0; Q(2,5) = 0; Q(2,6) = 0; Q(2,7) = 0; Q(2,8) = 0;
    Q(3,0) =  0; Q(3,1) = 0; Q(3,2) = 0; Q(3,3) = 0; Q(3,4) =  0; Q(3,5) = 0; Q(3,6) = 0; Q(3,7) = 0; Q(3,8) = 0;
    Q(4,0) =  0; Q(4,1) = 0; Q(4,2) = 0; Q(4,3) = 0; Q(4,4) = 50; Q(4,5) = 0; Q(4,6) = 0; Q(4,7) = 0; Q(4,8) = 0;
    Q(5,0) =  0; Q(5,1) = 0; Q(5,2) = 0; Q(5,3) = 0; Q(5,4) =  0; Q(5,5) = 0; Q(5,6) = 0; Q(5,7) = 0; Q(5,8) = 0;
    Q(6,0) =  0; Q(6,1) = 0; Q(6,2) = 0; Q(6,3) = 0; Q(6,4) =  0; Q(6,5) = 0; Q(6,6) = 0; Q(6,7) = 0; Q(6,8) = 0;
    Q(7,0) =  0; Q(7,1) = 0; Q(7,2) = 0; Q(7,3) = 0; Q(7,4) =  0; Q(7,5) = 0; Q(7,6) = 0; Q(7,7) = 0; Q(7,8) = 0;
    Q(8,0) =  0; Q(8,1) = 0; Q(8,2) = 0; Q(8,3) = 0; Q(8,4) =  0; Q(8,5) = 0; Q(8,6) = 0; Q(8,7) = 0; Q(8,8) = 0;


	matrix_type R(m,m);
	R(0,0) = 1; R(0,1) = 0; R(0,2) = 0;
	R(1,0) = 0; R(1,1) = 1; R(1,2) = 0;
	R(2,0) = 0; R(2,1) = 0; R(2,2) = 1;


	matrix_type expect_X(n,n);
	expect_X(0,0) =  1e+2*5.194221256889423; expect_X(0,1) =  1e+2*0.017093107843102; expect_X(0,2) =  1e+2*0.597943238405564; expect_X(0,3) =  1e+2*0.026730739405916; expect_X(0,4) = -1e+2*0.313928655928748; expect_X(0,5) = -1e+2*0.004690645261773; expect_X(0,6) = 1e+2*0; expect_X(0,7) = -1e+2*0.052623035300436; expect_X(0,8) = -1e+2*0.057714706758673;
	expect_X(1,0) =  1e+2*0.017093107843102; expect_X(1,1) =  1e+2*6.062940897213983; expect_X(1,2) =  1e+2*0.026906043422892; expect_X(1,3) =  1e+2*0.001396884300817; expect_X(1,4) = -1e+2*0.163784382994394; expect_X(1,5) =  1e+2*0.012779657048731; expect_X(1,6) = 1e+2*0; expect_X(1,7) =  1e+2*0.143430702460303; expect_X(1,8) =  1e+2*0.164521508295555;
	expect_X(2,0) =  1e+2*0.597943238405564; expect_X(2,1) =  1e+2*0.026906043422892; expect_X(2,2) =  1e+2*0.094662957940239; expect_X(2,3) =  1e+2*0.004250513282368; expect_X(2,4) = -1e+2*0.044326303597952; expect_X(2,5) = -1e+2*0.000662424480368; expect_X(2,6) = 1e+2*0; expect_X(2,7) = -1e+2*0.007428328837057; expect_X(2,8) = -1e+2*0.008137995705764;
	expect_X(3,0) =  1e+2*0.026730739405916; expect_X(3,1) =  1e+2*0.001396884300817; expect_X(3,2) =  1e+2*0.004250513282368; expect_X(3,3) =  1e+2*0.000190925476076; expect_X(3,4) = -1e+2*0.001984142852581; expect_X(3,5) = -1e+2*0.000028773557839; expect_X(3,6) = 1e+2*0; expect_X(3,7) = -1e+2*0.000322384797280; expect_X(3,8) = -1e+2*0.000351053981308;
	expect_X(4,0) = -1e+2*0.313928655928748; expect_X(4,1) = -1e+2*0.163784382994394; expect_X(4,2) = -1e+2*0.044326303597952; expect_X(4,3) = -1e+2*0.001984142852581; expect_X(4,4) =  1e+2*0.527988071584382; expect_X(4,5) =  1e+2*0.000111836392674; expect_X(4,6) = 1e+2*0; expect_X(4,7) =  1e+2*0.001410747036388; expect_X(4,8) =  1e+2*0.002405965955478;
	expect_X(5,0) = -1e+2-0.004690645261773; expect_X(5,1) =  1e+2*0.012779657048731; expect_X(5,2) = -1e+2*0.000662424480368; expect_X(5,3) = -1e+2*0.000028773557839; expect_X(5,4) =  1e+2*0.000111836392674; expect_X(5,5) =  1e+2*0.000041855087295; expect_X(5,6) = 1e+2*0; expect_X(5,7) =  1e+2*0.000476344734388; expect_X(5,8) =  1e+2*0.000583502334040;
	expect_X(6,0) =  1e+2*0.000000000000000; expect_X(6,1) =  1e+2*0.000000000000000; expect_X(6,2) =  1e+2*0.000000000000000; expect_X(6,3) =  1e+2*0.000000000000000; expect_X(6,4) =  1e+2*0.000000000000000; expect_X(6,5) =  1e+2*0.000000000000000; expect_X(6,6) = 1e+2*0; expect_X(6,7) =  1e+2*0.000000000000000; expect_X(6,8) =  1e+2*0.000000000000000;
	expect_X(7,0) = -1e+2*0.052623035300436; expect_X(7,1) =  1e+2*0.143430702460303; expect_X(7,2) = -1e+2*0.007428328837057; expect_X(7,3) = -1e+2*0.000322384797280; expect_X(7,4) =  1e+2*0.001410747036388; expect_X(7,5) =  1e+2*0.000476344734388; expect_X(7,6) = 1e+2*0; expect_X(7,7) =  1e+2*0.005438180536963; expect_X(7,8) =  1e+2*0.006767021988094;
	expect_X(8,0) = -1e+2*0.057714706758673; expect_X(8,1) =  1e+2*0.164521508295555; expect_X(8,2) = -1e+2*0.008137995705764; expect_X(8,3) = -1e+2*0.000351053981308; expect_X(8,4) =  1e+2*0.002405965955478; expect_X(8,5) =  1e+2*0.000583502334040; expect_X(8,6) = 1e+2*0; expect_X(8,7) =  1e+2*0.006767021988094; expect_X(8,8) =  1e+2*0.009074537095327;


	vector_type expect_l(n);
	expect_l(0) = -0.057714706758673;
	expect_l(1) =  0.164521508295555;
	expect_l(2) = -0.008137995705764;
	expect_l(3) = -0.000351053981308;
	expect_l(4) =  0.002405965955478;
	expect_l(5) =  0.000583502334040;
	expect_l(6) =  0.000000000000000;
	expect_l(7) =  0.006767021988094;
	expect_l(8) =  0.009074537095327;


	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.150278082884243; expect_G(0,1) =  0.143143688106951; expect_G(0,2) =  0.018203456199910; expect_G(0,3) =  0.000827154029396; expect_G(0,4) = -0.010056993372684; expect_G(0,5) =  0.000361963320200; expect_G(0,6) = 0; expect_G(0,7) =  0.004381429591722; expect_G(0,8) =  0.007080465330660;
	expect_G(1,0) =  0.593621020214211; expect_G(1,1) = -0.948632309320427; expect_G(1,2) =  0.080390245541794; expect_G(1,3) =  0.003513069380213; expect_G(1,4) = -0.036070591820291; expect_G(1,5) = -0.003806125042661; expect_G(1,6) = 0; expect_G(1,7) = -0.044545703037779; expect_G(1,8) = -0.062095480352244;
	expect_G(2,0) = -4.304428233551964; expect_G(2,1) =  0.015954154527100; expect_G(2,2) = -0.544493658368384; expect_G(2,3) = -0.024361609342909; expect_G(2,4) =  0.278050871512213; expect_G(2,5) =  0.003861907638793; expect_G(2,6) = 0; expect_G(2,7) =  0.042595101037027; expect_G(2,8) =  0.041926823496933;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "Solution X = " << std::scientific << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << std::scientific << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_mathematica_1 )
{
	DCS_DEBUG_TRACE("Test Case: Test Mathematica #1");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n(3);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) = 1; A(0,1) = -1; A(0,2) = 1;
	A(1,0) = 0; A(1,1) =  1; A(1,2) = 1;
	A(2,0) = 0; A(2,1) =  0; A(2,2) = 1;

	matrix_type B(n,m);
	B(0,0) = 1; B(0,1) = 0;
	B(1,0) = 1; B(1,1) = 0;
	B(2,0) = 0; B(2,1) = 1;

	matrix_type Q(n,n);
	Q(0,0) = 10; Q(0,1) = 0; Q(0,2) = 0;
	Q(1,0) =  0; Q(1,1) = 1; Q(1,2) = 0;
	Q(2,0) =  0; Q(2,1) = 0; Q(2,2) = 0.1;

	matrix_type R(m,m);
	R(0,0) = 10; R(0,1) = 0;
	R(1,0) =  0; R(1,1) = 0.1;

	matrix_type S;
	S = ublas::zero_matrix<value_type>(n,m);

	matrix_type E;
	E = ublas::identity_matrix<value_type>(n);

	matrix_type X;

	matrix_type expect_X(n,n);
	expect_X(0,0) = 42.2835365723965; expect_X(0,1) = - 68.5246566154061; expect_X(0,2) = - 3.9478276675312;
	expect_X(1,0) =-68.5246566154061; expect_X(1,1) =  154.0430907758175; expect_X(1,2) =  16.0017331334676;
	expect_X(2,0) =- 3.9478276675312; expect_X(2,1) =   16.0017331334676; expect_X(2,2) =   8.3319663731597;

	vector_type expect_l(n);
	expect_l(0) = complex_value_type(0.014790878209943,+0.065327030440924);
	expect_l(1) = complex_value_type(0.014790878209943,-0.065327030440924);
	expect_l(2) = complex_value_type(0.507913088596327,0);

	matrix_type expect_G(m,n);
	expect_G(0,0) = -0.395758347061307; expect_G(0,1) = 1.599377754500512; expect_G(0,2) = 0.810607779840516;
	expect_G(1,0) =  0.097558030819426; expect_G(1,1) = 0.079555884622979; expect_G(1,2) = 1.258885747544583;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R, S, E);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "S = " << S );
	DCS_DEBUG_TRACE( "E = " << E );
	DCS_DEBUG_TRACE( "Solution X = " << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues l = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_mathematica_2 )
{
	DCS_DEBUG_TRACE("Test Case: Test Mathematica #2");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n(2);
	const std::size_t m(1);

	matrix_type A(n,n);
	A(0,0) = -1.0; A(0,1) =  0.5;
	A(1,0) = -0.5; A(1,1) = -2.0;

	matrix_type B(n,m);
	B(0,0) = 1.0;
	B(1,0) = 0.5;

	matrix_type Q(n,n);
	Q(0,0) = 1; Q(0,1) = 0;
	Q(1,0) = 0; Q(1,1) = 2;

	matrix_type R(m,m);
	R(0,0) = 1;

	matrix_type S;
	S = ublas::zero_matrix<value_type>(n,m);

	matrix_type E;
	E = ublas::identity_matrix<value_type>(n);

	matrix_type X;

	matrix_type expect_X(n,n);
	expect_X(0,0) = 1.925056693914408; expect_X(0,1) =  2.610775872464338;
	expect_X(1,0) = 2.610775872464338; expect_X(1,1) = 31.230349238719803;

	vector_type expect_l(n);
	expect_l(0) = complex_value_type(-0.384777684926709,+0.143417474092307);
	expect_l(1) = complex_value_type(-0.384777684926709,-0.143417474092307);

	matrix_type expect_G(m,n);
	expect_G(0,0) = -0.925056693914411; expect_G(0,1) = -2.610775872464342;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R, S, E);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "S = " << S );
	DCS_DEBUG_TRACE( "E = " << E );
	DCS_DEBUG_TRACE( "Solution X = " << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues l = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


DCS_TEST_DEF( test_mathematica_3 )
{
	DCS_DEBUG_TRACE("Test Case: Test Mathematica #3");

	typedef double value_type;
	typedef ::std::complex<value_type> complex_value_type;
	typedef ublas::matrix<value_type> matrix_type;
	typedef ublas::vector<complex_value_type> vector_type;

	const std::size_t n(3);
	const std::size_t m(2);

	matrix_type A(n,n);
	A(0,0) = 0; A(0,1) =  1.000; A(0,2) =   0.0;
	A(1,0) = 0; A(1,1) = -0.010; A(1,2) =   0.3;
	A(2,0) = 0; A(2,1) = -0.003; A(2,2) = -10.0;

	matrix_type B(n,m);
	B(0,0) = 0.0; B(0,1) =  0;
	B(1,0) = 0.0; B(1,1) = -1;
	B(2,0) = 0.1; B(2,1) =  0;

	matrix_type Q(n,n);
	Q(0,0) = 1; Q(0,1) = 0; Q(0,2) = 0;
	Q(1,0) = 0; Q(1,1) = 1; Q(1,2) = 0;
	Q(2,0) = 0; Q(2,1) = 0; Q(2,2) = 1;

	matrix_type R(m,m);
	R(0,0) = 1; R(0,1) = 0;
	R(1,0) = 0; R(1,1) = 1;

	matrix_type S(n,m);
	S(0,0) =  0.10; S(0,1) = 0.0;
	S(1,0) = -0.20; S(1,1) = 0.0;
	S(2,0) =  0.15; S(2,1) = 0.2;

	matrix_type E;
	E = ublas::identity_matrix<value_type>(n);

	matrix_type X;

	matrix_type expect_X(n,n);
	expect_X(0,0) = 1.0e+3*0.000999899234593; expect_X(0,1) =  1.0e+3*0.000002166454981; expect_X(0,2) =  1.0e+3*0.009899384711326;
	expect_X(1,0) = 1.0e+3*0.000002166454981; expect_X(1,1) =  1.0e+3*0.001981954398444; expect_X(1,2) = -1.0e+3*0.017830763208183;
	expect_X(2,0) = 1.0e+3*0.009899384711326; expect_X(2,1) = -1.0e+3*0.017830763208183; expect_X(2,2) =  1.0e+3*9.930660768134448;


	vector_type expect_l(n);
	expect_l(0) = complex_value_type(-0.004926321073849,+0.012040888479689);
	expect_l(1) = complex_value_type(-0.004926321073849,-0.012040888479689);
	expect_l(2) = complex_value_type(-0.099826902789823, 0.000000000000000);


	matrix_type expect_G(m,n);
	expect_G(0,0) =  0.001007654067580; expect_G(0,1) = -0.021664549812126; expect_G(0,2) = -98.993847113263143;
	expect_G(1,0) = -0.000602532388965; expect_G(1,1) =  0.000935743736164; expect_G(1,2) = - 0.733959482458431;

	dcs_control::dare_solver<value_type> solver;

	solver.solve(A, B, Q, R, S, E);

	DCS_DEBUG_TRACE( "A = " << A );
	DCS_DEBUG_TRACE( "B = " << B );
	DCS_DEBUG_TRACE( "Q = " << Q );
	DCS_DEBUG_TRACE( "R = " << R );
	DCS_DEBUG_TRACE( "S = " << S );
	DCS_DEBUG_TRACE( "E = " << E );
	DCS_DEBUG_TRACE( "Solution X = " << solver.solution() );
	DCS_DEBUG_TRACE( "Gain G = " << solver.gain() );
	DCS_DEBUG_TRACE( "Closed-loop eigenvalues l = " << std::scientific << solver.eigenvalues() );
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.solution(), expect_X, n, n, tol);
	DCS_TEST_CHECK_MATRIX_CLOSE(solver.gain(), expect_G, m, n, tol);
	DCS_TEST_CHECK_VECTOR_CLOSE(solver.eigenvalues(), expect_l, n, tol);
}


int main()
{
	// Tests labeled with 'darex' keyword are taken from:
	//  'DAREX -- A Collection of Benchmark Examples for Discrete-Time
	//   Algebraic Riccati Equations'
	//  (http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.125.6784)
	// Tests labeled with 'mathematica' keyword are taken from the Wolfram
	// Mathematica 8 documentation
	// (http://reference.wolfram.com/mathematica/ref/DiscreteRiccatiSolve.html)
	//
	// All tests has been validated with MATLAB 2009b

	DCS_TEST_BEGIN();

	DCS_TEST_DO( test_1 );
	DCS_TEST_DO( test_2 );
	DCS_TEST_DO( test_3 );
	DCS_TEST_DO( test_darex_ex_1_1 );
	DCS_TEST_DO( test_darex_ex_1_2 );
	DCS_TEST_DO( test_darex_ex_1_3 );
	DCS_TEST_DO( test_darex_ex_1_4 );
	DCS_TEST_DO( test_darex_ex_1_5 );
	DCS_TEST_DO( test_darex_ex_1_6 );
	DCS_TEST_DO( test_darex_ex_1_7 );
	DCS_TEST_DO( test_darex_ex_1_8 );
	DCS_TEST_DO( test_darex_ex_1_9 );
//	DCS_TEST_DO( test_darex_ex_1_10 ); //FIXME: don't work!!!
	DCS_TEST_DO( test_mathematica_1 );
	DCS_TEST_DO( test_mathematica_2 );
	DCS_TEST_DO( test_mathematica_3 );

	DCS_TEST_END();
}
