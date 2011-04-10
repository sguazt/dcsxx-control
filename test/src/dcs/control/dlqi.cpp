#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <complex>
#include <cstddef>
#include <dcs/control/design/dlqi.hpp>
#include <dcs/test.hpp>
#include <iostream>


namespace ublas = boost::numeric::ublas;
namespace dcs_ctrl = dcs::control;


const double tol = 1.0e-5;


DCS_TEST_DEF( free_func_ill_cond )
{
	DCS_DEBUG_TRACE("Test Case: free function - ill-conditioned problem");

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


	matrix_type expect_K(nx,nx+1); // state-feedback optimal gain
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
	DCS_TEST_CHECK_MATRIX_CLOSE( expect_K, K, nx, nx+1, tol );
}


int main()
{
	// Tests labeled with 'mathematica' keyword are taken from the Wolfram
	// Mathematica 8 documentation
	// (http://reference.wolfram.com/mathematica/ref/LQRegulatorGains.html)
	//
	// All tests has been validated with MATLAB 2009b
 
	DCS_TEST_BEGIN();

	DCS_TEST_DO( free_func_ill_cond );

	DCS_TEST_END();
}
