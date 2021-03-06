/**
 * \file test/src/dcs/control/controllability.cpp
 *
 * \brief Test suite for system controllability.
 *
 * Copyright (C) 2009-2010  Distributed Computing System (DCS) Group, Computer
 * Science Department - University of Piemonte Orientale, Alessandria (Italy).
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
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 */

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublasx/operation/rank.hpp>
#include <dcs/debug.hpp>
#include <dcs/test.hpp>
#include <dcs/control/analysis/controllability.hpp>


namespace ublas = boost::numeric::ublas;
namespace ublasx = boost::numeric::ublasx;


const double tol = 1.0e-5;


DCS_TEST_DEF( state_controllability_matrix )
{
	DCS_TEST_CASE("State Controllability Matrix");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	const std::size_t n = 2;
	const std::size_t m = 1;

	matrix_type A(n,n);
	A(0,0) =  1; A(0,1) = 1;
	A(1,0) = -1; A(1,1) = 2;

	matrix_type B(n,m);
	B(0,0) = 1;
	B(1,0) = 2;


	matrix_type expect_C(n, n*m);
	expect_C(0,0) = 1; expect_C(0,1) = 3;
	expect_C(1,0) = 2; expect_C(1,1) = 3;

	matrix_type C = dcs::control::make_controllability_matrix(A, B);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("P = [B AB ... A^{" << (n-1) << "}B] = " << C);
	DCS_DEBUG_TRACE("Number of uncontrollable states " << (n-ublasx::rank(C)));
	DCS_TEST_CHECK_MATRIX_CLOSE(C, expect_C, n, n*m, tol);
}


DCS_TEST_DEF( state_controllability_check )
{
	DCS_TEST_CASE("State Controllability Check");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	// State controllable system
	{
		const std::size_t n = 2;
		const std::size_t m = 1;

		matrix_type A(n,n);
		A(0,0) =  1; A(0,1) = 1;
		A(1,0) = -1; A(1,1) = 2;

		matrix_type B(n,m);
		B(0,0) = 1;
		B(1,0) = 2;


		bool ctrb = dcs::control::is_controllable(A, B);
		bool expect_ctrb = true;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("Is State Controllable = " << std::boolalpha << ctrb);
		DCS_TEST_CHECK(ctrb == expect_ctrb);
	}

	// Non state controllable system
	{
		const std::size_t n = 2;
		const std::size_t m = 2;

		matrix_type A(n,n);
		A(0,0) = 1; A(0,1) =  1;
		A(1,0) = 4; A(1,1) = -2;

		matrix_type B(n,m);
		B(0,0) = 1; B(0,1) = -1;
		B(1,0) = 1; B(1,1) = -1;


		bool ctrb = dcs::control::is_controllable(A, B);
		bool expect_ctrb = false;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("Is State Controllable = " << std::boolalpha << ctrb);
		DCS_TEST_CHECK(ctrb == expect_ctrb);
	}
}


DCS_TEST_DEF( output_controllability_matrix )
{
	DCS_TEST_CASE("Output Controllability Matrix");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	const std::size_t n = 2;
	const std::size_t m = 1;
	const std::size_t p = 1;

	matrix_type A(n,n);
	A(0,0) =  1; A(0,1) = 1;
	A(1,0) = -1; A(1,1) = 2;

	matrix_type B(n,m);
	B(0,0) = 1;
	B(1,0) = 2;

	matrix_type C(p,n);
	C(0,0) = 1; C(0,1) = -1;

	matrix_type D(p,m);
	D(0,0) = 1;


	matrix_type expect_P(p, (n+1)*m);
	expect_P(0,0) = -1; expect_P(0,1) = 0; expect_P(0,2) = 1;

	matrix_type P = dcs::control::make_output_controllability_matrix(A, B, C, D);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("B = " << B);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("D = " << D);
	DCS_DEBUG_TRACE("P = [CB CAB ... CA^{" << (n-1) << "}B D] = " << P);
	DCS_TEST_CHECK_MATRIX_CLOSE(P, expect_P, m, (n+1)*m, tol);
}


DCS_TEST_DEF( output_controllability_check )
{
	DCS_TEST_CASE("Output Controllability Check");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	// Output controllable system
	{
		const std::size_t n = 2;
		const std::size_t m = 1;
		const std::size_t p = 2;

		matrix_type A(n,n);
		A(0,0) =  1; A(0,1) = 1;
		A(1,0) = -1; A(1,1) = 2;

		matrix_type B(n,m);
		B(0,0) = 1;
		B(1,0) = 2;

		matrix_type C(p,n);
		C(0,0) = 1; C(0,1) = -1;
		C(1,0) = 2; C(1,1) = -3;

		matrix_type D(p,m);
		D(0,0) = 1;
		D(1,0) = 2;


		bool ctrb = dcs::control::is_output_controllable(A, B, C, D);
		bool expect_ctrb = true;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("C = " << C);
		DCS_DEBUG_TRACE("D = " << D);
		DCS_DEBUG_TRACE("Is Output Controllable = " << std::boolalpha << ctrb);
		DCS_TEST_CHECK(ctrb == expect_ctrb);
	}

	// Non output controllable system
	{
		const std::size_t n = 2;
		const std::size_t m = 1;
		const std::size_t p = 2;

		matrix_type A(n,n);
		A(0,0) =  0; A(0,1) =  1;
		A(1,0) = -1; A(1,1) = -2;

		matrix_type B(n,m);
		B(0,0) =  1;
		B(1,0) = -1;

		matrix_type C(p,n);
		C(0,0) = 1; C(0,1) = 0;
		C(1,0) = 1; C(1,1) = 1;

		matrix_type D(p,m);
		D(0,0) = 0;
		D(1,0) = 0;


		bool ctrb = dcs::control::is_output_controllable(A, B, C, D);
		bool expect_ctrb = false;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("C = " << C);
		DCS_DEBUG_TRACE("D = " << D);
		DCS_DEBUG_TRACE("Is Output Controllable = " << std::boolalpha << ctrb);
		DCS_TEST_CHECK(ctrb == expect_ctrb);
	}
}


DCS_TEST_DEF( controllable_decomposition )
{
	DCS_TEST_CASE("Controllable Decomposition");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	// Output controllable system
	{
//		const std::size_t n = 2;
//		const std::size_t m = 2;
//		const std::size_t p = 2;
		const std::size_t n = 3;
		const std::size_t m = 2;
		const std::size_t p = 2;

		matrix_type A(n,n);
//		A(0,0) = 1; A(0,1) =  1;
//		A(1,0) = 4; A(1,1) = -2;
		A(0,0) = 1; A(0,1) = 0; A(0,2) = 0;
		A(1,0) = 0; A(1,1) = 2; A(1,2) = 0;
		A(2,0) = 0; A(2,1) = 0; A(2,2) = 3;

		matrix_type B(n,m);
//		B(0,0) = 1; B(0,1) = -1;
//		B(1,0) = 1; B(1,1) = -1;
		B(0,0) = 0; B(0,1) = 0;
		B(1,0) = 4; B(1,1) = 5;
		B(2,0) = 0; B(2,1) = 0;

		matrix_type C(p,n);
//		C(0,0) = 1; C(0,1) = 0;
//		C(1,0) = 0; C(1,1) = 1;
		C(0,0) = 6; C(0,1) = 0; C(0,2) = 0;
		C(1,0) = 0; C(1,1) = 0; C(1,2) = 7;


		dcs::control::controllable_decomposition<value_type> ctrbf;

		ctrbf.decompose(A, B, 1.7764e-15);
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("C = " << C);
		DCS_DEBUG_TRACE("Ab = " << ctrbf.A_bar());
		DCS_DEBUG_TRACE("Bb = " << ctrbf.B_bar());
		DCS_DEBUG_TRACE("Cb = " << ctrbf.C_bar(C));
		DCS_DEBUG_TRACE("T = " << ctrbf.T());
		DCS_DEBUG_TRACE("k = " << ctrbf.k());
//		DCS_DEBUG_TRACE("Is Output Controllable = " << std::boolalpha << ctrb);
//		DCS_TEST_CHECK(ctrb == expect_ctrb);
	}
}


int main()
{
	DCS_TEST_SUITE("DCS Control :: Controllability");

	DCS_TEST_BEGIN();

	DCS_TEST_DO( state_controllability_matrix );
	DCS_TEST_DO( state_controllability_check );
	DCS_TEST_DO( output_controllability_matrix );
	DCS_TEST_DO( output_controllability_check );
	DCS_TEST_DO( controllable_decomposition );

	DCS_TEST_END();
}
