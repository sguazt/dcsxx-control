/**
 * \file test/src/dcs/control/observability.cpp
 *
 * \brief Test suite for system observability.
 *
 * Copyright (C) 2009-2011  Distributed Computing System (DCS) Group, Computer
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
#include <dcs/control/analysis/observability.hpp>


namespace ublas = boost::numeric::ublas;
namespace ublasx = boost::numeric::ublasx;


const double tol = 1.0e-5;


DCS_TEST_DEF( observability_matrix )
{
	DCS_TEST_CASE("Observability Matrix");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	const std::size_t n = 3;
	const std::size_t m = 1;

	matrix_type A(n,n);
	A(0,0) = -1; A(0,1) =  2; A(0,2) = -1;
	A(1,0) = -1; A(1,1) = -2; A(1,2) =  0;
	A(2,0) = -3; A(2,1) = -3; A(2,2) =  2;

	matrix_type C(m,n);
	C(0,0) = 1; C(0,1) = 0; C(0,2) = 0;


	matrix_type expect_O(n*m, n);
	expect_O(0,0) =  1; expect_O(0,1) =  0; expect_O(0,2) =  0;
	expect_O(1,0) = -1; expect_O(1,1) =  2; expect_O(1,2) = -1;
	expect_O(2,0) =  2; expect_O(2,1) = -3; expect_O(2,2) = -1;

	matrix_type O = dcs::control::make_observability_matrix(A, C);
	DCS_DEBUG_TRACE("A = " << A);
	DCS_DEBUG_TRACE("C = " << C);
	DCS_DEBUG_TRACE("O = [C CA ... CA^{" << (n-1) << "}] = " << O);
	DCS_DEBUG_TRACE("Number of unobservable states: " << (n-ublasx::rank(O)) );
	DCS_TEST_CHECK_MATRIX_CLOSE(O, expect_O, n*m, n, tol);
}


DCS_TEST_DEF( observability_check )
{
	DCS_TEST_CASE("Observability Check");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	// observable system
	{
		const std::size_t n = 3;
		const std::size_t m = 1;

		matrix_type A(n,n);
		A(0,0) =  1; A(0,1) = -3; A(0,2) = -1;
		A(1,0) = -1; A(1,1) = -2; A(1,2) =  2;
		A(1,0) = -1; A(1,1) =  2; A(1,2) = -2;

		matrix_type C(m,n);
		C(0,0) = 1; C(0,1) = 1; C(0,2) = 0;


		bool obsv = dcs::control::is_observable(A, C);
		bool expect_obsv = true;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("C = " << C);
		DCS_DEBUG_TRACE("Is Observable = " << std::boolalpha << obsv);
		DCS_TEST_CHECK(obsv == expect_obsv);
	}

	// Non observable system
	{
		const std::size_t n = 3;
		const std::size_t m = 2;

		matrix_type A(n,n);
		A(0,0) = -3; A(0,1) =  1; A(0,2) =  0;
		A(1,0) =  0; A(1,1) = -3; A(1,2) =  1;
		A(2,0) =  0; A(2,1) =  0; A(2,2) = -3;

		matrix_type C(m,n);
		C(0,0) = 0; C(0,1) = 1; C(0,2) =  1;
		C(1,0) = 0; C(1,1) = 2; C(1,2) = -1;


		bool obsv = dcs::control::is_observable(A, C);
		bool expect_obsv = false;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("C = " << C);
		DCS_DEBUG_TRACE("Is Observable = " << std::boolalpha << obsv);
		DCS_TEST_CHECK(obsv == expect_obsv);
	}
}


int main()
{
	DCS_TEST_SUITE("DCS Control :: Observability");

	DCS_TEST_BEGIN();

	DCS_TEST_DO( observability_matrix );
	DCS_TEST_DO( observability_check );

	DCS_TEST_END();
}
