/**
 * \file test/src/dcs/control/stabilizability.cpp
 *
 * \brief Test suite for system stabilizability.
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
 * \author Marco Guazzone, &lt;marco.guazzone@mfn.unipmn.it&gt;
 */

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublasx/operation/rank.hpp>
#include <dcs/debug.hpp>
#include <dcs/test.hpp>
#include <dcs/control/analysis/stabilizability.hpp>


namespace ublas = boost::numeric::ublas;
namespace ublasx = boost::numeric::ublasx;


const double tol = 1.0e-5;


DCS_TEST_DEF( stabilizability_check )
{
	DCS_DEBUG_TRACE("Test Case: Stabilizability Check");

	typedef double value_type;
	typedef ublas::matrix<value_type> matrix_type;

	// Stabilizable and controllable system
	{
		const std::size_t n = 2;
		const std::size_t m = 1;

		matrix_type A(n,n);
		A(0,0) =  1; A(0,1) = 1;
		A(1,0) = -1; A(1,1) = 2;

		matrix_type B(n,m);
		B(0,0) = 1;
		B(1,0) = 2;


		bool stbl = dcs::control::is_stabilizable(A, B, true);
		bool expect_stbl = true;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("Is Stabilizable = " << std::boolalpha << stbl);
		DCS_TEST_CHECK(stbl == expect_stbl);
	}

	// Stabilizable but not controllable system
	{
		const std::size_t n = 3;
		const std::size_t m = 2;

		matrix_type A(n,n);
		A(0,0) = -1; A(0,1) = 3; A(0,2) = 0;
		A(1,0) =  0; A(1,1) = 1; A(1,2) = 1;
		A(1,0) =  0; A(1,1) = 2; A(1,2) = 0;

		matrix_type B(n,m);
		B(0,0) = 1; B(0,1) = 2.0/3.0;
		B(1,0) = 0; B(1,1) = 1.0/3.0;
		B(2,0) = 0; B(2,1) = 1.0/3.0;


		bool stbl = dcs::control::is_stabilizable(A, B, true);
		bool expect_stbl = true;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("Is Stabilizable = " << std::boolalpha << stbl);
		DCS_TEST_CHECK(stbl == expect_stbl);
	}

	// Non stabilizable system
	{
		const std::size_t n = 2;
		const std::size_t m = 2;

		matrix_type A(n,n);
		A(0,0) = 1; A(0,1) =  1;
		A(1,0) = 4; A(1,1) = -2;

		matrix_type B(n,m);
		B(0,0) = 1; B(0,1) = -1;
		B(1,0) = 1; B(1,1) = -1;


		bool stbl = dcs::control::is_stabilizable(A, B, true);
		bool expect_stbl = false;
		DCS_DEBUG_TRACE("A = " << A);
		DCS_DEBUG_TRACE("B = " << B);
		DCS_DEBUG_TRACE("Is Stabilizable = " << std::boolalpha << stbl);
		DCS_TEST_CHECK(stbl == expect_stbl);
	}
}


int main()
{
	DCS_TEST_BEGIN();

	DCS_TEST_DO( stabilizability_check );

	DCS_TEST_END();
}
