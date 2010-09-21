/**
 * \file dcs/control/pid_controller.cpp
 *
 * \brief Test suite for PID controllers.
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

#include <cstddef>
#include <dcs/control/design/pid_controller.hpp>
#include <dcs/math/la/container/dense_matrix.hpp>
#include <dcs/math/la/container/dense_vector.hpp>
#include <dcs/math/la/operation/io.hpp>
#include <dcs/debug.hpp>
#include <dcs/test.hpp>
#include <vector>


DCS_TEST_DEF( test_single_pid_controller )
{
	DCS_DEBUG_TRACE("Test Case: Single-loop PID Controller");

	typedef double real_type;
	typedef ::std::size_t size_type;


	real_type Kp(5);
	real_type Ki(10);
	real_type Kd(0.25);
	real_type ts(1);

	dcs::control::pid_controller<real_type> pid(Kp, Ki, Kd, ts);

	real_type errors[] = {
			-0.08, -0.06, -0.05,  0.04, -0.38, 0.48,  0.34, 0.76,
			-0.07,  0.18, -0.07, -0.31,  0.06, 0.69, -0.32, 0.36
		};
	real_type signals_ok[] = {
			-1.20, -1.70, -2.15, -1.28, -7.31,  2.12,  4.57, 14.41,
			 9.24, 12.56, 10.49,  6.19,  8.79, 18.91, 10.25, 17.67
		};

	size_type n = sizeof(errors)/sizeof(errors[0]);
	for (size_type i = 0; i < n; ++i)
	{
		real_type e = errors[i];
		real_type u_ok = signals_ok[i];
		real_type u = pid.control(e);

		DCS_DEBUG_TRACE("e = " << e << ", u = " << u << " ==> " << u_ok);

		DCS_TEST_CHECK_CLOSE(u, u_ok, 1e-2);
	}
}


DCS_TEST_DEF( test_single_pid_controller_nonunit_step_time )
{
	DCS_DEBUG_TRACE("Test Case: Single-loop PID Controller with non-unit step-time");

	typedef double real_type;
	typedef ::std::size_t size_type;


	real_type Kp(5);
	real_type Ki(10);
	real_type Kd(0.25);
	real_type ts(0.2);

	dcs::control::pid_controller<real_type> pid(Kp, Ki, Kd, ts);

	real_type errors[] = {
			-0.08, -0.06, -0.05,  0.04, -0.38, 0.48,  0.34, 0.76,
			-0.07,  0.18, -0.07, -0.31,  0.06, 0.69, -0.32, 0.36
		};
	real_type signals_ok[] = {
			-0.56, -0.56, -0.62,  0.01, -3.49, 3.38,  2.11, 6.43,
			 0.57,  3.53,  1.52, -0.29,  2.44, 7.30, -0.44, 5.79
		};

	size_type n = sizeof(errors)/sizeof(errors[0]);
	for (size_type i = 0; i < n; ++i)
	{
		real_type e = errors[i];
		real_type u_ok = signals_ok[i];
		real_type u = pid.control(e);

		DCS_DEBUG_TRACE("e = " << e << ", u = " << u << " ==> " << u_ok);

		DCS_TEST_CHECK_CLOSE(u, u_ok, 1e-2);
	}
}


DCS_TEST_DEF( test_multiloop_pid_controller )
{
	DCS_DEBUG_TRACE("Test Case: Multi-loop PID Controller");

	typedef double real_type;
	typedef ::std::size_t size_type;
	typedef ::dcs::math::la::dense_vector<real_type> vector_type;


	vector_type Kp(3);
	Kp(0) = 4; Kp(1) = 3; Kp(2) = 2;
	vector_type Ki(3);
	Ki(0) = 0.4; Ki(1) = 0.3; Ki(2) = 0.2;
	vector_type Kd(3);
	Kd(0) = 0.02; Kd(1) = 0.03; Kd(2) = 0.04;
	real_type ts(1);

	dcs::control::multiloop_pid_controller<vector_type,real_type> pid(Kp, Ki, Kd, ts);

	vector_type v(3);

	::std::vector<vector_type> errors;

	v(0) = -0.08; v(1) =  0.04; v(2) = -0.01;
	errors.push_back(v);
	v(0) = -0.06; v(1) =  0.02; v(2) = -0.02;
	errors.push_back(v);
	v(0) = -0.05; v(1) = -0.02; v(2) = -0.01;
	errors.push_back(v);
	v(0) =  0.04; v(1) = -0.01; v(2) =  0.04;
	errors.push_back(v);
	v(0) = -0.38; v(1) =  0.09; v(2) =  0.02;
	errors.push_back(v);
	v(0) =  0.48; v(1) =  0.12; v(2) = -0.08;
	errors.push_back(v);
	v(0) =  0.34; v(1) =  0.04; v(2) =  0.05;
	errors.push_back(v);
	v(0) =  0.76; v(1) =  0.50; v(2) =  0.24;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.33; v(2) =  0.32;
	errors.push_back(v);
	v(0) =  0.18; v(1) =  0.20; v(2) =  0.55;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.02; v(2) =  0.83;
	errors.push_back(v);
	v(0) = -0.31; v(1) = -0.04; v(2) =  0.60;
	errors.push_back(v);
	v(0) =  0.06; v(1) = -0.09; v(2) =  0.35;
	errors.push_back(v);
	v(0) =  0.69; v(1) = -1.00; v(2) = -0.30;
	errors.push_back(v);
	v(0) = -0.32; v(1) = -0.02; v(2) = -0.05;
	errors.push_back(v);
	v(0) =  0.36; v(1) =  0.40; v(2) = -0.45;
	errors.push_back(v);

	::std::vector<vector_type> signals_ok;

	v(0) = -0.35; v(1) =  0.13; v(2) =	-0.02;
	signals_ok.push_back(v);
	v(0) = -0.30; v(1) =  0.08; v(2) = -0.05;
	signals_ok.push_back(v);
	v(0) = -0.28; v(1) = -0.05; v(2) = -0.03;
	signals_ok.push_back(v);
	v(0) =  0.10; v(1) = -0.02; v(2) =  0.08;
	signals_ok.push_back(v);
	v(0) = -1.74; v(1) =  0.31; v(2) =  0.04;
	signals_ok.push_back(v);
	v(0) =  1.92; v(1) =  0.43; v(2) = -0.18;
	signals_ok.push_back(v);
	v(0) =  1.47; v(1) =  0.20; v(2) =  0.10;
	signals_ok.push_back(v);
	v(0) =  3.47; v(1) =  1.75; v(2) =  0.53;
	signals_ok.push_back(v);
	v(0) =  0.10; v(1) =  1.32; v(2) =  0.75;
	signals_ok.push_back(v);
	v(0) =  1.19; v(1) =  0.99; v(2) =  1.33;
	signals_ok.push_back(v);
	v(0) =  0.15; v(1) =  0.45; v(2) =  2.06;
	signals_ok.push_back(v);
	v(0) = -0.93; v(1) =  0.27; v(2) =  1.70;
	signals_ok.push_back(v);
	v(0) =  0.58; v(1) =  0.09; v(2) =  1.27;
	signals_ok.push_back(v);
	v(0) =  3.38; v(1) = -2.97; v(2) = -0.11;
	signals_ok.push_back(v);
	v(0) = -0.82; v(1) =  0.02; v(2) =  0.42;
	signals_ok.push_back(v);
	v(0) =  2.08; v(1) =  1.39; v(2) = -0.50;
	signals_ok.push_back(v);


	size_type nrows = errors.size();
	for (size_type r = 0; r < nrows; ++r)
	{
		vector_type e = errors[r];
		vector_type u_ok = signals_ok[r];
		vector_type u = pid.control(e);

		DCS_DEBUG_TRACE("e = " << e << ", u = " << u << " ==> " << u_ok);

		size_type ncols = errors[r].size();
		for (size_type c = 0; c < ncols; ++c)
		{
			DCS_TEST_CHECK_CLOSE(u[c], u_ok[c], 1e-2);
		}
	}
}


DCS_TEST_DEF( test_multiloop_pid_controller_nonunit_step_time )
{
	DCS_DEBUG_TRACE("Test Case: Multi-loop PID Controller with non-unit step-time");

	typedef double real_type;
	typedef ::std::size_t size_type;
	typedef ::dcs::math::la::dense_vector<real_type> vector_type;


	vector_type Kp(3);
	Kp(0) = 4; Kp(1) = 3; Kp(2) = 2;
	vector_type Ki(3);
	Ki(0) = 0.4; Ki(1) = 0.3; Ki(2) = 0.2;
	vector_type Kd(3);
	Kd(0) = 0.02; Kd(1) = 0.03; Kd(2) = 0.04;
	real_type ts(0.2);

	dcs::control::multiloop_pid_controller<vector_type,real_type> pid(Kp, Ki, Kd, ts);

	vector_type v(3);

	::std::vector<vector_type> errors;

	v(0) = -0.08; v(1) =  0.04; v(2) = -0.01;
	errors.push_back(v);
	v(0) = -0.06; v(1) =  0.02; v(2) = -0.02;
	errors.push_back(v);
	v(0) = -0.05; v(1) = -0.02; v(2) = -0.01;
	errors.push_back(v);
	v(0) =  0.04; v(1) = -0.01; v(2) =  0.04;
	errors.push_back(v);
	v(0) = -0.38; v(1) =  0.09; v(2) =  0.02;
	errors.push_back(v);
	v(0) =  0.48; v(1) =  0.12; v(2) = -0.08;
	errors.push_back(v);
	v(0) =  0.34; v(1) =  0.04; v(2) =  0.05;
	errors.push_back(v);
	v(0) =  0.76; v(1) =  0.50; v(2) =  0.24;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.33; v(2) =  0.32;
	errors.push_back(v);
	v(0) =  0.18; v(1) =  0.20; v(2) =  0.55;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.02; v(2) =  0.83;
	errors.push_back(v);
	v(0) = -0.31; v(1) = -0.04; v(2) =  0.60;
	errors.push_back(v);
	v(0) =  0.06; v(1) = -0.09; v(2) =  0.35;
	errors.push_back(v);
	v(0) =  0.69; v(1) = -1.00; v(2) = -0.30;
	errors.push_back(v);
	v(0) = -0.32; v(1) = -0.02; v(2) = -0.05;
	errors.push_back(v);
	v(0) =  0.36; v(1) =  0.40; v(2) = -0.45;
	errors.push_back(v);

	::std::vector<vector_type> signals_ok;

	v(0) = -0.33; v(1) = 0.12; v(2) = -0.02;
	signals_ok.push_back(v);
	v(0) = -0.25; v(1) = 0.06; v(2) = -0.04;
	signals_ok.push_back(v);
	v(0) = -0.21; v(1) = -0.06; v(2) = -0.02;
	signals_ok.push_back(v);
	v(0) = 0.16; v(1) = -0.03; v(2) = 0.09;
	signals_ok.push_back(v);
	v(0) = -1.6; v(1) = 0.29; v(2) = 0.04;
	signals_ok.push_back(v);
	v(0) = 2; v(1) = 0.38; v(2) = -0.18;
	signals_ok.push_back(v);
	v(0) = 1.37; v(1) = 0.12; v(2) = 0.13;
	signals_ok.push_back(v);
	v(0) = 3.17; v(1) = 1.62; v(2) = 0.53;
	signals_ok.push_back(v);
	v(0) = -0.28; v(1) = 1.03; v(2) = 0.68;
	signals_ok.push_back(v);
	v(0) = 0.84; v(1) = 0.66; v(2) = 1.19;
	signals_ok.push_back(v);
	v(0) = -0.22; v(1) = 0.11; v(2) = 1.79;
	signals_ok.push_back(v);
	v(0) = -1.2; v(1) = -0.05; v(2) = 1.26;
	signals_ok.push_back(v);
	v(0) = 0.34; v(1) = -0.21; v(2) = 0.77;
	signals_ok.push_back(v);
	v(0) = 2.95; v(1) = -3.12; v(2) = -0.63;
	signals_ok.push_back(v);
	v(0) = -1.28; v(1) = 0.1; v(2) = 0.05;
	signals_ok.push_back(v);
	v(0) = 1.63; v(1) = 1.3; v(2) = -0.9;
	signals_ok.push_back(v);

	size_type nrows = errors.size();
	for (size_type r = 0; r < nrows; ++r)
	{
		vector_type e = errors[r];
		vector_type u_ok = signals_ok[r];
		vector_type u = pid.control(e);

		DCS_DEBUG_TRACE("e = " << e << ", u = " << u << " ==> " << u_ok);

		size_type ncols = errors[r].size();
		for (size_type c = 0; c < ncols; ++c)
		{
			DCS_TEST_CHECK_CLOSE(u[c], u_ok[c], 1e-2);
		}
	}
}


DCS_TEST_DEF( test_mimo_pid_controller )
{
	DCS_DEBUG_TRACE("Test Case: MIMO PID Controller");

	typedef double real_type;
	typedef ::std::size_t size_type;
	typedef ::dcs::math::la::dense_vector<real_type> vector_type;
	typedef ::dcs::math::la::dense_matrix<real_type> matrix_type;


	matrix_type Kp(3,3);
	Kp(0,0) = 4.0; Kp(0,1) = 3.0; Kp(0,2) = 2.0;
	Kp(1,0) = 1.0; Kp(1,1) = 1.3; Kp(1,2) = 0.4;
	Kp(2,0) = 2.2; Kp(2,1) = 1.2; Kp(2,2) = 3.8;
	matrix_type Ki(3,3);
	Ki(0,0) = 0.4; Ki(0,1) = 0.30; Ki(0,2) = 0.2;
	Ki(1,0) = 0.9; Ki(1,1) = 0.63; Ki(1,2) = 1.0;
	Ki(2,0) = 2.2; Ki(2,1) = 0.90; Ki(2,2) = 0.1;
	matrix_type Kd(3,3);
	Kd(0,0) = 0.02; Kd(0,1) = 0.03; Kd(0,2) = 0.04;
	Kd(1,0) = 0.10; Kd(1,1) = 0.05; Kd(1,2) = 0.02;
	Kd(2,0) = 0.40; Kd(2,1) = 0.03; Kd(2,2) = 0.60;

	real_type ts(1);

	dcs::control::mimo_pid_controller<vector_type,matrix_type,real_type> pid(Kp, Ki, Kd, ts);

	vector_type v(3);

	::std::vector<vector_type> errors;

	v(0) = -0.08; v(1) =  0.04; v(2) = -0.01;
	errors.push_back(v);
	v(0) = -0.06; v(1) =  0.02; v(2) = -0.02;
	errors.push_back(v);
	v(0) = -0.05; v(1) = -0.02; v(2) = -0.01;
	errors.push_back(v);
	v(0) =  0.04; v(1) = -0.01; v(2) =  0.04;
	errors.push_back(v);
	v(0) = -0.38; v(1) =  0.09; v(2) =  0.02;
	errors.push_back(v);
	v(0) =  0.48; v(1) =  0.12; v(2) = -0.08;
	errors.push_back(v);
	v(0) =  0.34; v(1) =  0.04; v(2) =  0.05;
	errors.push_back(v);
	v(0) =  0.76; v(1) =  0.50; v(2) =  0.24;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.33; v(2) =  0.32;
	errors.push_back(v);
	v(0) =  0.18; v(1) =  0.20; v(2) =  0.55;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.02; v(2) =  0.83;
	errors.push_back(v);
	v(0) = -0.31; v(1) = -0.04; v(2) =  0.60;
	errors.push_back(v);
	v(0) =  0.06; v(1) = -0.09; v(2) =  0.35;
	errors.push_back(v);
	v(0) =  0.69; v(1) = -1.00; v(2) = -0.30;
	errors.push_back(v);
	v(0) = -0.32; v(1) = -0.02; v(2) = -0.05;
	errors.push_back(v);
	v(0) =  0.36; v(1) =  0.40; v(2) = -0.45;
	errors.push_back(v);

	::std::vector<vector_type> signals_ok;

	v(0) = -0.24; v(1) = -0.09; v(2) = -0.31;
	signals_ok.push_back(v);
	v(0) = -0.26; v(1) = -0.16; v(2) = -0.44;
	signals_ok.push_back(v);
	v(0) = -0.35; v(1) = -0.27; v(2) = -0.55;
	signals_ok.push_back(v);
	v(0) = 0.16; v(1) = -0.06; v(2) = -0.01;
	signals_ok.push_back(v);
	v(0) = -1.39; v(1) = -0.67; v(2) = -1.89;
	signals_ok.push_back(v);
	v(0) = 2.17; v(1) = 0.74; v(2) = 1.28;
	signals_ok.push_back(v);
	v(0) = 1.78; v(1) = 0.82; v(2) = 1.89;
	signals_ok.push_back(v);
	v(0) = 5.75; v(1) = 3.24; v(2) = 6.51;
	signals_ok.push_back(v);
	v(0) = 2.17; v(1) = 2.53; v(2) = 4.38;
	signals_ok.push_back(v);
	v(0) = 3.51; v(1) = 3.65; v(2) = 6.8;
	signals_ok.push_back(v);
	v(0) = 2.66; v(1) = 4.01; v(2) = 6.87;
	signals_ok.push_back(v);
	v(0) = 1.03; v(1) = 3.89; v(2) = 4.44;
	signals_ok.push_back(v);
	v(0) = 1.94; v(1) = 4.5; v(2) = 4.57;
	signals_ok.push_back(v);
	v(0) = 0.31; v(1) = 3.36; v(2) = 2.82;
	signals_ok.push_back(v);
	v(0) = -0.38; v(1) = 3.32; v(2) = 1.93;
	signals_ok.push_back(v);
	v(0) = 2.97; v(1) = 4.64; v(2) = 3.79;
	signals_ok.push_back(v);


	size_type nrows = errors.size();
	for (size_type r = 0; r < nrows; ++r)
	{
		vector_type e = errors[r];
		vector_type u_ok = signals_ok[r];
		vector_type u = pid.control(e);

		DCS_DEBUG_TRACE("e = " << e << ", u = " << u << " ==> " << u_ok);

		size_type ncols = errors[r].size();
		for (size_type c = 0; c < ncols; ++c)
		{
			DCS_TEST_CHECK_CLOSE(u[c], u_ok[c], 1e-2);
		}
	}
}


DCS_TEST_DEF( test_mimo_pid_controller_nonunit_step_time )
{
	DCS_DEBUG_TRACE("Test Case: MIMO PID Controller with non-unit step-time");

	typedef double real_type;
	typedef ::std::size_t size_type;
	typedef ::dcs::math::la::dense_vector<real_type> vector_type;
	typedef ::dcs::math::la::dense_matrix<real_type> matrix_type;


	matrix_type Kp(3,3);
	Kp(0,0) = 4.0; Kp(0,1) = 3.0; Kp(0,2) = 2.0;
	Kp(1,0) = 1.0; Kp(1,1) = 1.3; Kp(1,2) = 0.4;
	Kp(2,0) = 2.2; Kp(2,1) = 1.2; Kp(2,2) = 3.8;
	matrix_type Ki(3,3);
	Ki(0,0) = 0.4; Ki(0,1) = 0.30; Ki(0,2) = 0.2;
	Ki(1,0) = 0.9; Ki(1,1) = 0.63; Ki(1,2) = 1.0;
	Ki(2,0) = 2.2; Ki(2,1) = 0.90; Ki(2,2) = 0.1;
	matrix_type Kd(3,3);
	Kd(0,0) = 0.02; Kd(0,1) = 0.03; Kd(0,2) = 0.04;
	Kd(1,0) = 0.10; Kd(1,1) = 0.05; Kd(1,2) = 0.02;
	Kd(2,0) = 0.40; Kd(2,1) = 0.03; Kd(2,2) = 0.60;

	real_type ts(0.2);

	dcs::control::mimo_pid_controller<vector_type,matrix_type,real_type> pid(Kp, Ki, Kd, ts);

	vector_type v(3);

	::std::vector<vector_type> errors;

	v(0) = -0.08; v(1) =  0.04; v(2) = -0.01;
	errors.push_back(v);
	v(0) = -0.06; v(1) =  0.02; v(2) = -0.02;
	errors.push_back(v);
	v(0) = -0.05; v(1) = -0.02; v(2) = -0.01;
	errors.push_back(v);
	v(0) =  0.04; v(1) = -0.01; v(2) =  0.04;
	errors.push_back(v);
	v(0) = -0.38; v(1) =  0.09; v(2) =  0.02;
	errors.push_back(v);
	v(0) =  0.48; v(1) =  0.12; v(2) = -0.08;
	errors.push_back(v);
	v(0) =  0.34; v(1) =  0.04; v(2) =  0.05;
	errors.push_back(v);
	v(0) =  0.76; v(1) =  0.50; v(2) =  0.24;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.33; v(2) =  0.32;
	errors.push_back(v);
	v(0) =  0.18; v(1) =  0.20; v(2) =  0.55;
	errors.push_back(v);
	v(0) = -0.07; v(1) =  0.02; v(2) =  0.83;
	errors.push_back(v);
	v(0) = -0.31; v(1) = -0.04; v(2) =  0.60;
	errors.push_back(v);
	v(0) =  0.06; v(1) = -0.09; v(2) =  0.35;
	errors.push_back(v);
	v(0) =  0.69; v(1) = -1.00; v(2) = -0.30;
	errors.push_back(v);
	v(0) = -0.32; v(1) = -0.02; v(2) = -0.05;
	errors.push_back(v);
	v(0) =  0.36; v(1) =  0.40; v(2) = -0.45;
	errors.push_back(v);

	::std::vector<vector_type> signals_ok;

	v(0) = -0.22; v(1) = -0.04; v(2) = -0.19;
	signals_ok.push_back(v);
	v(0) = -0.23; v(1) = -0.06; v(2) = -0.23;
	signals_ok.push_back(v);
	v(0) = -0.3; v(1) = -0.12; v(2) = -0.21;
	signals_ok.push_back(v);
	v(0) = 0.22; v(1) = 0.07; v(2) = 0.5;
	signals_ok.push_back(v);
	v(0) = -1.28; v(1) = -0.52; v(2) = -1.75;
	signals_ok.push_back(v);
	v(0) = 2.2; v(1) = 1.04; v(2) = 2.34;
	signals_ok.push_back(v);
	v(0) = 1.62; v(1) = 0.42; v(2) = 1.26;
	signals_ok.push_back(v);
	v(0) = 5.31; v(1) = 2.18; v(2) = 5.27;
	signals_ok.push_back(v);
	v(0) = 1.42; v(1) = 0.46; v(2) = 0.65;
	signals_ok.push_back(v);
	v(0) = 2.69; v(1) = 1.37; v(2) = 4.66;
	signals_ok.push_back(v);
	v(0) = 1.69; v(1) = 0.9; v(2) = 4.09;
	signals_ok.push_back(v);
	v(0) = 0; v(1) = 0.53; v(2) = 1;
	signals_ok.push_back(v);
	v(0) = 0.9; v(1) = 1.11; v(2) = 1.98;
	signals_ok.push_back(v);
	v(0) = -0.81; v(1) = 0.11; v(2) = -0.89;
	signals_ok.push_back(v);
	v(0) = -1.14; v(1) = 0.15; v(2) = -1.43;
	signals_ok.push_back(v);
	v(0) = 2.03; v(1) = 1.88; v(2) = 0.62;
	signals_ok.push_back(v);



	size_type nrows = errors.size();
	for (size_type r = 0; r < nrows; ++r)
	{
		vector_type e = errors[r];
		vector_type u_ok = signals_ok[r];
		vector_type u = pid.control(e);

		DCS_DEBUG_TRACE("e = " << e << ", u = " << u << " ==> " << u_ok);

		size_type ncols = errors[r].size();
		for (size_type c = 0; c < ncols; ++c)
		{
			DCS_TEST_CHECK_CLOSE(u[c], u_ok[c], 1e-2);
		}
	}
}


int main()
{
	DCS_TEST_SUITE("Test Suite for DCS++ Control System Library");

	DCS_TEST_BEGIN();

	DCS_TEST_DO( test_single_pid_controller );

	DCS_TEST_DO( test_single_pid_controller_nonunit_step_time );

	DCS_TEST_DO( test_multiloop_pid_controller );

	DCS_TEST_DO( test_multiloop_pid_controller_nonunit_step_time );

	DCS_TEST_DO( test_mimo_pid_controller );

	DCS_TEST_DO( test_mimo_pid_controller_nonunit_step_time );

	DCS_TEST_END();
}
