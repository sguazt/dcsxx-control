/**
 * \file dcs/control/bindings/detail/fortran.hpp
 *
 * \brief Definitions for interfacing with Fortran routines.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright (C) 2012       Marco Guazzone (marco.guazzone@gmail.com)
 *                          [Distributed Computing System (DCS) Group,
 *                           Computer Science Institute,
 *                           Department of Science and Technological Innovation,
 *                           University of Piemonte Orientale,
 *                           Alessandria (Italy)]
 *
 * This file is part of dcsxx-control (below referred to as "this program").
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DCS_CONTROL_BINDINGS_DETAIL_FORTRAN_HPP
#define DCS_CONTROL_BINDINGS_DETAIL_FORTRAN_HPP


#include <complex>


extern "C" {

// "g77" or "gfortran" or mkl_intel_lp64
//#undef DCS_CONTROL_BINDINGS_DETAIL_BIND_FORTRAN_INTEGER_8
// clapack or "gfortran -fdefault-integer-8" or mkl_intel_ilp64
//#define DCS_CONTROL_BINDINGS_DETAIL_BIND_FORTRAN_INTEGER_8

//#ifndef DCS_CONTROL_BINDINGS_DETAIL_BIND_FORTRAN_INTEGER_8
////From BOOST-Numeric_Bindings
//typedef int fortran_int_t;
//typedef unsigned int fortran_bool_t;
////From LAPACK
////typedef int fortran_int_t;
////typedef fortran_int_t fortran_bool_t;
//#else // DCS_CONTROL_BINDINGS_DETAIL_BIND_FORTRAN_INTEGER_8
//typedef std::ptrdiff_t fortran_int_t;
//typedef std::size_t fortran_bool_t;
//#endif // DCS_CONTROL_BINDINGS_DETAIL_BIND_FORTRAN_INTEGER_8

typedef char dcs_control_bindings_fortran_char;
typedef int dcs_control_bindings_fortran_int;
typedef unsigned int dcs_control_bindings_fortran_logical;
typedef float dcs_control_bindings_fortran_single_real;
typedef double dcs_control_bindings_fortran_double_real;
typedef std::complex<float> dcs_control_bindings_fortran_single_complex;
typedef std::complex<double> dcs_control_bindings_fortran_double_complex;
typedef dcs_control_bindings_fortran_logical (*dcs_control_bindings_fortran_function)(...);

} // extern "C"

#define DCS_CONTROL_BINDINGS_FORTRAN_FUNCTION_NAME(id, ID) id##_

#endif // DCS_CONTROL_BINDINGS_DETAIL_FORTRAN_HPP
