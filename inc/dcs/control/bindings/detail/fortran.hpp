/**
 * \file dcs/control/bindings/detail/fortran.hpp
 *
 * \brief Definitions for interfacing with Fortran routines.
 *
 * \author Marco Guazzone (marco.guazzone@gmail.com)
 *
 * <hr/>
 *
 * Copyright 2012 Marco Guazzone (marco.guazzone@gmail.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
