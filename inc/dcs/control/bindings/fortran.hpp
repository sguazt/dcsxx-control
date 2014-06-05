/**
 * \file dcs/control/bindings/fortran.hpp
 *
 * \brief Bindings for generic Fortran routines.
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

#ifndef DCS_CONTROL_BINDINGS_FORTRAN_HPP
#define DCS_CONTROL_BINDINGS_FORTRAN_HPP


#include <dcs/control/bindings/detail/fortran.hpp>


namespace dcs { namespace control { namespace bindings {

/// Binding for the Fortran CHARACTER type
typedef dcs_control_bindings_fortran_char fortran_char;
/// Binding for the Fortran COMPLEX double precision type
typedef dcs_control_bindings_fortran_double_complex fortran_double_complex;
/// Binding for the Fortran COMPLEX single precision type
typedef dcs_control_bindings_fortran_single_complex fortran_single_complex;
/// Binding for the Fortran REAL double precision type
typedef dcs_control_bindings_fortran_double_real fortran_double_real;
/// Binding for the Fortran REAL single precision type
typedef dcs_control_bindings_fortran_single_real fortran_single_real;
/// Binding for the Fortran INTEGER type
typedef dcs_control_bindings_fortran_int fortran_int;
/// Binding for the Fortran LOGICAL type
typedef dcs_control_bindings_fortran_logical fortran_logical;
/// Binding for the Fortran external function type
typedef dcs_control_bindings_fortran_function fortran_function;

}}} //Namespace dcs::control::bindings

#endif // DCS_CONTROL_BINDINGS_FORTRAN_HPP
