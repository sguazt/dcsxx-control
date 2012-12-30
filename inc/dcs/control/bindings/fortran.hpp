/**
 * \file dcs/control/bindings/fortran.hpp
 *
 * \brief Bindings for generic Fortran routines.
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
