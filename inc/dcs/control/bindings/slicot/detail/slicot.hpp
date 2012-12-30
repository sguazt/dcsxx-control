/**
 * \file dcs/control/bindings/slicot/detail/slicot.hpp
 *
 * \brief Declarations of SLICOT Fortran routines.
 *
 * The subroutine library SLICOT provides Fortran 77 implementations of
 * numerical algorithms for computations in systems and control theory. Based on
 * numerical linear algebra routines from BLAS and LAPACK libraries, SLICOT
 * provides methods for the design and analysis of control systems.
 *
 * The SLICOT library is Copyright 2002-2012 NICONET e.V.
 * For more information, see <a href="http://www.slicot.org>SLICOT web site</a>.
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

#ifndef DCS_CONTROL_BINDINGS_SLICOT_DETAIL_SLICOT_HPP
#define DCS_CONTROL_BINDINGS_SLICOT_DETAIL_SLICOT_HPP


#include <dcs/control/bindings/detail/fortran.hpp>


#define DCS_CONTROL_BINDINGS_SLICOT_SB02OD DCS_CONTROL_BINDINGS_FORTRAN_FUNCTION_NAME(sb02od, SB02OD)
#define DCS_CONTROL_BINDINGS_SLICOT_SG02AD DCS_CONTROL_BINDINGS_FORTRAN_FUNCTION_NAME(sg02ad, SG02AD)


extern "C" {

void DCS_CONTROL_BINDINGS_SLICOT_SB02OD(const dcs_control_bindings_fortran_char* dico,
										const dcs_control_bindings_fortran_char* jobb,
										const dcs_control_bindings_fortran_char* fact,
										const dcs_control_bindings_fortran_char* uplo,
										const dcs_control_bindings_fortran_char* jobl,
										const dcs_control_bindings_fortran_char* sort,
										const ::dcs_control_bindings_fortran_int* n,
										const ::dcs_control_bindings_fortran_int* m,
										const ::dcs_control_bindings_fortran_int* p,
										const dcs_control_bindings_fortran_double_real* a,
										const ::dcs_control_bindings_fortran_int* lda,
										const dcs_control_bindings_fortran_double_real* b,
										const ::dcs_control_bindings_fortran_int* ldb,
										const dcs_control_bindings_fortran_double_real* q,
										const ::dcs_control_bindings_fortran_int* ldq,
										const dcs_control_bindings_fortran_double_real* r,
										const ::dcs_control_bindings_fortran_int* ldr,
										const dcs_control_bindings_fortran_double_real* l,
										const ::dcs_control_bindings_fortran_int* ldl,
										dcs_control_bindings_fortran_double_real* rcond,
										dcs_control_bindings_fortran_double_real* x,
										const ::dcs_control_bindings_fortran_int* ldx,
										dcs_control_bindings_fortran_double_real* alfar,
										dcs_control_bindings_fortran_double_real* alfai,
										dcs_control_bindings_fortran_double_real* beta,
										dcs_control_bindings_fortran_double_real* s,
										const ::dcs_control_bindings_fortran_int* lds,
										dcs_control_bindings_fortran_double_real* t,
										const ::dcs_control_bindings_fortran_int* ldt,
										dcs_control_bindings_fortran_double_real* u,
										const ::dcs_control_bindings_fortran_int* ldu,
										const dcs_control_bindings_fortran_double_real* tol,
										::dcs_control_bindings_fortran_int* iwork,
										dcs_control_bindings_fortran_double_real* dwork,
										const ::dcs_control_bindings_fortran_int* ldwork,
										::dcs_control_bindings_fortran_logical* bwork,
										::dcs_control_bindings_fortran_int* info);

void DCS_CONTROL_BINDINGS_SLICOT_SG02AD(const dcs_control_bindings_fortran_char* dico,
										const dcs_control_bindings_fortran_char* jobb,
										const dcs_control_bindings_fortran_char* fact,
										const dcs_control_bindings_fortran_char* uplo,
										const dcs_control_bindings_fortran_char* jobl,
										const dcs_control_bindings_fortran_char* scal,
										const dcs_control_bindings_fortran_char* sort,
										const dcs_control_bindings_fortran_char* acc,
										const ::dcs_control_bindings_fortran_int* n,
										const ::dcs_control_bindings_fortran_int* m,
										const ::dcs_control_bindings_fortran_int* p,
										const dcs_control_bindings_fortran_double_real* a,
										const ::dcs_control_bindings_fortran_int* lda,
										const dcs_control_bindings_fortran_double_real* e,
										const ::dcs_control_bindings_fortran_int* lde,
										const dcs_control_bindings_fortran_double_real* b,
										const ::dcs_control_bindings_fortran_int* ldb,
										const dcs_control_bindings_fortran_double_real* q,
										const ::dcs_control_bindings_fortran_int* ldq,
										const dcs_control_bindings_fortran_double_real* r,
										const ::dcs_control_bindings_fortran_int* ldr,
										const dcs_control_bindings_fortran_double_real* l,
										const ::dcs_control_bindings_fortran_int* ldl,
										dcs_control_bindings_fortran_double_real* rcond,
										dcs_control_bindings_fortran_double_real* x,
										const ::dcs_control_bindings_fortran_int* ldx,
										dcs_control_bindings_fortran_double_real* alfar,
										dcs_control_bindings_fortran_double_real* alfai,
										dcs_control_bindings_fortran_double_real* beta,
										dcs_control_bindings_fortran_double_real* s,
										const ::dcs_control_bindings_fortran_int* lds,
										dcs_control_bindings_fortran_double_real* t,
										const ::dcs_control_bindings_fortran_int* ldt,
										dcs_control_bindings_fortran_double_real* u,
										const ::dcs_control_bindings_fortran_int* ldu,
										const dcs_control_bindings_fortran_double_real* tol,
										::dcs_control_bindings_fortran_int* iwork,
										dcs_control_bindings_fortran_double_real* dwork,
										const ::dcs_control_bindings_fortran_int* ldwork,
										::dcs_control_bindings_fortran_logical* bwork,
										::dcs_control_bindings_fortran_int* iwarn,
										::dcs_control_bindings_fortran_int* info);
} // extern "C"


#endif // DCS_CONTROL_BINDINGS_SLICOT_DETAIL_SLICOT_HPP
