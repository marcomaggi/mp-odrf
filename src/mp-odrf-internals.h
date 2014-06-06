/*
  Part of: Multiple Precision One Dimensional Root Finding
  Contents: private header file
  Date: Fri Jun  6, 2014

  Abstract

	This header  file is for internal  use.  It must be  included in
	all the C source files.

  Copyright (C) 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

  This program is  free software: you can redistribute  it and/or modify
  it under the  terms of the GNU General Public  License as published by
  the Free Software Foundation, either  version 3 of the License, or (at
  your option) any later version.

  This program  is distributed in the  hope that it will  be useful, but
  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
  General Public License for more details.

  You  should have received  a copy  of the  GNU General  Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MPODRF_INTERNALS_H
#define MPODRF_INTERNALS_H 1


/** --------------------------------------------------------------------
 ** Headers.
 ** ----------------------------------------------------------------- */

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include "mp-odrf.h"


/** --------------------------------------------------------------------
 ** Preprocessor macros.
 ** ----------------------------------------------------------------- */

#define SAFE_FUNC_CALL(EE, RETVAL, F, X, Y)			\
  do {								\
    RETVAL = MP_ODRF_MPFR_FN_EVAL(F, Y, X, EE);			\
    if (MP_ODRF_ERROR == RETVAL) {				\
      ;								\
    } else if (! mpfr_number_p(Y)) {				\
      *EE	= &mp_odrf_function_value_is_not_finite;	\
      RETVAL	= MP_ODRF_ERROR;				\
    }								\
  } while (0);


/** --------------------------------------------------------------------
 ** Constants.
 ** ----------------------------------------------------------------- */

#define GSL_DBL_EPSILON        2.2204460492503131e-16


/** --------------------------------------------------------------------
 ** Private data structures.
 ** ----------------------------------------------------------------- */

mp_odrf_private_decl const mp_odrf_error_t mp_odrf_error_no_memory;
mp_odrf_private_decl const mp_odrf_error_t mp_odrf_error_invalid_bracket_interval;
mp_odrf_private_decl const mp_odrf_error_t mp_odrf_error_relative_tolerance_is_negative;
mp_odrf_private_decl const mp_odrf_error_t mp_odrf_error_absolute_tolerance_is_negative;
mp_odrf_private_decl const mp_odrf_error_t mp_odrf_error_lower_bound_larger_than_upper_bound;
mp_odrf_private_decl const mp_odrf_error_t mp_odrf_error_endpoints_do_not_straddle;
mp_odrf_private_decl const mp_odrf_error_t mp_odrf_function_value_is_not_finite;


/** --------------------------------------------------------------------
 ** Functions.
 ** ----------------------------------------------------------------- */




/** --------------------------------------------------------------------
 ** Done.
 ** ----------------------------------------------------------------- */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* MPODRF_INTERNALS_H */

/* end of file */
