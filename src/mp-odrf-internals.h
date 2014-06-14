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

#define SAFE_FUNC_CALL(RETVAL, F, X, Y)				\
  do {								\
    RETVAL = MP_ODRF_MPFR_FN_EVAL(F, Y, X);			\
    if (MP_ODRF_ERROR == RETVAL) {				\
      ;								\
    } else if (! mpfr_number_p(Y)) {				\
      RETVAL	= MP_ODRF_ERROR_FUNCTION_VALUE_IS_NOT_FINITE;	\
    }								\
  } while (0);


/** --------------------------------------------------------------------
 ** Constants.
 ** ----------------------------------------------------------------- */

#define GSL_DBL_EPSILON        2.2204460492503131e-16


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
