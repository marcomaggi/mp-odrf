/*
  Part of: Multiple Precision One-Dimensional Root-Finding
  Contents: definitions for error handling
  Date: Fri Jun  6, 2014

  Abstract

	This file  defines the  error messages  associated to  the error
	codes.

  Copyright (C) 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

  This program is  free software: you can redistribute  it and/or modify
  it under the  terms of the GNU General Public  License as published by
  the Free Software Foundation, either version  3 of the License, or (at
  your option) any later version.

  This program  is distributed in the  hope that it will  be useful, but
  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See the  GNU
  General Public License for more details.

  You should  have received  a copy  of the  GNU General  Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/** --------------------------------------------------------------------
 ** Headers.
 ** ----------------------------------------------------------------- */

#include "mp-odrf-internals.h"


/** --------------------------------------------------------------------
 ** Definitions.
 ** ----------------------------------------------------------------- */

const char *
mp_odrf_strerror (mp_odrf_code_t code)
{
  switch (code) {
  case MP_ODRF_OK:
    return "no error";
  case MP_ODRF_ERROR_NO_MEMORY_FOR_STATE_STRUCT:
    return "failed to allocate space for root solver state";
  case MP_ODRF_ERROR_INVALID_BRACKET_INTERVAL:
    return "invalid bracket interval (lower > upper)";
  case MP_ODRF_ERROR_RELATIVE_TOLERANCE_IS_NEGATIVE:
    return "relative tolerance is negative";
  case MP_ODRF_ERROR_ABSOLUTE_TOLERANCE_IS_NEGATIVE:
    return "absolute tolerance is negative";
  case MP_ODRF_ERROR_LOWER_BOUND_LARGER_THAN_UPPER_BOUND:
    return "lower bound larger than upper bound";
  case MP_ODRF_ERROR_ENDPOINTS_DO_NOT_STRADDLE:
    return "endpoints do not straddle y=0";
  case MP_ODRF_ERROR_FUNCTION_VALUE_IS_NOT_FINITE:
    return "function value is not finite";
  case MP_ODRF_ERROR_DERIVATIVE_IS_ZERO:
    return "derivative is zero";
  case MP_ODRF_ERROR_FUNCTION_OR_DERIVATIVE_VALUE_INVALID:
    return "function or derivative value is not finite or not a number";
  default:
    return "unknown or invalid error code";
  }
}

/* end of file */
