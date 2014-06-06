/*
  Part of: Multiple Precision One-Dimensional Root-Finding
  Contents: static definitions of error descriptors
  Date: Fri Jun  6, 2014

  Abstract

	This file defines  the error descriptor structs used  by the API
	functions to describe exceptions.

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

const mp_odrf_error_t	mp_odrf_error_no_memory = {
  .code		= MP_ODRF_NO_MEMORY_FOR_STATE_STRUCT,
  .description	= "failed to allocate space for root solver state"
};
const mp_odrf_error_t	mp_odrf_error_invalid_bracket_interval = {
  .code		= MP_ODRF_INVALID_BRACKET_INTERVAL,
  .description	= "invalid bracket interval (lower > upper)"
};
const mp_odrf_error_t	mp_odrf_error_relative_tolerance_is_negative = {
  .code		= MP_ODRF_RELATIVE_TOLERANCE_IS_NEGATIVE,
  .description	= "relative tolerance is negative"
};
const mp_odrf_error_t	mp_odrf_error_absolute_tolerance_is_negative = {
  .code		= MP_ODRF_ABSOLUTE_TOLERANCE_IS_NEGATIVE,
  .description	= "absolute tolerance is negative"
};
const mp_odrf_error_t	mp_odrf_error_lower_bound_larger_than_upper_bound = {
  .code		= MP_ODRF_LOWER_BOUND_LARGER_THAN_UPPER_BOUND,
  .description	= "lower bound larger than upper bound"
};
const mp_odrf_error_t	mp_odrf_error_endpoints_do_not_straddle = {
  .code		= MP_ODRF_ERROR_ENDPOINTS_DO_NOT_STRADDLE,
  .description	= "endpoints do not straddle y=0"
};
const mp_odrf_error_t	mp_odrf_error_function_value_is_not_finite = {
  .code		= MP_ODRF_ERROR_FUNCTION_VALUE_IS_NOT_FINITE,
  .description	= "function value is not finite"
};
const mp_odrf_error_t	mp_odrf_error_derivative_is_zero = {
  .code		= MP_ODRF_ERROR_DERIVATIVE_IS_ZERO,
  .description	= "derivative is zero"
};
const mp_odrf_error_t	mp_odrf_error_function_or_derivative_value_invalid = {
  .code		= MP_ODRF_ERROR_FUNCTION_OR_DERIVATIVE_VALUE_INVALID,
  .description	= "function or derivative value is not finite or not a number"
};

/* end of file */
