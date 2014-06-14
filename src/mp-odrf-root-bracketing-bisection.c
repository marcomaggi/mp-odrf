/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: root bracketing bisection algorithm
   Date: Sat Mar  7, 2009

   Abstract

	This  module  implements  bisection  root  bracketing  algorithm
	driver.

   Copyright (c) 2009, 2014 Marco Maggi <marco.maggi-ipsu@poste.it>
   Copyright (c)  1996, 1997, 1998,  1999, 2000, 2007  Reid Priedhorsky,
   Brian Gough.

   This program is free software:  you can redistribute it and/or modify
   it under the terms of the  GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or (at
   your option) any later version.

   This program is  distributed in the hope that it  will be useful, but
   WITHOUT  ANY   WARRANTY;  without   even  the  implied   warranty  of
   MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE.  See  the GNU
   General Public License for more details.

   You should  have received  a copy of  the GNU General  Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/** --------------------------------------------------------------------
 ** Headers.
 ** ----------------------------------------------------------------- */

#include "mp-odrf-internals.h"

typedef struct {
  mpfr_t	y_lower;
  mpfr_t	y_upper;
} bisection_state_t;


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: functions.
 ** ----------------------------------------------------------------- */

static void
bisection_init (void * driver_state)
{
  bisection_state_t *	state = driver_state;
  mpfr_init(state->y_lower);
  mpfr_init(state->y_upper);
}
static void
bisection_final (void * driver_state)
{
  bisection_state_t *	state = driver_state;
  mpfr_clear(state->y_lower);
  mpfr_clear(state->y_upper);
}
static mp_odrf_code_t
bisection_set (void * driver_state, mp_odrf_mpfr_function_t * f,
	       mpfr_ptr root, mpfr_ptr x_lower, mpfr_ptr x_upper)
{
  mp_odrf_code_t	retval = MP_ODRF_OK;
  {
    mpfr_t	tmp;
    mpfr_init(tmp);
    {
      /* tmp = (x_lower + x_upper) / 2 */
      mpfr_add(tmp, x_lower, x_upper, GMP_RNDN);
      mpfr_mul_d(root, tmp, 0.5, GMP_RNDN);
    }
    mpfr_clear(tmp);
  }
  {
    bisection_state_t *	state = driver_state;
    int			clo, cup;
    SAFE_FUNC_CALL(retval, f, x_lower, state->y_lower);
    if (MP_ODRF_OK == retval) {
      SAFE_FUNC_CALL(retval, f, x_upper, state->y_upper);
      if (MP_ODRF_OK == retval) {
	clo = mpfr_cmp_si(state->y_lower, 0);
	cup = mpfr_cmp_si(state->y_upper, 0);
	if (((clo < 0) && (cup < 0)) ||
	    ((clo > 0) && (cup > 0))) {
	  retval = MP_ODRF_ERROR_ENDPOINTS_DO_NOT_STRADDLE;
	}
      }
    }
  }
  return retval;
}
static mp_odrf_code_t
bisection_iterate (void * driver_state, mp_odrf_mpfr_function_t * f,
		   mpfr_t root, mpfr_t x_lower, mpfr_t x_upper)
{
  mp_odrf_code_t	retval	= MP_ODRF_OK;
  bisection_state_t *	state	= driver_state;
  if (mpfr_zero_p(state->y_lower)) {
    mpfr_set(root,    x_lower, GMP_RNDN);
    mpfr_set(x_upper, x_lower, GMP_RNDN);
  } else if (mpfr_zero_p(state->y_upper)) {
    mpfr_set(root,    x_upper, GMP_RNDN);
    mpfr_set(x_lower, x_upper, GMP_RNDN);
  } else {
    mpfr_t	x_bisect, y_bisect, tmp;
    int		clo, cbi;
    mpfr_init(y_bisect);
    mpfr_init(x_bisect);
    mpfr_init(tmp);
    {
      mpfr_add(tmp, x_lower, x_upper, GMP_RNDN);
      mpfr_mul_d(x_bisect, tmp, 0.5, GMP_RNDN);
      SAFE_FUNC_CALL(retval, f, x_bisect, y_bisect);
      if (MP_ODRF_OK == retval) {
	if (mpfr_zero_p(y_bisect)) {
	  mpfr_set(root,    x_bisect, GMP_RNDN);
	  mpfr_set(x_lower, x_bisect, GMP_RNDN);
	  mpfr_set(x_upper, x_bisect, GMP_RNDN);
	} else {
	  /* Discard the half of the  interval which doesn't contain the
	     root. */
	  clo = mpfr_cmp_si(state->y_lower, 0);
	  cbi = mpfr_cmp_si(y_bisect,       0);
	  if (((clo > 0) && (cbi < 0)) ||
	      ((clo < 0) && (cbi > 0))) {
	    mpfr_add(tmp, x_lower, x_bisect, GMP_RNDN);
	    mpfr_mul_d(root, tmp, 0.5, GMP_RNDN);
	    mpfr_set(x_upper, x_bisect, GMP_RNDN);
	    mpfr_set(state->y_upper, y_bisect, GMP_RNDN);
	  } else {
	    mpfr_add(tmp, x_bisect, x_upper, GMP_RNDN);
	    mpfr_mul_d(root, tmp, 0.5, GMP_RNDN);
	    mpfr_set(x_lower, x_bisect, GMP_RNDN);
	    mpfr_set(state->y_lower, y_bisect, GMP_RNDN);
	  }
	}
      }
    }
    mpfr_clear(tmp);
    mpfr_clear(x_bisect);
    mpfr_clear(y_bisect);
  }
  return retval;
}


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: struct definition.
 ** ----------------------------------------------------------------- */

static const mp_odrf_mpfr_root_fsolver_driver_t bisection_driver = {
  .name			= "bisection",
  .driver_state_size	= sizeof(bisection_state_t),
  .init			= bisection_init,
  .final		= bisection_final,
  .set			= bisection_set,
  .iterate		= bisection_iterate
};

const mp_odrf_mpfr_root_fsolver_driver_t * \
  mp_odrf_mpfr_root_fsolver_bisection = &bisection_driver;

/* end of file */
