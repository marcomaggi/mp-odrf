/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: root bracketing falsepos algorithm
   Date: Sat Mar  7, 2009

   Abstract

	This  module  implements   falsepos  root  bracketing  algorithm
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
} falsepos_state_t;


/** --------------------------------------------------------------------
 ** Falsepos root bracketing driver: functions.
 ** ----------------------------------------------------------------- */

static void
falsepos_init (void * driver_state)
/* This is equal to "bisection_init()". */
{
  falsepos_state_t *	state = driver_state;
  mpfr_init(state->y_lower);
  mpfr_init(state->y_upper);
}
static void
falsepos_final (void * driver_state)
/* This is equal to "bisection_final()". */
{
  falsepos_state_t *	state = driver_state;
  mpfr_clear(state->y_lower);
  mpfr_clear(state->y_upper);
}
static mp_odrf_operation_code_t
falsepos_set (void * driver_state, mp_odrf_mpfr_function_t * f,
	       mpfr_ptr root, mpfr_ptr x_lower, mpfr_ptr x_upper,
	       const mp_odrf_error_t ** EE)
/* This is equal to "bisection_set()". */
{
  mp_odrf_operation_code_t	retval = MP_ODRF_OK;
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
    falsepos_state_t *	state = driver_state;
    int			clo, cup;
    SAFE_FUNC_CALL(EE, retval, f, x_lower, state->y_lower);
    if (MP_ODRF_OK == retval) {
      SAFE_FUNC_CALL(EE, retval, f, x_upper, state->y_upper);
      if (MP_ODRF_OK == retval) {
	clo = mpfr_cmp_si(state->y_lower, 0);
	cup = mpfr_cmp_si(state->y_upper, 0);
	if (((clo < 0) && (cup < 0)) ||
	    ((clo > 0) && (cup > 0))) {
	  *EE = &mp_odrf_error_endpoints_do_not_straddle;
	  retval = MP_ODRF_ERROR;
	}
      }
    }
  }
  return retval;
}
static mp_odrf_operation_code_t
falsepos_iterate (void * driver_state, mp_odrf_mpfr_function_t * f,
		  mpfr_t root, mpfr_t x_lower, mpfr_t x_upper,
		  const mp_odrf_error_t ** EE)
{
  mp_odrf_operation_code_t	retval	= MP_ODRF_OK;
  falsepos_state_t *		state	= driver_state;
  if (mpfr_zero_p(state->y_lower)) {
    mpfr_set(root,    x_lower, GMP_RNDN);
    mpfr_set(x_upper, x_lower, GMP_RNDN);
    retval = MP_ODRF_OK;
  } else if (mpfr_zero_p(state->y_upper)) {
    mpfr_set(root,    x_upper, GMP_RNDN);
    mpfr_set(x_lower, x_upper, GMP_RNDN);
    retval = MP_ODRF_OK;
  } else {
    mpfr_t	tmp1, tmp2, tmp3;
    mpfr_t	x_linear, y_linear;
    mpfr_t	x_bisect, y_bisect;
    int		clow, clin, cbis;
    mpfr_init(tmp1);
    mpfr_init(tmp2);
    mpfr_init(tmp3);
    mpfr_init(x_linear);
    mpfr_init(y_linear);
    mpfr_init(x_bisect);
    mpfr_init(y_bisect);
    {
      /* Draw  a line  between f(*lower_bound)  and f(*upper_bound)  and
	 note where  it crosses the X  axis; that's where we  will split
	 the interval. */
      mpfr_sub(tmp1,     x_lower,        x_upper,        GMP_RNDN);
      mpfr_sub(tmp2,     state->y_lower, state->y_upper, GMP_RNDN);
      mpfr_div(tmp3,     tmp1,           tmp2,           GMP_RNDN);
      mpfr_mul(tmp1,     state->y_upper, tmp3,           GMP_RNDN);
      mpfr_sub(x_linear, x_upper,        tmp1,           GMP_RNDN);

      SAFE_FUNC_CALL(EE, retval, f, x_linear, y_linear);
      if (MP_ODRF_OK != retval) {
	goto end;
      }
      if (mpfr_zero_p(y_linear)) {
	mpfr_set(root,    x_linear, GMP_RNDN);
	mpfr_set(x_lower, x_linear, GMP_RNDN);
	mpfr_set(x_upper, x_linear, GMP_RNDN);
	goto end;
      }
      /* Discard  the half  of the  interval which  doesn't contain  the
	 root. */
      mpfr_set(root, x_linear, GMP_RNDN);
      clow = mpfr_cmp_si(state->y_lower, 0);
      clin = mpfr_cmp_si(y_linear,       0);
      if (((clow > 0) && (clin < 0)) ||
	  ((clow < 0) && (clin > 0))) {
	mpfr_set(x_upper,        x_linear, GMP_RNDN);
	mpfr_set(state->y_upper, y_linear, GMP_RNDN);
	mpfr_sub(tmp3, x_linear, x_lower, GMP_RNDN);
      } else {
	mpfr_set(x_lower,        x_linear, GMP_RNDN);
	mpfr_set(state->y_lower, y_linear, GMP_RNDN);
	mpfr_sub(tmp3, x_upper, x_linear, GMP_RNDN);
      }
      mpfr_sub(tmp1, x_upper, x_lower, GMP_RNDN);
      mpfr_mul_d(tmp2, tmp1, 0.5, GMP_RNDN);
      if (mpfr_less_p(tmp3, tmp2))
	goto end;
      mpfr_add(tmp1, x_lower, x_upper, GMP_RNDN);
      mpfr_mul_d(x_bisect, tmp1, 0.5, GMP_RNDN);
      SAFE_FUNC_CALL(EE, retval, f, x_bisect, y_bisect);
      if (MP_ODRF_OK != retval) {
	goto end;
      }
      clow = mpfr_cmp_si(state->y_lower, 0);
      cbis = mpfr_cmp_si(y_bisect,       0);
      if (((clow > 0) && (cbis < 0)) ||
	  ((clow < 0) && (cbis > 0))) {
	mpfr_set(x_upper,        x_bisect, GMP_RNDN);
	mpfr_set(state->y_upper, y_bisect, GMP_RNDN);
	if (mpfr_greater_p(root, x_bisect)) {
	  mpfr_add(tmp1, x_lower, x_bisect, GMP_RNDN);
	  mpfr_mul_d(root, tmp1, 0.5, GMP_RNDN);
	}
      } else {
	mpfr_set(x_lower,        x_bisect, GMP_RNDN);
	mpfr_set(state->y_lower, y_bisect, GMP_RNDN);
	if (mpfr_less_p(root, x_bisect)) {
	  mpfr_add(tmp1, x_bisect, x_upper, GMP_RNDN);
	  mpfr_mul_d(root, tmp1, 0.5, GMP_RNDN);
	}
      }
    }
  end:
    mpfr_clear(tmp1);
    mpfr_clear(tmp2);
    mpfr_clear(tmp3);
    mpfr_clear(x_linear);
    mpfr_clear(y_linear);
    mpfr_clear(x_bisect);
    mpfr_clear(y_bisect);
  }
  return retval;
}


/** --------------------------------------------------------------------
 ** Falsepos root bracketing driver: struct definition.
 ** ----------------------------------------------------------------- */

static const mp_odrf_mpfr_root_fsolver_driver_t falsepos_driver = {
  .name			= "falsepos",
  .driver_state_size	= sizeof(falsepos_state_t),
  .init			= falsepos_init,
  .final		= falsepos_final,
  .set			= falsepos_set,
  .iterate		= falsepos_iterate
};

const mp_odrf_mpfr_root_fsolver_driver_t * \
  mp_odrf_mpfr_root_fsolver_falsepos = &falsepos_driver;

/* end of file */
