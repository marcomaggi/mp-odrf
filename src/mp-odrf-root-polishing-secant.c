/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: root polishing secant algorithm
   Date: Sat Mar  7, 2009

   Abstract

	This module implements secant root polishing algorithm driver.

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
  mpfr_t	f, df;
} secant_state_t;


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: functions.
 ** ----------------------------------------------------------------- */

static void
secant_init (void * driver_state)
{
  secant_state_t *	state = driver_state;
  mpfr_init(state->f);
  mpfr_init(state->df);
}
static void
secant_final (void * driver_state)
{
  secant_state_t *	state = driver_state;
  mpfr_clear(state->f);
  mpfr_clear(state->df);
}
mp_odrf_operation_code_t
secant_set (void * driver_state,
	    mp_odrf_mpfr_function_fdf_t * FDF,
	    mpfr_ptr initial_guess,
	    const mp_odrf_error_t ** EE)
{
  secant_state_t *		state	= driver_state;
  return MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(FDF, initial_guess, state->f, state->df, EE);
}
mp_odrf_operation_code_t
secant_iterate (void * driver_state,
		mp_odrf_mpfr_function_fdf_t * FDF,
		mpfr_ptr root,
		const mp_odrf_error_t ** EE)
{
  mp_odrf_operation_code_t	retval	= MP_ODRF_OK;
  secant_state_t *		state	= driver_state;
  if (mpfr_zero_p(state->df)) {
    retval = MP_ODRF_ERROR;
    *EE    = &mp_odrf_error_derivative_is_zero;
  } else {
    mpfr_t	F_new, DF_new, X_new;
    mpfr_t	deltaF, deltaX;
    mpfr_init(deltaF);
    mpfr_init(deltaX);
    mpfr_init(X_new);
    mpfr_init(F_new);
    mpfr_init(DF_new);
    {
      mpfr_div(deltaF, state->f, state->df, GMP_RNDN);
      mpfr_sub(X_new, root, deltaF, GMP_RNDN);
      /* F_new = F(X_new) */
      retval = MP_ODRF_MPFR_FN_FDF_EVAL_F(FDF, F_new, X_new, EE);
      if (MP_ODRF_OK == retval) {
	/* Compute the incremental ratio of F. */
	{
	  /* deltaF = F_new - state->f */
	  mpfr_sub(deltaF, F_new, state->f, GMP_RNDN);
	  /* deltaX = X_new - root */
	  mpfr_sub(deltaX, X_new, root, GMP_RNDN);
	  /* DF_new = deltaF / deltaX */
	  mpfr_div(DF_new, deltaF, deltaX, GMP_RNDN);
	}
	mpfr_set(root,       X_new, GMP_RNDN);
	mpfr_set(state->f,   F_new, GMP_RNDN);
	mpfr_set(state->df, DF_new, GMP_RNDN);
	if ((!mpfr_number_p(F_new)) || (!mpfr_number_p(DF_new))) {
	  retval = MP_ODRF_ERROR;
	  *EE    = &mp_odrf_error_function_or_derivative_value_invalid;
	}
      }
    }
    mpfr_clear(deltaF);
    mpfr_clear(deltaX);
    mpfr_clear(X_new);
    mpfr_clear(F_new);
    mpfr_clear(DF_new);
  }
  return retval;
}


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: struct definition.
 ** ----------------------------------------------------------------- */

static const mp_odrf_mpfr_root_fdfsolver_driver_t secant_driver = {
  .name			= "secant",
  .driver_state_size	= sizeof(secant_state_t),
  .init			= secant_init,
  .final		= secant_final,
  .set			= secant_set,
  .iterate		= secant_iterate
};

const mp_odrf_mpfr_root_fdfsolver_driver_t * \
  mp_odrf_mpfr_root_fdfsolver_secant = &secant_driver;

/* end of file */
