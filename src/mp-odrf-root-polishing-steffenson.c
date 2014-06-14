/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: root polishing steffenson algorithm
   Date: Sat Mar  7, 2009

   Abstract

	This module implements steffenson root polishing algorithm driver.

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
  mpfr_t	x;
  mpfr_t	x_1;
  mpfr_t	x_2;
  int		count;
} steffenson_state_t;


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: functions.
 ** ----------------------------------------------------------------- */

static void
steffenson_init (void * driver_state)
{
  steffenson_state_t *	state = driver_state;
  mpfr_init(state->f);
  mpfr_init(state->df);
  mpfr_init(state->x);
  mpfr_init(state->x_1);
  mpfr_init(state->x_2);
}
static void
steffenson_final (void * driver_state)
{
  steffenson_state_t *	state = driver_state;
  mpfr_clear(state->f);
  mpfr_clear(state->df);
  mpfr_clear(state->x);
  mpfr_clear(state->x_1);
  mpfr_clear(state->x_2);
}
mp_odrf_code_t
steffenson_set (void * driver_state,
		mp_odrf_mpfr_function_fdf_t * FDF, mpfr_ptr initial_guess)
{
  mp_odrf_code_t	retval	= MP_ODRF_OK;
  steffenson_state_t *	state	= driver_state;
  retval = MP_ODRF_MPFR_FN_FDF_EVAL_F(FDF, state->f, initial_guess);
  if (MP_ODRF_OK == retval) {
    retval = MP_ODRF_MPFR_FN_FDF_EVAL_DF(FDF, state->df, initial_guess);
    if (MP_ODRF_OK == retval) {
      mpfr_set(state->x, initial_guess, GMP_RNDN);
      mpfr_set_si(state->x_1, 0, GMP_RNDN);
      mpfr_set_si(state->x_2, 0, GMP_RNDN);
      state->count = 1;
    }
  }
  return retval;
}
mp_odrf_code_t
steffenson_iterate (void * driver_state,
		    mp_odrf_mpfr_function_fdf_t * FDF, mpfr_ptr root)
{
  mp_odrf_code_t	retval	= MP_ODRF_OK;
  steffenson_state_t *	state	= driver_state;
  if (mpfr_zero_p(state->df)) {
    retval = MP_ODRF_ERROR_DERIVATIVE_IS_ZERO;
  } else {
    mpfr_t	X_new, F_new, DF_new;
    mpfr_t	tmp1, tmp2;
    mpfr_init(X_new);
    mpfr_init(F_new);
    mpfr_init(DF_new);
    mpfr_init(tmp1);
    mpfr_init(tmp2);
    {
      mpfr_div(tmp1, state->f, state->df, GMP_RNDN);
      mpfr_sub(X_new, state->x, tmp1, GMP_RNDN);
      retval = MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(FDF, DF_new, F_new, X_new);
      if (MP_ODRF_OK == retval) {
	mpfr_set(state->x_2, state->x_1, GMP_RNDN);
	mpfr_set(state->x_1, state->x,   GMP_RNDN);
	mpfr_set(state->x,   X_new,      GMP_RNDN);
	mpfr_set(state->f,   F_new,      GMP_RNDN);
	mpfr_set(state->df,  DF_new,     GMP_RNDN);
	if (!mpfr_number_p(F_new)) {
	  retval = MP_ODRF_ERROR_FUNCTION_OR_DERIVATIVE_VALUE_INVALID;
	} else {
	  if (state->count < 3) {
	    mpfr_set(root, X_new, GMP_RNDN);
	    state->count++;
	  } else {
	    mpfr_t	u, v;
	    mpfr_init(u);
	    mpfr_init(v);
	    {
	      mpfr_sub(u, state->x, state->x_1, GMP_RNDN);
	      mpfr_mul_si(tmp1, state->x, 2, GMP_RNDN);
	      mpfr_sub(tmp2, X_new, tmp1, GMP_RNDN);
	      mpfr_add(v, tmp2, state->x_1, GMP_RNDN);
	      if (mpfr_zero_p(v))
		mpfr_set(root, X_new, GMP_RNDN);  /* avoid division by zero */
	      else {
		mpfr_div(tmp1, u, v, GMP_RNDN);
		mpfr_mul(tmp2, u, tmp1, GMP_RNDN);
		mpfr_sub(root, state->x_1, tmp2, GMP_RNDN);  /* accelerated value */
	      }
	    }
	    mpfr_clear(u);
	    mpfr_clear(v);
	  }
	  if (!mpfr_number_p(DF_new)) {
	    retval = MP_ODRF_ERROR_FUNCTION_OR_DERIVATIVE_VALUE_INVALID;
	  }
	}
      }
    }
    mpfr_clear(X_new);
    mpfr_clear(F_new);
    mpfr_clear(DF_new);
    mpfr_clear(tmp1);
    mpfr_clear(tmp2);
  }
  return retval;
}


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: struct definition.
 ** ----------------------------------------------------------------- */

static const mp_odrf_mpfr_root_fdfsolver_driver_t steffenson_driver = {
  .name			= "steffenson",
  .driver_state_size	= sizeof(steffenson_state_t),
  .init			= steffenson_init,
  .final		= steffenson_final,
  .set			= steffenson_set,
  .iterate		= steffenson_iterate
};

const mp_odrf_mpfr_root_fdfsolver_driver_t * \
  mp_odrf_mpfr_root_fdfsolver_steffenson = &steffenson_driver;

/* end of file */
