/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: root polishing newton algorithm
   Date: Sat Mar  7, 2009

   Abstract

	This module implements newton root polishing algorithm driver.

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
} newton_state_t;


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: functions.
 ** ----------------------------------------------------------------- */

static void
newton_init (void * driver_state)
{
  newton_state_t *	state = driver_state;
  mpfr_init(state->f);
  mpfr_init(state->df);
}
static void
newton_final (void * driver_state)
{
  newton_state_t *	state = driver_state;
  mpfr_clear(state->f);
  mpfr_clear(state->df);
}
mp_odrf_operation_code_t
newton_set (void * driver_state,
	    mp_odrf_mpfr_function_fdf_t * FDF,
	    mpfr_ptr initial_guess,
	    const mp_odrf_error_t ** EE)
{
  newton_state_t *		state	= driver_state;
  return MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(FDF, initial_guess, state->f, state->df, EE);
}
mp_odrf_operation_code_t
newton_iterate (void * driver_state,
		mp_odrf_mpfr_function_fdf_t * FDF,
		mpfr_ptr root,
		const mp_odrf_error_t ** EE)
{
  mp_odrf_operation_code_t	retval	= MP_ODRF_OK;
  newton_state_t *		state	= driver_state;
  mpfr_t			tmp1, tmp2;
  mpfr_init(tmp1);
  mpfr_init(tmp2);
  {
    if (mpfr_zero_p(state->df)) {
      retval = MP_ODRF_ERROR;
      *EE    = &mp_odrf_error_derivative_is_zero;
    } else {
      mpfr_div(tmp1, state->f, state->df, GMP_RNDN);
      mpfr_sub(tmp2, root, tmp1, GMP_RNDN);
      mpfr_set(root, tmp2, GMP_RNDN);
      retval = MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(FDF, root, state->f, state->df, EE);
      if (MP_ODRF_OK == retval) {
	if ((!mpfr_number_p(state->f)) || (!mpfr_number_p(state->df))) {
	  retval = MP_ODRF_ERROR;
	  *EE    = &mp_odrf_error_function_or_derivative_value_invalid;
	}
      }
    }
  }
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  return retval;
}


/** --------------------------------------------------------------------
 ** Bisection root bracketing driver: struct definition.
 ** ----------------------------------------------------------------- */

static const mp_odrf_mpfr_root_fdfsolver_driver_t newton_driver = {
  .name			= "newton",
  .driver_state_size	= sizeof(newton_state_t),
  .init			= newton_init,
  .final		= newton_final,
  .set			= newton_set,
  .iterate		= newton_iterate
};

const mp_odrf_mpfr_root_fdfsolver_driver_t * \
  mp_odrf_mpfr_root_fdfsolver_newton = &newton_driver;

/* end of file */
