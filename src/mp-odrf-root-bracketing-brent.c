/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: one dimensional root finding
   Date: Sat Mar  7, 2009

   Abstract

	This module implements brent root bracketing algorithm driver.

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
  mpfr_t	a, b, c, d, e;
  mpfr_t	fa, fb, fc;
} brent_state_t;


/** --------------------------------------------------------------------
 ** Bracketing algorithm: brent.
 ** ----------------------------------------------------------------- */

static void
brent_init (void * driver_state)
{
  brent_state_t *	state = driver_state;
  mpfr_init(state->a);
  mpfr_init(state->b);
  mpfr_init(state->c);
  mpfr_init(state->d);
  mpfr_init(state->e);
  mpfr_init(state->fa);
  mpfr_init(state->fb);
  mpfr_init(state->fc);
}
static void
brent_final (void * driver_state)
{
  brent_state_t *	state = driver_state;
  mpfr_clear(state->a);
  mpfr_clear(state->b);
  mpfr_clear(state->c);
  mpfr_clear(state->d);
  mpfr_clear(state->e);
  mpfr_clear(state->fa);
  mpfr_clear(state->fb);
  mpfr_clear(state->fc);
}
static mp_odrf_code_t
brent_set (void * driver_state, mp_odrf_mpfr_function_t * f,
	   mpfr_ptr root, mpfr_ptr x_lower, mpfr_ptr x_upper)
{
  mp_odrf_code_t	retval = MP_ODRF_OK;
  brent_state_t *	state = driver_state;
  mpfr_t		y_lower, y_upper;
  mpfr_t		tmp1;
  int			clo, cup;
  mpfr_init(y_lower);
  mpfr_init(y_upper);
  mpfr_init(tmp1);
  {
    mpfr_add(tmp1, x_lower, x_upper, GMP_RNDN);
    mpfr_mul_d(root, tmp1, 0.5, GMP_RNDN);
    SAFE_FUNC_CALL(retval, f, x_lower, y_lower);
    if (MP_ODRF_OK == retval) {
      SAFE_FUNC_CALL(retval, f, x_upper, y_upper);
      if (MP_ODRF_OK == retval) {
	mpfr_set(state->a,  x_lower, GMP_RNDN);
	mpfr_set(state->fa, y_lower, GMP_RNDN);
	mpfr_set(state->b,  x_upper, GMP_RNDN);
	mpfr_set(state->fb, y_upper, GMP_RNDN);
	mpfr_set(state->c,  x_upper, GMP_RNDN);
	mpfr_set(state->fc, y_upper, GMP_RNDN);
	mpfr_sub(state->d, x_upper, x_lower, GMP_RNDN);
	mpfr_sub(state->e, x_upper, x_lower, GMP_RNDN);
	clo = mpfr_cmp_si(y_lower, 0);
	cup = mpfr_cmp_si(y_upper, 0);
	if (((clo < 0) && (cup < 0)) ||
	    ((clo > 0) && (cup > 0))) {
	  retval = MP_ODRF_ERROR_ENDPOINTS_DO_NOT_STRADDLE;
	}
      }
    }
  }
  mpfr_clear(y_lower);
  mpfr_clear(y_upper);
  mpfr_clear(tmp1);
  return retval;
}
static mp_odrf_code_t
brent_iterate (void * driver_state, mp_odrf_mpfr_function_t * f,
	       mpfr_t root, mpfr_t x_lower, mpfr_t x_upper)
{
  mp_odrf_code_t	retval	= MP_ODRF_OK;
  brent_state_t *	state = driver_state;
  mpfr_t		tol, m;
  mpfr_t		tmp1, tmp2, tmp3;
  int			ac_equal = 0;
#define A	state->a
#define B	state->b
#define C	state->c
#define D	state->d
#define E	state->e
#define FA	state->fa
#define FB	state->fb
#define FC	state->fc
  mpfr_init(tol);
  mpfr_init(m);
  mpfr_init(tmp1);
  mpfr_init(tmp2);
  mpfr_init(tmp3);
  {
    if (((mpfr_cmp_si(FB, 0) < 0) && (mpfr_cmp_si(FC, 0) < 0)) ||
	((mpfr_cmp_si(FB, 0) > 0) && (mpfr_cmp_si(FC, 0) > 0))) {
      ac_equal = 1;
      mpfr_set(C,    A, GMP_RNDN);
      mpfr_set(FC,  FA, GMP_RNDN);
      mpfr_sub(D, B, A, GMP_RNDN);
      mpfr_sub(E, B, A, GMP_RNDN);
    }
    mpfr_abs(tmp1, FC, GMP_RNDN);
    mpfr_abs(tmp2, FB, GMP_RNDN);
    if (mpfr_less_p(tmp1, tmp2)) {
      ac_equal = 1;
      mpfr_set(A,   B, GMP_RNDN);
      mpfr_set(B,   C, GMP_RNDN);
      mpfr_set(C,   A, GMP_RNDN);
      mpfr_set(FA, FB, GMP_RNDN);
      mpfr_set(FB, FC, GMP_RNDN);
      mpfr_set(FC, FA, GMP_RNDN);
    }
    if (mpfr_zero_p(FB)) {
      mpfr_set(root,    B, GMP_RNDN);
      mpfr_set(x_lower, B, GMP_RNDN);
      mpfr_set(x_upper, B, GMP_RNDN);
      goto end;
    }
    mpfr_abs(tmp1, B, GMP_RNDN);
    mpfr_mul_d(tol, tmp1, 0.5 * GSL_DBL_EPSILON, GMP_RNDN); /* FIXME usage of GSL_DBL_EPSILON */
    mpfr_sub(tmp1, C, B, GMP_RNDN);
    mpfr_mul_d(m, tmp1, 0.5, GMP_RNDN);
    mpfr_abs(tmp1, m, GMP_RNDN);
    if (mpfr_lessequal_p(tmp1, tol)) {
      mpfr_set(root, B, GMP_RNDN);
      if (mpfr_less_p(B, C)) {
	mpfr_set(x_lower, B, GMP_RNDN);
	mpfr_set(x_upper, C, GMP_RNDN);
      } else {
	mpfr_set(x_lower, C, GMP_RNDN);
	mpfr_set(x_upper, B, GMP_RNDN);
      }
      goto end;
    }
    mpfr_abs(tmp1,  E, GMP_RNDN);
    mpfr_abs(tmp2, FA, GMP_RNDN);
    mpfr_abs(tmp3, FB, GMP_RNDN);
    if (mpfr_less_p(tmp1, tol) || mpfr_lessequal_p(tmp2, tmp3)) {
      mpfr_set(D, m, GMP_RNDN);            /* use bisection */
      mpfr_set(E, m, GMP_RNDN);
    } else {
      mpfr_t p, q, r, s;   /* use inverse cubic interpolation */
      mpfr_init(p);
      mpfr_init(q);
      mpfr_init(r);
      mpfr_init(s);
      {
	mpfr_div(s, FB, FA, GMP_RNDN);
	if (ac_equal) {
	  mpfr_mul(tmp1, m, s, GMP_RNDN);
	  mpfr_mul_d(p, tmp1, 0.5, GMP_RNDN);
	  mpfr_sub_si(tmp1, s, 1, GMP_RNDN);
	  mpfr_neg(q, tmp1, GMP_RNDN);
	} else {
	  mpfr_div(q, FA, FC, GMP_RNDN);
	  mpfr_div(r, FB, FC, GMP_RNDN);

	  mpfr_sub(tmp1, B, A, GMP_RNDN);	/* tmp1 = b - a */
	  mpfr_sub_si(tmp2, r, 1, GMP_RNDN);	/* tmp2 = r - 1 */
	  mpfr_mul(tmp3, tmp1, tmp2, GMP_RNDN);	/* tmp3 = (b - a) * (r - 1) */
	  mpfr_sub(tmp1, q, r, GMP_RNDN);	/* tmp1 = q - r */
	  mpfr_mul(tmp2, q, tmp1, GMP_RNDN);	/* tmp2 = q * (q - r) */
	  mpfr_mul(tmp1, tmp2, m, GMP_RNDN);	/* tmp1 = m * q * (q - r) */
	  mpfr_mul_si(tmp2, tmp1, 2, GMP_RNDN);	/* tmp2 = 2 * m * q * (q - r) */
	  mpfr_sub(tmp1, tmp2, tmp3, GMP_RNDN); /* tmp1 = [2 * m * q * (q - r)] -
						   [(b - a) * (r - 1)] */
	  mpfr_mul(p, tmp1, s, GMP_RNDN);	/* p = s * {[2 * m * q * (q - r)] -
						   [(b - a) * (r - 1)]} */

	  mpfr_sub_si(tmp1, q, 1, GMP_RNDN);	/* tmp1 = q - 1 */
	  mpfr_sub_si(tmp2, r, 1, GMP_RNDN);	/* tmp2 = r - 1 */
	  mpfr_mul(tmp3, tmp1, tmp2, GMP_RNDN); /* tmp3 = (q - 1) * (r - 1) */
	  mpfr_sub_si(tmp1, s, 1, GMP_RNDN);	/* tmp1 = (s - 1) */
	  mpfr_mul(q, tmp3, tmp1, GMP_RNDN);	/* q = (q - 1) * (r - 1) * (s - 1) */
	}
	if (mpfr_cmp_si(p, 0) > 0) {
	  mpfr_neg(q, q, GMP_RNDN);
	} else {
	  mpfr_neg(p, p, GMP_RNDN);
	}
	mpfr_mul_si(tmp1, p, 2, GMP_RNDN);
	{
	  mpfr_t	tmp4;
	  mpfr_init(tmp4);
	  {
	    mpfr_mul(tmp2, tol, q, GMP_RNDN);		/* tmp2 = tol * q */
	    mpfr_abs(tmp3, tmp2, GMP_RNDN);		/* tmp3 = fabs(tol * q) */
	    mpfr_mul(tmp2, m, q, GMP_RNDN);		/* tmp2 = m * q */
	    mpfr_mul_si(tmp4, tmp2, 3, GMP_RNDN);	/* tmp4 = 3 * m * q */
	    mpfr_sub(tmp2, tmp4, tmp3, GMP_RNDN);	/* tmp2 = 3*m*1 - fabs(tol*q) */
	    mpfr_mul(tmp3, E, q, GMP_RNDN);		/* tmp3 = e * q */
	    mpfr_abs(tmp4, tmp3, GMP_RNDN);		/* tmp4 = fabs(e * q) */
	    mpfr_min(tmp3, tmp2, tmp4, GMP_RNDN);	/* tmp3 = min(3*m*1 - fabs(tol*q),
							   fabs(e * q)) */
	  }
	  mpfr_clear(tmp4);
	}
	if (mpfr_less_p(tmp1, tmp3)) {
	  mpfr_set(E, D, GMP_RNDN);
	  mpfr_div(D, p, q, GMP_RNDN);
	} else {
	  /* interpolation failed, fall back to bisection */
	  mpfr_set(D, m, GMP_RNDN);
	  mpfr_set(E, m, GMP_RNDN);
	}
      }
      mpfr_clear(p);
      mpfr_clear(q);
      mpfr_clear(r);
      mpfr_clear(s);
    }
    mpfr_set(A,   B, GMP_RNDN);
    mpfr_set(FA, FB, GMP_RNDN);
    mpfr_abs(tmp1, D, GMP_RNDN);
    if (mpfr_greater_p(tmp1, tol)) {
      mpfr_add(tmp1, B, D, GMP_RNDN);
      mpfr_set(B, tmp1, GMP_RNDN);
    } else {
      if (mpfr_cmp_si(m, 0) > 0) {
	mpfr_add(tmp1, B, tol, GMP_RNDN);
	mpfr_set(B, tmp1, GMP_RNDN);
      } else {
	mpfr_sub(tmp1, B, tol, GMP_RNDN);
	mpfr_set(B, tmp1, GMP_RNDN);
      }
    }
    SAFE_FUNC_CALL(retval, f, B, FB);
    if (MP_ODRF_OK != retval) {
      goto end;
    }
    /* Update  the  best  estimate  of  the  root  and  bounds  on  each
       iteration */
    mpfr_set(root, B, GMP_RNDN);

    if (((mpfr_cmp_si(FB, 0) < 0) && (mpfr_cmp_si(FC, 0) < 0)) ||
	((mpfr_cmp_si(FB, 0) > 0) && (mpfr_cmp_si(FC, 0) > 0))) {
      mpfr_set(C, A, GMP_RNDN);
    }
    if (mpfr_less_p(B, C)) {
      mpfr_set(x_lower, B, GMP_RNDN);
      mpfr_set(x_upper, C, GMP_RNDN);
    } else {
      mpfr_set(x_lower, C, GMP_RNDN);
      mpfr_set(x_upper, B, GMP_RNDN);
    }
  }
 end:
  mpfr_clear(tol);
  mpfr_clear(m);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  mpfr_clear(tmp3);
  return retval;
}


/** --------------------------------------------------------------------
 ** Brent root bracketing driver: struct definition.
 ** ----------------------------------------------------------------- */

static const mp_odrf_mpfr_root_fsolver_driver_t brent_driver = {
  .name			= "brent",
  .driver_state_size	= sizeof(brent_state_t),
  .init			= brent_init,
  .final		= brent_final,
  .set			= brent_set,
  .iterate		= brent_iterate
};

const mp_odrf_mpfr_root_fsolver_driver_t * mp_odrf_mpfr_root_fsolver_brent = &brent_driver;

/* end of file */
