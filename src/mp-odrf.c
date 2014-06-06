/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: one dimensional root finding
   Date: Sat Mar  7, 2009

   Abstract

	This module implements  one dimensional root-finding algorithms.
	For the original  GSL code look in the GSL  source tree, "roots"
	directory.

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

#define DEBUGGING		0
#define MP_ODRF_MPFR_USE_INLINE	1

#include "mp-odrf-internals.h"


/** --------------------------------------------------------------------
 ** Bracketing algorithm: falsepos.
 ** ----------------------------------------------------------------- */

typedef bisection_state_t	falsepos_state_t;

#define falsepos_init	bisection_init
#define falsepos_final	bisection_final
#define falsepos_set	bisection_set

static int
falsepos_iterate (void * vstate, mp_odrf_mpfr_function * f,
		  mpfr_t root, mpfr_t x_lower, mpfr_t x_upper)
{
  /*
    falsepos_state_t * state = (falsepos_state_t *) vstate;
    double x_linear, f_linear;
    double x_bisect, f_bisect;
    double x_left  = *x_lower;
    double x_right = *x_upper;
    double f_lower = state->f_lower;
    double f_upper = state->f_upper;
    double w;
    if (f_lower == 0.0) {
      *root = x_left ;
      *x_upper = x_left;
      return GSL_SUCCESS;
    }
    if (f_upper == 0.0) {
      *root = x_right ;
      *x_lower = x_right;
      return GSL_SUCCESS;
    }
    x_linear = x_right - (f_upper * (x_left - x_right) / (f_lower - f_upper));
    SAFE_FUNC_CALL (f, x_linear, &f_linear);
    if (f_linear == 0.0) {
      *root = x_linear;
      *x_lower = x_linear;
      *x_upper = x_linear;
      return GSL_SUCCESS;
    }
    if ((f_lower > 0.0 && f_linear < 0.0) ||
        (f_lower < 0.0 && f_linear > 0.0)) {
      *root = x_linear ;
      *x_upper = x_linear;
      state->f_upper = f_linear;
      w = x_linear - x_left ;
    } else {
      *root = x_linear;
      *x_lower = x_linear;
      state->f_lower = f_linear;
      w = x_right - x_linear;
    }
    if (w < 0.5 * (x_right - x_left))
      return GSL_SUCCESS;
    x_bisect = 0.5 * (x_left + x_right);
    SAFE_FUNC_CALL (f, x_bisect, &f_bisect);
    if ((f_lower > 0.0 && f_bisect < 0.0) ||
        (f_lower < 0.0 && f_bisect > 0.0)) {
      *x_upper = x_bisect;
      state->f_upper = f_bisect;
      if (*root > x_bisect)
        *root = 0.5 * (x_left + x_bisect) ;
    } else {
      *x_lower = x_bisect;
      state->f_lower = f_bisect;
      if (*root < x_bisect)
        *root = 0.5 * (x_bisect + x_right) ;
    }
    return GSL_SUCCESS;
  */
  falsepos_state_t * state = vstate;
  if (mpfr_zero_p(state->y_lower)) {
    mpfr_set(root,    x_lower, GMP_RNDN);
    mpfr_set(x_upper, x_lower, GMP_RNDN);
    return GSL_SUCCESS;
  } else if (mpfr_zero_p(state->y_upper)) {
    mpfr_set(root,    x_upper, GMP_RNDN);
    mpfr_set(x_lower, x_upper, GMP_RNDN);
    return GSL_SUCCESS;
  } else {
    mpfr_t	tmp1, tmp2, tmp3;
    mpfr_t	x_linear, y_linear;
    mpfr_t	x_bisect, y_bisect;
    int		e, retval = GSL_SUCCESS;
    int		clow, clin, cbis;
    mpfr_init(tmp1);
    mpfr_init(tmp2);
    mpfr_init(tmp3);
    mpfr_init(x_linear);
    mpfr_init(y_linear);
    mpfr_init(x_bisect);
    mpfr_init(y_bisect);
    {
      /* Draw  a line  between f(*lower_bound)  and  f(*upper_bound) and
	 note where  it crosses the X  axis; that's where  we will split
	 the interval. */
      mpfr_sub(tmp1,     x_lower,        x_upper,        GMP_RNDN);
      mpfr_sub(tmp2,     state->y_lower, state->y_upper, GMP_RNDN);
      mpfr_div(tmp3,     tmp1,           tmp2,           GMP_RNDN);
      mpfr_mul(tmp1,     state->y_upper, tmp3,           GMP_RNDN);
      mpfr_sub(x_linear, x_upper,        tmp1,           GMP_RNDN);

      SAFE_FUNC_CALL(e, f, x_linear, y_linear);
      if (GSL_SUCCESS != e) {
	retval = e;
	goto end;
      }
      if (mpfr_zero_p(y_linear)) {
	mpfr_set(root,    x_linear, GMP_RNDN);
	mpfr_set(x_lower, x_linear, GMP_RNDN);
	mpfr_set(x_upper, x_linear, GMP_RNDN);
	goto end;
      }
      /* Discard  the half  of the  interval which  doesn't  contain the
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
      if (mpfr_less_p(tmp3, tmp2)) goto end;
      mpfr_add(tmp1, x_lower, x_upper, GMP_RNDN);
      mpfr_mul_d(x_bisect, tmp1, 0.5, GMP_RNDN);
      SAFE_FUNC_CALL(e, f, x_bisect, y_bisect);
      if (GSL_SUCCESS != e) {
	retval = e;
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
    return retval;
  }
}

static const mp_odrf_mpfr_root_fsolver_type falsepos_type = {
  .name		= "falsepos",
  .size		= sizeof(falsepos_state_t),
  .init		= falsepos_init,
  .final	= falsepos_final,
  .set		= falsepos_set,
  .iterate	= falsepos_iterate
};

const mp_odrf_mpfr_root_fsolver_type * mp_odrf_mpfr_root_fsolver_falsepos = &falsepos_type;


/** --------------------------------------------------------------------
 ** Bracketing algorithm: brent.
 ** ----------------------------------------------------------------- */

typedef struct {
  mpfr_t a, b, c, d, e;
  mpfr_t fa, fb, fc;
} brent_state_t;

static void
brent_init (void * vstate)
{
  brent_state_t *	state = vstate;
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
brent_final (void * vstate)
{
  brent_state_t *	state = vstate;
  mpfr_clear(state->a);
  mpfr_clear(state->b);
  mpfr_clear(state->c);
  mpfr_clear(state->d);
  mpfr_clear(state->e);
  mpfr_clear(state->fa);
  mpfr_clear(state->fb);
  mpfr_clear(state->fc);
}
static int
brent_set (void * vstate, mp_odrf_mpfr_function * f,
	   mpfr_t root, mpfr_t x_lower, mpfr_t x_upper)
{
  /*
    brent_state_t * state = (brent_state_t *) vstate;
    double f_lower, f_upper;
    *root = 0.5 * (x_lower + x_upper);
    SAFE_FUNC_CALL (f, x_lower, &f_lower);
    SAFE_FUNC_CALL (f, x_upper, &f_upper);
    state->a  = x_lower;
    state->fa = f_lower;
    state->b  = x_upper;
    state->fb = f_upper;
    state->c  = x_upper;
    state->fc = f_upper;
    state->d  = x_upper - x_lower;
    state->e  = x_upper - x_lower;
    if ((f_lower < 0.0 && f_upper < 0.0) ||
        (f_lower > 0.0 && f_upper > 0.0))
      GSL_ERROR ("endpoints do not straddle y=0", GSL_EINVAL);
    return GSL_SUCCESS;
  */
  brent_state_t * state = vstate;
  mpfr_t	y_lower, y_upper;
  mpfr_t	tmp1;
  int		e, clo, cup, retval = GSL_SUCCESS;
  mpfr_init(y_lower);
  mpfr_init(y_upper);
  mpfr_init(tmp1);
  {
    mpfr_add(tmp1, x_lower, x_upper, GMP_RNDN);
    mpfr_mul_d(root, tmp1, 0.5, GMP_RNDN);
    SAFE_FUNC_CALL(e, f, x_lower, y_lower);
    if (GSL_SUCCESS != e) {
      retval = e;
      goto end;
    }
    SAFE_FUNC_CALL(e, f, x_upper, y_upper);
    if (GSL_SUCCESS != e) {
      retval = e;
      goto end;
    }
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
	retval = GSL_EINVAL;
	goto end;
    }
  }
  mpfr_clear(y_lower);
  mpfr_clear(y_upper);
  mpfr_clear(tmp1);
 end:
  if (GSL_SUCCESS == retval)
    return retval;
  else
    GSL_ERROR("endpoints do not straddle y=0", retval);
}
static int
brent_iterate (void * vstate, mp_odrf_mpfr_function * f,
	       mpfr_t root, mpfr_t x_lower, mpfr_t x_upper)
{
  /*
    brent_state_t * state = (brent_state_t *) vstate;
    double tol, m;
    int ac_equal = 0;
    double a = state->a, b = state->b, c = state->c;
    double fa = state->fa, fb = state->fb, fc = state->fc;
    double d = state->d, e = state->e;
    if ((fb < 0 && fc < 0) ||
        (fb > 0 && fc > 0)) {
      ac_equal = 1;
      c = a;
      fc = fa;
      d = b - a;
      e = b - a;
    }
    if (fabs (fc) < fabs (fb)) {
      ac_equal = 1;
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol = 0.5 * GSL_DBL_EPSILON * fabs (b);
    m = 0.5 * (c - b);
    if (fb == 0) {
      *root = b;
      *x_lower = b;
      *x_upper = b;
      return GSL_SUCCESS;
    }
    if (fabs (m) <= tol) {
      *root = b;
      if (b < c) {
        *x_lower = b;
        *x_upper = c;
      } else {
        *x_lower = c;
        *x_upper = b;
      }
      return GSL_SUCCESS;
    }
    if (fabs (e) < tol ||
        fabs (fa) <= fabs (fb)) {
      d = m;
      e = m;
    } else {
      double p, q, r;
      double s = fb / fa;
      if (ac_equal) {
        p = 2 * m * s;
        q = 1 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2 * m * q * (q - r) - (b - a) * (r - 1));
        q = (q - 1) * (r - 1) * (s - 1);
      }
      if (p > 0) {
        q = -q;
      } else {
        p = -p;
      }
      if (2 * p < GSL_MIN (3 * m * q - fabs (tol * q), fabs (e * q))) {
        e = d;
        d = p / q;
      } else {
        d = m;
        e = m;
      }
    }
    a = b;
    fa = fb;
    if (fabs (d) > tol) {
      b += d;
    } else {
      b += (m > 0 ? +tol : -tol);
    }
    SAFE_FUNC_CALL (f, b, &fb);
    state->a = a;
    state->b = b;
    state->c = c;
    state->d = d;
    state->e = e;
    state->fa = fa;
    state->fb = fb;
    state->fc = fc;
    *root = b;
    if ((fb < 0 && fc < 0) ||
        (fb > 0 && fc > 0)) {
      c = a;
    }
    if (b < c) {
      *x_lower = b;
      *x_upper = c;
    } else {
      *x_lower = c;
      *x_upper = b;
    }
    return GSL_SUCCESS;
  */
  brent_state_t *	state = vstate;
  mpfr_t	tol, m;
  mpfr_t	tmp1, tmp2, tmp3;
  int		e, ac_equal = 0, retval = GSL_SUCCESS;
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
    SAFE_FUNC_CALL(e, f, B, FB);
    if (GSL_SUCCESS != e) {
      retval = e;
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

static const mp_odrf_mpfr_root_fsolver_type brent_type = {
  .name		= "brent",
  .size		= sizeof(brent_state_t),
  .init		= brent_init,
  .final	= brent_final,
  .set		= brent_set,
  .iterate	= brent_iterate
};

const mp_odrf_mpfr_root_fsolver_type * mp_odrf_mpfr_root_fsolver_brent = &brent_type;


/** --------------------------------------------------------------------
 ** Root polishing algorithm: newton.
 ** ----------------------------------------------------------------- */

typedef struct {
  mpfr_t	f, df;
} newton_state_t;

static void
newton_init (void * vstate)
{
  newton_state_t *	state = vstate;
  mpfr_init(state->f);
  mpfr_init(state->df);
}
static void
newton_final (void * vstate)
{
  newton_state_t *	state = vstate;
  mpfr_clear(state->f);
  mpfr_clear(state->df);
}
static int
newton_set (void * vstate, mp_odrf_mpfr_function_fdf * fdf, mpfr_t initial_guess)
{
  /*
    newton_state_t * state = (newton_state_t *) vstate;
    const double x = *root;
    state->f  = GSL_FN_FDF_EVAL_F  (fdf, x);
    state->df = GSL_FN_FDF_EVAL_DF (fdf, x);
    return GSL_SUCCESS;
  */
  newton_state_t * state = vstate;
  return MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(fdf, initial_guess, state->f, state->df);
}
static int
newton_iterate (void * vstate, mp_odrf_mpfr_function_fdf * fdf, mpfr_t root)
{
  /*
    newton_state_t * state = (newton_state_t *) vstate;
    double root_new, f_new, df_new;
    if (state->df == 0.0)
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    root_new = *root - (state->f / state->df);
    *root = root_new ;
    GSL_FN_FDF_EVAL_F_DF(fdf, root_new, &f_new, &df_new);
    state->f = f_new;
    state->df = df_new;
    if (!gsl_finite(f_new))
      GSL_ERROR ("function value is not finite", GSL_EBADFUNC);
    if (!gsl_finite (df_new))
      GSL_ERROR ("derivative value is not finite", GSL_EBADFUNC);
    return GSL_SUCCESS;
  */
  newton_state_t * state = vstate;
  mpfr_t	tmp1, tmp2;
  int		e, retval = GSL_SUCCESS;
  mpfr_init(tmp1);
  mpfr_init(tmp2);
  {
    if (mpfr_zero_p(state->df)) {
      retval = GSL_EZERODIV;
      goto end;
    }
    mpfr_div(tmp1, state->f, state->df, GMP_RNDN);
    mpfr_sub(tmp2, root, tmp1, GMP_RNDN);
    mpfr_set(root, tmp2, GMP_RNDN);
    e = GSL_FN_FDF_EVAL_F_DF(fdf, root, state->f, state->df);
    if (GSL_SUCCESS != e) {
      retval = e;
      goto end;
    }
    if ((!mpfr_number_p(state->f)) ||
	(!mpfr_number_p(state->df))) {
      retval = GSL_EBADFUNC;
      goto end;
    }
  }
 end:
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  switch (retval) {
  case GSL_SUCCESS:
    return retval;
  case GSL_EZERODIV:
    GSL_ERROR("derivative is zero", retval);
  case GSL_EBADFUNC:
    GSL_ERROR("function or derivative value is not finite or not a number", retval);
  default:
    GSL_ERROR("error evaluating the function", retval);
  }
}

static const mp_odrf_mpfr_root_fdfsolver_type newton_type = {
  .name		= "newton",
  .size		= sizeof(newton_state_t),
  .init		= newton_init,
  .final	= newton_final,
  .set		= newton_set,
  .iterate	= newton_iterate
};

const mp_odrf_mpfr_root_fdfsolver_type * mp_odrf_mpfr_root_fdfsolver_newton = &newton_type;


/** --------------------------------------------------------------------
 ** Root polishing algorithm: secant.
 ** ----------------------------------------------------------------- */

typedef newton_state_t	secant_state_t;
#define secant_init	newton_init
#define secant_final	newton_final
#define secant_set	newton_set

static int
secant_iterate (void * vstate, mp_odrf_mpfr_function_fdf * fdf, mpfr_t root)
{
  /*
    secant_state_t * state = (secant_state_t *) vstate;
    const double x = *root ;
    const double f = state->f;
    const double df = state->df;
    double x_new, f_new, df_new;
    if (state->df == 0.0)
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    x_new = x - (f / df);
    f_new = GSL_FN_FDF_EVAL_F(fdf, x_new);
    df_new = (f_new - f) / (x_new - x);
    *root = x_new;
    state->f  = f_new;
    state->df = df_new;
    if (!gsl_finite (f_new))
      GSL_ERROR ("function value is not finite", GSL_EBADFUNC);
    if (!gsl_finite (df_new))
      GSL_ERROR ("derivative value is not finite", GSL_EBADFUNC);
    return GSL_SUCCESS;
  */
  secant_state_t * state = vstate;
  if (mpfr_zero_p(state->df))
    GSL_ERROR("derivative is zero", GSL_EZERODIV);
  mpfr_t	f_new, df_new;
  mpfr_t	tmp1, tmp2, x_new;
  int		e, retval = GSL_SUCCESS;
  mpfr_init(tmp1);
  mpfr_init(tmp2);
  mpfr_init(x_new);
  mpfr_init(f_new);
  mpfr_init(df_new);
  {
    mpfr_div(tmp1, state->f, state->df, GMP_RNDN);
    mpfr_sub(x_new, root, tmp1, GMP_RNDN);
    e = MP_ODRF_MPFR_FN_FDF_EVAL_F(fdf, f_new, x_new);
    if (GSL_SUCCESS != e) {
      retval = e;
      goto end;
    }
    mpfr_sub(tmp1, f_new, state->f, GMP_RNDN);
    mpfr_sub(tmp2, x_new, root, GMP_RNDN);
    mpfr_div(df_new, tmp1, tmp2, GMP_RNDN);
    mpfr_set(root,       x_new, GMP_RNDN);
    mpfr_set(state->f,   f_new, GMP_RNDN);
    mpfr_set(state->df, df_new, GMP_RNDN);

    if ((!mpfr_number_p(f_new)) ||
	(!mpfr_number_p(df_new))) {
      retval = GSL_EBADFUNC;
      goto end;
    }
  }
 end:
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  mpfr_clear(x_new);
  mpfr_clear(f_new);
  mpfr_clear(df_new);
  switch (retval) {
  case GSL_EBADFUNC:
    GSL_ERROR("function or derivative value is not finite", retval);
    break;
  case GSL_SUCCESS:
    return retval;
  default:
    GSL_ERROR("error computing function", retval);
  }
}

static const mp_odrf_mpfr_root_fdfsolver_type secant_type = {
  .name		= "secant",
  .size		= sizeof(secant_state_t),
  .init		= secant_init,
  .final	= secant_final,
  .set		= secant_set,
  .iterate	= secant_iterate
};

const mp_odrf_mpfr_root_fdfsolver_type * mp_odrf_mpfr_root_fdfsolver_secant = &secant_type;


/** --------------------------------------------------------------------
 ** Root polishing algorithm: steffenson.
 ** ----------------------------------------------------------------- */

typedef struct {
  mpfr_t f, df;
  mpfr_t x;
  mpfr_t x_1;
  mpfr_t x_2;
  int count;
} steffenson_state_t;

static void
steffenson_init (void * vstate)
{
  steffenson_state_t *	state = vstate;
  mpfr_init(state->f);
  mpfr_init(state->df);
  mpfr_init(state->x);
  mpfr_init(state->x_1);
  mpfr_init(state->x_2);
}
static void
steffenson_final (void * vstate)
{
  steffenson_state_t *	state = vstate;
  mpfr_clear(state->f);
  mpfr_clear(state->df);
  mpfr_clear(state->x);
  mpfr_clear(state->x_1);
  mpfr_clear(state->x_2);
}
static int
steffenson_set (void * vstate, mp_odrf_mpfr_function_fdf * fdf, mpfr_t initial_guess)
{
  /*
    steffenson_state_t * state = (steffenson_state_t *) vstate;
    const double x = *root ;
    state->f = GSL_FN_FDF_EVAL_F (fdf, x);
    state->df = GSL_FN_FDF_EVAL_DF (fdf, x) ;
    state->x = x;
    state->x_1 = 0.0;
    state->x_2 = 0.0;
    state->count = 1;
    return GSL_SUCCESS;
  */
  steffenson_state_t * state = vstate;
  int		e;
  e = MP_ODRF_MPFR_FN_FDF_EVAL_F(fdf, state->f, initial_guess);
  if (GSL_SUCCESS != e) return e;
  e = MP_ODRF_MPFR_FN_FDF_EVAL_DF (fdf, state->df, initial_guess);
  if (GSL_SUCCESS != e) return e;
  mpfr_set(state->x, initial_guess, GMP_RNDN);
  mpfr_set_si(state->x_1, 0, GMP_RNDN);
  mpfr_set_si(state->x_2, 0, GMP_RNDN);
  state->count = 1;
  return GSL_SUCCESS;
}
static int
steffenson_iterate (void * vstate, mp_odrf_mpfr_function_fdf * fdf, mpfr_t guess)
{
  /*
    steffenson_state_t * state = (steffenson_state_t *) vstate;
    double x_new, f_new, df_new;
    double x_1 = state->x_1;
    double x   = state->x;
    if (state->df == 0.0)
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    x_new = x - (state->f / state->df);
    GSL_FN_FDF_EVAL_F_DF(fdf, x_new, &f_new, &df_new);
    state->x_2 = x_1;
    state->x_1 = x;
    state->x = x_new;
    state->f = f_new;
    state->df = df_new;
    if (!gsl_finite (f_new))
      GSL_ERROR ("function value is not finite", GSL_EBADFUNC);
    if (state->count < 3) {
      *root = x_new ;
      state->count++ ;
    } else {
      double u = (x - x_1) ;
      double v = (x_new - 2 * x + x_1);
      if (v == 0)
        *root = x_new;
      else
        *root = x_1 - u * u / v ;
    }
    if (!gsl_finite (df_new))
      GSL_ERROR ("derivative value is not finite", GSL_EBADFUNC);
    return GSL_SUCCESS;
  */
  steffenson_state_t * state = (steffenson_state_t *) vstate;
  if (mpfr_zero_p(state->df))
    GSL_ERROR("derivative is zero", GSL_EZERODIV);
  mpfr_t	x_new, f_new, df_new;
  mpfr_t	tmp1, tmp2;
  int		e, retval = GSL_SUCCESS;
  mpfr_init(x_new);
  mpfr_init(f_new);
  mpfr_init(df_new);
  mpfr_init(tmp1);
  mpfr_init(tmp2);
  {
    mpfr_div(tmp1, state->f, state->df, GMP_RNDN);
    mpfr_sub(x_new, state->x, tmp1, GMP_RNDN);
    e = MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(fdf, x_new, f_new, df_new);
    if (GSL_SUCCESS != e) {
      retval = e;
      goto end;
    }
    mpfr_set(state->x_2, state->x_1, GMP_RNDN);
    mpfr_set(state->x_1, state->x,   GMP_RNDN);
    mpfr_set(state->x,   x_new,      GMP_RNDN);
    mpfr_set(state->f,   f_new,      GMP_RNDN);
    mpfr_set(state->df,  df_new,     GMP_RNDN);
    if (!mpfr_number_p(f_new)) {
      retval = GSL_EBADFUNC;
      goto end;
    }
    if (state->count < 3) {
      mpfr_set(guess, x_new, GMP_RNDN);
      state->count++;
    } else {
      mpfr_t	u, v;
      mpfr_init(u);
      mpfr_init(v);
      {
	mpfr_sub(u, state->x, state->x_1, GMP_RNDN);
	mpfr_mul_si(tmp1, state->x, 2, GMP_RNDN);
	mpfr_sub(tmp2, x_new, tmp1, GMP_RNDN);
	mpfr_add(v, tmp2, state->x_1, GMP_RNDN);
	if (mpfr_zero_p(v))
	  mpfr_set(guess, x_new, GMP_RNDN);  /* avoid division by zero */
	else {
	  mpfr_div(tmp1, u, v, GMP_RNDN);
	  mpfr_mul(tmp2, u, tmp1, GMP_RNDN);
	  mpfr_sub(guess, state->x_1, tmp2, GMP_RNDN);  /* accelerated value */
	}
      }
      mpfr_clear(u);
      mpfr_clear(v);
    }
    if (!mpfr_number_p(df_new)) {
      retval = GSL_EBADFUNC;
      goto end;
    }
  }
 end:
  mpfr_clear(x_new);
  mpfr_clear(f_new);
  mpfr_clear(df_new);
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  switch (retval) {
  case GSL_EBADFUNC:
    GSL_ERROR("function or derivative value is not finite", retval);
    break;
  case GSL_SUCCESS:
    return retval;
  default:
    GSL_ERROR("error computing function", retval);
  }
}

static const mp_odrf_mpfr_root_fdfsolver_type steffenson_type = {
  .name		= "steffenson",
  .size		= sizeof(steffenson_state_t),
  .init		= steffenson_init,
  .final	= steffenson_final,
  .set		= steffenson_set,
  .iterate	= steffenson_iterate
};

const mp_odrf_mpfr_root_fdfsolver_type * mp_odrf_mpfr_root_fdfsolver_steffenson = &steffenson_type;

/* end of file */
