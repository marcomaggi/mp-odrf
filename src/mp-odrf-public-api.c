/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: public API
   Date: Sat Mar  7, 2009

   Abstract

	This  module  implements  the  public  API  of  one  dimensional
	root-finding algorithms.  For the original  GSL code look in the
	GSL source tree, "roots" directory.

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


/** --------------------------------------------------------------------
 ** Root bracketing solver API.
 ** ----------------------------------------------------------------- */

mp_odrf_mpfr_root_fsolver_t *
mp_odrf_mpfr_root_fsolver_alloc (const mp_odrf_mpfr_root_fsolver_driver_t * T)
/* Allocate and initialise a new root bracketing state struct to use the
   selected algorithm driver. */
{
  mp_odrf_mpfr_root_fsolver_t * S = malloc(sizeof(mp_odrf_mpfr_root_fsolver_t));
  if (NULL != S) {
    S->driver_state = malloc(T->driver_state_size);
    if (NULL != S->driver_state) {
      T->init(S->driver_state);
      S->driver		= T;
      S->function	= NULL;
      mpfr_init(S->root);
      mpfr_init(S->x_lower);
      mpfr_init(S->x_upper);
    } else {
      free(S);
      S = NULL;
    }
  }
  return S;
}
void
mp_odrf_mpfr_root_fsolver_free (mp_odrf_mpfr_root_fsolver_t * S)
/* Finalise and release a root bracketing state struct. */
{
  mpfr_clear(S->root);
  mpfr_clear(S->x_lower);
  mpfr_clear(S->x_upper);
  S->driver->final(S->driver_state);
  free(S->driver_state);
  free(S);
}
mp_odrf_code_t
mp_odrf_mpfr_root_fsolver_set (mp_odrf_mpfr_root_fsolver_t * S,
			       mp_odrf_mpfr_function_t * F,
			       mpfr_t x_lower, mpfr_t x_upper)
/* Select the  math function to be  searched for roots for  a given root
   bracketing state struct.  Also selects the search bracket. */
{
  mp_odrf_code_t	retval = MP_ODRF_OK;
  if (mpfr_greater_p(x_lower, x_upper)) {
    retval = MP_ODRF_ERROR_INVALID_BRACKET_INTERVAL;
  } else {
    S->function = F;
    {
      mpfr_t	tmp;
      mpfr_init(tmp);
      {
	/* s->root = 0.5 * (x_lower + x_upper); */
	mpfr_add(tmp, x_lower, x_upper, GMP_RNDN);
	mpfr_mul_d(S->root, tmp, 0.5, GMP_RNDN);
      }
      mpfr_clear(tmp);
    }
    mpfr_set(S->x_lower, x_lower, GMP_RNDD);
    mpfr_set(S->x_upper, x_upper, GMP_RNDU);
    retval = (S->driver->set)(S->driver_state, S->function, S->root, x_lower, x_upper);
  }
  return retval;
}
mp_odrf_code_t
mp_odrf_mpfr_root_fsolver_iterate (mp_odrf_mpfr_root_fsolver_t * S)
/* Perform a search iteration for a root bracketing state struct. */
{
  return (S->driver->iterate) (S->driver_state, S->function, S->root,
			       S->x_lower, S->x_upper);
}
const char *
mp_odrf_mpfr_root_fsolver_name (const mp_odrf_mpfr_root_fsolver_t * S)
/* Return the name of the algorithms. */
{
  return S->driver->name;
}
mpfr_ptr
mp_odrf_mpfr_root_fsolver_root (const mp_odrf_mpfr_root_fsolver_t * S)
/* Return the current estimate solution. */
{
  return (mpfr_ptr)S->root;
}
mpfr_ptr
mp_odrf_mpfr_root_fsolver_x_lower (const mp_odrf_mpfr_root_fsolver_t * S)
/* Return the bracket's lower bound. */
{
  return (mpfr_ptr)S->x_lower;
}
mpfr_ptr
mp_odrf_mpfr_root_fsolver_x_upper (const mp_odrf_mpfr_root_fsolver_t * S)
/* Return the bracket's upper bound. */
{
  return (mpfr_ptr)S->x_upper;
}


/** --------------------------------------------------------------------
 ** Root polishing solver API.
 ** ----------------------------------------------------------------- */

mp_odrf_mpfr_root_fdfsolver_t *
mp_odrf_mpfr_root_fdfsolver_alloc (const mp_odrf_mpfr_root_fdfsolver_driver_t * T)
/* Allocate and initialise a new root  polishing state struct to use the
   selected algorithm driver. */
{
  mp_odrf_mpfr_root_fdfsolver_t * S = malloc(sizeof(mp_odrf_mpfr_root_fdfsolver_t));
  if (NULL != S) {
    S->driver_state = malloc(T->driver_state_size);
    if (NULL != S->driver_state) {
      T->init(S->driver_state);
      S->driver	= T;
      S->fdf	= NULL;
      mpfr_init(S->root);
    } else {
      free(S);
      S = NULL;
    }
  }
  return S;
}
void
mp_odrf_mpfr_root_fdfsolver_free (mp_odrf_mpfr_root_fdfsolver_t * S)
/* Finalise and release a root polishing state struct. */
{
  mpfr_clear(S->root);
  S->driver->final(S->driver_state);
  free(S->driver_state);
  free(S);
}
mp_odrf_code_t
mp_odrf_mpfr_root_fdfsolver_set (mp_odrf_mpfr_root_fdfsolver_t * S,
				 mp_odrf_mpfr_function_fdf_t * F,
				 mpfr_t root)
/* Select the  math function to be  searched for roots for  a given root
   polishing state struct.  Also selects the initial solution guess. */
{
  S->fdf = F;
  mpfr_set(S->root, root, GMP_RNDN);
  return (S->driver->set)(S->driver_state, S->fdf, S->root);
}
mp_odrf_code_t
mp_odrf_mpfr_root_fdfsolver_iterate (mp_odrf_mpfr_root_fdfsolver_t * S)
/* Perform a search iteration for a root polishing state struct. */
{
  return (S->driver->iterate) (S->driver_state, S->fdf, S->root);
}
const char *
mp_odrf_mpfr_root_fdfsolver_name (const mp_odrf_mpfr_root_fdfsolver_t * S)
/* Return the name of the algorithms. */
{
  return S->driver->name;
}
mpfr_ptr
mp_odrf_mpfr_root_fdfsolver_root (const mp_odrf_mpfr_root_fdfsolver_t * S)
/* Return the current estimate solution. */
{
  return (mpfr_ptr)S->root;
}


/** --------------------------------------------------------------------
 ** Convergence tests API.
 ** ----------------------------------------------------------------- */

mp_odrf_code_t
mp_odrf_mpfr_root_test_interval (mpfr_ptr x_lower, mpfr_ptr x_upper,
				 mpfr_ptr epsabs,  mpfr_ptr epsrel)
{
  mp_odrf_code_t	retval = MP_ODRF_CONTINUE;
  mpfr_t		abs_lower;
  mpfr_t		abs_upper;
  mpfr_t		min_abs;
  mpfr_t		tolerance;
  mpfr_t		tmp;
  int			clo, cup;
  if (mpfr_cmp_si(epsrel, 0) < 0) {
    retval	= MP_ODRF_ERROR_RELATIVE_TOLERANCE_IS_NEGATIVE;
  } else if (mpfr_cmp_si(epsabs, 0) < 0) {
    retval	= MP_ODRF_ERROR_ABSOLUTE_TOLERANCE_IS_NEGATIVE;
  } else if (mpfr_greater_p(x_lower, x_upper)) {
    retval	= MP_ODRF_ERROR_LOWER_BOUND_LARGER_THAN_UPPER_BOUND;
  } else {
    mpfr_init(tmp);
    mpfr_init(abs_lower);
    mpfr_init(abs_upper);
    mpfr_init(min_abs);
    mpfr_init(tolerance);
    {
      mpfr_abs(abs_lower, x_lower, GMP_RNDN);
      mpfr_abs(abs_upper, x_upper, GMP_RNDN);
      clo = mpfr_cmp_si(x_lower, 0);
      cup = mpfr_cmp_si(x_upper, 0);
      if (((clo > 0) && (cup > 0)) ||
	  ((clo < 0) && (cup < 0))) {
	mpfr_min(tmp, abs_lower, abs_upper, GMP_RNDN);
	mpfr_set(min_abs, tmp, GMP_RNDN);
      } else {
	mpfr_set_d(min_abs, 0.0, GMP_RNDN);
      }
      mpfr_mul(tmp, epsrel, min_abs, GMP_RNDN);
      mpfr_add(tolerance, epsabs, tmp, GMP_RNDN);
      mpfr_sub(tmp, x_upper, x_lower, GMP_RNDN);
      mpfr_abs(tmp, tmp, GMP_RNDN);
      if (mpfr_less_p(tmp, tolerance))
	retval = MP_ODRF_OK;
    }
    mpfr_clear(tmp);
    mpfr_clear(abs_lower);
    mpfr_clear(abs_upper);
    mpfr_clear(min_abs);
    mpfr_clear(tolerance);
  }
  return retval;
}
mp_odrf_code_t
mp_odrf_mpfr_root_test_delta (mpfr_ptr x1, mpfr_ptr x0,
			      mpfr_ptr epsabs, mpfr_ptr epsrel)
{
  mp_odrf_code_t	retval = MP_ODRF_CONTINUE;
  if (mpfr_cmp_si(epsrel, 0) < 0) {
    retval	= MP_ODRF_ERROR_RELATIVE_TOLERANCE_IS_NEGATIVE;
  } else if (mpfr_cmp_si(epsabs, 0) < 0) {
    retval	= MP_ODRF_ERROR_ABSOLUTE_TOLERANCE_IS_NEGATIVE;
  } else if (0 == mpfr_cmp(x1, x0)) {
    retval	= MP_ODRF_OK;
  } else {
    mpfr_t			tolerance, tmp, tmp1;
    mpfr_init(tmp);
    mpfr_init(tmp1);
    mpfr_init(tolerance);
    {
      mpfr_abs(tmp, x1, GMP_RNDN);
      mpfr_mul(tmp1, epsrel, tmp, GMP_RNDN);
      mpfr_add(tolerance, epsabs, tmp1, GMP_RNDN);
      mpfr_sub(tmp, x1, x0, GMP_RNDN);
      mpfr_abs(tmp, tmp, GMP_RNDN);
      if (mpfr_less_p(tmp, tolerance))
	retval = MP_ODRF_OK;
    }
    mpfr_clear(tmp);
    mpfr_clear(tmp1);
    mpfr_clear(tolerance);
  }
  return retval;
}
mp_odrf_code_t
mp_odrf_mpfr_root_test_residual (mpfr_ptr f, mpfr_ptr epsabs)
{
  mp_odrf_code_t	retval = MP_ODRF_CONTINUE;
  if (mpfr_cmp_si(epsabs, 0) < 0) {
    retval	= MP_ODRF_ERROR_ABSOLUTE_TOLERANCE_IS_NEGATIVE;
  } else {
    mpfr_t	tmp;
    mpfr_init(tmp);
    mpfr_abs(tmp, f, GMP_RNDN);
    if (mpfr_less_p(tmp, epsabs)) {
      retval = MP_ODRF_OK;
    }
  }
  return retval;
}


/** --------------------------------------------------------------------
 ** Comparison API.
 ** ----------------------------------------------------------------- */

int
mp_odrf_mpfr_fcmp (mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr epsilon)
/* This  is  a  MPFR-isation  of the  following  function  for  "double"
   numbers.  The original code itself is Copyright (C) 2002 Gert Van den
   Eynde, the  algorithm is  by D.E.  Knuth  "Seminumerical Algorithms",
   Section 4.2.2 (3rd Edition).

  int
  gsl_fcmp (const double x1, const double x2, const double epsilon)
  {
    int exponent;
    double delta, difference;
    double max = (fabs (x1) > fabs (x2)) ? x1 : x2;
    frexp (max, &exponent);
    delta = ldexp (epsilon, exponent);
    difference = x1 - x2;
    if (difference > delta)		return 1;
    else if (difference < -delta)	return -1;
    else				return 0;
  }
*/
{
  mpfr_t	delta, difference;
  mpfr_srcptr	max;
  int		retval;
  mp_prec_t	prec = mpfr_get_prec(epsilon);
  mpfr_init2(delta,prec);
  mpfr_init2(difference,prec);
  {
    max = (0 < mpfr_cmpabs(a, b))? a : b;
    if (mpfr_zero_p(max)) {
      retval = 0;
      goto end;
    } else {
      mpfr_mul_2si(delta, epsilon, mpfr_get_exp(max), GMP_RNDD);
      mpfr_sub(difference, a, b, GMP_RNDU);
      if (mpfr_greater_p(difference, delta)) {
	retval = 1;
	goto end;
      }
      mpfr_neg(delta, delta, GMP_RNDN);
      if (mpfr_less_p(difference, delta))
	retval = -1;
      else
	retval = 0;
    }
  }
 end:
  mpfr_clear(delta);
  mpfr_clear(difference);
  return retval;
}
int
mp_odrf_mpfr_absdiff_equal_p (mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr epsilon)
{
  mpfr_t	x;
  int		retval=1;
  mpfr_init(x);
  {
    mpfr_sub(x, a, b, mpfr_greater_p(a, b)? GMP_RNDU : GMP_RNDD);
    if (0 <= mpfr_cmpabs(x, epsilon))
      retval = 0;
  }
  mpfr_clear(x);
  return retval;
}
int
mp_odrf_mpfr_reldiff_equal_p (mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr epsilon)
{
  mpfr_t	x;
  int		retval=1;
  mpfr_init(x);
  {
    mpfr_reldiff(x, a, b, GMP_RNDN);
    if (0 <= mpfr_cmpabs(x, epsilon))
      retval = 0;
  }
  mpfr_clear(x);
  return retval;
}

/* end of file */
