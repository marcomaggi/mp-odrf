/*
   Part of: MPGSL
   Contents: tests for the one-dimensional root finding
   Date: Mon Mar  9, 2009

   Abstract

	Root polishing algorithms.

   Copyright (c) 2009 Marco Maggi <marcomaggi@gna.org>

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

#include <mpgsl.h>
#include <debug.h>
#include <test.h>

typedef struct {
  const mpgsl_root_fdfsolver_type *	type;
  double			initial_guess;
  mpgsl_function_fun_t *	function;
  mpgsl_function_fun_t *	derivative;
  mpgsl_function_deriv_fun_t *	function_and_derivative;
  const char *			description;
} polish_meta_data_tag_t;
typedef polish_meta_data_tag_t *	polish_meta_data_t;

int VERBOSE = 0;


/** --------------------------------------------------------------------
 ** Math functions.
 ** ----------------------------------------------------------------- */

static int
sine_function (mpfr_t y, mpfr_t x, void * p MPGSL_UNUSED)
{
  mpfr_sin(y, x, GMP_RNDN);
  return GSL_SUCCESS;
}
static int
cosine_function (mpfr_t y, mpfr_t x, void * p MPGSL_UNUSED)
{
  mpfr_cos(y, x, GMP_RNDN);
  return GSL_SUCCESS;
}
static int
sine_and_cosine_function (mpfr_t x, void * p MPGSL_UNUSED, mpfr_t y, mpfr_t dy)
{
  mpfr_sin(y,  x, GMP_RNDN);
  mpfr_cos(dy, x, GMP_RNDN);
  return GSL_SUCCESS;
}


/** --------------------------------------------------------------------
 ** Test with delta criterion.
 ** ----------------------------------------------------------------- */

static void
test_with_delta_criterion (polish_meta_data_t data)
{
  mpgsl_root_fdfsolver * solver;
  mpfr_t	guess;
  mpfr_t	epsabs, epsrel, x1;
  mpfr_ptr	result;
  int		e;
  mpgsl_function_fdf	F = {
    .f		= data->function,
    .df		= data->derivative,
    .fdf	= data->function_and_derivative,
    .params	= NULL
  };
  start("delta criterion", data->description);
  solver = mpgsl_root_fdfsolver_alloc(data->type);
  if (NULL == solver) {
    perror("error initialising solver");
    exit(EXIT_FAILURE);
  }
  report("(%s) ", mpgsl_root_fdfsolver_name(solver));
  mpfr_init(guess);
  mpfr_init(epsabs);
  mpfr_init(epsrel);
  mpfr_init(x1);
  {
    mpfr_set_d(guess, data->initial_guess, GMP_RNDN);
    mpfr_set_d(epsabs, 1e-6, GMP_RNDN);
    mpfr_set_d(epsrel, 0.0001, GMP_RNDN);
    mpfr_set_d(x1, data->initial_guess, GMP_RNDN);

    debug("setting");
    e = mpgsl_root_fdfsolver_set(solver, &F, guess);
    validate(GSL_SUCCESS == e, "error setting: %s", gsl_strerror(e));
    if (GSL_SUCCESS != e) goto end;
    if (VERBOSE)
      mpfr_fprintf(stderr, "\n- initial guess\t%30Rf\n",
		   mpgsl_root_fdfsolver_root(solver));
    debug("starting iteration");
    do {
      debug("iteration");
      e = mpgsl_root_fdfsolver_iterate(solver);
      validate(GSL_SUCCESS == e, "error iterating: %s", gsl_strerror(e));
      if (GSL_SUCCESS != e) goto end;

      debug("testing");
      if (VERBOSE)
	mpfr_fprintf(stderr, "- current values: x1 = %Rf, x0 = %Rf\n",
		     x1, mpgsl_root_fdfsolver_root(solver));
      e = mpgsl_root_test_delta(x1, mpgsl_root_fdfsolver_root(solver),
				epsabs, epsrel);
      switch (e) {
      case GSL_SUCCESS:
	goto solved;
      case GSL_CONTINUE:
	mpfr_set(x1, mpgsl_root_fdfsolver_root(solver), GMP_RNDN);
	break;
      default:
	error(gsl_strerror(e));
	goto end;
      }
    } while (GSL_CONTINUE == e);
 solved:
    result = mpgsl_root_fdfsolver_root(solver);
    if (VERBOSE)
      mpfr_fprintf(stderr, "- result %30Rf\n", result);
    validate_expected_got(0.0, result);
  }
 end:
  mpfr_clear(x1);
  mpfr_clear(epsrel);
  mpfr_clear(epsabs);
  mpfr_clear(guess);
  mpgsl_root_fdfsolver_free(solver);
  fine();
}


/** --------------------------------------------------------------------
 ** Test with residual criterion.
 ** ----------------------------------------------------------------- */

static void
test_with_residual_criterion (polish_meta_data_t data)
{
  mpgsl_root_fdfsolver * solver;
  mpfr_t	guess;
  mpfr_t	epsabs, residual;
  mpfr_ptr	result;
  int		e;
  mpgsl_function_fdf	F = {
    .f		= data->function,
    .df		= data->derivative,
    .fdf	= data->function_and_derivative,
    .params	= NULL
  };
  start("residual criterion", data->description);
  solver = mpgsl_root_fdfsolver_alloc(data->type);
  if (NULL == solver) {
    perror("error initialising solver");
    exit(EXIT_FAILURE);
  }
  report("(%s) ", mpgsl_root_fdfsolver_name(solver));
  mpfr_init(guess);
  mpfr_init(epsabs);
  mpfr_init(residual);
  {
    mpfr_set_d(guess, data->initial_guess, GMP_RNDN);
    mpfr_set_d(epsabs, 1e-6, GMP_RNDN);

    debug("setting");
    e = mpgsl_root_fdfsolver_set(solver, &F, guess);
    validate(GSL_SUCCESS == e, "error setting: %s", gsl_strerror(e));
    if (GSL_SUCCESS != e) goto end;
    if (VERBOSE)
      mpfr_fprintf(stderr, "\n- initial guess\t%30Rf\n",
		   mpgsl_root_fdfsolver_root(solver));
    debug("starting iteration");
    do {
      debug("iteration");
      e = mpgsl_root_fdfsolver_iterate(solver);
      validate(GSL_SUCCESS == e, "error iterating: %s", gsl_strerror(e));
      if (GSL_SUCCESS != e) goto end;

      debug("testing");
      MPGSL_FN_FDF_EVAL_F(&F,residual,mpgsl_root_fdfsolver_root(solver));
      if (VERBOSE) {
	mpfr_fprintf(stderr, "- current guess\t%30Rf\n",
		     mpgsl_root_fdfsolver_root(solver));
	mpfr_fprintf(stderr, "- current residual: %Rf\n", residual);
      }
      e = mpgsl_root_test_residual(residual, epsabs);
      switch (e) {
      case GSL_SUCCESS:
	goto solved;
      case GSL_CONTINUE:
	break;
      default:
	error(gsl_strerror(e));
	goto end;
      }
    } while (GSL_CONTINUE == e);
 solved:
    result = mpgsl_root_fdfsolver_root(solver);
    if (VERBOSE)
      mpfr_fprintf(stderr, "- result %30Rf\n", result);
    validate_expected_got(0.0, result);
  }
 end:
  mpfr_clear(residual);
  mpfr_clear(epsabs);
  mpfr_clear(guess);
  mpgsl_root_fdfsolver_free(solver);
  fine();
}


/** --------------------------------------------------------------------
 ** Main.
 ** ----------------------------------------------------------------- */

static void
subdoit (polish_meta_data_t data)
{
  data->initial_guess	= -1.0;
  data->description	= "leftist initial guess";
  test_with_delta_criterion(data);
  data->initial_guess	= +1.0;
  data->description	= "rightist initial guess";
  test_with_delta_criterion(data);

  data->initial_guess	= -1.0;
  data->description	= "leftist initial guess";
  test_with_residual_criterion(data);
  data->initial_guess	= +1.0;
  data->description	= "rightist initial guess";
  test_with_residual_criterion(data);
}
static void
doit (polish_meta_data_t data)
{
  subtitle("zero of sine function");
  data->function		= sine_function;
  data->derivative		= cosine_function;
  data->function_and_derivative	= sine_and_cosine_function;
  subdoit(data);
}
int
main (void)
{
  polish_meta_data_tag_t	data;
  const char *	s;

  s = getenv("verbose");
  if (s && 0 == strcmp("yes", s))
    VERBOSE=1;

  title("one dimensional root finding, newton algorithm");
  data.type = mpgsl_root_fdfsolver_newton;
  doit(&data);

  title("one dimensional root finding, secant algorithm");
  data.type		= mpgsl_root_fdfsolver_secant;
  doit(&data);

  title("one dimensional root finding, steffenson algorithm");
  data.type		= mpgsl_root_fdfsolver_steffenson;
  doit(&data);

  exit(EXIT_SUCCESS);
}

/* end of file */
