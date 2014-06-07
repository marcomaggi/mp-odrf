/*
   Part of: Multiple Precision One-Dimensional Root-Finding
   Contents: tests for the one-dimensional root finding
   Date: Mon Mar  9, 2009

   Abstract

	Tests for root polishing algorithms.

   Copyright (c) 2009, 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

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
#include <mp-odrf.h>
#include <debug.h>
#include <test.h>

typedef struct {
  const mp_odrf_mpfr_root_fdfsolver_driver_t *	driver;
  double				initial_guess;
  mp_odrf_mpfr_wrapped_f_t *		function;
  mp_odrf_mpfr_wrapped_f_t *		derivative;
  mp_odrf_mpfr_wrapped_fdf_t *		function_and_derivative;
  const char *				description;
} polish_meta_data_tag_t;
typedef polish_meta_data_tag_t *	polish_meta_data_t;

int VERBOSE = 0;

/* ------------------------------------------------------------------ */

/* Solve the problem using the  preset polishing algorithm, for the sine
   function.  The  argument DATA  must be  already initialised  with the
   selected algorithm driver. */
static void doit (polish_meta_data_t data);

/* Solve  the problem  using  the preset  polishing  algorithm and  math
   function.  The  argument DATA  must be  already initialised  with the
   selected  algorithm  driver  and  the selected  math  function;  this
   function will reuse it for multiple initial root guesses and multiple
   convergence tests. */
static void subdoit (polish_meta_data_t data);

static void test_with_delta_criterion (polish_meta_data_t data);
static void test_with_residual_criterion (polish_meta_data_t data);

/* Trigonometric sine function wrapped to  be used by the root polishing
   algorithm.  This is the target function;  we know that the root is at
   zero. */
static mp_odrf_mpfr_wrapped_f_t		sine_function;

/* Trigonometric  cosine  function  wrapped  to  be  used  by  the  root
   polishing  algorithm.    This  is   the  derivative  of   the  target
   function. */
static mp_odrf_mpfr_wrapped_f_t		cosine_function;

/* Trigonometric sine  and cosine  functions wrapped to  be used  by the
   root  polishing algorithm.   This function  computes both  the target
   function and its derivative. */
static mp_odrf_mpfr_wrapped_fdf_t	sine_and_cosine_function;


/** --------------------------------------------------------------------
 ** Main.
 ** ----------------------------------------------------------------- */

int
main (void)
{
  polish_meta_data_tag_t	data;
  const char *	s;

  s = getenv("VERBOSE");
  if (s && 0 == strcmp("yes", s))
    VERBOSE=1;

  title("one dimensional root finding, newton algorithm");
  data.driver = mp_odrf_mpfr_root_fdfsolver_newton;
  doit(&data);

  title("one dimensional root finding, secant algorithm");
  data.driver = mp_odrf_mpfr_root_fdfsolver_secant;
  doit(&data);

  title("one dimensional root finding, steffenson algorithm");
  data.driver = mp_odrf_mpfr_root_fdfsolver_steffenson;
  doit(&data);

  exit(EXIT_SUCCESS);
}


/** --------------------------------------------------------------------
 ** Solve with algorithms.
 ** ----------------------------------------------------------------- */

static void
doit (polish_meta_data_t data)
/* Solve the problem using the  preset polishing algorithm, for the sine
   function.  The  argument DATA  must be  already initialised  with the
   selected algorithm driver. */
{
  subtitle("zero of sine function");
  data->function		= sine_function;
  data->derivative		= cosine_function;
  data->function_and_derivative	= sine_and_cosine_function;
  subdoit(data);
}
static void
subdoit (polish_meta_data_t data)
/* Solve  the problem  using  the preset  polishing  algorithm and  math
   function.  The  argument DATA  must be  already initialised  with the
   selected  algorithm  driver  and  the selected  math  function;  this
   function will reuse it for multiple initial root guesses and multiple
   convergence tests. */
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


/** --------------------------------------------------------------------
 ** Test with delta criterion.
 ** ----------------------------------------------------------------- */

static void
test_with_delta_criterion (polish_meta_data_t data)
{
  mp_odrf_mpfr_root_fdfsolver_t *	solver;
  mpfr_t			guess;
  mpfr_t			epsabs, epsrel, x1;
  mpfr_ptr			result;
  const mp_odrf_error_t *	E;
  mp_odrf_operation_code_t	rv;
  mp_odrf_mpfr_function_fdf_t	FDF = {
    .f		= data->function,
    .df		= data->derivative,
    .fdf	= data->function_and_derivative,
    .params	= NULL
  };
  start("delta criterion", data->description);
  E = mp_odrf_mpfr_root_fdfsolver_alloc(&solver, data->driver);
  if (NULL == solver) {
    perror("error initialising solver");
    exit(EXIT_FAILURE);
  }
  report("(%s) ", mp_odrf_mpfr_root_fdfsolver_name(solver));
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
    rv = mp_odrf_mpfr_root_fdfsolver_set(solver, &FDF, guess, &E);
    validate(MP_ODRF_OK == rv, "error setting: %s", E->description);
    if (MP_ODRF_OK != rv) goto end;
    if (VERBOSE)
      mpfr_fprintf(stderr, "\n- initial guess\t%30Rf\n",
		   mp_odrf_mpfr_root_fdfsolver_root(solver));
    debug("starting iteration");
    do {
      debug("iteration");
      rv = mp_odrf_mpfr_root_fdfsolver_iterate(solver, &E);
      validate(MP_ODRF_OK == rv, "error iterating: %s", E->description);
      if (MP_ODRF_OK != rv) goto end;

      debug("testing");
      if (VERBOSE)
	mpfr_fprintf(stderr, "- current values: x1 = %Rf, x0 = %Rf\n",
		     x1, mp_odrf_mpfr_root_fdfsolver_root(solver));
      rv = mp_odrf_mpfr_root_test_delta(x1, mp_odrf_mpfr_root_fdfsolver_root(solver),
					epsabs, epsrel, &E);
      switch (rv) {
      case MP_ODRF_OK:
	goto solved;
      case MP_ODRF_CONTINUE:
	mpfr_set(x1, mp_odrf_mpfr_root_fdfsolver_root(solver), GMP_RNDN);
	break;
      default:
	error(E->description);
	goto end;
      }
    } while (MP_ODRF_CONTINUE == rv);
  solved:
    result = mp_odrf_mpfr_root_fdfsolver_root(solver);
    if (VERBOSE)
      mpfr_fprintf(stderr, "- result %30Rf\n", result);
    validate_expected_got(0.0, result);
  }
 end:
  mpfr_clear(x1);
  mpfr_clear(epsrel);
  mpfr_clear(epsabs);
  mpfr_clear(guess);
  mp_odrf_mpfr_root_fdfsolver_free(solver);
  fine();
}


/** --------------------------------------------------------------------
 ** Test with residual criterion.
 ** ----------------------------------------------------------------- */

static void
test_with_residual_criterion (polish_meta_data_t data)
{
  mp_odrf_mpfr_root_fdfsolver_t *	solver;
  mpfr_t			guess;
  mpfr_t			epsabs, residual;
  mpfr_ptr			result;
  const mp_odrf_error_t *	E;
  mp_odrf_operation_code_t	rv;
  mp_odrf_mpfr_function_fdf_t	FDF = {
    .f		= data->function,
    .df		= data->derivative,
    .fdf	= data->function_and_derivative,
    .params	= NULL
  };
  start("residual criterion", data->description);
  E = mp_odrf_mpfr_root_fdfsolver_alloc(&solver, data->driver);
  if (NULL == solver) {
    perror("error initialising solver");
    exit(EXIT_FAILURE);
  }
  report("(%s) ", mp_odrf_mpfr_root_fdfsolver_name(solver));
  mpfr_init(guess);
  mpfr_init(epsabs);
  mpfr_init(residual);
  {
    mpfr_set_d(guess, data->initial_guess, GMP_RNDN);
    mpfr_set_d(epsabs, 1e-6, GMP_RNDN);

    debug("setting");
    rv = mp_odrf_mpfr_root_fdfsolver_set(solver, &FDF, guess, &E);
    validate(MP_ODRF_OK == rv, "error setting: %s", E->description);
    if (MP_ODRF_OK != rv) goto end;
    if (VERBOSE)
      mpfr_fprintf(stderr, "\n- initial guess\t%30Rf\n",
		   mp_odrf_mpfr_root_fdfsolver_root(solver));
    debug("starting iteration");
    do {
      debug("iteration");
      rv = mp_odrf_mpfr_root_fdfsolver_iterate(solver, &E);
      validate(MP_ODRF_OK == rv, "error iterating: %s", E->description);
      if (MP_ODRF_OK != rv) goto end;

      debug("testing");
      MP_ODRF_MPFR_FN_FDF_EVAL_F(&FDF, residual,
				 mp_odrf_mpfr_root_fdfsolver_root(solver), &E);
      if (VERBOSE) {
	mpfr_fprintf(stderr, "- current guess\t%30Rf\n",
		     mp_odrf_mpfr_root_fdfsolver_root(solver));
	mpfr_fprintf(stderr, "- current residual: %Rf\n", residual);
      }
      rv = mp_odrf_mpfr_root_test_residual(residual, epsabs, &E);
      switch (rv) {
      case MP_ODRF_OK:
	goto solved;
      case MP_ODRF_CONTINUE:
	break;
      default:
	error(E->description);
	goto end;
      }
    } while (MP_ODRF_CONTINUE == rv);
 solved:
    result = mp_odrf_mpfr_root_fdfsolver_root(solver);
    if (VERBOSE)
      mpfr_fprintf(stderr, "- result %30Rf\n", result);
    validate_expected_got(0.0, result);
  }
 end:
  mpfr_clear(residual);
  mpfr_clear(epsabs);
  mpfr_clear(guess);
  mp_odrf_mpfr_root_fdfsolver_free(solver);
  fine();
}


/** --------------------------------------------------------------------
 ** Math functions.
 ** ----------------------------------------------------------------- */

/* We know  that the root is  at zero.  So  we will test the  result for
   zero. */

static mp_odrf_operation_code_t
sine_function (mpfr_t y, mpfr_t x, void * p MP_ODRF_UNUSED,
	       const mp_odrf_error_t ** E MP_ODRF_UNUSED)
{
  mpfr_sin(y, x, GMP_RNDN);
  return MP_ODRF_OK;
}
static mp_odrf_operation_code_t
cosine_function (mpfr_t y, mpfr_t x, void * p MP_ODRF_UNUSED,
		 const mp_odrf_error_t ** E MP_ODRF_UNUSED)
{
  mpfr_cos(y, x, GMP_RNDN);
  return MP_ODRF_OK;
}
static mp_odrf_operation_code_t
sine_and_cosine_function (mpfr_t x, void * p MP_ODRF_UNUSED, mpfr_t y, mpfr_t dy,
			  const mp_odrf_error_t ** E MP_ODRF_UNUSED)
{
  mpfr_sin(y,  x, GMP_RNDN);
  mpfr_cos(dy, x, GMP_RNDN);
  return MP_ODRF_OK;
}

/* end of file */
