/*
   Part of: Multiple Precision One-Dimensional Root-Finding
   Contents: tests for the one-dimensional root finding
   Date: Sun Mar  8, 2009

   Abstract

	Tests for bracketing algorithms.

   Copyright (c) 2009, 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

   This program is free software:  you can redistribute it and/or modify
   it under the terms of the  GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or (at
   your option) any later version.

   This program is  distributed in the hope that it  will be useful, but
   WITHOUT  ANY   WARRANTY;  without   even  the  implied   warranty  of
   MERCHANTABILITY  or FITNESS FOR  A PARTICULAR  PURPOSE.  See  the GNU
   General Public License for more details.

   You should  have received a  copy of  the GNU General  Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/** --------------------------------------------------------------------
 ** Headers.
 ** ----------------------------------------------------------------- */

#define DEBUGGING		0
#include <mp-odrf.h>
#include <test.h>
#include <debug.h>

typedef struct {
  const mp_odrf_mpfr_root_fsolver_driver_t *	driver;
  double			x_lower;
  double			x_upper;
  mp_odrf_mpfr_function_fun_t *	function;
  const char *			description;
} bracket_meta_data_tag_t;
typedef bracket_meta_data_tag_t *	bracket_meta_data_t;

int VERBOSE=0;

/* Solve the problems with all  the root bracketing algorithms, for both
   the sine and minus sine functions.  The argument DATA must be already
   initialised with  the selected  algorithm driver; this  function will
   reuse it for the 2 math functions. */
static void doit (bracket_meta_data_t data);

/* Solve  the problem  using the  preset bracketing  algorithm and  math
   function.  The  argument DATA  must be  already initialised  with the
   selected  algorithm  driver  and  the selected  math  function;  this
   function will reuse it for multiple initial brackets. */
static void subdoit (bracket_meta_data_t data);

/* Trigonometric sine and minus  trigonometric sine functions wrapped to
   be used by the root bracketing  algorithms.  We know that the root it
   at zero. */
static int sine_function       (mpfr_t y, mpfr_t x, void * p MP_ODRF_UNUSED);
static int minus_sine_function (mpfr_t y, mpfr_t x, void * p MP_ODRF_UNUSED);

static void test_with_interval_criterion (bracket_meta_data_t data);
static void test_with_delta_criterion    (bracket_meta_data_t data);
static void test_with_residual_criterion (bracket_meta_data_t data);

static mp_odrf_mpfr_function_fun_t	sine_function;
static mp_odrf_mpfr_function_fun_t	minus_sine_function;


/** --------------------------------------------------------------------
 ** Main.
 ** ----------------------------------------------------------------- */

int
main (void)
{
  bracket_meta_data_tag_t	data;
  const char *			s;

  s = getenv("verbose");
  if (s && 0 == strcmp("yes", s))
    VERBOSE=1;

  title("one dimensional root finding, bisection algorithm");
  data.driver = mp_odrf_mpfr_root_fsolver_bisection;
  doit(&data);

#if 0
  title("one dimensional root finding, falsepos algorithm");
  data.driver = mp_odrf_mpfr_root_fsolver_falsepos;
  doit(&data);

  title("one dimensional root finding, brent algorithm");
  data.driver = mp_odrf_mpfr_root_fsolver_brent;
  doit(&data);
#endif

  exit(EXIT_SUCCESS);
}


/** --------------------------------------------------------------------
 ** Solve with algorithms.
 ** ----------------------------------------------------------------- */

static void
doit (bracket_meta_data_t data)
/* Solve the problem using the preset bracketing algorithm, for both the
   sine and  minus sine  functions.  The argument  DATA must  be already
   initialised with  the selected  algorithm driver; this  function will
   reuse it for the 2 math functions. */
{
  subtitle("zero of sine function");
  data->function = sine_function;
  subdoit(data);
  subtitle("zero of minus sine function");
  data->function = minus_sine_function;
  subdoit(data);
}
static void
subdoit (bracket_meta_data_t data)
/* Solve  the problem  using the  preset bracketing  algorithm and  math
   function.  The  argument DATA  must be  already initialised  with the
   selected  algorithm  driver  and  the selected  math  function;  this
   function will reuse it for multiple initial brackets. */
{
  data->x_lower		= -1.0;
  data->x_upper		= +1.0;
  data->description	= "symmetric initial interval";
  test_with_interval_criterion(data);
  data->x_lower		= -1.0;
  data->x_upper		= +0.5;
  data->description	= "leftist initial interval";
  test_with_interval_criterion(data);
  data->x_lower		= -0.5;
  data->x_upper		= +1.0;
  data->description	= "rightist initial interval";
  test_with_interval_criterion(data);

  data->x_lower		= -1.0;
  data->x_upper		= +1.0;
  data->description	= "symmetric initial delta";
  test_with_delta_criterion(data);
  data->x_lower		= -1.0;
  data->x_upper		= +0.5;
  data->description	= "leftist initial delta";
  test_with_delta_criterion(data);
  data->x_lower		= -0.5;
  data->x_upper		= +1.0;
  data->description	= "rightist initial delta";
  test_with_delta_criterion(data);

  data->x_lower		= -1.0;
  data->x_upper		= +1.0;
  data->description	= "symmetric initial residual";
  test_with_residual_criterion(data);
  data->x_lower		= -1.0;
  data->x_upper		= +0.5;
  data->description	= "leftist initial residual";
  test_with_residual_criterion(data);
  data->x_lower		= -0.5;
  data->x_upper		= +1.0;
  data->description	= "rightist initial residual";
  test_with_residual_criterion(data);
}

/** --------------------------------------------------------------------
 ** Test with interval criterion.
 ** ----------------------------------------------------------------- */

static void
test_with_interval_criterion (bracket_meta_data_t data)
{
  mp_odrf_mpfr_root_fsolver_t * solver;
  mpfr_t			x_lower, x_upper;
  mpfr_t			epsabs, epsrel;
  mpfr_ptr			result;
  const mp_odrf_error_t *	E;
  mp_odrf_operation_code_t	rv;
  mp_odrf_mpfr_function_t	F = {
    .function	= data->function,
    .params	= NULL
  };
  start("interval criterion", data->description);
  E = mp_odrf_mpfr_root_fsolver_alloc(&solver, data->driver);
  if (NULL == solver) {
    perror("error initialising solver");
    exit(EXIT_FAILURE);
  }
  report("(%s) ", mp_odrf_mpfr_root_fsolver_name(solver));
  mpfr_init(x_lower);
  mpfr_init(x_upper);
  mpfr_init(epsabs);
  mpfr_init(epsrel);
  {
    mpfr_set_d(x_lower, data->x_lower, GMP_RNDN);
    mpfr_set_d(x_upper, data->x_upper, GMP_RNDN);
    mpfr_set_d(epsabs, 1e-6, GMP_RNDN);
    mpfr_set_d(epsrel, 0.0001, GMP_RNDN);

    rv = mp_odrf_mpfr_root_fsolver_set(solver, &F, x_lower, x_upper, &E);
    validate(MP_ODRF_OK == rv, "error setting: %s", E->description);
    if (MP_ODRF_OK != rv) goto end;
    if (VERBOSE)
      mpfr_fprintf(stderr, "\n- start interval\t[%30Rf, %30Rf]\n",
		   mp_odrf_mpfr_root_fsolver_x_lower(solver),
		   mp_odrf_mpfr_root_fsolver_x_upper(solver));
    do {
      rv = mp_odrf_mpfr_root_fsolver_iterate(solver, &E);
      validate(MP_ODRF_OK == rv, "error iterating: %s", E->description);
      if (MP_ODRF_OK != rv) goto end;

      if (1 == VERBOSE) {
	mpfr_fprintf(stderr, "- current interval\t[%30Rf, %30Rf]\n",
		     mp_odrf_mpfr_root_fsolver_x_lower(solver),
		     mp_odrf_mpfr_root_fsolver_x_upper(solver));
      }
      rv = mp_odrf_mpfr_root_test_interval (mp_odrf_mpfr_root_fsolver_x_lower(solver),
					    mp_odrf_mpfr_root_fsolver_x_upper(solver),
					    epsabs, epsrel, &E);
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
    result = mp_odrf_mpfr_root_fsolver_root(solver);
    if (VERBOSE)
      mpfr_fprintf(stderr, "- result %30Rf\n", result);
    validate_expected_got(0.0, result);
  }
 end:
  mpfr_clear(epsrel);
  mpfr_clear(epsabs);
  mpfr_clear(x_upper);
  mpfr_clear(x_lower);
  mp_odrf_mpfr_root_fsolver_free(solver);
  fine();
}


/** --------------------------------------------------------------------
 ** Test with delta criterion.
 ** ----------------------------------------------------------------- */

static void
test_with_delta_criterion (bracket_meta_data_t data)
{
  mp_odrf_mpfr_root_fsolver_t *	solver;
  mpfr_t			x_lower, x_upper;
  mpfr_t			epsabs, epsrel, x1;
  mpfr_ptr			result;
  const mp_odrf_error_t *	E;
  mp_odrf_operation_code_t	rv;
  mp_odrf_mpfr_function_t	F = {
    .function	= data->function,
    .params	= NULL
  };
  start("delta criterion", data->description);
  E = mp_odrf_mpfr_root_fsolver_alloc(&solver, data->driver);
  if (NULL == solver) {
    perror("error initialising solver");
    exit(EXIT_FAILURE);
  }
  report("(%s) ", mp_odrf_mpfr_root_fsolver_name(solver));
  mpfr_init(x_lower);
  mpfr_init(x_upper);
  mpfr_init(epsabs);
  mpfr_init(epsrel);
  mpfr_init(x1);
  {
    mpfr_set_d(x_lower, data->x_lower, GMP_RNDN);
    mpfr_set_d(x_upper, data->x_upper, GMP_RNDN);
    mpfr_set_d(epsabs, 1e-6, GMP_RNDN);
    mpfr_set_d(epsrel, 0.0001, GMP_RNDN);
    mpfr_set_d(x1, data->x_lower, GMP_RNDN);

    rv = mp_odrf_mpfr_root_fsolver_set(solver, &F, x_lower, x_upper, &E);
    validate(MP_ODRF_OK == rv, "error setting: %s", E->description);
    if (MP_ODRF_OK != rv) goto end;
    if (VERBOSE)
      mpfr_fprintf(stderr, "\n- start interval\t[%30Rf, %30Rf]\n",
		   mp_odrf_mpfr_root_fsolver_x_lower(solver),
		   mp_odrf_mpfr_root_fsolver_x_upper(solver));
    do {
      rv = mp_odrf_mpfr_root_fsolver_iterate(solver, &E);
      validate(MP_ODRF_OK == rv, "error iterating: %s", E->description);
      if (MP_ODRF_OK != rv) goto end;

      if (VERBOSE) {
	mpfr_fprintf(stderr, "- current interval\t[%30Rf, %30Rf]\n",
		     mp_odrf_mpfr_root_fsolver_x_lower(solver),
		     mp_odrf_mpfr_root_fsolver_x_upper(solver));
	mpfr_fprintf(stderr, "- current values: x1 = %Rf, x2 = %Rf\n",
		     x1, mp_odrf_mpfr_root_fsolver_root(solver));
      }
      rv = mp_odrf_mpfr_root_test_delta(x1, mp_odrf_mpfr_root_fsolver_root(solver),
					epsabs, epsrel, &E);
      switch (rv) {
      case MP_ODRF_OK:
	goto solved;
      case MP_ODRF_CONTINUE:
	mpfr_set(x1, mp_odrf_mpfr_root_fsolver_root(solver), GMP_RNDN);
	break;
      default:
	error(E->description);
	goto end;
      }
    } while (MP_ODRF_CONTINUE == rv);
 solved:
    result = mp_odrf_mpfr_root_fsolver_root(solver);
    if (VERBOSE)
      mpfr_fprintf(stderr, "- result %30Rf\n", result);
    validate_expected_got(0.0, result);
  }
 end:
  mpfr_clear(x1);
  mpfr_clear(epsrel);
  mpfr_clear(epsabs);
  mpfr_clear(x_upper);
  mpfr_clear(x_lower);
  mp_odrf_mpfr_root_fsolver_free(solver);
  fine();
}


/** --------------------------------------------------------------------
 ** Test with residual criterion.
 ** ----------------------------------------------------------------- */

static void
test_with_residual_criterion (bracket_meta_data_t data)
{
  mp_odrf_mpfr_root_fsolver_t *	solver;
  mpfr_t			x_lower, x_upper;
  mpfr_t			epsabs, residual;
  mpfr_ptr			result;
  const mp_odrf_error_t *	E;
  mp_odrf_operation_code_t	rv;
  mp_odrf_mpfr_function_t	F = {
    .function	= data->function,
    .params	= NULL
  };
  start("residual criterion", data->description);
  E = mp_odrf_mpfr_root_fsolver_alloc(&solver, data->driver);
  if (NULL == solver) {
    perror("error initialising solver");
    exit(EXIT_FAILURE);
  }
  report("(%s) ", mp_odrf_mpfr_root_fsolver_name(solver));
  mpfr_init(x_lower);
  mpfr_init(x_upper);
  mpfr_init(epsabs);
  mpfr_init(residual);
  {
    mpfr_set_d(x_lower, data->x_lower, GMP_RNDN);
    mpfr_set_d(x_upper, data->x_upper, GMP_RNDN);
    mpfr_set_d(epsabs, 1e-6, GMP_RNDN);

    rv = mp_odrf_mpfr_root_fsolver_set(solver, &F, x_lower, x_upper, &E);
    validate(MP_ODRF_OK == rv, "error setting: %s", E->description);
    if (MP_ODRF_OK != rv) goto end;
    if (VERBOSE)
      mpfr_fprintf(stderr, "\n- start interval\t[%30Rf, %30Rf]\n",
		   mp_odrf_mpfr_root_fsolver_x_lower(solver),
		   mp_odrf_mpfr_root_fsolver_x_upper(solver));
    do {
      rv = mp_odrf_mpfr_root_fsolver_iterate(solver, &E);
      validate(MP_ODRF_OK == rv, "error iterating: %s", E->description);
      if (MP_ODRF_OK != rv) goto end;

      MP_ODRF_MPFR_FN_EVAL(&F,residual,mp_odrf_mpfr_root_fsolver_root(solver));
      if (VERBOSE) {
	mpfr_fprintf(stderr, "- current interval\t[%30Rf, %30Rf]\n",
		     mp_odrf_mpfr_root_fsolver_x_lower(solver),
		     mp_odrf_mpfr_root_fsolver_x_upper(solver));
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
    result = mp_odrf_mpfr_root_fsolver_root(solver);
    if (VERBOSE)
      mpfr_fprintf(stderr, "- result %30Rf\n", result);
    validate_expected_got(0.0, result);
  }
 end:
  mpfr_clear(residual);
  mpfr_clear(epsabs);
  mpfr_clear(x_upper);
  mpfr_clear(x_lower);
  mp_odrf_mpfr_root_fsolver_free(solver);
  fine();
}


/** --------------------------------------------------------------------
 ** Math functions.
 ** ----------------------------------------------------------------- */

/* We know  that the root is  at zero.  So  we will test the  result for
   zero. */

static int
sine_function (mpfr_t y, mpfr_t x, void * p MP_ODRF_UNUSED)
{
  mpfr_sin(y, x, GMP_RNDN);
  return MP_ODRF_OK;
}
static int
minus_sine_function (mpfr_t y, mpfr_t x, void * p MP_ODRF_UNUSED)
{
  mpfr_sin(y, x, GMP_RNDN);
  mpfr_neg(y, y, GMP_RNDN);
  return MP_ODRF_OK;
}

/* end of file */
