/*
  Part of: Multiple Precision One Dimensional Root Finding
  Contents: documentation example
  Date: Sun Jun 15, 2014

  Abstract

	Documentation example of root  polishing problem with Steffenson
	algorithm.

  Copyright (C) 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

  This program is  free software: you can redistribute  it and/or modify
  it under the  terms of the GNU General Public  License as published by
  the Free Software Foundation, either version  3 of the License, or (at
  your option) any later version.

  This program  is distributed in the  hope that it will  be useful, but
  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See the  GNU
  General Public License for more details.

  You should  have received  a copy  of the  GNU General  Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mp-odrf.h>

static int
sine_function (mpfr_ptr Y, mpfr_ptr X, void * params_ MP_ODRF_UNUSED)
{
  mpfr_sin(Y, X, GMP_RNDN);
  return MP_ODRF_OK;
}
static int
cosine_function (mpfr_ptr Y, mpfr_ptr X, void * params_ MP_ODRF_UNUSED)
{
  mpfr_cos(Y, X, GMP_RNDN);
  return MP_ODRF_OK;
}
static int
sine_and_cosine_function (mpfr_ptr DY, mpfr_ptr Y,
			  mpfr_ptr X, void * params_ MP_ODRF_UNUSED)
{
  mpfr_sin(Y,  X, GMP_RNDN);
  mpfr_cos(DY, X, GMP_RNDN);
  return MP_ODRF_OK;
}

int
main (void)
{
  mp_odrf_mpfr_root_fdfsolver_t *	solver;
  mp_odrf_mpfr_function_fdf_t		FDF = {
    .f      = sine_function,
    .df     = cosine_function,
    .fdf    = sine_and_cosine_function,
    .params = NULL
  };
  mpfr_t		epsabs, residual, initial_guess;
  int			rv, iteration_count = 0;

  printf("*** One-dimensional root finding:\n\
\tsteffenson algorithm,\n\
\tresidual stop criterion.\n");

  solver = mp_odrf_mpfr_root_fdfsolver_alloc(mp_odrf_mpfr_root_fdfsolver_steffenson);
  if (NULL == solver) {
    fprintf(stderr, "error: %s\n",
	    mp_odrf_strerror(MP_ODRF_ERROR_NO_MEMORY_FOR_STATE_STRUCT));
    exit(EXIT_FAILURE);
  }

  mpfr_init(initial_guess);
  mpfr_init(epsabs);
  mpfr_init(residual);
  {
    mpfr_set_d(initial_guess, -1.0, GMP_RNDN);
    mpfr_set_d(epsabs,        1e-6, GMP_RNDN);

    rv = mp_odrf_mpfr_root_fdfsolver_set(solver, &FDF, initial_guess);
    if (MP_ODRF_OK != rv) {
      fprintf(stderr, "error setting up: %s\n", mp_odrf_strerror(rv));
      goto end;
    }

    mpfr_printf("iteration %2d: x = %-32Rf\n",
                iteration_count++,
                mp_odrf_mpfr_root_fdfsolver_root(solver));

    do {
      rv = mp_odrf_mpfr_root_fdfsolver_iterate(solver);
      if (MP_ODRF_OK != rv) {
        fprintf(stderr, "error iterating: %s\n", mp_odrf_strerror(rv));
        goto end;
      }

      mpfr_printf("iteration %2d: x = %-32Rf\n",
                  iteration_count++,
                  mp_odrf_mpfr_root_fdfsolver_root(solver));

      MP_ODRF_MPFR_FN_FDF_EVAL_F(&FDF, residual,
				 mp_odrf_mpfr_root_fdfsolver_root(solver));

      rv = mp_odrf_mpfr_root_test_residual(residual, epsabs);
      switch (rv) {
      case MP_ODRF_OK:
        goto solved;
      case MP_ODRF_CONTINUE:
        break;
      default:
        fprintf(stderr, "error testing: %s\n", mp_odrf_strerror(rv));
        goto end;
      }
    } while (MP_ODRF_OK == rv);
  solved:
    mpfr_printf("result = %Rf\n", mp_odrf_mpfr_root_fdfsolver_root(solver));
  }
 end:
  mpfr_clear(epsabs);
  mpfr_clear(initial_guess);
  mp_odrf_mpfr_root_fdfsolver_free(solver);
  exit(EXIT_SUCCESS);
}

/* end of file */
