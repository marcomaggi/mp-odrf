/*
  Part of: Multiple Precision One Dimensional Root Finding
  Contents: documentation example
  Date: Sun Jun 15, 2014

  Abstract

	Documentation example of root  bracketing problem with bisection
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

int
main (void)
{
  mp_odrf_mpfr_root_fsolver_t *	solver;
  mp_odrf_mpfr_function_t	F = {
    .function = sine_function,
    .params   = NULL
  };
  mpfr_t	x_lower, x_upper;
  mpfr_t	epsabs, epsrel;
  int		rv, iteration_count = 0;

  printf("*** One-dimensional root finding:\n\
\tbisection algorithm,\n\
\tinterval stop criterion.\n");

  solver = mp_odrf_mpfr_root_fsolver_alloc(mp_odrf_mpfr_root_fsolver_bisection);
  if (NULL == solver) {
    fprintf(stderr, "error: %s\n",
	    mp_odrf_strerror(MP_ODRF_ERROR_NO_MEMORY_FOR_STATE_STRUCT));
    exit(EXIT_FAILURE);
  }

  mpfr_init(x_lower);
  mpfr_init(x_upper);
  mpfr_init(epsabs);
  mpfr_init(epsrel);
  {
    mpfr_set_d(x_lower, -1.0, GMP_RNDN);
    mpfr_set_d(x_upper, +0.5, GMP_RNDN);
    mpfr_set_d(epsabs,  1e-6, GMP_RNDN);
    mpfr_set_d(epsrel,  1e-3, GMP_RNDN);

    rv = mp_odrf_mpfr_root_fsolver_set(solver, &F, x_lower, x_upper);
    if (MP_ODRF_OK !=rv) {
      fprintf(stderr, "error setting up: %s\n",
	      mp_odrf_strerror(rv));
      goto end;
    }

    mpfr_printf("iteration %2d: [%-26Rf, %-26Rf]\n",
                iteration_count++,
                mp_odrf_mpfr_root_fsolver_x_lower(solver),
                mp_odrf_mpfr_root_fsolver_x_upper(solver));

    do {
      rv = mp_odrf_mpfr_root_fsolver_iterate(solver);
      if (MP_ODRF_OK != rv) {
        fprintf(stderr, "error iterating: %s\n",
		mp_odrf_strerror(rv));
        goto end;
      }

      mpfr_printf("iteration %2d: [%-26Rf, %-26Rf]\n",
                  iteration_count++,
                  mp_odrf_mpfr_root_fsolver_x_lower(solver),
                  mp_odrf_mpfr_root_fsolver_x_upper(solver));

      rv = mp_odrf_mpfr_root_test_interval(mp_odrf_mpfr_root_fsolver_x_lower(solver),
					   mp_odrf_mpfr_root_fsolver_x_upper(solver),
					   epsabs, epsrel);
      switch (rv) {
      case MP_ODRF_OK:
        goto solved;
      case MP_ODRF_CONTINUE:
        break;
      default:
        fprintf(stderr, "error testing: %s\n", mp_odrf_strerror(rv));
        goto end;
      }
    } while (MP_ODRF_CONTINUE == rv);
  solved:
    mpfr_printf("result = %Rf\n", mp_odrf_mpfr_root_fsolver_root(solver));
  }
 end:
  mpfr_clear(epsrel);
  mpfr_clear(epsabs);
  mpfr_clear(x_upper);
  mpfr_clear(x_lower);
  mp_odrf_mpfr_root_fsolver_free(solver);
  exit(EXIT_SUCCESS);
}

/* end of file */
