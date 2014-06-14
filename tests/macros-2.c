/*
  Part of: Multiple Precision One Dimensional Root Finding
  Contents: tests for macros
  Date: Fri Jun 13, 2014

  Abstract



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

#include <mp-odrf.h>

typedef struct {
  mpfr_t        ell;
} parameters_t;

static void
parameters_init (parameters_t * params)
{
  mpfr_init(params->ell);
  mpfr_set_d(params->ell, 2.0, GMP_RNDN);
}
static void
parameters_final (parameters_t * params)
{
  mpfr_clear(params->ell);
}

static int
function (mpfr_ptr Y, mpfr_ptr X, void * params_)
{
  mpfr_t        tmp;
  parameters_t *params = params_;

  mpfr_init(tmp);
  {
    /* exp(ell * x) */
    mpfr_mul(tmp, params->ell, X, GMP_RNDN);
    mpfr_exp(Y, tmp, GMP_RNDN);
  }
  mpfr_clear(tmp);
  return MP_ODRF_OK;
}

static int
derivative (mpfr_ptr DY, mpfr_ptr X, void * params_)
{
  mpfr_t        tmp1, tmp2;
  parameters_t *params = params_;

  mpfr_init(tmp1);
  mpfr_init(tmp2);
  {
    /* ell * exp(ell * x) */
    mpfr_mul(tmp1, params->ell, X, GMP_RNDN);
    mpfr_exp(tmp2, tmp1, GMP_RNDN);
    mpfr_mul(DY, params->ell, tmp2, GMP_RNDN);
  }
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  return MP_ODRF_OK;
}

static int
function_and_derivative (mpfr_ptr DY, mpfr_ptr Y,
			 mpfr_ptr X, void * params_)
{
  mpfr_t        tmp;
  parameters_t *params = params_;

  mpfr_init(tmp);
  {
    /*
      y  = exp(ell * x)
      dy = ell * exp(ell * x) = ell * y
    */
    mpfr_mul(tmp, params->ell, X, GMP_RNDN);
    mpfr_exp(Y, tmp, GMP_RNDN);
    mpfr_mul(DY, params->ell, Y, GMP_RNDN);
  }
  mpfr_clear(tmp);
  return MP_ODRF_OK;
}

int
main (int argc MP_ODRF_UNUSED, char * argv [])
{
  parameters_t                  params;
  mpfr_t                        DY, Y, X;
  mp_odrf_mpfr_function_fdf_t   FDF = {
    .f          = function,
    .df         = derivative,
    .fdf        = function_and_derivative,
    .params     = &params
  };
  double                        x = 2.4;
  int                           rv;

  parameters_init(&params);
  mpfr_init(DY);
  mpfr_init(Y);
  mpfr_init(X);
  {
    mpfr_set_d(X, 2.4, GMP_RNDN);
    rv = MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(&FDF, DY, Y, X);
    mpfr_fprintf(stderr,
                 "%s: y = %20RNf, dy = %20RNf, should be %f and %f\n",
                 argv[0], Y, DY,
                 (exp(2.0 * x)),
                 (2.0 * exp(2.0 * x)));
    rv = MP_ODRF_MPFR_FN_FDF_EVAL_F(&FDF, Y, X);
    mpfr_fprintf(stderr, "%s: y  = %20RNf\n", argv[0], Y);
    rv = MP_ODRF_MPFR_FN_FDF_EVAL_DF(&FDF, DY, X);
    mpfr_fprintf(stderr, "%s: dy = %20RNf\n", argv[0], DY);
  }
  mpfr_clear(X);
  mpfr_clear(Y);
  mpfr_clear(DY);
  parameters_final(&params);
  exit(EXIT_SUCCESS);
}

/* end of file */
