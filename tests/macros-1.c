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
  mpfr_t        a;
  mpfr_t        b;
  mpfr_t        c;
} parameters_t;

static void
parameters_init (parameters_t * params)
{
  mpfr_init(params->a);
  mpfr_init(params->b);
  mpfr_init(params->c);

  mpfr_set_d(params->a, 3.0, GMP_RNDN);
  mpfr_set_d(params->b, 2.0, GMP_RNDN);
  mpfr_set_d(params->c, 1.0, GMP_RNDN);
}
static void
parameters_final (parameters_t * params)
{
  mpfr_clear(params->a);
  mpfr_clear(params->b);
  mpfr_clear(params->c);
}

static int
function (mpfr_ptr Y, mpfr_ptr X, void * params_)
{
  mpfr_t        tmp1;
  mpfr_t        tmp2;
  parameters_t *params = params_;

  mpfr_init(tmp1);
  mpfr_init(tmp2);
  {
    /* (a * x + b) * x + c */
    mpfr_mul(tmp1, params->a, X, GMP_RNDN);
    mpfr_add(tmp2, tmp1, params->b, GMP_RNDN);
    mpfr_mul(tmp1, tmp2, X, GMP_RNDN);
    mpfr_add(Y, tmp1, params->c, GMP_RNDN);
  }
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  return MP_ODRF_OK;
}

static int
derivative (mpfr_ptr DY, mpfr_ptr X, void * params_)
{
  mpfr_t        tmp1;
  mpfr_t        tmp2;
  parameters_t *params = params_;

  mpfr_init(tmp1);
  mpfr_init(tmp2);
  {
    /* 2 * a * x + b */
    mpfr_mul_d(tmp1, params->a, 2.0, GMP_RNDN);
    mpfr_mul(tmp2, tmp1, X, GMP_RNDN);
    mpfr_add(DY, tmp2, params->b, GMP_RNDN);
  }
  mpfr_clear(tmp1);
  mpfr_clear(tmp2);
  return MP_ODRF_OK;
}

static int
function_and_derivative (mpfr_ptr DY, mpfr_ptr Y,
			 mpfr_ptr X, void * params_)
{
  int	rv;
  rv = function(Y, X, params_);
  if (MP_ODRF_OK == rv) {
     rv = derivative(DY, X, params_);
  }
  return rv;
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
  int                           rv;
  double                        x = 2.4;

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
                 (3.0 * x * x + 2.0 * x + 1.0),
                 (2.0 * 3.0 * x + 2.0));
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
