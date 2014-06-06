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

