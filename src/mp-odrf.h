/*
  Part of: Multiple Precision One Dimensional Root Finding
  Contents: public header file
  Date: Fri Jun  6, 2014

  Abstract

	This is the public header file.   It must be included in all the
	C sources making use of the MP ODRF API.

  Copyright (C) 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

  In the original "roots/gsl_roots.h" the definitions are: Copyright (C)
  1996, 1997, 1998, 1999, 2000, 2007 Reid Priedhorsky, Brian Gough.

  This program is  free software: you can redistribute  it and/or modify
  it under the  terms of the GNU General Public  License as published by
  the Free Software Foundation, either  version 3 of the License, or (at
  your option) any later version.

  This program  is distributed in the  hope that it will  be useful, but
  WITHOUT   ANY  WARRANTY;   without  even   the  implied   warranty  of
  MERCHANTABILITY  or FITNESS  FOR A  PARTICULAR PURPOSE.   See  the GNU
  General Public License for more details.

  You  should have received  a copy  of the  GNU General  Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef MP_ODRF_H
#define MP_ODRF_H 1


/** --------------------------------------------------------------------
 ** Headers.
 ** ----------------------------------------------------------------- */

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mpfr.h>


/** --------------------------------------------------------------------
 ** Preprocessor macros.
 ** ----------------------------------------------------------------- */

/* The macro MP_ODRF_UNUSED indicates that a function, function argument
   or variable may potentially be unused. Usage examples:

   static int unused_function (char arg) MP_ODRF_UNUSED;
   int foo (char unused_argument MP_ODRF_UNUSED);
   int unused_variable MP_ODRF_UNUSED;
*/
#ifdef __GNUC__
#  define MP_ODRF_UNUSED		__attribute__((unused))
#else
#  define MP_ODRF_UNUSED		/* empty */
#endif

#ifndef __GNUC__
#  define __attribute__(...)	/* empty */
#endif

/* I found  the following chunk on  the Net.  (Marco Maggi;  Sun Feb 26,
   2012) */
#if defined _WIN32 || defined __CYGWIN__
#  ifdef BUILDING_DLL
#    ifdef __GNUC__
#      define mp_odrf_decl		__attribute__((dllexport))
#    else
#      define mp_odrf_decl		__declspec(dllexport)
#    endif
#  else
#    ifdef __GNUC__
#      define mp_odrf_decl		__attribute__((dllimport))
#    else
#      define mp_odrf_decl		__declspec(dllimport)
#    endif
#  endif
#  define mp_odrf_private_decl		extern
#else
#  if __GNUC__ >= 4
#    define mp_odrf_decl		__attribute__((visibility ("default")))
#    define mp_odrf_private_decl	__attribute__((visibility ("hidden")))
#  else
#    define mp_odrf_decl		extern
#    define mp_odrf_private_decl	extern
#  endif
#endif


/** --------------------------------------------------------------------
 ** Error handling.
 ** ----------------------------------------------------------------- */

typedef enum {
  MP_ODRF_CONTINUE					=   1,
  MP_ODRF_OK						=   0,
  MP_ODRF_ERROR						=  -1,
  MP_ODRF_ERROR_NO_MEMORY_FOR_STATE_STRUCT		=  -2,
  MP_ODRF_ERROR_INVALID_BRACKET_INTERVAL		=  -3,
  MP_ODRF_ERROR_RELATIVE_TOLERANCE_IS_NEGATIVE		=  -4,
  MP_ODRF_ERROR_ABSOLUTE_TOLERANCE_IS_NEGATIVE		=  -5,
  MP_ODRF_ERROR_LOWER_BOUND_LARGER_THAN_UPPER_BOUND	=  -6,
  MP_ODRF_ERROR_ENDPOINTS_DO_NOT_STRADDLE		=  -7,
  MP_ODRF_ERROR_FUNCTION_VALUE_IS_NOT_FINITE		=  -8,
  MP_ODRF_ERROR_DERIVATIVE_IS_ZERO			=  -9,
  MP_ODRF_ERROR_FUNCTION_OR_DERIVATIVE_VALUE_INVALID	= -10,
  MP_ODRF_NEGATIVE_ERROR_CODE_LIMIT			= -11
} mp_odrf_code_t;

mp_odrf_decl const char * mp_odrf_strerror (mp_odrf_code_t code);


/** --------------------------------------------------------------------
 ** MPFR type definitions: wrapped math functions.
 ** ----------------------------------------------------------------- */

typedef int mp_odrf_mpfr_wrapped_f_t	(mpfr_ptr y, mpfr_ptr x, void * params);
typedef int mp_odrf_mpfr_wrapped_fdf_t	(mpfr_ptr dy,
					 mpfr_ptr y, mpfr_ptr x, void * params);

typedef struct {
  mp_odrf_mpfr_wrapped_f_t *		function;
  void *				params;
} mp_odrf_mpfr_function_t;

typedef struct {
  mp_odrf_mpfr_wrapped_f_t *	f;
  mp_odrf_mpfr_wrapped_f_t *	df;
  mp_odrf_mpfr_wrapped_fdf_t *	fdf;
  void *			params;
} mp_odrf_mpfr_function_fdf_t;

/* Given  the target  math function  wrapped in  a structure  F of  type
   "mp_odrf_mpfr_function_t":   compute  its   value  in   the  abscissa
   "mpfr_ptr" X and store the result in the ordinate "mpfr_ptr" Y.  Both
   X and Y must have been already initialised.

   If  an  error  occurs  computing   the  function:  a  suitable  error
   code must be returned. */
#define MP_ODRF_MPFR_FN_EVAL(F,Y,X)			\
  (((F)->function)((Y), (X), (F)->params))

/* Given  the target  math  function  and its  derivative  wrapped in  a
   structure  FDF  of  type "mp_odrf_mpfr_function_fdf_t":  compute  the
   function's value in the abscissa "mpfr_ptr" X and store the result in
   the  ordinate "mpfr_ptr"  Y.  Both  X and  Y must  have been  already
   initialised.

   If an  error occurs:  a suitable  error code must be returned. */
#define MP_ODRF_MPFR_FN_FDF_EVAL_F(FDF,Y,X)		\
  (((FDF)->f) ((Y), (X), (FDF)->params))

/* Given  the target  math  function  and its  derivative  wrapped in  a
   structure  FDF  of  type "mp_odrf_mpfr_function_fdf_t":  compute  the
   derivetive's value in the abscissa  "mpfr_ptr" X and store the result
   in  the  "mpfr_ptr"  Y.   Both  X   and  Y  must  have  been  already
   initialised.

   If an  error occurs:  a suitable  error code must be returned. */
#define MP_ODRF_MPFR_FN_FDF_EVAL_DF(FDF,DY,X)	\
  (((FDF)->df) ((DY), (X), (FDF)->params))

/* Given  the target  math  function  and its  derivative  wrapped in  a
   structure FDF of type "mp_odrf_mpfr_function_fdf_t": compute both the
   function's  value   and  the  derivative's  value   in  the  abscissa
   "mpfr_ptr" X; store  the function's value in  the ordinate "mpfr_ptr"
   Y; store the  derivative's value in the ordinate  "mpfr_ptr" DY.  All
   of X, Y and DY must have been already initialised.

   If an  error occurs:  a suitable  error code must be returned. */
#define MP_ODRF_MPFR_FN_FDF_EVAL_F_DF(FDF,DY,Y,X)	\
  (((FDF)->fdf)((DY), (Y), (X), (FDF)->params))


/** --------------------------------------------------------------------
 ** MPFR type definitions: algorithm drivers.
 ** ----------------------------------------------------------------- */

/* Prototype of function used to compute  the value of the user supplied
   math  function to  be searched  for roots.   It is  used by  the root
   bracketing algorithm drivers, the client code should never use it. */
typedef mp_odrf_code_t \
  mp_odrf_mpfr_roots_f_fun_t (void * driver_state,
			      mp_odrf_mpfr_function_t * F,
			      mpfr_ptr root,
			      mpfr_ptr x_lower, mpfr_ptr x_upper);

/* Prototype of function used to compute  the value of the user supplied
   math function, and derivative, to be  searched for roots.  It is used
   by the root polishing algorithm drivers, the client code should never
   use it. */
typedef mp_odrf_code_t \
  mp_odrf_mpfr_roots_fdf_fun_t (void * driver_state,
				mp_odrf_mpfr_function_fdf_t * FDF,
				mpfr_ptr root);

/* Prototype of function  used to initalise the state  of a root-finding
   problem.   It is  used by  the algorithm's  drivers, the  client code
   should never use it. */
typedef void mp_odrf_mpfr_roots_init_fun_t	(void * driver_state);

/* Prototype of  function used to  finalise the state of  a root-finding
   problem.   It is  used by  the algorithm's  drivers, the  client code
   should never use it. */
typedef void mp_odrf_mpfr_roots_final_fun_t	(void * driver_state);

/* Driver  for  root  bracketing  algorithms.   The  library  statically
   allocates  and  initalises  an  instance  of  this  struct  for  each
   implemented algorithm. */
typedef struct {
  const char *				name;
  size_t				driver_state_size;
  mp_odrf_mpfr_roots_init_fun_t *	init;
  mp_odrf_mpfr_roots_final_fun_t *	final;
  mp_odrf_mpfr_roots_f_fun_t *		set;
  mp_odrf_mpfr_roots_f_fun_t *		iterate;
} mp_odrf_mpfr_root_fsolver_driver_t;

/* Driver  for  root  polishing   algorithms.   The  library  statically
   allocates  and  initalises  an  instance  of  this  struct  for  each
   implemented algorithm. */
typedef struct {
  const char *				name;
  size_t				driver_state_size;
  mp_odrf_mpfr_roots_init_fun_t *	init;
  mp_odrf_mpfr_roots_final_fun_t *	final;
  mp_odrf_mpfr_roots_fdf_fun_t *	set;
  mp_odrf_mpfr_roots_fdf_fun_t *	iterate;
} mp_odrf_mpfr_root_fdfsolver_driver_t;


/** --------------------------------------------------------------------
 ** MPFR global values: predefined root-finding algorithms.
 ** ----------------------------------------------------------------- */

/* Root bracketing algorithms. */
mp_odrf_decl const mp_odrf_mpfr_root_fsolver_driver_t  * mp_odrf_mpfr_root_fsolver_bisection;
mp_odrf_decl const mp_odrf_mpfr_root_fsolver_driver_t  * mp_odrf_mpfr_root_fsolver_brent;
mp_odrf_decl const mp_odrf_mpfr_root_fsolver_driver_t  * mp_odrf_mpfr_root_fsolver_falsepos;

/* Root polishing algorithms. */
mp_odrf_decl const mp_odrf_mpfr_root_fdfsolver_driver_t  * mp_odrf_mpfr_root_fdfsolver_newton;
mp_odrf_decl const mp_odrf_mpfr_root_fdfsolver_driver_t  * mp_odrf_mpfr_root_fdfsolver_secant;
mp_odrf_decl const mp_odrf_mpfr_root_fdfsolver_driver_t  * mp_odrf_mpfr_root_fdfsolver_steffenson;


/** --------------------------------------------------------------------
 ** MPFR functions: root bracketing problems.
 ** ----------------------------------------------------------------- */

/* Root-finding computation  state using a bracketing  algorithm.  Every
   time we want to solver a root-finding problem we allocate an instance
   of this struct. */
typedef struct {
  const mp_odrf_mpfr_root_fsolver_driver_t * driver;
  mp_odrf_mpfr_function_t *		function;
  mpfr_t				root;
  mpfr_t				x_lower;
  mpfr_t				x_upper;
  void *				driver_state;
} mp_odrf_mpfr_root_fsolver_t;

/* Allocate and initialise a new root bracketing state struct to use the
   selected algorithm driver. */
mp_odrf_decl mp_odrf_mpfr_root_fsolver_t * \
  mp_odrf_mpfr_root_fsolver_alloc (const mp_odrf_mpfr_root_fsolver_driver_t * T);

/* Finalise and release a root bracketing state struct. */
mp_odrf_decl void mp_odrf_mpfr_root_fsolver_free (mp_odrf_mpfr_root_fsolver_t * S);

/* Select the  math function to be  searched for roots for  a given root
   bracketing state struct.  Also selects the search bracket. */
mp_odrf_decl mp_odrf_code_t \
  mp_odrf_mpfr_root_fsolver_set	(mp_odrf_mpfr_root_fsolver_t * S,
				 mp_odrf_mpfr_function_t * f,
				 mpfr_t x_lower, mpfr_t x_upper);

/* Perform a search iteration for a root bracketing state struct. */
mp_odrf_decl mp_odrf_code_t \
  mp_odrf_mpfr_root_fsolver_iterate (mp_odrf_mpfr_root_fsolver_t * S);

/* Inspect the current state of a root bracketing problem. */
mp_odrf_decl const char * mp_odrf_mpfr_root_fsolver_name (const mp_odrf_mpfr_root_fsolver_t * S);
mp_odrf_decl mpfr_ptr mp_odrf_mpfr_root_fsolver_root    (const mp_odrf_mpfr_root_fsolver_t * S);
mp_odrf_decl mpfr_ptr mp_odrf_mpfr_root_fsolver_x_lower (const mp_odrf_mpfr_root_fsolver_t * S);
mp_odrf_decl mpfr_ptr mp_odrf_mpfr_root_fsolver_x_upper (const mp_odrf_mpfr_root_fsolver_t * S);


/** --------------------------------------------------------------------
 ** MPFR functions: root polishing problems.
 ** ----------------------------------------------------------------- */

/* Root-finding computation  state using  a polishing  algorithm.  Every
   time we want to solver a root-finding problem we allocate an instance
   of this struct. */
typedef struct {
  const mp_odrf_mpfr_root_fdfsolver_driver_t * driver;
  mp_odrf_mpfr_function_fdf_t *		fdf;
  mpfr_t				root;
  void *				driver_state;
} mp_odrf_mpfr_root_fdfsolver_t;

/* Allocate and initialise a new root  polishing state struct to use the
   selected algorithm driver. */
mp_odrf_decl mp_odrf_mpfr_root_fdfsolver_t * \
  mp_odrf_mpfr_root_fdfsolver_alloc (const mp_odrf_mpfr_root_fdfsolver_driver_t * T);

/* Finalise and release a root polishing state struct. */
mp_odrf_decl void mp_odrf_mpfr_root_fdfsolver_free (mp_odrf_mpfr_root_fdfsolver_t * S);

/* Select the  math function to be  searched for roots for  a given root
   polishing state struct.  Also selects the initial solution guess. */
mp_odrf_decl mp_odrf_code_t \
  mp_odrf_mpfr_root_fdfsolver_set	(mp_odrf_mpfr_root_fdfsolver_t * S,
					 mp_odrf_mpfr_function_fdf_t * fdf,
					 mpfr_t root);

/* Perform a search iteration for a root polishing state struct. */
mp_odrf_decl mp_odrf_code_t \
  mp_odrf_mpfr_root_fdfsolver_iterate	(mp_odrf_mpfr_root_fdfsolver_t * S);

/* Inspect the current state of a root polishing problem. */
mp_odrf_decl const char * mp_odrf_mpfr_root_fdfsolver_name (const mp_odrf_mpfr_root_fdfsolver_t * S);
mp_odrf_decl mpfr_ptr mp_odrf_mpfr_root_fdfsolver_root (const mp_odrf_mpfr_root_fdfsolver_t * S);


/** --------------------------------------------------------------------
 ** MPFR functions: convergence tests.
 ** ----------------------------------------------------------------- */

mp_odrf_decl mp_odrf_code_t \
  mp_odrf_mpfr_root_test_interval	(mpfr_ptr x_lower, mpfr_ptr x_upper,
					 mpfr_ptr epsabs,  mpfr_ptr epsrel);

mp_odrf_decl mp_odrf_code_t \
  mp_odrf_mpfr_root_test_delta		(mpfr_ptr x1,     mpfr_ptr x0,
					 mpfr_ptr epsabs, mpfr_ptr epsrel);

mp_odrf_decl mp_odrf_code_t \
  mp_odrf_mpfr_root_test_residual	(mpfr_ptr f, mpfr_ptr epsabs);


/** --------------------------------------------------------------------
 ** MPFR functions: comparison.
 ** ----------------------------------------------------------------- */

mp_odrf_decl int mp_odrf_mpfr_fcmp (mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr epsilon);
mp_odrf_decl int mp_odrf_mpfr_absdiff_equal_p (mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr epsilon);
mp_odrf_decl int mp_odrf_mpfr_reldiff_equal_p (mpfr_srcptr a, mpfr_srcptr b, mpfr_srcptr epsilon);


/** --------------------------------------------------------------------
 ** Version functions.
 ** ----------------------------------------------------------------- */

mp_odrf_decl const char *	mp_odrf_version_string		(void);
mp_odrf_decl int	mp_odrf_version_interface_current	(void);
mp_odrf_decl int	mp_odrf_version_interface_revision	(void);
mp_odrf_decl int	mp_odrf_version_interface_age		(void);


/** --------------------------------------------------------------------
 ** Done.
 ** ----------------------------------------------------------------- */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* MP_ODRF_H */

/* end of file */
