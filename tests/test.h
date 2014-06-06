/*
   Part of: Multiple Precision One-Dimensional Root-Finding
   Contents: test header file
   Date: Sun Jun  8, 2008

   Abstract



   Copyright (c) 2008, 2009, 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

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

#ifndef TEST_H
#define TEST_H 1


/** ------------------------------------------------------------
 ** Headers.
 ** ----------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <mpfr.h>


/** ------------------------------------------------------------
 ** Definitions.
 ** ----------------------------------------------------------*/

/*
  "test_result == 1": an error occurred
  "test_result == 0": success
 */
int		test_result	= 0;
const char *	test_name	= NULL;

/* ------------------------------------------------------------ */

#define title(...)						\
  do {								\
    fprintf(stderr,"\n****** ");				\
    fprintf(stderr,__VA_ARGS__);				\
    fprintf(stderr,"\n");					\
  } while(0);

#define subtitle(...)						\
  do {								\
    fprintf(stderr,"\n*** ");					\
    fprintf(stderr,__VA_ARGS__);				\
    fprintf(stderr,"\n");					\
  } while(0);

/* ------------------------------------------------------------------ */

#if (DEBUGGING == 1)
#  define start(NAME,DESCR)	{				\
    test_result=0;						\
    test_name=(NAME);						\
    report("%s: %-50s... \n", test_name, DESCR);		\
    }
#else
#  define start(NAME,DESCR)	{				\
    test_result=0;						\
    test_name=(NAME);						\
    report("%s: %-50s... ", test_name, DESCR);			\
  }
#endif

#define fine()		{					\
    if (test_result) {						\
      report("%s (%s): %s\n", test_name, __func__,"*** FAILED");\
    } else {							\
      mpfr_fprintf(stderr, "fine\n");				\
    }								\
  }

/* ------------------------------------------------------------------ */

#define report(...)		mpfr_fprintf(stderr,__VA_ARGS__)

#define error(...)						\
  do {								\
    report(__VA_ARGS__);					\
    test_result=1;						\
  } while(0);


#define validate(EXPR,...)	{				\
    if (!(EXPR)) {						\
      report("\n\terror: ");					\
      report(__VA_ARGS__);					\
      report("\n");						\
      test_result=1;						\
    }								\
  }

#define error_if_true(EXPR,...)		\
  do {					\
    if (EXPR) {				\
      report("\n\terror: ");		\
      report(__VA_ARGS__);		\
      report("\n");			\
      test_result=1;			\
    }					\
  } while(0)

#define error_if_false(EXPR,...)	\
  do {					\
    if (!(EXPR)) {			\
      report("\n\terror: ");		\
      report(__VA_ARGS__);		\
      report("\n");			\
      test_result=1;			\
    }					\
  } while(0)

#define validate_expected_got(EXPECTED,GOT)			\
  do {								\
    mpfr_t	expected_result, result_tolerance;		\
    mpfr_init(expected_result);					\
    mpfr_init(result_tolerance);				\
    mpfr_set_d(expected_result, EXPECTED, GMP_RNDN);		\
    mpfr_set_d(result_tolerance, 1.0, GMP_RNDN);		\
    validate(mp_odrf_mpfr_absdiff_equal_p(expected_result, GOT,	\
				     result_tolerance),		\
	     "expected %.20RNf, got %.20RNf",			\
	     expected_result, GOT);				\
    mpfr_clear(expected_result);				\
    mpfr_clear(result_tolerance);				\
  } while(0);

#define validate_expected_got_mpfr(EXPECTED,GOT)		\
  do {								\
    mpfr_t	result_tolerance;				\
    mpfr_init(result_tolerance);				\
    mpfr_set_d(result_tolerance, 1e-6, GMP_RNDN);		\
    validate((!mpfr_nan_p(GOT)) &&				\
	     (mp_odrf_mpfr_absdiff_equal_p(EXPECTED, GOT,	\
				    result_tolerance)),		\
	     "expected %.20RNf, got %.20RNf",			\
	     EXPECTED, GOT);					\
    mpfr_clear(result_tolerance);				\
  } while(0);

#endif /* TEST_H */

/* end of file */
