/*
  Part of: Multiple Precision One Dimensional Root Finding
  Contents: public header file
  Date: Fri Jun  6, 2014

  Abstract

	This is the public header file.   It must be included in all the
	C sources making use of the MP ODRF API.

  Copyright (C) 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

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
#  define mp_odrf_private_decl	extern
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
 ** Constants.
 ** ----------------------------------------------------------------- */




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
