/*
   Part of: Multiple Precision One Dimensional Root Finding
   Contents: root bracketing bisection algorithm
   Date: Sat Mar  7, 2009

   Abstract

	This  module  implements  bisection  root  bracketing  algorithm
	driver.

   Copyright (c) 2009, 2014 Marco Maggi <marco.maggi-ipsu@poste.it>
   Copyright (c)  1996, 1997, 1998,  1999, 2000, 2007  Reid Priedhorsky,
   Brian Gough.

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

#include "mp-odrf.h"

typedef struct {
  const mp_odrf_mpfr_root_fsolver_driver_t *	driver;
  double			x_lower;
  double			x_upper;
  mp_odrf_mpfr_function_fun_t *	function;
  const char *			description;
} bracket_meta_data_t;




/* end of file */
