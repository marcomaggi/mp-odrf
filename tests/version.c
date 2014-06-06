/*
  Part of: Multiple Precision One Dimensional Root Finding
  Contents: test for version functions
  Date: Fri Jun  6, 2014

  Abstract

	Test file for version functions.

  Copyright (C) 2014 Marco Maggi <marco.maggi-ipsu@poste.it>

  See the COPYING file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <mp-odrf.h>

int
main (int argc, const char *const argv[])
{
  printf("version number string: %s\n", mp_odrf_version_string());
  printf("libtool version number: %d:%d:%d\n",
	 mp_odrf_version_interface_current(),
	 mp_odrf_version_interface_revision(),
	 mp_odrf_version_interface_age());
  exit(EXIT_SUCCESS);
}

/* end of file */
