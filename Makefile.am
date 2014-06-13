## Process this file with automake to produce Makefile.in

ACLOCAL_AMFLAGS		= -I meta/autotools
AUTOMAKE_OPTIONS	= foreign
EXTRA_DIST		= INSTALL configure.sh
dist_doc_DATA		= README COPYING
AM_CFLAGS		= -Wall

## --------------------------------------------------------------------

EXTRA_DIST		+= meta/scripts/libmpodrf.pc.in

pkgconfigdir		= $(libdir)/pkgconfig
nodist_pkgconfig_DATA	= meta/scripts/libmpodrf.pc

## --------------------------------------------------------------------

AM_MAKEINFOFLAGS	= --no-split

info_TEXINFOS		= doc/mp-odrf.texi
doc_mp_odrf_TEXINFOS	= \
	doc/macros.texi		\
	doc/gpl-3.0.texi	\
	doc/fdl-1.3.texi

#page
#### libraries

mp_odrf_CURRENT		= @mp_odrf_VERSION_INTERFACE_CURRENT@
mp_odrf_REVISION	= @mp_odrf_VERSION_INTERFACE_REVISION@
mp_odrf_AGE		= @mp_odrf_VERSION_INTERFACE_AGE@

include_HEADERS		= src/mp-odrf.h

lib_LTLIBRARIES		= libmpodrf.la
libmpodrf_la_LDFLAGS	= -version-info $(mp_odrf_CURRENT):$(mp_odrf_REVISION):$(mp_odrf_AGE)
libmpodrf_la_SOURCES	= \
	src/mp-odrf-internals.h				\
	src/mp-odrf-version.c				\
	src/mp-odrf-errors.c				\
	src/mp-odrf-public-api.c			\
	src/mp-odrf-root-bracketing-bisection.c		\
	src/mp-odrf-root-bracketing-falsepos.c		\
	src/mp-odrf-root-bracketing-brent.c		\
	src/mp-odrf-root-polishing-newton.c		\
	src/mp-odrf-root-polishing-secant.c		\
	src/mp-odrf-root-polishing-steffenson.c

#page
#### tests

check_PROGRAMS	= \
	tests/version		\
	tests/root-bracketing	\
	tests/root-polishing	\
	tests/macros-0		\
	tests/macros-1		\
	tests/macros-2

TESTS				= $(check_PROGRAMS)

mp_odrf_tests_cppflags		= -I$(top_srcdir)/src -I$(top_srcdir)/tests
mp_odrf_tests_ldadd		= libmpodrf.la
mp_odrf_tests_sources		= tests/test.h tests/debug.h

tests_version_CPPFLAGS		= $(mp_odrf_tests_cppflags)
tests_version_LDADD		= $(mp_odrf_tests_ldadd)
tests_version_SOURCES		= tests/version.c $(mp_odrf_tests_sources)

tests_root_bracketing_CPPFLAGS	= $(mp_odrf_tests_cppflags)
tests_root_bracketing_LDADD	= $(mp_odrf_tests_ldadd)
tests_root_bracketing_SOURCES	= tests/root-bracketing.c $(mp_odrf_tests_sources)

tests_root_polishing_CPPFLAGS	= $(mp_odrf_tests_cppflags)
tests_root_polishing_LDADD	= $(mp_odrf_tests_ldadd)
tests_root_polishing_SOURCES	= tests/root-polishing.c $(mp_odrf_tests_sources)

tests_macros_0_CPPFLAGS		= $(mp_odrf_tests_cppflags)
tests_macros_0_LDADD		= $(mp_odrf_tests_ldadd)
tests_macros_0_SOURCES		= tests/macros-0.c $(mp_odrf_tests_sources)

tests_macros_1_CPPFLAGS		= $(mp_odrf_tests_cppflags)
tests_macros_1_LDADD		= $(mp_odrf_tests_ldadd)
tests_macros_1_SOURCES		= tests/macros-1.c $(mp_odrf_tests_sources)

tests_macros_2_CPPFLAGS		= $(mp_odrf_tests_cppflags)
tests_macros_2_LDADD		= $(mp_odrf_tests_ldadd)
tests_macros_2_SOURCES		= tests/macros-2.c $(mp_odrf_tests_sources)

## --------------------------------------------------------------------

installcheck-local: $(TESTS)
	@for f in $(TESTS); do $(builddir)/$$f; done

### end of file
