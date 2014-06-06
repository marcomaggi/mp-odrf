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
	src/mp-odrf-internals.h		\
	src/mp-odrf-version.c

#page
#### tests

check_PROGRAMS	= \
	tests/version

TESTS			= $(check_PROGRAMS)

mp_odrf_tests_cppflags	= -I$(top_srcdir)/src
mp_odrf_tests_ldadd	= libmpodrf.la

tests_version_CPPFLAGS	= $(mp_odrf_tests_cppflags)
tests_version_LDADD	= $(mp_odrf_tests_ldadd)

installcheck-local: $(TESTS)
	@for f in $(TESTS); do $(builddir)/$$f; done

### end of file
