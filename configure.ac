dnl @configure_input@
dnl

AC_PREREQ([2.69])
AC_INIT([Multiple Precision One-Dimensional Root-Finding],[0.1d0],
  [marco.maggi-ipsu@poste.it],
  [mp-odrf],[https://github.com/marcomaggi/mp-odrf/])
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_MACRO_DIR([meta/autotools])
AC_CONFIG_AUX_DIR([meta/autotools])
AC_CANONICAL_BUILD
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([foreign dist-xz no-dist-gzip subdir-objects])
AM_MAINTAINER_MODE

AM_PROG_AR
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_MKDIR_P

LT_PREREQ([2.4])
LT_INIT

dnl page
#### libraries interface version

dnl This is the version stored in the pkg-config data file.
AC_SUBST([MP_ODRF_PKG_CONFIG_VERSION],[0.1])

MM_LIBTOOL_LIBRARY_VERSIONS([mp_odrf],0,0,0)

dnl page
#### basic system inspection

AC_LANG([C])
dnl AC_DEFINE([_GNU_SOURCE],[1],[enable everything GNU])
AC_PROG_CC_C99
AM_PROG_AS
AC_PROG_CC
AM_PROG_CC_C_O
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([assert.h stdint.h strings.h errno.h stddef.h])
AC_CHECK_SIZEOF(void *)
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

AC_CACHE_SAVE

dnl page
#### checks for library functions

AC_FUNC_MALLOC
AC_FUNC_MEMCMP
dnl AC_CHECK_FUNCS([memmove memset strerror])

AC_CACHE_SAVE

dnl page
#### external libraries

AC_CHECK_HEADERS([gmp.h],,
  [AC_MSG_ERROR([Cannot find libgmp.])])
AC_CHECK_LIB([gmp],[__gmpn_add],,
  [AC_MSG_ERROR([Cannot find libgmp.])])

AC_CHECK_HEADERS([mpfr.h],,
  [AC_MSG_FAILURE([missing mpfr.h header],[2])])
AC_CHECK_LIB([mpfr],[mpfr_init],,
  [AC_MSG_FAILURE([missing or unusable MPFR library],[2])])

AC_CACHE_SAVE

dnl page
#### operating system specific configuration

dnl case "$target_os" in
dnl   *linux*)
dnl     ;;
dnl   *bsd*)
dnl     ;;
dnl   *cygwin*)
dnl     ;;
dnl   *darwin*)
dnl     cctemplate_OS_DARWIN=yes
dnl     ;;
dnl esac
dnl AM_CONDITIONAL([WANT_LINUX], [test "x$cctemplate_OS_LINUX"  = xyes])
dnl AM_CONDITIONAL([WANT_BSD],   [test "x$cctemplate_OS_BSD"    = xyes])
dnl AM_CONDITIONAL([WANT_CYGWIN],[test "x$cctemplate_OS_CYGWIN" = xyes])
dnl AM_CONDITIONAL([WANT_DARWIN],[test "x$cctemplate_OS_DARWIN" = xyes])

dnl page
#### finish

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile]
   [meta/scripts/libmpodrf.pc]
   [meta/slackware/slack-desc])
AC_OUTPUT

### end of file
# Local Variables:
# mode: autoconf
# page-delimiter: "^dnl page"
# End:
