#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([VOISE], [1.0.0rc10], [p.guio@ucl.ac.uk])
AC_CONFIG_SRCDIR([matlab/divideVD.m])
dnl AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_TARGET

AM_INIT_AUTOMAKE([1.11 dist-bzip2])

dnl Avoid generating makefile rules to rebuild maintainer-only files by
dnl default.  Maintainers may override this default and generate these
dnl makefile rules using the `--enable-maintainer-mode' configure option.
AM_MAINTAINER_MODE

AC_MSG_NOTICE
AC_MSG_NOTICE([Configuring $PACKAGE_STRING for $target])
AC_MSG_NOTICE


# Checks for programs.

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

dnl Generate Makefiles
AC_CONFIG_FILES([Makefile matlab/Makefile share/Makefile])

AC_CONFIG_FILES([package-webVOISE.sh:package-webVOISE.sh.in], [chmod +x package-webVOISE.sh])

AC_CONFIG_FILES([matlab/start_VOISE.m:matlab/start_VOISE.m.in], [])

AC_CONFIG_FILES([matlab/webVOISE.m:matlab/webVOISE.m.in], [])

AC_OUTPUT
