#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.60)
AC_INIT(DFLib, 0.01, russo@bogodyn.org)
AM_INIT_AUTOMAKE(DFLib,0.01)
AC_CONFIG_SRCDIR([DF_Abstract_Point.hpp])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CHECK_LIB(proj,pj_init,,AC_MSG_ERROR([DFLib requires proj.4 libraries.]))

DFLIB_CHECK_GDAL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
