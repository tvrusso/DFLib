#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.72])
AC_INIT([DFLib],[1.00],[russo@bogodyn.org])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([DF_Abstract_Point.hpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX([c++ CC g++ clang++  cxx])
AC_PROG_CC([cc gcc clang])
AC_PROG_RANLIB

AC_DISABLE_SHARED
LT_INIT

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_C_RESTRICT

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CHECK_LIB(proj,proj_create_argv,,AC_MSG_ERROR([DFLib requires proj libraries.]))

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
