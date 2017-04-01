#-*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(pocketkaldi, 0.0.1, ling032x@gmail.com)
AC_CONFIG_AUX_DIR([scripts])

AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH([kaldi],
[AS_HELP_STRING([--with-kaidi=DIR], [the root directory of Kaldi])],
[if test x$with_kaldi = xyes; then
  AC_MSG_ERROR([--with-kaidi=DIR: path of Kaldi should be specified])
fi
# Extracts ATLASLIBS from KALDI_ROOT\src\kaldi.mk
eval with_kaldi=$with_kaldi
echo "Kaldi path is: " $with_kaldi
atlaslibs=`cat $with_kaldi/src/kaldi.mk | grep -oP '(?<=^ATLASLIBS = )(.*)$'`
if test "x$atlaslibs" = "x"; then
  AC_MSG_ERROR([--with-kaidi=DIR: unable to find ATLASLIBS line in kaldi_dir/src/kaldi.mk])
fi
AC_SUBST([ATLASLIBS], [$atlaslibs])
AC_SUBST([KALDI_ROOT], [$with_kaldi])])

AC_ARG_ENABLE([tools],
[AS_HELP_STRING([--enable-tools], [enable tools for model training])],
[if test x$with_kaldi = xno || test x$with_kaldi = x; then
  AC_MSG_ERROR([--enable-tools: path of Kaldi should be specified (using --with-kaidi=DIR)])
fi
AM_CONDITIONAL(ENABLE_TOOLS, true)], [AM_CONDITIONAL(ENABLE_TOOLS, false)])

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
