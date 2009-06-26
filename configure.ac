AC_INIT(RNAcode, 0.2beta, wash@tbi.univie.ac.at, RNAcode)
AM_INIT_AUTOMAKE

AC_PROG_CC


if test "$GCC" = yes; then
  AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
  AC_DEFINE(UNUSED,, only for gcc)
fi
      
AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)


AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h strings.h unistd.h)
AC_C_CONST
AC_TYPE_SIZE_T
AC_INLINE

AC_CHECK_FUNCS(strdup strstr strchr erand48)

AC_PROG_CXX

AC_PROG_RANLIB

AC_PROG_INSTALL

AC_CONFIG_HEADERS(config.h) 

AC_OUTPUT([Makefile phyml/Makefile levmar/Makefile seqgen/Makefile src/Makefile librna/Makefile])


