# Based on the configure.ac file from the gsl R package.

AC_INIT([vines], [0.1.0])

AC_PATH_PROG([GSL_CONFIG], [gsl-config])
if test -z "${GSL_CONFIG}"; then
   AC_MSG_ERROR([could not find GSL configuration tool])
fi

AC_MSG_CHECKING([for required GSL version])
AS_VERSION_COMPARE(`${GSL_CONFIG} --version`, [1.14],
                   AC_MSG_ERROR([could not find required version of GSL]),
                   AC_MSG_RESULT([yes]),
                   AC_MSG_RESULT([yes]))

GSL_CFLAGS=`${GSL_CONFIG} --cflags`
GSL_LIBS=`${GSL_CONFIG} --libs`

AC_SUBST([GSL_CFLAGS])
AC_SUBST([GSL_LIBS])

AC_CONFIG_FILES([src/Makevars])

AC_OUTPUT()
