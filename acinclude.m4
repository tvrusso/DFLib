AC_DEFUN([DFLIB_CHECK_GDAL],
[
use_gdal=no
#
# Important: DO NOT use "use_gdal" as the variable here, because AC_CHECK_PROG
# will do nothing if the variable is already set!
#
AC_CHECK_PROG(found_gdal_config, [gdal-config], yes, no)
if test "${found_gdal_config}" = "yes"; then
   save_cppflags="$CPPFLAGS" 
   save_libs="$LIBS" 
   save_ldflags="$LDFLAGS" 

   GDAL_BIN="gdal-config"
   CPPFLAGS="$CPPFLAGS `${GDAL_BIN} --cflags`" 
#
# This is an annoyance: gdal-config --libs outputs both LDFLAGS and LIBS 
# stuff.  AC_CHECK_LIB puts the -l in if it works, and we only want the LDFLAGS
#   LIBS="$LIBS `${GDAL_BIN} --libs`"
# Remove the -lgdal from what gdal-config --libs returns, because AC_CHECK_LIB
# will put it into LIBS for us.
#
   LDFLAGS="$LDFLAGS `${GDAL_BIN} --libs | sed -e 's/-lgdal[^ ]*//'`"
   AC_CHECK_HEADERS(gdal.h, [AC_CHECK_LIB(gdal, GDALAllRegister,
                    [use_gdal="yes"
                     LIBS="$LIBS -lgdal"
                     AC_DEFINE(HAVE_LIBGDAL, , 
                      [Define to 1 if you have the `gdal' library (-lgdal).])],
                    [CPPFLAGS=${save_cppflags}
                     LDFLAGS=${save_ldflags}
                     LIBS=${save_libs}])])
else
   AC_MSG_NOTICE([Cannot find gdal-config:  Checking standard locations.])
   AC_CHECK_HEADERS(gdal.h, [AC_CHECK_LIB(gdal, GDALAllRegister,
                    [use_gdal="yes"
                     LIBS="$LIBS -lgdal"
                     AC_DEFINE(HAVE_LIBGDAL, , 
                      [Define to 1 if you have the `gdal' library (-lgdal).])],)])
fi
]
)
