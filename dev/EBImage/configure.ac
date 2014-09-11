MagickConfig='/vol/vipdata/users/kpk09/imagemagick/bin/Magick-config'; export MagickConfig
WandConfig='/vol/vipdata/users/kpk09/imagemagick/bin/Wand-config'; export WandConfig

AC_INIT([EBImage], 2.2.0)

AC_CHECK_PROG(magickfound, [$MagickConfig --cppflags], yes, no)
if test "$magickfound" = "yes"; then
    echo "* found ImageMagick version `$MagickConfig --version`"

    ## check version number
    $MagickConfig --version | head -n 1 | awk '{
            x[1] = substr($1,1,1);
            x[2] = substr($1,3,1);
            x[3] = substr($1,5,1);
            if ( x[1] < 6 ) exit(1);
            if ( x[1] > 6 ) exit(0);
            if ( x[2] < 3 ) exit(2);
            if ( x[2] > 3 ) exit(0);
            if ( x[3] < 7 ) exit(3);
            exit(0);
          }'
    ext=$?
    if [[ $ext -ge 1 ]] ; then
      AC_MSG_ERROR("Requires ImageMagick version >= 6.3.7")
    fi

    ## set up include and library paths
    echo "* found ImageMagick includes in `$MagickConfig --cppflags`"
    COMPILERFLAGS="`$MagickConfig --cflags`"
    COMPILERFLAGS="${COMPILERFLAGS} `$MagickConfig --cppflags`"
    CFLAGS="-Wall ${COMPILERFLAGS} ${CFLAGS}"
    CXXFLAGS="-Wall ${COMPILERFLAGS} ${CXXFLAGS}"
    LINKERFLAGS="`$MagickConfig --ldflags`"
    LINKERFLAGS="${LINKERFLAGS} `$MagickConfig --libs`"
    LINKERFLAGS="${LINKERFLAGS} `$WandConfig --libs`"
    LIBS="${LINKERFLAGS} ${LIBS} ${LFLAGS}"

    ## remove -lgomp, which causes a seg fault on openSUSE 11.2|3 when dynamically linked together with -lpthread
    LIBS=`echo $LIBS | awk '{gsub(/\-lgomp/,"")}; 1'`
else
    AC_MSG_ERROR("ImageMagick or ImageMagick includes not found!")
fi

PKGCONF=pkg-config
AC_CHECK_PROG(pkgfound, [pkg-config --version], yes, no)
if test "$pkgfound" = "no"; then
  PKGCONF=pkgconfig
  AC_CHECK_PROG(pkgfound, [pkgconfig --version], yes, no)
fi

if test "$pkgfound" = "yes"; then
    if `$PKGCONF --cflags gtk+-2.0 > /dev/null` ; then
        echo "* compiling with GTK+ support"
        COMPILERFLAGS="`$PKGCONF --cflags gtk+-2.0`"
        CFLAGS="-DUSE_GTK -DGLIB_GETTEXT ${COMPILERFLAGS} ${CFLAGS}"
        ## DO NOT ADD GTK keys to c++, that code does not use it
        ## CXXFLAGS="-DUSE_GTK -DGLIB_GETTEXT ${COMPILERFLAGS} ${CXXFLAGS}"
        LINKERFLAGS="`$PKGCONF --libs gtk+-2.0`"
        LIBS="${LINKERFLAGS} ${LIBS} ${LFLAGS}"
        LFLAGS="${LIBS}"
    else
        echo "* GTK+2.0 not found: compiling without GTK+ support"
    fi
else
    echo "* pkg-config/pkgconfig not found: compiling without GTK+ support"
fi

if test -z "${R_HOME}"; then
    AC_MSG_ERROR("Cannot determine R_HOME...")
fi

echo -----------------------------------------------------

AC_SUBST(CFLAGS)
AC_SUBST(CXXFLAGS)
AC_SUBST(LIBS)

AC_OUTPUT(src/Makevars)
