dnl Process this file with autoconf to produce a configure script.
AC_PREREQ([2.59])
AC_INIT(libMems/Match.h)
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR([m4])

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------

GENERIC_LIBRARY_NAME=libMems

#release versioning
GENERIC_MAJOR_VERSION=1
GENERIC_MINOR_VERSION=6
GENERIC_MICRO_VERSION=0

#API version (often = GENERIC_MAJOR_VERSION.GENERIC_MINOR_VERSION)
GENERIC_API_VERSION=1.6
AC_SUBST(GENERIC_API_VERSION)

#shared library versioning
GENERIC_LIBRARY_VERSION=1:0:0
#                       | | |
#                +------+ | +---+
#                |        |     |
#             current:revision:age
#                |        |     |
#                |        |     +- increment if interfaces have been added
#                |        |        set to zero if interfaces have been removed
#                                  or changed
#                |        +- increment if source code has changed
#                |           set to zero if current is incremented
#                +- increment if interfaces have been added, removed or changed



dnl --------------------------------
dnl Package name and version number
dnl --------------------------------

AC_SUBST(GENERIC_LIBRARY_VERSION)

PACKAGE=$GENERIC_LIBRARY_NAME
AC_SUBST(GENERIC_LIBRARY_NAME)

GENERIC_VERSION=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION.$GENERIC_MICRO_VERSION
GENERIC_RELEASE=$GENERIC_MAJOR_VERSION.$GENERIC_MINOR_VERSION
AC_SUBST(GENERIC_RELEASE)
AC_SUBST(GENERIC_VERSION)

VERSION=$GENERIC_VERSION

AM_INIT_AUTOMAKE($PACKAGE, $VERSION, no-define)

dnl Override default O2
CFLAGS=${CFLAGS-""}
CXXFLAGS=${CXXFLAGS-""}

AC_PREFIX_DEFAULT(/usr/local)

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S 
AM_PROG_LIBTOOL
AC_SYS_LARGEFILE

dnl Checks for header files.
AC_HEADER_STDC

dnl Check what compiler we're using
AM_CONDITIONAL(ICC, test x$CXX = xicc )
EXTRA_CXX_FLAGS=""
if( test x$CC == "xgcc" ) then
	EXTRA_CXX_FLAGS="-Wno-deprecated"
fi
AC_SUBST(EXTRA_CXX_FLAGS)

dnl Allow debugging compilation
AC_ARG_ENABLE(debug,
[  --enable-debug    Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

dnl Get location of Boost
AX_BOOST_BASE(1.34.0)
AC_CXX_NAMESPACES
AX_BOOST_FILESYSTEM
AX_BOOST_PROGRAM_OPTIONS
AX_BOOST_IOSTREAMS
AX_BOOST_SYSTEM

dnl Get location of libGenome Headers
PKG_CHECK_MODULES(DEPS, libGenome-1.3 >= 1.3.1  libMUSCLE-3.7 >= 1.0.0)
AC_SUBST(DEPS_CFLAGS)

dnl Check for OpenMP
#AX_OPENMP()
AC_SUBST(OPENMP_CFLAGS)
dnl OpenMP checker only defines for C when compiling both C and C++
OPENMP_CXXFLAGS=$OPENMP_CFLAGS
AC_SUBST(OPENMP_CXXFLAGS)

dnl ensure portability for OS X with these checks
AC_CHECK_HEADERS(sys/types.h)
AC_CHECK_HEADERS(sys/aio.h)
AC_CHECK_HEADERS(aio.h)
AC_CHECK_HEADERS(features.h)

dnl certain parts of the library need async io and threads
AC_CHECK_LIB(pthread, pthread_mutex_unlock)
AC_CHECK_LIB(rt, aio_write)
AC_SUBST(DEPS_LIBS)

dnl Make doxygen docs
DX_INIT_DOXYGEN( "libMems", "projects/libMems.doxygen", "doc" )

AM_CONFIG_HEADER(config.h)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
dnl AC_C_BIGENDIAN
AC_HEADER_TIME

dnl Checks for library functions.
AC_PROG_GCC_TRADITIONAL

dnl SAVE_LIBRARY_VERSION
AC_SUBST(LIBTOOL_VERSION_INFO)

AC_OUTPUT(Makefile libMems/Makefile libMems-1.6.pc )
#doc/html/Makefile  doc/man/Makefile  doc/man/man3/Makefile)
