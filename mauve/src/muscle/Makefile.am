ACLOCAL_AMFLAGS = -I m4
include doxygen.am
MOSTLYCLEANFILES = $(DX_CLEANFILES)

EXTRA_DIST = \
projects/libMUSCLE.vcproj \
projects/muscle.vcproj \
projects/libmuscle.doxygen


pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libMUSCLE-@GENERIC_API_VERSION@.pc

SUBDIRS = libMUSCLE  

