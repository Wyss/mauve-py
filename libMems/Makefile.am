ACLOCAL_AMFLAGS = -I m4
include doxygen.am
MOSTLYCLEANFILES = $(DX_CLEANFILES)

EXTRA_DIST = \
projects/libMems.doxygen \
projects/libMems.sln \
projects/libMems.vcproj

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libMems-@GENERIC_API_VERSION@.pc

SUBDIRS = libMems 

