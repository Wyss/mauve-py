ACLOCAL_AMFLAGS = -I m4
include doxygen.am
MOSTLYCLEANFILES = $(DX_CLEANFILES)

EXTRA_DIST = \
projects/libGenome.sln \
projects/libGenome.vcproj \
projects/libGenome.nsi \
projects/libGenome.doxygen

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = libGenome-@GENERIC_API_VERSION@.pc

SUBDIRS = libGenome   

