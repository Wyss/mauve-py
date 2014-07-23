AM_CFLAGS = -DUSE_POSIX_AIO

DMSML_H = \
asyncio.h alinuxaio.h aPOSIXaio.h \
alibc.h awin32aio.h buffer.h \
util.h sorting.h dmsort.h \
timing.h sml.h

DMSML_SRC = \
asyncio.c alinuxaio.c aPOSIXaio.c \
alibc.c awin32aio.c buffer.c \
util.c sorting.c dmsort.c \
timing.c sml.c

library_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)/dmSML

library_include_HEADERS = $(DMSML_H)

noinst_LTLIBRARIES = libdmSML.la
libdmSML_la_SOURCES = $(DMSML_SRC)

INCLUDES = -I$(top_srcdir) $(DEPS_CFLAGS)
