AM_CFLAGS = -I$(top_srcdir)/src/lib -Wall

AM_CFLAGS += $(GLIB_CFLAGS)

# The header files to include
libspop2_ladir = $(prefix)/spop2

# The libraries to build
lib_LTLIBRARIES =  libspop2.la

# The files to add to the library and to the source distribution
libspop2_la_HEADERS = gamma.h ran_gen.h uthash.h spop2.h

libspop2_la_SOURCES = \
	                 $(libspop2_la_HEADERS) \
				     ran_gen.c \
				     gamma.c \
					 spop2.c

libspop2_la_LIBADD = $(GLIB_LIBS)
