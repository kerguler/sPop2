AM_CFLAGS = -I$(top_srcdir)/src/lib -Wall

AM_CFLAGS += $(GLIB_CFLAGS)

# The libraries to build
lib_LTLIBRARIES =  libspop2.la libspop2py.la

# The header files to include
libspop2_ladir = $(includedir)/spop2
libspop2py_ladir = $(includedir)/spop2

# The files to add to the library and to the source distribution
# ---
libspop2_la_HEADERS = uthash.h spop2.h

libspop2_la_SOURCES = \
	                 $(libspop2_la_HEADERS) \
				     ran_gen.c \
				     gamma.c \
				     qnta.c \
					 spop2.c

libspop2_la_LIBADD = $(GLIB_LIBS)

# ---
libspop2py_la_HEADERS = uthash.h spop2.h

libspop2py_la_SOURCES = \
	                 $(libspop2py_la_HEADERS) \
				     ran_gen.c \
				     gamma.c \
				     qnta.c \
					 spop2.c

libspop2py_la_LIBADD = $(GLIB_LIBS)
libspop2py_la_LDFLAGS = -lgsl -lgslcblas -module -avoid-version -shared
