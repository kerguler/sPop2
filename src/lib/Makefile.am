AM_CFLAGS = -I$(top_srcdir)/src/lib -Wall

AM_CFLAGS += $(GLIB_CFLAGS)
AM_CFLAGS += -lgsl
AM_CFLAGS += -lm

# The libraries to build
noinst_LIBRARIES = libspop2.a

# The files to add to the library and to the source distribution
libspop2_a_SOURCES = \
	                 $(libspop2_a_HEADERS) \
				     ran_gen.c \
				     gamma.c \
					 spop2.c
