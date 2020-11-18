AM_CFLAGS = -I$(top_srcdir)/src/bin -Wall

AM_CFLAGS += $(GLIB_CFLAGS)
AM_CFLAGS += -lgsl
AM_CFLAGS += -lm

bin_PROGRAMS = survival
survival_SOURCES = survival.c

survival_LDADD = $(top_builddir)/src/lib/libspop2.a $(AM_LDFLAGS)
survival_LDADD += $(GLIB_LIBS)
