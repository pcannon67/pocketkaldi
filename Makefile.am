ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = gnu subdir-objects foreign
AM_CXXFLAGS = -I$(top_srcdir)/src -DMODEL_DIR=\"$(pkgdatadir)/\" -fno-rtti \
              -fno-exceptions -fPIC

bin_PROGRAMS = main
main_SOURCES = src/main.c
